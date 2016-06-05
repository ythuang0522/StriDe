//-----------------------------------------------
// Copyright 2016 National Chung Cheng University
// Written by Yao-Ting Huang & Ping-Yeh Chen
// Released under the GPL
//-----------------------------------------------
//
// PacBioHybridCorrectionProcess - Hybrid correction using FM-index walk for PacBio reads
//

#include "PacBioHybridCorrectionProcess.h"
#include "CorrectionThresholds.h"
#include "HashMap.h"
#include <iomanip>
#include "SAIPBHybridCTree.h"
#include "SAIPBSelfCTree.h"
#include "Timer.h"
#include "multiple_alignment.h"
#include "LongReadOverlap.h"

using namespace std;


PacBioHybridCorrectionProcess::PacBioHybridCorrectionProcess(const PacBioHybridCorrectionParameters params):m_params(params)
{
}

PacBioHybridCorrectionProcess::~PacBioHybridCorrectionProcess()
{
}

// PacBio Hybrid Correction by Ya, v20160216.
PacBioHybridCorrectionResult PacBioHybridCorrectionProcess::PBHybridCorrection(const SequenceWorkItem& workItem)
{
	PacBioHybridCorrectionResult result;
	
	// std::cout << workItem.read.id <<"\n";
	std::vector<SeedFeature> seedVec, pacbioCorrectedStrs;
	std::string readSeq = workItem.read.seq.toString();
	// seedVec = seedingByDynamicKmer(readSeq);
	// seedVec = seedingByDynamicKmer_v2(readSeq);
	seedVec = seedingByDynamicKmer_v3(readSeq);
	
	if(seedVec.size() >= 2)
	{
		result.correctedLen += seedVec.at(0).seedLength;
		pacbioCorrectedStrs.push_back(seedVec.at(0));
	}
	else
	{
		result.merge = false;
		return result;
	}

	// FMWalk for each pair of seeds
	for(size_t targetSeed = 1 ; targetSeed < seedVec.size() ; targetSeed++)
	{
		SeedFeature preTarget = seedVec.at(targetSeed-1);
		SeedFeature source = pacbioCorrectedStrs.back();
		SeedFeature target = seedVec.at(targetSeed);
		int dis_between_src_target = target.seedStartPos - preTarget.seedEndPos - 1;
		// we need raw pacbio string between source and target seed and 10 bp around
		// to do global alignmet with fmwalk result.
		string strBetweenSrcTarget = readSeq.substr(preTarget.seedEndPos+1-10, dis_between_src_target+20);
		int FMWalkReturnType;
		FMWalkResult FMWResult;
		FMWalkReturnType = extendBetweenSeeds(source, target, strBetweenSrcTarget, dis_between_src_target, &FMWResult, targetSeed);

		// debug
		// if(targetSeed < seedVec.size())
		// {
			// string kmer = seedVec.at(targetSeed-1).seedStr;
			// int fwdKmerFreqs = BWTAlgorithms::countSequenceOccurrencesSingleStrand(kmer, m_params.indices);
			// int rvcKmerFreqs = BWTAlgorithms::countSequenceOccurrencesSingleStrand(reverseComplement(kmer), m_params.indices);
			// int kmerFreqs = fwdKmerFreqs+rvcKmerFreqs;
			// std::cout << ">" << seedVec.at(targetSeed-1).seedStartPos << "."
			// << kmerFreqs << endl
			// << seedVec.at(targetSeed-1).seedStr << endl;
			// kmer = seedVec.at(targetSeed).seedStr;
			// fwdKmerFreqs = BWTAlgorithms::countSequenceOccurrencesSingleStrand(kmer, m_params.indices);
			// rvcKmerFreqs = BWTAlgorithms::countSequenceOccurrencesSingleStrand(reverseComplement(kmer), m_params.indices);
			// kmerFreqs = fwdKmerFreqs+rvcKmerFreqs;
			// std::cout << ">" << seedVec.at(targetSeed).seedStartPos << "."
			// << kmerFreqs << endl
			// << seedVec.at(targetSeed).seedStr << endl;
			// cout << FMWalkReturnType << "----\n";
		// }
		
		// debug
		// std::cout << targetSeed << " " << preTarget.seedStartPos << " " << FMWalkReturnType << "\n";
		
		// debug
		// if(FMWResult.mergedSeq.length() != 0)
		//	if(targetSeed == 34)
		//		std::cout << ">" << targetSeed << "." << seedVec.at(targetSeed-1).seedStartPos << "." << FMWResult.mergedSeq.substr(source.seedLength - seedVec.at(targetSeed-1).seedLength).length() << "\n" << FMWResult.mergedSeq.substr(source.seedLength - seedVec.at(targetSeed-1).seedLength) << "\n";
		
		// debug
		// if(targetSeed == 4)
		// {
			// std::cout << ">" << targetSeed << "." << seedVec.at(targetSeed-1).seedStartPos << "." << readSeq.substr(seedVec.at(targetSeed-1).seedStartPos,seedVec.at(targetSeed).seedEndPos-seedVec.at(targetSeed-1).seedStartPos+1).length() << "\n" << readSeq.substr(seedVec.at(targetSeed-1).seedStartPos,seedVec.at(targetSeed).seedEndPos-seedVec.at(targetSeed-1).seedStartPos+1) << "\n";
			// std::cout << ">" << targetSeed << "." << seedVec.at(targetSeed-1).seedStartPos << "." << mergedseq.substr(source.seedLength - seedVec.at(targetSeed-1).seedLength).length() << "\n" << mergedseq.substr(source.seedLength - seedVec.at(targetSeed-1).seedLength) << "\n";
			// std::cout << seedVec.at(targetSeed-1).seedStr << "\n";
			// std::cout << target.seedStr << "\n";
		// }
		
		// debug
		// if(targetSeed == 10 || targetSeed == 9)
		// {
			// for(int i=0;i<=(target.seedLength-21);i++)
			// {
				// string kmer = "GAAAGCAGCGCAGCAACCCAGAT";//target.seedStr.substr(i, 21);
				// int fwdKmerFreqs = BWTAlgorithms::countSequenceOccurrencesSingleStrand(kmer, m_params.indices);
				// int rvcKmerFreqs = BWTAlgorithms::countSequenceOccurrencesSingleStrand(reverseComplement(kmer), m_params.indices);
				// int kmerFreqs = fwdKmerFreqs+rvcKmerFreqs;
				// std::cout << kmer << " " << kmerFreqs << ".21" << endl;
			// }
		// }
		
		// record corrected pacbio reads string
		// FMWalk success
		if(FMWalkReturnType == 1)
		{
			size_t gainPos = source.seedLength;
			// have gain ground
			if(FMWResult.mergedSeq.length() > gainPos)
			{
				string gainStr = FMWResult.mergedSeq.substr(gainPos);
				pacbioCorrectedStrs.back().append(gainStr);
				pacbioCorrectedStrs.back().isRepeat = target.isRepeat;
				pacbioCorrectedStrs.back().isPBSeed = target.isPBSeed;
				pacbioCorrectedStrs.back().isPBSeed = target.isNextRepeat;
				pacbioCorrectedStrs.back().startBestKmerSize = target.startBestKmerSize;
				pacbioCorrectedStrs.back().endBestKmerSize = target.endBestKmerSize;
				pacbioCorrectedStrs.back().seedEndPos = target.seedEndPos;
				pacbioCorrectedStrs.back().seedStartPos = target.seedStartPos;

				result.correctedLen += gainStr.length();
			}
		}
		// FMWalk failure: 
		// 1. high error 
		// 2. exceed leaves
		// 3. exceed depth
		else
		{
			//std::cout << FMWalkReturnType << ":\t" << source.seedStr << "\t" << target.seedStr << "\t" << dis_between_src_target << "\n";	
			pacbioCorrectedStrs.push_back(target);
			result.correctedLen += target.seedLength;
		}
		
		// output information
		result.totalWalkNum++;			
		result.seedDis += dis_between_src_target;
		if(FMWalkReturnType == 1)
			result.correctedNum++;
		// else if(FMWalkReturnType == -1)
			// result.highErrorNum++;
		// else if(FMWalkReturnType == -2)
			// result.exceedDepthNum++;
		// else if(FMWalkReturnType == -3)
			// result.exceedLeaveNum++;
	}

	result.totalSeedNum = seedVec.size();
	result.totalReadsLen = readSeq.length();
	result.merge = true;
	for(size_t result_count = 0 ; result_count < pacbioCorrectedStrs.size() ; result_count++)
		result.correctedPacbioStrs.push_back(pacbioCorrectedStrs[result_count].seedStr);
	
	return result;
}

// PacBio Hybrid Correction Seeding By Dynamic Kmer, v20160217 by Ya.
// find seeds by dynamic kmers, which is suitable for the PacBio hybrid error correction,
// where repeat regions require large kmers and error-prone regions require small kmers.
std::vector<SeedFeature> PacBioHybridCorrectionProcess::seedingByDynamicKmer(const string& readSeq)
{
	std::vector<SeedFeature> seedVec;
	std::vector<int> seedEndPosVec;
	size_t maxKmerSize = m_params.kmerLength;
	size_t minKmerSize = m_params.minKmerLength;
	size_t kmerThreshold = m_params.seedKmerThreshold;
	
	// prevention of short reads
	if(readSeq.length() < maxKmerSize)
		return seedVec;
	
	// computing maximum interval of variously consecutive matches length (dk figure)
	std::vector<int> maxSeedInterval;
	for(size_t i = 0 ; i <= maxKmerSize ; i++)
		maxSeedInterval.push_back(2*3.8649*pow(2.7183,0.1239*i));
	
	// set dynamic kmer as largest kmer initially, 
	// which will reduce size when no seeds can be found within a maximum interval.
	size_t dynamicKmerSize = maxKmerSize;
	size_t dynamicKmerThreshold = kmerThreshold;
		
	// search for solid kmers and group consecutive solids kmers into one seed
	// reduce kmer size if no seeds can be found within maxSeedInterval
	for(size_t i = 0 ; i+dynamicKmerSize <= readSeq.length() ; i++)
	{
		string kmer = readSeq.substr(i, dynamicKmerSize);
		size_t fwdKmerFreqs = BWTAlgorithms::countSequenceOccurrencesSingleStrand(kmer, m_params.indices);
		size_t rvcKmerFreqs = BWTAlgorithms::countSequenceOccurrencesSingleStrand(reverseComplement(kmer), m_params.indices);
		size_t kmerFreqs = fwdKmerFreqs+rvcKmerFreqs;
		
		// std::cout << i << ":\t" << kmer << "\t" << kmerFreqs << "\n";
		// finding the seed
		if(kmerFreqs >= dynamicKmerThreshold)
		{
			size_t seedStartPos = i;
			size_t maxKmerFreq = kmerFreqs;
			
			// group consecutive solid kmers into one seed if possible
			for(i++ ; i+dynamicKmerSize <= readSeq.length() ; i++)
			{
				kmer = readSeq.substr(i, dynamicKmerSize);
				fwdKmerFreqs = BWTAlgorithms::countSequenceOccurrencesSingleStrand(kmer, m_params.indices);
				rvcKmerFreqs = BWTAlgorithms::countSequenceOccurrencesSingleStrand(reverseComplement(kmer), m_params.indices);
				kmerFreqs = fwdKmerFreqs+rvcKmerFreqs;
				maxKmerFreq = std::max(maxKmerFreq,kmerFreqs);

				// std::cout << i << ":\t" << kmer << "\t" << kmerFreqs << "\n";

				if(kmerFreqs < dynamicKmerThreshold)
					break;
			}
			
			// small-sized seed has 50% chance of errors in C elegans, 
			// skip only if there is another seed nearby 30bp
			if( (i-seedStartPos) < 2 && !seedVec.empty() && i-seedVec.back().seedEndPos <= 30) continue;

			size_t seedEndPos = i+dynamicKmerSize-2;

			// repeat seeds are less accurate at boundary and should be trimmed 
			if(maxKmerFreq > m_params.coverage*4)
				trimRepeatSeed(readSeq, m_params.coverage, seedStartPos, seedEndPos);

			// super repeat seeds with frequency > 2000 are troublesome, often lead to -3 but no good solution so far, mark first
			bool isSuperRepeat = maxKmerFreq > m_params.coverage*15?true:false;
			SeedFeature newSeed(seedStartPos, readSeq.substr(seedStartPos, seedEndPos-seedStartPos+1), isSuperRepeat, dynamicKmerSize, kmerThreshold*4);

			// skip low-complexity sequencing errors of PacBio
			// bool isShortAndHighFreq = i-seedStartPos <= 2 && maxKmerFreq > 80;
			if(!isLowComplexity(newSeed.seedStr, 0.9))
				seedVec.push_back(newSeed);
			
			seedEndPosVec.push_back(seedEndPos);
			
			// debug
			// std::cout << ">" << seedStartPos << "." << newSeed.seedLength << ":" << maxKmerFreq << "\t" << dynamicKmerSize << "\t" << dynamicKmerThreshold <<  "\n" << newSeed.seedStr << "\n";
			
			// jump to the index after new seed
			i = seedEndPos;
			dynamicKmerSize = maxKmerSize;
			dynamicKmerThreshold = kmerThreshold;
		}
		
		// reduce kmer size if no seeds can be found within a maximum interval
		int prevSeedEndPos = seedEndPosVec.empty() ? 0 : seedEndPosVec.back()+1;
		int distToPrevSeed = i + 1 - prevSeedEndPos;
		if(distToPrevSeed >= (int)maxSeedInterval.at(dynamicKmerSize))
		{
			// reducing the kmer size
			if(dynamicKmerSize > minKmerSize)
			{
				i = prevSeedEndPos;
				dynamicKmerSize -= 2;
				dynamicKmerThreshold +=2;
			}
			else
			{
				// collect seeds from PacBio reads if distToPrevSeed is too long due to sequencing gaps
				prevSeedEndPos = seedVec.empty() ? 0 : seedVec.back().seedEndPos+1;
				distToPrevSeed = i + 1 - prevSeedEndPos;
				
				/* In large sequencing gaps (>7kb), no seeds can be found in Illumina index
					If no seed is found within PBSearchDepth, 
					seeds from PB index will be serached instead*/
				if(distToPrevSeed >= m_params.PBSearchDepth)
				{	
					seedingByPacBio(readSeq, seedVec, seedEndPosVec, prevSeedEndPos);

				}
				// no seeds can be found even using minKmerSize, 
				// so we reset start position of searching to give up the maximum interval of minKmerSize.
				seedEndPosVec.push_back(i);
				dynamicKmerSize = maxKmerSize;
				dynamicKmerThreshold = kmerThreshold;

			}
		}
	}
	
	return seedVec;
}

// PacBio Hybrid Correction Seeding By Dynamic Kmer v2, v20160511 by Ya.
// find seeds by dynamic kmers, which is suitable for the PacBio hybrid error correction,
// where repeat regions require large kmers and error-prone regions require small kmers.
std::vector<SeedFeature> PacBioHybridCorrectionProcess::seedingByDynamicKmer_v2(const string& readSeq)
{
	std::vector<std::vector<size_t> > PBKmerFreqsVec;
	std::vector<SeedFeature> seedVec;
	std::vector<int> seedEndPosVec;
	size_t maxKmerSize = m_params.kmerLength;
	size_t minKmerSize = m_params.minKmerLength;
	size_t kmerThreshold = m_params.seedKmerThreshold;
	
	// prevention of short reads
	if(readSeq.length() < maxKmerSize)
		return seedVec;
	
	// computing maximum interval of variously consecutive matches length (dk figure)
	std::vector<int> maxSeedInterval;
	for(size_t i = 0 ; i <= maxKmerSize ; i++)
		maxSeedInterval.push_back(2*3.8649*pow(2.7183,0.1239*i));
	
	int numValidPos = calculateKmerFreqsEachPBPos(readSeq, PBKmerFreqsVec);
	// only 24% of pacbio reads have more than two seeds
	// if number of valid pos (high kmer freqs in the pacbio read) is smaller than two
	// we end the seeding function as early as possible
	if(numValidPos < 2)
		return seedVec;
	
	// set dynamic kmer as largest kmer initially, 
	// which will reduce size when no seeds can be found within a maximum interval.
	size_t dynamicKmerSize = maxKmerSize;
	//size_t dynamicKmerThreshold = kmerThreshold;

	// search for solid kmers and group consecutive solids kmers into one seed
	// reduce kmer size if no seeds can be found within maxSeedInterval
	for(size_t i = 0 ; i+dynamicKmerSize <= readSeq.length() ; i++)
	{
		size_t PBKmerSize = minKmerSize+PBKmerFreqsVec.at(i).size()-1;
		if(PBKmerSize >= dynamicKmerSize)
		{
			size_t seedStartPos = i;
			size_t maxKmerFreq = PBKmerFreqsVec.at(i).back();

			// group consecutive solid kmers into one seed if possible
			for(i++ ; i+dynamicKmerSize <= readSeq.length() ; i++)
			{
				PBKmerSize = minKmerSize+PBKmerFreqsVec.at(i).size()-1;
				if(PBKmerSize >= dynamicKmerSize)
					maxKmerFreq = std::max(maxKmerFreq, PBKmerFreqsVec.at(i).back());
				else
				{
					PBKmerSize = minKmerSize+PBKmerFreqsVec.at(i-1).size()-1;
					break;
				}
			}
			
			// small-sized seed has 50% chance of errors in C elegans, 
			// skip only if there is another seed nearby 30bp
			if( (i-seedStartPos) < 2 && !seedVec.empty() && i-seedVec.back().seedEndPos <= 30)
				continue;

			size_t seedEndPos = i+PBKmerSize-2;

			// repeat seeds are less accurate at boundary and should be trimmed 
			if(maxKmerFreq > m_params.coverage*4)
				trimRepeatSeed(readSeq, m_params.coverage, seedStartPos, seedEndPos);

			// super repeat seeds with frequency > 2000 are troublesome, often lead to -3 but no good solution so far, mark first
			bool isSuperRepeat = maxKmerFreq > m_params.coverage*15?true:false;
			SeedFeature newSeed(seedStartPos, readSeq.substr(seedStartPos, seedEndPos-seedStartPos+1), isSuperRepeat, PBKmerSize, kmerThreshold*4);

			// skip low-complexity sequencing errors of PacBio
			// bool isShortAndHighFreq = i-seedStartPos <= 2 && maxKmerFreq > 80;
			if(!isLowComplexity(newSeed.seedStr, 0.9))
				seedVec.push_back(newSeed);
			
			seedEndPosVec.push_back(seedEndPos);
			
			// debug
			// std::cout << ">" << seedStartPos << "." << newSeed.seedLength << ":" << maxKmerFreq << "\t" << dynamicKmerSize << "\t" << dynamicKmerThreshold <<  "\n" << newSeed.seedStr << "\n";
			
			// jump to the index after new seed
			i = seedEndPos;
			dynamicKmerSize = maxKmerSize;
			//dynamicKmerThreshold = kmerThreshold;
		}
		
		// reduce kmer size if no seeds can be found within a maximum interval
		int prevSeedEndPos = seedEndPosVec.empty() ? 0 : seedEndPosVec.back()+1;
		int distToPrevSeed = i + 1 - prevSeedEndPos;
		if(distToPrevSeed >= (int)maxSeedInterval.at(dynamicKmerSize))
		{
			// reducing the kmer size
			if(dynamicKmerSize > minKmerSize)
			{
				i = prevSeedEndPos;
				dynamicKmerSize -= 2;
				//dynamicKmerThreshold +=2;
			}
			else
			{
				// collect seeds from PacBio reads if distToPrevSeed is too long due to sequencing gaps
				prevSeedEndPos = seedVec.empty() ? 0 : seedVec.back().seedEndPos+1;
				distToPrevSeed = i + 1 - prevSeedEndPos;
				
				// In large sequencing gaps (>7kb), no seeds can be found in Illumina index
				// If no seed is found within PBSearchDepth, 
				// seeds from PB index will be serached instead
				if(distToPrevSeed >= m_params.PBSearchDepth)
					seedingByPacBio(readSeq, seedVec, seedEndPosVec, prevSeedEndPos);
				
				// no seeds can be found even using minKmerSize, 
				// so we reset start position of searching to give up the maximum interval of minKmerSize.
				seedEndPosVec.push_back(i);
				dynamicKmerSize = maxKmerSize;
				//dynamicKmerThreshold = kmerThreshold;
			}
		}
	}
	
	return seedVec;
}

int PacBioHybridCorrectionProcess::calculateKmerFreqsEachPBPos(const std::string& readSeq, std::vector<std::vector<size_t> >& PBKmerFreqsVec)
{
	size_t maxKmerSize = m_params.kmerLength;
	size_t minKmerSize = m_params.minKmerLength;
	size_t kmerThreshold = m_params.seedKmerThreshold;
	assert(kmerThreshold - 2*(maxKmerSize-minKmerSize) > 0);
	
	PBKmerFreqsVec.clear();
	PBKmerFreqsVec.resize(readSeq.length());
	
	int numValidPos = 0;
	for(size_t pos = 0 ; pos+maxKmerSize<=readSeq.length() ; pos++)
	{
		string kmer = readSeq.substr(pos, minKmerSize);
		BWTInterval fwdInterval = BWTAlgorithms::findInterval(m_params.indices.pRBWT, reverse(kmer));
		BWTInterval rvcInterval = BWTAlgorithms::findInterval(m_params.indices.pBWT, reverseComplement(kmer));
		size_t kmerFreqs = fwdInterval.size() + rvcInterval.size();
		size_t dynamicKmerThreshold = kmerThreshold;
		
		if(kmerFreqs < dynamicKmerThreshold)
			continue;
	
		PBKmerFreqsVec.at(pos).push_back(kmerFreqs);
		numValidPos++;
		for(size_t extendPos=pos+minKmerSize; extendPos<pos+maxKmerSize && extendPos<readSeq.length() ; extendPos++)
		{
			char b = readSeq.at(extendPos);
			char rcb;
			switch(b)
			{
				case 'A': rcb='T'; break;
				case 'T': rcb='A'; break;
				case 'C': rcb='G'; break;
				case 'G': rcb='C'; break;
			}
			BWTAlgorithms::updateInterval(fwdInterval, b, m_params.indices.pRBWT);
			BWTAlgorithms::updateInterval(rvcInterval, rcb, m_params.indices.pBWT);
			kmerFreqs = fwdInterval.size() + rvcInterval.size();
			dynamicKmerThreshold -= 2;
			if(kmerFreqs >= dynamicKmerThreshold)
				PBKmerFreqsVec.at(pos).push_back(kmerFreqs);
			else
				break;
		}
	}
	
	return numValidPos;
}

// PacBio Hybrid Correction Seeding By Dynamic Kmer v3, v20160517 by Ya.
// find seeds by dynamic kmers, which is suitable for the PacBio hybrid error correction,
// where repeat regions require large kmers and error-prone regions require small kmers.
std::vector<SeedFeature> PacBioHybridCorrectionProcess::seedingByDynamicKmer_v3(const string& readSeq)
{
	std::vector<SeedFeature> seedVec;
	std::vector<int> seedEndPosVec;
	size_t maxKmerSize=m_params.kmerLength;
	size_t minKmerSize=m_params.minKmerLength;
	
	// YTH
	// int prevRepeatEnd = -10000;
	// size_t prevRepeatFreq = 0;
	
	// prevention of short reads
	if(readSeq.length()<=maxKmerSize)
		return seedVec;
	
	// computing kmer threshold of various kmer size
	// kmer size	kmer threshold*(short reads coverage/100)+3
	// 21	28*(short reads coverage/100)+3
	// 31	21*(short reads coverage/100)+3
	// 41	15*(short reads coverage/100)+3
	// 51	10*(short reads coverage/100)+3
	// 61	6*(short reads coverage/100)+3
	// 71	3*(short reads coverage/100)+3
	// 81	1*(short reads coverage/100)+3
	// 91	0*(short reads coverage/100)+3
	// 101	0*(short reads coverage/100)+3
	// 201	0*(short reads coverage/100)+3
	// regression: kmer threshold=(0.005*(kmer size)^2-0.96*kmer size+45.955)*(short reads coverage/100)+3
	std::vector<float> kmerThreshold;
	kmerThreshold.resize(201+1,3);
	for(size_t kmerSize=0 ; kmerSize<=91 ; kmerSize++)
	{
		float kmerThresholdValue=(0.005*pow(kmerSize,2)-0.96*kmerSize+45.955)*((float)m_params.coverage/100);
		kmerThreshold.at(kmerSize)+=kmerThresholdValue;
	}
	
	// search for solid kmers as seeds
	for(size_t pos=0 ; pos+minKmerSize<readSeq.length() ; pos++)
	{
		string kmer=readSeq.substr(pos,minKmerSize);
		BWTInterval fwdInterval=BWTAlgorithms::findInterval(m_params.indices.pRBWT, reverse(kmer));
		BWTInterval rvcInterval=BWTAlgorithms::findInterval(m_params.indices.pBWT, reverseComplement(kmer));
		size_t kmerFreqs=fwdInterval.size()+rvcInterval.size();
		size_t dynamicKmerSize=minKmerSize;
		size_t dynamicKmerThreshold=kmerThreshold.at(minKmerSize);
		
		if(kmerFreqs<dynamicKmerThreshold)
		{
			// In large sequencing gaps (>7kb), no seeds can be found in Illumina index
			// If no seed is found within PBSearchDepth, 
			// seeds from PB index will be serached instead
			int prevSeedEndPos=seedEndPosVec.empty()?0:seedEndPosVec.back()+1;
			int distToPrevSeed=pos+1-prevSeedEndPos;

			// skip repeat regions
			// bool isWithinRepeat = pos - prevRepeatEnd < m_params.PBSearchDepth;
			if(/*!isWithinRepeat &&*/ distToPrevSeed>=m_params.PBSearchDepth)
			{
				// if(seedingByPacBio(readSeq,seedVec,seedEndPosVec,prevSeedEndPos))
				if(seedingByPacBio_v2(readSeq,seedVec,seedEndPosVec,prevSeedEndPos))
					;//cout << seedVec.back().seedStr << " " << seedVec.back().seedStartPos << " " << BWTAlgorithms::countSequenceOccurrencesSingleStrand(seedVec.back().seedStr, m_params.PBindices)+BWTAlgorithms::countSequenceOccurrencesSingleStrand(reverseComplement(seedVec.back().seedStr), m_params.PBindices) << endl;
				else
					seedEndPosVec.push_back(pos);
				//cout << seedEndPosVec.size() << endl;
				pos=seedEndPosVec.back();
			}
			continue;
		}
		
		size_t seedStartPos=pos;
		size_t maxKmerFreq=kmerFreqs;
		
		// search for longest solid kmer as one seed if possible
		for(pos=pos+minKmerSize ; pos<readSeq.length() ; pos++)
		{
			char b=readSeq.at(pos);
			char rcb;
			switch(b)
			{
				case 'A': rcb='T'; break;
				case 'T': rcb='A'; break;
				case 'C': rcb='G'; break;
				case 'G': rcb='C'; break;
			}
			BWTAlgorithms::updateInterval(fwdInterval,b,m_params.indices.pRBWT);
			BWTAlgorithms::updateInterval(rvcInterval,rcb,m_params.indices.pBWT);
			kmerFreqs = fwdInterval.isValid()?fwdInterval.size():0 + rvcInterval.isValid()?rvcInterval.size():0;
			
			dynamicKmerSize++;			
			assert(dynamicKmerSize<=kmerThreshold.size());
			dynamicKmerThreshold=kmerThreshold.at(dynamicKmerSize);
			
			if(kmerFreqs>=dynamicKmerThreshold)
				maxKmerFreq=kmerFreqs;
			else
			{
				dynamicKmerSize--;
				dynamicKmerThreshold=kmerThreshold.at(dynamicKmerSize);
				break;
			}
		}
		
		// small-sized seed has 50% chance of errors in C elegans, 
		// skip only if there is another seed nearby 30bp
		//if( (pos-1-minKmerSize-seedStartPos) < 2 && !seedVec.empty() && pos-1-minKmerSize-seedVec.back().seedEndPos <= 30)
		//if(!seedVec.empty() && pos-1-minKmerSize-seedVec.back().seedEndPos <= 30)
		//	continue;

		size_t seedEndPos=pos-1;

		// repeat seeds are less accurate at boundary and should be trimmed 
		if(maxKmerFreq >= m_params.coverage*4)
		{
			trimRepeatSeed(readSeq,m_params.coverage,seedStartPos,seedEndPos);

		}
		
		// bool isTooCloseToPrevSeed = seedEndPos - prevRepeatEnd <= prevRepeatFreq*0.02+100;
		// if(isTooCloseToPrevSeed)
			// continue;
				
		// super repeat seeds with frequency > 2000 are troublesome, often lead to -3 but no good solution so far, mark first
		bool isSuperRepeat=maxKmerFreq>m_params.coverage*15?true:false;
		SeedFeature newSeed(seedStartPos, readSeq.substr(seedStartPos, seedEndPos-seedStartPos+1), isSuperRepeat, dynamicKmerSize, m_params.PBcoverage/2);
		newSeed.estimateBestKmerSize(m_params.PBindices.pBWT);

		// Some high-GC error seeds have unexpected large freq
		size_t PBkmerFreq = BWTAlgorithms::countSequenceOccurrences(newSeed.seedStr, m_params.PBindices);
		if(PBkmerFreq > m_params.PBcoverage/4)
			continue;

		// skip low-complexity sequencing errors of PacBio
		if(!isLowComplexity(newSeed.seedStr,0.9))
		{
			// std::cout << ">" << seedStartPos << ":" << maxKmerFreq << ":" << readSeq.substr(seedStartPos,seedEndPos-seedStartPos+1) << "\t" 
			// << prevRepeatEnd << ":" << PBkmerFreq << "\n";

			seedVec.push_back(newSeed);
		}
			
		seedEndPosVec.push_back(seedEndPos);
		// jump to the index after new seed
		pos=seedEndPos;
	}
	
	return seedVec;
}

bool PacBioHybridCorrectionProcess::seedingByPacBio(const string& readSeq, std::vector<SeedFeature>& seedVec, 	std::vector<int>& seedEndPosVec, size_t prevEndPos)
{
	size_t dynamicKmerSize = 17;
	size_t dynamicKmerThreshold = 8;
	
	// search for solid kmers and group consecutive solids kmers into one seed
	// reduce kmer size if no seeds can be found within maxSeedInterval
	for(int i = prevEndPos ; i+dynamicKmerSize <= readSeq.length() ; i++)
	{
		string kmer = readSeq.substr(i, dynamicKmerSize);
		size_t fwdKmerFreqs = BWTAlgorithms::countSequenceOccurrencesSingleStrand(kmer, m_params.PBindices);
		size_t rvcKmerFreqs = BWTAlgorithms::countSequenceOccurrencesSingleStrand(reverseComplement(kmer), m_params.PBindices);
		size_t kmerFreqs = fwdKmerFreqs+rvcKmerFreqs;
		
		// std::cout << i << "::\t" << kmer << "\t" << kmerFreqs << "\n";
		// finding the seed
		if(kmerFreqs >= dynamicKmerThreshold)
		{
			size_t seedStartPos = i;
			size_t maxKmerFreq = kmerFreqs;
			
			// group consecutive solid kmers into one seed if possible
			for(i++ ; i+dynamicKmerSize <= readSeq.length() ; i++)
			{
				kmer = readSeq.substr(i, dynamicKmerSize);
				fwdKmerFreqs = BWTAlgorithms::countSequenceOccurrencesSingleStrand(kmer, m_params.PBindices);
				rvcKmerFreqs = BWTAlgorithms::countSequenceOccurrencesSingleStrand(reverseComplement(kmer), m_params.PBindices);
				kmerFreqs = fwdKmerFreqs+rvcKmerFreqs;
				maxKmerFreq = std::max(maxKmerFreq,kmerFreqs);

				  // std::cout << i << "::\t" << kmer << "\t" << kmerFreqs << "\n";

				if(kmerFreqs < dynamicKmerThreshold)
					break;
			}
			
			//skip repeat seeds
			if(maxKmerFreq > 200) continue;

			// require sufficient seed length in repeats
			if(i-seedStartPos <= 4){
				// std::cout << readSeq.substr(seedStartPos, i+dynamicKmerSize-2-seedStartPos+1) << "\n";
				continue;
			}

			size_t seedEndPos = i+dynamicKmerSize-2;

			// repeat seeds are less accurate at boundary and should be trimmed 
			// if(maxKmerFreq > m_params.coverage*4)
				// trimRepeatSeed(readSeq, m_params.coverage, seedStartPos, seedEndPos);

			// super repeat seeds with frequency > 2000 are troublesome, often lead to -3 but no good solution so far, mark first
			bool isSuperRepeat = maxKmerFreq > m_params.coverage*4?true:false;
			SeedFeature newSeed(seedStartPos, readSeq.substr(seedStartPos, seedEndPos-seedStartPos+1), isSuperRepeat, dynamicKmerSize, m_params.coverage*4);

			// skip low-complexity sequencing errors of PacBio
			// bool isShortAndHighFreq = i-seedStartPos <= 2 && maxKmerFreq > 80;
			if(!isLowComplexity(newSeed.seedStr, 0.9))
			{
				// debug
				 // std::cout << ">>" << seedStartPos << "." << newSeed.seedLength << ":" << maxKmerFreq << "\t" << dynamicKmerSize << "\t" << dynamicKmerThreshold <<  "\n" << newSeed.seedStr << "\n";

				newSeed.isPBSeed = true;
				seedVec.push_back(newSeed);
				seedEndPosVec.push_back(seedEndPos);
				return true;
			}
		}
		
		// reduce kmer size if no seeds can be found within a maximum interval
		// int prevSeedEndPos = seedEndPosVec.empty() ? 0 : seedEndPosVec.back()+1;
		int prevSeedEndPos = seedVec.empty() ? 0 : seedVec.back().seedEndPos+1;
		int distToPrevSeed = i + 1 - prevSeedEndPos;
		if(distToPrevSeed >= 1000)
		{
			return false;
		}
	}
	
	return false;
}

bool PacBioHybridCorrectionProcess::seedingByPacBio_v2(const string& readSeq, std::vector<SeedFeature>& seedVec, 	std::vector<int>& seedEndPosVec, size_t prevEndPos)
{	
	// computing kmer threshold of various kmer size
	// kmer size	kmer threshold*(PB reads coverage/60)+5
	// 17	8*(PB reads coverage/60)+5
	// 27	7*(PB reads coverage/60)+5
	// 37	6*(PB reads coverage/60)+5
	// 47	5*(PB reads coverage/60)+5
	// 57	4*(PB reads coverage/60)+5
	// 67	3*(PB reads coverage/60)+5
	// 77	2*(PB reads coverage/60)+5
	// 87	1*(PB reads coverage/60)+5
	// 97	0*(PB reads coverage/60)+5
	// regression: kmer threshold=(-0.1*(kmer size)+9.7)*(PB reads coverage/60)+5
	std::vector<float> kmerThreshold;
	kmerThreshold.resize(97+1,5);
	for(size_t kmerSize=0 ; kmerSize<=97 ; kmerSize++)
	{
		float kmerThresholdValue=(-0.1*kmerSize+9.7)*((float)m_params.PBcoverage/60);
		kmerThreshold.at(kmerSize)+=kmerThresholdValue;
	}
	
	// search for solid kmers as seeds
	for(size_t pos=prevEndPos ; pos+m_params.PBKmerLength<readSeq.length() && pos-prevEndPos <= m_params.PBSearchDepth ; pos++)
	{
		size_t dynamicKmerSize = m_params.PBKmerLength;

		string kmer=readSeq.substr(pos,dynamicKmerSize);
		BWTInterval fwdInterval=BWTAlgorithms::findInterval(m_params.PBindices.pRBWT, reverse(kmer));
		BWTInterval rvcInterval=BWTAlgorithms::findInterval(m_params.PBindices.pBWT, reverseComplement(kmer));
		size_t kmerFreqs=fwdInterval.size()+rvcInterval.size();
		size_t dynamicKmerThreshold=kmerThreshold.at(dynamicKmerSize);
		
		if(kmerFreqs<dynamicKmerThreshold)
			continue;

		size_t seedStartPos=pos;
		size_t maxKmerFreq=kmerFreqs;
		
		// search for longest solid kmer as one seed if possible
		for(pos=pos+dynamicKmerSize ; pos+dynamicKmerSize<readSeq.length() ; pos++)
		{
			char b=readSeq.at(pos);
			char rcb;
			switch(b)
			{
				case 'A': rcb='T'; break;
				case 'T': rcb='A'; break;
				case 'C': rcb='G'; break;
				case 'G': rcb='C'; break;
			}
			BWTAlgorithms::updateInterval(fwdInterval,b,m_params.PBindices.pRBWT);
			BWTAlgorithms::updateInterval(rvcInterval,rcb,m_params.PBindices.pBWT);
			kmerFreqs=fwdInterval.isValid()?fwdInterval.size():0 + rvcInterval.isValid()?rvcInterval.size():0;
			
			dynamicKmerSize++;
			// assert(dynamicKmerSize<=kmerThreshold.size());
			if(dynamicKmerSize>=kmerThreshold.size()) break;
			
			dynamicKmerThreshold=kmerThreshold.at(dynamicKmerSize);
			
			if(kmerFreqs>=dynamicKmerThreshold)
				maxKmerFreq=kmerFreqs;
			else
			{
				dynamicKmerSize--;
				dynamicKmerThreshold=kmerThreshold.at(dynamicKmerSize);
				break;
			}
		}

		//skip repeat seeds
		if(maxKmerFreq >= m_params.PBcoverage*2)
			continue;

		// require sufficient seed length in repeats
		if(maxKmerFreq >= m_params.PBcoverage && dynamicKmerSize-m_params.PBKmerLength <= 4)
			// std::cout << readSeq.substr(seedStartPos, i+dynamicKmerSize-2-seedStartPos+1) << "\n";
			continue;

		size_t seedEndPos = pos-1;

		// repeat seeds are less accurate at boundary and should be trimmed 
		// if(maxKmerFreq > m_params.coverage*4)
			// trimRepeatSeed(readSeq, m_params.coverage, seedStartPos, seedEndPos);

		// super repeat seeds with frequency > 2000 are troublesome, often lead to -3 but no good solution so far, mark first
		bool isSuperRepeat = maxKmerFreq >= m_params.PBcoverage?true:false;
		SeedFeature newSeed(seedStartPos, readSeq.substr(seedStartPos, seedEndPos-seedStartPos+1), isSuperRepeat, dynamicKmerSize, m_params.PBcoverage/2);
		newSeed.estimateBestKmerSize(m_params.PBindices.pBWT);

		// skip low-complexity sequencing errors of PacBio
		// bool isShortAndHighFreq = i-seedStartPos <= 2 && maxKmerFreq > 80;
		if(!isLowComplexity(newSeed.seedStr,0.8))
		{
			// debug
			// std::cout << ">>" << seedStartPos << "." << newSeed.seedLength << ":" << maxKmerFreq << "\t" << dynamicKmerSize << "\t" << dynamicKmerThreshold <<  "\n" << newSeed.seedStr << "\n";
			newSeed.isPBSeed = true;
			seedVec.push_back(newSeed);
			seedEndPosVec.push_back(seedEndPos);
			return true;
		}
	}
	
	return false;
}

// return is low complexity of seq, 20160525 by Ya.
bool PacBioHybridCorrectionProcess::isLowComplexity(std::string& seq ,const float& ratioThreshold)
{
	size_t seqLen=seq.length();
	size_t countG=0;
	size_t countC=0;
	size_t countT=0;
	size_t countA=0;

	for (size_t i=0;i<seqLen;i++)
	{
		switch(seq[i])
		{
			case 'A': countA ++ ;break;
			case 'T': countT ++ ;break;
			case 'C': countC ++ ;break;
			case 'G': countG ++ ;break;
			default:  assert(false);
		}
	}

	// GCratio = (float)(countG+countC)/seqLen ;

	if( ((float)countA/seqLen>=ratioThreshold) ||
		((float)countT/seqLen>=ratioThreshold) ||
		((float)countC/seqLen>=ratioThreshold) ||
		((float)countG/seqLen>=ratioThreshold) ||
		(countA == 0) ||
		(countT == 0) ||
		(countC == 0) ||
		(countG == 0) )
		return true;

	return false;
}

int PacBioHybridCorrectionProcess::extendBetweenSeeds(SeedFeature& source, SeedFeature& target, string& strBetweenSrcTarget, int dis_between_src_target, FMWalkResult* FMWResult, int debugTargetSeed)
{	
	int FMWalkReturnType = -2;
	int PrevFMWalkReturnType = 0;	
	
	// compute minOverlap from min of seedLength and maxOverlap
	int minOverlap = std::min(source.seedLength, target.seedLength);
		minOverlap = std::min(minOverlap, m_params.maxOverlap);

	int initMinOverlap = minOverlap;
	
	FMWalkParameters FMWParams;
	FMWParams.indices = m_params.indices;
	FMWParams.maxOverlap = m_params.maxOverlap;
	FMWParams.SAThreshold = m_params.FMWKmerThreshold;
	FMWParams.disBetweenSrcTarget = dis_between_src_target;
	FMWParams.maxLeaves = m_params.maxLeaves;
	FMWParams.coverage = m_params.coverage;
	
	// debug
	// if(debugTargetSeed == 2){
		// FMWParams.debugMode = true;
		// std::cout << "\nfmwalk id: " << debugTargetSeed << ", dis_between_src_target length: " << dis_between_src_target << "\t"<< minOverlap << ".----\n";
	// }
	
	bool isSequencingGap = false;
	bool isSeedfromPB = source.isPBSeed || target.isPBSeed || source.isNextRepeat;
	
	// Correction by FM-index extension from source to target with iterative minOverlap reduction
	while( (FMWalkReturnType == -1 || FMWalkReturnType == -2) && 
			(minOverlap >= m_params.minKmerLength) && !isSeedfromPB)
	{
		FMWParams.sourceSeed = source.seedStr;
		FMWParams.targetSeed = target.seedStr;
		FMWParams.strBetweenSrcTarget = strBetweenSrcTarget;
		FMWParams.minOverlap = minOverlap;
		SAIntervalPBHybridCTree SAITree(FMWParams);
		FMWalkReturnType = SAITree.mergeTwoSeeds(*FMWResult);
		
		// if(target.seedLength == 27)
			// std::cout << FMWalkReturnType << ":\t" << minOverlap <<"\n";

		// Correction by FM-index extension from target to source
		if(FMWalkReturnType > 0)
		{
			FMWalkResult FMWResult2;
			FMWParams.sourceSeed = reverseComplement(target.seedStr);
			FMWParams.targetSeed = reverseComplement(source.seedStr);
			FMWParams.strBetweenSrcTarget = reverseComplement(strBetweenSrcTarget);
			SAIntervalPBHybridCTree SAITree2(FMWParams);
			FMWalkReturnType = SAITree2.mergeTwoSeeds(FMWResult2);

			//if(dis_between_src_target==236)
			//	std::cout << FMWalkReturnType << ":\t" << (*mergedseq).length() << "\t" << mergedseq2.length() << "\n";

			if((*FMWResult).mergedSeq.length() == FMWResult2.mergedSeq.length())
			{
				// select the extended seq with higher alignment score
				if((*FMWResult).alnScore < FMWResult2.alnScore)
				{
					FMWResult2.mergedSeq = reverseComplement(FMWResult2.mergedSeq);
					(*FMWResult) = FMWResult2;
				}
				break;
			}
			// Give up the extension if 2nd extension does not match the 1st one.
			if(FMWalkReturnType > 0)
				FMWalkReturnType = -4;
		}
		
		/*** Extension Failed 1) try smaller minOverlap 2) try correction by MSA for sequencing gap ***/
		// Sequencing gaps in short reads often lead to -1 or -2
		// But over-reduced minOverlap size also lead to -2 due to error seeds
		if( (FMWalkReturnType == -2 && minOverlap >= m_params.kmerLength) || 
			(FMWalkReturnType == -1 && minOverlap == initMinOverlap) )
		// if(minOverlap==initMinOverlap&&(FMWalkReturnType==-2||FMWalkReturnType==-1))
			isSequencingGap = true;
		
		// The FMWalkReturnType may switch from -1 to -3 due to over-reduced minOverlap
		// We have to skip further reduction and reset the FMWalkReturnType
		if(FMWalkReturnType == -3 && PrevFMWalkReturnType == -1)
		{
			// the current -3 state is incorrect, reset to prev status
			FMWalkReturnType = PrevFMWalkReturnType;
			break;
		}
		
		PrevFMWalkReturnType = FMWalkReturnType;
		minOverlap--;	
		
		// FMWalkReturnType -2 leads to misassembly in C elegans
		// FMWalkReturnType < 0 generates bette r assembly in E coli
		if(source.isRepeat && minOverlap < m_params.kmerLength-1) break;
	}

	// std::cout << FMWalkReturnType << "\t" << target.seedStr << "\t" << source.isNextRepeat << target.isRepeat << source.isPBSeed << target.isPBSeed << "\n";
	
	// try correction by MSA using FM-index of low quality long reads for sequencing gaps
	if( (FMWalkReturnType == -2 || FMWalkReturnType == -1) && !source.isRepeat && !target.isRepeat 
		&& (isSequencingGap || isSeedfromPB) )
	{
		// mimic the overlap correction process
		// const std::string query = source.seedStr.substr(source.seedLength - m_params.PBKmerLength)+strBetweenSrcTarget.substr(10, dis_between_src_target)+target.seedStr;

		// switch to best kmer in PB index
		const std::string query = source.seedStr.substr(source.seedLength - source.endBestKmerSize)+strBetweenSrcTarget.substr(10, dis_between_src_target)+target.seedStr;
		size_t sourceKmerLength = source.endBestKmerSize;
		size_t targetKmerLength = target.endBestKmerSize;

		// std::cout << minOverlap << "\t" << query.length() << "\n" << query << "\n";
		
		// Self correction by aligning all reads having seeds with query
		MultipleAlignment maquery = LongReadOverlap::buildMultipleAlignment(query,
													sourceKmerLength, //m_params.PBKmerLength, 
													targetKmerLength, //m_params.PBKmerLength,
													query.length()/10, 
													0.73,	// alignment identity < 0.7 are often false positive repeats
													m_params.PBcoverage,
													m_params.PBindices);

		// skip insufficient number of overlapping reads for correction
		if(maquery.getNumRows() <= 3)
			return FMWalkReturnType;
		
		// maquery.print(120);
		std::string consensus = maquery.calculateBaseConsensus(100000, -1);
		FMWResult->mergedSeq = source.seedStr + consensus.substr(m_params.PBKmerLength);
		// std::cout << ">" << consensus.length() <<"\n" << consensus << "\n";
		
		return 1;
		
		/* Self-correction by FM-index extenion is suitable for high-coverage PB sequencing
		
		FMWParams.lowCoverageHighErrorMode = true;
		FMWParams.minOverlap = m_params.minKmerLength;
		FMWParams.maxOverlap = m_params.kmerLength;
		FMWParams.sourceSeed = source.seedStr;
		FMWParams.targetSeed = target.seedStr;
		SAIntervalPBHybridCTree SAITree(FMWParams);

		const double maxRatio = 1.1;
		const double minRatio = 0.9;
		const int minOffSet = 30;	//PB159615_16774.fa contains large indels > 30bp
		const size_t extendKmerSize = 15;
		const size_t srcKmerSize = 17;
		
		SAIPBSelfCorrectTree SAITree(m_params.PBindices.pBWT, m_params.PBindices.pRBWT, strBetweenSrcTarget, 2);
			
		// this occurs when one end is repeat requiring larger extendKmerSize while the other requires small extendKmerSize
		// the source str should be the longest one
		// const int srcMaxLength = maxRatio*(dis_between_src_target+minOffSet) + source.seedLength + extendKmerSize;
		// size_t sourceFreq = SAITree.addHashBySingleSeed(source.seedStr, source.endBestKmerSize, extendKmerSize, srcMaxLength, m_params.isFirst);
		
		std::string srcStr = source.seedStr.substr(source.seedStr.length()-srcKmerSize);
		size_t srcMaxLength = maxRatio*(dis_between_src_target+minOffSet) + srcStr.length() + extendKmerSize;
		// size_t sourceFreq = SAITree.addHashBySingleSeed(srcStr, srcKmerSize, extendKmerSize, srcMaxLength, true);

		// srcStr = source.seedStr.substr(source.seedStr.length()-srcKmerSize*3.5, srcKmerSize);
		// SAITree.addHashBySingleSeed(srcStr, srcKmerSize, extendKmerSize, srcMaxLength, true);		
		srcStr = source.seedStr.substr(source.seedStr.length()-srcKmerSize*2, srcKmerSize);
		SAITree.addHashBySingleSeed(srcStr, srcKmerSize, extendKmerSize, srcMaxLength, true);
		srcStr = source.seedStr.substr(source.seedStr.length()-srcKmerSize*3, srcKmerSize);
		SAITree.addHashBySingleSeed(srcStr, srcKmerSize, extendKmerSize, srcMaxLength, true);
		srcStr = source.seedStr.substr(source.seedStr.length()-srcKmerSize*1.5, srcKmerSize);
		SAITree.addHashBySingleSeed(srcStr, srcKmerSize, extendKmerSize, srcMaxLength, true);		
		
		srcStr = source.seedStr.substr(source.seedStr.length()-srcKmerSize);
		
		// Collect local kmer frequency from target upto targetMaxLength
		std::string rvcTargetStr = reverseComplement(target.seedStr);
		const int targetMaxLength = maxRatio*(dis_between_src_target+minOffSet) + rvcTargetStr.length() + srcKmerSize;
		size_t expectedLength = dis_between_src_target + rvcTargetStr.length();
		assert(rvcTargetStr.length()>=extendKmerSize);
		size_t targetFreq = SAITree.addHashBySingleSeed(rvcTargetStr, srcKmerSize, extendKmerSize, targetMaxLength, true, expectedLength);

		int srcMinLength = minRatio*(dis_between_src_target-minOffSet) + srcStr.length() + extendKmerSize;
		if(srcMinLength < 0) srcMinLength = 0;
		expectedLength = srcStr.length() + dis_between_src_target + target.seedLength;
		
		std::string pbseq;
		FMWalkReturnType = SAITree.mergeTwoSeedsUsingHash(srcStr, target.seedStr, pbseq, extendKmerSize, m_params.maxLeaves,
														  srcMinLength, srcMaxLength, expectedLength);	

		// debug
		// std::cout << dis_between_src_target << "\t" << FMWalkReturnType << "\t"<< source.seedStartPos << "\t" << target.seedStr << "\n";
		if(!pbseq.empty())
			*mergedseq = source.seedStr+pbseq.substr(srcKmerSize);
		*/
	}
	
	// if(FMWalkReturnType==-2)
	// {
		// if(debugTargetSeed == 3){
			// FMWParams.debugMode = true;
			// std::cout << "\nfmwalk id: " << debugTargetSeed << ", dis_between_src_target length: " << dis_between_src_target << "\t"<< minOverlap << ".----\n";
		// }
		// FMWParams.lowCoverageHighErrorMode = true;
		// FMWParams.minOverlap = m_params.minKmerLength;
		// FMWParams.maxOverlap = m_params.kmerLength;
		// FMWParams.sourceSeed = source.seedStr;
		// FMWParams.targetSeed = target.seedStr;
		// SAIntervalPBHybridCTree SAITree(FMWParams);
		// FMWalkReturnType = SAITree.mergeTwoSeeds(*FMWResult);
		// if(FMWalkReturnType < 0)
		// {
			// FMWParams.sourceSeed = reverseComplement(target.seedStr);
			// FMWParams.targetSeed = reverseComplement(source.seedStr);
			// SAIntervalPBHybridCTree SAITree(FMWParams);
			// FMWalkReturnType = SAITree.mergeTwoSeeds(*mergedseq);
			// (*mergedseq) = reverseComplement((*mergedseq));
		// }
	// }
	
	return FMWalkReturnType;
}


// boundary of repeat seeds are less reliable
// trim the base with lower frequency in comparison with others
void PacBioHybridCorrectionProcess::trimRepeatSeed(const string& readSeq, size_t coverage, size_t& seedStartPos, size_t& seedEndPos)
{
	// initially set to max unisnged int value
	size_t newSeedStartPos = (size_t)-1;
	size_t newSeedEndPos = (size_t)-1;
	
	size_t kmerSize = m_params.kmerLength;
	
	const int minRepeatFreq = coverage;
	double minFreqDiff = 0.5;
	
	std::string kmer = readSeq.substr(seedStartPos, kmerSize);
	int initKmerFreq = BWTAlgorithms::countSequenceOccurrences(kmer, m_params.indices);
	int prevKmerFreq = initKmerFreq;
	size_t startKmerFreq=0, endKmerFreq=0;


	// first kmer is a repeat
	if(initKmerFreq > minRepeatFreq)
	{
		newSeedStartPos = seedStartPos;
		startKmerFreq = initKmerFreq;
	}
	
	
	// identify breakpoints of large freq difference between two kmers	
	for(size_t i=seedStartPos+1 ; i+kmerSize-1 <= seedEndPos; i++)
	{
		kmer = readSeq.substr(i, kmerSize);
		int currKmerFreq = BWTAlgorithms::countSequenceOccurrences(kmer, m_params.indices);

		// std::cout << i << ": " << kmer << "\t" << currKmerFreq << "\n";

		// error kmers within repeats often lead to large freq diff
		bool isLargeFreqDiff = (currKmerFreq - prevKmerFreq)/(double)currKmerFreq > minFreqDiff;
		
		// PB36993_4517.fa, TTATGTAAGGAGTATTTGAT
		// some error kmers are with moderate frequency and no freq diff can be observed
		// pick up the first repeat kmer as starting point
		bool isRepeatKmer = (newSeedStartPos == (size_t)-1) && (currKmerFreq >= (int)minRepeatFreq);
		if(isLargeFreqDiff || isRepeatKmer)
		{
			// capture the highest repeat start pos
			bool isBetterRepeatKmer = (startKmerFreq!=0 && currKmerFreq > (int)startKmerFreq);
			if(newSeedStartPos == (size_t)-1 || isBetterRepeatKmer)
			{
				newSeedStartPos = i;
				startKmerFreq = currKmerFreq;
			}
		}
			
		// repeat end is reached
		// PB36993_4517.fa, AGGCTTGTCTGTAATCGGG
		isLargeFreqDiff = (prevKmerFreq - currKmerFreq)/(double)prevKmerFreq > minFreqDiff;
		if(isLargeFreqDiff /*|| currKmerFreq < minFreqDiff*/)
		{
			// do not enter unless start pos was found
			// if(newSeedStartPos != (size_t)-1)
			// {
				newSeedEndPos = i + kmerSize -2;
				endKmerFreq = prevKmerFreq;
				break;
			// }
		}
			
		prevKmerFreq = currKmerFreq;
	}
	
	if(newSeedStartPos == (size_t)-1)
	{
		newSeedStartPos = seedStartPos;
		startKmerFreq = initKmerFreq;
	}
		
	if(newSeedEndPos == (size_t)-1)
	{
		newSeedEndPos = seedEndPos;
		endKmerFreq = prevKmerFreq;
	}
/*
	std::string kmer = readSeq.substr(seedStartPos, kmerSize);
	int currKmerFreq = BWTAlgorithms::countSequenceOccurrences(kmer, m_params.indices);

	// std::cout << i << ": " << kmer << "\t" << currKmerFreq << "\n";

	bool isNonRepeatKmer = currKmerFreq < (int) coverage;
	if(isNonRepeatKmer)
	{
		newSeedStartPos = seedStartPos+1;
	}
			
	// repeat end is reached
	kmer = readSeq.substr(seedEndPos-kmerSize+1, kmerSize);
	currKmerFreq = BWTAlgorithms::countSequenceOccurrences(kmer, m_params.indices);
	isNonRepeatKmer = currKmerFreq < (int) coverage;
	if(isNonRepeatKmer)
	{
		newSeedEndPos = seedEndPos-1;
	}
				
	if(newSeedStartPos == (size_t)-1)
	{
		newSeedStartPos = seedStartPos;
	}
		
	if(newSeedEndPos == (size_t)-1)
	{
		newSeedEndPos = seedEndPos;
	}
	
	// std::cout << newSeedStartPos << "\t" << newSeedEndPos << "\n";
*/
	seedStartPos = newSeedStartPos;
	seedEndPos = newSeedEndPos;
	return;
}

//
//
PacBioHybridCorrectionPostProcess::PacBioHybridCorrectionPostProcess(std::ostream* pCorrectedWriter, std::ostream* pDiscardWriter, const PacBioHybridCorrectionParameters params):
	m_pCorrectedWriter(pCorrectedWriter),
	m_pDiscardWriter(pDiscardWriter),
	m_params(params),
	m_totalReadsLen(0),
	m_correctedLen(0),
	m_totalSeedNum(0),
	m_totalWalkNum(0),
	m_correctedNum(0),
	m_highErrorNum(0),
	m_exceedDepthNum(0),
	m_exceedLeaveNum(0),
	m_seedDis(0)
{
}

//
PacBioHybridCorrectionPostProcess::~PacBioHybridCorrectionPostProcess()
{	
	if(m_totalWalkNum>0 && m_totalReadsLen>0)
	{
		std::cout << std::endl;
		std::cout << "totalReadsLen: " << m_totalReadsLen << ", ";
		std::cout << "correctedLen: " << m_correctedLen << ", ratio: " 
			<< (float)(m_correctedLen)/m_totalReadsLen << "%." << std::endl;
		std::cout << "totalSeedNum: " << m_totalSeedNum << "." << std::endl;
		std::cout << "totalWalkNum: " << m_totalWalkNum << ", ";
		std::cout << "correctedNum: " << m_correctedNum << ", ratio: " 
			<< (float)(m_correctedNum*100)/m_totalWalkNum << "%." << std::endl;
		std::cout << "seedDis: " << (float)(m_seedDis)/m_totalWalkNum << "." << std::endl;
		//std::cout << "highErrorNum: " << m_highErrorNum << ", ratio: " << (float)(m_highErrorNum*100)/m_totalWalkNum << "%." << std::endl;
		//std::cout << "exceedDepthNum: " << m_exceedDepthNum << ", ratio: " << (float)(m_exceedDepthNum*100)/m_totalWalkNum << "%." << std::endl;
		//std::cout << "exceedLeaveNum: " << m_exceedLeaveNum << ", ratio: " << (float)(m_exceedLeaveNum*100)/m_totalWalkNum << "%." << std::endl;
	}
}


// Writting results for kmerize and validate
void PacBioHybridCorrectionPostProcess::process(const SequenceWorkItem& item, const PacBioHybridCorrectionResult& result)
{
	if(result.merge)
	{
		m_totalReadsLen += result.totalReadsLen;
		m_correctedLen += result.correctedLen;
		m_totalSeedNum += result.totalSeedNum;
		m_totalWalkNum += result.totalWalkNum;
		m_correctedNum += result.correctedNum;
		m_highErrorNum += result.highErrorNum;
		m_exceedDepthNum += result.exceedDepthNum;
		m_exceedLeaveNum += result.exceedLeaveNum;
		m_seedDis += result.seedDis;
		
		for(size_t i = 0 ; i < result.correctedPacbioStrs.size() ; i++)
		{
			SeqItem mergeRecord;
			std::stringstream ss;
			ss << item.read.id << "_" << i << "_" << result.correctedPacbioStrs[i].toString().length();
			mergeRecord.id = ss.str();
			mergeRecord.seq = result.correctedPacbioStrs[i];
			mergeRecord.write(*m_pCorrectedWriter);
		}
	}
	else
	{
		// write into discard.fa
		SeqItem mergeRecord;
		mergeRecord.id = item.read.id;
		mergeRecord.seq = item.read.seq;
		mergeRecord.write(*m_pDiscardWriter);
	}
}
