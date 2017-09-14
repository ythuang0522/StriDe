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
#include "ShortReadOverlapTree.h"
#include "PacBioSelfCorrectionProcess.h"

using namespace std;

PacBioHybridCorrectionProcess::PacBioHybridCorrectionProcess(const PacBioHybridCorrectionParameters params):m_params(params)
{
}

PacBioHybridCorrectionProcess::~PacBioHybridCorrectionProcess()
{
}

// PacBio Hybrid Correction by Ya, v20170602.
PacBioHybridCorrectionResult PacBioHybridCorrectionProcess::PBHybridCorrection(SequenceWorkItem& workItem)
{
	PacBioHybridCorrectionResult result;
	
	// std::cout << workItem.read.id << endl;
	std::vector<SeedFeature> seedVec;
	SeedFeature pacbioCorrectedStrs;
	std::vector<bool> isPBPosCorrectedByHybridCorrection;
	std::string readSeq = workItem.read.seq.toString();
	seedVec = dynamicSeedingFromSR(readSeq);
	seedVec = filterErrorSRSeeds(seedVec);
	
	if(seedVec.size() >= 2)
	{
		result.correctedLen += seedVec.at(0).seedLength;
		pacbioCorrectedStrs = seedVec.at(0);
		// record position of PBHC read whether is corrected.
		// for(int i=0; i<seedVec.at(0).seedStartPos; i++)
			// isPBPosCorrectedByHybridCorrection.push_back(false);
		for(int i=0; i<seedVec.at(0).seedLength; i++)
			isPBPosCorrectedByHybridCorrection.push_back(true);
	}
	else
		// calling ChengWei's PB Self Correction if no seed.
		return PBSelfCorrection(workItem);
	
	// FMWalk for each pair of seeds
	for(size_t numFMWalk = 1 ; numFMWalk < seedVec.size() ; numFMWalk++)
	{
		FMWalkResult FMWResult;
		SeedFeature seedSource = pacbioCorrectedStrs,
			seedTarget = seedVec.at(numFMWalk);

		// in order to escape the error seed target potentially, 
		// we try to use the next seed as a new target.
		int tryNext = 0;
		for(int nextFMWalk = numFMWalk + tryNext; 
		// tryNext <= N, N means the maximum try times.
		tryNext <= 3 && 
		(FMWResult.typeFMWalkResult==-1||FMWResult.typeFMWalkResult==-2) && 
		nextFMWalk < seedVec.size(); 
		tryNext++, nextFMWalk = numFMWalk + tryNext)
		{
			seedTarget = seedVec.at(nextFMWalk);
			extendBetweenSeeds(readSeq, seedSource, seedTarget, FMWResult);
			// cout << numFMWalk << ": " <<  FMWResult.typeFMWalkResult <<"\n";
		}
		if(FMWResult.typeFMWalkResult > 0)
			numFMWalk += (tryNext - 1);
		seedTarget = seedVec.at(numFMWalk);
		
		// debug:
		// string seed = seedVec.at(numFMWalk-1).seedStr;
		// int fwdSF = BWTAlgorithms::countSequenceOccurrencesSingleStrand(seed, m_params.indices);
		// int rvcSF = BWTAlgorithms::countSequenceOccurrencesSingleStrand(reverseComplement(seed), m_params.indices);
		// int SRSF = fwdSF+rvcSF;
		// fwdSF = BWTAlgorithms::countSequenceOccurrencesSingleStrand(seed, m_params.PBindices);
		// rvcSF = BWTAlgorithms::countSequenceOccurrencesSingleStrand(reverseComplement(seed), m_params.PBindices);
		// int PBSF = fwdSF+rvcSF;
		// std::cout << ">" << numFMWalk-1
		// << ", pos:" << seedVec.at(numFMWalk-1).seedStartPos
		// << ", len:" << seedVec.at(numFMWalk-1).seedLength
		// << ", isPB:" << (seedVec.at(numFMWalk-1).isPBSeed?"yes":"no")
		// << ", SRSF:" << SRSF
		// << ", PBSF:" << PBSF
		// << ", GC&tandem%:" << (int)(GCAndTandemRatio(seed)*100) << "%"
		// << ", FMWRT:" << FMWResult.typeFMWalkResult << ".\n"
		// << seed << endl;
		// if(FMWResult.typeFMWalkResult==1)
		// {
			// cout << "raw length: "
			// << seedTarget.seedStartPos-seedSource.seedEndPos-1
			// << ", corrected length: "
			// << FMWResult.mergedSeq.substr(seedSource.seedLength).length()-seedTarget.seedLength
			// << "\n"
			// << readSeq.substr(seedSource.seedEndPos+1,seedTarget.seedEndPos-seedSource.seedEndPos)
			// << "\n"
			// << FMWResult.mergedSeq.substr(seedSource.seedLength)
			// << "\n";
		// }
		
		// record corrected pacbio reads string--
		// result.seedDis += target.seedStartPos-source.seedEndPos-1;
		// A. FMWalk success:
		if(FMWResult.typeFMWalkResult == 1)
		{
			assert(FMWResult.mergedSeq.length() > seedSource.seedLength);
			string extendedStr = FMWResult.mergedSeq.substr(seedSource.seedLength);
			pacbioCorrectedStrs.append(extendedStr);
			pacbioCorrectedStrs.seedStartPos = seedTarget.seedStartPos;
			pacbioCorrectedStrs.seedEndPos = seedTarget.seedEndPos;
			pacbioCorrectedStrs.isRepeat = seedTarget.isRepeat;
			pacbioCorrectedStrs.isPBSeed = seedTarget.isPBSeed;
			pacbioCorrectedStrs.isNextRepeat = seedTarget.isNextRepeat;
			pacbioCorrectedStrs.startBestKmerSize = seedTarget.startBestKmerSize;
			pacbioCorrectedStrs.endBestKmerSize = seedTarget.endBestKmerSize;
			result.correctedLen += extendedStr.length();
			
			// record position of PBHC read whether is corrected.
			for(int i=0; i<extendedStr.length(); i++)
				isPBPosCorrectedByHybridCorrection.push_back(true);
		}
		// B. FMWalk failure: 
		// 1. high error 
		// 2. exceed leaves
		// 3. exceed depth
		else
		{
			string extendedStr = readSeq.substr(seedSource.seedEndPos+1,seedTarget.seedEndPos-seedSource.seedEndPos);
			pacbioCorrectedStrs.append(extendedStr);
			pacbioCorrectedStrs.seedStartPos = seedTarget.seedStartPos;
			pacbioCorrectedStrs.seedEndPos = seedTarget.seedEndPos;
			pacbioCorrectedStrs.isRepeat = seedTarget.isRepeat;
			pacbioCorrectedStrs.isPBSeed = seedTarget.isPBSeed;
			pacbioCorrectedStrs.isNextRepeat = seedTarget.isNextRepeat;
			pacbioCorrectedStrs.startBestKmerSize = seedTarget.startBestKmerSize;
			pacbioCorrectedStrs.endBestKmerSize = seedTarget.endBestKmerSize;
			result.correctedLen += extendedStr.length();
			
			// record position of PBHC read whether is corrected.
			for(int i=0; i<(seedTarget.seedStartPos-seedSource.seedEndPos-1); i++)
				isPBPosCorrectedByHybridCorrection.push_back(false);
			for(int i=0; i<seedTarget.seedLength; i++)
				isPBPosCorrectedByHybridCorrection.push_back(true);
		}
		
		// output information
		result.totalWalkNum++;
		if(FMWResult.typeFMWalkResult == 1)
			result.correctedNum++;
		// else if(FMWalkReturnType == -1)
			// result.highErrorNum++;
		// else if(FMWalkReturnType == -2)
			// result.exceedDepthNum++;
		// else if(FMWalkReturnType == -3)
			// result.exceedLeaveNum++;
	}
	
	// result.strPBHC =
		// (seedVec.front().seedStartPos==0?"":readSeq.substr(0,seedVec.front().seedStartPos))+
		// pacbioCorrectedStrs.seedStr+
		// (seedVec.back().seedEndPos==(readSeq.length()-1)?"":readSeq.substr(seedVec.back().seedEndPos+1));
	// record position of PBHC read whether is corrected.
	// for(int i=seedVec.back().seedEndPos+1; i<readSeq.length(); i++)
			// isPBPosCorrectedByHybridCorrection.push_back(false);

	result.strPBHC = pacbioCorrectedStrs.seedStr;
	result.totalSeedNum = seedVec.size();
	result.totalReadsLen = readSeq.length();
	result.merge = true;
	
	// return result;
	
	assert(result.strPBHC.length()==isPBPosCorrectedByHybridCorrection.size());
	
	workItem.read.seqPBHC = result.strPBHC;
	workItem.read.isPBPosCorrectedByHybridCorrection = isPBPosCorrectedByHybridCorrection;
	// calling ChengWei's PB Self Correction to solve FMWalk failed in PB Hybrid Correction.
	return PBSelfCorrection(workItem);
}

// calling ChengWei's PB Self Correction to solve FMWalk failed in PB Hybrid Correction.
PacBioHybridCorrectionResult PacBioHybridCorrectionProcess::PBSelfCorrection(const SequenceWorkItem& workItem)
{
	PacBioSelfCorrectionParameters selfECParams;
	selfECParams.indices = m_params.PBindices;
	selfECParams.kmerLength = 17;
	selfECParams.maxLeaves = 32;
	selfECParams.minKmerLength = 13;
    selfECParams.idmerLength = 9;
    selfECParams.ErrorRate = 0.15;
	selfECParams.FMWKmerThreshold = 3;
	selfECParams.numOfNextTarget = 1;
	selfECParams.collectedSeeds = 5;
    selfECParams.PBcoverage = m_params.PBcoverage;
	selfECParams.isSplit = false;
	selfECParams.isFirst = false;
    selfECParams.DebugExtend = false;
    selfECParams.DebugSeed = false;
	selfECParams.maxSeedInterval = 500;
	PacBioSelfCorrectionProcess processPBSC(selfECParams);
	PacBioSelfCorrectionResult resultPBSC=processPBSC.PBSelfCorrection(workItem);
	PacBioHybridCorrectionResult resultPBHC;
	resultPBHC.merge=resultPBSC.merge;
	resultPBHC.strPBHC=resultPBSC.strPBSC;
	return resultPBHC;
}

// dynamic seeding from short reads, v20160517 by Ya.
// identify seeds by dynamic kmers from short reads, 
// which is suitable for the PacBio hybrid error correction,
// where repeat regions require large kmers and error-prone regions require small kmers.
std::vector<SeedFeature> PacBioHybridCorrectionProcess::dynamicSeedingFromSR(const string& readSeq)
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
	std::vector<float> kmerThresholdVec;
	kmerThresholdVec.resize(201+1,3);
	for(size_t kmerSize=0 ; kmerSize<=91 ; kmerSize++)
	{
		float kmerThresholdValue=(0.005*pow(kmerSize,2)-0.96*kmerSize+45.955)*((float)m_params.coverage/100);
		kmerThresholdVec.at(kmerSize)+=kmerThresholdValue;
	}
	
	// search for solid kmers as seeds
	for(size_t pos=0 ; pos+minKmerSize<=readSeq.length() ; pos++)
	{
		string kmer=readSeq.substr(pos,minKmerSize);
		BWTInterval fwdInterval=BWTAlgorithms::findInterval(m_params.indices.pRBWT, reverse(kmer));
		BWTInterval rvcInterval=BWTAlgorithms::findInterval(m_params.indices.pBWT, reverseComplement(kmer));
		size_t kmerFreqs = (fwdInterval.isValid()?fwdInterval.size():0) + (rvcInterval.isValid()?rvcInterval.size():0);
		size_t dynamicKmerSize=minKmerSize;
		size_t dynamicKmerThreshold=kmerThresholdVec.at(minKmerSize);

		// std::cout << pos << ":" << kmerFreqs << ":" << dynamicKmerSize << ":" 
			// << BWTAlgorithms::countSequenceOccurrencesSingleStrand(kmer, m_params.indices) << ":" 
			// << BWTAlgorithms::countSequenceOccurrencesSingleStrand(reverseComplement(kmer), m_params.indices) << "\n";
			
		if(kmerFreqs<dynamicKmerThreshold)
		{
			// In large sequencing gaps (>7kb), no seeds can be found in Illumina index
			// If no seed is found within PBSearchDepth, 
			// seeds from PB index will be serached instead
			// int prevSeedEndPos=seedEndPosVec.empty()?0:(seedEndPosVec.back()+1);
			// int distToPrevSeed=pos+1-prevSeedEndPos;

			// skip repeat regions
			// bool isWithinRepeat = pos - prevRepeatEnd < m_params.PBSearchDepth;
			// if(/*!isWithinRepeat &&*/ distToPrevSeed>=m_params.PBSearchDepth)
			// {
				// if(seedingByPacBio(readSeq,seedVec,seedEndPosVec,prevSeedEndPos))
				// if(dynamicSeedingFromPB(readSeq,seedVec,seedEndPosVec,prevSeedEndPos))
					// ;//cout << seedVec.back().seedStr << " " << seedVec.back().seedStartPos << " " << BWTAlgorithms::countSequenceOccurrencesSingleStrand(seedVec.back().seedStr, m_params.PBindices)+BWTAlgorithms::countSequenceOccurrencesSingleStrand(reverseComplement(seedVec.back().seedStr), m_params.PBindices) << endl;
				// else
					// seedEndPosVec.push_back(pos);
				// cout << seedEndPosVec.size() << endl;
				// pos=seedEndPosVec.back();
			// }
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
			kmerFreqs = (fwdInterval.isValid()?fwdInterval.size():0) + (rvcInterval.isValid()?rvcInterval.size():0);
			
			dynamicKmerSize++;			
			assert(dynamicKmerSize<=kmerThresholdVec.size());
			dynamicKmerThreshold=kmerThresholdVec.at(dynamicKmerSize);
			
			if(kmerFreqs>=dynamicKmerThreshold)
				maxKmerFreq=kmerFreqs;
			else
			{
				dynamicKmerSize--;
				dynamicKmerThreshold=kmerThresholdVec.at(dynamicKmerSize);
				break;
			}
		}
		
		// cout << dynamicKmerSize << "\n";
		
		// small-sized seed has 50% chance of errors in C elegans, 
		// skip only if there is another seed nearby 30bp
		//if( (pos-1-minKmerSize-seedStartPos) < 2 && !seedVec.empty() && pos-1-minKmerSize-seedVec.back().seedEndPos <= 30)
		//if(!seedVec.empty() && pos-1-minKmerSize-seedVec.back().seedEndPos <= 30)
		//	continue;
	
		size_t seedEndPos=pos-1;
		
		// cut bases off the end of a seed, 
		// if below a threshold kmer frequency.
		for(size_t i=seedEndPos-minKmerSize+1 ; i>seedStartPos ; i--)
		{
			string kmer=readSeq.substr(i,minKmerSize);
			BWTInterval fwdInterval=BWTAlgorithms::findInterval(m_params.indices.pRBWT, reverse(kmer));
			BWTInterval rvcInterval=BWTAlgorithms::findInterval(m_params.indices.pBWT, reverseComplement(kmer));
			size_t kmerFreqs = (fwdInterval.isValid()?fwdInterval.size():0) + (rvcInterval.isValid()?rvcInterval.size():0);
			size_t kmerThreshold=kmerThresholdVec.at(minKmerSize);
			if(kmerFreqs>=kmerThreshold)
			{
				seedEndPos=i+minKmerSize-1;
				break;
			}
			dynamicKmerSize--;
		}
		
		// if seed frequency in PB great tha in SR, 
		// we skipped it.
		if(dynamicKmerSize==minKmerSize && 
			BWTAlgorithms::countSequenceOccurrences(readSeq.substr(seedStartPos, seedEndPos-seedStartPos+1), m_params.PBindices)>maxKmerFreq)
		{
			// pos=seedStartPos;
			continue;
		}

		// repeat seeds are less accurate at boundary and should be trimmed 
		if(maxKmerFreq >= m_params.coverage*4)
			trimRepeatSeed(readSeq,m_params.coverage,seedStartPos,seedEndPos);
		
		// bool isTooCloseToPrevSeed = seedEndPos - prevRepeatEnd <= prevRepeatFreq*0.02+100;
		// if(isTooCloseToPrevSeed)
			// continue;
				
		// super repeat seeds with frequency > 2000 are troublesome, often lead to -3 but no good solution so far, mark first
		bool isSuperRepeat = maxKmerFreq>(m_params.coverage*15)?true:false;
		SeedFeature newSeed(seedStartPos, readSeq.substr(seedStartPos, seedEndPos-seedStartPos+1), isSuperRepeat, dynamicKmerSize, m_params.PBcoverage/2);
		// newSeed.estimateBestKmerSize(m_params.PBindices.pBWT);

		// Some high-GC error seeds have unexpected large freq
		// size_t PBkmerFreq = BWTAlgorithms::countSequenceOccurrences(newSeed.seedStr, m_params.PBindices);
		// if(PBkmerFreq > m_params.PBcoverage/4)
			// continue;

		// skip low-complexity sequencing errors of PacBio
		if(!isLowComplexity(newSeed.seedStr,0.9))
		{
			// std::cout << "SR> " << newSeed.seedStartPos << ": " 
				// << newSeed.seedLength << ": "
				// << maxKmerFreq << ": "
				// << GCAndTandemRatio(newSeed.seedStr) << ": "
				// << newSeed.seedStr << "\n";
			seedVec.push_back(newSeed);
		}
			seedEndPosVec.push_back(seedEndPos);
			
		// jump to the index after new seed
		pos=seedEndPos;
	}
	
	return seedVec;
}

// dynamic seeding from pacbio reads, v20160517 by Ya.
// identify seeds by dynamic kmers from pacbio reads, 
// which is suitable for the PacBio hybrid error correction,
// where repeat regions require large kmers and error-prone regions require small kmers.
bool PacBioHybridCorrectionProcess::dynamicSeedingFromPB(const string& readSeq, std::vector<SeedFeature>& seedVec, std::vector<int>& seedEndPosVec, size_t prevEndPos)
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
		size_t kmerFreqs = (fwdInterval.isValid()?fwdInterval.size():0) + (rvcInterval.isValid()?rvcInterval.size():0);
		size_t dynamicKmerThreshold=kmerThreshold.at(dynamicKmerSize);
		
		if(kmerFreqs<dynamicKmerThreshold)
			continue;
		
		// if the seed is in repeat region, we increase the PB seed threshold.
		bool isRepeatRegion = 
		(seedVec.size()!=0 && 
		BWTAlgorithms::countSequenceOccurrences(seedVec.back().seedStr, m_params.indices)>m_params.coverage*4)
		?true:false;
		if(isRepeatRegion && kmerFreqs<dynamicKmerThreshold*2)
			continue;
		// cout << kmer << " " << dynamicKmerThreshold << " " << dynamicKmerThreshold*2 << endl;

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
			kmerFreqs = (fwdInterval.isValid()?fwdInterval.size():0) + (rvcInterval.isValid()?rvcInterval.size():0);
			
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

		// skip repeat seeds
		if(maxKmerFreq >= m_params.PBcoverage*2)
			continue;
		
		// skip repeat seeds, 2nd.
		if(isRepeatRegion && kmerFreqs<dynamicKmerThreshold*2)
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
		// skip GC and tandem ratio >= 0.9 seeds
		if(!isLowComplexity(newSeed.seedStr,0.8) && GCAndTandemRatio(newSeed.seedStr)<0.9)
		{
			// std::cout << "PB> " << newSeed.seedStartPos << ": " 
				// << newSeed.seedLength << ": "
				// << maxKmerFreq << ": "
				// << GCAndTandemRatio(newSeed.seedStr) << ": "
				// << newSeed.seedStr << "\n";
			newSeed.isPBSeed = true;
			seedVec.push_back(newSeed);
			seedEndPosVec.push_back(seedEndPos);
			return true;
		}
	}
	
	return false;
}

// filter error seeds of short reads, v20161108 by Ya.
// seed frequency: average k-mer frequency in a seed.
// if seed frequency in short reads FM-index is below all average k-mer frequency 
// and is between high frequency seeds , the seed will be deleted.
std::vector<SeedFeature> PacBioHybridCorrectionProcess::filterErrorSRSeeds(std::vector<SeedFeature>& seedVec)
{
	if(seedVec.size() < 3)
		return seedVec;
	
	std::vector<SeedFeature> newSeedVec;
	std::vector<float> seedFreqsVec;
	
	// caculate all average k-mer frequency
	float allSeedsAvgKmerFreqs=0;
	int allSeedsKmerFreqs=0, count=0;
	for(int i=0 ; i<seedVec.size() ; i++)
	{
		// cout << i << "-- " << seedVec.at(i).seedStartPos << "\n";
		
		// seeds from short reads.
		if(!seedVec.at(i).isPBSeed)
		{
			int j=0, seedKmerFreqs=0;
			for(; j<=seedVec.at(i).seedLength-m_params.minKmerLength ; j++)
			{
				string kmer = seedVec.at(i).seedStr.substr(j, m_params.minKmerLength);
				int kmerFreqs = 
				BWTAlgorithms::countSequenceOccurrencesSingleStrand(kmer, m_params.indices) +
				BWTAlgorithms::countSequenceOccurrencesSingleStrand(reverseComplement(kmer), m_params.indices);
				seedKmerFreqs += kmerFreqs;
				allSeedsKmerFreqs += kmerFreqs;
				count++;
				// cout << kmerFreqs << " ";
			}
			// seedFreqsVec.push_back((float)seedKmerFreqs/j);
			// seedFreqsVec is a prefix sum array of kmer freqs for computing interval average kmer freq later
			if(i > 0)
				seedFreqsVec.push_back((float)seedKmerFreqs/j + seedFreqsVec.at(i-1));
			else
				seedFreqsVec.push_back((float)seedKmerFreqs/j);
			// cout << (float)seedKmerFreqs/j << " ";
		}
		// seeds from pacbio reads.
		else
			seedFreqsVec.push_back(0);
		// cout << endl;
	}
	allSeedsAvgKmerFreqs = (float)allSeedsKmerFreqs/count;
	
	int SRSearchDepth = 0.8*m_params.PBSearchDepth;
	// cout << SRSearchDepth << endl;
	// cout << seedFreqsVec.size() << endl;
	newSeedVec.push_back(seedVec.at(0));
	for(int i=1 ; i<(seedVec.size()-1) ; i++)
	{
		bool isErrorSeed=false;
		
		// compute average kmer freq within a window centered at i
		const size_t window = 10;
		// special case of short reads
		if(seedFreqsVec.size() <= window)
			allSeedsAvgKmerFreqs = (float)seedFreqsVec.back()/seedFreqsVec.size();
		// left boundary case < window/2
		else if(i <= window/2)
			allSeedsAvgKmerFreqs = (float)seedFreqsVec.at(window-1)/window;
		// right boundary case > window/2
		else if(i >= seedFreqsVec.size()-window/2)
			allSeedsAvgKmerFreqs = (float)(seedFreqsVec.back()-seedFreqsVec.at(seedFreqsVec.size()-window-1))/window;
		// main case: window centered at i
		else
			allSeedsAvgKmerFreqs = (float)(seedFreqsVec.at(i+window/2-1)-seedFreqsVec.at(i-window/2))/window;
						
		float currSeedFreq = seedFreqsVec.at(i)-seedFreqsVec.at(i-1) ;
		// cout << i << ": " << currSeedFreq << "-> " << allSeedsAvgKmerFreqs << "\n";
		
		if(!seedVec.at(i).isPBSeed && currSeedFreq<allSeedsAvgKmerFreqs/2)
		{
			bool isPreSeedHighFreqs=false, isSufSeedHighFreqs=false;
			for(int j=0 ; 
			j<seedVec.size() &&
			((int)seedVec.at(j).seedStartPos-(int)seedVec.at(i).seedStartPos)<=SRSearchDepth;
			j++)
			{
				float SeedJFreq = seedFreqsVec.at(j)-(j>0?seedFreqsVec.at(j-1):0);
				if(!seedVec.at(j).isPBSeed && SeedJFreq > allSeedsAvgKmerFreqs*1.2)
				{
					// pre seeds
					if(j<i &&
						((int)seedVec.at(i).seedStartPos-(int)seedVec.at(j).seedStartPos)<=SRSearchDepth)
					{
						isPreSeedHighFreqs=true;
						// cout << j << " ";
					}
					// suf seeds
					else if(j>i)
					{
						isSufSeedHighFreqs=true;
						// cout << j << " ";
					}
				}
			}
			if(isPreSeedHighFreqs && isSufSeedHighFreqs)
				isErrorSeed=true;
			
			// if(isErrorSeed)
				// cout << "." << i << seedVec.at(i).seedStr << "\n";
		}
		// cout << endl;
		
		if(!isErrorSeed)
			newSeedVec.push_back(seedVec.at(i));
	}
	newSeedVec.push_back(seedVec.at(seedVec.size()-1));
	
	return newSeedVec;
}

void PacBioHybridCorrectionProcess::extendBetweenSeeds(std::string& readSeq, SeedFeature& seedSource, SeedFeature& seedTarget, FMWalkResult& FMWResult)
{
	// FMWalk params initialized
	FMWalkParameters FMWParams;
	// if(seedTarget.seedStr == "ATTAATCAATAAAATTTACGATTATT")
		// FMWParams.debugMode=true;
	FMWParams.indices = m_params.indices;
	FMWParams.maxOverlap = m_params.maxOverlap;
	FMWParams.SAThreshold = m_params.FMWKmerThreshold;
	FMWParams.maxLeaves = m_params.maxLeaves;
	
	// compute minOverlap from min of seedLength and maxOverlap
	FMWParams.minOverlap = std::min(seedSource.seedLength, seedTarget.seedLength);
	FMWParams.minOverlap = std::min(FMWParams.minOverlap, m_params.maxOverlap);
	// fix a bug, corrected PB pos may more short than raw PB.
	FMWParams.minOverlap = std::min(FMWParams.minOverlap, (int)seedSource.seedEndPos+1);
	
	FMWParams.strSourceSeed = seedSource.seedStr;
	FMWParams.strTargetSeed = seedTarget.seedStr;
	FMWParams.disBetweenSrcTarget = seedTarget.seedStartPos-seedSource.seedEndPos-1;
	FMWParams.rawPBStrBetweenSrcTargetWith2Minoverlap = 
		readSeq.substr(seedSource.seedEndPos+1-FMWParams.minOverlap, FMWParams.disBetweenSrcTarget+2*FMWParams.minOverlap);
	
	// FMWalk 1st: Correction by FM-index extension from source to target
	SAIntervalPBHybridCTree SAITree(FMWParams);
	SAITree.mergeTwoSeeds(FMWResult);
	if(FMWResult.typeFMWalkResult==1)
		return;
	
	int iniTypeFMWalkResult = FMWResult.typeFMWalkResult;
	
	// FMWalk 2nd: Correction by FM-index extension from target to source
	if(FMWResult.typeFMWalkResult==-1||FMWResult.typeFMWalkResult==-2)
	{
		FMWParams.strSourceSeed = reverseComplement(seedTarget.seedStr);
		FMWParams.strTargetSeed = reverseComplement(seedSource.seedStr);
		FMWParams.rawPBStrBetweenSrcTargetWith2Minoverlap = reverseComplement(FMWParams.rawPBStrBetweenSrcTargetWith2Minoverlap);
		SAIntervalPBHybridCTree SAITree2(FMWParams);
		SAITree2.mergeTwoSeeds(FMWResult);
		if(FMWResult.typeFMWalkResult==1)
		{
			FMWResult.mergedSeq=reverseComplement(FMWResult.mergedSeq);
			return;
		}
	}
	
	return;

	// FMWalk 3rd: MSA Correction
	const std::string query = 
		readSeq.substr(seedSource.seedEndPos-m_params.PBKmerLength+1,
			seedTarget.seedEndPos-seedSource.seedEndPos+m_params.PBKmerLength);

	// Self correction by aligning all reads having seeds with query
	MultipleAlignment maquery = LongReadOverlap::buildMultipleAlignment(query,
												m_params.PBKmerLength, 
												m_params.PBKmerLength,
												query.length()/10, 
												0.73,	// alignment identity < 0.7 are often false positive repeats
												m_params.PBcoverage,
												m_params.PBindices);
	
	// skip insufficient number of overlapping reads for correction
	if(maquery.getNumRows() <= 3)
	{
		FMWResult.typeFMWalkResult = iniTypeFMWalkResult;
		return;
	}
	
	std::string consensus = maquery.calculateBaseConsensus(100000, -1);
	FMWResult.mergedSeq = seedSource.seedStr + consensus.substr(m_params.PBKmerLength);
	FMWResult.typeFMWalkResult = 1;
	return;
}

// 3rd correction:
// try correction by MSA using FM-index of low quality long reads for sequencing gaps
int PacBioHybridCorrectionProcess::MSACorrection(FMWalkParameters FMWParams, std::string& readSeq, SeedFeature& source, SeedFeature& target, bool& isSequencingGap, FMWalkResult& FMWResult)
{
	// mimic the overlap correction process
	// const std::string query = source.seedStr.substr(source.seedLength - m_params.PBKmerLength)+strBetweenSrcTarget.substr(10, dis_between_src_target)+target.seedStr;

	// switch to best kmer in PB index
	const std::string query = 
		readSeq.substr(source.seedEndPos-source.endBestKmerSize+1,
			target.seedEndPos-source.seedEndPos-1+source.endBestKmerSize);
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
		return 0;
	// if(target.seedStartPos == 14013)
		// cout << "\n" << query << "\n" << maquery.getNumRows() << "\n"
		// << sourceKmerLength << " " << targetKmerLength << "\n";
		// maquery.print(120);
	std::string consensus = maquery.calculateBaseConsensus(100000, -1);
	FMWResult.mergedSeq = source.seedStr + consensus.substr(m_params.PBKmerLength);
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

// return is low complexity of seq, v20160525 by Ya.
bool PacBioHybridCorrectionProcess::isLowComplexity(std::string& seq, const float& ratioThreshold)
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

// return GC and tandem nucleotides ratio of seq, v20161011 by Ya.
float PacBioHybridCorrectionProcess::GCAndTandemRatio(std::string& seq)
{
	size_t seqLen=seq.length();
	size_t countIsolatedAT=0;

	for(size_t i=0 ; i<seqLen ; i++)
	{
		if(seq[i]=='A'||seq[i]=='T')
		{
			if(i==0 && seq[i]!=seq[i+1])
				countIsolatedAT++;
			else if(i==(seqLen-1) && seq[i]!=seq[i-1])
				countIsolatedAT++;
			else if(seq[i]!=seq[i+1] && seq[i]!=seq[i-1])
				countIsolatedAT++;
		}
	}

	return (float)(seqLen-countIsolatedAT)/seqLen ;
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
			<< (float)(m_correctedLen*100)/m_totalReadsLen << "%." << std::endl;
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
		
		// for(size_t i = 0 ; i < result.correctedPacbioStrs.size() ; i++)
		{
			SeqItem mergeRecord;
			std::stringstream ss;
			ss << item.read.id << "_" << item.read.seq.length() << "_" << result.strPBHC.toString().length();
			mergeRecord.id = ss.str();
			mergeRecord.seq = result.strPBHC;
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