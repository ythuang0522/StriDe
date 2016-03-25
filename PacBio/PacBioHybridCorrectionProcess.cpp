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
#include "RollingPBSelfCTree.h"
#include "Timer.h"
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
	
	std::vector<SeedFeature> seedVec, pacbioCorrectedStrs;
	std::string readSeq = workItem.read.seq.toString();
		
	seedVec = seedingByDynamicKmer(readSeq);

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
		SeedFeature source = pacbioCorrectedStrs.back();
		SeedFeature target = seedVec.at(targetSeed);
		int dis_between_src_target = target.seedStartPos - seedVec.at(targetSeed-1).seedStartPos - seedVec.at(targetSeed-1).seedLength;
		int FMWalkReturnType;
		string mergedseq;
		FMWalkReturnType = extendBetweenSeeds(source, target, dis_between_src_target, &mergedseq);
		
		// record corrected pacbio reads string
		// FMWalk success
		if(FMWalkReturnType == 1)
		{
			size_t gainPos = source.seedLength;
			// have gain ground
			if(mergedseq.length() > gainPos)
			{
				string gainStr = mergedseq.substr(gainPos);
				pacbioCorrectedStrs.back().append(gainStr);
				result.correctedLen += gainStr.length();
			}
		}
		// FMWalk failure: 
		// 1. high error 
		// 2. exceed leaves
		// 3. exceed depth
		else
		{
			pacbioCorrectedStrs.push_back(target);
			result.correctedLen += target.seedLength;
		}
		
		// output information
		result.totalWalkNum++;			
		result.seedDis += dis_between_src_target;
		if(FMWalkReturnType == 1)
			result.correctedNum++;
		/*else if(FMWalkReturnType == -1)
			result.highErrorNum++;
		else if(FMWalkReturnType == -2)
			result.exceedDepthNum++;
		else if(FMWalkReturnType == -3)
			result.exceedLeaveNum++;*/
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
std::vector<SeedFeature> PacBioHybridCorrectionProcess::seedingByDynamicKmer(const string readSeq)
{
	std::vector<SeedFeature> seedVec;
	std::vector<int> seedEndPosVec;
	int maxKmerSize = m_params.kmerLength;
	int minKmerSize = m_params.minKmerLength;
	int kmerThreshold = m_params.seedKmerThreshold;
	
	// prevention of short reads
	if(readSeq.length() < maxKmerSize)
		return seedVec;
	
	// computing maximum interval of variously consecutive matches length (dk figure)
	std::vector<int> maxSeedInterval;
	for(int i = 0 ; i <= maxKmerSize ; i++)
		maxSeedInterval.push_back(2*3.8649*pow(2.7183,0.1239*i));
	
	// set dynamic kmer as largest kmer initially, 
	// which will reduce size when no seeds can be found within a maximum interval.
	int dynamicKmerSize = maxKmerSize;
		
	// search for solid kmers and group consecutive solids kmers into one seed
	// reduce kmer size if no seeds can be found within maxSeedInterval
	for(int i = 0 ; i+dynamicKmerSize <= readSeq.length() ; i++)
	{
		string kmer = readSeq.substr(i, dynamicKmerSize);
		int fwdKmerFreqs = BWTAlgorithms::countSequenceOccurrencesSingleStrand(kmer, m_params.indices);
		int rvcKmerFreqs = BWTAlgorithms::countSequenceOccurrencesSingleStrand(reverseComplement(kmer), m_params.indices);
		int kmerFreqs = fwdKmerFreqs+rvcKmerFreqs;
		
		// finding the seed
		if(kmerFreqs >= kmerThreshold)
		{
			int seedStartPos = i;
			int maxKmerFreq = kmerFreqs;
			
			// group consecutive solid kmers into one seed if possible
			for(i++ ; i+dynamicKmerSize <= readSeq.length() ; i++)
			{
				kmer = readSeq.substr(i, dynamicKmerSize);
				fwdKmerFreqs = BWTAlgorithms::countSequenceOccurrencesSingleStrand(kmer, m_params.indices);
				rvcKmerFreqs = BWTAlgorithms::countSequenceOccurrencesSingleStrand(reverseComplement(kmer), m_params.indices);
				kmerFreqs = fwdKmerFreqs+rvcKmerFreqs;
				maxKmerFreq = std::max(maxKmerFreq,kmerFreqs);
				if(kmerFreqs < kmerThreshold)
					break;
			}
			
			int seedEndPos = i+dynamicKmerSize-2;
			SeedFeature newSeed(seedStartPos, readSeq.substr(seedStartPos, seedEndPos-seedStartPos+1), true, dynamicKmerSize, kmerThreshold*4);
			seedVec.push_back(newSeed);
			seedEndPosVec.push_back(seedEndPos);
			//std::cout << ">" << seedStartPos << "." << newSeed.seedLength << "\n" << newSeed.seedStr << "\n";
			
			// jump to the index after new seed
			i = seedEndPos;
			dynamicKmerSize = maxKmerSize;
		}
		
		// reduce kmer size if no seeds can be found within a maximum interval
		int prevSeedEndPos = seedEndPosVec.empty() ? 0 : seedEndPosVec.back()+1;
		int distToPrevSeed = i + 1 - prevSeedEndPos;
		if(distToPrevSeed >= maxSeedInterval.at(dynamicKmerSize))
		{
			// reducing the kmer size
			if(dynamicKmerSize > minKmerSize)
			{
				i = prevSeedEndPos;
				dynamicKmerSize -= 2;
			}
			// no seeds can be found even using minKmerSize, 
			// so we reset start position of searching to give up the maximum interval of minKmerSize.
			else
			{
				seedEndPosVec.push_back(i);
				dynamicKmerSize = maxKmerSize;
			}
		}
	}
	
	return seedVec;
}

int PacBioHybridCorrectionProcess::extendBetweenSeeds(SeedFeature source, SeedFeature target, int dis_between_src_target, string* mergedseq)
{
	int FMWalkReturnType;
	int minOverlap = std::min(source.seedLength, target.seedLength);
		minOverlap = std::min(minOverlap, m_params.maxOverlap);
	PBHCFMWalkParameters PBHCFMWParams;
	PBHCFMWParams.indices = m_params.indices;
	PBHCFMWParams.maxOverlap = m_params.maxOverlap;
	PBHCFMWParams.SAThreshold = m_params.FMWKmerThreshold;
	PBHCFMWParams.disBetweenSrcTarget = dis_between_src_target;
	PBHCFMWParams.maxLeaves = m_params.maxLeaves;
	
	do
	{
		PBHCFMWParams.sourceSeed = source.seedStr;
		PBHCFMWParams.targetSeed = target.seedStr;
		PBHCFMWParams.minOverlap = minOverlap;
		SAIntervalPBHybridCTree SAITree(PBHCFMWParams);
		FMWalkReturnType = SAITree.mergeTwoReads(*mergedseq);
		if(FMWalkReturnType > 0)
		{
			assert((*mergedseq).empty() != true);
			string mergedseq2;
			PBHCFMWParams.sourceSeed = reverseComplement(target.seedStr);
			PBHCFMWParams.targetSeed = reverseComplement(source.seedStr);
			SAIntervalPBHybridCTree SAITree2(PBHCFMWParams);
			FMWalkReturnType = SAITree2.mergeTwoReads(mergedseq2);
			if((*mergedseq).length() == mergedseq2.length())
				break;
			else if(FMWalkReturnType > 0)
				FMWalkReturnType = -4;
		}
		minOverlap--;		
	}while(FMWalkReturnType < 0 && minOverlap >= m_params.minKmerLength);
	
	if(FMWalkReturnType < 0)
	{
		PBHCFMWParams.lowCoverageHighErrorMode = true;
		PBHCFMWParams.minOverlap = m_params.minKmerLength;
		PBHCFMWParams.maxOverlap = m_params.kmerLength;
		PBHCFMWParams.sourceSeed = source.seedStr;
		PBHCFMWParams.targetSeed = target.seedStr;
		SAIntervalPBHybridCTree SAITree(PBHCFMWParams);
		FMWalkReturnType = SAITree.mergeTwoReads(*mergedseq);
		if(FMWalkReturnType < 0)
		{
			PBHCFMWParams.sourceSeed = reverseComplement(target.seedStr);
			PBHCFMWParams.targetSeed = reverseComplement(source.seedStr);
			SAIntervalPBHybridCTree SAITree(PBHCFMWParams);
			FMWalkReturnType = SAITree.mergeTwoReads(*mergedseq);
			(*mergedseq) = reverseComplement((*mergedseq));
		}
	}
	
	return FMWalkReturnType;
}

//
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