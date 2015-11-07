///-----------------------------------------------
// Copyright 2015 National Chung Cheng University
// Written by Yao-Ting Huang
// Released under the GPL
//-----------------------------------------------
//
// PacBioCorrectionProcess.cpp - Self-correction or hybrid correction using FM-index walk for PacBio reads
//
#include "PacBioCorrectionProcess.h"
#include "CorrectionThresholds.h"
#include "HashMap.h"
#include <iomanip>
#include "SAIPBSelfCTree.h"
#include "SAIPBHybridCTree.h"
using namespace std;


PacBioCorrectionProcess::PacBioCorrectionProcess(const PacBioCorrectionParameters params) : m_params(params)
{
}

PacBioCorrectionProcess::~PacBioCorrectionProcess()
{

}


// PacBio Self Correction by Ya, v20151001.
PacBioCorrectionResult PacBioCorrectionProcess::PBSelfCorrection(const SequenceWorkItem& workItem)
{	
	PacBioCorrectionResult result;
	
	std::vector<std::pair<int, std::string> > seeds, pacbioCorrectedStrs;
	std::string readSeq = workItem.read.seq.toString();
	
	// find seeds
	seeds = searchingSeedsUsingSolidKmer(readSeq);	
	result.totalSeedNum = seeds.size();
	
	// initialize corrected pacbio reads string
	if(seeds.size() >= 1)
	{
		result.correctedLen += seeds[0].second.length();
		pacbioCorrectedStrs.push_back(seeds[0]);
	}
	else
	{
		result.merge = false;
		return result;
	}
	
	for(size_t targetSeed = 1 ; targetSeed < seeds.size() ; targetSeed++)
	{
		size_t largeKmerSize = m_params.kmerLength;
		size_t smallKmerSize = m_params.minKmerLength;

		int FMWalkReturnType = -1;
		std::string mergedseq;
		std::pair<int,std::string> source = pacbioCorrectedStrs.back();
		
		// Multiple targets will be tested for FM-index walk from source to target, until m_params.downward times.
		for(int nextTargetSeed = 0 ; nextTargetSeed < m_params.downward && targetSeed + nextTargetSeed < seeds.size() ; nextTargetSeed++)
		{
			std::cout << "======= " << result.totalWalkNum << " =======\n";

			/********* PingYe's implementation  ***********/
			// // <distance between targets, read(ATCG)>
			// std::vector<std::pair<int, std::string> > targets;
			
			// // Multple targets are collected into targets vector, but for what ???
			// for(size_t collectedSeed = 0, currentTarget = targetSeed + nextTargetSeed ; 
			// collectedSeed < (size_t) m_params.collectedSeeds && (currentTarget + collectedSeed) < seeds.size() ; collectedSeed++)
			// {
			// int dis_between_targets = seeds[currentTarget+collectedSeed].first - seeds[targetSeed-1].first - seeds[targetSeed-1].second.length();
			// std::string target_str = seeds[currentTarget+collectedSeed].second;
			// targets.push_back(make_pair(dis_between_targets, target_str));
			// }
			

			// // FM-index extension for correction
			// SAIPBSelfCorrectTree SAITree(m_params.indices.pBWT, m_params.indices.pRBWT, m_params.FMWKmerThreshold);

			// SAITree.addHashFromPairedSeed(source.second, targets, smallKmerSize);

			// int maxLength = 1.4*(targets[0].first+40) + source.second.length() + smallKmerSize;
			// int minLength = 0.6*(targets[0].first-40) + source.second.length() + smallKmerSize;
			// size_t expectedLength = source.second.length() + targets[0].first + smallKmerSize;
			
			// FMWalkReturnType = SAITree.mergeTwoSeedsUsingHash(source.second, targets[0].second, mergedseq, smallKmerSize, m_params.maxLeaves,
			// minLength, maxLength, expectedLength);

			// std::cout << pacbioCorrectedStrs.size()-1 << ": " //<< reverseComplement((source.second).substr(source.second.length()-minOverlap)) << " " << reverseComplement((target.second).substr(0, minOverlap)) << " "
			// //<< (source.second).substr(source.second.length()-minOverlap) << " " << (target.second).substr(0, minOverlap) << " "
			// << source.first << "-" << source.first+source.second.length()-1 <<  ":" << source.second.length() << ", " 
			// <<	seeds[targetSeed+nextTargetSeed].first << "-" << seeds[targetSeed+nextTargetSeed].first+seeds[targetSeed+nextTargetSeed].second.length()-1 <<  ":" << seeds[targetSeed+nextTargetSeed].second.length() << ", dis: "
			// << targets[0].first << ". " << FMWalkReturnType << ".\n";
			
			
			/********* YT's implementation  ***********/					
			size_t currTargetIndex = targetSeed+nextTargetSeed;
			std::string targetStr = seeds[currTargetIndex].second;
			int dis_between_src_target = seeds[currTargetIndex].first - seeds[targetSeed-1].first - seeds[targetSeed-1].second.length();
			
			// Estimate upper/lower/expected bounds of search depth
			int maxLength = 1.2*(dis_between_src_target+20) + source.second.length() + smallKmerSize;
			int minLength = 0.8*(dis_between_src_target-20) + source.second.length() + smallKmerSize;
			if(minLength<0) minLength = 0;
			
			size_t expectedLength = dis_between_src_target + source.second.length() + smallKmerSize;

			// std::cout << "before build\t" << expectedLength << "\t" << maxLength << "\t" << minLength << "\n";
			SAIPBSelfCorrectTree SAITree(m_params.indices.pBWT, m_params.indices.pRBWT, m_params.FMWKmerThreshold);
			
			// Collect local kmer frequency from source (1st time) and left seeds
			const int seedMoveStepSize = 2;
			const size_t maxkmerCount = 2;
			int seedIdx=source.second.length()-largeKmerSize - (nextTargetSeed*seedMoveStepSize*maxkmerCount);
			size_t kmerCount = 0;
			for(; seedIdx >= 0 && kmerCount < maxkmerCount; seedIdx-=seedMoveStepSize, kmerCount++)
			{
				std::string srcStr = source.second.substr(seedIdx, largeKmerSize);
				SAITree.addHashFromSingleSeed(srcStr, smallKmerSize, maxLength);
			}

			// Collect local kmer frequency from right targets
			kmerCount = 0;			
			for(int seedIdx=targetStr.length()-largeKmerSize; seedIdx >= 0 && kmerCount < 2; seedIdx-=2, kmerCount++)				
			{
				std::string destStr = targetStr.substr(seedIdx, largeKmerSize);			
				SAITree.addHashFromSingleSeed(reverseComplement(destStr), smallKmerSize, maxLength);
			}
			
			// Extension using local kmer hash with smaller kmer size = m_params.minKmerLength
			FMWalkReturnType = SAITree.mergeTwoSeedsUsingHash(source.second, targetStr, mergedseq, smallKmerSize, m_params.maxLeaves,
			minLength, maxLength, expectedLength);
			
			std::cout << pacbioCorrectedStrs.size()-1 << ": " //<< reverseComplement((source.second).substr(source.second.length()-minOverlap)) << " " << reverseComplement((target.second).substr(0, minOverlap)) << " "
			//<< (source.second).substr(source.second.length()-minOverlap) << " " << (target.second).substr(0, minOverlap) << " "
			<< source.first << "-" << source.first+source.second.length()-1 <<  ":" << source.second.length() << ", " 
			<<	seeds[targetSeed+nextTargetSeed].first << "-" << seeds[targetSeed+nextTargetSeed].first+seeds[targetSeed+nextTargetSeed].second.length()-1 <<  ":" << seeds[targetSeed+nextTargetSeed].second.length() << ", dis: "
			<< dis_between_src_target << ". " << FMWalkReturnType << ".\n";
			

			// FMWalk success
			if(FMWalkReturnType > 0)
			{				
				size_t gainPos = source.second.length();
				assert(mergedseq.length() > gainPos);
				std::string gainStr = mergedseq.substr(gainPos);
				pacbioCorrectedStrs.back().second += gainStr;
				result.correctedLen += gainStr.length();
				result.correctedNum++;
				result.seedDis += seeds[targetSeed+nextTargetSeed].first - source.first - source.second.length();
				targetSeed = targetSeed + nextTargetSeed;
				break;
			}
			
			
			if(FMWalkReturnType==-3)
			{
				smallKmerSize++;
				// nextTargetSeed--;
			}
		}
		
		// All multiple targets failure: 
		// 1. high error 
		// 2. exceed leaves
		// 3. exceed depth
		// Then what??
		if(FMWalkReturnType < 0)
		{
			// seed distance
			result.seedDis += seeds[targetSeed].first - source.first - source.second.length();
			result.correctedLen += seeds[targetSeed].second.length();
			// if(round != 1)
			// {
			// not cut off
			pacbioCorrectedStrs.back().second += readSeq.substr(source.first + source.second.length(), 
			seeds[targetSeed].first + seeds[targetSeed].second.length() - source.first - source.second.length());
			// }
			
			// else if(round == 1)
			// {
			// cut off
			// pacbioCorrectedStrs.push_back(seeds[targetSeed]);
			// }
			
			// statistics of FM extension
			if(FMWalkReturnType == -1)
			result.highErrorNum++;
			else if(FMWalkReturnType == -2)
			result.exceedDepthNum++;
			else if(FMWalkReturnType == -3)
			result.exceedLeaveNum++;
		}
		result.totalWalkNum++;
	}
	
	result.merge = true;
	result.totalReadsLen = readSeq.length();
	for(size_t result_count = 0 ; result_count < pacbioCorrectedStrs.size() ; result_count++)
	result.correctedPacbioStrs.push_back(pacbioCorrectedStrs[result_count].second);
	return result;
}

std::vector<std::pair<int, std::string> > PacBioCorrectionProcess::searchingSeedsUsingSolidKmer(const std::string readSeq)
{
	std::vector<std::pair<int, std::string> > seeds;
	int kmerLen = m_params.kmerLength, 
	kmerThreshold = m_params.seedKmerThreshold, 
	readLen = readSeq.length();

	if(readLen >= kmerLen)
	{
		for(int i = 0 ; i+kmerLen <= readLen ; i++)
		{
			std::string kmer = readSeq.substr(i, kmerLen);
			// int kmerFreqs = 0;
			// kmerFreqs += BWTAlgorithms::countSequenceOccurrencesSingleStrand(kmer, m_params.indices);
			// kmerFreqs += BWTAlgorithms::countSequenceOccurrencesSingleStrand(reverseComplement(kmer), m_params.indices);
			int kmerFreqs = BWTAlgorithms::countSequenceOccurrences(kmer, m_params.indices);
			
			// std::cout << i << ": " << kmerFreqs << "\n";
			
			if(kmerFreqs >= kmerThreshold)
			{
				int seedStartPos = i, 
				seedLen = 0;
				
				// Group consecutive solid kmers into one seed if possible
				for(i++ ; i+kmerLen <= readLen ; i++)
				{
					// std::cout << i << ": " << kmerFreqs << "\n";
					kmer = readSeq.substr(i, kmerLen);
					// kmerFreqs = 0;
					// kmerFreqs += BWTAlgorithms::countSequenceOccurrencesSingleStrand(kmer, m_params.indices);
					// kmerFreqs += BWTAlgorithms::countSequenceOccurrencesSingleStrand(reverseComplement(kmer), m_params.indices);
					kmerFreqs = BWTAlgorithms::countSequenceOccurrences(kmer, m_params.indices);;
					// std::cout << kmer << ": " << kmerFreqs << "\n";
					
					if(kmerFreqs >= kmerThreshold)
					seedLen++;
					else
					break;
				}
				
				seeds.push_back(make_pair(seedStartPos, readSeq.substr(seedStartPos, seedLen+kmerLen)));
				i = i - 2 + kmerLen;
			}
		}
	}
	
	return seeds;
}

// PacBio Hybrid Correction by Ya, v20150617.
PacBioCorrectionResult PacBioCorrectionProcess::PBHybridCorrection(const SequenceWorkItem& workItem)
{
	PacBioCorrectionResult result;
	
	// get parameters
	string readSeq = workItem.read.seq.toString();
	string corReadSeq = readSeq;
	vector<string> pacbioCorrectedStrs;	
	vector<pair<int,string> > seeds;
	
	for(int round = 3 ; round > 0 ; round--)
	{
		// initialize
		result.correctedLen = 0;
		pacbioCorrectedStrs.clear();
		seeds.clear();
		
		// find seeds
		seeds = findSeedsUsingDynamicKmerLen(corReadSeq);
		
		// initialize corrected pacbio reads string
		if(seeds.size() >= 1)
		{
			result.correctedLen += seeds[0].second.length();
			pacbioCorrectedStrs.push_back(seeds[0].second);
		}
		else
		{
			result.merge = false;
			return result;
		}
		
		// FMWalk for each pair of seeds
		for(size_t i = 1 ; i < seeds.size() ; i++)
		{			
			int needWalkLen = (seeds[i].first - seeds[i-1].first - seeds[i-1].second.length());
			int FMWalkReturnType = -1, minOverlap;
			string mergedseq;
			pair<int,string> source = seeds[i-1];
			pair<int,string> target = seeds[i];
			
			if(source.second.length() >= (size_t) m_params.maxOverlap && target.second.length() >= (size_t) m_params.maxOverlap)
			minOverlap = m_params.maxOverlap-2;
			else
			{
				if(source.second.length() <= target.second.length())
				minOverlap = source.second.length();
				else
				minOverlap = target.second.length();
			}
			
			FMWalkReturnType = solveHighError(source, target, minOverlap, needWalkLen, &mergedseq);
			
			// record corrected pacbio reads string
			// FMWalk success
			if(FMWalkReturnType == 1)
			{
				size_t gainPos = source.second.length();
				// have gain ground
				if(mergedseq.length() > gainPos)
				{
					string gainStr = mergedseq.substr(gainPos);
					pacbioCorrectedStrs.back() += gainStr;
					if(round == 1)
					result.correctedLen += gainStr.length();
				}
			}
			// FMWalk failure: 
			// 1. high error 
			// 2. exceed leaves
			// 3. exceed depth
			else
			{
				if(round != 1)
				{
					// not cut off
					int startPos = seeds[i-1].first + seeds[i-1].second.length();
					int distanceBetweenSeeds = seeds[i].first + seeds[i].second.length() - seeds[i-1].first - seeds[i-1].second.length();
					pacbioCorrectedStrs.back() += corReadSeq.substr(startPos, distanceBetweenSeeds);
				}
				else if(round == 1)
				{
					// cut off
					pacbioCorrectedStrs.push_back(seeds[i].second);
					result.correctedLen += seeds[i].second.length();
				}
			}
			// output information
			if(round == 3)
			{
				result.totalSeedNum = seeds.size();
				result.totalWalkNum++;				
				if(FMWalkReturnType == 1)
				result.correctedNum++;
				else if(FMWalkReturnType == -1)
				result.highErrorNum++;
				else if(FMWalkReturnType == -2)
				result.exceedDepthNum++;
				else if(FMWalkReturnType == -3)
				result.exceedLeaveNum++;
			}
		}
		
		assert(pacbioCorrectedStrs.size() != 0);
		corReadSeq = pacbioCorrectedStrs.back();
	}
	result.merge = true;
	result.totalReadsLen = readSeq.length();
	for(size_t result_count = 0 ; result_count < pacbioCorrectedStrs.size() ; result_count++)
	result.correctedPacbioStrs.push_back(pacbioCorrectedStrs[result_count]);
	return result;
}

vector<pair<int,string> > PacBioCorrectionProcess::findSeedsUsingDynamicKmerLen(const string readSeq)
{
	vector<pair<int,string> > seeds;
	int kmerLen, iniKmerLen, minKmerLen, kmerThreshold;
	
	kmerLen = iniKmerLen = m_params.kmerLength;
	minKmerLen = m_params.minKmerLength;
	kmerThreshold = m_params.seedKmerThreshold;
	
	int readLen = readSeq.length();
	if(readLen >= iniKmerLen)
	{
		bool start = false;
		int newStartPos = -1, newStartPos2 = -1, seedStartPos = 0, walkDistance = 0;
		
		// starting search seeds
		for(int i = 0 ; i+kmerLen <= readLen ; i++)
		{
			string kmer = readSeq.substr(i, kmerLen);
			
			int kmerFreqs = 0;
			kmerFreqs += BWTAlgorithms::countSequenceOccurrencesSingleStrand(kmer, m_params.indices);
			kmerFreqs += BWTAlgorithms::countSequenceOccurrencesSingleStrand(reverseComplement(kmer), m_params.indices);
			
			walkDistance++;
			if(kmerFreqs >= kmerThreshold)	// find the seed.
			{
				// debug
				//cout << round << "-" << kmerLen << ": " << i << ", " << kmerFreqs << ", ";
				//cout << endl << kmer << endl;
				
				start = true;
				seedStartPos = i;
				
				// Until not contiguous.
				for(i++ ; i+kmerLen <= readLen ; i++)
				{
					kmer = readSeq.substr(i, kmerLen);
					kmerFreqs = 0;
					kmerFreqs += BWTAlgorithms::countSequenceOccurrencesSingleStrand(kmer, m_params.indices);
					kmerFreqs += BWTAlgorithms::countSequenceOccurrencesSingleStrand(reverseComplement(kmer), m_params.indices);
					if(kmerFreqs < kmerThreshold)
					break;
				}
				// debug
				//cout << i << ", " << kmerFreqs << ".\n";
				
				seeds.push_back(make_pair(seedStartPos, readSeq.substr(seedStartPos, kmerLen + i - seedStartPos - 1)));
				
				kmerLen = iniKmerLen;
				i = seeds.back().first + seeds.back().second.length() - 1;
				newStartPos2 = i;
				walkDistance = 0;
			}
			else if(walkDistance >= m_params.seedWalkDistance[kmerLen])	// It's too long to walk, so using dynamic kmer.
			{
				walkDistance = 0;
				if(kmerLen <= minKmerLen) // It has not searched the seed using dynamic kmer.
				{
					newStartPos = i;
					newStartPos2 = i;
					kmerLen = iniKmerLen;
				}
				else
				{
					kmerLen -= 2;
					if(start == false)	// It's too long to walk from starting point.
					i = newStartPos;
					else
					i = newStartPos2;
				}
			}
		}
	}
	
	return seeds;
}

int PacBioCorrectionProcess::doubleFMWalkForPacbio(pair<int,string> firstSeed, pair<int,string> secondSeed, int minOverlap, int needWalkLen, string* mergedseq)
{
	assert(minOverlap <= (int)firstSeed.second.length() && minOverlap <= (int)secondSeed.second.length());

	int FMWalkReturnType;
	
	SAIntervalPBHybridCTree SAITree(&firstSeed.second, minOverlap, m_params.maxOverlap, needWalkLen, m_params.maxLeaves,
	m_params.indices.pBWT, m_params.indices.pRBWT, secondSeed.second, m_params.FMWKmerThreshold);
	FMWalkReturnType = SAITree.mergeTwoReads(*mergedseq);
	
	if(FMWalkReturnType < 0)
	return FMWalkReturnType;

	assert((*mergedseq).empty() != true);
	
	string mergedseq2;
	string firstSeq = reverseComplement(firstSeed.second);
	string secondSeq = reverseComplement(secondSeed.second);
	SAIntervalPBHybridCTree SAITree2(&secondSeq, minOverlap, m_params.maxOverlap, needWalkLen, m_params.maxLeaves,
	m_params.indices.pBWT, m_params.indices.pRBWT, firstSeq, m_params.FMWKmerThreshold);
	FMWalkReturnType = SAITree2.mergeTwoReads(mergedseq2);

	if((*mergedseq).length() == mergedseq2.length())
	return FMWalkReturnType;
	else if(FMWalkReturnType > 0)
	return -4;
	else if(FMWalkReturnType < 0)
	return FMWalkReturnType;
	else
	assert(false);
}

int PacBioCorrectionProcess::solveHighError(pair<int,string> firstSeed, pair<int,string> secondSeed, int minOverlap, int needWalkLen, string* mergedseq)
{
	int FMWalkReturnType;
	int minOverlapTmp = minOverlap;
	
	do
	{
		FMWalkReturnType = doubleFMWalkForPacbio(firstSeed, secondSeed, minOverlapTmp, needWalkLen, mergedseq);
		//minOverlapTmp--;
		//minOverlapTmp-=2;
		minOverlapTmp=(minOverlapTmp*2)/3;
		
	}while(FMWalkReturnType != 1 && minOverlapTmp >= m_params.minKmerLength);
	
	return FMWalkReturnType;
}


//
//
//
PacBioCorrectionPostProcess::PacBioCorrectionPostProcess(std::ostream* pCorrectedWriter,
std::ostream* pDiscardWriter,
const PacBioCorrectionParameters params) :
m_pCorrectedWriter(pCorrectedWriter),
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
PacBioCorrectionPostProcess::~PacBioCorrectionPostProcess()
{	
	if(m_params.algorithm == PBC_SELF || m_params.algorithm == PBC_HYBRID)
	{
		std::cout << std::endl;
		std::cout << "totalReadsLen: " << m_totalReadsLen << ", ";
		std::cout << "correctedLen: " << m_correctedLen << ", ratio: " << (float)(m_correctedLen)/m_totalReadsLen << "." << std::endl;
		std::cout << "totalSeedNum: " << m_totalSeedNum << "." << std::endl;
		std::cout << "totalWalkNum: " << m_totalWalkNum << ", ";
		std::cout << "correctedNum: " << m_correctedNum << ", ratio: " << (float)(m_correctedNum*100)/m_totalWalkNum << "%." << std::endl;
		std::cout << "highErrorNum: " << m_highErrorNum << ", ratio: " << (float)(m_highErrorNum*100)/m_totalWalkNum << "%." << std::endl;
		std::cout << "exceedDepthNum: " << m_exceedDepthNum << ", ratio: " << (float)(m_exceedDepthNum*100)/m_totalWalkNum << "%." << std::endl;
		std::cout << "exceedLeaveNum: " << m_exceedLeaveNum << ", ratio: " << (float)(m_exceedLeaveNum*100)/m_totalWalkNum << "%." << std::endl;
		
		if(m_params.algorithm == PBC_SELF)
		std::cout << "disBetweenSeeds: " << m_seedDis/m_totalWalkNum << std::endl << std::endl;
	}
}


// Writting results for kmerize and validate
void PacBioCorrectionPostProcess::process(const SequenceWorkItem& item, const PacBioCorrectionResult& result)
{
	if(m_params.algorithm == PBC_SELF || m_params.algorithm == PBC_HYBRID)
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
		
		if (result.merge)
		{
			//cout << result.correctSequence.toString();
			/*SeqItem mergeRecord;
			stringstream ss;
			ss << item.read.id << "_before_len:" << result.correctSequence.toString().length();
			mergeRecord.id = ss.str();
			mergeRecord.seq = result.correctSequence;
			mergeRecord.write(*m_pCorrectedWriter);*/
			
			for(size_t i = 0 ; i < result.correctedPacbioStrs.size() ; i++)
			{
				SeqItem mergeRecord2;
				std::stringstream ss2;
				ss2 << item.read.id << "_" << i << "_" << result.correctedPacbioStrs[i].toString().length();
				mergeRecord2.id = ss2.str();
				mergeRecord2.seq = result.correctedPacbioStrs[i];
				mergeRecord2.write(*m_pCorrectedWriter);
			}
		}
	}
}

