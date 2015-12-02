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
#include "RollingPBSelfCTree.h"

#include "Timer.h"
using namespace std;


PacBioCorrectionProcess::PacBioCorrectionProcess(const PacBioCorrectionParameters params) : m_params(params)
{
}

PacBioCorrectionProcess::~PacBioCorrectionProcess()
{

}


// PacBio Self Correction by Ya and YTH, v20151202.
// 1. Identify highly-accurate seeds within PacBio reads
// 2. For each pair of seeds, perform kmer extension using local kmer frequency collected by FM-index extension
PacBioCorrectionResult PacBioCorrectionProcess::PBSelfCorrection(const SequenceWorkItem& workItem)
{	
	PacBioCorrectionResult result;
	
	std::vector<SeedFeature> seedVec, pacbioCorrectedStrs;
	std::string readSeq = workItem.read.seq.toString();
		
	// find seeds using large kmer, which can not be less than 17bp.
	seedVec = searchingSeedsUsingSolidKmer(readSeq);	
	result.totalSeedNum = seedVec.size();
	
	std::cout << workItem.read.id << "\n";

	// push the first seed into pacbioCorrectedStrs, which will be popped later as source seed
	if(seedVec.size() >= 2)
	{
		result.correctedLen += seedVec[0].seedStr.length();
		pacbioCorrectedStrs.push_back(seedVec[0]);
	}
	else
	{
		// give up reads with less than 2 seeds
		result.merge = false;
		return result;
	}
	
	// for each pair of seeds, perform kmer extension using local kmer frequency collected by LF mapping
	for(size_t targetSeed = 1 ; targetSeed < seedVec.size() ; targetSeed++)
	{
		// large kmer is used for collecting overlapping reads having common source or target seedVec
		size_t largeKmerSize = m_params.kmerLength;
		
		// small kmer is used for extension using kmer hashtable from overlapping reads
		size_t smallKmerSize = m_params.minKmerLength;
		
		// number of trials for extension to the same target seed
		size_t numOfTrials = 0;
		
		int FMWalkReturnType = 0, prevFMWalkReturnType = 0;
		SeedFeature source = pacbioCorrectedStrs.back();
		
		// Multiple targets will be tested for FM-index walk from source to target, until m_params.numOfNextTarget times.
		for(int nextTargetSeed = 0 ; nextTargetSeed < (int)m_params.numOfNextTarget && targetSeed + nextTargetSeed < seedVec.size() ; nextTargetSeed++)
		{
			// std::cout << "======= " << result.totalWalkNum << " =======\n";

			size_t currTargetIndex = targetSeed+nextTargetSeed;
			SeedFeature target =  seedVec[currTargetIndex];

			// kmer size should be adjusted according to seed sizes
			if( numOfTrials == 0 && ((source.seedLength < largeKmerSize) || (target.seedLength < largeKmerSize) ))
			{
				largeKmerSize = std::min(source.seedLength, target.seedLength);
				smallKmerSize = largeKmerSize - 2;
			}
			
			// Estimate distance between source and target, but this may over-estimate due to more insertion errors
			// Note that source seed has been updated and no long stands for the original seed, which is seedVec[targetSeed-1]
			int dis_between_src_target = target.seedStartPos - seedVec[targetSeed-1].seedStartPos - seedVec[targetSeed-1].seedStr.length();
			
			// skip seeds with large distance in between for speedup
			if( (source.seedLength < largeKmerSize) || (target.seedLength < largeKmerSize) 
				|| (dis_between_src_target >= (int)m_params.maxSeedInterval) )
			{
				// std::cout << pacbioCorrectedStrs.size()-1 << ": " 
				// << source.seedStartPos << "-" << source.seedStartPos+source.seedLength-1 <<  ":" << source.seedLength << ", " 
				// <<	target.seedStartPos << "-" << target.seedStartPos+target.seedLength-1 <<  ":" << target.seedLength 
				// << ", dis: " << dis_between_src_target << ":" << FMWalkReturnType << ".\n";
				break;
			}
			// speed up by starting from smaller kmer when distance is large, in which majority of walks failed
			else if(m_params.isFirst && dis_between_src_target >= 70 && prevFMWalkReturnType != -3)
			{
				largeKmerSize -= 2;
				smallKmerSize -= 2;
			}
		
			std::string mergedseq;
			FMWalkReturnType = extendBetweenSeeds(source, target, mergedseq, largeKmerSize, smallKmerSize, dis_between_src_target);

			if(FMWalkReturnType > 0)
			{
				// FMWalk success
				size_t extendStartPos = source.seedLength;
				std::string extendedStr = mergedseq.substr(extendStartPos);
				
				// append extended string into last corrected seed string and update related seed attributes
				pacbioCorrectedStrs.back().append(extendedStr);

				// result statistics
				result.correctedLen += extendedStr.length();
				result.correctedNum++;
				result.seedDis += dis_between_src_target;
				
				// jump to nextTargetSeed+1 if more than one target was tried and succeeded
				targetSeed = targetSeed + nextTargetSeed;
				break;
			}
			else
			{
				// return <0: give up this source seed
				// return 0: retry the same target
				// return >0: move on to next target
				int ActionFlag = FailureActions(largeKmerSize, smallKmerSize, numOfTrials, source.seedLength, target.seedLength, FMWalkReturnType, prevFMWalkReturnType);
				if(ActionFlag <0)
					break;
				else if(ActionFlag == 0)
					nextTargetSeed--;
			}
			
			prevFMWalkReturnType = FMWalkReturnType;
		}// end of next target seed
		
		// All targets failure: 
		// 0: seed inter-distance too large
		// -1: kmer extension failed at later stage 
		// -4: kmer extension failed at early stage 
		// -2: exceed depth
		// -3: exceed leaves
		if(FMWalkReturnType <= 0)
		{
			// push seedVec[targetSeed] into results, which will become new source in the next iteration
			result.seedDis += seedVec[targetSeed].seedStartPos - seedVec[targetSeed-1].seedStartPos - seedVec[targetSeed-1].seedStr.length();
			result.correctedLen += seedVec[targetSeed].seedStr.length();

			// retain uncorrected parts of reads
			if(!m_params.isSplit)
			{
				size_t startPos = seedVec[targetSeed-1].seedStartPos + seedVec[targetSeed-1].seedStr.length();
				size_t extendedLen = seedVec[targetSeed].seedStartPos + seedVec[targetSeed].seedStr.length() - startPos;

				pacbioCorrectedStrs.back().append(readSeq.substr(startPos,extendedLen));
			}
			else
			{
				// split original read into seeds and discard uncorrected parts of reads
				pacbioCorrectedStrs.push_back(seedVec[targetSeed]);
			}
			
			// statistics of FM extension
			if(FMWalkReturnType == -1 || FMWalkReturnType == -4)
				result.highErrorNum++;
			else if(FMWalkReturnType == -2)
				result.exceedDepthNum++;
			else if(FMWalkReturnType == -3)
				result.exceedLeaveNum++;
			
		}
		
		result.totalWalkNum++;
	}// end of each target seed
	
	result.merge = true;
	result.totalReadsLen = readSeq.length();
	for(size_t result_count = 0 ; result_count < pacbioCorrectedStrs.size() ; result_count++)
		result.correctedPacbioStrs.push_back(pacbioCorrectedStrs[result_count].seedStr);
	
	return result;
}

std::vector<SeedFeature> PacBioCorrectionProcess::searchingSeedsUsingSolidKmer(const std::string readSeq, size_t contaminatedCutoff)
{
	std::vector<SeedFeature> seedVec;
	const size_t smallKmerSize = m_params.kmerLength;
	const size_t largeKmerSize = m_params.isFirst?smallKmerSize:smallKmerSize+6;
	
	size_t dynamicKmerSize = largeKmerSize;
	
	size_t kmerThreshold = m_params.seedKmerThreshold;

	// prevention of short reads
	if(readSeq.length() < largeKmerSize) return seedVec;

	// search for solid kmers and group consecutive solids kmers into one seed
	// reduce kmerThreshold if no seeds can be found within 500bp
	for(size_t i = 0 ; i+dynamicKmerSize <= readSeq.length() ; i++)
	{
		std::string kmer = readSeq.substr(i, dynamicKmerSize);

		size_t fwdKmerFreqs = BWTAlgorithms::countSequenceOccurrencesSingleStrand(kmer, m_params.indices);
		size_t rvcKmerFreqs = BWTAlgorithms::countSequenceOccurrencesSingleStrand(reverseComplement(kmer), m_params.indices);
		size_t kmerFreqs = fwdKmerFreqs+rvcKmerFreqs;
		
		// std::cout << i << ": " << kmer << "\t" << kmerFreqs <<":" << fwdKmerFreqs << ":" << rvcKmerFreqs << "\n";
		
		if(kmerFreqs >= kmerThreshold && fwdKmerFreqs>=1 && rvcKmerFreqs>=1)
		{
			// skip repeat kmers which may bring noisy reads
			if(m_params.isFirst && kmerFreqs >= 40)
				continue;
			
			// skip low-complexity seeds, majority of high-GC seeds are accurate
			float GCRatio = 0;
			if(isLowComplexity(kmer, GCRatio) /*|| GCRatio > 0.7*/)
			{
				// std::cout << "Problematic Seed: " << GCRatio << "\n";
				continue;
			}
			
			size_t seedStartPos = i, 
			seedLen = 0;
			
			// Group consecutive solid kmers into one seed if possible
			size_t maxKmerFreq = kmerFreqs;
			for(i++ ; i+dynamicKmerSize <= readSeq.length(); i++)
			{
				kmer = readSeq.substr(i, dynamicKmerSize);
				fwdKmerFreqs = BWTAlgorithms::countSequenceOccurrencesSingleStrand(kmer, m_params.indices);
				rvcKmerFreqs = BWTAlgorithms::countSequenceOccurrencesSingleStrand(reverseComplement(kmer), m_params.indices);
				kmerFreqs = fwdKmerFreqs + rvcKmerFreqs;
		
				// std::cout << i << ": " << kmer << "\t" << kmerFreqs <<":" << fwdKmerFreqs << ":" << rvcKmerFreqs << "\n";

				// skip repeat kmers which may bring noisy reads
				if(m_params.isFirst && kmerFreqs >= 40)
					break;

				// contaminated seeds lead to large kmer freq	
				maxKmerFreq = std::max(maxKmerFreq,kmerFreqs);

				if(kmerFreqs >= (size_t) kmerThreshold && fwdKmerFreqs>=1 && rvcKmerFreqs>=1)
					seedLen++;
				else
					break;
			}

			size_t seedEndPos = seedStartPos+seedLen+dynamicKmerSize-1;
			
			// skip contaminated seeds for the 1st round only
			// the 2nd round may rule out seeds within large repeats
			if(maxKmerFreq > contaminatedCutoff && m_params.isFirst)
				continue;
							
			// push concatenated seeds into seed vector
			// seedVec.push_back(make_pair(seedStartPos, readSeq.substr(seedStartPos, seedEndPos-seedStartPos+1)));
			SeedFeature newSeed(seedStartPos, readSeq.substr(seedStartPos, seedEndPos-seedStartPos+1));
			seedVec.push_back(newSeed);
			
			// std::cout << seedLen+dynamicKmerSize << "\t" << readSeq.substr(seedStartPos, seedEndPos-seedStartPos+1) << "\n";

			// jump to the index after new seed
			i = i + dynamicKmerSize - 2;
			
			// reset kmer size and threshold back to strict criterion
			kmerThreshold = m_params.seedKmerThreshold;
			dynamicKmerSize = largeKmerSize;
		}// end of sufficient kmerThreshold

		// reduce kmerThreshold if no seeds can be found within 500bp
		size_t prevSeedEndPos = seedVec.empty()? 0 : seedVec.back().seedStartPos + seedVec.back().seedStr.length();
		int distToPrevSeed = i + 1 - prevSeedEndPos;
		
		/*  It happened when kmer become larger in 2nd round
			The regions with ultra-high errors are difficult to contain seeds and were mostly not corrected in the first round
			reduce kmer size and threshold if no seeds can be found within maxSeedInterval*/
		if(!m_params.isFirst && distToPrevSeed >= (int)m_params.maxSeedInterval)
		{
			if(kmerThreshold > 10)
			{
				kmerThreshold = 10;
				i = prevSeedEndPos;
			}
			
			if( dynamicKmerSize > 17)
			{
				dynamicKmerSize = 17;
				i = prevSeedEndPos;
			}
		}
			
	}// end of for
	
	return seedVec;
}

// Perform FMindex extension between source and target seeds
// Return FMWalkReturnType
int PacBioCorrectionProcess::extendBetweenSeeds(SeedFeature& source, SeedFeature& target, std::string& mergedseq, 
											size_t largeKmerSize, size_t smallKmerSize, size_t dis_between_src_target)
{
	SAIPBSelfCorrectTree SAITree(m_params.indices.pBWT, m_params.indices.pRBWT, m_params.FMWKmerThreshold);

	// Collect local kmer frequency from source upto srcMaxLength
	const int srcMaxLength = 1.2*(dis_between_src_target+20) + source.seedLength + smallKmerSize;
	SAITree.addHashBySingleSeed(source.seedStr, largeKmerSize, smallKmerSize, srcMaxLength, m_params.isFirst);

	// Collect local kmer frequency from target upto targetMaxLength
	std::string rvcTargetStr = reverseComplement(target.seedStr);
	const int targetMaxLength = 1.2*(dis_between_src_target+20) + rvcTargetStr.length() + smallKmerSize;
	size_t expectedLength = dis_between_src_target + rvcTargetStr.length();
	
	SAITree.addHashBySingleSeed(rvcTargetStr, largeKmerSize, smallKmerSize, targetMaxLength, m_params.isFirst, expectedLength);


	// Estimate upper/lower/expected bounds of search depth
	int srcMinLength = 0.8*(dis_between_src_target-20) + source.seedLength + smallKmerSize;
	if(srcMinLength < 0) srcMinLength = 0;
	expectedLength = source.seedLength + dis_between_src_target + smallKmerSize;
	
	int FMWalkReturnType = SAITree.mergeTwoSeedsUsingHash(source.seedStr, target.seedStr, mergedseq, smallKmerSize, m_params.maxLeaves,
													  srcMinLength, srcMaxLength, expectedLength);
	
	// std::cout << source.seedStartPos << "-" << source.seedStartPos+source.seedLength-1 <<  ":" << source.seedLength << ", " 
	// <<	target.seedStartPos << "-" << target.seedStartPos+target.seedLength-1 <<  ":" << target.seedLength 
	// << ", dis: " << dis_between_src_target << ", " << expectedLength << ", " << srcMaxLength << ", "
	// << FMWalkReturnType << ".\n";
	
	return FMWalkReturnType;
}
// refine seed interval using larger kmer
std::pair<size_t, size_t> PacBioCorrectionProcess::refineRepeatSeed(const std::string readSeq, size_t repeatKmerSize, size_t seedStartPos, size_t seedEndPos)
{
	size_t newSeedStartPos = (size_t)-1;
	size_t newSeedEndPos = (size_t)-1;
	
	for(size_t i=seedStartPos ; i+repeatKmerSize-1 <= seedEndPos; i++)
	{
		std::string kmer = readSeq.substr(i, repeatKmerSize);
		size_t fwdKmerFreqs = BWTAlgorithms::countSequenceOccurrencesSingleStrand(kmer, m_params.indices);
		size_t rvcKmerFreqs = BWTAlgorithms::countSequenceOccurrencesSingleStrand(reverseComplement(kmer), m_params.indices);
		size_t kmerFreqs = fwdKmerFreqs + rvcKmerFreqs;

		std::cout << i << ": " << kmer << "\t" << kmerFreqs <<":" << fwdKmerFreqs << ":" << rvcKmerFreqs << "\n";

		if(kmerFreqs >= (size_t) m_params.seedKmerThreshold && fwdKmerFreqs>=1 && rvcKmerFreqs>=1)
		{
			// only capture the first seed
			if(newSeedStartPos == (size_t)-1)
				newSeedStartPos = i;
				
			if(newSeedStartPos != (size_t)-1)
				newSeedEndPos = i;
		}
		// out of seeds
		else if(newSeedStartPos != (size_t)-1)
			break;
	}
	
	std::cout << newSeedStartPos << "\t" << newSeedEndPos << "\n";
	return std::make_pair(newSeedStartPos, newSeedEndPos+repeatKmerSize);
}

// return <0: give up and break
// return 0: retry the same target
// return >0: continue to next target
int PacBioCorrectionProcess::FailureActions(size_t& largeKmerSize, size_t& smallKmerSize, size_t& numOfTrials, 
								size_t sourceStrLen, size_t targetStrLen, int FMWalkReturnType, int prevFMWalkReturnType)
{
	// extension failed due to insufficient kmers, reduce large and small kmer sizes
	if(FMWalkReturnType==-1 || FMWalkReturnType==-4)
	{
		// kmers have been enlarged due to repeats, shrink will lead to infinite loop
		if(prevFMWalkReturnType==-3 )
			return -1;
		
		largeKmerSize -= 2;
		smallKmerSize -= 2;

		// don't aggressively reduce kmer in the 1st found where most kmers are errors
		if(m_params.isFirst && (largeKmerSize < 15 || smallKmerSize < 9) )
			return -1;

			
		return 0;
	}
	
	// increase smallKmerSize for reducing repeats
	else if(FMWalkReturnType==-3)
	{
		if(prevFMWalkReturnType==-4 || prevFMWalkReturnType==-1)
			return -1;

		// exponential growth is required in super large repeats. Otherwise speed is too slow
		smallKmerSize += pow(2, numOfTrials+1);
		largeKmerSize += pow(2, numOfTrials+1);

		if(sourceStrLen < largeKmerSize || targetStrLen < largeKmerSize)
			return -1;
		
		numOfTrials++;
		return 0;
	}
	else if(FMWalkReturnType==-2)
	{
		// probably chimera, need more observations
		// largeKmerSize = m_params.kmerLength;
		// smallKmerSize = largeKmerSize - 2;
		numOfTrials++;
		return 1;
	}
	
	return 1;
}

bool PacBioCorrectionProcess::isLowComplexity (std::string seq, float & GCratio, float threshold)
{
	size_t seqLen = seq.length();
	size_t countG =0 ;
	size_t countC =0 ;
	size_t countT =0 ;
	size_t countA =0 ;

	for (size_t i=0; i<seqLen; i++)
	{
		switch(seq[i]){
			case 'A': countA ++ ;break;
			case 'T': countT ++ ;break;
			case 'C': countC ++ ;break;
			case 'G': countG ++ ;break;
			default:  assert(false);
		}
	}

	GCratio = (float)(countG+countC)/seqLen ;

	if (  ((float) countA/seqLen >= threshold ) || ((float) countT/seqLen >= threshold)
			|| ((float) countC/seqLen >= threshold ) || ((float) countG/seqLen >= threshold) )
		return true;

	return false;

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
PacBioCorrectionPostProcess::~PacBioCorrectionPostProcess()
{	
	if(m_params.algorithm == PBC_SELF || m_params.algorithm == PBC_HYBRID)
	{
		if(m_totalWalkNum>0 && m_totalReadsLen>0)
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
}


// Writting results for kmerize and validate
void PacBioCorrectionPostProcess::process(const SequenceWorkItem& item, const PacBioCorrectionResult& result)
{
	if(m_params.algorithm == PBC_SELF || m_params.algorithm == PBC_HYBRID)
	{		
		if (result.merge)
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
		else
		{
			// write into discard.fa
			SeqItem mergeRecord2;
			mergeRecord2.id = item.read.id;
			mergeRecord2.seq = item.read.seq;
			mergeRecord2.write(*m_pDiscardWriter);
		}
	}
}

