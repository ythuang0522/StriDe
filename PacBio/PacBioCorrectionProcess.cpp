///-----------------------------------------------
// Copyright 2015 National Chung Cheng University
// Written by Yao-Ting Huang
// Released under the GPL
//-----------------------------------------------
//
// PacBioCorrectionProcess.cpp - Self-correction using FM-index walk for PacBio reads
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
	// std::cout << workItem.read.id << "\n";

	PacBioCorrectionResult result;
	
	std::vector<SeedFeature> seedVec, pacbioCorrectedStrs;
	std::string readSeq = workItem.read.seq.toString();
		
	// find seeds using fixed or dynamic kmers depending on 1st round or not
	if(m_params.isFirst)
		// seeding by fixed kmer size suitable for the 1st round of error correction where most regions are errors
		// error seeds within repeat regions are ruled out in particular
		seedVec = seedingByFixedKmer(readSeq);	
	else
		// find seeds by dynamic kmers, which is suitable for the 2nd round of error correction
		// where repeat regions require large kmers and error-prone regions require small kmers
		seedVec = seedingByDynamicKmer(readSeq);
		
	result.totalSeedNum = seedVec.size();

	// push the first seed into pacbioCorrectedStrs, which will be popped later as source seed
	if(seedVec.size() >= 2)
	{
		result.correctedLen += seedVec.at(0).seedStr.length();
		// if(m_params.isSplit)
			pacbioCorrectedStrs.push_back(seedVec.at(0));
		// else
			// pacbioCorrectedStrs.push_back(readSeq.substr(0, seedVec.at(0).seedEndPos+1));
		
		// if(!m_params.isSplit)
			// pacbioCorrectedStrs.back().seedStr.reserve(readSeq.length()*1.5);
	}
	else
	{
		// give up reads with less than 2 seeds
		result.merge = false;
		return result;
	}
	
	if(m_params.isFirst)
	{
		// reserve sufficient str length for fast append
		pacbioCorrectedStrs.back().seedStr.reserve(readSeq.length());
		initCorrect(readSeq, seedVec, pacbioCorrectedStrs, result);
	}
	else
		realCorrect(readSeq, seedVec, pacbioCorrectedStrs, result);
		
	
	result.merge = true;
	result.totalReadsLen = readSeq.length();
	for(size_t result_count = 0 ; result_count < pacbioCorrectedStrs.size() ; result_count++)
		result.correctedPacbioStrs.push_back(pacbioCorrectedStrs[result_count].seedStr);
	
	return result;
}

void PacBioCorrectionProcess::initCorrect(std::string& readSeq, std::vector<SeedFeature>& seedVec, std::vector<SeedFeature>& pacbioCorrectedStrs, PacBioCorrectionResult& result)
{
	// for each pair of seeds, perform kmer extension using local kmer frequency
	for(size_t targetSeed = 1 ; targetSeed < seedVec.size() ; targetSeed++)
	{				
		// number of trials of extension to the same target seed
		size_t numOfTrials = 0;
		
		int FMWalkReturnType = 0, prevFMWalkReturnType = 0;
		
		// source is increasing because no split in the 1st round, this is slow, better replace with pointer
		SeedFeature source = pacbioCorrectedStrs.back();
		SeedFeature target =  seedVec.at(targetSeed);

		// extension kmer is used for extension using local kmer hashtable collected from overlapping reads
		// default: smaller beset kmer size from both seeds -2
		size_t extendKmerSize = std::min(source.endBestKmerSize, seedVec.at(targetSeed).startBestKmerSize) - 2;
					
		// Multiple targets will be tested for FM-index walk from source to target, if target is error seed, until m_params.numOfNextTarget times.
		for(int nextTargetSeed = 0 ; nextTargetSeed < (int)m_params.numOfNextTarget && targetSeed + nextTargetSeed < seedVec.size() ; nextTargetSeed++)
		{
			// std::cout << "======= " << result.totalWalkNum << " =======\n";
			
			// extendKmerSize should be large for repeat seeds
			if((source.isRepeat || target.isRepeat) )
			{
				extendKmerSize = std::min((int)source.seedLength, (int)target.seedLength);
				if(extendKmerSize > m_params.kmerLength+2) 
						extendKmerSize = m_params.kmerLength+2;
			}

			// Estimate distance between source and target, but this may over-estimate due to more insertion errors
			// Note that source seed has been updated and no long stands for the original seed, which is seedVec[targetSeed-1]
			int dis_between_src_target = target.seedStartPos - seedVec.at(targetSeed-1).seedStartPos - seedVec.at(targetSeed-1).seedStr.length();
		
			// skip seeds with large distance in between for speedup
			if(dis_between_src_target >= (int)m_params.maxSeedInterval) 
				break;

			// PB159615_16774.fa require smaller kmer
			// if(dis_between_src_target >= 50 && !source.isRepeat && !target.isRepeat && extendKmerSize>11)
				// extendKmerSize -= 2;
		
			// PB36993_4517.fa, 185
			// Only one path is found but is wrong 
			// both seeds are correct but are tandem repeats which appear more than once in different copies
			// Skip extension if both seeds are repeat and dist is large and kmer freq is large
			if((source.isRepeat && target.isRepeat) && dis_between_src_target>=70 && (source.endKmerFreq>40 || target.startKmerFreq>40))
				break;

			// extension using local kmer hashtable collected from overlapping reads
			std::string mergedseq;
			FMWalkReturnType = extendBetweenSeeds(source, target, readSeq, mergedseq, extendKmerSize, dis_between_src_target);

			if(FMWalkReturnType > 0)
			{
				// FMWalk success
				// size_t extendStartPos = source.seedLength;
				// std::string extendedStr = mergedseq.substr(extendStartPos);
				
				// append extended string into last corrected seed string and update related seed attributes
				// pacbioCorrectedStrs.back().append(extendedStr);
				pacbioCorrectedStrs.back().append(mergedseq);
				// the last seed will become new source and should be updated
				pacbioCorrectedStrs.back().endBestKmerSize = target.endBestKmerSize;
				pacbioCorrectedStrs.back().isRepeat = target.isRepeat;

				// result statistics
				// result.correctedLen += extendedStr.length();
				result.correctedLen += mergedseq.length();
				result.correctedNum++;
				result.seedDis += dis_between_src_target;
				
				// jump to nextTargetSeed+1 if more than one target was tried and succeeded
				targetSeed = targetSeed + nextTargetSeed;
				break;
			}
			else
			{
				// return <0: give up this source seed
				int ActionFlag = FMWalkFailedActions(extendKmerSize, numOfTrials, source, target, FMWalkReturnType, prevFMWalkReturnType);
				if(ActionFlag <0)
					break;
				// return 0: retry the same target
				else if(ActionFlag == 0)
					nextTargetSeed--;
				// return >0: move on to next target
				else
				{
					// target =  targetSeed+nextTargetSeed+1<seedVec.size()?seedVec[targetSeed+nextTargetSeed+1]:target;
					// extendKmerSize = std::min(source.endBestKmerSize, target.startBestKmerSize) - 2;
				}
			}
			
			prevFMWalkReturnType = FMWalkReturnType;
		}// end of next target seed
		
		// All targets failure: 
		// 0: seed inter-distance too large
		// -1: kmer extension failed at later stage close to target
		// -4: kmer extension failed at early stage close to source
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
				pacbioCorrectedStrs.back().endBestKmerSize = seedVec.at(targetSeed).endBestKmerSize;
				pacbioCorrectedStrs.back().isRepeat = seedVec.at(targetSeed).isRepeat;
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
}

void PacBioCorrectionProcess::realCorrect(std::string& readSeq, std::vector<SeedFeature>& seedVec, std::vector<SeedFeature>& pacbioCorrectedStrs, PacBioCorrectionResult& result)
{
	// for each pair of seeds, perform kmer extension using local kmer frequency
	for(size_t targetSeed = 1 ; targetSeed < seedVec.size() ; targetSeed++)
	{				
		// number of trials of extension to the same target seed
		size_t numOfTrials = 0;
		
		int FMWalkReturnType = 0, prevFMWalkReturnType = 0;
		SeedFeature source = pacbioCorrectedStrs.back();
		
		//PB36993_4517, GTAATAAGGAATATCTCAATTTT is a repeat seed with large frequency leading to -2
		// re-estimate the best kmer size if source becomes larger
		if(source.endKmerFreq>90)
			source.estimateBestKmerSize(m_params.indices.pBWT);
		
		SeedFeature target =  seedVec.at(targetSeed);

		// small kmer is used for extension using kmer hashtable from overlapping reads
		size_t extendKmerSize = std::min(source.endBestKmerSize, seedVec.at(targetSeed).startBestKmerSize) - 2;

		// Multiple targets will be tested for FM-index walk from source to target, if target is error seed, until m_params.numOfNextTarget times.
		for(int nextTargetSeed = 0 ; nextTargetSeed < (int)m_params.numOfNextTarget && targetSeed + nextTargetSeed < seedVec.size() ; nextTargetSeed++)
		{
			// std::cout << "======= " << result.totalWalkNum << " =======\n";
			
			// Estimate distance between source and target, but this may over-estimate due to more insertion errors
			// Note that source seed has been updated and no long stands for the original seed, which is seedVec[targetSeed-1]
			int dis_between_src_target = target.seedStartPos - seedVec.at(targetSeed-1).seedStartPos - seedVec.at(targetSeed-1).seedStr.length();
					
			// skip seeds with large distance in between for speedup
			if(dis_between_src_target >= (int)m_params.maxSeedInterval) 
				break;

			// PB159615_16774.fa require smaller kmer
			// else if(m_params.isFirst && dis_between_src_target >= 50 && !source.isRepeat && !target.isRepeat && extendKmerSize>11)
				// extendKmerSize -= 2;
		

			std::string mergedseq;
			FMWalkReturnType = extendBetweenSeeds(source, target, readSeq, mergedseq, extendKmerSize, dis_between_src_target);

			if(FMWalkReturnType > 0)
			{
				// FMWalk success
				// size_t extendStartPos = source.seedLength;
				// std::string extendedStr = mergedseq.substr(extendStartPos);
				
				// append extended string into last corrected seed string and update related seed attributes
				// pacbioCorrectedStrs.back().append(extendedStr);
				pacbioCorrectedStrs.back().append(mergedseq);
				
				// the last seed will become new source and should be updated
				pacbioCorrectedStrs.back().endBestKmerSize = target.endBestKmerSize;
				pacbioCorrectedStrs.back().isRepeat = target.isRepeat;

				// result statistics
				// result.correctedLen += extendedStr.length();
				result.correctedLen += mergedseq.length();
				result.correctedNum++;
				result.seedDis += dis_between_src_target;
				
				// jump to nextTargetSeed+1 if more than one target was tried and succeeded
				targetSeed = targetSeed + nextTargetSeed;
				break;
			}
			else
			{
				// return <0: give up this source seed
				int ActionFlag = FMWalkFailedActions(extendKmerSize, numOfTrials, source, target, FMWalkReturnType, prevFMWalkReturnType);
				if(ActionFlag <0)
					break;
				// return 0: retry the same target
				else if(ActionFlag == 0)
					nextTargetSeed--;
				// return >0: move on to next target
				else
				{
					// target =  targetSeed+nextTargetSeed+1<seedVec.size()?seedVec[targetSeed+nextTargetSeed+1]:target;
					// extendKmerSize = std::min(source.endBestKmerSize, target.startBestKmerSize) - 2;
				}
			}
			
			prevFMWalkReturnType = FMWalkReturnType;
		}// end of next target seed
		
		// All targets failure: 
		// 0: seed inter-distance too large
		// -1: kmer extension failed at later stage close to target
		// -4: kmer extension failed at early stage close to source
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
				pacbioCorrectedStrs.back().endBestKmerSize = seedVec.at(targetSeed).endBestKmerSize;
				pacbioCorrectedStrs.back().isRepeat = seedVec.at(targetSeed).isRepeat;
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
}

// seeding by fixed kmer size which is suitable for the 1st round of error correction where most regions are errors
// error seeds within repeat regions are ruled out in particular
std::vector<SeedFeature> PacBioCorrectionProcess::seedingByFixedKmer(const std::string readSeq, size_t contaminatedCutoff)
{
	std::vector<SeedFeature> seedVec;
	const size_t kmerSize = m_params.kmerLength;
	size_t kmerThreshold = m_params.seedKmerThreshold;

	// prevention of short reads
	if(readSeq.length() < kmerSize) return seedVec;

	// search for solid kmers and group consecutive solids kmers into one seed
	// but remove seeds within repeats
	for(size_t i = 0 ; i+kmerSize <= readSeq.length() ; i++)
	{
		std::string kmer = readSeq.substr(i, kmerSize);

		size_t fwdKmerFreqs = BWTAlgorithms::countSequenceOccurrencesSingleStrand(kmer, m_params.indices);
		size_t rvcKmerFreqs = BWTAlgorithms::countSequenceOccurrencesSingleStrand(reverseComplement(kmer), m_params.indices);
		size_t kmerFreqs = fwdKmerFreqs+rvcKmerFreqs;
		
		// std::cout << i << ": " << kmer << "\t" << kmerFreqs <<":" << fwdKmerFreqs << ":" << rvcKmerFreqs << "\n";
		
		if(kmerFreqs >= kmerThreshold && fwdKmerFreqs>=1 && rvcKmerFreqs>=1)
		{	
			// skip low-complexity seeds, e.g., AAAAAAAATAAAA, which are often error seeds
			float GCRatio = 0;
			if(isLowComplexity(kmer, GCRatio))
			{
				bool isPrevSeedCloseToRepeat = !seedVec.empty() && !seedVec.back().isRepeat 
										// remove previous seed if it's too close to the current repeat within kmerSize
										&& i - seedVec.back().seedEndPos < kmerSize 
										// ensure that the kmer frequency difference is also large
										&& seedVec.back().seedLength-kmerSize<=3; 
				
				if(isPrevSeedCloseToRepeat)	
					seedVec.pop_back();
				continue;
			}
			size_t seedStartPos = i, 
			seedLen = 0;
			
			// Group consecutive solid kmers into one seed if possible
			size_t maxKmerFreq = kmerFreqs;
			for(i++ ; i+kmerSize <= readSeq.length(); i++)
			{
				kmer = readSeq.substr(i, kmerSize);
				fwdKmerFreqs = BWTAlgorithms::countSequenceOccurrencesSingleStrand(kmer, m_params.indices);
				rvcKmerFreqs = BWTAlgorithms::countSequenceOccurrencesSingleStrand(reverseComplement(kmer), m_params.indices);
				kmerFreqs = fwdKmerFreqs + rvcKmerFreqs;
		
				// std::cout << i << ": " << kmer << "\t" << kmerFreqs <<":" << fwdKmerFreqs << ":" << rvcKmerFreqs << "\n";

				maxKmerFreq = std::max(maxKmerFreq,kmerFreqs);

				if(isLowComplexity(kmer, GCRatio)) break;
				
				if(kmerFreqs >= (size_t) kmerThreshold && fwdKmerFreqs>=1 && rvcKmerFreqs>=1)
					seedLen++;
				else
					break;
			}

			size_t seedEndPos = seedStartPos + seedLen + kmerSize-1;

			// skip contaminated seeds
			if(maxKmerFreq > contaminatedCutoff)
				continue;

			// this is a repeat seed
			if(maxKmerFreq > 40)	//40 should be determined using kmerDistribution
			{
				// For seeds within repeats, error seeds will still exceed the cutoff, e.g., 12 11 15 60 65 70 20 19..
				// refine the exact boundary by finding the segments with highest kmer frequency 60 65 70
				std::pair<size_t, size_t> kmerFreqPair = refineRepeatSeed(readSeq, seedStartPos, seedEndPos);

				// PB135123_7431.fa, error seeds close to repeats are also with high frequency
				bool isPrevSeedCloseToRepeat = !seedVec.empty() && !seedVec.back().isRepeat 
										// remove previous seed if it's too close to the current repeat within kmerSize
										&& seedStartPos - seedVec.back().seedEndPos < kmerSize 
										// ensure that the kmer frequency difference is also large
										&& (std::abs((int)seedVec.back().endKmerFreq-(int)kmerFreqPair.first)>40); 
				
				if(isPrevSeedCloseToRepeat)	
					seedVec.pop_back();

				// PB135123_7431.fa TGTAATCAGGCTGAAAA
				bool isPrevSeedBetweenRepeat = seedVec.size()>=2 && !seedVec.back().isRepeat // previous seed is not repeat
										// but previous previous seed and current seed are both repeats
										&& seedVec.at(seedVec.size()-2).isRepeat
										// and the distance between two repeat seeds are not large
										&& seedStartPos - seedVec.at(seedVec.size()-2).seedEndPos < 80; 
										
				//PB135123_7431.fa: CGGCTTTCTTTAATGAT
				bool isPrevSeedWithinLargeRepeat = seedVec.size()>=3 && !seedVec.back().isRepeat // previous seed is not repeat
										// but both previous two seeds are repeats
										&& seedVec.at(seedVec.size()-2).isRepeat && seedVec.at(seedVec.size()-3).isRepeat
										// and the distance between two repeat seeds are not large
										&& seedStartPos - seedVec.at(seedVec.size()-2).seedEndPos < 200; 


				if(isPrevSeedBetweenRepeat || isPrevSeedWithinLargeRepeat)	
					seedVec.pop_back();

				// PB135123_7431.fa: AAATTTGAAGAGACTCA, CCAGGGTATCTAAATCCTGTTT
				bool isPrevTwoSeedWithinLargeRepeat = seedVec.size()>=4 && !seedVec.back().isRepeat && !seedVec.at(seedVec.size()-2).isRepeat // previous two seed are not repeat
										// but previous seed is repeats
										&& seedVec.at(seedVec.size()-3).isRepeat 
										// and the distance between two repeat seeds are not large
										&& seedStartPos - seedVec.at(seedVec.size()-3).seedEndPos < 200 
										&& (seedVec.back().seedLength-kmerSize<=3 || seedVec.at(seedVec.size()-2).seedLength-kmerSize<=3);

				if(isPrevTwoSeedWithinLargeRepeat){
					// std::cout << "hahaha" << seedVec.back().isSmall() << "\t" << seedVec.back().seedStr << "\n";
				
					seedVec.pop_back();
					seedVec.pop_back();
				}
										
				// if(isPrevSeedCloseToRepeat || isPrevSeedBetweenRepeat || isPrevSeedWithinLargeRepeat)	
					// seedVec.pop_back();
				
				/** Unsolved cases
				PB38933_7552.fa
				1414: GTTCAGCGGAAATTTTC 1:1:0
				1415: TTCAGCGGAAATTTTCC 20:11:9
				1416: TCAGCGGAAATTTTCCA 1:1:0
				Seed:   17      TTCAGCGGAAATTTTCC
				*/
				
				SeedFeature newSeed(seedStartPos, readSeq.substr(seedStartPos, seedEndPos-seedStartPos+1), true, kmerSize, kmerThreshold*4);
				newSeed.estimateBestKmerSize(m_params.indices.pBWT);
				seedVec.push_back(newSeed);				
			}
			else 
			{
				bool isCloseToPrevRepeatSeed = !seedVec.empty() && seedVec.back().isRepeat 
												&& (seedStartPos - seedVec.back().seedEndPos <= kmerSize);
									
									
				// most error seeds are not too large or with lower frequency
				// PB135123_7431.fa: AAAACTTCGCAGTGAAC is not error but discarded
				// && seedEndPos+1-seedStartPos-kmerSize < 7;
												
				// if(seedVec.empty() || !seedVec.back().isRepeat || (seedStartPos - seedVec.back().seedEndPos > kmerSize) )
				if(seedVec.empty() || !isCloseToPrevRepeatSeed)
				{
					// push concatenated seeds into seed vector
					SeedFeature newSeed(seedStartPos, readSeq.substr(seedStartPos, seedEndPos-seedStartPos+1), false, kmerSize, kmerThreshold*4);
					// newSeed.setBestKmerSize(kmerSize);
					newSeed.estimateBestKmerSize(m_params.indices.pBWT);
					seedVec.push_back(newSeed);
				}
			
				// debug:: don't push if previous seed is repeat and too close
				// if(!seedVec.empty() && seedVec.back().isRepeat && (seedStartPos - seedVec.back().seedEndPos <= kmerSize) )
					// std::cout << "hehe\t" << readSeq.substr(seedStartPos, seedEndPos-seedStartPos+1) << "\n";
			}

			// std::cout << "Seed:\t"<< seedVec.back().seedLength << "\t" << seedVec.back().seedStr << "\n";
			
			// restart from end of previous repeat seed because multiple repeat seeds may be within the same region
			i=seedEndPos;

		}// end of sufficient kmerThreshold
			
	}// end of for
	
	return seedVec;
}

// find seeds by dynamic kmers, which is suitable for the 2nd round of error correction
// where repeat regions require large kmers and error-prone regions require small kmers
std::vector<SeedFeature> PacBioCorrectionProcess::seedingByDynamicKmer(const std::string readSeq)
{
	std::vector<SeedFeature> seedVec;
	const size_t smallKmerSize = m_params.kmerLength;
	const size_t largeKmerSize = smallKmerSize+6;
	
	// set dynamic kmer as largest kmer initially, which will reduce size when no seeds can be found within a maximum interval
	size_t dynamicKmerSize = largeKmerSize;
	size_t kmerThreshold = m_params.seedKmerThreshold;

	// prevention of short reads
	if(readSeq.length() < dynamicKmerSize) return seedVec;

	// search for solid kmers and group consecutive solids kmers into one seed
	// reduce kmer size and kmerThreshold if no seeds can be found within maxSeedInterval
	for(size_t i = 0 ; i+dynamicKmerSize <= readSeq.length() ; i++)
	{
		std::string kmer = readSeq.substr(i, dynamicKmerSize);

		size_t fwdKmerFreqs = BWTAlgorithms::countSequenceOccurrencesSingleStrand(kmer, m_params.indices);
		size_t rvcKmerFreqs = BWTAlgorithms::countSequenceOccurrencesSingleStrand(reverseComplement(kmer), m_params.indices);
		size_t kmerFreqs = fwdKmerFreqs+rvcKmerFreqs;
		
		// std::cout << i << ": " << kmer << "\t" << kmerFreqs <<":" << fwdKmerFreqs << ":" << rvcKmerFreqs << "\n";
		
		if(kmerFreqs >= kmerThreshold && fwdKmerFreqs>=1 && rvcKmerFreqs>=1)
		{	
			// skip low-complexity seeds, e.g., AAAAAAAATAAAA
			float GCRatio = 0;
			if(isLowComplexity(kmer, GCRatio)) continue;
			
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

				if(isLowComplexity(kmer, GCRatio)) break;

				maxKmerFreq = std::max(maxKmerFreq,kmerFreqs);

				if(kmerFreqs >= (size_t) kmerThreshold && fwdKmerFreqs>=1 && rvcKmerFreqs>=1)
					seedLen++;
				else
					break;
					
				
			}

			size_t seedEndPos = seedStartPos+seedLen+dynamicKmerSize-1;
																	
			// refine repeat seed only if it's small size
			if(dynamicKmerSize < largeKmerSize && maxKmerFreq > 40)
			{
				refineRepeatSeed(readSeq, seedStartPos, seedEndPos);
			}
			
			if(maxKmerFreq > 80)
			{
				// push concatenated seeds into seed vector
				SeedFeature newSeed(seedStartPos, readSeq.substr(seedStartPos, seedEndPos-seedStartPos+1), true, dynamicKmerSize, kmerThreshold+10);
				// newSeed.setBestKmerSize(dynamicKmerSize);
				newSeed.estimateBestKmerSize(m_params.indices.pBWT);
				seedVec.push_back(newSeed);
			}
			else
			{
				// push concatenated seeds into seed vector
				SeedFeature newSeed(seedStartPos, readSeq.substr(seedStartPos, seedEndPos-seedStartPos+1), false, dynamicKmerSize, kmerThreshold+10);
				// newSeed.setBestKmerSize(dynamicKmerSize);
				newSeed.estimateBestKmerSize(m_params.indices.pBWT);
				seedVec.push_back(newSeed);
			}
			
			// std::cout << "Seed:\t " << seedVec.back().seedLength << "\t" << seedVec.back().seedStr << "\n";

			// jump to the index after new seed
			i = i + dynamicKmerSize - 2;
			
			// reset kmer size and threshold back to strict criterion
			kmerThreshold = m_params.seedKmerThreshold;
			dynamicKmerSize = largeKmerSize;
		}// end of sufficient kmerThreshold

		/* reduce kmerThreshold or kmer size if no seeds can be found within maxSeedInterval
		It happened when kmer become larger in 2nd round
		The regions with ultra-high errors are difficult to contain seeds and were mostly not corrected in the first round
		reduce kmer size and threshold if no seeds can be found within maxSeedInterval*/
		size_t prevSeedEndPos = seedVec.empty()? 0 : seedVec.back().seedStartPos + seedVec.back().seedStr.length();
		int distToPrevSeed = i + 1 - prevSeedEndPos;
		
		if(distToPrevSeed >= (int)m_params.maxSeedInterval)
		{
			if(kmerThreshold > 10)
			{
				kmerThreshold = 10;
				i = prevSeedEndPos;
			}
			
			// PB36993_4517, CATTTCAAAATCACCATA is a wrong seed within large tandem repeats
			// previous seed is a correct repeat seed
			if( dynamicKmerSize > 17 && !seedVec.empty() && !seedVec.back().isRepeat)
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
int PacBioCorrectionProcess::extendBetweenSeeds(SeedFeature& source, SeedFeature& target, std::string& rawSeq, std::string& mergedseq, 
												size_t extendKmerSize, size_t dis_between_src_target)
{
	const double maxRatio = 1.1;
	const double minRatio = 0.9;
	const int minOffSet = 30;	//PB159615_16774.fa contains large indels > 30bp

	size_t srcKmerSize = std::max(source.endBestKmerSize, extendKmerSize);

	std::string rawSubseq = rawSeq.substr(target.seedStartPos+1-dis_between_src_target-srcKmerSize, target.seedEndPos-source.seedStartPos+1);	
	SAIPBSelfCorrectTree SAITree(m_params.indices.pBWT, m_params.indices.pRBWT, rawSubseq, m_params.FMWKmerThreshold);
		
	// this occurs when one end is repeat requiring larger extendKmerSize while the other requires small extendKmerSize
	// the source str should be the longest one
	// const int srcMaxLength = maxRatio*(dis_between_src_target+minOffSet) + source.seedLength + extendKmerSize;
	// size_t sourceFreq = SAITree.addHashBySingleSeed(source.seedStr, source.endBestKmerSize, extendKmerSize, srcMaxLength, m_params.isFirst);
	
	//size_t srcKmerSize = std::max(source.endBestKmerSize, extendKmerSize);
	std::string srcStr = source.seedStr.substr(source.seedStr.length()-srcKmerSize);
	const int srcMaxLength = maxRatio*(dis_between_src_target+minOffSet) + srcStr.length() + extendKmerSize;
	size_t sourceFreq = SAITree.addHashBySingleSeed(srcStr, source.endBestKmerSize, extendKmerSize, srcMaxLength, m_params.isFirst);
	
	// Collect local kmer frequency from target upto targetMaxLength
	std::string rvcTargetStr = reverseComplement(target.seedStr);
	const int targetMaxLength = maxRatio*(dis_between_src_target+minOffSet) + rvcTargetStr.length() + extendKmerSize;
	size_t expectedLength = dis_between_src_target + rvcTargetStr.length();
	assert(rvcTargetStr.length()>=extendKmerSize);
	size_t targetFreq = SAITree.addHashBySingleSeed(rvcTargetStr, target.startBestKmerSize, extendKmerSize, targetMaxLength, m_params.isFirst, expectedLength);

	// Estimate upper/lower/expected bounds of search depth
	// int srcMinLength = minRatio*(dis_between_src_target-minOffSet) + source.seedLength + extendKmerSize;
	// if(srcMinLength < 0) srcMinLength = 0;
	// expectedLength = source.seedLength + dis_between_src_target + target.seedLength;	
	// int FMWalkReturnType = SAITree.mergeTwoSeedsUsingHash(source.seedStr, target.seedStr, mergedseq, extendKmerSize, m_params.maxLeaves,
													  // srcMinLength, srcMaxLength, expectedLength);
	
	int srcMinLength = minRatio*(dis_between_src_target-minOffSet) + srcStr.length() + extendKmerSize;
	if(srcMinLength < 0) srcMinLength = 0;
	expectedLength = srcStr.length() + dis_between_src_target + target.seedLength;
	int FMWalkReturnType = SAITree.mergeTwoSeedsUsingHash(srcStr, target.seedStr, mergedseq, extendKmerSize, m_params.maxLeaves,
													  srcMinLength, srcMaxLength, expectedLength);
	
	// retain only extended portion
	if(!mergedseq.empty())
		mergedseq = mergedseq.substr(srcStr.length());
	
	// std::cout << source.seedStartPos << "-" << source.seedStartPos+source.seedLength-1 <<  ":" << source.seedLength << ", " 
	// <<	target.seedStartPos << "-" << target.seedStartPos+target.seedLength-1 <<  ":" << target.seedLength 
	// << ", dis: " << dis_between_src_target << ", " << expectedLength << ", " << srcMaxLength << ", "
	// << FMWalkReturnType << ".\n";
	
	// repeat seeds also lead to -1 or -4, don't reduce kmer and give up this seed
	// PB80779_18480.fa 1,596,115-1,598,135 lead to wrong extension due to over-reduction of repeat source seed
	if( !m_params.isFirst && ((FMWalkReturnType==-1 && sourceFreq > (size_t) m_params.seedKmerThreshold*3) || 
		(FMWalkReturnType==-4 && targetFreq > (size_t) m_params.seedKmerThreshold*3)) )
			FMWalkReturnType = -2;

		
	return FMWalkReturnType;
}
// refine seed interval using larger kmer
std::pair<size_t, size_t> PacBioCorrectionProcess::refineRepeatSeed(const std::string readSeq, size_t& seedStartPos, size_t& seedEndPos)
{
	// initially set to max unisnged int value
	size_t newSeedStartPos = (size_t)-1;
	size_t newSeedEndPos = (size_t)-1;
	size_t startKmerFreq=0, endKmerFreq=0;
	
	const int minRepeatFreq = 40, minFreqDiff = 30;
	
	size_t kmerSize = m_params.kmerLength;
	
	std::string kmer = readSeq.substr(seedStartPos, kmerSize);
	int initKmerFreq = BWTAlgorithms::countSequenceOccurrences(kmer, m_params.indices);
	int prevKmerFreq = initKmerFreq;

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
		bool isLargeFreqDiff = currKmerFreq - prevKmerFreq > minFreqDiff;
		
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
		if(prevKmerFreq - currKmerFreq > minFreqDiff /*|| currKmerFreq < minFreqDiff*/)
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
	
	// std::cout << newSeedStartPos << "\t" << newSeedEndPos << "\n";

	seedStartPos = newSeedStartPos;
	seedEndPos = newSeedEndPos;
	return std::make_pair(startKmerFreq, endKmerFreq);
}

// return <0: give up and break
// return 0: retry the same target
// return >0: continue to next target
int PacBioCorrectionProcess::FMWalkFailedActions(size_t& extendKmerSize, size_t& numOfTrials, 
								SeedFeature& source, SeedFeature& target, int FMWalkReturnType, int prevFMWalkReturnType)
{
	numOfTrials++;
	// extension failed due to insufficient kmers, reduce large and small kmer sizes
	if(FMWalkReturnType==-1 || FMWalkReturnType==-4)
	{
		// kmers have been enlarged due to repeats, shrink will lead to infinite loop
		if(prevFMWalkReturnType==-3 )
			return -1;
		
		// PB36993_4517.fa, AGGCTTGTCTGTAATCGGG
		if(m_params.isFirst /*&& (source.isRepeat || target.isRepeat)*/)
			return -1;		
		
		source.endBestKmerSize -=2;
		
		target.startBestKmerSize -=2;
		
		extendKmerSize -= 2;

		// std::cout << source.endBestKmerSize << "\t" << target.startBestKmerSize << "\n";

		
		// don't aggressively reduce kmer in the 1st found where most kmers are errors
		if(m_params.isFirst && (source.endBestKmerSize < 15 || target.startBestKmerSize < 15))
			return -1;
			
		if(source.endBestKmerSize < 11 || target.startBestKmerSize < 11 || extendKmerSize < 9)
			return -1;
			
		return 0;
	}
	
	// increase extendKmerSize for reducing repeats
	else if(FMWalkReturnType==-3)
	{
		if(prevFMWalkReturnType==-4 || prevFMWalkReturnType==-1)
			return -1;

		// exponential growth is required in super large repeats. Otherwise speed is too slow
		source.endBestKmerSize += pow(2, numOfTrials+1);
		target.startBestKmerSize += pow(2, numOfTrials+1);
		extendKmerSize += pow(2, numOfTrials+1);
				
		// bug: PB7017_14581_0_14647.fa
		// extendKmerSize is less than seedLength , dunno why
		if(source.seedLength < source.endBestKmerSize || target.seedLength < target.startBestKmerSize ||
			source.seedLength < extendKmerSize || target.seedLength < extendKmerSize )
			return -1;
		
		return 0;
	}
	else if(FMWalkReturnType==-2)
	{
		// probably chimera, need more observations
		// largeKmerSize = m_params.kmerLength;
		// extendKmerSize = largeKmerSize - 2;
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
		std::cout << "disBetweenSeeds: " << m_seedDis/m_totalWalkNum << std::endl << std::endl;
	}
}


// Writting results for kmerize and validate
void PacBioCorrectionPostProcess::process(const SequenceWorkItem& item, const PacBioCorrectionResult& result)
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

/***************************/
/*** Seed Feature Body *****/
/***************************/
SeedFeature::SeedFeature(size_t startPos, std::string str, bool repeat, size_t kmerSize, size_t repeatCutoff)
	:seedStartPos(startPos), seedStr(str), isRepeat(repeat), freqUpperBound(repeatCutoff), freqLowerBound(repeatCutoff/2), 
	minKmerSize(17), isPBSeed(false), isNextRepeat(false), stepSize(1)
{
	seedEndPos = seedStartPos + seedStr.length() -1;
	seedLength = seedStr.length();
	startBestKmerSize = endBestKmerSize = kmerSize<=seedLength?kmerSize:seedLength;
}

// append current seed string with extendedStr
void SeedFeature::append(std::string extendedStr)
{
	seedStr += extendedStr;
	seedLength += extendedStr.length();
	seedStartPos += extendedStr.length();
	seedEndPos += extendedStr.length();
}

void SeedFeature::setBestKmerSize(size_t kmerSize)
{
	startBestKmerSize = endBestKmerSize = kmerSize;
}

void SeedFeature::estimateBestKmerSize(const BWT* pBWT)
{			
	std::string kmerStr = seedStr.substr(0, startBestKmerSize);
	startKmerFreq = BWTAlgorithms::countSequenceOccurrences(kmerStr, pBWT);

	if(startKmerFreq > freqUpperBound)
		increaseStartKmerSize(pBWT);
	else if(startKmerFreq < freqLowerBound)
		decreaseStartKmerSize(pBWT);
		
	kmerStr = seedStr.substr(seedLength-endBestKmerSize);
	endKmerFreq = BWTAlgorithms::countSequenceOccurrences(kmerStr, pBWT);

	if(endKmerFreq > freqUpperBound)
		increaseEndKmerSize(pBWT);
	else if(endKmerFreq < freqLowerBound)
		decreaseEndKmerSize(pBWT);
	
}
	
//estimate kmer size
void SeedFeature::increaseStartKmerSize(const BWT* pBWT)
{
	while(startKmerFreq > freqUpperBound && startBestKmerSize <= seedLength - stepSize)
	{
		startBestKmerSize+=stepSize;
		std::string kmerStr = seedStr.substr(0, startBestKmerSize);
		startKmerFreq = BWTAlgorithms::countSequenceOccurrences(kmerStr, pBWT);
	}
	
	// over increase kmer size
	if(startKmerFreq < freqLowerBound)
	{
		startBestKmerSize-=stepSize;
		std::string kmerStr = seedStr.substr(0, startBestKmerSize);
		startKmerFreq = BWTAlgorithms::countSequenceOccurrences(kmerStr, pBWT);
	}
}

void SeedFeature::decreaseStartKmerSize(const BWT* pBWT)
{
	while(startKmerFreq < freqLowerBound && startBestKmerSize > minKmerSize)
	{
		startBestKmerSize-=stepSize;
		std::string kmerStr = seedStr.substr(0, startBestKmerSize);
		startKmerFreq = BWTAlgorithms::countSequenceOccurrences(kmerStr, pBWT);
	}

	// over reduce kmer size
	if(startKmerFreq>freqUpperBound)
	{
		startBestKmerSize+=stepSize;
		std::string kmerStr = seedStr.substr(0, startBestKmerSize);
		startKmerFreq = BWTAlgorithms::countSequenceOccurrences(kmerStr, pBWT);
	}
}

//estimate kmer size
void SeedFeature::increaseEndKmerSize(const BWT* pBWT)
{
	while(endKmerFreq > freqUpperBound && endBestKmerSize <= seedLength - stepSize)
	{
		endBestKmerSize+=stepSize;
		assert(seedLength >= endBestKmerSize);
		std::string kmerStr = seedStr.substr(seedLength - endBestKmerSize);
		endKmerFreq = BWTAlgorithms::countSequenceOccurrences(kmerStr, pBWT);
	}
	
	if(endKmerFreq < freqLowerBound)
	{
		endBestKmerSize-=stepSize;
		std::string kmerStr = seedStr.substr(seedLength - endBestKmerSize);
		endKmerFreq = BWTAlgorithms::countSequenceOccurrences(kmerStr, pBWT);
	}
}

void SeedFeature::decreaseEndKmerSize(const BWT* pBWT)
{
	while(endKmerFreq < freqLowerBound && endBestKmerSize > minKmerSize)
	{
		endBestKmerSize -= stepSize;
		std::string kmerStr = seedStr.substr(seedLength - endBestKmerSize);
		endKmerFreq = BWTAlgorithms::countSequenceOccurrences(kmerStr, pBWT);
	}
	
	if(endKmerFreq > freqUpperBound)
	{
		endBestKmerSize += stepSize;
		std::string kmerStr = seedStr.substr(seedLength - endBestKmerSize);
		endKmerFreq = BWTAlgorithms::countSequenceOccurrences(kmerStr, pBWT);
	}
}
