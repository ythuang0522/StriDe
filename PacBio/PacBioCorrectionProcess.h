//----------------------------------------------
// Copyright 2015 National Chung Cheng University
// Written by Yao-Ting Huang
// Released under the GPL
//-----------------------------------------------
//
// PacBioCorrectionProcess - Self-correction or hybrid correction using FM-index walk for PacBio reads
//
//

#ifndef PacBioCorrectionProcess_H
#define PacBioCorrectionProcess_H

#include "HashMap.h"
#include "Util.h"
#include "SequenceProcessFramework.h"
#include "SequenceWorkItem.h"
#include "Metrics.h"
#include "BWTIndexSet.h"
#include "SampledSuffixArray.h"
#include "BWTAlgorithms.h"
#include "KmerDistribution.h"

struct SeedFeature
{
	public:
		SeedFeature(size_t startPos, std::string str, bool repeat=false):seedStartPos(startPos), seedStr(str), isRepeat(repeat)
		{
			seedEndPos = seedStartPos + seedStr.length() -1;
			seedLength = seedStr.length();
		}
		
		// append current seed string with extendedStr
		void append(std::string extendedStr)
		{
			seedStr += extendedStr;
			seedLength += extendedStr.length();
			seedEndPos += extendedStr.length();
		}
		
		size_t seedStartPos;
		size_t seedEndPos;
		size_t seedLength;
		std::string seedStr;
		size_t startKmerFreq;
		size_t endKmerFreq;
		bool isRepeat;
};

enum PacBioCorrectionAlgorithm
{
	PBC_SELF,	// PacBio self correction
	PBC_HYBRID	// PacBio Hybrid Correction
};


// Parameter object for the error corrector
struct PacBioCorrectionParameters
{
	PacBioCorrectionAlgorithm algorithm;
	BWTIndexSet indices;

	int numKmerRounds;
	int kmerLength;

	// tree search parameters
	int maxLeaves;
	int minOverlap;
	int maxOverlap;
	

	// PACBIO
	int minKmerLength;
	int FMWKmerThreshold;
	int seedKmerThreshold;
	int numOfNextTarget;
	int collectedSeeds;
	std::vector<int> seedWalkDistance;

	bool isSplit;
	bool isFirst;
	size_t maxSeedInterval;
	
	KmerDistribution kd;
};


class PacBioCorrectionResult
{
public:
	PacBioCorrectionResult()
	: merge(false),
	totalReadsLen(0),
	correctedLen(0),
	totalSeedNum(0),
	totalWalkNum(0),
	correctedNum(0),
	highErrorNum(0),
	exceedDepthNum(0),
	exceedLeaveNum(0),
	seedDis(0) {}

	DNAString correctSequence;
	
	bool merge;
	
	size_t kmerLength;

	// PacBio reads correction by Ya, v20151001.
	std::vector<DNAString> correctedPacbioStrs;
	int64_t totalReadsLen;
	int64_t correctedLen;
	int64_t totalSeedNum;
	int64_t totalWalkNum;
	int64_t correctedNum;
	int64_t highErrorNum;
	int64_t exceedDepthNum;
	int64_t exceedLeaveNum;
	int64_t seedDis;
};

//
class PacBioCorrectionProcess
{
public:
	PacBioCorrectionProcess(const PacBioCorrectionParameters params);
	~PacBioCorrectionProcess();

	// PacBio correction by Ya, v20150305.
	PacBioCorrectionResult PBSelfCorrection(const SequenceWorkItem& workItem);
	PacBioCorrectionResult PBHybridCorrection(const SequenceWorkItem& workItem);

	PacBioCorrectionResult process(const SequenceWorkItem& workItem)
	{
		switch(m_params.algorithm)
		{
		case PBC_SELF:
			{
				return PBSelfCorrection(workItem);
				break;
			}
		case PBC_HYBRID:
			{
				return PBHybridCorrection(workItem);
				break;
			}

		default:
			{
				std::cout << "Unsupported algorithm\n";
				assert(false);
			}
		}
		PacBioCorrectionResult result;
		return result;
	}		

private:

	// PacBio correction by Ya, v20150305.
	// std::vector<std::pair<int, std::string> > searchingSeedsUsingSolidKmer(const std::string readSeq, size_t contaminatedCutoff=128);
	std::vector<SeedFeature> searchingSeedsUsingSolidKmer(const std::string readSeq, size_t contaminatedCutoff=128);

	std::pair<size_t, size_t> refineRepeatSeed(const std::string readSeq, size_t repeatKmerSize, size_t seedStartPos, size_t seedEndPos);

	// return complexity of seq, default: 0.9
	bool  isLowComplexity (std::string seq, float & GCratio, float threshold=0.7);

	// return <0: give up and break
	// return 0: retry the same target
	// return >0: continue to next target
	int  FailureActions (size_t& largeKmerSize, size_t& smallKmerSize, size_t& numOfTrials, 
						size_t sourceStrLen, size_t targetStrLen, int FMWalkReturnType, int prevFMWalkReturnType);


	// Perform FMindex extension between source and target seeds
	// Return FMWalkReturnType
	int extendBetweenSeeds(SeedFeature& source, SeedFeature& target, std::string& mergedseq,
						size_t largeKmerSize, size_t smallKmerSize, size_t dis_between_src_target);
	
	std::vector<std::pair<int,std::string> > findSeedsUsingDynamicKmerLen(const std::string readSeq);
	int doubleFMWalkForPacbio(std::pair<int,std::string> firstSeed, std::pair<int,std::string> secondSeed, int minOverlap, int needWalkLen, std::string* mergedseq);
	int solveHighError(std::pair<int,std::string> firstSeed, std::pair<int,std::string> secondSeed, int minOverlap, int needWalkLen, std::string* mergedseq);
	
	PacBioCorrectionParameters m_params;

};

// Write the results from the overlap step to an ASQG file
class PacBioCorrectionPostProcess
{
public:
	PacBioCorrectionPostProcess(std::ostream* pCorrectedWriter,
	std::ostream* pDiscardWriter,
	const PacBioCorrectionParameters params);

	~PacBioCorrectionPostProcess();

	void process(const SequenceWorkItem& item, const PacBioCorrectionResult& result);
	// void process(const SequenceWorkItemPair& itemPair, const PacBioCorrectionResult& result);

private:

	std::ostream* m_pCorrectedWriter;
	std::ostream* m_pDiscardWriter;
	PacBioCorrectionParameters m_params;
	
	int64_t m_totalReadsLen;
	int64_t m_correctedLen;
	int64_t m_totalSeedNum;
	int64_t m_totalWalkNum;
	int64_t m_correctedNum;
	int64_t m_highErrorNum;
	int64_t m_exceedDepthNum;
	int64_t m_exceedLeaveNum;
	int64_t m_seedDis;

};


#endif
