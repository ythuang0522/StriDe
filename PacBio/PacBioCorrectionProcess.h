//----------------------------------------------
// Copyright 2015 National Chung Cheng University
// Written by Yao-Ting Huang
// Released under the GPL
//-----------------------------------------------
//
// PacBioCorrectionProcess - Self-correction using FM-index walk for PacBio reads
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
		SeedFeature(size_t startPos, std::string str, bool repeat, size_t kmerSize, size_t repeatCutoff);

		SeedFeature(){};
		
		// append current seed string with extendedStr
		void append(std::string extendedStr);
		void setBestKmerSize(size_t kmerSize);		
		void estimateBestKmerSize(const BWT* pBWT);
		bool inline isSmall(){return seedLength<=17?true:false ;}
		
		size_t seedStartPos;
		size_t seedEndPos;
		size_t seedLength;
		std::string seedStr;
		bool isRepeat;
		bool isPBSeed;
		bool isNextRepeat;
		
		// estimated by calling estimateBestKmerSize
		size_t startBestKmerSize;
		size_t endBestKmerSize;
		size_t startKmerFreq;
		size_t endKmerFreq;
		
	private:
		size_t freqUpperBound;
		size_t freqLowerBound;
		size_t minKmerSize;
		size_t stepSize;
		//estimate kmer size
		void increaseStartKmerSize(const BWT* pBWT);
		void decreaseStartKmerSize(const BWT* pBWT);

		//estimate kmer size
		void increaseEndKmerSize(const BWT* pBWT);
		void decreaseEndKmerSize(const BWT* pBWT);
};

enum PacBioCorrectionAlgorithm
{
	PBC_SELF	// PacBio self correction
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
	
	PacBioCorrectionResult process(const SequenceWorkItem& workItem)
	{
		switch(m_params.algorithm)
		{
		case PBC_SELF:
			{
				return PBSelfCorrection(workItem);
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

	// PacBio correction by Yao-Ting Huang, v20151208
	std::vector<SeedFeature> seedingByFixedKmer(const std::string readSeq, size_t contaminatedCutoff=256);
	
	std::vector<SeedFeature> seedingByDynamicKmer(const std::string readSeq);

	void initCorrect(std::string& readSeq, std::vector<SeedFeature>& seeds, std::vector<SeedFeature>& pacbioCorrectedStrs, PacBioCorrectionResult& result);
	
	void realCorrect(std::string& readSeq, std::vector<SeedFeature>& seeds, std::vector<SeedFeature>& pacbioCorrectedStrs, PacBioCorrectionResult& result);
	
	// kmers around repeat seeds are often error seeds, split the repeat regions into high-confident seeds
	// return kmer freq of beginning and ending kmers
	std::pair<size_t, size_t> refineRepeatSeed(const std::string readSeq, size_t& seedStartPos, size_t& seedEndPos);

	// return complexity of seq, default: 0.9
	bool  isLowComplexity (std::string seq, float & GCratio, float threshold=0.7);

	// return <0: give up and break
	// return 0: retry the same target
	// return >0: continue to next target
	int  FMWalkFailedActions (size_t& smallKmerSize, size_t& numOfTrials, 
						SeedFeature& sourceStrLen, SeedFeature& target, int FMWalkReturnType, int prevFMWalkReturnType);


	// Perform FMindex extension between source and target seeds
	// Return FMWalkReturnType
	int extendBetweenSeeds(SeedFeature& source, SeedFeature& target, std::string& rawSeq, std::string& mergedseq,
							size_t smallKmerSize, size_t dis_between_src_target);
	
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
