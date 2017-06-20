//----------------------------------------------
// Copyright 2015 National Chung Cheng University
// Written by Yao-Ting Huang
// Released under the GPL
//-----------------------------------------------
//
// PacBioSelfCorrectionProcess - Self-correction using FM-index walk for PacBio reads
//

#ifndef PacBioSelfCorrectionProcess_H
#define PacBioSelfCorrectionProcess_H

#include "HashMap.h"
#include "Util.h"
#include "SequenceProcessFramework.h"
#include "SequenceWorkItem.h"
#include "Metrics.h"
#include "BWTIndexSet.h"
#include "SampledSuffixArray.h"
#include "BWTAlgorithms.h"
#include "KmerDistribution.h"
#include "LongReadCorrectByOverlap.h"
#include "SeedFeature.h"

// Parameter object for the error corrector
struct PacBioSelfCorrectionParameters
{
	BWTIndexSet indices;

	int numKmerRounds;
	int kmerLength;

	// tree search parameters
	int maxLeaves;
	int minOverlap;
	int maxOverlap;
    int idmerLength;
    

	// PACBIO
	int minKmerLength;
	int FMWKmerThreshold;
	int seedKmerThreshold;
	int numOfNextTarget;
	int collectedSeeds;
    double ErrorRate;
	bool isSplit;
	bool isFirst;
	size_t maxSeedInterval;
	size_t PBcoverage;
	KmerDistribution kd;
    
    bool DebugExtend;
    bool DebugSeed;    
};




class PacBioSelfCorrectionResult
{
public:
	PacBioSelfCorrectionResult()
	: merge(false),
	totalReadsLen(0),
	correctedLen(0),
	totalSeedNum(0),
	totalWalkNum(0),
	correctedNum(0),
	highErrorNum(0),
	exceedDepthNum(0),
	exceedLeaveNum(0),
    DPNum(0),
	seedDis(0),
    Timer_Seed(0),
    Timer_FM(0),
    Timer_DP(0) {}

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
    int64_t DPNum;
	int64_t seedDis;
    double Timer_Seed;
    double Timer_FM;
    double Timer_DP;
};

//
class PacBioSelfCorrectionProcess
{
public:
	PacBioSelfCorrectionProcess(const PacBioSelfCorrectionParameters params);
	~PacBioSelfCorrectionProcess();

	// PacBio correction by Ya, v20150305.
	PacBioSelfCorrectionResult PBSelfCorrection(const SequenceWorkItem& workItem);
	PacBioSelfCorrectionResult PBSelfCorrectionUsedByPBHybridCorrection(const SequenceWorkItem& workItem);
	
	PacBioSelfCorrectionResult process(const SequenceWorkItem& workItem)
	{
		return PBSelfCorrection(workItem);
		PacBioSelfCorrectionResult result;
		return result;
	}		

private:
    FMextendParameters FMextendParameter();
    void separatebykmer(std::string readid,std::string readSeq,size_t kmerSize);
	// PacBio correction by Yao-Ting Huang, v20151208
    

    std::vector<SeedFeature> hybridSeedingFromPB(const std::string& readSeq);
    
	void initCorrect(std::string& readSeq, std::vector<SeedFeature>& seeds, std::vector<SeedFeature>& pacbioCorrectedStrs, PacBioSelfCorrectionResult& result);
	
	void realCorrect(std::string& readSeq, std::vector<SeedFeature>& seeds, std::vector<SeedFeature>& pacbioCorrectedStrs, PacBioSelfCorrectionResult& result);
	int checkseedcorrect(std::vector<SeedFeature> seeds,std::string currseed,size_t currseedStartpos);
	// kmers around repeat seeds are often error seeds, split the repeat regions into high-confident seeds
	// return kmer freq of beginning and ending kmers
	std::pair<size_t, size_t> refineRepeatSeed(const std::string readSeq, size_t& seedStartPos, size_t& seedEndPos,size_t normal_freqs);

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
							size_t smallKmerSize, size_t dis_between_src_target,FMextendParameters FMextendParameter);
                            

	std::pair<size_t,size_t> alnscore;
	PacBioSelfCorrectionParameters m_params;
    size_t m_repeat_distance = 40;
    std::string m_readid;
    double m_total_FMtime;
    double m_total_DPtime;

};

// Write the results from the overlap step to an ASQG file
class PacBioSelfCorrectionPostProcess
{
public:
	PacBioSelfCorrectionPostProcess(std::ostream* pCorrectedWriter,
	std::ostream* pDiscardWriter,
	const PacBioSelfCorrectionParameters params);

	~PacBioSelfCorrectionPostProcess();

	void process(const SequenceWorkItem& item, const PacBioSelfCorrectionResult& result);
	// void process(const SequenceWorkItemPair& itemPair, const PacBioSelfCorrectionResult& result);

private:

	std::ostream* m_pCorrectedWriter;
	std::ostream* m_pDiscardWriter;
	PacBioSelfCorrectionParameters m_params;
	
	int64_t m_totalReadsLen;
	int64_t m_correctedLen;
	int64_t m_totalSeedNum;
	int64_t m_totalWalkNum;
	int64_t m_correctedNum;
	int64_t m_highErrorNum;
	int64_t m_exceedDepthNum;
	int64_t m_exceedLeaveNum;
    int64_t m_DPNum;
	int64_t m_seedDis;
    double m_Timer_Seed;
    double m_Timer_FM;
    double m_Timer_DP;

};


#endif
