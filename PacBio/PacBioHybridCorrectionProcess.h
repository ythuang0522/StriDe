//----------------------------------------------
// Copyright 2016 National Chung Cheng University
// Written by Yao-Ting Huang & Ping-Yeh Chen
// Released under the GPL
//-----------------------------------------------
//
// PacBioHybridCorrectionProcess - Hybrid correction using FM-index walk for PacBio reads
//

#ifndef PacBioHybridCorrectionProcess_H
#define PacBioHybridCorrectionProcess_H

#include "HashMap.h"
#include "Util.h"
#include "SequenceProcessFramework.h"
#include "SequenceWorkItem.h"
#include "Metrics.h"
#include "BWTIndexSet.h"
#include "SampledSuffixArray.h"
#include "BWTAlgorithms.h"
#include "KmerDistribution.h"
#include "PacBioCorrectionProcess.h"
#include "SAIPBHybridCTree.h"

// Parameter object for the error corrector
struct PacBioHybridCorrectionParameters
{
	// FM-index of high-quality short reads
	BWTIndexSet indices;

	// FM-index of low-quality long reads
	BWTIndexSet PBindices;
	
	int kmerLength;

	// tree search parameters
	int maxLeaves;
	int minOverlap;
	int maxOverlap;

	// PACBIO
	int minKmerLength;
	int FMWKmerThreshold;
	int seedKmerThreshold;

	size_t coverage;	// coverage of high-quality short reads	
	
	size_t PBKmerLength;	// kmer length used in PBself correction
	size_t PBcoverage;	// coverage of low-quality short reads	
	size_t PBSearchDepth;
	// KmerDistribution kd;
};

class PacBioHybridCorrectionResult
{
public:

	PacBioHybridCorrectionResult():
	merge(false),
	totalReadsLen(0),
	correctedLen(0),
	totalSeedNum(0),
	totalWalkNum(0),
	correctedNum(0),
	highErrorNum(0),
	exceedDepthNum(0),
	exceedLeaveNum(0),
	seedDis(0){}

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
class PacBioHybridCorrectionProcess
{
public:

	PacBioHybridCorrectionProcess(const PacBioHybridCorrectionParameters params);
	~PacBioHybridCorrectionProcess();

	// PacBio correction by Ya, v20150305.
	PacBioHybridCorrectionResult PBSelfCorrection(const SequenceWorkItem& workItem);
	PacBioHybridCorrectionResult PBHybridCorrection(const SequenceWorkItem& workItem);
	
	PacBioHybridCorrectionResult process(const SequenceWorkItem& workItem)
	{
		return PBHybridCorrection(workItem);
		PacBioHybridCorrectionResult result;
		return result;
	}		

private:

	std::vector<SeedFeature> seedingByDynamicKmer(const std::string& readSeq);
	std::vector<SeedFeature> seedingByDynamicKmer_v2(const std::string& readSeq);
	std::vector<SeedFeature> seedingByDynamicKmer_v3(const std::string& readSeq);
	int calculateKmerFreqsEachPBPos(const std::string& readSeq, std::vector<std::vector<size_t> >& PBKmerFreqsVec);
	int extendBetweenSeeds(SeedFeature& source, SeedFeature& target, std::string& strBetweenSrcTarget, int dis_between_src_target, FMWalkResult* FMWResult, int debugTargetSeed);
	void trimRepeatSeed(const std::string& readSeq, size_t coverage, size_t& seedStartPos, size_t& seedEndPos);
	bool seedingByPacBio(const std::string& readSeq, std::vector<SeedFeature>& seedVec, 	std::vector<int>& seedEndPosVec, size_t prevEndPos);
	bool seedingByPacBio_v2(const std::string& readSeq, std::vector<SeedFeature>& seedVec, 	std::vector<int>& seedEndPosVec, size_t prevEndPos);
	bool isLowComplexity(std::string& seq ,const float& ratioThreshold);
	PacBioHybridCorrectionParameters m_params;
};

// Write the results from the overlap step to an ASQG file
class PacBioHybridCorrectionPostProcess
{
public:

	PacBioHybridCorrectionPostProcess(std::ostream* pCorrectedWriter,
	std::ostream* pDiscardWriter,
	const PacBioHybridCorrectionParameters params);

	~PacBioHybridCorrectionPostProcess();

	void process(const SequenceWorkItem& item, const PacBioHybridCorrectionResult& result);
	// void process(const SequenceWorkItemPair& itemPair, const PacBioHybridCorrectionResult& result);

private:

	std::ostream* m_pCorrectedWriter;
	std::ostream* m_pDiscardWriter;
	PacBioHybridCorrectionParameters m_params;
	
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
