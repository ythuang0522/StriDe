//----------------------------------------------
// Copyright 2016 National Chung Cheng University
// Written by Yao-Ting Huang & Ping-Yeh Chen
// Released under the GPL
//-----------------------------------------------

#ifndef SAIPBHybridCTree_H
#define SAIPBHybridCTree_H

#include <list>
#include "BWT.h"
#include "BWTAlgorithms.h"
#include "SAINode.h"
#include "HashMap.h"

// Parameter object for the FM-index Walk of PacBio Hybrid Correction
struct PBHCFMWalkParameters
{
	BWTIndexSet indices;
	std::string sourceSeed;
	std::string strBetweenSrcTarget;
	std::string targetSeed;
    int disBetweenSrcTarget;
	int maxLeaves;
	int minOverlap;
	int maxOverlap;
    int SAThreshold = 3;
	bool kmerMode = false;
	bool lowCoverageHighErrorMode = false;
	bool debugMode = false;
};

class SAIntervalPBHybridCTree
{
    public:
        SAIntervalPBHybridCTree(PBHCFMWalkParameters parameters);

        ~SAIntervalPBHybridCTree();

        //return the merged string
        //bool mergeTwoReads(StringVector & mergeReads);
        int mergeTwoSeeds(std::string &mergedseq);
		size_t getKmerCoverage(){return m_maxKmerCoverage;};
		size_t getMaxUsedLeaves(){return m_maxUsedLeaves;};

        // Print all the strings represented by the tree
        void printAll();

    private:

        //
        // Functions
        //
        void extendLeaves();
        void attempToExtend(STNodePtrList &newLeaves);
        void refineSAInterval(size_t newKmer);
        std::vector<std::pair<std::string, BWTIntervalPair> > getFMIndexExtensions(SAIntervalNode* pNode);

        // Check if the leaves can be extended no further
        bool isTerminated(SAIntervalNodeResultVector& results);
        size_t calculateKmerCoverage (const std::string & seq , size_t kmerLength , const BWT* pBWT);

        void removeLeavesByRepeatKmer();
		int findTheBestPath(SAIntervalNodeResultVector results, std::string &mergedseq);

        //
        // Data
        //
        const std::string* m_pSourceSeed;
		std::string m_strBetweenSrcTarget;
        std::string m_targetSeed;
        size_t m_minOverlap;
		size_t m_maxOverlap;
        int m_MaxLength;
		int m_MinLength;
        size_t m_MaxLeaves;
        const BWT* m_pBWT;
        const BWT* m_pRBWT;
        size_t m_min_SA_threshold;
		int m_disBetweenSrcTarget;
		int m_expectedLength;

        SAIntervalNode* m_pRootNode;
        STNodePtrList m_leaves;
        size_t m_currentLength;
		size_t m_currentKmerSize;
		size_t m_maxKmerCoverage;
		size_t m_maxUsedLeaves;

        BWTInterval m_fwdTerminatedInterval;   //in rBWT
        BWTInterval m_rvcTerminatedInterval;   //in BWT
		
		size_t m_beginningIntervalSize;
		size_t m_terminatedIntervalSize;
		
        bool m_kmerMode;
		bool m_lowCoverageHighErrorMode;
		bool m_repeatMode=false;
		bool m_debugMode;
};

#endif
