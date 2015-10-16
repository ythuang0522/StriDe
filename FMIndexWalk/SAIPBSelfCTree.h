//----------------------------------------------
// Copyright 2015 National Chung Cheng University
// Written by Yao-Ting Huang & Ping-Yeh Chen
// Released under the GPL
//-----------------------------------------------

#ifndef SAIPBSelfCTree_H
#define SAIPBSelfCTree_H

#include <list>
#include "BWT.h"
#include "BWTAlgorithms.h"
#include "SAINode.h"
#include "HashMap.h"

struct END
{
	std::string endingKmer;
	int maxLength;
	int minLength;
    BWTInterval fwdTerminatedInterval;   //in rBWT
    BWTInterval rvcTerminatedInterval;   //in BWT
};

class SAIntervalTreeForPBGap
{
    public: 
		SAIntervalTreeForPBGap(const std::string* pQuery,
                       size_t minOverlap,
                       std::vector<std::pair<int, std::string> > targets,
                       size_t MaxLeaves,
                       const BWT* pBWT,
                       const BWT* pRBWT,
                       size_t SA_threshold=3);
		
        ~SAIntervalTreeForPBGap();

        //return the merged string
		int mergeTwoSeedsUsingHash(std::string &mergedseq);
		int buildHash(std::vector<END> ENDs);
		size_t getKmerCoverage(){return m_maxKmerCoverage;};
		size_t getMaxUsedLeaves(){return m_maxUsedLeaves;};

        // Print all the strings represented by the tree
        void printAll();

    private:

        //
        // Functions
        //
		void attempToExtendUsingHash(STNodePtrList &newLeaves);
        void refineSAInterval(size_t newKmer);
		std::vector<std::pair<std::string, BWTIntervalPair> > getFMIndexExtensionsForHash(SAIntervalNode* pNode);

        // Check if the leaves can be extended no further
		bool isTerminated(SAIntervalNodeResultVector& results, END target);
		bool isTerminated(SAIntervalNodeResultVector& results);

        //
        // Data
        //
        const std::string* m_pQuery;
        size_t m_minOverlap;
		size_t m_maxOverlap;
        int m_MaxLength;
		int m_MinLength;
        size_t m_MaxLeaves;
        const BWT* m_pBWT;
        const BWT* m_pRBWT;
        std::string m_secondread;
        size_t m_min_SA_threshold;
        bool m_kmerMode;
		SparseHashMap<std::string, std::pair<long long int, long long int> > kmerHash;
		bool debug;
		int m_PBLen;
		std::vector<END> ENDs;

        SAIntervalNode* m_pRootNode;
        STNodePtrList m_leaves;
        DenseHashMap<std::string, size_t, StringHasher> m_KmerIndexMap;
        int m_currentLength;
		size_t m_currentKmerSize;
		size_t m_maxKmerCoverage;
		size_t m_maxUsedLeaves;

        BWTInterval m_fwdTerminatedInterval;   //in rBWT
        BWTInterval m_rvcTerminatedInterval;   //in BWT
};

#endif
