//----------------------------------------------
// Copyright 2015 National Chung Cheng University
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

class SAIntervalPBHybridCTree
{
    public:
        SAIntervalPBHybridCTree(const std::string* pQuery,
                       size_t minOverlap,
					   size_t maxOverlap,
                       int MaxLength,
                       size_t MaxLeaves,
                       const BWT* pBWT,
                       const BWT* pRBWT,
                       std::string secondread,
                       size_t SA_threshold=3,
                       bool KmerMode=false);

        ~SAIntervalPBHybridCTree();

        //return the merged string
        //bool mergeTwoReads(StringVector & mergeReads);
        int mergeTwoReads(std::string &mergedseq);
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

        SAIntervalNode* m_pRootNode;
        STNodePtrList m_leaves;
        DenseHashMap<std::string, size_t, StringHasher> m_KmerIndexMap;
        size_t m_currentLength;
		size_t m_currentKmerSize;
		size_t m_maxKmerCoverage;
		size_t m_maxUsedLeaves;

        BWTInterval m_fwdTerminatedInterval;   //in rBWT
        BWTInterval m_rvcTerminatedInterval;   //in BWT
};

#endif
