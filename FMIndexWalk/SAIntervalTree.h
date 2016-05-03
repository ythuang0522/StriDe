//----------------------------------------------
// Copyright 2015 National Chung Cheng University
// Written by Yao-Ting Huang
// Released under the GPL
//-----------------------------------------------

//
// Re-written from Jared Simpson's StringThreaderNode and StringThreader class
// The search tree represents a traversal through implicit FM-index graph
//
#ifndef SAINTERVALTREE_H
#define SAINTERVALTREE_H

#include <list>
#include "BWT.h"
#include "BWTAlgorithms.h"
#include "SAINode.h"

class SAIntervalTree
{
    public:
        SAIntervalTree(const std::string* pQuery,
                       size_t minOverlap,
					   size_t maxOverlap,
                       size_t MaxLength,
                       size_t MaxLeaves,
					   BWTIndexSet indices,
                       std::string secondread,
                       size_t SA_threshold=3,
                       bool KmerMode=false);
		
		//for validation purpose
		SAIntervalTree(const std::string* pQuery,
                       size_t minOverlap,
					   size_t maxOverlap,
                       size_t MaxLength,
                       size_t MaxLeaves,
					   BWTIndexSet indices,
                       size_t SA_threshold=3,
                       bool KmerMode=false);

        ~SAIntervalTree();

        //return the merged string
        int mergeTwoReads(std::string &mergedseq);
		
		// validate if each kmer in the read is suuported by at least m_minOverlap overlap 
		int validate(std::string &mergedseq);
		size_t getKmerCoverage(){return m_maxKmerCoverage;};
		size_t getMaxUsedLeaves(){return m_maxUsedLeaves;};
		bool isBubbleCollapsed(){return m_isBubbleCollapsed;}

        // Print all the strings represented by the tree
        void printAll();

    private:

        //
        // Functions
        //
        void extendLeaves();
        void attempToExtend(STNodePtrList &newLeaves);
        void refineSAInterval(size_t newKmerSize);
        std::vector<std::pair<std::string, BWTIntervalPair> > getFMIndexExtensions(SAIntervalNode* pNode);

        // Check if the leaves can be extended no further
        bool isTerminated(SAIntervalNodeResultVector& results);
        bool isTwoReadsOverlap(std::string & mergedseq);
        size_t calculateKmerCoverage (const std::string & seq , size_t kmerLength , const BWT* pBWT);
		bool replaceLowFreqKmer (std::string & seq , size_t kmerLength);

        void removeLeavesByRepeatKmer();

        //
        // Data
        //
        const std::string* m_pQuery;
        size_t m_minOverlap;
		size_t m_maxOverlap;
        size_t m_MaxLength;
        size_t m_MaxLeaves;
		BWTIndexSet m_indices;
        std::string m_secondread;
        size_t m_min_SA_threshold;
        bool m_kmerMode;

        SAIntervalNode* m_pRootNode;
        STNodePtrList m_leaves;

        size_t m_currentLength;
		size_t m_currentKmerSize;
		size_t m_maxKmerCoverage;
		size_t m_maxUsedLeaves;
		bool m_isBubbleCollapsed;

        BWTInterval m_fwdTerminatedInterval;   //in rBWT
        BWTInterval m_rvcTerminatedInterval;   //in BWT
};

#endif
