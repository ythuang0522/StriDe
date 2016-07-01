//----------------------------------------------
// Copyright 2015 National Chung Cheng University
// Written by Yao-Ting Huang
// Released under the GPL
//-----------------------------------------------

//
// Re-written from Jared Simpson's StringThreaderNode and StringThreader class
// The search tree represents a traversal through implicit FM-index graph
//
#ifndef PBOVERLAPTREE_H
#define PBOVERLAPTREE_H

#include <list>
#include "BWT.h"
#include "BWTAlgorithms.h"
#include "SAINode.h"
#include "IntervalTree.h"

class PBOverlapTree
{
    public:
         PBOverlapTree(const std::string& query,
						size_t srcKmerSize,
                       size_t minOverlap,
					   size_t maxIndelSize,
					   const BWT* pBWT, 
					   const BWT* pRBWT, 
                       double errorRate=0.7,	
					   size_t maxLeaves=256,
					   size_t seedSize=11, size_t seedDist=1, size_t repeatFreq=256);
		
        ~PBOverlapTree();

        // extend all leaves one base pair during overlap computation
        int extendOverlap(std::vector<std::string>& results);
		
		// return emptiness of leaves
		inline bool isEmpty(){return m_leaves.empty();};
		
		// return size of leaves
		inline size_t size(){return m_leaves.size();};
		
		// return size of seed
		inline size_t getSeedSize(){return m_seedSize;};

		// return size of seed
		inline size_t getCurrentLength(){return m_currentLength;};
		
        // Print all the strings represented by the tree
        void printAll();

    private:

        //
        // Functions
        //
        void extendLeaves();

        std::vector<std::pair<std::string, BWTIntervalPair> > getLeftFMIndexExtensions(SAIOverlapNode* pNode);

		// prone the leaves without seeds in proximity
		bool PrunedBySeedSupport();

        // Check if the leaves reach $
        bool isTerminated(std::vector<std::string>& results);
		bool isOverlapAcceptable(SAIOverlapNode* currNode);
		bool isSupportedByNewSeed(SAIOverlapNode* currNode, size_t smallSeedIdx, size_t largeSeedIdx);
		double computeErrorRate(const SAIOverlapNode* currNode);
	
        //
        // Data
        //
        const std::string m_Query;
        size_t m_minOverlap;
		size_t m_maxIndelSize;
		const BWT* m_pBWT;
		const BWT* m_pRBWT;
		double m_errorRate;

		// Optional parameters
        size_t m_maxLeaves;
		size_t m_seedSize;
		size_t m_seedDist;
		size_t m_repeatFreq;

        SONodePtrList m_leaves;

        SAIOverlapNode* m_pRootNode;
        SONodePtrList m_RootNodes;

        size_t m_currentLength;

		// std::vector<BWTInterval> m_TerminatedIntervals;  
        std::vector<TreeInterval<size_t> > m_TerminatedIntervals;  
		IntervalTree<size_t> BWTIntervalTree;
		
};

#endif
