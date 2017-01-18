//----------------------------------------------
// Copyright 2016 National Chung Cheng University
// Written by Yao-Ting Huang & Ping-Yeh Chen
// Released under the GPL
//-----------------------------------------------
//
// Re-written from Jared Simpson's StringThreaderNode and StringThreader class
// The search tree represents a traversal through implicit FM-index graph
//
#ifndef ShortReadOverlapTree_H
#define ShortReadOverlapTree_H

#include <list>
#include "BWT.h"
#include "BWTAlgorithms.h"
#include "SAINode.h"
#include "IntervalTree.h"
#include "SAIPBHybridCTree.h"

class ShortReadOverlapTree
{
    public:
         ShortReadOverlapTree(FMWalkParameters& parameters);
		
        ~ShortReadOverlapTree();

        // extend all leaves one base pair during overlap computation
        int extendOverlap(FMWalkResult &FMWResult);
		
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
		void extendLeaves_v2();
		void attempToExtend(SONode2PtrList &newLeaves);
		void refineSAInterval(size_t newKmer);
		int findTheBestPath(SAIntervalNodeResultVector results, FMWalkResult &FMWResult);
		
        std::vector<std::pair<std::string, BWTIntervalPair> > getFMIndexExtensions(SAIOverlapNode2* pNode);

		// prone the leaves without seeds in proximity
		bool PrunedBySeedSupport();

        // Check if the leaves reach $
        bool isTerminated(SAIntervalNodeResultVector& results);
		bool isOverlapAcceptable(SAIOverlapNode2* currNode);
		bool isSupportedByNewSeed(SAIOverlapNode2* currNode, size_t smallSeedIdx, size_t largeSeedIdx);
		double computeErrorRate(const SAIOverlapNode2* currNode);
	
        //
        // Data
        //
		const std::string m_sourceSeed;
		const std::string m_rawPBStrBetweenSrcTargetWith2Minoverlap;
		const std::string m_targetSeed;		
		int m_disBetweenSrcTarget;
        size_t m_minOverlap;
		size_t m_maxOverlap;
		size_t m_maxIndelSize;
		const BWT* m_pBWT;
		const BWT* m_pRBWT;

		// Optional parameters
        size_t m_min_SA_threshold;
		double m_errorRate;
        size_t m_maxLeaves;
		size_t m_seedSize;
		size_t m_repeatFreq;
		
        int m_maxLength;
		int m_minLength;
        BWTInterval m_fwdTerminatedInterval;   //in rBWT
        BWTInterval m_rvcTerminatedInterval;   //in BWT
		
        SONode2PtrList m_leaves;

        SAIOverlapNode2* m_pRootNode;
        SONode2PtrList m_RootNodes;

        size_t m_currentLength;
		size_t m_currentKmerSize;

        std::vector<TreeInterval<size_t> > m_fwdIntervals;
		std::vector<TreeInterval<size_t> > m_rvcIntervals;
		IntervalTree<size_t> fwdIntervalTree;
		IntervalTree<size_t> rvcIntervalTree;
};

#endif