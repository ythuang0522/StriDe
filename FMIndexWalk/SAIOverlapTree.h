//----------------------------------------------
// Copyright 2015 National Chung Cheng University
// Written by Yao-Ting Huang
// Released under the GPL
//-----------------------------------------------

//
// Re-written from Jared Simpson's StringThreaderNode and StringThreader class
// The search tree represents a traversal through implicit FM-index graph
//
#ifndef SAIOVERLAPTREE_H
#define SAIOVERLAPTREE_H

#include <list>
#include "BWT.h"
#include "BWTAlgorithms.h"
#include "SAINode.h"
#include "OverlapBlock.h"


class SAIOverlapTree
{
    public:
         SAIOverlapTree(const std::string& query,
                       size_t minOverlap,
					   size_t maxIndelSize,
					   const BWT* pBWT, 
					   const BWT* pRBWT, 
					   AlignFlags af,
                       size_t maxLeaves=256,
					   size_t seedSize=17, size_t seedDist=1, size_t repeatFreq=256);
		
        ~SAIOverlapTree();

        // extend all leaves one base pair during overlap computation
        int extendOverlapOneBase(std::vector<OverlapBlock>& results);

        // terminate leaves which exceeds m_Query length
        bool terminateContainedBlocks(std::vector<OverlapBlock>& results);
		
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

		// add new roots for tolerating init errors
		void addNewRootNodes();

        void attempToExtend(SONodePtrList &newLeaves);

        std::vector<std::pair<std::string, BWTIntervalPair> > getLeftFMIndexExtensions(SAIOverlapNode* pNode);
        std::vector<std::pair<std::string, BWTIntervalPair> > getRightFMIndexExtensions(SAIOverlapNode* pNode);

		// prone the leaves without seeds in proximity
		bool PrunedBySeedSupport();

        // Check if the leaves reach $
        bool isTerminated(std::vector<OverlapBlock>& results);
		bool isOverlapAcceptable(SAIOverlapNode* currNode);
		bool isSupportedByNewSeed(SAIOverlapNode* currNode, size_t largeSeedIdx);
		double computeErrorRate(SAIOverlapNode* currNode);

		size_t computeRightExtensionLength(SAIOverlapNode* currNode);

		// extend to left or right extreme upto length
		std::list<BWTIntervalPair> extendToLeftExtreme(BWTIntervalPair& currIntervalPair, int length, bool& isLeftSubstring);
		std::list<BWTIntervalPair> extendToRightExtreme(BWTIntervalPair& currIntervalPair, int length);
	
		// extend to left or right extreme upto length and collect pre-terminated reads
		std::list<BWTIntervalPair> collectToRightExtreme(BWTIntervalPair& currIntervalPair, int length, std::list<BWTIntervalPair>& terminatedReads);

        //
        // Data
        //
        const std::string m_Query;
        size_t m_minOverlap;
		size_t m_maxIndelSize;
		double m_errorRate;
		const BWT* m_pBWT;
		const BWT* m_pRBWT;
		AlignFlags m_af;

		// Optional parameters
        size_t m_maxLeaves;
		size_t m_seedSize;
		size_t m_seedDist;
		size_t m_repeatFreq;

        SONodePtrList m_leaves;

        // SAIOverlapNode* m_pRootNode;
        SONodePtrList m_RootNodes;

        size_t m_currentLength;

        std::vector<BWTInterval> m_TerminatedIntervals;  
};

#endif
