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

class SAIPBSelfCorrectTree
{
    public: 
		SAIPBSelfCorrectTree(const BWT* pBWT, const BWT* pRBWT, size_t min_SA_threshold=3);
		
        ~SAIPBSelfCorrectTree();

		// Collect kmers from overlapping reads via right extension of single seed
		void addHashFromSingleSeedUsingFMExtension(std::string& seedStr, size_t hashKmerSize, size_t maxLength);
		void addHashFromSingleSeedUsingLFMapping(std::string& seedStr, size_t hashKmerSize, size_t maxLength);

		// Collect kmers from overlapping reads via FM-index walk between two seeds
		int addHashFromPairedSeed(std::string seedStr, std::vector<std::pair<int, std::string> >  &targets, size_t hashKmerSize);
		
        // find feasible extension using kmers collected by addHashFromSingleSeed or addHashFromPairedSeed
		int mergeTwoSeedsUsingHash(const std::string &src, const std::string &dest, std::string &mergedseq, 
									size_t minOverlap, size_t maxLeaves, size_t minLength, size_t maxLength, size_t expectedLength);
		
		
        // Print all the strings represented by the tree
        void printAll();

    private:

        //
        // Functions
        //
		// create root node using src
		void initializeSearchTree(std::string src, size_t hashKmerSize);
		void initializeTerminalIntervals(std::string dest, size_t hashKmerSize);

		void attempToExtendUsingHash(STNodePtrList &newLeaves, size_t hashKmerSize);
        void refineSAInterval(size_t newKmer);
		std::vector<std::pair<std::string, BWTIntervalPair> > getFMIndexRightExtensions(SAIntervalNode* pNode, const size_t IntervalSizeCutoff);

        // Check if the leaves can be extended no further
		bool isTerminated(SAIntervalNodeResultVector& results, END target);
		bool isTerminated(SAIntervalNodeResultVector& results);

        //
        // Data
        //
        
        const BWT* m_pBWT;
        const BWT* m_pRBWT;
        
        size_t m_min_SA_threshold;
		
		int m_PBLen;

		SparseHashMap<std::string, std::pair<long long int, long long int> > kmerHash;

		std::vector<END> ENDs;

        SAIntervalNode* m_pRootNode;
        STNodePtrList m_leaves;


        BWTInterval m_fwdTerminatedInterval;   //in rBWT
        BWTInterval m_rvcTerminatedInterval;   //in BWT
		
		bool debug;
};

#endif
