//----------------------------------------------
// Copyright 2015 National Chung Cheng University
// Written by Yao-Ting Huang & Ping-Yeh Chen
// Released under the GPL
//-----------------------------------------------

#ifndef RollingPBSelfCTree_H
#define RollingPBSelfCTree_H

#include <list>
#include "BWT.h"
#include "BWTAlgorithms.h"
#include "SAIPBSelfCTree.h"
#include "RollingNode.h"
#include "HashMap.h"


// Return a hash key for a Kmer
// The real hash function is implemented using rolling hash
struct RollingHasher
{
    size_t operator()(const size_t& key) const { return key; }
};


class RollingPBSelfCTree
{
    public: 
		RollingPBSelfCTree(const BWT* pBWT, const BWT* pRBWT, 
					size_t smallKmerSize, size_t min_SA_threshold=3, int m_maxLeavesAllowed=64);
		
        ~RollingPBSelfCTree();

		// Collect kmers from overlapping reads via right extension of single seed
		// void addHashFromSingleSeedUsingFMExtension(std::string& seedStr, size_t hashKmerSize, size_t maxLength);
		
		// LF-mapping of each SA index in the interval independently using loop instead of BFS tree expansion
		// Contaminated reads often are simple repeats C* or T* with large freq
		// Return false if the freq is way too large
		bool addHashBySingleSeed(std::string& seedStr, size_t largeKmerSize, size_t smallKmerSize, size_t maxLength, int expectedLength=-1);

		// Collect kmers from overlapping reads via FM-index walk between two seeds
		// int addHashFromPairedSeed(std::string seedStr, std::vector<std::pair<int, std::string> >  &targets, size_t hashKmerSize);
		
        // find feasible extension using kmers collected by addHashFromSingleSeed or addHashFromPairedSeed
		int mergeTwoSeedsUsingHash(const std::string &src, const std::string &dest, std::string &mergedseq, 
								size_t minLength, size_t maxLength, size_t expectedLength);
		
		
        // Print all the strings represented by the tree
        void printAll();
		void printLeaves(size_t hashKmerSize);

    private:

        //
        // Functions
        //
		// create root node using src
		void initializeSearchTree(std::string src, size_t hashKmerSize);
		void initializeTerminalHashValues(std::string dest, size_t hashKmerSize);

		void attempToExtendUsingRollingHash(RollingNodePtrList &newLeaves, size_t hashKmerSize);
		
        // Check if the leaves can be extended no further
		bool isTerminated(SAIntervalNodeResultVector& results);

        //
        // Data
        //
        
        const BWT* m_pBWT;
        const BWT* m_pRBWT;
        
		size_t m_smallKmerSize;
        size_t m_min_SA_threshold;
		size_t m_maxLeavesAllowed;
		
		int m_expectedLength;
		int m_currentLength;
		int m_seedLength;
		
		// std::vector<END> ENDs;

        RollingNode* m_pRootNode;
        RollingNodePtrList m_leaves;

        size_t m_fwdTerminatedValue;   //in rBWT
        size_t m_rvcTerminatedValue;   //in BWT

		
		// DenseHashMap<size_t, KmerFeatures*, RollingHasher> kmerHash;
		// typedef DenseHashMap<size_t, KmerFeatures*>::iterator kmerHashiter;

		HashMap<size_t, KmerFeatures*, RollingHasher> kmerHash;
		typedef HashMap<size_t, KmerFeatures*>::iterator kmerHashiter;

		// m_rollinghasher can be only one instance in each tree due to mapping by random number generator
		CyclicHash<size_t, unsigned char> *m_rollinghasher;
};


#endif
