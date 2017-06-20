//----------------------------------------------
// Copyright 2016 National Chung Cheng University
// Written by Yao-Ting Huang & Ping-Yeh Chen
// Released under the GPL
//-----------------------------------------------
//
// Re-written from Jared Simpson's StringThreaderNode and StringThreader class
// The search tree represents a traversal through implicit FM-index graph
//
#ifndef OverlapTree_H
#define OverlapTree_H

#include <list>
#include "BWT.h"
#include "BWTAlgorithms.h"
#include "SAINode.h"
#include "IntervalTree.h"
#include "HashtableSearch.h"


struct FMWalkResult2
{
	std::string mergedSeq;
	int alnScore;
	double kmerFreq;
};
struct FMextendParameters
{
    
    public:
    FMextendParameters(BWTIndexSet indices, int idmerLength, int maxLeaves,int minKmerLength,size_t PBcoverage, double ErrorRate,bool Debug):
    indices(indices),
    idmerLength(idmerLength),
    maxLeaves(maxLeaves),
    minKmerLength(minKmerLength),
    PBcoverage(PBcoverage),
    ErrorRate(ErrorRate),
    Debug(Debug)
    {};

    FMextendParameters(){};
    BWTIndexSet indices;
    int idmerLength;
    bool Debug;
    double ErrorRate;
    int maxLeaves;
    int minKmerLength;
    size_t PBcoverage;
};

class LongReadSelfCorrectByOverlap
{
    public:
    LongReadSelfCorrectByOverlap();
         LongReadSelfCorrectByOverlap(
                const std::string& sourceSeed,
				const std::string& strBetweenSrcTarget,
				const std::string& targetSeed,
				int m_disBetweenSrcTarget,
                size_t initkmersize,
                size_t maxOverlap,
                const FMextendParameters params,
				size_t m_min_SA_threshold = 3,
				double errorRate = 0.25,	
				size_t repeatFreq = 256,
                size_t localSimilarlykmerSize = 100
               );
		
        ~LongReadSelfCorrectByOverlap();

        // extend all leaves one base pair during overlap computation
        int extendOverlap(FMWalkResult2 &FMWResult);
		
		// return emptiness of leaves
		inline bool isEmpty(){return m_leaves.empty();};
		
		// return size of leaves
		inline size_t size(){return m_leaves.size();};
		
		// return size of seed
		inline size_t getSeedSize(){return m_seedSize;};

		// return size of seed
		inline size_t getCurrentLength(){return m_currentLength;};
		

        

        size_t SelectFreqsOfrange(size_t LowerBound,size_t UpBound,SONode3PtrList &newLeaves);
        

        std::pair<size_t,size_t> alnscore;
    private:

        //
        // Functions
        //
        
         void initialRootNode(std::string beginningkmer);
        void buildOverlapbyFMindex(std::string beginningkmer);
        
        

        SONode3PtrList extendLeaves();

        char transToRvc(char b);

		void attempToExtend(SONode3PtrList &newLeaves);
        void updateLeaves(SONode3PtrList &newLeaves,std::vector< std::pair<std::string, BWTIntervalPair> > &extensions,SAIOverlapNode3* pNode);
        
		void refineSAInterval(size_t newKmer);
        
		int findTheBestPath(SAIntervalNodeResultVector results, FMWalkResult2 &FMWResult);

		
        std::vector<std::pair<std::string, BWTIntervalPair> > getFMIndexExtensions(SAIOverlapNode3* pNode);
     

		// prone the leaves without seeds in proximity
		bool PrunedBySeedSupport(SONode3PtrList &newLeaves);
        //Check if need reduce kmer size
        bool isInsufficientFreqs(SONode3PtrList &newLeaves);
        // Check if the leaves reach $
        bool isTerminated(SAIntervalNodeResultVector& results);
		bool isOverlapAcceptable(SAIOverlapNode3* currNode);
		bool isSupportedByNewSeed(SAIOverlapNode3* currNode, size_t smallSeedIdx, size_t largeSeedIdx);
		double computeErrorRate(SAIOverlapNode3* currNode);
	
        //
        // Data
        //
		const std::string m_sourceSeed;
		const std::string m_strBetweenSrcTarget;
		const std::string m_targetSeed;	
        const bool m_isDebug;
		int m_disBetweenSrcTarget;
        size_t m_minOverlap;
		size_t m_maxOverlap;
        size_t m_initkmersize;
		size_t m_maxIndelSize;
        BWTIndexSet m_BWTindices;
		const BWT* m_pBWT;
		const BWT* m_pRBWT;
        const size_t m_PBcoverage;
        std::vector<double> freqsOfKmerSize;
		// Optional parameters
        size_t m_min_SA_threshold;

		double m_errorRate;
        double m_PacBioErrorRate;
        size_t m_maxLeaves;
		size_t m_seedSize;
		size_t m_repeatFreq;
        size_t m_localSimilarlykmerSize;
        


        
		std::string m_query;
		
        int m_maxLength;
		int m_minLength;
        std::vector<BWTInterval> m_fwdTerminatedInterval;   //in rBWT
        std::vector<BWTInterval> m_rvcTerminatedInterval;   //in BWT
		
        SONode3PtrList m_leaves;

        SAIOverlapNode3* m_pRootNode;
        SONode3PtrList m_RootNodes;
        
        

        size_t m_currentLength;
		size_t m_currentKmerSize;
        
        std::vector<TreeInterval<size_t> > m_fwdIntervals;
		std::vector<TreeInterval<size_t> > m_rvcIntervals;
		IntervalTree<size_t> fwdIntervalTree;
		IntervalTree<size_t> rvcIntervalTree;
        

        // HashtableSearch *hashIndex;

        size_t RemainedMaxLength;
        
};

#endif