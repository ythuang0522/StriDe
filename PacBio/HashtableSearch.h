//----------------------------------------------
// Copyright 2015 National Chung Cheng University
// Written by Yao-Ting Huang & Ping-Yeh Chen
// Released under the GPL
//-----------------------------------------------

#ifndef HashtableSearch_H
#define HashtableSearch_H

#include <list>

#include "BWT.h"
#include "BWTAlgorithms.h"
#include "SAINode.h"
#include "HashMap.h"




struct END2
{
	std::string endingKmer;
	int maxLength;
	int minLength;
    BWTInterval fwdTerminatedInterval;   //in rBWT
    BWTInterval rvcTerminatedInterval;   //in BWT
};

// class for storing various features of each kmer in the hashtable, including pos, freq, ...
class KmerFeatures2{
	private:
		std::vector<long long int> m_sumOfFreq;
		std::vector<long long int> m_sumOfPos;
		long long int m_intervalSize;
		long long int m_totalFreq;
		long long int m_totalSum;
		double m_maxAvgFreq;
		bool m_isVisited;

    public: 
		// maxIntervalSize is overestimated to be the max length
		KmerFeatures2(long long int pos, size_t maxIntervalSize, size_t intervalSize=35):
			m_intervalSize(intervalSize), 
			m_totalFreq(0),m_totalSum(0), m_maxAvgFreq(0), m_isVisited(false)
			{
				m_sumOfFreq.resize(maxIntervalSize/intervalSize + 1);
				m_sumOfPos.resize(maxIntervalSize/intervalSize + 1);
				add(pos);
			}
			
		~KmerFeatures2(){
		}

		// increment frequence and compute approximate pos
		void add(long long int pos)
		{
			m_totalFreq++;
			m_totalSum+=pos;
			
			// compute sampled index
			int index = pos/m_intervalSize;
			// extension may exceed the expected length
           
			if(index < 0)
				index = 0;
			else if(index > (int) m_sumOfFreq.size()-1)
				index = m_sumOfFreq.size()-1;

			
			m_sumOfFreq[index]++;
			m_sumOfPos[index]+=pos;
		}
		
		long long int getTotalFreq()
		{
			return m_totalFreq;
		}

		long long int getTotalSum()
		{
			return m_totalSum;
		}

		// convert pos into index in the m_sumOfPos
		long long int getSumOfFreq(long long int pos)
		{
			int index = pos/m_intervalSize;
           
			long long int sumOfFreq = m_sumOfFreq[index];
			if(index > 0)
				sumOfFreq += m_sumOfFreq[index-1];
			if(index < (int) m_sumOfFreq.size()-1)
				sumOfFreq += m_sumOfFreq[index+1];
				
			return sumOfFreq;
		}
		
		// convert pos into index in the m_sumOfPos
		long long int getSumOfPos(long long int pos)
		{
			int index = pos/m_intervalSize;
			long long int sumOfPos = m_sumOfPos[index];
			if(index > 0)
				sumOfPos += m_sumOfPos[index-1];
			if(index < (int) m_sumOfPos.size()-1)
				sumOfPos += m_sumOfPos[index+1];
				
			return sumOfPos;
		}
		
		void printFreq()
		{
			for(size_t i=0; i<m_sumOfFreq.size(); i++)
			{
				std::cout << i <<":" << m_sumOfFreq[i] << "\t" << m_sumOfPos[i] <<"\n";
			}
		}
		
		inline void setVisited(bool visited)
		{
			m_isVisited = visited;
		}
		
		inline bool isVisited()
		{
			return m_isVisited;
		}
		
		inline void setMaxAvgFreq(double newAvgFreq)
		{
			m_maxAvgFreq = newAvgFreq;
		}
		
		inline double getMaxAvgFreq()
		{
			return m_maxAvgFreq;
		}
};

class HashtableSearch
{
    public: 
		HashtableSearch(const BWT* pBWT, const BWT* pRBWT);
		
        ~HashtableSearch();

		// Collect kmers from overlapping reads via right extension of single seed
		// void addHashFromSingleSeedUsingFMExtension(std::string& seedStr, size_t hashKmerSize, size_t maxLength);
		
		// LF-mapping of each SA index in the interval independently using loop instead of BFS tree expansion
		// Contaminated reads often are simple repeats C* or T* with large freq
		// Return false if the freq is way too large
		size_t addHashBySingleSeed(std::string& seedStr, size_t largeKmerSize, size_t smallKmerSize, size_t maxLength, bool skipRepeat, int expectedLength=-1,bool istarget = 0);


        size_t hashkmerfreqs(std::string fwdkmer,size_t kemrposition);
    private:

        //
        // Functions
        //
		// create root node using src
		
		void insertKmerToHash(std::string& insertedKmer, size_t seedStrLen, size_t currentLength, size_t smallKmerSize, size_t maxLength, int expectedLength);


		
       
        //
        // Data
        //
        
		int findBestExtSeq(SAIntervalNodeResultVector& results, const std::string &src, const std::string &dest, std::string &mergedseq,size_t hashKmerSize, size_t expectedLength, size_t maxLength, size_t maxLeaves);
		
        const BWT* m_pBWT;
        const BWT* m_pRBWT;

        size_t m_min_SA_threshold;
		size_t m_maxLeavesAllowed;
		
		int m_expectedLength;
		int m_currentLength;
		int m_seedLength;
		
		std::vector<END2> ENDs;

        SAIntervalNode* m_pRootNode;
        STNodePtrList m_leaves;

        BWTInterval m_fwdTerminatedInterval;   //in rBWT
        BWTInterval m_rvcTerminatedInterval;   //in BWT
		
		// multiple means for repeats leading to multinomial occurrences
		DenseHashMap<std::string, KmerFeatures2*, StringHasher> kmerHash;
		typedef DenseHashMap<std::string, KmerFeatures2*> :: iterator kmerHashiter;
		
		bool m_isSourceRepeat;
		bool m_isTargetRepeat;
		size_t m_sourceKmerSize;
		size_t m_targetKmerSize;
		
		bool m_isLargeLeaveRemoved;
		size_t m_prevExtPos;
		// bool debug;
};

#endif
