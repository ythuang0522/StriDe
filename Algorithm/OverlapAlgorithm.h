//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// OverlapAlgorithm - This class implements all the logic
// for finding and outputting overlaps for sequence reads
//
#ifndef OVERLAPALGORITHM_H
#define OVERLAPALGORITHM_H

#include "BWT.h"
#include "OverlapBlock.h"
#include "SearchSeed.h"
#include "BWTAlgorithms.h"
#include "Util.h"

enum OverlapMode
{
	OM_OVERLAP,
	OM_FULLREAD
};

struct OverlapResult
{
	OverlapResult() : isSubstring(false), searchAborted(false) {}
	bool isSubstring;
	bool searchAborted;
};

class OverlapAlgorithm
{

	// a wrapper of BWTIntervalPair with additional overlap/error information
	struct BWTOverlapInfo{
		#define m_localRange 11
		BWTOverlapInfo(){
			overlapLength=0;
			mismatch=0;
			insertion=0;
			deletion=0;

			// number of local errors
			m_localErrors=0;

			head=0;
			tail=m_localRange-1;
			memset(m_history,0,sizeof(int)*m_localRange);
			m_localDeletion=m_localInsertion=m_lastInsertion=0;
		};
		
		BWTIntervalPair pair;
		int overlapLength;
		int mismatch;
		int insertion;
		int deletion;
					
		// offset to dignal in banded DP
		int diagonalOffset;
		
		double getErrorRate()
		{
			if(overlapLength>0)
			return (double)(getTotalErrors())/overlapLength;
			else
			return 0;
		}		
		
		int getTotalErrors()
		{
			return mismatch+insertion+deletion;
		}
		
		int getLocalRange()
		{
			return m_localRange;
		}

		double getLocalErrorRate()
		{
			return (double) m_localErrors/m_localRange;
		}

		int getLocalErrors()
		{
			return m_localErrors;
		}
		//update local match or mismatch/indels after a mismatch/indel
		// void updateLocalError(int error)
		// {
			// assert(error >= 0);
			// m_localErrors -= dequeue();
			// m_localErrors += error;
			// enqueue(error);
		// }

		void updateLocalError(int error)
		{
			assert(error >= 0);
			int pop=dequeue();
			if(pop==2)
				m_localInsertion--;
			else if(pop==3)
				m_localDeletion--;
			
			m_localErrors -= pop>0?1:0;
			m_localErrors+=error;
			enqueue(error);
		}

		void updateLocalInsertion(int error)
		{
			assert(error >= 0);
			int pop=dequeue();
			if(pop==2)
				m_localInsertion--;
			else if(pop==3)
				m_localDeletion--;
			m_localErrors -= pop>0?1:0;
			m_localErrors++;
			m_localInsertion++;
			m_lastInsertion=error;
			enqueue(2);
		}
		
		void updateLocalDeletion(int error)
		{
			assert(error >= 0);
			int pop=dequeue();
			if(pop==2)
				m_localInsertion--;
			else if(pop==3)
				m_localDeletion--;

			m_localErrors -= pop>0?1:0;
			m_localErrors++;
			m_localDeletion++;
			enqueue(3);
		}
		
		inline bool isLocalIndel()
		{
			return m_localDeletion>0 || m_localInsertion>0;
		}

		inline int getLastInsertion()
		{
			return m_lastInsertion;
		}
		
		private:
			// int m_localRange;
			int m_localErrors;
			int m_localInsertion;
			int m_localDeletion;
			int m_lastInsertion;
			
			
			//circular queue implementation for storing error history
			int head;
			int tail;
			int m_history[m_localRange];
			void enqueue(int data){
				m_history[tail] = data;
				tail = (tail+1)%(m_localRange);
			}
			 
			 
			int dequeue(){
				int temp;
				temp = m_history[head];
				head = (head+1)%(m_localRange);
				return temp;
			}
	};

public:
	
	//exact overlap constructor
	OverlapAlgorithm(const BWT* pBWT, const BWT* pRevBWT,
	SuffixArray* pFwdSAI, SuffixArray* pRevSAI,
	ReadInfoTable* pQueryRIT, ReadInfoTable* pTargetRIT) : 
	m_pBWT(pBWT), 
	m_pRevBWT(pRevBWT),
	m_pFwdSAI(pFwdSAI),
	m_pRevSAI(pRevSAI),
	m_pQueryRIT(pQueryRIT),
	m_pTargetRIT(pTargetRIT),
	// default transitive reduction and exact overlap
	m_bIrreducible(true),
	m_exactModeOverlap(true),
	m_exactModeIrreducible(true)
	{

	}
	~OverlapAlgorithm()
	{
	}
	
	// Inexact overlap constructor
	OverlapAlgorithm(const BWT* pBWT, const BWT* pRevBWT,
	SuffixArray* pFwdSAI, SuffixArray* pRevSAI,
	ReadInfoTable* pQueryRIT, ReadInfoTable* pTargetRIT, double errorRate, int maxIndels, std::string algorithm) : 
	m_pBWT(pBWT), 
	m_pRevBWT(pRevBWT),
	m_pFwdSAI(pFwdSAI),
	m_pRevSAI(pRevSAI),
	m_pQueryRIT(pQueryRIT),
	m_pTargetRIT(pTargetRIT),                                         
	m_errorRate(errorRate),
	m_maxIndels(maxIndels),
	m_algorithm(algorithm),
	m_bIrreducible(false),
	m_exactModeOverlap(false),
	m_exactModeIrreducible(false)
	{

	}

	OverlapAlgorithm(const BWT* pBWT, const BWT* pRevBWT,  
	double er, int seedLen, int seedStride, bool irrOnly, int maxSeeds = -1) : 
	m_pBWT(pBWT), 
	m_pRevBWT(pRevBWT),
	m_errorRate(er),
	m_seedLength(seedLen),
	m_seedStride(seedStride),
	m_bIrreducible(irrOnly),
	m_exactModeOverlap(false),
	m_exactModeIrreducible(false),
	m_maxSeeds(maxSeeds)
	{

	}
	
	// Perform the overlap
	// This function is threaded so everything must be const
	OverlapResult overlapRead(const SeqRecord& read, int minOverlap, OverlapBlockList* pOutList) const;
	
	// Perform an irreducible overlap
	OverlapResult overlapReadExact(const SeqRecord& read, int minOverlap, OverlapBlockList* pOBOut) const;

	// Perform an inexact overlap using FM-index walk
    OverlapResult overlapReadInexactFMWalk(const SeqRecord& read, int minOverlap, OverlapBlockList* pOBOut) const;

	// Find duplicate blocks for this read
	OverlapResult alignReadDuplicate(const SeqRecord& read, OverlapBlockList* pOBOut) const;

	// Perform an inexact overlap
	OverlapResult overlapReadInexact(const SeqRecord& read, int minOverlap, OverlapBlockList* pOBOut) const;

	// Write the result of an overlap to an ASQG file
	void writeResultASQG(std::ostream& writer, const SeqRecord& read, const OverlapResult& result) const;

	// Write all the overlap blocks pList to the filehandle
	void writeOverlapBlocks(std::ostream& writer, size_t readIdx, bool isSubstring, const OverlapBlockList* pList) const;

	// Build the forward history structures for the blocks
	void buildForwardHistory(OverlapBlockList* pList) const;

	// Set flag to use exact-match algorithms only
	void setExactModeOverlap(bool b) { m_exactModeOverlap = b; }
	void setExactModeIrreducible(bool b) { m_exactModeIrreducible = b; }

	//
	const BWT* getBWT() const { return m_pBWT; }
	const BWT* getRBWT() const { return m_pRevBWT; }
	const SuffixArray* getFwdSAI() const { return m_pFwdSAI; }
	const SuffixArray* getRevSAI() const { return m_pRevSAI; }
	const ReadInfoTable* getQueryRIT() const { return m_pQueryRIT; }
	const ReadInfoTable* getTargetRIT() const { return m_pTargetRIT; }
	
	
private:

	// Calculate the ranges in pBWT that contain a prefix of at least minOverlap basepairs that
	// overlaps with a suffix of w.
	void findOverlapBlocksExact(const std::string& w, const BWT* pBWT, const BWT* pRevBWT, 
								const AlignFlags& af, const int minOverlap, OverlapBlockList* pOverlapList, 
								OverlapBlockList* pContainList, OverlapResult& result) const;

	// Implement inexact overlap using approximate banded DP via FM-index walk
	bool findOverlapBlocksInexact(const std::string& w, const BWT* pBWT, const BWT* pRevBWT, 
								const AlignFlags& af, const int minOverlap, OverlapBlockList* pOverlapList, 
								OverlapBlockList* pContainList, OverlapResult& result) const;

	// Implement inexact overlap using locality-sensitive backward search via FM-index walk
	bool findOverlapBlocksInexactFMIndexWalk(const std::string& w, const BWT* pBWT, const BWT* pRevBWT, 
								const AlignFlags& af, const int minOverlap, OverlapBlockList* pOverlapList, 
								OverlapBlockList* pContainList, OverlapResult& result) const;

	bool TrimOBLInterval(OverlapBlockList* pOverlapList, int MaxInterval) const;

	// Initialize the OverlapInfoList wrt char w
	void initOverlapInfoList(std::vector<BWTOverlapInfo>& BWTOverlapInfoVec, const std::string& w, size_t idx, const BWT* pBWT, const BWT* pRevBWT) const;

	// Expand OverlapInfoList wrt new char w, the error rate should be under control
	void expandOverlapInfoList(BWTOverlapInfo& currentOverlap, std::vector<BWTOverlapInfo>& BWTExpanedVector, const std::string& w, size_t idx, const BWT* pBWT) const;

	void expandOverlapInfoListByMismatch(BWTOverlapInfo& currOverlapInfo, std::vector<BWTOverlapInfo>& BWTExpanedVector, const std::string& w, size_t idx, const BWT* pBWT) const;

	// Expand OverlapInfoList wrt w with insertion
	void expandOverlapInfoListByInsertion(BWTOverlapInfo& currOverlapInfo, std::vector<BWTOverlapInfo>& BWTExpanedVector, const std::string& w, size_t idx, const BWT* pBWT) const;

	void expandOverlapInfoListByDeletion(BWTOverlapInfo& currOverlapInfo, std::vector<BWTOverlapInfo>& BWTExpanedVector, const std::string& w, size_t idx, const BWT* pBWT) const;

	// Calculate the terminal extension for the contained blocks to make the intervals consistent
	void terminateContainedBlocks(const std::string& w, const AlignFlags& af, BWTOverlapInfo& currentOverlap, const BWT* pBWT,
							const BWT* pRevBWT, OverlapBlockList* pContainList, OverlapResult& result) const;

	// Calculate the terminal extension for any overlap block
	void terminateOverlapBlocks(const AlignFlags& af, BWTOverlapInfo& currentOverlap, const BWT* pBWT, OverlapBlockList* pOverlapList) const;
							
	//                    
	// Irreducible-only processing algorithms
	//
	// Reduce the block list pOBList by removing blocks that correspond to transitive edges
	void computeIrreducibleBlocks(const BWT* pBWT, const BWT* pRevBWT, 
	OverlapBlockList* pOBList, OverlapBlockList* pOBFinal) const;
	
	// these recursive functions do the actual work of computing the irreducible blocks
	void _processIrreducibleBlocksExact(const BWT* pBWT, const BWT* pRevBWT, 
	OverlapBlockList& obList, OverlapBlockList* pOBFinal) const;


	void _processIrreducibleBlocksExactIterative(const BWT* pBWT, 
	const BWT* pRevBWT, 
	OverlapBlockList& inList, 
	OverlapBlockList* pOBFinal) const;
	//
	void _processIrreducibleBlocksInexact(const BWT* pBWT, const BWT* pRevBWT, 
	OverlapBlockList& obList, OverlapBlockList* pOBFinal) const;

	// Update the overlap block list with a righthand extension to b, removing ranges that become invalid
	void updateOverlapBlockRangesRight(const BWT* pBWT, const BWT* pRevBWT, 
	OverlapBlockList& obList, char b) const;
	
	//                                  
	void extendActiveBlocksRight(const BWT* pBWT, const BWT* pRevBWT, 
	OverlapBlockList& activeList, 
	OverlapBlockList& terminalList,
	OverlapBlockList& containedList) const;

	double calculateBlockErrorRate(const OverlapBlock& terminalBlock, const OverlapBlock& otherBlock) const;
	bool isBlockSubstring(OverlapBlock& terminalBlock, const OverlapBlockList& blockList, double maxER) const;

	// Data
	const BWT* m_pBWT;
	const BWT* m_pRevBWT;
	//
	
	//Direct ASQG Write
	SuffixArray* m_pFwdSAI; 
	SuffixArray* m_pRevSAI;
	ReadInfoTable* m_pQueryRIT;
	ReadInfoTable* m_pTargetRIT;

	// DenseHashSet<std::string,StringHasher>* m_SuperRepeat;

	double m_errorRate;
	size_t m_maxIndels;
	
	std::string m_algorithm;

	int m_seedLength;
	int m_seedStride;
	bool m_bIrreducible;
	bool m_exactModeOverlap;
	bool m_exactModeIrreducible;
	
	// Optional parameter to limit the amount of branching that is performed
	int m_maxSeeds; 
	
};

#endif
