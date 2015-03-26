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
    public:
		
		//exact overlap constructor
        OverlapAlgorithm(const BWT* pBWT, const BWT* pRevBWT,
                         SuffixArray* pFwdSAI, SuffixArray* pRevSAI,
						 ReadInfoTable* pQueryRIT,ReadInfoTable* pTargetRIT) : 
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
										 	//Maintain a hashset for reducing edges of repeat vertices,by YTH
											m_SuperRepeat=new DenseHashSet<std::string,StringHasher>();
											m_SuperRepeat->set_empty_key("");
										}
		~OverlapAlgorithm()
		{
			delete m_SuperRepeat;
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
										 	//Maintain a hashset for reducing edges of repeat vertices,by YTH
											m_SuperRepeat=new DenseHashSet<std::string,StringHasher>();
											m_SuperRepeat->set_empty_key("");
										}

		
        // Perform the overlap
        // This function is threaded so everything must be const
        OverlapResult overlapRead(const SeqRecord& read, int minOverlap, OverlapBlockList* pOutList) const;
    
        // Perform an irreducible overlap
        OverlapResult overlapReadExact(const SeqRecord& read, int minOverlap, OverlapBlockList* pOBOut) const;

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
		bool isSuperRepeatVertex(std::string vstr) const { 
			DenseHashSet<std::string,StringHasher>::const_iterator found;
			#pragma omp critical (SuperRepeat)
			{
				found  = m_SuperRepeat->find(vstr);
			}
			if(found!= m_SuperRepeat->end()) return true;
			else return false;
			
		}
		
    private:

        // Calculate the ranges in pBWT that contain a prefix of at least minOverlap basepairs that
        // overlaps with a suffix of w.
        void findOverlapBlocksExact(const std::string& w, const BWT* pBWT, const BWT* pRevBWT, 
                                    const AlignFlags& af, const int minOverlap, OverlapBlockList* pOBTemp, 
                                    OverlapBlockList* pOBFinal, OverlapResult& result) const;

        // Same as above while allowing mismatches
        bool findOverlapBlocksInexact(const std::string& w, const BWT* pBWT, const BWT* pRevBWT, 
                                      const AlignFlags& af, const int minOverlap, OverlapBlockList* pOBList, 
                                      OverlapBlockList* pOBFinal, OverlapResult& result) const;

		bool TrimOBLInterval(OverlapBlockList* pOverlapList, int MaxInterval) const;

        //
        inline bool extendSeedExactRight(SearchSeed& seed, const std::string& w, const BWT* pBWT, const BWT* pRevBWT) const;
        inline bool extendSeedExactLeft(SearchSeed& seed, const std::string& w, const BWT* pBWT, const BWT* pRevBWT) const;

        inline void branchSeedRight(const SearchSeed& seed, const std::string& w, const BWT* pBWT, const BWT* pRevBWT, SearchSeedQueue* pQueue) const;
        inline void branchSeedLeft(const SearchSeed& seed, const std::string& w, const BWT* pBWT, const BWT* pRevBWT, SearchSeedQueue* pQueue) const;

        //
        inline void extendSeedInexactRight(SearchSeed& seed, const std::string& w, const BWT* pBWT, const BWT* pRevBWT, 
                                           SearchSeedVector* pOutVector) const;

        //
        inline void extendSeedInexactLeft(SearchSeed& seed, const std::string& w, const BWT* pBWT, const BWT* pRevBWT,
                                          SearchSeedVector* pOutVector) const;


        //
        inline void calculateSeedParameters(const std::string& w, const int minOverlap, int& seed_length, int& seed_stride) const;
        
        //
        inline int createSearchSeeds(const std::string& w, const BWT* pBWT, 
                                     const BWT* pRevBWT, int seed_length, int seed_stride, 
                                     SearchSeedVector* pOutVector) const;

        //
        inline void extendSeedsExactRightQueue(const std::string& w, const BWT* pBWT, const BWT* pRevBWT, 
                                                 ExtendDirection dir, const SearchSeedVector* pInVector, 
                                                 SearchSeedQueue* pOutQueue) const;

        //
        inline void extendSeedsExactRight(const std::string& w, const BWT* pBWT, const BWT* pRevBWT, 
                                                 ExtendDirection dir, const SearchSeedVector* pInVector, 
                                                 SearchSeedVector* pOutVector) const;
        
        // Calculate the terminal extension for the contained blocks to make the intervals consistent
        void terminateContainedBlocks(OverlapBlockList& containedBlocks) const;

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

		DenseHashSet<std::string,StringHasher>* m_SuperRepeat;

		double m_errorRate;
        int m_seedLength;
        int m_seedStride;
        bool m_bIrreducible;
        bool m_exactModeOverlap;
        bool m_exactModeIrreducible;
        
        // Optional parameter to limit the amount of branching that is performed
        int m_maxSeeds; 
};

#endif
