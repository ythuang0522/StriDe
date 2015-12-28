//----------------------------------------------
// Copyright 2015 National Chung Cheng University
// Written by Yao-Ting Huang, revised from Simpson's overlap structure
// Released under the GPL
//-----------------------------------------------

#include "OverlapAlgorithm.h"
#include "ASQG.h"
#include <tr1/unordered_set>
#include <math.h>
#include <SAIOverlapTree.h>

// Collect the complete set of overlaps in pOBOut
static const AlignFlags sufPreAF(false, false, false);
static const AlignFlags prePreAF(false, true, true);
static const AlignFlags sufSufAF(true, false, true);
static const AlignFlags preSufAF(true, true, false);

//#define DEBUGOVERLAP 1

// Perform the overlap
OverlapResult OverlapAlgorithm::overlapRead(const SeqRecord& read, int minOverlap, OverlapBlockList* pOutList) const
{

	OverlapResult r;
	if(static_cast<int>(read.seq.length()) < minOverlap)
	return r;

	if(!m_exactModeOverlap)
	{
		if(m_algorithm=="LSSF")
			r = overlapReadInexactFMWalk(read, minOverlap, pOutList);
		else if(m_algorithm=="ADPF")
			r = overlapReadInexact(read, minOverlap, pOutList);
		else
		{
			std::cout << "Unknown algorithm!!\n";
			assert(false);
		}
	}
	else
		r = overlapReadExact(read, minOverlap, pOutList);
	return r;
}

//
OverlapResult OverlapAlgorithm::overlapReadInexact(const SeqRecord& read, int minOverlap, OverlapBlockList* pOBOut) const
{
	OverlapResult result;
	// The complete set of overlap blocks are collected in obWorkingList
	// The filtered set (containing only irreducible overlaps) are placed into pOBOut
	// by calculateIrreducibleHits
	OverlapBlockList obWorkingList;
	std::string seq = read.seq.toString();

	// We store the various overlap blocks using a number of lists, one for the containments
	// in the forward and reverse index and one for each set of overlap blocks
	OverlapBlockList oblFwdContain;
	OverlapBlockList oblRevContain;
	
	OverlapBlockList oblSuffixFwd;
	OverlapBlockList oblSuffixRev;
	OverlapBlockList oblPrefixFwd;
	OverlapBlockList oblPrefixRev;

	// std::cout << read.id << "\n"; 
	
	// Match the suffix of seq to prefixes
	findOverlapBlocksInexact(seq, m_pBWT, m_pRevBWT, sufPreAF, minOverlap, &oblSuffixFwd, &oblFwdContain, result);
	if(result.isSubstring) return result;
	
	// getchar();
	// std::cout << complement(seq) << "\n";
	findOverlapBlocksInexact(complement(seq), m_pRevBWT, m_pBWT, prePreAF, minOverlap, &oblSuffixRev, &oblRevContain, result);
	if(result.isSubstring) return result;
	
	// Match the prefix of seq to suffixes
	// getchar();
	// std::cout << reverseComplement(seq) << "\n";
	findOverlapBlocksInexact(reverseComplement(seq), m_pBWT, m_pRevBWT, sufSufAF, minOverlap, &oblPrefixFwd, &oblFwdContain, result);
	if(result.isSubstring) return result;
	
	// getchar();
	// std::cout << reverse(seq) << "\n";
	findOverlapBlocksInexact(reverse(seq), m_pRevBWT, m_pBWT, preSufAF, minOverlap, &oblPrefixRev, &oblRevContain, result);
	if(result.isSubstring) return result;
	
	// Remove submaximal blocks for each block list including fully contained blocks
	// Copy the containment blocks into the prefix/suffix lists
	oblSuffixFwd.insert(oblSuffixFwd.end(), oblFwdContain.begin(), oblFwdContain.end());
	oblPrefixFwd.insert(oblPrefixFwd.end(), oblFwdContain.begin(), oblFwdContain.end());
	oblSuffixRev.insert(oblSuffixRev.end(), oblRevContain.begin(), oblRevContain.end());
	oblPrefixRev.insert(oblPrefixRev.end(), oblRevContain.begin(), oblRevContain.end());

	//Trim the OB list
	TrimOBLInterval(&oblSuffixFwd, seq.length());
	TrimOBLInterval(&oblSuffixRev, seq.length());
	TrimOBLInterval(&oblPrefixFwd, seq.length());
	TrimOBLInterval(&oblPrefixRev, seq.length());

	// Perform the submaximal filter, 
	// bug: Error in resolveOverlap: Overlap blocks with same length do not the have same coordinates
	// std::cout << oblSuffixFwd.size() << "\n";
	removeSubMaximalBlocks(&oblSuffixFwd, m_pBWT, m_pRevBWT);
	if(containSubstringBlocks(&oblSuffixFwd, seq.length()))
		result.isSubstring=true;
		
	removeSubMaximalBlocks(&oblPrefixFwd, m_pBWT, m_pRevBWT);
	if(containSubstringBlocks(&oblPrefixFwd, seq.length()))
		result.isSubstring=true;
		
	removeSubMaximalBlocks(&oblSuffixRev, m_pRevBWT, m_pBWT);
	if(containSubstringBlocks(&oblSuffixRev, seq.length()))
		result.isSubstring=true;
		
	removeSubMaximalBlocks(&oblPrefixRev, m_pRevBWT, m_pBWT);
	if(containSubstringBlocks(&oblPrefixRev, seq.length()))
		result.isSubstring=true;

	if(result.isSubstring) return result;
	
	// Remove the contain blocks from the suffix/prefix lists for transitive reduction
	// However, transitive reduction algorithm can't be applied for indel overlap
	removeContainmentBlocks(seq.length(), &oblSuffixFwd);
	removeContainmentBlocks(seq.length(), &oblPrefixFwd);
	removeContainmentBlocks(seq.length(), &oblSuffixRev);
	removeContainmentBlocks(seq.length(), &oblPrefixRev);

	
	// Move the containments to the output list first
	pOBOut->splice(pOBOut->end(), oblFwdContain);
	pOBOut->splice(pOBOut->end(), oblRevContain);

	// Filter out transitive overlap blocks if requested
	if(m_bIrreducible)
	{
		computeIrreducibleBlocks(m_pBWT, m_pRevBWT, &oblSuffixFwd, pOBOut);
		computeIrreducibleBlocks(m_pBWT, m_pRevBWT, &oblPrefixFwd, pOBOut);
	}
	else
	{
		pOBOut->splice(pOBOut->end(), oblSuffixFwd);
		pOBOut->splice(pOBOut->end(), oblPrefixFwd);
	}

	return result;
}

OverlapResult OverlapAlgorithm::overlapReadInexactFMWalk(const SeqRecord& read, int minOverlap, OverlapBlockList* pOBOut) const
{
	OverlapResult result;
	// The complete set of overlap blocks are collected in obWorkingList
	// The filtered set (containing only irreducible overlaps) are placed into pOBOut
	// by calculateIrreducibleHits
	OverlapBlockList obWorkingList;
	std::string seq = read.seq.toString();

	// We store the various overlap blocks using a number of lists, one for the containments
	// in the forward and reverse index and one for each set of overlap blocks
	OverlapBlockList oblFwdContain;
	OverlapBlockList oblRevContain;
	
	OverlapBlockList oblSuffixFwd;
	OverlapBlockList oblSuffixRev;
	OverlapBlockList oblPrefixFwd;
	OverlapBlockList oblPrefixRev;

	// std::cout << read.id << "\n"; 
	
	// skip short reads
	if(seq.length() < (size_t) minOverlap)
		return result;
	
	// Match the suffix of seq to prefixes
	findOverlapBlocksInexactFMIndexWalk(seq, m_pBWT, m_pRevBWT, sufPreAF, minOverlap, &oblSuffixFwd, &oblFwdContain, result);
	if(result.isSubstring) return result;
	
	// getchar();
	// std::cout << complement(seq) << "\n";
	findOverlapBlocksInexactFMIndexWalk(complement(seq), m_pRevBWT, m_pBWT, prePreAF, minOverlap, &oblSuffixRev, &oblRevContain, result);
	if(result.isSubstring) return result;
	
	// Match the prefix of seq to suffixes
	// getchar();
	// std::cout << reverseComplement(seq) << "\n";
	findOverlapBlocksInexactFMIndexWalk(reverseComplement(seq), m_pBWT, m_pRevBWT, sufSufAF, minOverlap, &oblPrefixFwd, &oblFwdContain, result);
	if(result.isSubstring) return result;
	
	// getchar();
	// std::cout << reverse(seq) << "\n";
	findOverlapBlocksInexactFMIndexWalk(reverse(seq), m_pRevBWT, m_pBWT, preSufAF, minOverlap, &oblPrefixRev, &oblRevContain, result);
	if(result.isSubstring) return result;
	
	// Remove submaximal blocks for each block list including fully contained blocks
	// Copy the containment blocks into the prefix/suffix lists
	oblSuffixFwd.insert(oblSuffixFwd.end(), oblFwdContain.begin(), oblFwdContain.end());
	oblPrefixFwd.insert(oblPrefixFwd.end(), oblFwdContain.begin(), oblFwdContain.end());
	oblSuffixRev.insert(oblSuffixRev.end(), oblRevContain.begin(), oblRevContain.end());
	oblPrefixRev.insert(oblPrefixRev.end(), oblRevContain.begin(), oblRevContain.end());

	//Trim the OB list
	// TrimOBLInterval(&oblSuffixFwd, seq.length());
	// TrimOBLInterval(&oblSuffixRev, seq.length());
	// TrimOBLInterval(&oblPrefixFwd, seq.length());
	// TrimOBLInterval(&oblPrefixRev, seq.length());

	// Perform the submaximal filter, 
	removeSubMaximalBlocks(&oblSuffixFwd, m_pBWT, m_pRevBWT);
	if(containSubstringBlocks(&oblSuffixFwd, seq.length()))
		result.isSubstring=true;
		
	removeSubMaximalBlocks(&oblPrefixFwd, m_pBWT, m_pRevBWT);
	if(containSubstringBlocks(&oblPrefixFwd, seq.length()))
		result.isSubstring=true;
		
	removeSubMaximalBlocks(&oblSuffixRev, m_pRevBWT, m_pBWT);
	if(containSubstringBlocks(&oblSuffixRev, seq.length()))
		result.isSubstring=true;
		
	removeSubMaximalBlocks(&oblPrefixRev, m_pRevBWT, m_pBWT);
	if(containSubstringBlocks(&oblPrefixRev, seq.length()))
		result.isSubstring=true;

	if(result.isSubstring) return result;
	

	// Join the suffix and prefix lists
	oblSuffixFwd.splice(oblSuffixFwd.end(), oblSuffixRev);
	oblPrefixFwd.splice(oblPrefixFwd.end(), oblPrefixRev);
	oblPrefixFwd.splice(oblPrefixFwd.end(), oblSuffixFwd);
	
	pOBOut->splice(pOBOut->end(), oblPrefixFwd);

	return result;
}

//
OverlapResult OverlapAlgorithm::alignReadDuplicate(const SeqRecord& read, OverlapBlockList* pOBOut) const
{
	OverlapResult result;
	OverlapBlockList obWorkingList;
	std::string seq = read.seq.toString();
	int readLength = seq.length();

	findOverlapBlocksInexact(seq, m_pBWT, m_pRevBWT, sufPreAF, readLength, &obWorkingList, pOBOut, result);
	findOverlapBlocksInexact(complement(seq), m_pRevBWT, m_pBWT, prePreAF, readLength, &obWorkingList, pOBOut, result);
	return result;
}


// Construct the set of blocks describing irreducible overlaps with READ
// and write the blocks to pOBOut
OverlapResult OverlapAlgorithm::overlapReadExact(const SeqRecord& read, int minOverlap, OverlapBlockList* pOBOut) const
{	
	OverlapResult result;
	// The complete set of overlap blocks are collected in obWorkingList
	// The filtered set (containing only irreducible overlaps) are placed into pOBOut
	// by calculateIrreducibleHits
	OverlapBlockList obWorkingList;
	std::string seq = read.seq.toString();

	// We store the various overlap blocks using a number of lists, one for the containments
	// in the forward and reverse index and one for each set of overlap blocks
	OverlapBlockList oblFwdContain;
	OverlapBlockList oblRevContain;
	
	OverlapBlockList oblSuffixFwd;
	OverlapBlockList oblSuffixRev;
	OverlapBlockList oblPrefixFwd;
	OverlapBlockList oblPrefixRev;

	// Match the suffix of seq to prefixes
	findOverlapBlocksExact(seq, m_pBWT, m_pRevBWT, sufPreAF, minOverlap, &oblSuffixFwd, &oblFwdContain, result);
	findOverlapBlocksExact(complement(seq), m_pRevBWT, m_pBWT, prePreAF, minOverlap, &oblSuffixRev, &oblRevContain, result);

	// Match the prefix of seq to suffixes
	findOverlapBlocksExact(reverseComplement(seq), m_pBWT, m_pRevBWT, sufSufAF, minOverlap, &oblPrefixFwd, &oblFwdContain, result);
	findOverlapBlocksExact(reverse(seq), m_pRevBWT, m_pBWT, preSufAF, minOverlap, &oblPrefixRev, &oblRevContain, result);

	//Trim the OB list
	TrimOBLInterval(&oblSuffixFwd, seq.length());
	TrimOBLInterval(&oblSuffixRev, seq.length());
	TrimOBLInterval(&oblPrefixFwd, seq.length());
	TrimOBLInterval(&oblPrefixRev, seq.length());


	// Remove submaximal blocks for each block list including fully contained blocks
	// Copy the containment blocks into the prefix/suffix lists
	oblSuffixFwd.insert(oblSuffixFwd.end(), oblFwdContain.begin(), oblFwdContain.end());
	oblPrefixFwd.insert(oblPrefixFwd.end(), oblFwdContain.begin(), oblFwdContain.end());
	oblSuffixRev.insert(oblSuffixRev.end(), oblRevContain.begin(), oblRevContain.end());
	oblPrefixRev.insert(oblPrefixRev.end(), oblRevContain.begin(), oblRevContain.end());
	
	// Perform the submaximal filter
	removeSubMaximalBlocks(&oblSuffixFwd, m_pBWT, m_pRevBWT);
	removeSubMaximalBlocks(&oblPrefixFwd, m_pBWT, m_pRevBWT);
	removeSubMaximalBlocks(&oblSuffixRev, m_pRevBWT, m_pBWT);
	removeSubMaximalBlocks(&oblPrefixRev, m_pRevBWT, m_pBWT);
	
	// Remove the contain blocks from the suffix/prefix lists
	removeContainmentBlocks(seq.length(), &oblSuffixFwd);
	removeContainmentBlocks(seq.length(), &oblPrefixFwd);
	removeContainmentBlocks(seq.length(), &oblSuffixRev);
	removeContainmentBlocks(seq.length(), &oblPrefixRev);

	// Join the suffix and prefix lists
	oblSuffixFwd.splice(oblSuffixFwd.end(), oblSuffixRev);
	oblPrefixFwd.splice(oblPrefixFwd.end(), oblPrefixRev);

	// Move the containments to the output list
	pOBOut->splice(pOBOut->end(), oblFwdContain);
	pOBOut->splice(pOBOut->end(), oblRevContain);

	// Filter out transitive overlap blocks if requested
	if(m_bIrreducible)
	{
		computeIrreducibleBlocks(m_pBWT, m_pRevBWT, &oblSuffixFwd, pOBOut);
		computeIrreducibleBlocks(m_pBWT, m_pRevBWT, &oblPrefixFwd, pOBOut);
	}
	else
	{
		pOBOut->splice(pOBOut->end(), oblSuffixFwd);
		pOBOut->splice(pOBOut->end(), oblPrefixFwd);
	}

	return result;
}

// Limit OBList interval By Ya: 20141022
// Only retain edges from longest overlap block to short ones
bool OverlapAlgorithm::TrimOBLInterval(OverlapBlockList* pOverlapList, int readLength) const
{
	// The overlapLen of insertions are not in the right order.
	pOverlapList->sort(OverlapBlock::sortSizeDescending);

	int Interval = 0;
	bool isSuperRepeat=false;
	std::list <OverlapBlock>::iterator OB = pOverlapList->end();
	
	OB--;	//Get the longest Overlap Block
	int longestOverlap=OB->getOverlapLength();

	//Prune the list from longest overlap to shortest ones
	for(; OB != pOverlapList->begin() ; OB--)
	{
		assert(OB->ranges.interval[1].isValid());
		Interval += OB->ranges.interval[1].size();
		
		// readLength can be from kmer size to insert size
		// For kmerized reads, if it's high-error, interval size should be small
		// For insert-sized long reads, no need to consider overlap diff too large
		if(Interval >= 128 || (longestOverlap - OB->getOverlapLength()>=readLength*0.5) ) 
		{
			// std::cout << readLength <<  "\t" << Interval << "\t" << longestOverlap <<"\n";
			for(;; OB--)
			{
				if(OB != pOverlapList->begin())
				pOverlapList->erase(OB++);
				else	//first OB, erase without movement
				{
					pOverlapList->erase(OB);
					break;
				}
			}
			break;
		}
	}
	
	//std::cout << "Interval size: " << Interval << "\n";
	
	return isSuperRepeat;
}

// Write overlap results to an ASQG file
void OverlapAlgorithm::writeResultASQG(std::ostream& writer, const SeqRecord& read, const OverlapResult& result) const
{
	ASQG::VertexRecord record(read.id, read.seq.toString());
	record.setSubstringTag(result.isSubstring);
	record.write(writer);
}

// Write overlap blocks out to a file
void OverlapAlgorithm::writeOverlapBlocks(std::ostream& writer, size_t readIdx, bool isSubstring, const OverlapBlockList* pList) const
{
	// Write the header info
	size_t numBlocks = pList->size();
	writer << readIdx << " " << isSubstring << " " << numBlocks << " ";
	//std::cout  << readIdx << " " << isSubstring << " " << numBlocks << std::endl;
	//std::cout << "<Wrote> idx: " << count << " count: " << numBlocks << "\n";
	for(OverlapBlockList::const_iterator iter = pList->begin(); iter != pList->end(); ++iter)
	{
		
		//OverlapBlock ob = *iter;
		//std::cout << ob << std::endl;
		writer << *iter << " ";
	}
	writer << "\n";
}

// Calculate the ranges in pBWT that contain a prefix of at least minOverlap basepairs that
// overlaps with a suffix of w. The ranges are added to the pOBList
void OverlapAlgorithm::findOverlapBlocksExact(const std::string& w, const BWT* pBWT,
const BWT* pRevBWT, const AlignFlags& af, int minOverlap,
OverlapBlockList* pOverlapList, OverlapBlockList* pContainList, 
OverlapResult& result) const
{
	// The algorithm is as follows:
	// We perform a backwards search using the FM-index for the string w.
	// As we perform the search we collect the intervals 
	// of the significant prefixes (len >= minOverlap) that overlap w.
	BWTIntervalPair ranges;
	size_t l = w.length();
	int start = l - 1;
	BWTAlgorithms::initIntervalPair(ranges, w[start], pBWT, pRevBWT);
	
	// Collect the OverlapBlocks
	for(size_t i = start - 1; i >= 1; --i)
	{
		// Compute the range of the suffix w[i, l]
		BWTAlgorithms::updateBothL(ranges, w[i], pBWT);
		int overlapLen = l - i;
		if(overlapLen >= minOverlap)
		{
			// Calculate which of the prefixes that match w[i, l] are terminal
			// These are the proper prefixes (they are the start of a read)
			BWTIntervalPair probe = ranges;
			BWTAlgorithms::updateBothL(probe, '$', pBWT);
			
			// The probe interval contains the range of proper prefixes
			if(probe.interval[1].isValid())
			{
				assert(probe.interval[1].lower > 0);
				pOverlapList->push_back(OverlapBlock(probe, ranges, overlapLen, 0, af));
			}
		}
	}

	// Determine if this sequence is contained and should not be processed further
	BWTAlgorithms::updateBothL(ranges, w[0], pBWT);

	// Ranges now holds the interval for the full-length read
	// To handle containments, we output the overlapBlock to the final overlap block list
	// and it will be processed later
	// Two possible containment cases:
	// 1) This read is a substring of some other read
	// 2) This read is identical to some other read
	
	// Case 1 is indicated by the existance of a non-$ left or right hand extension
	// In this case we return no alignments for the string
	AlphaCount64 left_ext = BWTAlgorithms::getExtCount(ranges.interval[0], pBWT);
	AlphaCount64 right_ext = BWTAlgorithms::getExtCount(ranges.interval[1], pRevBWT);
	if(left_ext.hasDNAChar() || right_ext.hasDNAChar())
	{
		result.isSubstring = true;
	}
	else
	{
		BWTIntervalPair probe = ranges;
		BWTAlgorithms::updateBothL(probe, '$', pBWT);
		if(probe.isValid())
		{
			// terminate the contained block and add it to the contained list
			BWTAlgorithms::updateBothR(probe, '$', pRevBWT);
			assert(probe.isValid());
			pContainList->push_back(OverlapBlock(probe, ranges, w.length(), 0, af));
		}
	}

	return;
}

// Simulate Banded-DP overlap detection using FM-index walk
// Each walk corresponds to one possible DP path in the banded DP.
// The walks are expanded or reduced wrt error rate
// The starting and ending overlap regions are required to be highly accurate.
bool OverlapAlgorithm::findOverlapBlocksInexact(const std::string& w, const BWT* pBWT, 
const BWT* pRevBWT, const AlignFlags& af, int minOverlap,
OverlapBlockList* pOverlapList, OverlapBlockList* pContainList, 
OverlapResult& result) const
{
	size_t l = w.length();
	int start = l - 1;
		
	//create a vector of BWTOverlapInfo
	std::vector<BWTOverlapInfo> BWTOverlapInfoVec;
	
	// Initially, assume the 1st base may contain error and thus consider all 4 possible bases
	initOverlapInfoList(BWTOverlapInfoVec, w, start, pBWT, pRevBWT);
	
	// Collect SA intervals overlapping with w
	// Expand the list of intervals if error rate is still acceptable
	for(size_t i = start - 1; i >= 1; --i)
	{
		int overlapLen = l - i;

		assert(BWTOverlapInfoVec.size()>0);
		
		/*** debugging code ***/
		// std::cout << i << "\t"<< w[i] << "\t" << BWTOverlapInfoVec.size() << "\t" 
		// << BWTOverlapInfoVec.at(0).pair.interval[0].size() << "\t"<< overlapLen <<"\n";

		// for(i < 8952, size_t j=0; j<BWTOverlapInfoVec.size(); j++)
			// std::cout << BWTOverlapInfoVec.at(j).getTotalErrors() << "\t" << BWTOverlapInfoVec.at(j).mismatch
				// <<"\t" <<BWTOverlapInfoVec.at(j).insertion << "\t"<< BWTOverlapInfoVec.at(j).deletion << "\t"<< 
				// BWTOverlapInfoVec.at(j).getLocalErrors() << "\t"<< BWTOverlapInfoVec.at(j).pair.interval[0] <<"\n";
		// getchar();

		
		// Expand the list of SA intervals if w[i] is error or SNP or indel
		std::vector<BWTOverlapInfo> BWTExpanedVector;
		for(size_t idx=0; idx < BWTOverlapInfoVec.size(); idx++)
		{
			// The list may grow exponentially 
			// (1) leaving repeat overlap; 
			bool isTooManyIntervlas = BWTOverlapInfoVec.size()>128 ;
			bool isAnyLocalError = BWTOverlapInfoVec.at(idx).getLocalErrors()>0;

			// (2) false SNP/indels after long overlap leading to false intervals (e.g.,1000*0.05=50)
			// Prone the list according to the local error rate
			// 27      2049306 13321M  NM:i:30
			// 29      2050171 10239M  NM:i:9
			bool isLocalErrorRateNotAcceptable = BWTOverlapInfoVec.at(idx).getLocalErrorRate() > 0.5;
			
			if( (isTooManyIntervlas && isAnyLocalError) || isLocalErrorRateNotAcceptable) 
				continue;
			
			// (3) The interval is often a redundant insertion+deletion, not better than mismatch
			if(BWTOverlapInfoVec.at(idx).insertion>0 && BWTOverlapInfoVec.at(idx).deletion>0 && idx>0)
			{
				int prevDiagonal = BWTOverlapInfoVec.at(idx-1).diagonalOffset;
				int currDiagonal = BWTOverlapInfoVec.at(idx).diagonalOffset;
				if(prevDiagonal == currDiagonal+BWTOverlapInfoVec.at(idx).getLastInsertion() &&
					BWTOverlapInfoVec.at(idx-1).pair == BWTOverlapInfoVec.at(idx).pair && 
					BWTOverlapInfoVec.at(idx-1).getTotalErrors()<BWTOverlapInfoVec.at(idx).getTotalErrors() )
					{
					// std::cout << BWTOverlapInfoVec.at(idx-1).pair.interval[0] << "\t" << BWTOverlapInfoVec.at(idx-1).insertion << "\t"
								// << BWTOverlapInfoVec.at(idx-1).deletion << "\t" << BWTOverlapInfoVec.at(idx-1).mismatch << "\t"
								// << BWTOverlapInfoVec.at(idx).pair.interval[0] << "\t" << BWTOverlapInfoVec.at(idx).insertion<< "\t"
								// << BWTOverlapInfoVec.at(idx).deletion << "\t" << BWTOverlapInfoVec.at(idx).mismatch << "\n";
					}
				// getchar();
					// continue;
			}
			// insertion may reach the end of w earlier than others
			if((int)i+BWTOverlapInfoVec.at(idx).diagonalOffset == 0 
				&& BWTOverlapInfoVec.at(idx).deletion>BWTOverlapInfoVec.at(idx).insertion)
			{
				terminateContainedBlocks(w, af, BWTOverlapInfoVec.at(idx), pBWT, pRevBWT, pContainList, result);
				continue;
			}
			
			if((int)i+BWTOverlapInfoVec.at(idx).diagonalOffset < 1)
			{
				// std::cout << i << "\t" << BWTOverlapInfoVec.size() << "\t" 
				// << BWTOverlapInfoVec.at(i).pair.interval[0].size() << "\t"
				// << BWTOverlapInfoVec.at(i).insertion <<"\t"
				// << BWTOverlapInfoVec.at(i).deletion <<"\n";

				continue;
			}
			//Expand the interval wrt SNP or indel
			expandOverlapInfoList(BWTOverlapInfoVec.at(idx), BWTExpanedVector, w, i, pBWT);
		}
				
		if(BWTExpanedVector.empty())
			return true;
		
		// Push SA intervals reaching $ into output
		BWTOverlapInfoVec.clear();
		int m_minOverlap = (int)w.length()<minOverlap?w.length()*0.8:minOverlap;
		if(overlapLen >= m_minOverlap)
		// if(overlapLen >= minOverlap)
		{
			// check if any interval reaching $
			for(size_t idx=0; idx < BWTExpanedVector.size(); idx++)
			{
				// move the overlap block into pOverlapList if reaching terminal $
				terminateOverlapBlocks(af, BWTExpanedVector.at(idx), pBWT, pOverlapList);
			}
		}
		
		// Always copy each interval in BWTExpanedVector into BWTOverlapInfoVec
		BWTOverlapInfoVec = BWTExpanedVector;
		BWTExpanedVector.clear();
	}// end of BWT update using w[i]
	
	// Determine if this sequence is contained and should not be processed further
	for(size_t idx=0; idx < BWTOverlapInfoVec.size(); idx++)
	{
		terminateContainedBlocks(w, af, BWTOverlapInfoVec.at(idx), pBWT, pRevBWT, pContainList, result);
		if(result.isSubstring) return false;

	}
	return true;
}

// Initialize the OverlapInfoList wrt char w
void OverlapAlgorithm::initOverlapInfoList(std::vector<BWTOverlapInfo>& BWTOverlapInfoVec, const std::string& w, size_t idx, const BWT* pBWT, const BWT* pRevBWT) const
{
	// Initially, assume the 1st base may contain error and thus consider all 4 possible bases
	// w: AACCCTTTGGG
	//       CCTTTGGA
	for(size_t i = 0; i < 4; i++)
	{
		char b=ALPHABET[i];
		BWTOverlapInfo overlapinfo;
		overlapinfo.mismatch = (b==w[idx])?0:1;
		overlapinfo.overlapLength=1;
		overlapinfo.diagonalOffset=0;
		overlapinfo.updateLocalError(overlapinfo.mismatch);
		BWTAlgorithms::initIntervalPair(overlapinfo.pair, b, pBWT, pRevBWT);
		BWTOverlapInfoVec.push_back(overlapinfo);
	}

	// Initialize insertion candidates
	// w: AACCCTTTGGG
	//       CCTTTG--
	// for(size_t i = 1; i <= m_maxIndels; i++)
	// {
		// BWTOverlapInfo overlapinfo;
		// overlapinfo.insertion=i;
		// overlapinfo.overlapLength=1;
		// overlapinfo.diagonalOffset= -i;

		// // we do not consider expanding mismatches right after indels, which is too rare
		// if(idx>=i){
			// BWTAlgorithms::initIntervalPair(overlapinfo.pair, w[idx+overlapinfo.diagonalOffset], pBWT, pRevBWT);
			// BWTOverlapInfoVec.push_back(overlapinfo);
		// }
	// }
	
	// Deletion candidates do not need initialization as the initialization of 1st base has included
	// w: AACCCTTTGGG--
	//       CCTTTGGGAA
}

// Expand OverlapInfoList wrt new char w, the error rate should be under control
void OverlapAlgorithm::expandOverlapInfoList(BWTOverlapInfo& currOverlapInfo, std::vector<BWTOverlapInfo>& BWTExpanedVector, const std::string& w, size_t idx, const BWT* pBWT) const
{
	// At first, try extending to w only, don't expand if all reads in the list match w extension
	BWTIntervalPair probe = currOverlapInfo.pair;
	int64_t prevSAISize = probe.interval[0].size();

	// Compute new SA intervals using prefix b
	// if((int)idx+currOverlapInfo.diagonalOffset >= 1)
		// BWTAlgorithms::updateBothL(probe, w[idx + currOverlapInfo.diagonalOffset], pBWT);
	// else
	// {
		// //push back for later containment check
		// // BWTExpanedVector.push_back(currOverlapInfo);
		// // buggy, requiring other solutions
		// return;
	// }
	BWTAlgorithms::updateBothL(probe, w[idx + currOverlapInfo.diagonalOffset], pBWT);
	
	if(probe.isValid())
	{
		// Successful extension using b
		// Store the new BWTInterval into BWTExpanedVector
		BWTOverlapInfo newOverlapInfo = currOverlapInfo;
		newOverlapInfo.overlapLength = currOverlapInfo.overlapLength+1;
		newOverlapInfo.pair = probe;
		newOverlapInfo.updateLocalError(0);
		BWTExpanedVector.push_back(newOverlapInfo);
		// std::cout << b << "\t" << currInfo.mismatch << "\n";
		
		// Don't expand the list if w is a feasible extension for all reads in the list
		int64_t currentSize = probe.interval[0].size();
		if(currentSize == prevSAISize)
			return;
		
		assert(currentSize < prevSAISize);
		// Skip if the reduction is due to reaching read ends
		BWTIntervalPair endingProbe = currOverlapInfo.pair;
		BWTAlgorithms::updateBothL(endingProbe, '$', pBWT);
		if(endingProbe.isValid()) 
			currentSize+=endingProbe.interval[0].size();
			
		if(currentSize == prevSAISize)
			return;
	}

	// Otherwise, expand the list for mismatches, insertions, and deletions
	if(!currOverlapInfo.isLocalIndel())
		expandOverlapInfoListByMismatch(currOverlapInfo, BWTExpanedVector, w, idx, pBWT);
	
	if(currOverlapInfo.insertion < (int)m_maxIndels*2 && !currOverlapInfo.isLocalIndel())
		expandOverlapInfoListByInsertion(currOverlapInfo, BWTExpanedVector, w, idx, pBWT);
	
	if(currOverlapInfo.deletion < (int)m_maxIndels*2 && !currOverlapInfo.isLocalIndel())
		expandOverlapInfoListByDeletion(currOverlapInfo, BWTExpanedVector, w, idx, pBWT);

}

// Expand the list of SA intervals for error/SNP 
// w: AACCCTTTTGG
//       CCTTTTTG
//       CCTTTTAG 
void OverlapAlgorithm::expandOverlapInfoListByMismatch(BWTOverlapInfo& currOverlapInfo, std::vector<BWTOverlapInfo>& BWTExpanedVector, const std::string& w, size_t idx, const BWT* pBWT) const
{
	// (1) error rate is acceptable if overlap length is enough
	// (2) the overlap length is too short to reflect error rate, tolerate one error upto maximum possible length	
	size_t newTotalErrors = currOverlapInfo.getTotalErrors()+1;
	double newErrorRate = (double)newTotalErrors/(currOverlapInfo.overlapLength+1);
	
	// int maxOneErrorLength = m_errorRate>0?(double)1/m_errorRate:1;
	if(newErrorRate > m_errorRate && currOverlapInfo.overlapLength+1 >= 31) return;
	if(newTotalErrors > 1 && currOverlapInfo.overlapLength+1 < 31) return; 

	//update with other chars and increase the error count
	for(int idx2 = 0; idx2 < 4; ++idx2)
	{
		char b = ALPHABET[idx2];
		if(b==w[idx+currOverlapInfo.diagonalOffset]) continue;
		
		// Get the pre-updated SA interval
		BWTIntervalPair probe = currOverlapInfo.pair;
		// Compute new SA intervals using prefix b
		BWTAlgorithms::updateBothL(probe, b, pBWT);
		
		if(probe.isValid())
		{
			// Successful extension using b
			// Store the new BWTInterval into BWTExpanedVector
			BWTOverlapInfo newOverlapInfo = currOverlapInfo;
			// increase the error count
			newOverlapInfo.mismatch = currOverlapInfo.mismatch+1;
			newOverlapInfo.overlapLength = currOverlapInfo.overlapLength+1;
			newOverlapInfo.pair = probe;
			newOverlapInfo.updateLocalError(1);
			BWTExpanedVector.push_back(newOverlapInfo);
			// std::cout << b << "\t" << currInfo.mismatch << "\n";
		}
	}// end of for each base
}

// Expand OverlapInfoList wrt w with insertion
// w: AACCCTTTTGG
//       CCTTTT-G
//       CCTTTTG- duplicate?
// Expand the list for insertions upto m_maxIndels
// False indels will generate false intervals with lots of mismatches after long overlap
// This should be avoided by checking local errors after each gap
void OverlapAlgorithm::expandOverlapInfoListByInsertion(BWTOverlapInfo& currOverlapInfo, std::vector<BWTOverlapInfo>& BWTExpanedVector, const std::string& w, size_t idx, const BWT* pBWT) const
{
	for(size_t idx3=1; idx3<=m_maxIndels; idx3++)
	{
		size_t newTotalErrors = currOverlapInfo.getTotalErrors()+idx3;
		double newErrorRate = (double)newTotalErrors/(currOverlapInfo.overlapLength+idx3);
	
		if( newErrorRate > m_errorRate && currOverlapInfo.overlapLength+idx3 >= 31) return;
		if( newTotalErrors > 1 && currOverlapInfo.overlapLength+idx3 < 31) return; 
		
		// Get the pre-updated SA interval
		BWTIntervalPair probe = currOverlapInfo.pair;
		// Compute new SA intervals using prefix 
		// w[idx]+currOverlapInfo.diagonalOffset -> current diagonal w[] index for this list
		// - idx3 -> adjusted for insertion
		if((int)idx + currOverlapInfo.diagonalOffset-(int)idx3 > 0)
			BWTAlgorithms::updateBothL(probe, w[idx+currOverlapInfo.diagonalOffset-idx3], pBWT);
		else
			return;
		
		if(probe.isValid())
		{
			// std::cout << w[(int)idx + currOverlapInfo.diagonalOffset-idx3] << "\t" << m_maxIndels <<"\n";

			// Successful extension using b
			// Store the new BWTInterval into BWTExpanedVector
			BWTOverlapInfo newOverlapInfo = currOverlapInfo;
			newOverlapInfo.insertion = currOverlapInfo.insertion+idx3;
			newOverlapInfo.diagonalOffset = currOverlapInfo.diagonalOffset-idx3;
			newOverlapInfo.overlapLength = currOverlapInfo.overlapLength+idx3;
			newOverlapInfo.pair = probe;
			// newOverlapInfo.updateLocalError(idx3);
			newOverlapInfo.updateLocalInsertion(idx3);
			BWTExpanedVector.push_back(newOverlapInfo);
		}
	}
}

// Expand OverlapInfoList wrt w with deletions
// w: AACCCTTC--A
//       CCTTCTGA
// False indels will generate false intervals with lots of mismatches after long overlap
// This should be avoided by checking local errors

void OverlapAlgorithm::expandOverlapInfoListByDeletion(BWTOverlapInfo& currOverlapInfo, std::vector<BWTOverlapInfo>& BWTExpanedVector, const std::string& w, size_t idx, const BWT* pBWT) const
{
	// Recursively expand the interval list for deletions
	// w: AACCCTTC--A
	// q:    CCTTCTGA
	// q has to extend from A by two rounds of FM-index extensions in order to reach C
	
	// store the failed-to-match SA interval lists, the matched lists will go to BWTExpanedVector
	std::vector<BWTOverlapInfo> BWTDeletionVector;
	BWTDeletionVector.push_back(currOverlapInfo);
	
	for(size_t idx1=1; idx1<=m_maxIndels; idx1++)
	{
		// Store the failed-to-match lists at each round
		std::vector<BWTOverlapInfo> BWTDeletionOneRoundVector;

		// Step 1: expand each list with all possible extensions
		for(size_t idx2=0; idx2<BWTDeletionVector.size(); idx2++)
		{
			// std::cout << idx1 << ":" << idx2 <<":"<< BWTDeletionVector.at(idx2).pair.interval[0].size() <<"\n";
			std::vector<BWTOverlapInfo> BWTDeletionTmpVector;
			BWTOverlapInfo tmpOverlapInfo = BWTDeletionVector.at(idx2);

			size_t newTotalErrors = tmpOverlapInfo.getTotalErrors()+idx1;
			double newErrorRate = (double)newTotalErrors/(tmpOverlapInfo.overlapLength);

			if(newErrorRate > m_errorRate && tmpOverlapInfo.overlapLength+1 >= 31) continue;
			if(newTotalErrors > 1 && tmpOverlapInfo.overlapLength+1 < 31) continue; 

			//extend with all 4 chars and increase the error count
			for(int idx3 = 0; idx3 < 4; ++idx3)
			{
				char b = ALPHABET[idx3];
				// We must skip the extension of interval matching w
				// This interval is no better than earlier matching rule
				if(b==w[idx+tmpOverlapInfo.diagonalOffset]) continue;
				
				// Get the pre-updated SA interval
				BWTIntervalPair probe = tmpOverlapInfo.pair;
				// Compute new SA intervals using prefix b
				BWTAlgorithms::updateBothL(probe, b, pBWT);
				
				if(probe.isValid())
				{
					// Successful extension using b
					// Store the new BWTInterval into BWTExpanedVector
					BWTOverlapInfo newOverlapInfo=tmpOverlapInfo;
					// increase the error count
					newOverlapInfo.deletion = tmpOverlapInfo.deletion+1;
					newOverlapInfo.overlapLength = tmpOverlapInfo.overlapLength;
					newOverlapInfo.pair = probe;
					// newOverlapInfo.updateLocalError(1);
					BWTDeletionTmpVector.push_back(newOverlapInfo);
					// break;
				}
			}
		
			//Step 2: validate if the base after deletion match w, no mismatch following deletions is considered
			for(size_t idx4=0; idx4 < BWTDeletionTmpVector.size(); idx4++)
			{
				BWTOverlapInfo validateOverlapInfo = BWTDeletionTmpVector.at(idx4);
				BWTIntervalPair probe = validateOverlapInfo.pair;
				
				// Compute new SA intervals using prefix b
				int deletionIdx = (int)idx + validateOverlapInfo.diagonalOffset;
				if(deletionIdx > 0)
					BWTAlgorithms::updateBothL(probe, w[deletionIdx], pBWT);
				else
					continue;
				
				if(probe.isValid())
				{
					// Successful extension using b
					// Store the new BWTInterval into BWTExpanedVector
					BWTOverlapInfo newOverlapInfo=validateOverlapInfo;
					newOverlapInfo.overlapLength = validateOverlapInfo.overlapLength+1;
					newOverlapInfo.pair = probe;
					newOverlapInfo.updateLocalDeletion(idx1);
					BWTExpanedVector.push_back(newOverlapInfo);
				}
				else
					BWTDeletionOneRoundVector.push_back(validateOverlapInfo);
			}
			BWTDeletionTmpVector.clear();
		}
		BWTDeletionVector.clear();
		BWTDeletionVector = BWTDeletionOneRoundVector;
		BWTDeletionOneRoundVector.clear();
	}// end of indel size
}

// Calculate the single right extension to the '$' for each the contained blocks
// so that the interval ranges are consistent
void OverlapAlgorithm::terminateContainedBlocks(const std::string& w, const AlignFlags& af, 
							BWTOverlapInfo& currentOverlap, const BWT* pBWT,
							const BWT* pRevBWT, OverlapBlockList* pContainList, OverlapResult& result) const
{
	// std::cout << BWTOverlapInfoVec.at(idx).pair.interval[0] << "\t" << BWTOverlapInfoVec.at(idx).mismatch+1 << "\n";
	// getchar();
	
	// Require high-accurate ending overlap
	if( (currentOverlap.getErrorRate() > m_errorRate ) ||
		(currentOverlap.getLocalErrors() > 0) )
			return;
			
	//assume the ending char may be SNP
	// for(int idx2 = 0; idx2 < 4; ++idx2)
	// {
		// char b = ALPHABET[idx2];
		char b = w[0];
						
		// Get the pre-updated SA interval
		BWTIntervalPair ranges = currentOverlap.pair;
		// Compute new SA intervals using prefix b
		BWTAlgorithms::updateBothL(ranges, b, pBWT);
		
		// Case 1 is indicated by the existence of a non-$ left or right hand extension
		// In this case we return no alignments for the string
		AlphaCount64 left_ext = BWTAlgorithms::getExtCount(ranges.interval[0], pBWT);
		AlphaCount64 right_ext = BWTAlgorithms::getExtCount(ranges.interval[1], pRevBWT);
		if(left_ext.hasDNAChar() || right_ext.hasDNAChar())
		{
			result.isSubstring = true;
			// std::cout << "is substring\n";
			// getchar();
			return;
		}
		else
		{
			BWTIntervalPair probe = ranges;
			BWTAlgorithms::updateBothL(probe, '$', pBWT);
			if(probe.isValid())
			{
				// terminate the contained block and add it to the contained list
				BWTAlgorithms::updateBothR(probe, '$', pRevBWT);
				assert(probe.isValid());
				pContainList->push_back(OverlapBlock(probe, ranges, w.length(), 0, af));
				// std::cout << "has containment!!\n";
				// getchar();
			}
		}
	// }// end of each ACGT
	
}

// Calculate the terminal extension for any overlap block
void OverlapAlgorithm::terminateOverlapBlocks(const AlignFlags& af, BWTOverlapInfo& currentOverlap, const BWT* pBWT,
							OverlapBlockList* pOverlapList) const
{
	if(currentOverlap.getErrorRate() > m_errorRate) return;
	
	// Require high-accurate ending overlap
	if(currentOverlap.getLocalErrors() > 0) return;
	
	// Calculate which of the prefixes that match w[i, l] are terminal
	// These are the proper prefixes (they are the start of a read)
	BWTIntervalPair probe = currentOverlap.pair;
	BWTAlgorithms::updateBothL(probe, '$', pBWT);
	
	// The SA interval reach $
	if(probe.isValid())
	{
		assert(probe.interval[1].lower > 0);
		// std::cout << i << ":\t" << BWTExpanedVector.at(idx).mismatch << "\t"<<
					// BWTExpanedVector.at(idx).insertion << "\t"<< 
					// BWTExpanedVector.at(idx).deletion << "\t"<< overlapLen << "\n";
		// getchar();

		pOverlapList->push_back(OverlapBlock(probe, currentOverlap.pair, currentOverlap.overlapLength, 
								currentOverlap.getTotalErrors(), currentOverlap.insertion, 
								currentOverlap.deletion, af));
	}

}


// Implement inexact overlap using locality-sensitive backward search via FM-index walk
bool OverlapAlgorithm::findOverlapBlocksInexactFMIndexWalk(const std::string& w, const BWT* pBWT, const BWT* pRevBWT, 
							const AlignFlags& af, const int minOverlap, 
							OverlapBlockList* pOverlapList, OverlapBlockList* pContainList, 
							OverlapResult& result) const
{
	// int m_minOverlap = (int)w.length()<minOverlap?w.length()*0.8:minOverlap;

	SAIOverlapTree OverlapTree(w, minOverlap, m_maxIndels, pBWT, pRevBWT, af);

	// SAIOverlapTree has computed seedSize overlap during construction	
	// Extend base by base for overlaps
	std::vector<OverlapBlock> OBTmpResults;

	// The number of extensions consider max deletions
	// for(size_t i = wlength - OverlapTree.getSeedSize() + m_maxIndels; i >= 1; --i)
	while(OverlapTree.getCurrentLength() < w.length()+m_maxIndels)
	{
		if(OverlapTree.isEmpty()) break;

		// Extend one base for all intervals and return intervals terminated with $
		int flag = OverlapTree.extendOverlapOneBase(OBTmpResults);

		if(flag == -3)
		{
			std::cout << "Too many possible overlapping reads: " << "\t" << OverlapTree.size() << "\n"; 
			// getchar();
			return false;
		}
		
		// Store intervals with sufficient overlap reach $, if any.
		for(size_t j=0; j<OBTmpResults.size(); j++)
			pOverlapList->push_back(OBTmpResults.at(j));

		// Some containment reads may be reached earlier due to indels
		OBTmpResults.clear();
		if(OverlapTree.getCurrentLength() >= w.length()-m_maxIndels)
		{
			bool isSubstring = OverlapTree.terminateContainedBlocks(OBTmpResults);

			if(isSubstring)
			{
				result.isSubstring=true;
				return false;
			}
			else
			{
				// push contained blocks into containment list
				for(size_t j=0; j<OBTmpResults.size(); j++)
					pContainList->push_back(OBTmpResults.at(j));
				
				OBTmpResults.clear();
			}
		}

	}// end of BWT update using w[i]
	
	return true;	
}

// Calculate the irreducible blocks from the vector of OverlapBlocks
void OverlapAlgorithm::computeIrreducibleBlocks(const BWT* pBWT, const BWT* pRevBWT, 
OverlapBlockList* pOBList, 
OverlapBlockList* pOBFinal) const
{
	// processIrreducibleBlocks requires the pOBList to be sorted in descending order
	pOBList->sort(OverlapBlock::sortSizeDescending);
	// if(m_exactModeIrreducible)
	_processIrreducibleBlocksExactIterative(pBWT, pRevBWT, *pOBList, pOBFinal);
	// else
	// _processIrreducibleBlocksInexact(pBWT, pRevBWT, *pOBList, pOBFinal);
	pOBList->clear();
}

// iterate through obList and determine the overlaps that are irreducible. This function is recursive.
// The final overlap blocks corresponding to irreducible overlaps are written to pOBFinal.
// Invariant: the blocks are ordered in descending order of the overlap size so that the longest overlap is first.
// Invariant: each block corresponds to the same extension of the root sequence w.
void OverlapAlgorithm::_processIrreducibleBlocksExactIterative(const BWT* pBWT, const BWT* pRevBWT, 
OverlapBlockList& inList, 
OverlapBlockList* pOBFinal) const
{
	if(inList.empty())
	return;
	
	// We store the overlap blocks in groups of blocks that have the same right-extension.
	// When a branch is found, the groups are split based on the extension
	typedef std::list<OverlapBlockList> BlockGroups;

	BlockGroups blockGroups;
	blockGroups.push_back(inList);
	int numExtensions = 0;
	int numBranches = 0;
	while(!blockGroups.empty())
	{
		// Perform one extenion round for each group.
		// If the top-level block has ended, push the result
		// to the final list and remove the group from processing
		BlockGroups::iterator groupIter = blockGroups.begin();
		BlockGroups incomingGroups; // Branched blocks are placed here

		while(groupIter != blockGroups.end())
		{
			OverlapBlockList& currList = *groupIter;
			bool bEraseGroup = false;

			// Count the extensions in the top level (longest overlap) blocks first
			int topLen = currList.front().overlapLen;
			AlphaCount64 ext_count;
			OBLIter blockIter = currList.begin();
			while(blockIter != currList.end() && blockIter->overlapLen == topLen)
			{
				ext_count += blockIter->getCanonicalExtCount(pBWT, pRevBWT);
				++blockIter;
			}
			
			// Three cases:
			// 1) The top level block has ended as it contains the extension $. Output TLB and end.
			// 2) There is a singular unique extension base for all the blocks. Update the blocks and continue.
			// 3) There are multiple extension bases, split the block group and continue.
			// If some block other than the TLB ended, it must be contained within the TLB and it is not output
			// or considered further. 
			// Likewise if multiple distinct strings in the TLB ended, we only output the top one. The rest
			// must have the same sequence as the top one and are hence considered to be contained with the top element.
			if(ext_count.get('$') > 0)
			{
				// An irreducible overlap has been found. It is possible that there are two top level blocks
				// (one in the forward and reverse direction). Since we can't decide which one
				// contains the other at this point, we output hits to both. Under a fixed 
				// length string assumption one will be contained within the other and removed later.
				OBLIter tlbIter = currList.begin();
				while(tlbIter != currList.end() && tlbIter->overlapLen == topLen)
				{
					AlphaCount64 test_count = tlbIter->getCanonicalExtCount(pBWT, pRevBWT);
					// One of the topLen block does not terminate with $, 
					// implying the other topLen blocks with $ are substring
					if(test_count.get('$') == 0)
					{
						// Remove substring blocks from pOBFinal
						while(tlbIter != currList.begin())
						{
							pOBFinal->erase(--pOBFinal->end());
							tlbIter--;
						}
						// jump to split blocks
						goto SARightExtension;
						// std::cerr << "Error: substring read found during overlap computation.\n";
						// std::cerr << "Please run sga rmdup before sga overlap\n";
						// exit(EXIT_FAILURE);
					}
					
					// Perform the final right-update to make the block terminal
					OverlapBlock branched = *tlbIter;
					BWTAlgorithms::updateBothR(branched.ranges, '$', branched.getExtensionBWT(pBWT, pRevBWT));
					pOBFinal->push_back(branched);
#ifdef DEBUGOVERLAP
					std::cout << "[IE] TLB of length " << branched.overlapLen << " has ended\n";
					std::cout << "[IE]\tBlock data: " << branched << "\n";
#endif             
					++tlbIter;
				} 

				// Set the flag to erase this group, it is finished
				bEraseGroup = true;
			}
			else
			{
SARightExtension:
				// Count the extension for the rest of the blocks
				while(blockIter != currList.end())
				{
					ext_count += blockIter->getCanonicalExtCount(pBWT, pRevBWT);
					++blockIter;
				}

				if(ext_count.hasUniqueDNAChar())
				{
					// Update all the blocks using the unique extension character
					// This character is in the canonical representation wrt to the query
					char b = ext_count.getUniqueDNAChar();
					updateOverlapBlockRangesRight(pBWT, pRevBWT, currList, b);
					numExtensions++;
					bEraseGroup = false;
				}
				else
				{
					// the substring blocks with $ will be discarded
					for(size_t idx = 0; idx < DNA_ALPHABET_SIZE; ++idx)
					{
						char b = ALPHABET[idx];
						if(ext_count.get(b) > 0)
						{
							numBranches++;
							OverlapBlockList branched = currList;
							updateOverlapBlockRangesRight(pBWT, pRevBWT, branched, b);
							incomingGroups.push_back(branched);
							bEraseGroup = true;
						}
					}
				}
			}

			if(bEraseGroup)
			groupIter = blockGroups.erase(groupIter);
			else
			++groupIter;
		}

		// Splice in the newly branched blocks, if any
		blockGroups.splice(blockGroups.end(), incomingGroups);
	}
}

// Classify the blocks in obList as irreducible, transitive or substrings. The irreducible blocks are
// put into pOBFinal. The remaining are discarded.
// Invariant: the blocks are ordered in descending order of the overlap size so that the longest overlap is first.
void OverlapAlgorithm::_processIrreducibleBlocksInexact(const BWT* pBWT, const BWT* pRevBWT, 
OverlapBlockList& activeList, 
OverlapBlockList* pOBFinal) const
{
	if(activeList.empty())
	return;
	
	// The activeList contains all the blocks that are not yet right terminal
	// Count the extensions in the top level (longest) blocks first
	bool all_eliminated = false;
	while(!activeList.empty() && !all_eliminated)
	{
		// The terminalBlock list contains all the blocks that became right-terminal
		// in the current extension round.
		OverlapBlockList terminalList;
		OverlapBlockList potentialContainedList;

		// Perform a single round of extension, any terminal blocks
		// are moved to the terminated list
		extendActiveBlocksRight(pBWT, pRevBWT, activeList, terminalList, potentialContainedList);

		// Compare the blocks in the contained list against the other terminal and active blocks
		// If they are a substring match to any of these, discard them
		OverlapBlockList::iterator containedIter = potentialContainedList.begin();
		for(; containedIter != potentialContainedList.end(); ++containedIter)
		{
			if(!isBlockSubstring(*containedIter, terminalList, m_errorRate) && 
					!isBlockSubstring(*containedIter, activeList, m_errorRate))
			{
				// Not a substring, move to terminal list
				terminalList.push_back(*containedIter);
				//std::cout << "Contained block kept: " << containedIter->overlapLen << "\n";
			}
			else
			{
				//std::cout << "Contained block found and removed: " << containedIter->overlapLen << "\n";
			}
		}

		// Using the terminated blocks, mark as eliminated any active blocks
		// that form a valid overlap to the terminal block. These are transitive edges
		// We do not compare two terminal blocks, we don't consider these overlaps to be
		// transitive
		OverlapBlockList::iterator terminalIter = terminalList.begin();
		for(; terminalIter != terminalList.end(); ++terminalIter)
		{
#ifdef DEBUGOVERLAP
			std::cout << "[II] ***TLB of length " << terminalIter->overlapLen << " has ended\n";
#endif       
			all_eliminated = true;
			OverlapBlockList::iterator activeIter = activeList.begin();
			for(; activeIter != activeList.end(); ++activeIter)
			{
				if(activeIter->isEliminated)
				continue; // skip previously marked blocks
				
				// Two conditions must be met for a block to be transitive wrt terminal:
				// 1) It must have a strictly shorter overlap than the terminal block
				// 2) The error rate between the block and terminal must be less than the threshold
				double inferredErrorRate = calculateBlockErrorRate(*terminalIter, *activeIter);
				if(activeIter->overlapLen < terminalIter->overlapLen && 
						isErrorRateAcceptable(inferredErrorRate, m_errorRate))
				{
#ifdef DEBUGOVERLAP_2                            
					std::cout << "Marking block of length " << activeIter->overlapLen << " as eliminated\n";
#endif
					activeIter->isEliminated = true;
				}
				else
				{
					all_eliminated = false;
				}
			} 
			
			// Move this block to the final list if it has not been previously marked eliminated
			if(!terminalIter->isEliminated)
			{
#ifdef DEBUGOVERLAP
				std::cout << "[II] Adding block " << *terminalIter << " to final list\n";
				//std::cout << "  extension: " << terminalIter->forwardHistory << "\n";
#endif                
				pOBFinal->push_back(*terminalIter);
			}
		}
	}

	activeList.clear();
}

// Extend all the blocks in activeList by one base to the right
// Move all right-terminal blocks to the termainl list. If a block 
// is terminal and potentially contained by another block, add it to 
// containedList
void OverlapAlgorithm::extendActiveBlocksRight(const BWT* pBWT, const BWT* pRevBWT, 
OverlapBlockList& activeList, 
OverlapBlockList& terminalList,
OverlapBlockList& /*containedList*/) const
{
	OverlapBlockList::iterator iter = activeList.begin();
	OverlapBlockList::iterator next;
	while(iter != activeList.end())
	{
		next = iter;
		++next;

		// Check if block is terminal
		AlphaCount64 ext_count = iter->getCanonicalExtCount(pBWT, pRevBWT);
		if(ext_count.get('$') > 0)
		{
			// Only consider this block to be terminal irreducible if it has at least one extension
			// or else it is a substring block
			if(iter->forwardHistory.size() > 0)
			{
				OverlapBlock branched = *iter;
				BWTAlgorithms::updateBothR(branched.ranges, '$', branched.getExtensionBWT(pBWT, pRevBWT));
				terminalList.push_back(branched);
#ifdef DEBUGOVERLAP_2            
				std::cout << "Block of length " << iter->overlapLen << " moved to terminal\n";
#endif
			}
		}

		int curr_extension = iter->forwardHistory.size();

		// Perform the right extensions
		
		// Best case, there is only a single extension character
		// Handle this case specially so we don't need to copy the potentially
		// large OverlapBlock structure and its full history
		if(ext_count.hasUniqueDNAChar())
		{
			// Get the extension character with respect to the queried sequence
			char canonical_base = ext_count.getUniqueDNAChar();

			// Flip the base into the frame of reference for the block
			char block_base = iter->flags.isQueryComp() ? complement(canonical_base) : canonical_base;

			// Update the block using the base in its frame of reference
			BWTAlgorithms::updateBothR(iter->ranges, block_base, iter->getExtensionBWT(pBWT, pRevBWT));

			// Add the base to the history in the frame of reference of the query read
			// This is so the history is consistent when comparing between blocks from different strands
			iter->forwardHistory.add(curr_extension, canonical_base);
		}
		else
		{
			for(size_t idx = 0; idx < DNA_ALPHABET_SIZE; ++idx)
			{
				char canonical_base = ALPHABET[idx];
				char block_base = iter->flags.isQueryComp() ? complement(canonical_base) : canonical_base;
				if(ext_count.get(canonical_base) == 0)
				continue;

				// Branch the sequence. This involves copying the entire history which can be large
				// if the input sequences are very long. This could be avoided by using the SearchHistoyNode/Link
				// structure but branches are infrequent enough to not have a large impact
				OverlapBlock branched = *iter;
				BWTAlgorithms::updateBothR(branched.ranges, block_base, branched.getExtensionBWT(pBWT, pRevBWT));
				assert(branched.ranges.isValid());

				// Add the base in the canonical frame
				branched.forwardHistory.add(curr_extension, canonical_base);

				// Insert the new block after the iterator
				activeList.insert(iter, branched);
			}

			// Remove the original block, which has been superceded by the branches
			activeList.erase(iter);
		}

		iter = next; // this skips the newly-inserted blocks
	}
} 

// Return true if the terminalBlock is a substring of any member of blockList
bool OverlapAlgorithm::isBlockSubstring(OverlapBlock& terminalBlock, const OverlapBlockList& blockList, double maxER) const
{
	OverlapBlockList::const_iterator iter = blockList.begin();
	size_t right_extension_length = terminalBlock.forwardHistory.size();
	for(; iter != blockList.end(); ++iter)
	{
		if(terminalBlock.overlapLen == iter->overlapLen && 
				right_extension_length == iter->forwardHistory.size())
		{
			continue; // same length, cannot be a substring
		}
		
		// Calculate error rate between blocks
		double er = calculateBlockErrorRate(terminalBlock, *iter);
		if(isErrorRateAcceptable(er, maxER))
		return true;
	}
	return false;
}

// Calculate the error rate between two overlap blocks using their history
double OverlapAlgorithm::calculateBlockErrorRate(const OverlapBlock& terminalBlock, const OverlapBlock& otherBlock) const
{
	int back_max = std::min(terminalBlock.overlapLen, otherBlock.overlapLen);
	int backwards_diff = SearchHistoryVector::countDifferences(terminalBlock.backHistory, otherBlock.backHistory, back_max);

	// We compare the forward (right) extension only up to the last position of the terminated block's extension
	int forward_len = terminalBlock.forwardHistory.size();
	int forward_max = forward_len - 1;
	int forward_diff = SearchHistoryVector::countDifferences(terminalBlock.forwardHistory, otherBlock.forwardHistory, forward_max);

	// Calculate the length of the inferred overlap
	int trans_overlap_length = back_max + forward_len;
	double er = static_cast<double>(backwards_diff + forward_diff) / trans_overlap_length;
	
#ifdef DEBUGOVERLAP_2
	std::cout << "OL: " << terminalBlock.overlapLen << "\n";
	std::cout << "TLB BH: " << terminalBlock.backHistory << "\n";
	std::cout << "TB  BH: " << otherBlock.backHistory << "\n";
	std::cout << "TLB FH: " << terminalBlock.forwardHistory << "\n";
	std::cout << "TB  FH: " << otherBlock.forwardHistory << "\n";
	std::cout << "BM: " << back_max << " FM: " << forward_max << "\n";
	std::cout << "IOL: " << trans_overlap_length << " TD: " << (backwards_diff + forward_diff) << "\n";
	std::cout << "Block of length " << otherBlock.overlapLen << " has ier: " << er << "\n";
#endif
	return er;
}

// Update the overlap block list with a righthand extension to b, removing ranges that become invalid
void OverlapAlgorithm::updateOverlapBlockRangesRight(const BWT* pBWT, const BWT* pRevBWT, 
OverlapBlockList& obList, char canonical_base) const
{
	OverlapBlockList::iterator iter = obList.begin(); 
	while(iter != obList.end())
	{
		char relative_base = iter->flags.isQueryComp() ? complement(canonical_base) : canonical_base;
		BWTAlgorithms::updateBothR(iter->ranges, relative_base, iter->getExtensionBWT(pBWT, pRevBWT));
		// remove the block from the list if its no longer valid
		if(!iter->ranges.isValid())
		{
			iter = obList.erase(iter);
		}
		else
		{
			// Add the base to the extension history
			int currExtension = iter->forwardHistory.size();
			iter->forwardHistory.add(currExtension, canonical_base);
			++iter;
		}
	}
}

