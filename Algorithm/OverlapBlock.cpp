//-----------------------------------------------
// Copyright 2015 National Chung Cheng University
// Written by Yao-Ting Huang
// Revised from Jared Simpson's overlap block
// Implement new submaximal filtration for indel overlap
// New data structure for indel overlap
// Released under the GPL
//-----------------------------------------------

//-----------------------------------------------
//
// OverlapBlock - Data structures holding
// the result of the alignment of a sequence read
// to a BWT
// 
#include "OverlapBlock.h"
#include "BWTAlgorithms.h"

//#define DEBUG_RESOLVE 1

// 
OverlapBlock::OverlapBlock(BWTIntervalPair r,
                           BWTIntervalPair rawI,
                           int ol, 
                           int nd, 
                           const AlignFlags& af, 
                           const SearchHistoryVector& backHist) : ranges(r), 
                                                                  rawRanges(rawI),
                                                                  overlapLen(ol), 
                                                                  numDiff(nd),
                                                                  flags(af),
                                                                  isEliminated(false),
                                                                  backHistory(backHist)
{
    backHistory.normalize(af.isQueryComp());
}

// Return a pointer to the BWT that should be used to extend the block
// this is the opposite BWT that was used in the backwards search
const BWT* OverlapBlock::getExtensionBWT(const BWT* pBWT, const BWT* pRevBWT) const
{
    if(!flags.isTargetRev())
        return pRevBWT;
    else
        return pBWT;
}

// 
AlphaCount64 OverlapBlock::getCanonicalExtCount(const BWT* pBWT, const BWT* pRevBWT) const
{
    AlphaCount64 out = BWTAlgorithms::getExtCount(ranges.interval[1], getExtensionBWT(pBWT, pRevBWT));
    if(flags.isQueryComp())
        out.complement();
    return out;
}

// Returns 0 if the BWT used for the overlap step was the forward BWT
int OverlapBlock::getCanonicalIntervalIndex() const
{
    if(!flags.isTargetRev())
        return 0;
    else
        return 1;
}

//
BWTInterval OverlapBlock::getCanonicalInterval() const
{
    return ranges.interval[getCanonicalIntervalIndex()];
}

// Get the string corresponding to the overlap block. This is the string found
// during the backwards search
std::string OverlapBlock::getOverlapString(const std::string& original) const
{
    std::string transformed = backHistory.transform(original, flags.isQueryRev());
    // If the query was reversed, we take the first overlapLen (the search
    // was from the front of the sequence) otherwise we take the last overlapLen
    if(flags.isQueryRev())
        return transformed.substr(0, overlapLen);
    else
        return transformed.substr(transformed.length() - overlapLen);
}

// Get the full string corresponding to this overlapblock using the forward history
std::string OverlapBlock::getFullString(const std::string& original) const
{
    std::string str = getOverlapString(original);
    std::string history = forwardHistory.getBaseString();

    if(history.empty() && overlapLen != (int)original.size())
    {
        WARN_ONCE("getFullString() called on block with no history")
    }
/*
    std::cout << "OVERLAP: " << str << "\n";
    std::cout << "HIST: " << history << "\n";
    std::cout << "QREV: " << flags.isQueryRev() << "\n";
    std::cout << "RC: " << flags.isReverseComplement() << "\n";
    std::cout << "QC: " << flags.isQueryComp() << "\n";
*/
    if(!flags.isQueryRev())
    {
        str.append(history);
    }
    else
    {
        history = reverse(history);
        history.append(str);
        str.swap(history);
    }

    if(flags.isReverseComplement())
        str = reverseComplement(str);
    return str;
}


EdgeDir OverlapBlock::getEdgeDir() const
{
    if(flags.isQueryRev())
        return ED_ANTISENSE;
    else
        return ED_SENSE;
}

//
Overlap OverlapBlock::toOverlap(const std::string queryID, const std::string targetID, int queryLen, int targetLen) const
{
	// std::cout << "###" << queryLen << "\t" << overlapLen << "\t" << targetLen <<"\t"
			// << numInsertion << "\t" << numDeletion <<"\n";
    // Compute the sequence coordinates
    int s1 = queryLen - overlapLen;
    int e1 = s1 + overlapLen - 1;
    SeqCoord sc1(s1, e1, queryLen);
	
	// The start of the second hit must be zero by definition of a prefix/suffix match
	int s2 = 0; 
	// deletion insertion
	// AT--CC   GGAATT
	//  TAACCC   G--TTA
    int e2 = s2 + overlapLen - 1 - numInsertion + numDeletion;	
    SeqCoord sc2(s2, e2, targetLen);

    // The coordinates are always with respect to the read, so flip them if
    // we aligned to/from the reverse of the read
    if(flags.isQueryRev())
    {
        sc1.flip();
    }
    if(flags.isTargetRev())
    {
        sc2.flip();
    }
    bool isRC = flags.isReverseComplement();

	assert(sc1.isExtreme() && sc2.isExtreme());
    Overlap o(queryID, sc1, targetID, sc2, isRC, numDiff);
    return o;
}

//
std::string OverlapBlock::toCanonicalID() const
{
    std::stringstream ss;
    int ci = getCanonicalIntervalIndex();
    ss << "IDX-" << ranges.interval[ci].lower;
    return ss.str();
}


//
void printBlockList(const OverlapBlockList* pList)
{
    for(OverlapBlockList::const_iterator i = pList->begin(); i != pList->end(); ++i)
    {
        std::cout << "Block: " << *i << "\n";
    }
}

// remove submaximal overlap due to errors
void removeSubMaximalBlocks(OverlapBlockList* pList, const BWT* /*pBWT*/, const BWT* /*pRevBWT*/)
{
    // This algorithm removes any sub-maximal OverlapBlocks from pList
    // The list is sorted by the left coordinate and iterated through
    // if two adjacent blocks overlap they are split into maximal contiguous regions
    // with resolveOverlap. The resulting list is merged back into pList. This process
    // is repeated until each block in pList is a unique range
    // The bookkeeping in the intersecting case could be more efficient 
    // but the vast vast majority of the cases will not have overlapping 
    // blocks.
    pList->sort(OverlapBlock::sortIntervalLeft);
    OverlapBlockList::iterator iter = pList->begin();
    OverlapBlockList::iterator last = pList->end();
    --last;

    while(iter != pList->end())
    {
        OverlapBlockList::iterator next = iter;
        ++next;

        if(next == pList->end())
            break;

        // Check if iter and next overlap
        if(Interval::isIntersecting(iter->ranges.interval[0].lower, iter->ranges.interval[0].upper, 
                                    next->ranges.interval[0].lower, next->ranges.interval[0].upper))
        {
            // OverlapBlockList resolvedList = resolveOverlap(*iter, *next, pBWT, pRevBWT);
            OverlapBlockList resolvedList = resolveOverlap(*iter, *next);
            
            // Merge the new elements in and start back from the beginning of the list
            pList->erase(iter);
            pList->erase(next);
            pList->merge(resolvedList, OverlapBlock::sortIntervalLeft);
			
			// Restart iteration from the beginning
            iter = pList->begin();

            //std::cout << "After splice: \n";
            //printList(pList);
        }
        else
        {
            ++iter;
        }
    }
}

// identify substring blocks after submaximal removal
bool containSubstringBlocks(OverlapBlockList* pList, int querylength)
{
    OverlapBlockList::iterator iter = pList->begin();
    while(iter != pList->end())
    {
		// substring has overlap length > query length
		// After submaximal removal this substring block is still optimal
        if( (*iter).getOverlapLength()>querylength) return true;
        ++iter;
	}
	return false;
}

// A tracing interval maps an interval of an overlap block representing
// a smaller overlap to an overlap block with a larger overlap. This is
// used to determine which intervals are redundant and can be removed.
struct TracingInterval
{
    int64_t foundPosForward;
    int64_t sourcePosReverse;
    BWTInterval tracing;
    BWTIntervalPair updateRanges;
};

typedef std::list<TracingInterval> TracingIntervalList;


// Resolve redundant overlap due to indels
OverlapBlockList resolveOverlap(OverlapBlock& A, OverlapBlock& B)
{
    OverlapBlockList outList;

    // Check if A and B have the same overlap length, if so they must be 
    // identical blocks (resulting from different seeds) and we can remove one
    // if(/*A.overlapLen == B.overlapLen &&*/ A.ranges.interval[0].lower==0 || B.ranges.interval[0].lower==0)
    // {
        // if(A.ranges.interval[0].lower == B.ranges.interval[0].lower &&
           // A.ranges.interval[0].upper == B.ranges.interval[0].upper)
        // {
			// std::cout << A.ranges.interval[0].lower << "\t"  << A.numDiff << "\t" << B.numDiff << "\n";
            // // outList.push_back(A);
            // // return outList;
        // }
    // }

    // Priority: 1. lower error rate 2. longer overlap length
	// Empirical observation, this is often occurred in inexact overlap
    OverlapBlock* pBetter;
    OverlapBlock* pWorse;
    if(A.numDiff < B.numDiff || 
		(A.numDiff == B.numDiff && A.overlapLen > B.overlapLen))
    {
        pBetter = &A;
        pWorse = &B;
    }
    else
    {
        pBetter = &B;
        pWorse = &A;
    }

	// always retain the better block
    outList.push_back(*pBetter);
	
	//identify the duplicated interval
	BWTInterval dupOverlap;
	Interval::intersect(pBetter->ranges.interval[0].lower, pBetter->ranges.interval[0].upper, 
						pWorse->ranges.interval[0].lower, pWorse->ranges.interval[0].upper,
						dupOverlap.lower, dupOverlap.upper);
	
	assert(dupOverlap.isValid());
	
	// Remove the duplicated interval from the worse block
	// case 1: skip identical block
	// |------|	pBetter
	// |------|	pWorse
	if(pBetter->ranges.interval[0].size() != dupOverlap.size())
	{
		// case 2
		// |------|		pBetter
		//    |------|	pWorse
		if(pBetter->ranges.interval[0].lower < pWorse->ranges.interval[0].lower)
			pWorse->ranges.interval[0].lower += dupOverlap.size();

		// case 3
		//    |------|	pBetter
		// |------|		pWorse
		else
			pWorse->ranges.interval[0].upper -= dupOverlap.size();
			
		// Case 4: The worse block may be entirely contained and invalid and don't push
		//  |------|	pBetter
		//    |--|		pWorse
		if(pWorse->ranges.interval[0].isValid())
			outList.push_back(*pWorse);
	}
	
	// Sort the outlist by left coordinate
    outList.sort(OverlapBlock::sortIntervalLeft);
    return outList;
}

// Partition the overlap block list into two lists, 
// one for the containment overlaps and one for the proper overlaps
void partitionBlockList(int readLen, OverlapBlockList* pCompleteList, 
                        OverlapBlockList* pOverlapList, 
                        OverlapBlockList* pContainList)
{
    OverlapBlockList::iterator iter = pCompleteList->begin();
    while(iter != pCompleteList->end())
    {
        if(iter->overlapLen == readLen)
            pContainList->splice(pContainList->end(), *pCompleteList, iter++);
        else
            pOverlapList->splice(pOverlapList->end(), *pCompleteList, iter++);
    }
}

// Filter out full-length (containment) overlaps from the block list
void removeContainmentBlocks(int readLen, OverlapBlockList* pList)
{
    OverlapBlockList::iterator iter = pList->begin();
    while(iter != pList->end())
    {
        if(iter->overlapLen == readLen)
            iter = pList->erase(iter);
        else
            ++iter;
    }    
}

// 
MultiOverlap blockListToMultiOverlap(const SeqRecord& record, OverlapBlockList& blockList)
{
    std::string read_idx = record.id;
    std::string read_seq = record.seq.toString();
    std::string qual = record.qual;
    MultiOverlap out(read_idx, read_seq, qual);

    for(OverlapBlockList::iterator iter = blockList.begin(); iter != blockList.end(); ++iter)
    {
        std::string overlap_string = iter->getOverlapString(read_seq);

        // Compute the endpoints of the overlap
        int s1 = read_seq.length() - iter->overlapLen;
        int e1 = s1 + iter->overlapLen - 1;
        SeqCoord sc1(s1, e1, read_seq.length());

        int s2 = 0; // The start of the second hit must be zero by definition of a prefix/suffix match
        int e2 = s2 + iter->overlapLen - 1;
        SeqCoord sc2(s2, e2, overlap_string.length());

        // The coordinates are always with respect to the read, so flip them if
        // we aligned to/from the reverse of the read
        if(iter->flags.isQueryRev())
            sc1.flip();
        if(iter->flags.isTargetRev())
            sc2.flip();

        bool isRC = false; // since we transformed the original sequence, they are never RC
        if(sc1.isContained())
            continue; // skip containments

        // Add an overlap for each member of the block
        for(int64_t i = iter->ranges.interval[0].lower; i <= iter->ranges.interval[0].upper; ++i)
        {
            Overlap o(read_idx, sc1, makeIdxString(i), sc2, isRC, -1);
            out.add(overlap_string, o);
        }
    }
    return out;
}

// make an id string from a read index
std::string makeIdxString(int64_t idx)
{
    std::stringstream ss;
    ss << idx;
    return ss.str();
}

