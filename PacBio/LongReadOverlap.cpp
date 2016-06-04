///-----------------------------------------------
// Copyright 2016 National Chung Cheng University
// Written by Yao-Ting Huang
// Released under the GPL
//-----------------------------------------------
//
// LongReadOverlap - Overlap computation functions extended from KmerOverlap
// Identify long reads of common seeds and align
//
#include "KmerOverlaps.h"
#include "HashMap.h"
#include "BWTAlgorithms.h"
#include "LongReadOverlap.h"

//
MultipleAlignment LongReadOverlap::buildMultipleAlignment(const std::string& query,
                                                       size_t srcKmerLength,
													   size_t tarKmerLength,
                                                       size_t min_overlap,
                                                       double min_identity,
                                                       size_t coverage,
                                                       BWTIndexSet& indices)
{
    MultipleAlignment multiple_alignment;
    multiple_alignment.addBaseSequence("query", query, "");

	// forward overlap from source seed
    SequenceOverlapPairVector overlap_vector;
	retrieveMatches(query, srcKmerLength, min_overlap, min_identity, coverage, indices, false, overlap_vector);

    for(size_t i = 0; i < overlap_vector.size(); ++i)
        multiple_alignment.addOverlap("Src", overlap_vector[i].sequence[1], "", overlap_vector[i].overlap);

	size_t srcSize = overlap_vector.size();
	// reverse overlap from target seed
	retrieveMatches(query, tarKmerLength, min_overlap, min_identity, coverage, indices, true, overlap_vector);
	
    for(size_t i = srcSize; i < overlap_vector.size(); ++i)
        multiple_alignment.addOverlap("Tar", overlap_vector[i].sequence[1], "", overlap_vector[i].overlap);

    return multiple_alignment;
}

MultipleAlignment LongReadOverlap::endMultipleAlignment(const std::string& query,
                                                       size_t srcKmerLength,
													   size_t min_overlap,
                                                       double min_identity,
                                                       size_t coverage,
                                                       BWTIndexSet& indices)
{
	// forward overlap from source seed
    SequenceOverlapPairVector overlap_vector;
	retrieveMatches(query, srcKmerLength, min_overlap, min_identity, coverage, indices, false, overlap_vector);
	
    MultipleAlignment multiple_alignment;
    multiple_alignment.addBaseSequence("query", query, "");
    for(size_t i = 0; i < overlap_vector.size(); ++i)
        multiple_alignment.addOverlap("null", overlap_vector[i].sequence[1], "", overlap_vector[i].overlap);

    return multiple_alignment;
}

//
void LongReadOverlap::retrieveMatches(const std::string& query,
									size_t k,
									size_t min_overlap,
									double min_identity,
									size_t coverage,
									BWTIndexSet& indices,
									bool isRC,
									SequenceOverlapPairVector& overlap_vector)
{
    assert(indices.pBWT != NULL);
    assert(indices.pRBWT != NULL);
    assert(indices.pSSA != NULL);

	std::vector<std::string> ovlStr;

	// retrive overlap reads
	size_t maxLength = query.length()*1.1+20;
	retrieveStr(query, k, maxLength, indices, isRC, coverage, ovlStr);

	
    // Refine the matches by computing proper overlaps between the sequences
    // Use the overlaps that meet the thresholds to build a multiple alignment
    for(std::vector<std::string>::iterator iter = ovlStr.begin(); iter != ovlStr.end(); ++iter)
    {
        std::string match_sequence = *iter;
			
        // Ignore identical sequence from forward or backward extension
        if( (!isRC && match_sequence.substr(0,query.length()) == query) || 
			(isRC && match_sequence.length() >= query.length() && match_sequence.substr(match_sequence.length()-query.length()) == query))
            continue;

        // Compute the overlap. If the kmer match occurs a single time in each sequence we use
        // the banded extension overlap strategy. Otherwise we use the slow O(M*N) overlapper.
        SequenceOverlap overlap;

		// bandwidth not yet completely tested. < 200 are insufficient dunno why yet. 
        size_t bandwidth = query.length()*0.15+100; // 200;

		// banded global DP alignment, PB requires large mismatch penalty -8
        if(isRC)
			overlap = Overlapper::extendMatch(query, match_sequence, query.length()-k, match_sequence.length()-k, bandwidth, 1, -1, -8);
			// overlap = Overlapper::computeOverlapAffine(query, match_sequence, pacbio_params);
			// overlap = Overlapper::bandedAffineOverlap(query, match_sequence, query.length()-k, match_sequence.length()-k, bandwidth, pacbio_params);
			
		else
			overlap = Overlapper::extendMatch(query, match_sequence, 0, 0, bandwidth, 1, -1, -8);
			// overlap = Overlapper::computeOverlapAffine(query, match_sequence, pacbio_params);
			// overlap = Overlapper::bandedAffineOverlap(query, match_sequence, 0, 0, bandwidth, pacbio_params);
			

        bool bPassedOverlap = (size_t)overlap.getOverlapLength() >= (size_t) min_overlap;
        bool bPassedIdentity = overlap.getPercentIdentity() / 100 >= min_identity;

		// if(query.length() == 10942)
			// std::cout << ">" << overlap.getPercentIdentity() / 100 << ":" << overlap.getOverlapLength() << "\n" ;
					 // << match_sequence.length() << "\n";

        if(bPassedOverlap && bPassedIdentity)
        {
			// std::cout << ">" << overlap.getPercentIdentity() / 100 << ":" << overlap.getOverlapLength() << "\n" ;
			// << match_sequence << "\n";
            SequenceOverlapPair op;
            //op.sequence[0] = query;
            op.sequence[1] = match_sequence;
            op.overlap = overlap;
			op.is_reversed = false;
            overlap_vector.push_back(op);
        }
    }
}

// LF-mapping of each SA index in the interval independently using loop instead of BFS tree expansion
// Contaminated reads often are simple repeats C* or T* with large freq
// Give up this read if the freq is way too large
void LongReadOverlap::retrieveStr(const std::string& query, size_t seedSize, size_t maxLength, 
					BWTIndexSet& indices, bool isRC, size_t coverage, std::vector<std::string>& ovlStr)
{	
	std::string initKmer;
	BWTInterval fwdInterval, rvcInterval;
	size_t totalKmerFreq = 0;
	size_t seedOffSet = 0;
	
	// do
	// {
		if(isRC)
			initKmer = reverseComplement( query.substr(query.length()-seedSize-seedOffSet, seedSize) );
		else
			initKmer = query.substr(0+seedOffSet, seedSize);

		fwdInterval=BWTAlgorithms::findInterval(indices.pRBWT, reverse(initKmer));
		rvcInterval=BWTAlgorithms::findInterval(indices.pBWT, reverseComplement(initKmer));

		size_t kmerFreq = 0;
		kmerFreq += fwdInterval.isValid()?fwdInterval.size():0;
		kmerFreq += fwdInterval.isValid()?rvcInterval.size():0;
		totalKmerFreq += kmerFreq;
			
		// std::cout << seedOffSet << "\t" << initKmer << "\t" << reverseComplement(initKmer) << "\t" << fwdInterval.size() << "\t" << rvcInterval.size() << "\t" << kmerFreq << "\n";
		
		// skip repeat and low-complexity seeds
		if(kmerFreq >= coverage*2 || kmerFreq >= 128) return;
		
		// extend each SA index and collect kmers of smallKmerSize along the extension
		for(int64_t fwdRootIndex = fwdInterval.lower; 
			fwdInterval.isValid() && fwdRootIndex <= fwdInterval.upper 
			&& (fwdRootIndex-fwdInterval.lower < coverage); 
			fwdRootIndex++)
		{		
			// extract the string of fwdIndex via LF mapping
			std::string currStr = initKmer;
			currStr.reserve(maxLength);

			int64_t fwdIndex = fwdRootIndex;
			for(size_t currentLength = initKmer.length(); currentLength < maxLength; currentLength++)
			{
				char b = indices.pRBWT->getChar(fwdIndex);
				if(b == '$') break;
				// currStr = currStr + b;
				currStr.append(1,b);
				fwdIndex = indices.pRBWT->getPC(b) + indices.pRBWT->getOcc(b, fwdIndex - 1);			
			}
			
			if(isRC)
				ovlStr.push_back( reverseComplement(currStr) );
			else
				ovlStr.push_back(currStr);
			
			// if(isRC)
				// std::cout << query << "\n";
				// std::cout << maxLength << ":" << currStr.length() << "\n";
		}
		
		// LF-mapping of each rvc index	
		for(int64_t rvcRootIndex=rvcInterval.lower; 
			rvcRootIndex <= rvcInterval.upper && rvcInterval.isValid() && 
			(rvcRootIndex-rvcInterval.lower < coverage); 
			rvcRootIndex++)
		{
			std::string currStr = reverseComplement(initKmer);
			currStr.reserve(maxLength);

			int64_t rvcIndex = rvcRootIndex;
			
			for(size_t currentLength = initKmer.length(); currentLength < maxLength; currentLength++)
			{
				char b = indices.pBWT->getChar(rvcIndex);
				if(b == '$') break;
				currStr = b + currStr;	// in reverse complement, currStr is before b
				rvcIndex = indices.pBWT->getPC(b) + indices.pBWT->getOcc(b, rvcIndex - 1);
			}
			
			if(isRC)
				ovlStr.push_back(currStr);
			else
				ovlStr.push_back( reverseComplement(currStr) );


		}
	
		// seedOffSet += 10;
	// }
	// while(totalKmerFreq <= 30 && seedOffSet <= 10);

}
