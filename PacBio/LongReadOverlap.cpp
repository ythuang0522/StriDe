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
#include "PBOverlapTree.h"

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
	size_t srcSize = overlap_vector.size();
	
	// reverse overlap from target seed
	retrieveMatches(query, tarKmerLength, min_overlap, min_identity, coverage, indices, true, overlap_vector);
	
	// push into multiple alignment matrix
    for(size_t i = 0; i < srcSize; ++i)
        multiple_alignment.addOverlap("Src", overlap_vector[i].sequence[1], "", overlap_vector[i].overlap);
	
    for(size_t i = srcSize; i < overlap_vector.size(); ++i)
        multiple_alignment.addOverlap("Tar", overlap_vector[i].sequence[1], "", overlap_vector[i].overlap);
	
	// filter low-quality overlap
	// if(!overlap_vector.empty())
	// {
		// // sort overlap by identity descend
		// std::sort(overlap_vector.begin(), overlap_vector.end(), SequenceOverlapPair::sortByOverlapIdentityDesc);
		// filterOverlap(overlap_vector);
	// }

	// for(size_t i = 0; i < overlap_vector.size(); ++i)
        // multiple_alignment.addOverlap("NULL", overlap_vector[i].sequence[1], "", overlap_vector[i].overlap);

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

std::string LongReadOverlap::HeadTailSeedMSA(const std::string& query,
                                                       SeedFeature& head,
													   SeedFeature& tail,
                                                       size_t minOverlap,
                                                       double min_identity,
													   size_t maxIndelSize,
                                                       BWTIndexSet& indices)
{
	double errorRate = 0.6;
	SequenceOverlapPairVector overlap_vector;
	const std::string seq = query.substr(head.seedStartPos, tail.seedEndPos-head.seedStartPos+1);

	// std::cout << seedVec.front().seedStr << ":" << seedVec.back().seedStr << "\n";
	// Match the suffix of query to prefixes in FM-index
	std::vector<std::string> overlapSufFwdStrVector; 
	findOverlapInexact(seq, tail.endBestKmerSize, indices.pBWT, indices.pRBWT, minOverlap, errorRate, maxIndelSize, overlapSufFwdStrVector);
	for(size_t i = 0; i < overlapSufFwdStrVector.size(); i++)
		overlapStrAlign(seq, overlapSufFwdStrVector.at(i), minOverlap, min_identity, overlap_vector, seq.length()-17, overlapSufFwdStrVector.at(i).length()-17);
	
	// std::cout << complement(seq) << "\n";
	std::vector<std::string> overlapSufComStrVector;
	findOverlapInexact(complement(seq), tail.endBestKmerSize, indices.pRBWT, indices.pBWT, minOverlap, errorRate, maxIndelSize, overlapSufComStrVector);
	for(size_t i = 0; i < overlapSufComStrVector.size(); i++)
		overlapStrAlign(seq, complement(overlapSufComStrVector.at(i)), minOverlap, min_identity, overlap_vector, seq.length()-17, overlapSufComStrVector.at(i).length()-17);
		
	// Match the prefix of query to suffixes
	// std::cout << reverseComplement(query) << "\n";
	std::vector<std::string> overlapPreRvcStrVector;
	findOverlapInexact(reverseComplement(seq), head.endBestKmerSize, indices.pBWT, indices.pRBWT, minOverlap, errorRate, maxIndelSize, overlapPreRvcStrVector);
	for(size_t i = 0; i < overlapPreRvcStrVector.size(); i++)
		overlapStrAlign(seq, reverseComplement(overlapPreRvcStrVector.at(i)), minOverlap, min_identity, overlap_vector, 0, 0);

	// std::cout << reverse(seq) << "\n";
	std::vector<std::string> overlapPreRevStrVector;
	findOverlapInexact(reverse(seq), head.endBestKmerSize, indices.pRBWT, indices.pBWT, minOverlap, errorRate, maxIndelSize, overlapPreRevStrVector);
	for(size_t i = 0; i < overlapPreRevStrVector.size(); i++)
		overlapStrAlign(seq, reverse(overlapPreRevStrVector.at(i)), minOverlap, min_identity, overlap_vector, 0, 0);

    MultipleAlignment multiple_alignment;
    multiple_alignment.addBaseSequence("query", seq, "");
    for(size_t i = 0; i < overlap_vector.size(); ++i)
        multiple_alignment.addOverlap("null", overlap_vector[i].sequence[1], "", overlap_vector[i].overlap);

	// skip insufficient number of overlapping reads for correction
	if(multiple_alignment.getNumRows() < 7 )
	{
		std::cout << "MSA Failed: " << multiple_alignment.getNumRows() << "\n";
		return "";
	}
	
	std::string consensus;
	if(overlapSufFwdStrVector.size()+overlapSufComStrVector.size()<3)
		consensus = multiple_alignment.calculateBaseConsensus(100000, 5);
	else if(overlapPreRevStrVector.size()+overlapPreRvcStrVector.size()<3)
	{
		consensus = multiple_alignment.calculateBaseConsensus(100000, -1);
		consensus = consensus.substr(consensus.length()/2);
	}
	else 
		consensus = multiple_alignment.calculateBaseConsensus(100000, -1);
	
	// maquery.print(120);

	// return multiple_alignment;
	return consensus;
}

std::string LongReadOverlap::HybridMSA(const std::string& query,
                                               std::vector<SeedFeature>& seedVec1,
                                               std::vector<SeedFeature>& seedVec2,
											   size_t minOverlap,
											   double min_identity,
											   size_t maxIndelSize,
											   BWTIndexSet& indices)
{
	double errorRate = 0.6;
	SequenceOverlapPairVector overlap_vector;
	std::string seq;
	// beginning part
	for(size_t seedIdx=0; seedIdx < seedVec1.size(); seedIdx++)
	{
		seq = query.substr(seedVec1.at(seedIdx).seedStartPos, seedVec1.back().seedEndPos-seedVec1.front().seedStartPos+1);

		// Match the prefix of query to suffixes
		// std::cout << reverseComplement(query) << "\n";
		SeedFeature head = seedVec1.at(seedIdx);
		std::vector<std::string> overlapPreRvcStrVector;
		findOverlapInexact(reverseComplement(seq), head.endBestKmerSize, indices.pBWT, indices.pRBWT, minOverlap, errorRate, maxIndelSize, overlapPreRvcStrVector);
		for(size_t i = 0; i < overlapPreRvcStrVector.size(); i++)
			overlapStrAlign(seq, reverseComplement(overlapPreRvcStrVector.at(i)), minOverlap, min_identity, overlap_vector, 0, 0);

		// std::cout << reverse(seq) << "\n";
		std::vector<std::string> overlapPreRevStrVector;
		findOverlapInexact(reverse(seq), head.endBestKmerSize, indices.pRBWT, indices.pBWT, minOverlap, errorRate, maxIndelSize, overlapPreRevStrVector);
		for(size_t i = 0; i < overlapPreRevStrVector.size(); i++)
			overlapStrAlign(seq, reverse(overlapPreRevStrVector.at(i)), minOverlap, min_identity, overlap_vector, 0, 0);
		
		if(overlap_vector.size() > 5) break;
	}
	
	// middle part
	// for(size_t seedIdx=0; seedIdx < seedVec2.size(); seedIdx++)
	// {
		// if(seedVec2.at(seedIdx).seedStartPos < 6000) continue;
		// if(seedVec2.at(seedIdx).seedStartPos > 1200) break;
		
		// std::string initKmer = seedVec2.at(seedIdx).seedStr;
		// // std::string initKmer = seedVec.at(seedIdx).seedStr.substr(0,seedVec.at(seedIdx).startBestKmerSize);
		// BWTIntervalPair bip=BWTAlgorithms::findIntervalPair(indices.pBWT, indices.pRBWT, initKmer);
		// std::cout << initKmer << "\t" << bip.interval[0].size() << "\t"<< initKmer.length() << "\n";
		
		// std::string extStr;
		// size_t extPos;
		// size_t overlapCount = 0;

		// // extend each fwd SA index both forward and backward
		// for(int64_t fwdRootIndex = bip.interval[0].lower; 
			// bip.interval[0].isValid() && fwdRootIndex <= bip.interval[0].upper; 
			// fwdRootIndex++)
		// {		
			// std::string prefix = initKmer;
			// size_t maxLength = seedVec2.at(seedIdx).seedEndPos + 100;

			// // backward extend the string of fwdIndex via LF mapping
			// int64_t fwdIndex = fwdRootIndex;
			// for(size_t currentLength = initKmer.length(); currentLength < maxLength; currentLength++)
			// {
				// char b = indices.pBWT->getChar(fwdIndex);
				// if(b == '$') break;
				// prefix = b + prefix;
				// fwdIndex = indices.pBWT->getPC(b) + indices.pBWT->getOcc(b, fwdIndex - 1);			
			// }

			// // forward extend the string of revIndex via LF mapping			
			// maxLength = query.length() - seedVec2.at(seedIdx).seedEndPos + 100;

			// // Assume 100bp PB str is unique for speedup
			// const size_t uniqueLen = prefix.length()>=100?100:prefix.length();
			// BWTInterval revbit = BWTAlgorithms::findInterval(indices.pRBWT, reverse(prefix.substr(prefix.length()-uniqueLen)));
			// // BWTInterval revbit = BWTAlgorithms::findInterval(indices.pRBWT, reverse(prefix));

			// assert(revbit.size() == 1);
			// int64_t revIndex = revbit.lower;
			// // std::cout << revbit.size() << ":" << revIndex << ":" << fwdRootIndex << "\n";
			
			// // assert(revbit.size()==1 && revbit.isValid());
			// std::string suffix;
			// // suffix.reserve(maxLength);
			// for(size_t currentLength = 0; currentLength < maxLength; currentLength++)
			// {
				// char b = indices.pRBWT->getChar(revIndex);
				// if(b == '$') break;
				// suffix.append(1, b);
				// revIndex = indices.pRBWT->getPC(b) + indices.pRBWT->getOcc(b, revIndex - 1);			
			// }

			// extStr = prefix + suffix;
			// extPos = prefix.length()-initKmer.length();
			
			// if(extStr.length() < minOverlap) continue;
			
			// bool isFeasibleOverlap = overlapStrAlign(seq, extStr, minOverlap, min_identity, overlap_vector, seedVec2.at(seedIdx).seedStartPos+1, extPos+1);
		// }
		
		// // extend each rvc SA index both forward and backward
		// bip=BWTAlgorithms::findIntervalPair(indices.pBWT, indices.pRBWT, reverseComplement(initKmer));
		// std::cout << initKmer << "\t" << bip.interval[0].size() << "\n";

		// for(int64_t fwdRootIndex = bip.interval[0].lower; 
			// bip.interval[0].isValid() && fwdRootIndex <= bip.interval[0].upper; 
			// fwdRootIndex++)
		// {		
			// std::string prefix = reverseComplement(initKmer);
			// size_t maxLength = query.length() - seedVec2.at(seedIdx).seedStartPos + 100;

			// // backward extend the string of fwdIndex via LF mapping
			// int64_t fwdIndex = fwdRootIndex;
			// for(size_t currentLength = prefix.length(); currentLength < maxLength; currentLength++)
			// {
				// char b = indices.pBWT->getChar(fwdIndex);
				// if(b == '$') break;
				// prefix = b + prefix;
				// fwdIndex = indices.pBWT->getPC(b) + indices.pBWT->getOcc(b, fwdIndex - 1);			
			// }

			// // forward extend the string of revIndex via LF mapping			
			// maxLength = seedVec2.at(seedIdx).seedStartPos + 100;
			
			// // Assume 100bp PB str is unique for speedup
			// const size_t uniqueLen = prefix.length()>=100?100:prefix.length();
			// BWTInterval revbit = BWTAlgorithms::findInterval(indices.pRBWT, reverse(prefix.substr(prefix.length()-uniqueLen)));
			// // BWTInterval revbit = BWTAlgorithms::findInterval(indices.pRBWT, reverse(prefix));
			
			// int64_t revIndex = revbit.lower;
			// assert(revbit.size()==1 && revbit.isValid());
			// std::string suffix;

			// for(size_t currentLength = 0; currentLength < maxLength; currentLength++)
			// {
				// char b = indices.pRBWT->getChar(revIndex);
				// if(b == '$') break;
				// suffix.append(1, b);
				// revIndex = indices.pRBWT->getPC(b) + indices.pRBWT->getOcc(b, revIndex - 1);			
			// }

			// extStr = reverseComplement(prefix + suffix);
			// extPos = suffix.length();
			// if(extStr.length() < minOverlap) continue;

			// bool isFeasibleOverlap = overlapStrAlign(seq, extStr, minOverlap, min_identity, overlap_vector, seedVec2.at(seedIdx).seedStartPos+1, extPos+1);
		// }// end of rvc SA indices	
	// }
	
	// ending part
	// for(size_t seedIdx=seedVec1.size()-1; seedIdx >= 0; seedIdx--)
	// {
		// SeedFeature tail = seedVec1.at(seedIdx);
		// // std::cout << seedVec.front().seedStr << ":" << seedVec.back().seedStr << "\n";
		// // Match the suffix of query to prefixes in FM-index
		// std::vector<std::string> overlapSufFwdStrVector; 
		// findOverlapInexact(seq, tail.endBestKmerSize, indices.pBWT, indices.pRBWT, minOverlap, errorRate, maxIndelSize, overlapSufFwdStrVector);
		// for(size_t i = 0; i < overlapSufFwdStrVector.size(); i++)
			// overlapStrAlign(seq, overlapSufFwdStrVector.at(i), minOverlap, min_identity, overlap_vector, seq.length()-17, overlapSufFwdStrVector.at(i).length()-17);
		
		// // std::cout << complement(seq) << "\n";
		// std::vector<std::string> overlapSufComStrVector;
		// findOverlapInexact(complement(seq), tail.endBestKmerSize, indices.pRBWT, indices.pBWT, minOverlap, errorRate, maxIndelSize, overlapSufComStrVector);
		// for(size_t i = 0; i < overlapSufComStrVector.size(); i++)
			// overlapStrAlign(seq, complement(overlapSufComStrVector.at(i)), minOverlap, min_identity, overlap_vector, seq.length()-17, overlapSufComStrVector.at(i).length()-17);
	// }
	
    MultipleAlignment multiple_alignment;
    multiple_alignment.addBaseSequence("query", seq, "");
    for(size_t i = 0; i < overlap_vector.size(); ++i)
        multiple_alignment.addOverlap("null", overlap_vector[i].sequence[1], "", overlap_vector[i].overlap);

	// skip insufficient number of overlapping reads for correction
	if(multiple_alignment.getNumRows() < 7 )
	{
		std::cout << "MSA Failed: " << multiple_alignment.getNumRows() << "\n";
		return "";
	}
	
	std::string consensus;
	// if(overlapSufFwdStrVector.size()+overlapSufComStrVector.size()<3)
		consensus = multiple_alignment.calculateBaseConsensus(100000, 5);
	// else if(overlapPreRevStrVector.size()+overlapPreRvcStrVector.size()<3)
	// {
		// consensus = multiple_alignment.calculateBaseConsensus(100000, -1);
		// consensus = consensus.substr(consensus.length()/2);
	// }
	// else 
		// consensus = multiple_alignment.calculateBaseConsensus(100000, -1);
	
	// maquery.print(120);

	// return multiple_alignment;
	return consensus;
}

std::string LongReadOverlap::AllSeedMSA(const std::string& query,
                                                       std::vector<SeedFeature> seedVec,
                                                       size_t minOverlap,
                                                       double min_identity,
													   size_t /*maxIndelSize*/,
                                                       BWTIndexSet& indices)
{
	SequenceOverlapPairVector overlap_vector;

	for(size_t seedIdx=0; seedIdx < seedVec.size(); seedIdx++)
	{		
		std::string initKmer = seedVec.at(seedIdx).seedStr;
		BWTIntervalPair bip=BWTAlgorithms::findIntervalPair(indices.pBWT, indices.pRBWT, initKmer);
		// std::cout << initKmer << "\t" << bip.interval[0].size() << "\t"<< initKmer.length() << "\n";
		
		std::string extStr;
		size_t extPos;
		size_t overlapCount = 0;

		// extend each fwd SA index both forward and backward
		for(int64_t fwdRootIndex = bip.interval[0].lower; 
			bip.interval[0].isValid() && fwdRootIndex <= bip.interval[0].upper; 
			fwdRootIndex++)
		{		
			std::string prefix = initKmer;
			size_t maxLength = seedVec.at(seedIdx).seedEndPos + 100;
			std::string strbuffer;
			strbuffer.reserve(maxLength);

			// backward extend the string of fwdIndex via LF mapping
			int64_t fwdIndex = fwdRootIndex;
			for(size_t currentLength = initKmer.length(); currentLength < maxLength; currentLength++)
			{
				char b = indices.pBWT->getChar(fwdIndex);
				if(b == '$') break;
				// prefix = b + prefix;
				strbuffer.push_back(b);
				fwdIndex = indices.pBWT->getPC(b) + indices.pBWT->getOcc(b, fwdIndex - 1);			
			}
			prefix = reverse(strbuffer)+prefix;

			// forward extend the string of revIndex via LF mapping			
			maxLength = query.length() - seedVec.at(seedIdx).seedEndPos + 100;

			// Assume 200bp PB str is unique for speedup
			const size_t uniqueLen = prefix.length()>=200?200:prefix.length();
			BWTInterval revbit = BWTAlgorithms::findInterval(indices.pRBWT, reverse(prefix.substr(prefix.length()-uniqueLen)));
			if(revbit.size() != 1) continue;

			int64_t revIndex = revbit.lower;
			
			std::string suffix;
			suffix.reserve(maxLength);
			for(size_t currentLength = 0; currentLength < maxLength; currentLength++)
			{
				char b = indices.pRBWT->getChar(revIndex);
				if(b == '$') break;
				// suffix.append(1, b);
				suffix.push_back(b);
				revIndex = indices.pRBWT->getPC(b) + indices.pRBWT->getOcc(b, revIndex - 1);
			}

			extStr = prefix + suffix;
			extPos = prefix.length()-initKmer.length();
			
			if(extStr.length() < minOverlap) continue;
			
			bool isFeasibleOverlap = overlapStrAlign(query, extStr, minOverlap, min_identity, overlap_vector, seedVec.at(seedIdx).seedStartPos+1, extPos+1);
		}
		
		// extend each rvc SA index both forward and backward
		bip=BWTAlgorithms::findIntervalPair(indices.pBWT, indices.pRBWT, reverseComplement(initKmer));
		for(int64_t fwdRootIndex = bip.interval[0].lower; 
			bip.interval[0].isValid() && fwdRootIndex <= bip.interval[0].upper; 
			fwdRootIndex++)
		{		
			std::string prefix = reverseComplement(initKmer);
			size_t maxLength = query.length() - seedVec.at(seedIdx).seedStartPos + 100;
			std::string strbuffer;
			strbuffer.reserve(maxLength);
			
			// backward extend the string of fwdIndex via LF mapping
			int64_t fwdIndex = fwdRootIndex;
			for(size_t currentLength = prefix.length(); currentLength < maxLength; currentLength++)
			{
				char b = indices.pBWT->getChar(fwdIndex);
				if(b == '$') break;
				// prefix = b + prefix;
				strbuffer.push_back(b);
				fwdIndex = indices.pBWT->getPC(b) + indices.pBWT->getOcc(b, fwdIndex - 1);			
			}
			prefix = reverse(strbuffer)+prefix;

			// forward extend the string of revIndex via LF mapping			
			maxLength = seedVec.at(seedIdx).seedStartPos + 100;
			
			// Assume 200bp PB str is unique for speedup
			const size_t uniqueLen = prefix.length()>=200?200:prefix.length();
			BWTInterval revbit = BWTAlgorithms::findInterval(indices.pRBWT, reverse(prefix.substr(prefix.length()-uniqueLen)));
			
			int64_t revIndex = revbit.lower;
			if(revbit.size()!=1) continue;

			std::string suffix;
			suffix.reserve(maxLength);
			for(size_t currentLength = 0; currentLength < maxLength; currentLength++)
			{
				char b = indices.pRBWT->getChar(revIndex);
				if(b == '$') break;
				// suffix.append(1, b);
				suffix.push_back(b);
				revIndex = indices.pRBWT->getPC(b) + indices.pRBWT->getOcc(b, revIndex - 1);			
			}

			extStr = reverseComplement(prefix + suffix);
			extPos = suffix.length();
			if(extStr.length() < minOverlap) continue;

			bool isFeasibleOverlap = overlapStrAlign(query, extStr, minOverlap, min_identity, overlap_vector, seedVec.at(seedIdx).seedStartPos+1, extPos+1);
		}// end of rvc SA indices
		
		// sufficient overlap has been found
		// if(overlap_vector.size()>50) break;
	}// end of each seed


    MultipleAlignment multiple_alignment;
    multiple_alignment.addBaseSequence("query", query, "");
    for(size_t i = 0; i < overlap_vector.size(); ++i)
        multiple_alignment.addOverlap("null", overlap_vector[i].sequence[1], "", overlap_vector[i].overlap);

	// skip insufficient number of overlapping reads for correction
	if(multiple_alignment.getNumRows() < 10 )
	{
		// std::cout << "MSA Failed: " << multiple_alignment.getNumRows() << "\n";
		return "";
	}
	
	// std::cout << "MSA rows: " << multiple_alignment.getNumRows() << "\n";
	
	std::string consensus;
	consensus = multiple_alignment.calculateBaseConsensus(1000, 7);
	
	// multiple_alignment.print(120);

	// head and tail are not good enough
	return consensus.substr(20, consensus.length()-40);
}


bool LongReadOverlap::overlapStrAlign(const std::string& query, const std::string& extStr, size_t minOverlap, double min_identity,
					SequenceOverlapPairVector& overlap_vector, int queryPos, int extPos)
{
	// skip original fwd and rvc of query
	if( (extStr.substr(0,query.length()) == query) || 
		(extStr.length() >= query.length() && extStr.substr(extStr.length()-query.length()) == query))
		return false;

	// bandwidth not yet completely tested. < 700 leads to more errors. dunno why. 
    size_t bandwidth = std::min(700, (int)(query.length()*0.15+100)); 

	// banded semi-global DP alignment, PB requires large mismatch penalty -8
	SequenceOverlap overlap = Overlapper::bandedAlignment(query, extStr, queryPos, extPos, bandwidth, 1, -1, -8);

	bool bPassedOverlap = (size_t)overlap.getOverlapLength() >= (size_t) minOverlap;
	bool bPassedIdentity = overlap.getPercentIdentity() / 100 >= min_identity;

	// std::cout << ">" << overlap.getPercentIdentity() / 100 << ":" << overlap.getOverlapLength() << "\n"; 
	// std::cout << extStr << "\n";

	if(bPassedOverlap && bPassedIdentity)
	{
		// std::cout << ">" << overlap.getPercentIdentity() / 100 << ":" << overlap.getOverlapLength() << "\n" ;
		// std::cout << extStr << "\n";
		SequenceOverlapPair op;
		//op.sequence[0] = query;
		op.sequence[1] = extStr;
		op.overlap = overlap;
		op.is_reversed = false;
		overlap_vector.push_back(op);
		return true;
	}
	
	return false;
}

// retrieve overlapping reads sharing sufficient seeds with query
void LongReadOverlap::findOverlapInexact(const std::string& query, size_t srcKmerLength, 
				const BWT* pBWT, const BWT* pRBWT, 
				size_t minOverlap, double errorRate, size_t maxIndelSize,
				std::vector<std::string>& overlapStrVector)
{
	PBOverlapTree pbOvlTree(query, srcKmerLength, minOverlap, maxIndelSize, pBWT, pRBWT, errorRate);
	pbOvlTree.extendOverlap(overlapStrVector);	
	// pbOvlTree.printAll();
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
        size_t bandwidth = 200;//query.length()*0.15+100;

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

		// std::cout << ">" << overlap.getPercentIdentity() / 100 << ":" << overlap.getOverlapLength() << "\n" 
				 // << match_sequence << "\n";

        if(bPassedOverlap && bPassedIdentity)
        {
			// std::cout << ">" << overlap.getPercentIdentity() / 100 << ":" << overlap.getOverlapLength() << "\n" 
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
			&& (fwdRootIndex-fwdInterval.lower < (int)coverage); 
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
			(rvcRootIndex-rvcInterval.lower < (int)coverage); 
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

//Remove low-quality overlap, high-quality overlap often have high identities group at front
void LongReadOverlap::filterOverlap(SequenceOverlapPairVector& overlap_vector)
{
	SequenceOverlapPairVector new_overlap_vector;
	size_t count = 0;

	new_overlap_vector.push_back(overlap_vector.at(0));
	double bestIdentity = overlap_vector.at(0).overlap.getPercentIdentity();

	for(size_t i=1; i< overlap_vector.size(); i++)
	{
		if(bestIdentity-overlap_vector.at(i).overlap.getPercentIdentity() < 5)
		{
			new_overlap_vector.push_back(overlap_vector.at(i));
			count++;
		}

		// std::cout << overlap_vector.at(i).overlap.getPercentIdentity() << "\n";
	}
		
	overlap_vector.clear();
	overlap_vector = new_overlap_vector;
	// getchar();
}
