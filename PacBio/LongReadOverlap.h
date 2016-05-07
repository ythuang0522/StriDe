///-----------------------------------------------
// Copyright 2012 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// KmerOverlaps - Overlap computation functions
// seeded by exact kmer matches
//
#ifndef LONGREADOVERLAP_H
#define LONGREADOVERLAP_H

#include "multiple_alignment.h"
#include "BWTIndexSet.h"
#include "SampledSuffixArray.h"
#include "KmerOverlaps.h"

namespace LongReadOverlap
{
	// Build a multiple alignment for the query, based on initial exact k-matches
	MultipleAlignment buildMultipleAlignment(const std::string& query,
                                         size_t srcKmerLength,
										 size_t tarKmerLength,
                                         size_t min_overlap,
                                         double min_identity,
                                         size_t coverage,
                                         BWTIndexSet& indices);

	// Retrieve matches to the query sequence
	void retrieveMatches(const std::string& query,
                                          size_t k,
                                          size_t min_overlap,
                                          double min_identity,
                                          size_t coverage,
                                          BWTIndexSet& indices,
										  bool isRC,
										  SequenceOverlapPairVector& overlap_vector);

	void retrieveStr(const std::string& query, 
							size_t seedSize,  
							size_t maxLength, 
							BWTIndexSet& indices, 
							bool isRC,
							size_t coverage,
							std::vector<std::string>& ovlStr);
};

#endif
