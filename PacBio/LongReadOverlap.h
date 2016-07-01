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
#include "PacBioCorrectionProcess.h"


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

	MultipleAlignment endMultipleAlignment(const std::string& query,
										   size_t srcKmerLength,
										   size_t min_overlap,
										   double min_identity,
										   size_t coverage,
										   BWTIndexSet& indices);

	std::string HeadTailSeedMSA(const std::string& query,
								   SeedFeature& head,
								   SeedFeature& tail,
								   size_t min_overlap,
								   double min_identity,
								   size_t maxIndelSize,
								   BWTIndexSet& indices);

	std::string AllSeedMSA(const std::string& query,
									std::vector<SeedFeature> seedVec,
								   size_t min_overlap,
								   double min_identity,
								   size_t maxIndelSize,
								   BWTIndexSet& indices);
				
	std::string HybridMSA(const std::string& query,
							   std::vector<SeedFeature>& seedVec1,
							   std::vector<SeedFeature>& seedVec2,
							   size_t minOverlap,
							   double min_identity,
							   size_t maxIndelSize,
							   BWTIndexSet& indices);
				
	void findOverlapInexact(const std::string& query, size_t srcKmerLength, const BWT* pBWT, const BWT* pRevBWT, 
				size_t minOverlap, double errorRate,  size_t maxIndelSize,
				std::vector<std::string>& overlapStrVector);
				
	void overlapAlign(const std::string& query, std::vector<std::string>& ovlStrVec, size_t minOverlap, double min_identity,
					SequenceOverlapPairVector& overlap_vector, int seedPos, bool isFrontSeed);

	bool overlapStrAlign(const std::string& query, const std::string& extStr, size_t minOverlap, double min_identity,
					SequenceOverlapPairVector& overlap_vector, int queryPos, int extPos);
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
							
	void filterOverlap(SequenceOverlapPairVector& overlap_vector);
};

#endif
