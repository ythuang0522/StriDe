//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// SGUtils - Data structures/Functions related
// to building and manipulating string graphs
//
#ifndef SGUTIL_H
#define SGUTIL_H

#include "Bigraph.h"
#include "ASQG.h"

// typedefs
typedef Bigraph StringGraph;

namespace SGUtil
{
	// Main string graph loading function
	// The allowContainments flag forces the string graph to retain identical vertices
	// Vertices that are substrings of other vertices (SS flag = 1) are never kept
	StringGraph* loadASQG(const std::string& filename, const unsigned int minOverlap, bool allowContainments = false, size_t maxEdges = -1);
	StringGraph* loadASQG(const StringVector & filenameList, const unsigned int minOverlap, bool allowContainments = false, size_t maxEdges = -1,GraphColor c =GC_WHITE);

	//Parallel loading asqg by YTH
	StringGraph* loadASQGVertex(const std::string& filename, const unsigned int minOverlap, bool allowContainments = false, size_t maxEdges=-1);
	StringGraph* loadASQGEdge(std::string ASQGFileName, const unsigned int minOverlap, bool allowContainments, size_t maxEdges, StringGraph* pGraph);

	StringGraph* loadASQG_Parallel(const StringVector & filenameList, const unsigned int minOverlap, bool allowContainments = false, size_t maxEdges = -1,GraphColor c =GC_WHITE);
	StringGraph* loadASQG_EDGE_Parallel(const StringVector & filenameList, const unsigned int minOverlap, bool allowContainments = false, size_t maxEdges = -1,GraphColor c =GC_WHITE , StringGraph* pGraph=NULL);


	// Load a string graph from a fasta file.
	// Returns a graph where each sequence in the fasta is a vertex but there are no edges in the graph.
	StringGraph* loadFASTA(const std::string& filename);

};
#endif
