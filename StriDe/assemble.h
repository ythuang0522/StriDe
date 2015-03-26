//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// assemble - Assemble reads into contigs
//
#ifndef ASSEMBLE_H
#define ASSEMBLE_H
#include <getopt.h>
#include "config.h"
#include "Bigraph.h"
#include "SGVisitors.h"

// functions
int assembleMain(int argc, char** argv);
void parseAssembleOptions(int argc, char** argv);
int assemble();

void outputGraphAndFasta(StringGraph* pGraph , std::string  name , int phase = -1);
void graphTrimAndSmooth (StringGraph* pGraph, size_t trimLength = 400, bool bIsGapPrecent=false);
void RemoveVertexWithBothShortEdges (StringGraph* pGraph, size_t vertexLength, size_t overlapLength, BWT* pBWT = NULL, size_t kmerLength = 0, float threshold = 0 );
void RemoveSmallOverlapRatioEdges ( StringGraph* pGraph, size_t chimeraLength);

#endif
