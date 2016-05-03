//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// assemble - Assemble reads into contigs
//
#ifndef ASMLONG_H
#define ASMLONG_H
#include <getopt.h>
#include "config.h"
#include "Bigraph.h"
#include "SGVisitors.h"

// functions
int asmlongMain(int argc, char** argv);
void parseAsmLongOptions(int argc, char** argv);
int asmlong();

void sequentialTrimAndSmooth (StringGraph* pGraph, size_t trimLength = 10000, bool bIsGapPrecent=true);

#endif
