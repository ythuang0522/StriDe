//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// correct - Correct sequencing errors in reads using the FM-index
//
#ifndef FMWALK_H
#define FMWALK_H
#include <getopt.h>
#include "config.h"
#include "BWT.h"
#include "Match.h"
#include "BWTAlgorithms.h"
#include "KmerDistribution.h"

// functions

//
int FMindexWalkMain(int argc, char** argv);

// options
void parseFMWalkOptions(int argc, char** argv);

#endif
