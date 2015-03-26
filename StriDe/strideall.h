//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// subgraph - extract a subgraph from an assembly graph
//
#ifndef STRIDEALL_H
#define STRIDEALL_H
#include <getopt.h>
#include "config.h"

// functions
int StrideMain(int argc, char** argv);
void parseStrideOptions(int argc, char** argv);
void runStride(int argc, char** argv);

#endif
