//-----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// bwt2fa - Transform a bwt back into a set of sequences
//
#ifndef KMERFREQ_H
#define KMERFREQ_H
#include <getopt.h>
#include "config.h"

int kmerfreqMain(int argc, char** argv);
void parseKMERFREQOptions(int argc, char** argv);

#endif
