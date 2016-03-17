//-----------------------------------------------
// Copyright 2016 National Chung Cheng University
// Written by Yao-Ting Huang & Ping-Yeh Chen
// Released under the GPL
//-----------------------------------------------
//
// PacBioHybridCorrection - Hybrid correction using FM-index walk for PacBio reads
//

#ifndef PacBioHybridCorrection_H
#define PacBioHybridCorrection_H
#include <getopt.h>
#include "config.h"
#include "BWT.h"
#include "BWTAlgorithms.h"

// functions

//
int PacBioHybridCorrectionMain(int argc, char** argv);

// options
void parsePacBioHybridCorrectionOptions(int argc, char** argv);

#endif
