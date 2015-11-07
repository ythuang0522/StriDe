///-----------------------------------------------
// Copyright 2015 National Chung Cheng University
// Written by Yao-Ting Huang
// Released under the GPL
//-----------------------------------------------
//
// PacBioCorrectionProcess.cpp - Self-correction or hybrid correction using FM-index walk for PacBio reads
//

#ifndef PACBIOCORRECTION_H
#define PACBIOCORRECTION_H
#include <getopt.h>
#include "config.h"
#include "BWT.h"
#include "BWTAlgorithms.h"

// functions

//
int PacBioCorrectionMain(int argc, char** argv);

// options
void parsePacBioCorrectionOptions(int argc, char** argv);

#endif
