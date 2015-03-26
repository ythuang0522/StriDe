//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// sga - Main assembler driver program
//
#include <string>
#include <iostream>
#include "index.h"
#include "overlap.h"
#include "assemble.h"
#include "oview.h"
#include "preprocess.h"
#include "correct.h"
#include "subgraph.h"
#include "filter.h"
#include "fm-merge.h"
#include "kmerfreq.h"
#include "grep.h"
#include "FMIndexWalk.h"
#include "strideall.h"

#define PROGRAM_BIN "stride"
#define AUTHOR "Yao-Ting Huang"

static const char *STRIDE_VERSION_MESSAGE =
"Stride Assembler  " PACKAGE_VERSION "\n"
"Written by Yao-Ting Huang.\n"
"\n"
"Copyright 2014 National Chung Cheng University\n";

static const char *STRIDE_USAGE_MESSAGE =
"Program: " PACKAGE_NAME "\n"
"Version: " PACKAGE_VERSION "\n"
"Contact: " AUTHOR " [" PACKAGE_BUGREPORT "]\n"
"Usage: " PROGRAM_BIN " <command> [options]\n\n"
"All-in-one Commands:\n"
"      all	  Perform error correction, long-read generation, overlap computation, and assembly in one run\n"
"\nStep-by-step Commands:\n"
"      preprocess  filter and quality-trim reads\n"
"      index       build FM-index for a set of reads\n"
"      correct     correct sequencing errors in reads \n"
"      fmwalk      merge paired reads into long reads via FM-index walk\n"
"      filter      remove redundant reads from a data set\n"
"      overlap     compute overlaps between reads\n"
"      assemble    generate contigs from an assembly graph\n"
"\nOther Commands:\n"
"      merge	merge multiple BWT/FM-index files into a single index\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

int main(int argc, char** argv)
{
    if(argc <= 1)
    {
        std::cout << STRIDE_USAGE_MESSAGE;
        return 0;
    }
    else
    {
        std::string command(argv[1]);
        if(command == "help" || command == "--help")
        {
            std::cout << STRIDE_USAGE_MESSAGE;
            return 0;
        }
        else if(command == "version" || command == "--version")
        {
            std::cout << STRIDE_VERSION_MESSAGE;
            return 0;
        }

		if(command == "all")
            StrideMain(argc - 1, argv + 1);
        else if(command == "preprocess")
            preprocessMain(argc - 1, argv + 1);
        else if(command == "index")
            indexMain(argc - 1, argv + 1);
        else if(command == "filter")
            filterMain(argc - 1, argv + 1);
        else if(command == "fm-merge")
            FMMergeMain(argc - 1, argv + 1);
        else if(command == "overlap")
            overlapMain(argc - 1, argv + 1);
        else if(command == "correct")
            correctMain(argc - 1, argv + 1);
        else if(command == "assemble")
            assembleMain(argc - 1, argv + 1);
        else if(command == "subgraph")
            subgraphMain(argc - 1, argv + 1);
        else if(command == "oview")
            oviewMain(argc - 1, argv + 1);
        else if(command == "kmerfreq")
            kmerfreqMain(argc - 1, argv + 1);
        else if(command == "grep")
            grepMain(argc - 1, argv + 1);
        else if(command == "fmwalk")
            FMindexWalkMain(argc - 1, argv + 1);

        else
        {
            std::cerr << "Unrecognized command: " << command << "\n";
            return 1;
        }
    }

    return 0;
}
