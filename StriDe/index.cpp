//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// index - Build a BWT/FM-index for a set of reads
//
#include <iostream>
#include <fstream>
#include <algorithm>
#include "SGACommon.h"
#include "Util.h"
#include "index.h"
#include "SuffixArray.h"
#include "SeqReader.h"
#include "SACAInducedCopying.h"
#include "BWT.h"
#include "Timer.h"
#include "BWTCARopebwt.h"
#include "SampledSuffixArray.h"

//
// Getopt
//
#define SUBPROGRAM "index"

static const char *INDEX_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2009 Wellcome Trust Sanger Institute\n";

static const char *INDEX_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... READSFILE\n"
"Index the reads in READSFILE using a suffixarray/bwt\n"
"\n"
"  -v, --verbose                        display verbose output\n"
"      --help                           display this help and exit\n"
"  -a, --algorithm=STR                  BWT construction algorithm. STR can be:\n"
"                                       sais - induced sort algorithm, slower but works for very long sequences\n"
"                                       ropebwt - Li's ropebwt algorithm, suitable for short reads (<200bp) \n"
"                                       ropebwt2 - Li's ropebwt2 algorithm, suitable for short and long reads (default)\n"
"  -t, --threads=NUM                    use NUM threads to construct the index (default: 1)\n"
"  -p, --prefix=PREFIX                  write index to file using PREFIX instead of prefix of READSFILE\n"
"      --no-reverse                     suppress construction of the reverse BWT. Use this option when building the index\n"
"                                       for reads that will be error corrected using the k-mer corrector, which only needs the forward index\n"
"      --no-forward                     suppress construction of the forward BWT. Use this option when building the forward and reverse index separately\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
    static unsigned int verbose;
    static std::string readsFile;
    static std::string prefix;
    static std::string algorithm = "ropebwt2";
    static int numReadsPerBatch = 2000000;
    static int numThreads = 8;
    static bool bDiskAlgo = false;
    static bool bBuildReverse = true;
    static bool bBuildForward = true;
    static bool validate;
    static int gapArrayStorage = 4;
}

static const char* shortopts = "p:a:m:t:d:g:cv";

enum { OPT_HELP = 1, OPT_VERSION, OPT_NO_REVERSE,OPT_NO_FWD };

static const struct option longopts[] = {
    { "verbose",     no_argument,       NULL, 'v' },
    { "check",       no_argument,       NULL, 'c' },
    { "prefix",      required_argument, NULL, 'p' },
    { "threads",     required_argument, NULL, 't' },
    { "disk",        required_argument, NULL, 'd' },
    { "gap-array",   required_argument, NULL, 'g' },
    { "algorithm",   required_argument, NULL, 'a' },
    { "no-reverse",  no_argument,       NULL, OPT_NO_REVERSE },
    { "no-forward",  no_argument,       NULL, OPT_NO_FWD },
    { "help",        no_argument,       NULL, OPT_HELP },
    { "version",     no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

int indexMain(int argc, char** argv)
{
    Timer* pTimer = new Timer("Build FM index");

    parseIndexOptions(argc, argv);
    if(!opt::bDiskAlgo)
    {
        if(opt::algorithm == "sais")
            indexInMemorySAIS();
        else if(opt::algorithm == "ropebwt")
            indexInMemoryRopebwt();
        else if(opt::algorithm == "ropebwt2")
            indexInMemoryRopebwt2();
	}
    else
	std::cout << "Unknown BWT algorithm!\n";
 
	
	delete pTimer;
    return 0;
}


// Wrapper for using RopeBWT by Heng Li
void indexInMemoryRopebwt()
{
    std::cout << "Building index for " << opt::readsFile << " in memory using ropebwt\n";

    bool use_threads = opt::numThreads >= 4;
	std::string bwt_filename = opt::prefix + BWT_EXT;
	std::string rbwt_filename = opt::prefix + RBWT_EXT;
	BWT *pBWT, *pRBWT;

	#pragma omp parallel
	{
		#pragma omp single nowait
		{
			if(opt::bBuildForward)
			{
				BWTCA::runRopebwt(opt::readsFile, bwt_filename, use_threads, false);
				std::cout << "\t done bwt construction, generating .sai file\n";
				pBWT = new BWT(bwt_filename);
			}
		}
		#pragma omp single nowait
		{	
			if(opt::bBuildReverse)
			{
				BWTCA::runRopebwt(opt::readsFile, rbwt_filename, use_threads, true);
				std::cout << "\t done rbwt construction, generating .rsai file\n";
				pRBWT = new BWT(rbwt_filename);
			}
		}
	}
	
	//Construct forward SAI
	if(opt::bBuildForward)
	{
		std::string sai_filename = opt::prefix + SAI_EXT;
		SampledSuffixArray ssa;
		ssa.buildLexicoIndex(pBWT, opt::numThreads);
		ssa.writeLexicoIndex(sai_filename);
		delete pBWT;
	}
	
	//Construct reverse SAI
	if(opt::bBuildReverse)
	{
		std::string rsai_filename = opt::prefix + RSAI_EXT;
		SampledSuffixArray rssa;
		rssa.buildLexicoIndex(pRBWT, opt::numThreads);
		rssa.writeLexicoIndex(rsai_filename);
		delete pRBWT;
	}
}

// Wrapper for using RopeBWT2 by Heng Li, 
//the RLO/RLCO functions are disabled due to the requirement of Paired-End indices of original order
void indexInMemoryRopebwt2()
{
    std::cout << "Building index for " << opt::readsFile << " in memory using RopeBWT2\n";

	std::string bwt_filename = opt::prefix + BWT_EXT;
	std::string rbwt_filename = opt::prefix + RBWT_EXT;
	BWT* pBWT;
	BWT* pRBWT;
	#pragma omp parallel
	{
		#pragma omp single nowait
		{	
			if(opt::bBuildForward)
			{
				BWTCA::runRopebwt2(opt::readsFile, bwt_filename, opt::numThreads, false);
				std::cout << "\t done bwt construction, generating .sai file\n";
				pBWT = new BWT(bwt_filename);
			}
		}
		#pragma omp single nowait
		{	
			if(opt::bBuildReverse)
			{
				BWTCA::runRopebwt2(opt::readsFile, rbwt_filename, opt::numThreads, true);
				std::cout << "\t done rbwt construction, generating .rsai file\n";
				pRBWT = new BWT(rbwt_filename);
			}
		}
	}

	//Construct forward SAI
	if(opt::bBuildForward)
	{	
		std::string sai_filename = opt::prefix + SAI_EXT;
		SampledSuffixArray ssa;
		ssa.buildLexicoIndex(pBWT, opt::numThreads);
		ssa.writeLexicoIndex(sai_filename);
		delete pBWT;
	}

	//Construct reverse SAI
	if(opt::bBuildReverse)
	{
		std::string rsai_filename = opt::prefix + RSAI_EXT;
		SampledSuffixArray rssa;
		rssa.buildLexicoIndex(pRBWT, opt::numThreads);
		rssa.writeLexicoIndex(rsai_filename);
		delete pRBWT;
	}
}

//
void indexInMemorySAIS()
{
    std::cout << "Building index for " << opt::readsFile << " in memory using SAIS\n";

	if(opt::bBuildForward || opt::bBuildReverse)
    {
		// Parse the initial read table
		ReadTable* pRT = new ReadTable(opt::readsFile);

		// Create and write the suffix array for the forward reads
		if(opt::bBuildForward)
		{
			buildIndexForTable(opt::prefix, pRT, false);
		}

		if(opt::bBuildReverse)
		{
			// Reverse all the reads
			pRT->reverseAll();

			// Build the reverse suffix array
			buildIndexForTable(opt::prefix, pRT, true);
		}

		delete pRT;
	}
}


//
void buildIndexForTable(std::string prefix, const ReadTable* pRT, bool isReverse)
{
    // Create suffix array from read table
    SuffixArray* pSA = new SuffixArray(pRT, opt::numThreads);

    if(opt::validate)
    {
        std::cout << "Validating suffix array\n";
        pSA->validate(pRT);
    }

    std::string bwt_filename = prefix + (!isReverse ? BWT_EXT : RBWT_EXT);
    pSA->writeBWT(bwt_filename, pRT);

    std::string sufidx_filename = prefix + (!isReverse ? SAI_EXT : RSAI_EXT);
    pSA->writeIndex(sufidx_filename);

    delete pSA;
    pSA = NULL;
}

//
// Handle command line arguments
//
void parseIndexOptions(int argc, char** argv)
{
	optind=1;
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case 'p': arg >> opt::prefix; break;
            case '?': die = true; break;
            case 'c': opt::validate = true; break;
            case 'd': opt::bDiskAlgo = true; arg >> opt::numReadsPerBatch; break;
            case 't': arg >> opt::numThreads; break;
            case 'g': arg >> opt::gapArrayStorage; break;
            case 'a': arg >> opt::algorithm; break;
            case 'v': opt::verbose++; break;
            case OPT_NO_REVERSE: opt::bBuildReverse = false; break;
            case OPT_NO_FWD: opt::bBuildForward = false; break;
            case OPT_HELP:
                std::cout << INDEX_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << INDEX_VERSION_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }

    // Transform algorithm parameter to lower case
    std::transform(opt::algorithm.begin(), opt::algorithm.end(), opt::algorithm.begin(), ::tolower);

    if (argc - optind < 1)
    {
        std::cerr << SUBPROGRAM ": missing arguments\n";
        die = true;
    }
    else if (argc - optind > 1)
    {
        std::cerr << SUBPROGRAM ": too many arguments\n";
        die = true;
    }

    if(opt::gapArrayStorage != 4 && opt::gapArrayStorage != 8 &&
       opt::gapArrayStorage != 16 && opt::gapArrayStorage != 32)
    {
        std::cerr << SUBPROGRAM ": invalid argument, --gap-array,-g must be one of 4,8,16,32 (found: " << opt::gapArrayStorage << ")\n";
        die = true;
    }

    if(opt::numThreads <= 0)
    {
        std::cerr << SUBPROGRAM ": invalid number of threads: " << opt::numThreads << "\n";
        die = true;
    }

    if(opt::algorithm != "sais" && opt::algorithm != "bcr" && opt::algorithm != "ropebwt" && opt::algorithm != "ropebwt2")
    {
        std::cerr << SUBPROGRAM ": unrecognized algorithm string " << opt::algorithm << ". --algorithm must be sais, bcr or ropebwt\n";
        die = true;
    }

    if (die)
    {
        std::cout << "\n" << INDEX_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

    // Parse the input filenames
    opt::readsFile = argv[optind++];
    //if(opt::prefix.empty())
        opt::prefix = stripFilename(opt::readsFile);

    // Check if input file is empty
    size_t filesize = getFilesize(opt::readsFile);
    if(filesize == 0)
    {
        std::cerr << SUBPROGRAM ": input file is empty\n";
        exit(EXIT_FAILURE);
    }
}
