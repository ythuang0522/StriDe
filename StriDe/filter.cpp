//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// filter - remove reads from a data set based on various criteria
//
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include "Util.h"
#include "filter.h"
#include "SuffixArray.h"
#include "BWT.h"
#include "SGACommon.h"
#include "OverlapCommon.h"
#include "Timer.h"
#include "BWTAlgorithms.h"
#include "ASQG.h"
#include "gzstream.h"
#include "SequenceProcessFramework.h"
#include "QCProcess.h"
#include "BitVector.h"
#include "BWTCARopebwt.h"


// Defines
#define PROCESS_FILTER_SERIAL SequenceProcessFramework::processSequencesSerial<SequenceWorkItem, QCResult, \
                                                                               QCProcess, QCPostProcess>

#define PROCESS_FILTER_PARALLEL SequenceProcessFramework::processSequencesParallel<SequenceWorkItem, QCResult, \
                                                                                   QCProcess, QCPostProcess>

// Functions

//
// Getopt
//
#define SUBPROGRAM "filter"
static const char *FILTER_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"Revised by Yao-Ting Huang\n"
"Copyright 2014 National Chung Cheng University\n";

static const char *FILTER_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... READSFILE\n"
"Remove reads from a data set.\n"
"The currently available filters are removing exact-match duplicates\n"
"and removing reads with low-frequency k-mers.\n"
"Automatically rebuilds the FM-index without the discarded reads.\n"
"\n"
"      --help                           display this help and exit\n"
"      -v, --verbose                    display verbose output\n"
"      -p, --prefix=PREFIX              use PREFIX for the names of the index files (default: prefix of the input file)\n"
"      -o, --outfile=FILE               write the qc-passed reads to FILE (default: READSFILE.filter.pass.fa)\n"
"      -t, --threads=NUM                use NUM threads to compute the overlaps (default: 1)\n"
"      -d, --sample-rate=N              use occurrence array sample rate of N in the FM-index. Higher values use significantly\n"
"                                       less memory at the cost of higher runtime. This value must be a power of 2 (default: 128)\n"
"      --no-duplicate-check             turn off duplicate removal\n"
"      --substring-only                 when removing duplicates, only remove substring sequences, not full-length matches\n"
"      --no-kmer-check                  turn off the kmer check\n"
"      --homopolymer-check              check reads for hompolymer run length sequencing errors\n"
"      --low-complexity-check           filter out low complexity reads\n"
"\nK-mer filter options:\n"
"      -k, --kmer-size=N                The length of the kmer to use. (default: 31)\n"
"      -x, --kmer-threshold=N           Require at least N kmer coverage for each kmer in a read. (default: 3)\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* PROGRAM_IDENT =
PACKAGE_NAME "::" SUBPROGRAM;

namespace opt
{
    static unsigned int verbose;
    static int numThreads = 8;
    static std::string prefix;
    static std::string readsFile;
    static std::string outFile;
    static std::string discardFile;
    static int sampleRate = BWT::DEFAULT_SAMPLE_RATE_SMALL;

    static bool dupCheck = true;
    static bool substringOnly = false;
    static bool kmerCheck = false;
    static bool hpCheck = false;
    static bool lowComplexityCheck = false;

    static int kmerLength = 31;
    static int kmerThreshold = 3;
}

static const char* shortopts = "p:d:t:o:k:x:v";

enum { OPT_HELP = 1, OPT_VERSION, OPT_SUBSTRING_ONLY, OPT_NO_RMDUP, OPT_NO_KMER, OPT_CHECK_HPRUNS, OPT_CHECK_COMPLEXITY };

static const struct option longopts[] = {
    { "verbose",               no_argument,       NULL, 'v' },
    { "threads",               required_argument, NULL, 't' },
    { "outfile",               required_argument, NULL, 'o' },
    { "prefix",                required_argument, NULL, 'p' },
    { "sample-rate",           required_argument, NULL, 'd' },
    { "kmer-size",             required_argument, NULL, 'k' },
    { "kmer-threshold",        required_argument, NULL, 'x' },
    { "help",                  no_argument,       NULL, OPT_HELP },
    { "version",               no_argument,       NULL, OPT_VERSION },
    { "no-duplicate-check",    no_argument,       NULL, OPT_NO_RMDUP },
    { "no-kmer-check",         no_argument,       NULL, OPT_NO_KMER },
    { "homopolymer-check",     no_argument,       NULL, OPT_CHECK_HPRUNS },
    { "low-complexity-check",  no_argument,       NULL, OPT_CHECK_COMPLEXITY },
    { "substring-only",        no_argument,       NULL, OPT_SUBSTRING_ONLY },
    { NULL, 0, NULL, 0 }
};

//
// Main
//
int filterMain(int argc, char** argv)
{
    parseFilterOptions(argc, argv);
    Timer* pTimer = new Timer(PROGRAM_IDENT);


    BWT* pBWT = new BWT(opt::prefix + BWT_EXT, opt::sampleRate);
    BWT* pRBWT = new BWT(opt::prefix + RBWT_EXT, opt::sampleRate);
    //pBWT->printInfo();

    std::ostream* pWriter = createWriter(opt::outFile);
    std::ostream* pDiscardWriter = createWriter(opt::discardFile);
    QCPostProcess* pPostProcessor = new QCPostProcess(pWriter, pDiscardWriter);

    // If performing duplicate check, create a bitvector to record
    // which reads are duplicates
    BitVector* pSharedBV = NULL;
    if(opt::dupCheck)
        pSharedBV = new BitVector(pBWT->getNumStrings());

    // Set up QC parameters
    QCParameters params;
    params.pBWT = pBWT;
    params.pRevBWT = pRBWT;
    params.pSharedBV = pSharedBV;

    params.checkDuplicates = opt::dupCheck;
    params.substringOnly = opt::substringOnly;
    params.checkKmer = opt::kmerCheck;
    params.checkHPRuns = opt::hpCheck;
    params.checkDegenerate = opt::lowComplexityCheck;

    params.verbose = opt::verbose;

    params.kmerLength = opt::kmerLength;
    params.kmerThreshold = opt::kmerThreshold;

    params.hpKmerLength = 51;
    params.hpHardAcceptCount = 10;
    params.hpMinProportion = 0.1f;
    params.hpMinLength = 6;

    if(opt::numThreads <= 1)
    {
        // Serial mode
        QCProcess processor(params);
        PROCESS_FILTER_SERIAL(opt::readsFile, &processor, pPostProcessor);
    }
    else
    {
        // Parallel mode
        std::vector<QCProcess*> processorVector;
        for(int i = 0; i < opt::numThreads; ++i)
        {
            QCProcess* pProcessor = new QCProcess(params);
            processorVector.push_back(pProcessor);
        }

        PROCESS_FILTER_PARALLEL(opt::readsFile, processorVector, pPostProcessor);

        for(int i = 0; i < opt::numThreads; ++i)
            delete processorVector[i];
    }

    delete pPostProcessor;
    delete pWriter;
    delete pDiscardWriter;

    delete pBWT;
    delete pRBWT;

    if(pSharedBV != NULL)
        delete pSharedBV;

    std::cout << "RE-building index for " << opt::outFile << " in memory using ropebwt2\n";
    std::string prefix=stripFilename(opt::outFile);
        //BWT *pBWT, *pRBWT;
		#pragma omp parallel
		{
			#pragma omp single nowait
			{	
			    std::string bwt_filename = prefix + BWT_EXT;
				BWTCA::runRopebwt2(opt::outFile, bwt_filename, opt::numThreads, false);
				std::cout << "\t done bwt construction, generating .sai file\n";
				pBWT = new BWT(bwt_filename);
			}
			#pragma omp single nowait
			{	
				std::string rbwt_filename = prefix + RBWT_EXT;
				BWTCA::runRopebwt2(opt::outFile, rbwt_filename, opt::numThreads, true);
				std::cout << "\t done rbwt construction, generating .rsai file\n";
				pRBWT = new BWT(rbwt_filename);
			}
		}
        std::string sai_filename = prefix + SAI_EXT;
		SampledSuffixArray ssa;
        ssa.buildLexicoIndex(pBWT, opt::numThreads);
        ssa.writeLexicoIndex(sai_filename);
        delete pBWT;

        std::string rsai_filename = prefix + RSAI_EXT;
        SampledSuffixArray rssa;
        rssa.buildLexicoIndex(pRBWT, opt::numThreads);
        rssa.writeLexicoIndex(rsai_filename);
        delete pRBWT;

    // Cleanup
    delete pTimer;

    return 0;
}

//
// Handle command line arguments
//
void parseFilterOptions(int argc, char** argv)
{
	optind=1;
    std::string algo_str;
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case 'p': arg >> opt::prefix; break;
            case 'o': arg >> opt::outFile; break;
            case 't': arg >> opt::numThreads; break;
            case 'd': arg >> opt::sampleRate; break;
            case 'k': arg >> opt::kmerLength; break;
            case 'x': arg >> opt::kmerThreshold; break;
            case OPT_NO_RMDUP: opt::dupCheck = false; break;
            case OPT_NO_KMER: opt::kmerCheck = false; break;
            case OPT_CHECK_HPRUNS: opt::hpCheck = true; break;
            case OPT_CHECK_COMPLEXITY: opt::lowComplexityCheck = true; break;
            case OPT_SUBSTRING_ONLY: opt::substringOnly = true; break;
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
            case OPT_HELP:
                std::cout << FILTER_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << FILTER_VERSION_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }

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

    if(opt::numThreads <= 0)
    {
        std::cerr << SUBPROGRAM ": invalid number of threads: " << opt::numThreads << "\n";
        die = true;
    }

    if(opt::kmerLength <= 0)
    {
        std::cerr << SUBPROGRAM ": invalid kmer length: " << opt::kmerLength << ", must be greater than zero\n";
        die = true;
    }

    if(opt::kmerThreshold <= 0)
    {
        std::cerr << SUBPROGRAM ": invalid kmer threshold: " << opt::kmerThreshold << ", must be greater than zero\n";
        die = true;
    }

    if (die)
    {
        std::cout << "\n" << FILTER_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

    // Parse the input filenames
    opt::readsFile = argv[optind++];

    if(opt::prefix.empty())
    {
        opt::prefix = stripFilename(opt::readsFile);
    }

    if(opt::outFile.empty())
    {
        opt::outFile = opt::prefix + ".filter.pass.fa";
        opt::discardFile = opt::prefix + ".discard.fa";
    }
    else
    {
        opt::discardFile = stripFilename(opt::outFile) + ".discard.fa";
    }
}
