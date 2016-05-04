//-----------------------------------------------
// Copyright 2016 National Chung Cheng University
// Written by Yao-Ting Huang & Ping-Yeh Chen
// Released under the GPL
//-----------------------------------------------
//
// PacBioHybridCorrection - Hybrid correction using FM-index walk for PacBio reads
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include "Util.h"
#include "PacBioHybridCorrection.h"
#include "SuffixArray.h"
#include "BWT.h"
#include "SGACommon.h"
#include "OverlapCommon.h"
#include "Timer.h"
#include "BWTAlgorithms.h"
#include "ASQG.h"
#include "gzstream.h"
#include "SequenceProcessFramework.h"
#include "PacBioHybridCorrectionProcess.h"
#include "CorrectionThresholds.h"
#include "BWTIntervalCache.h"


//
// Getopt
//
#define SUBPROGRAM "PacBioHybridCorrection"
static const char *CORRECT_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Yao-Ting Huang & Ping-Yeh Chen.\n"
"\n"
"Copyright 2016 National Chung Cheng University\n";

static const char *CORRECT_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... READSFILE\n"
"Correct PacBio reads via FM-index walk\n"
"\n"
"      --help                           display this help and exit\n"
"      -c, --coverage=N                 coverage of high-quality short reds\n"
"      -v, --verbose                    display verbose output\n"
"      -p, --prefix=PREFIX              use PREFIX for the names of the high-quality index files \n"
"      -f, --lqprefix=PREFIX              use PREFIX for the names of the low-quality index files (default: prefix of the input file)\n"
"      -o, --outfile=FILE               write the corrected reads to FILE (default: READSFILE.ec.fa)\n"
"      -t, --threads=NUM                use NUM threads for the computation (default: 1)\n"
"\nPacBio correction parameters:\n"
"      -k, --kmer-size=N                The length of the kmer to use. (default: 31, recommend: 31 (PacBioH), 17 (PacBioS).)\n"
"      -s, --min-seed-size=N            The minimum length of the kmer to use. (default: 9, recommend: 21 (PacBioH), 9 (PacBioS).)\n"
"      -x, --kmer-threshold=N           Attempt to correct kmers that are seen less than N times. (default: 3)\n"
"      -y, --seed-kmer-threshold=N      Attempt to find kmers of seed that are seen large than N times. (default: 10)\n"
"      -L, --max-leaves=N               Number of maximum leaves in the search tree. (default: 32)\n"
"      -M, --max-overlap=N              the max overlap (default: -1, recommend: avg read length*0.9 (PacBioH).)\n"
"      --split                  		split the uncorrected reads (default: false)\n"

"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* PROGRAM_IDENT =
PACKAGE_NAME "::" SUBPROGRAM;

namespace opt
{
	static unsigned int verbose;
	static int numThreads = 1;
	static std::string prefix;
	static std::string readsFile;
	static std::string outFile;
	static std::string discardFile;
	static int sampleRate = BWT::DEFAULT_SAMPLE_RATE_SMALL;
	static int kmerLength = 31;
	static int kmerThreshold = 3;	// FM-index extension threshold
	static int maxLeaves = 256;
	static int minOverlap = 81;
	static int maxOverlap = -1;
	static int minSeedLength = 21;	// min seed size
	static int seedKmerThreshold = 10;	// min seed frequency threshold
	static int coverage = -1;	// coverage of high-quality short reads
	static std::string lqprefix;	// prefix of FM-index of low quality long reads

}

static const char* shortopts = "p:t:o:k:x:L:m:s:M:y:f:c:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "verbose",       no_argument,       NULL, 'v' },
	{ "threads",       required_argument, NULL, 't' },
	{ "outfile",       required_argument, NULL, 'o' },
	{ "prefix",        required_argument, NULL, 'p' },
	{ "kmer-size",     required_argument, NULL, 'k' },
	{ "kmer-threshold" ,required_argument, NULL, 'x' },
	{ "max-leaves",    required_argument, NULL, 'L' },
	{ "min-overlap"    ,required_argument, NULL, 'm' },
	{ "min-seed-size"  ,required_argument, NULL, 's' },
	{ "seed-kmer-threshold"   ,required_argument, NULL, 'y' },
	{ "coverage"        ,required_argument, NULL, 'c' },
	{ "lqprefix",       required_argument,  NULL, 'f' },
	{ "help",          no_argument,       NULL, OPT_HELP },
	{ "version",       no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

//
// Main
//
int PacBioHybridCorrectionMain(int argc, char** argv)
{
	parsePacBioHybridCorrectionOptions(argc, argv);

	// Set the error correction parameters
	PacBioHybridCorrectionParameters ecParams;
	
	// Load FM-index of high-quality short reads
	BWT *pBWT, *pRBWT;
	SampledSuffixArray* pSSA;

	// Load FM-index of low-quality long reads
	BWT *plqBWT, *plqRBWT;
	SampledSuffixArray* plqSSA;

	// decompression of FM-index from disks is CPU-bound
	#pragma omp parallel
	{
		// Load FM-index of high-quality short reads
		#pragma omp single nowait
		{	//Initialization of large BWT takes some time, pass the disk to next job
			std::cout << std::endl << "Loading BWT: " << opt::prefix + BWT_EXT << std::endl;
			pBWT = new BWT(opt::prefix + BWT_EXT, opt::sampleRate);
		}
		#pragma omp single nowait
		{
			std::cout << "Loading RBWT: " << opt::prefix + RBWT_EXT << std::endl;
			pRBWT = new BWT(opt::prefix + RBWT_EXT, opt::sampleRate);
		}
		#pragma omp single nowait
		{
			std::cout << "Loading Sampled Suffix Array: " << opt::prefix + SAI_EXT << std::endl;
			pSSA = new SampledSuffixArray(opt::prefix + SAI_EXT, SSA_FT_SAI);
		}

		// Load FM-index of low-quality long reads
		#pragma omp single nowait
		{	//Initialization of large BWT takes some time, pass the disk to next job
			std::cout << std::endl << "Loading BWT: " << opt::lqprefix + BWT_EXT << std::endl;
			plqBWT = new BWT(opt::lqprefix + BWT_EXT, opt::sampleRate);
		}
		#pragma omp single nowait
		{
			std::cout << "Loading RBWT: " << opt::lqprefix + RBWT_EXT << std::endl;
			plqRBWT = new BWT(opt::lqprefix + RBWT_EXT, opt::sampleRate);
		}
		#pragma omp single nowait
		{
			std::cout << "Loading Sampled Suffix Array: " << opt::lqprefix + SAI_EXT << std::endl;
			plqSSA = new SampledSuffixArray(opt::lqprefix + SAI_EXT, SSA_FT_SAI);
		}

	}

	// Sample 100000 kmer counts into KmerDistribution from reverse BWT 
	// Don't sample from forward BWT as Illumina reads are bad at the 3' end
	// ecParams.kd = BWTAlgorithms::sampleKmerCounts(opt::kmerLength, 100000, pBWT);
	// ecParams.kd.computeKDAttributes();
	// ecParams.kd.print(100);
	// const size_t RepeatKmerFreq = ecParams.kd.getCutoffForProportion(0.95); 
	// std::cout << "Median kmer frequency: " <<ecParams.kd.getMedian() << "\t Std: " <<  ecParams.kd.getSdv() 
				// <<"\t 95% kmer frequency: " << ecParams.kd.getCutoffForProportion(0.95)
				// << "\t Repeat frequency cutoff: " << ecParams.kd.getRepeatKmerCutoff() << "\n";
	
	BWTIndexSet indexSet;
	indexSet.pBWT = pBWT;
	indexSet.pRBWT = pRBWT;
	indexSet.pSSA = pSSA;
	ecParams.indices = indexSet;

	BWTIndexSet lqindexSet;
	lqindexSet.pBWT = plqBWT;
	lqindexSet.pRBWT = plqRBWT;
	lqindexSet.pSSA = plqSSA;
	ecParams.lq_indices = lqindexSet;

	// Open outfiles and start a timer
	std::ostream* pWriter = createWriter(opt::outFile);
	std::ostream* pDiscardWriter = (!opt::discardFile.empty() ? createWriter(opt::discardFile) : NULL);
	Timer* pTimer = new Timer(PROGRAM_IDENT);

	ecParams.kmerLength = opt::kmerLength;
	ecParams.maxLeaves = opt::maxLeaves;
	ecParams.minOverlap = opt::minOverlap;
	ecParams.maxOverlap = opt::maxOverlap;
	ecParams.minKmerLength = opt::minSeedLength;
	ecParams.seedKmerThreshold = opt::seedKmerThreshold;
	ecParams.FMWKmerThreshold = opt::kmerThreshold;
	ecParams.coverage = opt::coverage;
	
	std::cout << std::endl << "Correcting PacBio reads for " << opt::readsFile << " using--" << std::endl
		<< "number of threads:\t" << opt::numThreads << std::endl
		<< "max kmer size:\t" << ecParams.kmerLength << std::endl 
		<< "min kmer size:\t" << ecParams.minKmerLength << std::endl
		<< "seed kmer threshold:\t" << ecParams.seedKmerThreshold << std::endl
		<< "max distance of searching seed:\t2* tendency distance" << std::endl							
		<< "max overlap:\t" <<  ecParams.maxOverlap << std::endl 
		<< "max leaves:\t" << ecParams.maxLeaves  << std::endl
		<< "search depth:\t1.1~0.9* (length between two seeds +/- 10)" << std::endl
		<< "kmer threshold:\t" << ecParams.FMWKmerThreshold << std::endl << std::endl;

	// Setup post-processor
	PacBioHybridCorrectionPostProcess postProcessor(pWriter, pDiscardWriter, ecParams);

	if(opt::numThreads <= 1)
	{
		// Serial mode
		PacBioHybridCorrectionProcess processor(ecParams);

		SequenceProcessFramework::processSequencesSerial<SequenceWorkItem,
		PacBioHybridCorrectionResult,
		PacBioHybridCorrectionProcess,
		PacBioHybridCorrectionPostProcess>(opt::readsFile, &processor, &postProcessor);
	}
	else
	{
		// Parallel mode
		std::vector<PacBioHybridCorrectionProcess*> processorVector;
		for(int i = 0; i < opt::numThreads; ++i)
		{
			PacBioHybridCorrectionProcess* pProcessor = new PacBioHybridCorrectionProcess(ecParams);
			processorVector.push_back(pProcessor);
		}

		SequenceProcessFramework::processSequencesParallel<SequenceWorkItem,
		PacBioHybridCorrectionResult,
		PacBioHybridCorrectionProcess,
		PacBioHybridCorrectionPostProcess>(opt::readsFile, processorVector, &postProcessor);

		// SequenceProcessFramework::processSequencesParallelOpenMP<SequenceWorkItem,
		// PacBioCorrectionResult,
		// PacBioCorrectionProcess,
		// PacBioCorrectionPostProcess>(opt::readsFile, processorVector, &postProcessor);
		
		for(int i = 0; i < opt::numThreads; ++i)
			delete processorVector[i];
	}

	delete pBWT;
	if(pRBWT != NULL)
		delete pRBWT;
	if(pSSA != NULL)
		delete pSSA;
	delete pTimer;
	delete pWriter;
	if(pDiscardWriter != NULL)
		delete pDiscardWriter;
	
	return 0;
}


//
// Handle command line arguments
//
void parsePacBioHybridCorrectionOptions(int argc, char** argv)
{
	optind=1;	//reset getopt
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
			case 'k': arg >> opt::kmerLength; break;
			case 'x': arg >> opt::kmerThreshold; break;
			case '?': die = true; break;
			case 'v': opt::verbose++; break;
			case 'L': arg >> opt::maxLeaves; break;
			case 'm': arg >> opt::minOverlap; break;
			case 'M': arg >> opt::maxOverlap; break;
			case 's': arg >> opt::minSeedLength; break;
			case 'y': arg >> opt::seedKmerThreshold; break;
			case 'c': arg >> opt::coverage; break;
			case 'f': arg >> opt::lqprefix; break;

			case OPT_HELP:
				std::cout << CORRECT_USAGE_MESSAGE;
				exit(EXIT_SUCCESS);
			case OPT_VERSION:
				std::cout << CORRECT_VERSION_MESSAGE;
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

	if(opt::coverage <= 0)
	{
		std::cerr << SUBPROGRAM ": coverage is invalid or not provided: " << opt::coverage << ", must be greater than zero\n";
		die = true;
	}

	if(opt::prefix.empty())
	{
		std::cerr << SUBPROGRAM ": FM-index of high-quality short reads is not provided.\n";
		die = true;
	}
	
	if (die)
	{
		std::cout << "\n" << CORRECT_USAGE_MESSAGE;
		exit(EXIT_FAILURE);
	}

	// Parse the input filenames
	opt::readsFile = argv[optind++];

	if(opt::lqprefix.empty())
		opt::lqprefix = stripFilename(opt::readsFile);

	// Set the correction threshold
	if(opt::kmerThreshold <= 0)
	{
		std::cerr << "Invalid kmer support threshold: " << opt::kmerThreshold << "\n";
		exit(EXIT_FAILURE);
	}
	CorrectionThresholds::Instance().setBaseMinSupport(opt::kmerThreshold);

	std::string out_prefix = stripFilename(opt::readsFile);
	if(opt::outFile.empty())
		opt::outFile = out_prefix + ".PBHybridCor.fa";
	opt::discardFile = out_prefix + ".discard.fa";
}
