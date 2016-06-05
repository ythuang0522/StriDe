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
"Correct PacBio reads using hybrid FM-indices of short and long reads\n"
"\n"
"\nMandatory parameters:\n"
"      -p, --prefix=PREFIX              PREFIX of the names of the high-quality index files \n"
"      -f, --PBprefix=PREFIX            PREFIX of the names of the low-quality index files (default: prefix of the input file)\n"
"      -r, --readlen=NUM                Length of high-quality short reads\n"
"      -c, --coverage=N                 Coverage of high-quality short reds\n"
"      -C, --PBcoverage=N               Coverage of PacBio reads\n"
"\nPacBio correction parameters:\n"
"      -o, --outfile=FILE               Write the corrected reads to FILE (default: READSFILE.ec.fa)\n"
"      -t, --threads=NUM                NUM threads for the computation (default: 1)\n"
"      -K, --max-seed-size=N            Length of max seed size in high-quality sshort reads. (default: 31)\n"
"      -k, --min-seed-size=N            Length of min seed size in high-quality short reads. (default: 21)\n"
"      -x, --kmer-threshold=N           FM-index extension threshold. (default: 3)\n"
"      -L, --max-leaves=N               Number of maximum leaves in the search tree. (default: 256)\n"
"      -M, --max-overlap=N              the max overlap during extension (default: read length*0.9)\n"
"      -m, --min-overlap=N              the min overlap during extension (default: read length*0.8)\n"
"      -v, --verbose                    display verbose output\n"
"      --help                           display this help and exit\n"

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
	static int minOverlap = -1;
	static int maxOverlap = -1;
	static int minSeedLength = 21;	// min seed size
	static int seedKmerThreshold = 30;	// min seed frequency threshold
	static int coverage = -1;	// coverage of high-quality short reads
	static int readLen = -1;	// read length of high-quality short reads
	
	// PacBio Parameteres
	static std::string PBprefix;	// prefix of FM-index of low quality long reads
	static size_t PBKmerLength = 17;	// seed size in PacBio index
	static size_t PBcoverage = 60;		// coverage of PacBio
	static size_t PBSearchDepth = 1000;	// PB seed searh depth

}

static const char* shortopts = "p:t:o:K:x:L:m:k:M:f:r:c:C:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "threads",       required_argument, NULL, 't' },
	{ "outfile",       required_argument, NULL, 'o' },
	{ "prefix",        required_argument, NULL, 'p' },
	{ "max-seed-size", required_argument, NULL, 'K' },
	{ "kmer-threshold",required_argument, NULL, 'x' },
	{ "max-leaves",    required_argument, NULL, 'L' },
	{ "min-overlap",   required_argument, NULL, 'm' },
	{ "min-seed-size", required_argument, NULL, 'k' },
	{ "PBprefix",      required_argument,  NULL, 'f' },
	{ "readLen",       required_argument,  NULL, 'r' },
	{ "coverage",      required_argument, NULL, 'c' },
	{ "PBcoverage",    required_argument, NULL, 'C' },
	{ "verbose",       no_argument,       NULL, 'v' },
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

		if(!opt::PBprefix.empty())
		{
			// Load FM-index of low-quality long reads
			#pragma omp single nowait
			{		//Initialization of large BWT takes some time, pass the disk to next job
				std::cout << std::endl << "Loading BWT: " << opt::PBprefix + BWT_EXT << std::endl;
				plqBWT = new BWT(opt::PBprefix + BWT_EXT, opt::sampleRate);
			}
			#pragma omp single nowait
			{
				std::cout << "Loading RBWT: " << opt::PBprefix + RBWT_EXT << std::endl;
				plqRBWT = new BWT(opt::PBprefix + RBWT_EXT, opt::sampleRate);
			}
			#pragma omp single nowait
			{
				std::cout << "Loading Sampled Suffix Array: " << opt::PBprefix + SAI_EXT << std::endl;
				plqSSA = new SampledSuffixArray(opt::PBprefix + SAI_EXT, SSA_FT_SAI);
			}
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

	// PacBio index
	BWTIndexSet lqindexSet;
	lqindexSet.pBWT = plqBWT;
	lqindexSet.pRBWT = plqRBWT;
	lqindexSet.pSSA = plqSSA;
	ecParams.PBindices = lqindexSet;

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
	ecParams.PBKmerLength = opt::PBKmerLength;
	ecParams.PBcoverage = opt::PBcoverage;
	ecParams.PBSearchDepth = opt::PBSearchDepth;

	
	std::cout << std::endl << "Correcting PacBio reads for " << opt::readsFile << " using--" << std::endl
		<< "number of threads:\t" << opt::numThreads << std::endl
		<< "max seed size:\t" << ecParams.kmerLength << std::endl 
		<< "min seed size:\t" << ecParams.minKmerLength << std::endl
		// << "seed threshold:\t" << ecParams.seedKmerThreshold << std::endl
		<< "max overlap:\t" << ecParams.maxOverlap << std::endl 
		<< "min overlap:\t" << ecParams.minOverlap << std::endl 
		<< "max leaves:\t" << ecParams.maxLeaves << std::endl
		<< "FMW kmer threshold:\t" << ecParams.FMWKmerThreshold << std::endl
		<< "short reads coverage:\t" << ecParams.coverage << std::endl
		<< "PB kmer length:\t" << ecParams.PBKmerLength << std::endl
		<< "PB reads coverage:\t" << ecParams.PBcoverage << std::endl
		<< "PB search depth:\t" << ecParams.PBSearchDepth << std::endl
		<< std::endl;

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
			case 'K': arg >> opt::kmerLength; break;
			case 'x': arg >> opt::kmerThreshold; break;
			case '?': die = true; break;
			case 'v': opt::verbose++; break;
			case 'L': arg >> opt::maxLeaves; break;
			case 'm': arg >> opt::minOverlap; break;
			case 'M': arg >> opt::maxOverlap; break;
			case 'k': arg >> opt::minSeedLength; break;
			case 'c': arg >> opt::coverage; break;
			case 'C': arg >> opt::PBcoverage; break;
			case 'f': arg >> opt::PBprefix; break;
			case 'r': arg >> opt::readLen; break;

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
		std::cerr << "Warnning: coverage of high-quality short reads is invalid or not provided: " << opt::coverage << ", reset to 100 by default\n";
		opt::coverage = 100;
	}

	if(opt::readLen <= 0)
	{
		std::cerr << "Warnning: read length of high-quality short reads is invalid or not provided: " << opt::readLen << ", reset to 100 by default\n";
		opt::readLen = 100;
	}

	if(opt::minOverlap < 0)
		opt::minOverlap = opt::readLen * 0.8 +1;

	if(opt::maxOverlap < 0)
		opt::maxOverlap = opt::readLen * 0.9 +1;
	
	if(opt::prefix.empty())
	{
		std::cerr << SUBPROGRAM ": FM-index of short reads is not provided.\n";
		die = true;
	}
	
	if (die)
	{
		std::cout << "\n" << CORRECT_USAGE_MESSAGE;
		exit(EXIT_FAILURE);
	}

	// Parse the input filenames
	opt::readsFile = argv[optind++];

	if(opt::PBprefix.empty())
	{
		opt::PBprefix = stripFilename(opt::readsFile);
		std::cout << "Warning: FM-index of long reads is not provided. Search for index " << opt::PBprefix << ".bwt/sai/rbwt/rsai\n";
	}
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
