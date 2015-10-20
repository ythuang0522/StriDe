//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// correct - Correct sequencing errors in reads using the FM-index
//
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include "Util.h"
#include "FMIndexWalk.h"
#include "SuffixArray.h"
#include "BWT.h"
#include "SGACommon.h"
#include "OverlapCommon.h"
#include "Timer.h"
#include "BWTAlgorithms.h"
#include "ASQG.h"
#include "gzstream.h"
#include "SequenceProcessFramework.h"
#include "FMIndexWalkProcess.h"
#include "CorrectionThresholds.h"
#include "KmerDistribution.h"
#include "BWTIntervalCache.h"


//
// Getopt
//
#define SUBPROGRAM "FMIndexWalk"
static const char *CORRECT_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Yao-Ting Huang.\n"
"\n"
"Copyright 2015 National Chung Cheng University\n";

static const char *CORRECT_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... READSFILE\n"
"Merge paired end reads into long reads via FM-index walk\n"
"\n"
"      --help                           display this help and exit\n"
"      -v, --verbose                    display verbose output\n"
"      -p, --prefix=PREFIX              use PREFIX for the names of the index files (default: prefix of the input file)\n"
"      -o, --outfile=FILE               write the corrected reads to FILE (default: READSFILE.ec.fa)\n"
"      -t, --threads=NUM                use NUM threads for the computation (default: 1)\n"
"      -a, --algorithm=STR              specify the walking algorithm (hybrid (default), merge, kmerize, or validate)\n"
"                                       pacbioH: pacbio hybrid correction (using NGS reads to correct PB reads)\n"
"                                       pacbioS: pacbio self correction (using PB reads to correct PB reads)\n"
"\nAlgorithm parameters:\n"
"      -k, --kmer-size=N                The length of the kmer to use. (default: 31)\n"
"      -x, --kmer-threshold=N           Attempt to correct kmers that are seen less than N times. (default: 3)\n"
"      -L, --max-leaves=N               Number of maximum leaves in the search tree (default: 32)\n"
"      -I, --max-insertsize=N           the maximum insert size (i.e. search depth) (deault: 400)\n"
"      -m, --min-overlap=N              the min overlap (default: 81)\n"
"      -M, --max-overlap=N              the max overlap (default: avg read length*0.9)\n"
"\nPacBio correction parameters:\n"
"      -k, --kmer-size=N                The length of the kmer to use. (default: 31, recommend: 31 (PacBioH), 17 (PacBioS).)\n"
"      -s, --min-kmer-size=N            The minimum length of the kmer to use. (default: 9, recommend: 21 (PacBioH), 9 (PacBioS).)\n"
"      -x, --kmer-threshold=N           Attempt to correct kmers that are seen less than N times. (default: 3)\n"
"      -y, --seed-kmer-threshold=N      Attempt to find kmers of seed that are seen large than N times. (default: 10)\n"
"      -L, --max-leaves=N               Number of maximum leaves in the search tree. (default: 32)\n"
"      -M, --max-overlap=N              the max overlap (default: -1, recommend: avg read length*0.9 (PacBioH).)\n"
"      -d, --downward=N                 for each possible source, we consider N downward seeds as targets. (default: 3)\n"
"      -c, --collect=N                  for each possible source, we consider N downward seeds to collect reads. (default: 5)\n"

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
	static std::string discardFile = "kmerized.fa";
	static int sampleRate = BWT::DEFAULT_SAMPLE_RATE_SMALL;


	static int kmerLength = 31;
	static int kmerThreshold = 3;
	static bool bLearnKmerParams = false;

	static int maxLeaves=32;
	static int maxInsertSize=400;
	static int minOverlap=81;
	static int maxOverlap=-1;
	
	static int minKmerLength = 9;
	static int seedKmerThreshold = 10;
	static int downward = 3;
	static int collect = 5;

	static FMIndexWalkAlgorithm algorithm = FMW_HYBRID;
}

static const char* shortopts = "p:t:o:a:k:x:L:I:m:M:s:y:d:c:v";

enum { OPT_HELP = 1, OPT_VERSION, OPT_METRICS, OPT_DISCARD, OPT_LEARN };

static const struct option longopts[] = {
	{ "verbose",       no_argument,       NULL, 'v' },
	{ "threads",       required_argument, NULL, 't' },
	{ "outfile",       required_argument, NULL, 'o' },
	{ "prefix",        required_argument, NULL, 'p' },
	{ "algorithm",     required_argument, NULL, 'a' },
	{ "kmer-size",     required_argument, NULL, 'k' },
	{ "kmer-threshold" ,required_argument, NULL, 'x' },
	{ "max-leaves",    required_argument, NULL, 'L' },
	{ "max-insertsize" ,required_argument, NULL, 'I' },
	{ "min-overlap"    ,required_argument, NULL, 'm' },
	{ "min-kmer-size"  ,required_argument, NULL, 's' },
	{ "seed-kmer-threshold"   ,required_argument, NULL, 'y' },
	{ "downward"       ,required_argument, NULL, 'd' },
	{ "collect"        ,required_argument, NULL, 'c' },
	{ "learn",         no_argument,       NULL, OPT_LEARN },
	{ "discard",       no_argument,       NULL, OPT_DISCARD },
	{ "help",          no_argument,       NULL, OPT_HELP },
	{ "version",       no_argument,       NULL, OPT_VERSION },
	{ "metrics",       required_argument, NULL, OPT_METRICS },
	{ NULL, 0, NULL, 0 }
};

//
// Main
//
int FMindexWalkMain(int argc, char** argv)
{
	parseFMWalkOptions(argc, argv);

	// Set the error correction parameters
	FMIndexWalkParameters ecParams;
	BWT *pBWT, *pRBWT;
	SampledSuffixArray* pSSA;

	// Load indices
	#pragma omp parallel
	{
		#pragma omp single nowait
		{	//Initialization of large BWT takes some time, pass the disk to next job
			std::cout << std::endl << "Loading BWT: " << opt::prefix + BWT_EXT << "\n";
			pBWT = new BWT(opt::prefix + BWT_EXT, opt::sampleRate);
		}
		#pragma omp single nowait
		{
			std::cout << "Loading RBWT: " << opt::prefix + RBWT_EXT << "\n";
			pRBWT = new BWT(opt::prefix + RBWT_EXT, opt::sampleRate);
		}
		#pragma omp single nowait
		{
			std::cout << "Loading Sampled Suffix Array: " << opt::prefix + SAI_EXT << "\n";
			pSSA = new SampledSuffixArray(opt::prefix + SAI_EXT, SSA_FT_SAI);
		}
	}

	BWTIndexSet indexSet;
	indexSet.pBWT = pBWT;
	indexSet.pRBWT = pRBWT;
	indexSet.pSSA = pSSA;
	ecParams.indices = indexSet;

	// Sample 100000 kmer counts into KmerDistribution from reverse BWT 
	// Don't sample from forward BWT as Illumina reads are bad at the 3' end
	ecParams.kd = BWTAlgorithms::sampleKmerCounts(opt::minOverlap, 100000, pRBWT);
	ecParams.kd.computeKDAttributes();
	// const size_t RepeatKmerFreq = ecParams.kd.getCutoffForProportion(0.95); 
	std::cout << "Median kmer frequency: " <<ecParams.kd.getMedian() << "\t Std: " <<  ecParams.kd.getSdv() 
	<<"\t 95% kmer frequency: " << ecParams.kd.getCutoffForProportion(0.95)
	<< "\t Repeat frequency cutoff: " << ecParams.kd.getRepeatKmerCutoff() << "\n";
	
	// Open outfiles and start a timer
	std::ostream* pWriter = createWriter(opt::outFile);
	std::ostream* pDiscardWriter = (!opt::discardFile.empty() ? createWriter(opt::discardFile) : NULL);
	Timer* pTimer = new Timer(PROGRAM_IDENT);

	ecParams.algorithm = opt::algorithm;
	ecParams.kmerLength = opt::kmerLength;
	ecParams.printOverlaps = opt::verbose > 0;
	ecParams.maxLeaves = opt::maxLeaves;
	ecParams.maxInsertSize = opt::maxInsertSize;
	ecParams.minOverlap = opt::minOverlap;
	ecParams.maxOverlap = opt::maxOverlap;
	ecParams.minKmerLength = opt::minKmerLength;
	ecParams.seedKmerThreshold = opt::seedKmerThreshold;
	ecParams.FMWKmerThreshold = opt::kmerThreshold;
	ecParams.downward = opt::downward;
	ecParams.collectedSeeds = opt::collect;
	if(ecParams.algorithm == FMW_PACBIOSELF)
	{
		std::cout << std::endl << "Correcting PacBio reads for " << opt::readsFile << " using--" << std::endl
			   	<< "number of threads:\t" << opt::numThreads << std::endl
				<< "max kmer size:\t" << ecParams.kmerLength << std::endl 
				<< "min kmer size:\t" << ecParams.minKmerLength << std::endl
				<< "seed kmer threshold:\t" << ecParams.seedKmerThreshold << std::endl
				<< "FMW max leaves:\t" << ecParams.maxLeaves  << std::endl
				<< "FMW search depth:\t1.4~0.6* (length between two seeds +- 40)" << std::endl
				<< "FMW kmer threshold:\t" << ecParams.FMWKmerThreshold << std::endl
				<< "downward each pair of seeds:\t" << ecParams.downward << std::endl
				<< "number of target to be collected:\t" << ecParams.collectedSeeds << std::endl << std::endl;
	}
	else if(ecParams.algorithm == FMW_PACBIOHYB)
	{
		std::cout << std::endl << "Correcting PacBio reads for " << opt::readsFile << " using--" << std::endl
				<< "number of threads:\t" << opt::numThreads << std::endl
				<< "max kmer size:\t" << ecParams.kmerLength << std::endl 
				<< "min kmer size:\t" << ecParams.minKmerLength << std::endl
				<< "seed kmer threshold:\t" << ecParams.seedKmerThreshold << std::endl
				<< "max distance of searching seed:\t2* tendency distance" << std::endl							
				<< "FMW max overlap:\t" <<  ecParams.maxOverlap << std::endl 
				<< "FMW max leaves:\t" << ecParams.maxLeaves  << std::endl
				<< "FMW search depth:\t1.2~0.8* (length between two seeds +- 10)" << std::endl
				<< "FMW kmer threshold:\t" << ecParams.FMWKmerThreshold << std::endl << std::endl;
		
		// computing distance of various continuous matches length (dk)
		for(int i = 0 ; i <= ecParams.kmerLength ; i++)
		{
			if(i >= ecParams.minKmerLength && i <= ecParams.kmerLength)
				ecParams.seedWalkDistance.push_back(2*3.8649*pow(2.7183,0.1239*i));
			else
				ecParams.seedWalkDistance.push_back(0);
		}
	}
	else
	{
		std::cout << "Merge paired end reads into long reads for " << opt::readsFile << " using \n" 
		<< "min overlap=" <<  ecParams.minOverlap << "\t"
		<< "max overlap=" <<  ecParams.maxOverlap << "\t"
		<< "max leaves=" << opt::maxLeaves << "\t"
		<< "max Insert size=" << opt::maxInsertSize << "\t"
		<< "kmer size=" << opt::kmerLength << "\n\n";
	}
	
	// Setup post-processor
	FMIndexWalkPostProcess postProcessor(pWriter, pDiscardWriter, ecParams);

	if(opt::numThreads <= 1)
	{
		// Serial mode
		FMIndexWalkProcess processor(ecParams);

		if (ecParams.algorithm == FMW_HYBRID || ecParams.algorithm == FMW_MERGE )
		SequenceProcessFramework::processSequencesSerial<SequenceWorkItemPair,
		FMIndexWalkResult,
		FMIndexWalkProcess,
		FMIndexWalkPostProcess>(opt::readsFile, &processor, &postProcessor);

		else
		SequenceProcessFramework::processSequencesSerial<SequenceWorkItem,
		FMIndexWalkResult,
		FMIndexWalkProcess,
		FMIndexWalkPostProcess>(opt::readsFile, &processor, &postProcessor);
	}
	else
	{
		// Parallel mode
		std::vector<FMIndexWalkProcess*> processorVector;
		for(int i = 0; i < opt::numThreads; ++i)
		{
			FMIndexWalkProcess* pProcessor = new FMIndexWalkProcess(ecParams);
			processorVector.push_back(pProcessor);
		}

		if (ecParams.algorithm == FMW_HYBRID || ecParams.algorithm == FMW_MERGE )
		SequenceProcessFramework::processSequencesParallel<SequenceWorkItemPair,
		FMIndexWalkResult,
		FMIndexWalkProcess,
		FMIndexWalkPostProcess>(opt::readsFile, processorVector, &postProcessor);

		else
		SequenceProcessFramework::processSequencesParallel<SequenceWorkItem,
		FMIndexWalkResult,
		FMIndexWalkProcess,
		FMIndexWalkPostProcess>(opt::readsFile, processorVector, &postProcessor);

		for(int i = 0; i < opt::numThreads; ++i)
		{
			delete processorVector[i];
		}
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
void parseFMWalkOptions(int argc, char** argv)
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
		case 'a': arg >> algo_str; break;
		case 'k': arg >> opt::kmerLength; break;
		case 'x': arg >> opt::kmerThreshold; break;
		case '?': die = true; break;
		case 'v': opt::verbose++; break;
		case 'L': arg >> opt::maxLeaves; break;
		case 'I': arg >> opt::maxInsertSize; break;
		case 'm': arg >> opt::minOverlap; break;
		case 'M': arg >> opt::maxOverlap; break;
		case 's': arg >> opt::minKmerLength; break;
		case 'y': arg >> opt::seedKmerThreshold; break;
		case 'd': arg >> opt::downward; break;
		case 'c': arg >> opt::collect; break;
		case OPT_LEARN: opt::bLearnKmerParams = true; break;
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

	// Determine the correction algorithm to use
	if(!algo_str.empty())
	{
		if(algo_str == "merge")
		opt::algorithm = FMW_MERGE;
		else if(algo_str == "hybrid")
		opt::algorithm = FMW_HYBRID;
		else if(algo_str == "kmerize")
		opt::algorithm = FMW_KMERIZE;
		else if(algo_str == "validate")
		opt::algorithm = FMW_VALIDATE;
		else if(algo_str == "pacbioS")
		opt::algorithm = FMW_PACBIOSELF;
		else if(algo_str == "pacbioH")
		opt::algorithm = FMW_PACBIOHYB;
		else
		{
			std::cerr << SUBPROGRAM << ": unrecognized -a,--algorithm parameter: " << algo_str << "\n";
			die = true;
		}
	}

	if (die)
	{
		std::cout << "\n" << CORRECT_USAGE_MESSAGE;
		exit(EXIT_FAILURE);
	}

	// Parse the input filenames
	opt::readsFile = argv[optind++];

	if(opt::prefix.empty())
	{
		opt::prefix = stripFilename(opt::readsFile);
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
	{
		if (opt::algorithm == FMW_HYBRID || opt::algorithm ==  FMW_MERGE )
		opt::outFile = out_prefix + ".merge.fa";
		else if (opt::algorithm == FMW_PACBIOSELF)
			opt::outFile = out_prefix + ".PBSelfCor.fa";
		else if (opt::algorithm == FMW_PACBIOHYB)
			opt::outFile = out_prefix + ".PBHybridCor.fa";
		else
		opt::outFile = out_prefix + ".origin.fa";
	}


	if (opt::algorithm == FMW_KMERIZE || opt::algorithm == FMW_HYBRID || opt::algorithm == FMW_VALIDATE)
	{
		opt::discardFile.clear();
		opt::discardFile = out_prefix + ".kmerized.fa";

	}
}
