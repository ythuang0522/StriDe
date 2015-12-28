//-----------------------------------------------
// Copyright 2015 National Chung Cheng University
// Written by Yao-Ting Huang
// Released under the GPL
//-----------------------------------------------
//
// asmlong - convert read overlaps into contigs for long reads
//
#include "asmlong.h"
#include "assemble.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <map>
#include "Util.h"
#include "SGUtil.h"
#include "SGAlgorithms.h"
#include "SGVisitors.h"
#include "Timer.h"
#include "EncodedString.h"
#include "SuffixArray.h"
#include "BWT.h"
#include "SGACommon.h"
#include "BWTAlgorithms.h"
#include "BWTIntervalCache.h"
#include "SGSearch.h"

//
// Getopt
//
#define SUBPROGRAM "assemble"
static const char *ASSEMBLE_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Yao-Ting Huang\n"
"\n"
"Copyright 2015 National Chung Cheng University, Taiwan\n";

static const char *ASSEMBLE_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... ASQGFILE\n"
"Create contigs from the assembly graph ASQGFILE.\n"
"\nMandatory arguments:\n"
"      -i,  --insert-size              median size of the long read library\n"
"      -p, --prefix=NAME               prefix of FM-index of paired-end reads (bwt, rbwt, sai, rsai)\n"
"\nSimplify graph parameters:\n"
"      -x, --max-chimera=LEN            maximum chimera length (default: insert size*2)\n"
"      -v, --verbose                        display verbose output\n"
"      --help                           display this help and exit\n"
"\nOther minor control options:\n"
"      -o, --out-prefix=NAME            use NAME as the prefix of the output files (output files will be NAME-contigs.fa, etc)\n"
"      -m, --min-overlap=LEN            only use overlaps of at least LEN. This can be used to filter\n"
"          --transitive-reduction       remove transitive edges from the graph. Off by default.\n"
"          --max-edges=N                limit each vertex to a maximum of N edges. For highly repetitive regions\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";


namespace asmlongopt
{
	static unsigned int verbose;
	//bwt prefix
	static std::string prefix;
	static std::string asqgFile;
	static std::string outContigsFile;
	static std::string outVariantsFile;
	static std::string outGraphFile;
	static std::string outTrimFile;
	static std::string outSimpleBubblesFile;
	static std::string outFilePrefix = "StriDe" ;

	static unsigned int minOverlap=30;
	static size_t maxEdges = 512;

	// Bubble parameters
	static int maxIndelLength = 9;

	//
	static bool bExact = false;

	//FM index files
	BWTIndexSet indices;
	static BWT* pBWT =NULL;
    static BWT* pRBWT =NULL;
    static SampledSuffixArray* pSSA = NULL;

    //Visitor parameters
	static double minOverlapRatio=0.8;
	static size_t insertSize=0;

	static size_t maxChimeraLength = 0;

}

static const char* shortopts = "p:o:m:i:x:v";

enum { OPT_HELP = 1, OPT_VERSION, OPT_EXACT, OPT_MAXINDEL,OPT_MAXEDGES};

static const struct option longopts[] = {
	{ "verbose",               no_argument,       NULL, 'v' },
	{ "out-prefix",            required_argument, NULL, 'o' },
	{ "prefix",                required_argument, NULL, 'p' },
	{ "min-overlap",           required_argument, NULL, 'm' },
	{ "max-indel",             required_argument, NULL, OPT_MAXINDEL },
	{ "max-edges",             required_argument, NULL, OPT_MAXEDGES },
	{ "max-chimera",           required_argument, NULL, 'x' },
	{ "insert-size",           required_argument, NULL, 'i' },
	{ "exact",                 no_argument,       NULL, OPT_EXACT },
	{ "help",                  no_argument,       NULL, OPT_HELP },
	{ "version",               no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

//
// Main
//
int asmlongMain(int argc, char** argv)
{
	Timer* pTimer = new Timer("StriDe assembly");
	parseAsmLongOptions(argc, argv);

	std::cout << "\n#---  Parameters   ---#\n";
	std::cout << "Maximum Chimera Length : " << asmlongopt::maxChimeraLength << std::endl;
	std::cout << "Insert Size            : " << asmlongopt::insertSize << std::endl;

	asmlong();
	delete pTimer;

	return 0;
}

int asmlong()
{
	StringGraph* pGraph;
	#pragma omp parallel
	{
		#pragma omp single nowait
		{
			std::cout << "\n[ Loading string graph: " << asmlongopt::asqgFile <<  " ]\n";
			pGraph=SGUtil::loadASQGVertex(asmlongopt::asqgFile, asmlongopt::minOverlap, true, asmlongopt::maxEdges);
		}
		#pragma omp single nowait
		{
			std::cout << "[ Loading BWT ]\n";
			asmlongopt::pBWT = new BWT(asmlongopt::prefix + BWT_EXT, BWT::DEFAULT_SAMPLE_RATE_SMALL);
		}
		#pragma omp single nowait
		{
			std::cout << "[ Loading RBWT ]\n";
			asmlongopt::pRBWT = new BWT(asmlongopt::prefix + RBWT_EXT, BWT::DEFAULT_SAMPLE_RATE_SMALL);
		}
		#pragma omp single nowait
		{
			std::cout << "[ Loading SAI ]\n";
			asmlongopt::pSSA = new SampledSuffixArray(asmlongopt::prefix + SAI_EXT, SSA_FT_SAI);
		}
	}
    asmlongopt::indices.pBWT = asmlongopt::pBWT;
    asmlongopt::indices.pRBWT = asmlongopt::pRBWT;
    asmlongopt::indices.pSSA = asmlongopt::pSSA;
	
	pGraph=SGUtil::loadASQGEdge(asmlongopt::asqgFile, asmlongopt::minOverlap, true, asmlongopt::maxEdges, pGraph);

	// if(asmlongopt::bExact)
	pGraph->setExactMode(asmlongopt::bExact);
	//pGraph->printMemSize();

	// // Pre-assembly graph stats
	SGGraphStatsVisitor statsVisit;
	std::cout << "[Stats] Input graph:\n";
	pGraph->visitP(statsVisit);
	
	// int phase = 0 ;

	// Remove containments from the graph
	std::cout << "Removing contained vertices from graph\n";
	SGContainRemoveVisitor containVisit;
	while(pGraph->hasContainment())
		pGraph->visit(containVisit);

	/*---Remove Transitive Edges---*/
	std::cout << "Removing transitive edges\n";
	SGTransitiveReductionVisitor trVisit;
	pGraph->visit(trVisit);

	// Compact together unbranched chains of vertices
	std::cout << "Start to simplify unipaths ...\n";
    pGraph->simplify();
	std::cout << "[Stats] Simplified graph:\n";
	pGraph->visitP(statsVisit);


	/*** Pop Bubbles ***/
	std::cout << "\n[ Remove bubbles and tips ]\n";
	sequentialTrimAndSmooth(pGraph, asmlongopt::maxChimeraLength);
	pGraph->contigStats();
	pGraph->visit(statsVisit);
	// outputGraphAndFasta(pGraph,"",++phase);
	// pGraph->writeDot("debug.dot",0);
	// getchar();

    std::cout << "\n[ Remove chimeric edges with insufficient overlap or large-overlap difference from large vertices]\n";
	size_t minOverlapLen = asmlongopt::insertSize*asmlongopt::minOverlapRatio;
	SGRemoveByOverlapLenDiffVisitor srv1(1600, minOverlapLen, asmlongopt::insertSize/10, false);	
	pGraph->visitP(srv1);
	sequentialTrimAndSmooth(pGraph, asmlongopt::maxChimeraLength);

	// Rename requires extra memory which should
	// only be done after the string graph has been greatly simplified.
	pGraph->renameVertices("");

	std::cout << "\n[Stats] Final graph statistic:\n";
	pGraph->contigStats();
	pGraph->visitP(statsVisit);
	
	// Write the results
	std::cout << "\n<Printing the contig file> : " << asmlongopt::outContigsFile << " \n" << std::endl;
	SGFastaVisitor av(asmlongopt::outContigsFile);
	pGraph->visit(av);
	//std::cout << "<Printing the ASQG and dot files>\n" << std::endl;
	pGraph->writeASQG(asmlongopt::outGraphFile);
    pGraph->writeDot("StriDe-graph.dot",0);

	delete pGraph;
	return 0;

}

//
// Handle command line arguments
//
void parseAsmLongOptions(int argc, char** argv)
{
	optind=1;	//reset getopt

	// Set defaults
	bool die = false;
	for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
	{
		std::istringstream arg(optarg != NULL ? optarg : "");
		switch (c)
		{
		case 'o': arg >> asmlongopt::outFilePrefix; break;
		case 'p': arg >> asmlongopt::prefix; break;
		case 'm': arg >> asmlongopt::minOverlap; break;
		case '?': die = true; break;
		case 'v': asmlongopt::verbose++; break;
		case 'x': arg >> asmlongopt::maxChimeraLength; break;
        case 'i': arg >> asmlongopt::insertSize; break;
		case OPT_MAXEDGES: arg >> asmlongopt::maxEdges; break;
		case OPT_MAXINDEL: arg >> asmlongopt::maxIndelLength; break;
		case OPT_EXACT: asmlongopt::bExact = true; break;
		case OPT_HELP:
			std::cout << ASSEMBLE_USAGE_MESSAGE;
			exit(EXIT_SUCCESS);
		case OPT_VERSION:
			std::cout << ASSEMBLE_VERSION_MESSAGE;
			exit(EXIT_SUCCESS);

		}
	}

	// Build the output names
	asmlongopt::outContigsFile = asmlongopt::outFilePrefix + "-contigs.fa";
	asmlongopt::outVariantsFile = asmlongopt::outFilePrefix + "-variants.fa";
	asmlongopt::outGraphFile = asmlongopt::outFilePrefix + "-graph.asqg.gz";
	asmlongopt::outTrimFile  = asmlongopt::outFilePrefix + "-terminal.fa";
	asmlongopt::outSimpleBubblesFile = asmlongopt::outFilePrefix + "-simpeBubbles.fa";


	if (asmlongopt::prefix.empty() ||  asmlongopt::insertSize==0)
	{
		std::cerr << SUBPROGRAM ": missing mandatory arguments\n";
		die = true;
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

	if (die)
	{
		std::cout << "\n" << ASSEMBLE_USAGE_MESSAGE;
		exit(EXIT_FAILURE);
	}

	// Parse the input filename
	asmlongopt::asqgFile = argv[optind++];

	if (asmlongopt::maxChimeraLength ==0 )
	{
		asmlongopt::maxChimeraLength = asmlongopt::insertSize *2 ;
	}
}

void sequentialTrimAndSmooth (StringGraph* pGraph, size_t trimLength, bool bIsGapPrecent)
{
	pGraph->simplify();
	SGTrimVisitor trimVisit("",trimLength);
	SGSmoothingVisitor smoothingVisit(asmlongopt::maxIndelLength, asmlongopt::pBWT, bIsGapPrecent);

	if (pGraph->visit(trimVisit))
		pGraph->simplify();


	if (pGraph->visit(smoothingVisit))
	{
		pGraph->simplify();
		if (pGraph->visit(trimVisit))
			pGraph->simplify();

	}	
}
