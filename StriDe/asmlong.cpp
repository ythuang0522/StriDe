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
#include "SAIntervalTree.h"
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
"      -r, --read-length=LEN            original read length\n"
"      -i,  --insert-size               insert size of the paired-end library\n"
"      -p, --prefix=NAME               prefix of FM-index of paired-end reads (bwt, rbwt, sai, rsai)\n"
"\nSimplify graph parameters:\n"
"      -k, --kmer-size=N                The length of the kmer to use. (default: 31)\n"
"      -t, --kmer-threshold=N           filter average kmer frequency vertex less than N (default: 3)\n"
"      -x, --max-chimera=LEN            maximum chimera length (default: read length(R)*2 )\n"
"      -c, --credible-overlap=LEN       credible overlap length (default: 80) \n"
"      -T, --min-overlap-ratio=double      min-overlap-ratio in double (default: 0.8)\n"
"  -v, --verbose                        display verbose output\n"
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

	//kmer frequence parameters
	static size_t kmerLength = 31;
	static size_t kmerThreshold = 3;


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
	static size_t readLength = 0 ;
	static double minOverlapRatio=0.8;
	static size_t insertSize=0;

	static size_t credibleOverlapLength = 0;
	static size_t maxChimeraLength = 0;

}

static const char* shortopts = "k:t:p:o:m:i:r:T:x:c:v";

enum { OPT_HELP = 1, OPT_VERSION, OPT_EXACT, OPT_MAXINDEL,OPT_MAXEDGES};

static const struct option longopts[] = {
	{ "verbose",               no_argument,       NULL, 'v' },
	{ "out-prefix",            required_argument, NULL, 'o' },
	{ "prefix",                required_argument, NULL, 'p' },
	{ "min-overlap",           required_argument, NULL, 'm' },
	{ "min-overlap-ratio",           required_argument, NULL, 'T' },
	{ "max-indel",             required_argument, NULL, OPT_MAXINDEL },
	{ "max-edges",             required_argument, NULL, OPT_MAXEDGES },
	{ "kmer-length",           required_argument, NULL, 'k' },
	{ "kmer-threshold",        required_argument, NULL, 't' },
	{ "read-length",           required_argument, NULL, 'r' },
	{ "max-chimera",           required_argument, NULL, 'x' },
	{ "credible-overlap",      required_argument, NULL, 'c' },
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
	std::cout << "Kmer Size              : " << asmlongopt::kmerLength << std::endl;
	std::cout << "Kmer Threshold         : " << asmlongopt::kmerThreshold << std::endl;
	std::cout << "Original Read Length   : " << asmlongopt::readLength << std::endl;
	std::cout << "Maximum Chimera Length : " << asmlongopt::maxChimeraLength << std::endl;
	std::cout << "Reliable Overlap Length: " << asmlongopt::credibleOverlapLength <<std::endl;
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
	
	int phase = 0 ;

	// Remove containments from the graph
	std::cout << "Removing contained vertices from graph\n";
	SGContainRemoveVisitor containVisit;
	while(pGraph->hasContainment())
		pGraph->visit(containVisit);

	/*---Remove Transitive Edges---*/
	std::cout << "Removing transitive edges\n";
	SGTransitiveReductionVisitor trVisit;
	pGraph->visit(trVisit);
	/*---Remove Transitive Edges---*/
	// getchar();

	// Compact together unbranched chains of vertices
	std::cout << "Start to simplify unipaths ...\n";
    pGraph->simplify();
	std::cout << "[Stats] Simplified graph:\n";
	pGraph->visitP(statsVisit);


	/**********************Compute overlap raio and diff (for debug)**************
	std::ofstream ssol  ("simpleOverlapLength.histo", std::ofstream::out);
	std::map<size_t,int> simpleStats = pGraph->getCountMap() ;
	//Stats simple overlap length
	std::cout << "< Compute Stats of Overlap Length and Ratios>" << std::endl;
	//ssol << "Length\tCount" << std::endl;
	for (std::map<size_t,int> ::iterator iter=simpleStats.begin(); iter!=simpleStats.end(); ++iter)
		ssol << iter->first << "\t" << iter->second << std::endl;
	ssol.close();

	std::ofstream orStats  ("overlapRaito.stats", std::ofstream::out);
	std::vector < overlapRatioInfo > overlapRatioStats = pGraph->getOverlapRatioStats();
	for (std::vector < overlapRatioInfo >::iterator iter=overlapRatioStats.begin(); iter!=overlapRatioStats.end(); ++iter)
		orStats << (*iter).overlapLength  << "\t" << (*iter).originalLength << "\t"
		        << std::setiosflags(std::ios::fixed) << std::setprecision(3)
				<< (double)(*iter).overlapLength/(*iter).originalLength << "\t" << (*iter).ratioDiff << std::endl;

	orStats.close();
	**********************Compute overlap raio and diff***************************/
		
	// /*** Remove Illegal Kmer Edges ***/
	// SGRemoveIllegalKmerEdgeVisitor ikeVisit(asmlongopt::pBWT,asmlongopt::kmerLength,asmlongopt::kmerThreshold, asmlongopt::credibleOverlapLength);
	// std::cout << "\n[ Remove Illegal Kmer Edges due to kmerization ]" << std::endl;
	// pGraph->visitP(ikeVisit);
	// pGraph->simplify();

	// outputGraphAndFasta(pGraph, "debug", 0);
	// pGraph->writeDot("debug.dot",0);
	// getchar();
	/*** trim dead end vertices from small to large***/
	size_t trimLen = asmlongopt::kmerLength+1 , stepsize=asmlongopt::insertSize/5;
    while (trimLen < asmlongopt::insertSize)
    {
			SGTrimVisitor shortTrimVisit("",trimLen);
			std::cout << "[ Trimming short vertices (<" << trimLen  << ") ]\n";

			if(pGraph->visitP(shortTrimVisit))
				pGraph->simplify();

			 trimLen=trimLen+stepsize;
    }

	/*** Pop Bubbles ***/
	std::cout << "\n[ Remove bubbles and tips ]\n";
	graphTrimAndSmooth (pGraph, asmlongopt::maxChimeraLength);
	// outputGraphAndFasta(pGraph,"popBubbles",++phase);

	/*** Remove small chimeric vertices ***/
	std::cout << "\n[ Remove small chimera vertices ]\n";

	//remove edges with overlap <= 1st overlap length peak (from unmerged original reads)
	RemoveVertexWithBothShortEdges (pGraph, asmlongopt::readLength, pGraph->getMinOverlap());
	RemoveVertexWithBothShortEdges (pGraph, asmlongopt::readLength, asmlongopt::credibleOverlapLength);
	RemoveVertexWithBothShortEdges (pGraph, asmlongopt::insertSize, asmlongopt::credibleOverlapLength);
	RemoveVertexWithBothShortEdges (pGraph, asmlongopt::maxChimeraLength, asmlongopt::readLength-1);
    
	pGraph->contigStats();
	pGraph->visit(statsVisit);
	outputGraphAndFasta(pGraph,"",++phase);

	   /*** remove edges from 1st peak to 2nd peak in the overlap length distribution
	1st peak: asmlongopt::credibleOverlapLength , 2nd peak: asmlongopt::insertSize*asmlongopt::minOverlapRatio 
    Define large vertex as asmlongopt::maxChimeraLength*4
    be very careful and little by little for avoiding misassembly ***/
    std::cout << "\n[ Remove chimeric edges with insufficient overlap or large-overlap difference from large vertices]\n";
	size_t stepsize2 = (asmlongopt::insertSize*asmlongopt::minOverlapRatio-asmlongopt::credibleOverlapLength)/4;	//iterate 4 times
    for (size_t len = asmlongopt::credibleOverlapLength ; len <= asmlongopt::insertSize*asmlongopt::minOverlapRatio ; len+=stepsize2 )	
	{
		SGRemoveByOverlapLenDiffVisitor srv1(1600, len, (asmlongopt::insertSize*asmlongopt::minOverlapRatio)+asmlongopt::credibleOverlapLength-len);	
		if(pGraph->visitP(srv1)) 
            graphTrimAndSmooth(pGraph, asmlongopt::maxChimeraLength);
    }

	pGraph->contigStats(); 
	pGraph->visit(statsVisit);
	outputGraphAndFasta(pGraph,"",++phase);
	
	/** Incrementally remove chimera edges from small to large chimeric vertices using overlap ratio 
	  Be very careful and little by little for avoiding misassembly 
	  This visitor resolves chimeric with size roughly equal to read length+100 (first peak)
	  Resolve larger vertices lead to misassembly ***/
    std::cout << "\n[ Remove chimeric edges from small vertices using overlap ratios ]\n";
	SGLowOverlapRatioEdgeSweepVisitor sweepLowRatioEdgesVisit(1000000, asmlongopt::minOverlapRatio, asmlongopt::maxChimeraLength);
	if(pGraph->visitP(sweepLowRatioEdgesVisit))
        graphTrimAndSmooth (pGraph, asmlongopt::maxChimeraLength);

	RemoveVertexWithBothShortEdges (pGraph, asmlongopt::maxChimeraLength*4, asmlongopt::maxChimeraLength*0.8);

	graphTrimAndSmooth (pGraph, asmlongopt::maxChimeraLength*2);

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
		case 'k': arg >> asmlongopt::kmerLength; break;
		case 't': arg >> asmlongopt::kmerThreshold; break;
		case 'r': arg >> asmlongopt::readLength; break;
		case 'T': arg >> asmlongopt::minOverlapRatio; break;
		case 'x': arg >> asmlongopt::maxChimeraLength; break;
		case 'c': arg >> asmlongopt::credibleOverlapLength; break;
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


	if (asmlongopt::prefix.empty() || asmlongopt::readLength==0 || asmlongopt::insertSize==0)
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
	if (asmlongopt::credibleOverlapLength ==0 )
	{
		asmlongopt::credibleOverlapLength = asmlongopt::readLength * asmlongopt::minOverlapRatio ;
	}
}


