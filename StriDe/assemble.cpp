//-----------------------------------------------
// Copyright 2014 National Chung Cheng University
// Written by Yao-Ting Huang
// Released under the GPL
//-----------------------------------------------
//
// assemble - convert read overlaps into contigs
//
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <map>
#include "Util.h"
#include "assemble.h"
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
"Copyright 2014 National Chung Cheng University, Taiwan\n";

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
"      -l, --min-branch-length=LEN      remove terminal branches only if they are less than LEN bases in length (default: 200)\n"
"  -v, --verbose                        display verbose output\n"
"      --help                           display this help and exit\n"
"\nOther minor control options:\n"
"      -o, --out-prefix=NAME            use NAME as the prefix of the output files (output files will be NAME-contigs.fa, etc)\n"
"      -m, --min-overlap=LEN            only use overlaps of at least LEN. This can be used to filter\n"
"          --transitive-reduction       remove transitive edges from the graph. Off by default.\n"
"          --max-edges=N                limit each vertex to a maximum of N edges. For highly repetitive regions\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";


namespace opt
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
	static bool bExact = true;

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
	{ "min-branch-length",     required_argument, NULL, 'l' },
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
int assembleMain(int argc, char** argv)
{
	Timer* pTimer = new Timer("StriDe assembly");
	parseAssembleOptions(argc, argv);

	std::cout << "\n#---  Parameters   ---#\n";
	std::cout << "Kmer Size              : " << opt::kmerLength << std::endl;
	std::cout << "Kmer Threshold         : " << opt::kmerThreshold << std::endl;
	std::cout << "Original Read Length   : " << opt::readLength << std::endl;
	std::cout << "Maximum Chimera Length : " << opt::maxChimeraLength << std::endl;
	std::cout << "Reliable Overlap Length: " << opt::credibleOverlapLength <<std::endl;
	std::cout << "Insert Size            : " << opt::insertSize << std::endl;

	assemble();
	delete pTimer;

	return 0;
}

int assemble()
{
	StringGraph* pGraph;
	#pragma omp parallel
	{
		#pragma omp single nowait
		{
			std::cout << "\n[ Loading string graph: " << opt::asqgFile <<  " ]\n";
			pGraph=SGUtil::loadASQGVertex(opt::asqgFile, opt::minOverlap, true, opt::maxEdges);
		}
		#pragma omp single nowait
		{
			std::cout << "[ Loading BWT ]\n";
			opt::pBWT = new BWT(opt::prefix + BWT_EXT, BWT::DEFAULT_SAMPLE_RATE_SMALL);
		}
		#pragma omp single nowait
		{
			std::cout << "[ Loading RBWT ]\n";
			opt::pRBWT = new BWT(opt::prefix + RBWT_EXT, BWT::DEFAULT_SAMPLE_RATE_SMALL);
		}
		#pragma omp single nowait
		{
			std::cout << "[ Loading SAI ]\n";
			opt::pSSA = new SampledSuffixArray(opt::prefix + SAI_EXT, SSA_FT_SAI);
		}
	}
    opt::indices.pBWT = opt::pBWT;
    opt::indices.pRBWT = opt::pRBWT;
    opt::indices.pSSA = opt::pSSA;
	
	pGraph=SGUtil::loadASQGEdge(opt::asqgFile, opt::minOverlap, true, opt::maxEdges, pGraph);

	if(opt::bExact)
		pGraph->setExactMode(true);
	//pGraph->printMemSize();

	// // Pre-assembly graph stats
	SGGraphStatsVisitor statsVisit;
	std::cout << "[Stats] Input graph:\n";
	pGraph->visitP(statsVisit);
	
	int phase = 0 ;

	// Remove containments from the graph
	std::cout << "Removing contained vertices from graph\n";
	SGContainRemoveVisitor containVisit;
	if(pGraph->hasContainment())
		pGraph->visit(containVisit);

	/*---Remove Transitive Edges---*/
	//std::cout << "Removing transitive edges\n";
	//SGTransitiveReductionVisitor trVisit;
	//pGraph->visit(trVisit);
	/*---Remove Transitive Edges---*/

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
	SGRemoveIllegalKmerEdgeVisitor ikeVisit(opt::pBWT,opt::kmerLength,opt::kmerThreshold, opt::credibleOverlapLength);
	std::cout << "\n[ Remove Illegal Kmer Edges due to kmerization ]" << std::endl;
	pGraph->visitP(ikeVisit);
	pGraph->simplify();

	/*** trim dead end vertices from small to large***/
	size_t trimLen = opt::kmerLength+1 , stepsize=opt::insertSize/5;
    while (trimLen < opt::insertSize *3/2)
    {
			SGTrimVisitor shortTrimVisit("",trimLen);
			std::cout << "[ Trimming short vertices (<" << trimLen  << ") ]\n";

			if(pGraph->visitP(shortTrimVisit))
				pGraph->simplify();

			 trimLen=trimLen+stepsize;
    }

	/*** Pop Bubbles ***/
	std::cout << "\n[ Remove bubbles and tips ]\n";
	graphTrimAndSmooth (pGraph, opt::maxChimeraLength);
	// outputGraphAndFasta(pGraph,"popBubbles",++phase);

	/*** Remove small chimeric vertices ***/
	std::cout << "\n[ Remove small chimera vertices ]\n";
	for (size_t threshold=2; threshold<=opt::kmerThreshold; threshold++)
		RemoveVertexWithBothShortEdges (pGraph, opt::readLength, opt::credibleOverlapLength, opt::pBWT, opt::kmerLength, threshold);

	//remove edges with overlap <= 1st overlap length peak (from unmerged original reads)
	RemoveVertexWithBothShortEdges (pGraph, opt::readLength, pGraph->getMinOverlap());
	RemoveVertexWithBothShortEdges (pGraph, opt::readLength, opt::credibleOverlapLength);
	RemoveVertexWithBothShortEdges (pGraph, opt::insertSize, opt::credibleOverlapLength);
	RemoveVertexWithBothShortEdges (pGraph, opt::maxChimeraLength, opt::credibleOverlapLength);
    
	pGraph->contigStats();
	pGraph->visit(statsVisit);
	outputGraphAndFasta(pGraph,"",++phase);

	   /*** remove edges from 1st peak to 2nd peak in the overlap length distribution
	1st peak: opt::credibleOverlapLength , 2nd peak: opt::insertSize*opt::minOverlapRatio 
    Define large vertex as opt::maxChimeraLength*4
    be very careful and little by little for avoiding misassembly ***/
    std::cout << "\n[ Remove chimeric edges with insufficient overlap or large-overlap difference from large vertices]\n";
	size_t stepsize2 = (opt::insertSize*opt::minOverlapRatio-opt::credibleOverlapLength)/4;	//iterate 4 times
    for (size_t len = opt::credibleOverlapLength ; len <= opt::insertSize*opt::minOverlapRatio ; len+=stepsize2 )	
	{
		SGRemoveByOverlapLenDiffVisitor srv1(1600, len, (opt::insertSize*opt::minOverlapRatio)+opt::credibleOverlapLength-len);	
		if(pGraph->visitP(srv1)) 
            graphTrimAndSmooth(pGraph, opt::maxChimeraLength);
    }

    // fix min overlap length= insertsize*opt::minOverlapRatio, remove more edges by overlap diff from large to small diff
	// iterate 2 times
    for (size_t stepsize3 = opt::credibleOverlapLength/4; stepsize3 <= opt::credibleOverlapLength/2 ; stepsize3+=stepsize3 )
    {
        SGRemoveByOverlapLenDiffVisitor srv1(1600, /*opt::insertSize*opt::minOverlapRatio*/ 0, opt::credibleOverlapLength-stepsize3);
        if(pGraph->visitP(srv1)) 
            graphTrimAndSmooth(pGraph, opt::maxChimeraLength);
    }

	// SGRemoveByOverlapLenDiffVisitor srv1(opt::maxChimeraLength*3, /*opt::insertSize*opt::minOverlapRatio*/ 0, opt::credibleOverlapLength/3);
    // if(pGraph->visitP(srv1)) 
            // graphTrimAndSmooth(pGraph, opt::maxChimeraLength);

	RemoveVertexWithBothShortEdges (pGraph, opt::readLength+100, opt::readLength*0.9);

	pGraph->contigStats(); 
	pGraph->visit(statsVisit);
	outputGraphAndFasta(pGraph,"",++phase);

	
	/*** Remove edges according to PE links, which remove small repeat vertices most of the time ***/
	for (size_t minPELink = 1; minPELink <= 1; minPELink++ )
    {
		// 51 is more accurate while 31 produces longer contigs
		SGRemoveEdgeByPEVisitor REPEVistit(opt::indices, opt::insertSize, 51, minPELink);
		if(pGraph->visitP(REPEVistit))
			graphTrimAndSmooth(pGraph, opt::maxChimeraLength);
    }

	pGraph->contigStats();
	pGraph->visit(statsVisit);
	outputGraphAndFasta(pGraph,"",++phase);

	
	/** Incrementally remove chimera edges from small to large chimeric vertices using overlap ratio 
	  Be very careful and little by little for avoiding misassembly 
	  This visitor resolves chimeric with size roughly equal to read length+100 (first peak)
	  Resolve larger vertices lead to misassembly ***/
    std::cout << "\n[ Remove chimeric edges from small vertices using overlap ratios ]\n";
	size_t stepsize4 = 15;
	for (size_t len = opt::readLength ; len <= opt::readLength+100 ; len+=stepsize4 )
		RemoveSmallOverlapRatioEdges ( pGraph, len);
	
	// Rename requires extra memory which should
	// only be done after the string graph has been greatly simplified.
	pGraph->renameVertices("");

	/******* Re-join broken islands/tips due to high-GC errors ********/
	size_t min_size_of_islandtip=opt::maxChimeraLength;

    /***************** 1. Trim bad ends of island/tip *****************/
	SGFastaErosionVisitor eFAVisit (opt::pBWT, opt::kmerLength, opt::kmerThreshold, min_size_of_islandtip);
	pGraph->visitP(eFAVisit);
	
	pGraph->contigStats();
	pGraph->visitP(statsVisit);
	outputGraphAndFasta(pGraph,"", ++phase);
    
    /*** 2. Collect read IDs mapped to large island/tip with size > min_size_of_islandtip ***/
	ThreadSafeListVector tslv;
	tslv.resize(opt::pSSA->getNumberOfReads());
    SGIslandCollectVisitor sgicv(&tslv, opt::indices, opt::insertSize, 51, min_size_of_islandtip);
    pGraph->visitP(sgicv);
    
	/*** 3. Join islands/tips with PE support using FM-index walk (depth,leaves,minoverlap)=(150, 2000, 19) ***/
	SGJoinIslandVisitor sgjiv(100, 4000, opt::kmerLength/2+4, min_size_of_islandtip, &tslv, opt::indices, 3);
	pGraph->visitProgress(sgjiv);
	graphTrimAndSmooth (pGraph, opt::maxChimeraLength, false);

	std::cout << "\n[Stats] Final graph statistic:\n";
	pGraph->contigStats();
	pGraph->visitP(statsVisit);
	
	// Write the results
	std::cout << "\n<Printing the contig file> : " << opt::outContigsFile << " \n" << std::endl;
	SGFastaVisitor av(opt::outContigsFile);
	pGraph->visit(av);
	//std::cout << "<Printing the ASQG and dot files>\n" << std::endl;
	pGraph->writeASQG(opt::outGraphFile);
    pGraph->writeDot("StriDe-graph.dot",0);

	delete pGraph;
	return 0;

}

//
// Handle command line arguments
//
void parseAssembleOptions(int argc, char** argv)
{
	optind=1;	//reset getopt

	// Set defaults
	bool die = false;
	for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
	{
		std::istringstream arg(optarg != NULL ? optarg : "");
		switch (c)
		{
		case 'o': arg >> opt::outFilePrefix; break;
		case 'p': arg >> opt::prefix; break;
		case 'm': arg >> opt::minOverlap; break;
		case '?': die = true; break;
		case 'v': opt::verbose++; break;
		case 'k': arg >> opt::kmerLength; break;
		case 't': arg >> opt::kmerThreshold; break;
		case 'r': arg >> opt::readLength; break;
		case 'T': arg >> opt::minOverlapRatio; break;
		case 'x': arg >> opt::maxChimeraLength; break;
		case 'c': arg >> opt::credibleOverlapLength; break;
        case 'i': arg >> opt::insertSize; break;
		case OPT_MAXEDGES: arg >> opt::maxEdges; break;
		case OPT_MAXINDEL: arg >> opt::maxIndelLength; break;
		case OPT_EXACT: opt::bExact = true; break;
		case OPT_HELP:
			std::cout << ASSEMBLE_USAGE_MESSAGE;
			exit(EXIT_SUCCESS);
		case OPT_VERSION:
			std::cout << ASSEMBLE_VERSION_MESSAGE;
			exit(EXIT_SUCCESS);

		}
	}

	// Build the output names
	opt::outContigsFile = opt::outFilePrefix + "-contigs.fa";
	opt::outVariantsFile = opt::outFilePrefix + "-variants.fa";
	opt::outGraphFile = opt::outFilePrefix + "-graph.asqg.gz";
	opt::outTrimFile  = opt::outFilePrefix + "-terminal.fa";
	opt::outSimpleBubblesFile = opt::outFilePrefix + "-simpeBubbles.fa";


	if (opt::prefix.empty() || opt::readLength==0 || opt::insertSize==0)
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
	opt::asqgFile = argv[optind++];

	if (opt::maxChimeraLength ==0 )
	{
		opt::maxChimeraLength = opt::insertSize *2 ;
	}
	if (opt::credibleOverlapLength ==0 )
	{
		opt::credibleOverlapLength = opt::readLength * opt::minOverlapRatio ;
	}
}


void graphTrimAndSmooth (StringGraph* pGraph, size_t trimLength, bool bIsGapPrecent)
{
	pGraph->simplify();
	SGTrimVisitor trimVisit("",trimLength);
	SGSmoothingVisitor smoothingVisit(opt::maxIndelLength, opt::pBWT, bIsGapPrecent);
	//SGSimpleBubbleVisitor sbVisit(opt::pBWT,opt::kmerLength,opt::maxBubbleGapDivergence, opt::maxBubbleDivergence, opt::maxIndelLength);

	if (pGraph->visitP(trimVisit))
		pGraph->simplify();


	if (pGraph->visitP(smoothingVisit))
	{
		pGraph->simplify();
		if (pGraph->visitP(trimVisit))
			pGraph->simplify();

	}

	// SGDuplicateVisitor sgdup;
	// if(pGraph->visit(sgdup))
		// pGraph->simplify();
				
	// if (pGraph->visit(sbVisit))
	// {
		// pGraph->simplify();
		// if (pGraph->visitP(trimVisit))
			// pGraph->simplify();
	// }
	
}

void RemoveVertexWithBothShortEdges (StringGraph* pGraph ,size_t vertexLength ,size_t overlapLength, BWT* pBWT , size_t kmerLength, float threshold )
{
	if (pBWT !=NULL)
		assert (kmerLength>0 && threshold>0);

    /* debug
	std::stringstream filename;
	if (pBWT!=NULL)
		filename << opt::outFilePrefix  << "-bothShortEdgeVertexL" << vertexLength << "O" << overlapLength << "K" << kmerLength << "T" << threshold << ".fa";
	else
		filename << opt::outFilePrefix  <<  "-bothShortEdgeVertexL" << vertexLength << "O" <<overlapLength << ".fa";

	SGBothShortEdgesRemoveVisitor serVisit(vertexLength,overlapLength,filename.str(),pBWT,kmerLength,threshold);
    */
    SGBothShortEdgesRemoveVisitor serVisit(vertexLength,overlapLength,pBWT,kmerLength,threshold);
	if(pGraph->visitP(serVisit))
        graphTrimAndSmooth (pGraph, opt::maxChimeraLength);

}


void RemoveSmallOverlapRatioEdges ( StringGraph* pGraph, size_t chimeraLength)
{
	size_t overlapLength = chimeraLength * opt::minOverlapRatio ;

    /* debug
	std::stringstream sweepSS;
	sweepSS << "R" << opt::minOverlapRatio << "M" << overlapLength;
	std::string sweepFilename = opt::outFilePrefix + "-sweepLowOverlapRatioEdge" + sweepSS.str() + ".log";
	SGLowOverlapRatioEdgeSweepVisitor sweepLowRatioEdgesVisit(sweepFilename,opt::minOverlapRatio,overlapLength);
	*/
	SGLowOverlapRatioEdgeSweepVisitor sweepLowRatioEdgesVisit(chimeraLength, opt::minOverlapRatio, overlapLength);
	if(pGraph->visitP(sweepLowRatioEdgesVisit))
        graphTrimAndSmooth (pGraph, opt::maxChimeraLength);

}

void outputGraphAndFasta(StringGraph* pGraph , std::string  name , int phase)
{
	std::cout << "\n<Printing the fasta & ASQG file>" << std::endl;
	std::cout << "< name: " << name << " phase: " << phase << " >\n"<<  std::endl;

	std::string out = "phase";

	if (phase>0)
	{
		std::stringstream  ss ;
		ss << phase ;
		out = out + ss.str() + "-";
	}
	
	if(!name.empty()) out = out + name+ "-";

	SGFastaVisitor fastaVisit (out+"contigs.fa");
	pGraph->visit(fastaVisit);
	// pGraph->writeASQG(out+"graph.asqg.gz");
}
