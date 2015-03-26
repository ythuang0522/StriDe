//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// overlap - compute pairwise overlaps between reads
//
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <iomanip>
#include "Util.h"
#include "overlap.h"
#include "SuffixArray.h"
#include "BWT.h"
#include "SGACommon.h"
#include "OverlapCommon.h"
#include "Timer.h"
#include "BWTAlgorithms.h"
#include "ASQG.h"
#include "gzstream.h"
#include "SequenceProcessFramework.h"
#include "OverlapProcess.h"
#include "ReadInfoTable.h"
#include <sys/stat.h>

/*Tatsuki include */
#include <tr1/unordered_map>
#include <tr1/unordered_set>
#define HashMap std::tr1::unordered_map
#define HashSet std::tr1::unordered_set
typedef std::tr1::hash<std::string> StringHasher;
#include <sparsehash/dense_hash_set>
#define DenseHashSet google::dense_hash_set
typedef DenseHashSet<std::string, StringHasher> ReadHashSet;
typedef ReadHashSet::iterator ReadHashSetIter;
#include <time.h>
#include <omp.h>
/*Tatsuki include */





//
enum OutputType
{
	OT_ASQG,
	OT_RAW
};

// Functions
size_t computeHitsSerial(const std::string& prefix, const std::string& readsFile, const OverlapAlgorithm* pOverlapper, int minOverlap, StringVector& filenameVec, std::ostream* pASQGWriter);

size_t computeHitsParallel(int numThreads, const std::string& prefix, const std::string& readsFile, const OverlapAlgorithm* pOverlapper, int minOverlap, StringVector& filenameVec, std::ostream* pASQGWriter);

//
void convertHitsToASQG(const std::string& indexPrefix, const StringVector& hitsFilenames, std::ostream* pASQGWriter, DenseHashSet < std::string, StringHasher > *SuperRepeatVertices);

void convertHitsToASQGedgeParallel(const std::string hitsFilename, std::ostream* pASQGWriter ,std::ostream* pDiscardEdgeWriter, 
	SuffixArray* pFwdSAI ,SuffixArray* pRevSAI , ReadInfoTable* pQueryRIT , ReadInfoTable* pTargetRIT, DenseHashSet < std::string, StringHasher > *SuperRepeatVertices);

//
// Getopt
//
#define SUBPROGRAM "overlap"
static const char *OVERLAP_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2009 Wellcome Trust Sanger Institute\n";

static const char *OVERLAP_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... READSFILE\n"
"Compute pairwise overlap between all the sequences in READS\n"
"\n"
"      --help                           display this help and exit\n"
"      -v, --verbose                    display verbose output\n"
"      -t, --threads=NUM                use NUM worker threads to compute the overlaps (default: no threading)\n"
"      -e, --error-rate                 the maximum error rate allowed to consider two sequences aligned (default: exact matches only)\n"
"      -m, --min-overlap=LEN            minimum overlap required between two reads (default: 45)\n"
"      -f, --target-file=FILE           perform the overlap queries against the reads in FILE\n"
"      -p, --paired-overlap             output only paired overlaps\n"
"      -x, --exhaustive                 output all overlaps, including transitive edges\n"
"          --exact                      force the use of the exact-mode irreducible block algorithm. This is faster\n"
"                                       but requires that no substrings are present in the input set.\n"
"      -l, --seed-length=LEN            force the seed length to be LEN. By default, the seed length in the overlap step\n"
"                                       is calculated to guarantee all overlaps with --error-rate differences are found.\n"
"                                       This option removes the guarantee but will be (much) faster. As SGA can tolerate some\n"
"                                       missing edges, this option may be preferable for some data sets.\n"
"      -s, --seed-stride=LEN            force the seed stride to be LEN. This parameter will be ignored unless --seed-length\n"
"                                       is specified (see above). This parameter defaults to the same value as --seed-length\n"
"      -d, --sample-rate=N              sample the symbol counts every N symbols in the FM-index. Higher values use significantly\n"
"                                       less memory at the cost of higher runtime. This value must be a power of 2 (default: 128)\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* PROGRAM_IDENT =
PACKAGE_NAME "::" SUBPROGRAM;

namespace opt
{
	static unsigned int verbose;
	static int numThreads = 1;
	static OutputType outputType = OT_ASQG;
	static std::string readsFile;
	static std::string targetFile;
	static std::string outFile;
	
	static double errorRate = 0.0f;
	static unsigned int minOverlap = DEFAULT_MIN_OVERLAP;
	static int seedLength = 0;
	static int seedStride = 0;
	static int sampleRate = BWT::DEFAULT_SAMPLE_RATE_SMALL;
	static bool bIrreducibleOnly = true;
	static bool bExactIrreducible = false;
	static bool bIsPairedOverlapOnly  = false;
}

static const char* shortopts = "m:d:e:t:l:s:o:f:vixp";

enum { OPT_HELP = 1, OPT_VERSION, OPT_EXACT };

static const struct option longopts[] = {
	{ "verbose",     no_argument,       NULL, 'v' },
	{ "threads",     required_argument, NULL, 't' },
	{ "min-overlap", required_argument, NULL, 'm' },
	{ "sample-rate", required_argument, NULL, 'd' },
	{ "outfile",     required_argument, NULL, 'o' },
	{ "target-file", required_argument, NULL, 'f' },
	{ "prefix",      required_argument, NULL, 'p' },
	{ "error-rate",  required_argument, NULL, 'e' },
	{ "seed-length", required_argument, NULL, 'l' },
	{ "seed-stride", required_argument, NULL, 's' },
	{ "exhaustive",  no_argument,       NULL, 'x' },
	{ "paired-overlap",no_argument,     NULL, 'p' },
	{ "exact",       no_argument,       NULL, OPT_EXACT },
	{ "help",        no_argument,       NULL, OPT_HELP },
	{ "version",     no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

//
// Main
//
int overlapMain(int argc, char** argv)
{
	parseOverlapOptions(argc, argv);

	// Prepare the output ASQG file
	assert(opt::outputType == OT_ASQG);

	// Open output file
	std::ostream* pASQGWriter = createWriter(opt::outFile);

	// Build and write the ASQG header
	ASQG::HeaderRecord headerRecord;
	headerRecord.setOverlapTag(opt::minOverlap);
	headerRecord.setErrorRateTag(opt::errorRate);
	headerRecord.setInputFileTag(opt::readsFile);
	headerRecord.setContainmentTag(false); // containments are always present
	headerRecord.setTransitiveTag(!opt::bIrreducibleOnly);
	headerRecord.write(*pASQGWriter);

	// Compute the overlap hits
	StringVector hitsFilenames;
	//StringVector EdgeRecordFilenames(opt::numThreads);
	//StringVector discardEdgeRecordFilenames(opt::numThreads);
	
	// Determine which index files to use. If a target file was provided,
	// use the index of the target reads
	std::string indexPrefix;
	if(!opt::targetFile.empty())
		indexPrefix = stripFilename(opt::targetFile);
	else
		indexPrefix = stripFilename(opt::readsFile);

	BWT *pBWT, *pRBWT;
	SuffixArray *pFwdSAI, *pRevSAI;
	#pragma omp parallel
	{
		#pragma omp single nowait
		{
			pBWT = new BWT(indexPrefix + BWT_EXT, opt::sampleRate);
		}
		#pragma omp single nowait
		{
			pRBWT = new BWT(indexPrefix + RBWT_EXT, opt::sampleRate);
		}
		#pragma omp single nowait
		{
			pFwdSAI = new SuffixArray(indexPrefix + SAI_EXT);
		}
		#pragma omp single nowait
		{
			pRevSAI = new SuffixArray(indexPrefix + RSAI_EXT);
		}
	}
	// Load the ReadInfoTable for the queries to look up the ID and lengths of the hits
	ReadInfoTable* pQueryRIT = new ReadInfoTable(opt::readsFile);
	// If the target file is not the query file, load its ReadInfoTable
	ReadInfoTable* pTargetRIT;
	if(!opt::targetFile.empty() && opt::targetFile != opt::readsFile)
		pTargetRIT = new ReadInfoTable(opt::targetFile);
	else
		pTargetRIT = pQueryRIT;

	OverlapAlgorithm* pOverlapper = new OverlapAlgorithm(pBWT, pRBWT, pFwdSAI, pRevSAI, pQueryRIT, pTargetRIT);
	
	Timer* pTimer = new Timer(PROGRAM_IDENT);

	// Make a prefix for the hit edges files
	std::string outPrefix;
	outPrefix = stripFilename(opt::readsFile);
	if(!opt::targetFile.empty())
	{
		outPrefix.append(1, '.');
		outPrefix.append(stripFilename(opt::targetFile));
	}
	
	time_t now = time(NULL);	
	std::cout << "\n# start time of overlapping: " << asctime(localtime(&now))<<std::endl;
	
	if(opt::numThreads <= 1)
	{
		printf("[%s] starting serial-mode overlap computation\n", PROGRAM_IDENT);
		computeHitsSerial(outPrefix, opt::readsFile, pOverlapper, opt::minOverlap, hitsFilenames, pASQGWriter);
	}
	else
	{
		printf("[%s] starting parallel-mode overlap computation with %d threads\n", PROGRAM_IDENT, opt::numThreads);
		computeHitsParallel(opt::numThreads, outPrefix, opt::readsFile, pOverlapper, opt::minOverlap, hitsFilenames, pASQGWriter);
	}

	delete pOverlapper;
	delete pBWT; 
	delete pRBWT;
	delete pASQGWriter;
	delete pTimer;

	return 0;
}

// Compute the hits for each read in the input file without threading
// Return the number of reads processed
size_t computeHitsSerial(const std::string& prefix, const std::string& readsFile, const OverlapAlgorithm* pOverlapper, int minOverlap, 
						StringVector& filenameVec, std::ostream* pASQGWriter)
{
	std::string filename = prefix + "-thread0" +HITS_EXT + GZIP_EXT;
	filenameVec.push_back(filename);

	OverlapProcess processor(filename, pOverlapper, minOverlap);
	OverlapPostProcess postProcessor(pASQGWriter, pOverlapper);

	size_t numProcessed = 
	SequenceProcessFramework::processSequencesSerial<SequenceWorkItem,
	OverlapResult, 
	OverlapProcess, 
	OverlapPostProcess>(readsFile, &processor, &postProcessor);
	return numProcessed;
}

// Compute the hits for each read in the SeqReader file with threading
// The way this works is we create a vector of numThreads OverlapProcess pointers and 
// pass this to the SequenceProcessFramework which wraps the processes
// in threads and distributes the reads to each thread.
// The number of reads processsed is returned
size_t computeHitsParallel(int numThreads, const std::string& prefix, const std::string& readsFile, 
											const OverlapAlgorithm* pOverlapper, int minOverlap, 
											StringVector& filenameVec, std::ostream* pASQGWriter)
{
	//std::string filename = prefix + HITS_EXT + GZIP_EXT;

	std::vector<OverlapProcess*> processorVector;
	for(int i = 0; i < numThreads; ++i)
	{
		std::stringstream ss;
		ss << prefix << "-thread" << i << HITS_EXT << GZIP_EXT;
		std::string outfile = ss.str();
		filenameVec.push_back(outfile);
		OverlapProcess* pProcessor = new OverlapProcess(outfile, pOverlapper, minOverlap);
		processorVector.push_back(pProcessor);
	}

	//Remove previous edge files generated by larger threads, otherwise, subsequent assembly may load inconsistent edges files
	std::stringstream ss;
	size_t fileIdx=numThreads;
	ss << prefix << "-thread" << fileIdx << HITS_EXT << GZIP_EXT;
	std::string edgefile = ss.str();
	
	struct stat buffer;   
	while(stat (edgefile.c_str(), &buffer) == 0)
	{
		remove(edgefile.c_str());
		//search for next edge file
		ss.str("");
		fileIdx++;
		ss << prefix << "-thread" << fileIdx << HITS_EXT << GZIP_EXT;
		edgefile = ss.str();
	}
	/*
	std::istream* pistream = createReader(edgefile);
	while(pistream->peek()!=std::istream::traits_type::eof())
	{
		//std::cout << edgefile << std::endl;
		unlink(edgefile.c_str());
		//search for next edge file
		ss.str("");
		fileIdx++;
		ss << prefix << "-thread" << fileIdx << HITS_EXT << GZIP_EXT;
		edgefile = ss.str();
		pistream = createReader(edgefile);		
	}*/

	
	// The post processing is performed serially so only one post processor is created
	OverlapPostProcess postProcessor(pASQGWriter, pOverlapper);
	
	size_t numProcessed = 
	SequenceProcessFramework::processSequencesParallel<SequenceWorkItem,
	OverlapResult, 
	OverlapProcess, 
	OverlapPostProcess>(readsFile, processorVector, &postProcessor);

	for(int i = 0; i < numThreads; ++i)
		delete processorVector[i];

	return numProcessed;
}

//
// void convertHitsToASQG(const std::string& indexPrefix, const StringVector& hitsFilenames, std::ostream* pASQGWriter, 
										// DenseHashSet < std::string, StringHasher > *SuperRepeatVertices)
// {
	// // Load the suffix array index and the reverse suffix array index
	// // Note these are not the full suffix arrays
	// SuffixArray* pFwdSAI = new SuffixArray(indexPrefix + SAI_EXT);
	// SuffixArray* pRevSAI = new SuffixArray(indexPrefix + RSAI_EXT);

	// // Load the ReadInfoTable for the queries to look up the ID and lengths of the hits
	// ReadInfoTable* pQueryRIT = new ReadInfoTable(opt::readsFile);

	// // If the target file is not the query file, load its ReadInfoTable
	// ReadInfoTable* pTargetRIT;
	// if(!opt::targetFile.empty() && opt::targetFile != opt::readsFile)
	// pTargetRIT = new ReadInfoTable(opt::targetFile);
	// else
	// pTargetRIT = pQueryRIT;

	// bool bIsSelfCompare = pTargetRIT == pQueryRIT;
	
	// // Convert the hits to overlaps and write them to the asqg file as initial edges
	// for(StringVector::const_iterator iter = hitsFilenames.begin(); iter != hitsFilenames.end(); ++iter)
	// {
		// printf("[%s] parsing file %s\n", PROGRAM_IDENT, iter->c_str());
		// std::istream* pReader = createReader(*iter);
		
		// // Read each hit sequentially, converting it to an overlap
		// std::string line;
		// while(getline(*pReader, line))
		// {
			// size_t readIdx;
			// size_t totalEntries;
			// bool isSubstring;
			// OverlapVector ov;
			// OverlapCommon::parseHitsString(line, pQueryRIT, pTargetRIT, pFwdSAI, pRevSAI, bIsSelfCompare, readIdx, totalEntries, ov, isSubstring, SuperRepeatVertices);
			// for(OverlapVector::iterator iter = ov.begin(); iter != ov.end(); ++iter)
			// {
				// ASQG::EdgeRecord edgeRecord(*iter);
				// edgeRecord.write(*pASQGWriter);
			// }
		// }
		// delete pReader;

		// // delete the hits file
		// unlink(iter->c_str());
	// }

	// // Deallocate data
	// if(pTargetRIT != pQueryRIT)
	// delete pTargetRIT;
	// delete pFwdSAI;
	// delete pRevSAI;
	// delete pQueryRIT;
// }

/************************************************************************************************************************************************/

// void convertHitsToASQGedgeParallel(const std::string hitsFilename, std::ostream* pASQGWriter ,std::ostream* pDiscardEdgeWriter
	// ,SuffixArray* pFwdSAI ,SuffixArray* pRevSAI , ReadInfoTable* pQueryRIT , ReadInfoTable* pTargetRIT, DenseHashSet < std::string, StringHasher > *SuperRepeatVertices)
// {


	// bool bIsSelfCompare = pTargetRIT == pQueryRIT;

	// // Convert the hits to overlaps and write them to the asqg file as initial edges

	// printf("[%s] parsing file %s\n", PROGRAM_IDENT, hitsFilename.c_str());
	// std::istream* pReader = createReader(hitsFilename);

	// // paired-overlap mode

	// std::string line,line2;
	// while(getline(*pReader, line) && getline(*pReader, line2) )
	// {
		// size_t readIdx;
		// size_t totalEntries;
		// bool isSubstring;

		// OverlapVector ov1;
		// OverlapCommon::parseHitsString(line, pQueryRIT, pTargetRIT, pFwdSAI, pRevSAI, bIsSelfCompare, readIdx, totalEntries, ov1, isSubstring, SuperRepeatVertices);

		// OverlapVector ov2;
		// OverlapCommon::parseHitsString(line2, pQueryRIT, pTargetRIT, pFwdSAI, pRevSAI, bIsSelfCompare, readIdx, totalEntries, ov2, isSubstring, SuperRepeatVertices);

		// ReadHashSet ov2_set;
		// ReadHashSet intersection;
		// ov2_set.set_empty_key("");
		// intersection.set_empty_key("");
		
		// for(OverlapVector::iterator iter = ov2.begin(); iter != ov2.end(); ++iter)
		// {
			// Overlap r2=*iter;
			// if (!r2.match.isReverse) 
			// ov2_set.insert (r2.id[1].substr(0,r2.id[1].find_last_of("/")));
		// }			

		
		// for(OverlapVector::iterator iter = ov1.begin(); iter != ov1.end(); ++iter)
		// {
			// ASQG::EdgeRecord edgeRecord(*iter);
			// Overlap r1=*iter;
			// if (r1.match.isReverse) {  edgeRecord.write (*pDiscardEdgeWriter); continue; }
			// std::string mainName=r1.id[1].substr(0,r1.id[1].find_last_of("/"));
			// ReadHashSetIter RHSiter = ov2_set.find (mainName);
			// if (RHSiter != ov2_set.end()) 
			// {
				// intersection.insert (mainName);	
				// edgeRecord.write(*pASQGWriter);	
			// }
			// else edgeRecord.write (*pDiscardEdgeWriter);
		// }
		

		// for(OverlapVector::iterator iter = ov2.begin(); iter != ov2.end(); ++iter)
		// {
			// ASQG::EdgeRecord edgeRecord(*iter);
			// Overlap r2=*iter;
			// if (r2.match.isReverse) { edgeRecord.write (*pDiscardEdgeWriter); continue; }
			// ReadHashSetIter RHSiter = intersection.find (r2.id[1].substr(0,r2.id[1].find_last_of("/")));
			// if (RHSiter != intersection.end()) edgeRecord.write(*pASQGWriter);
			// else edgeRecord.write (*pDiscardEdgeWriter);
		// }
		
	// }
	

	// delete pReader;

	// // delete the hits file
	// unlink(hitsFilename.c_str());
	
// }



/*************************************************************************************************************************************************/
// 
// Handle command line arguments
//
void parseOverlapOptions(int argc, char** argv)
{
	optind=1;	//reset getopt
	bool die = false;
	for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
	{
		std::istringstream arg(optarg != NULL ? optarg : "");
		switch (c) 
		{
		case 'm': arg >> opt::minOverlap; break;
		case 'o': arg >> opt::outFile; break;
		case 'e': arg >> opt::errorRate; break;
		case 't': arg >> opt::numThreads; break;
		case 'l': arg >> opt::seedLength; break;
		case 's': arg >> opt::seedStride; break;
		case 'd': arg >> opt::sampleRate; break;
		case 'f': arg >> opt::targetFile; break;
		case OPT_EXACT: opt::bExactIrreducible = true; break;
		case 'x': opt::bIrreducibleOnly = false; break;
		case 'p': opt::bIsPairedOverlapOnly = true;  opt::bIrreducibleOnly = false; break;
		case '?': die = true; break;
		case 'v': opt::verbose++; break;
		case OPT_HELP:
			std::cout << OVERLAP_USAGE_MESSAGE;
			exit(EXIT_SUCCESS);
		case OPT_VERSION:
			std::cout << OVERLAP_VERSION_MESSAGE;
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

	if(!IS_POWER_OF_2(opt::sampleRate))
	{
		std::cerr << SUBPROGRAM ": invalid parameter to -d/--sample-rate, must be power of 2. got: " << opt::sampleRate << "\n";
		die = true;
	}

	if (die) 
	{
		std::cout << "\n" << OVERLAP_USAGE_MESSAGE;
		exit(EXIT_FAILURE);
	}

	// Validate parameters
	if(opt::errorRate <= 0)
	opt::errorRate = 0.0f;
	
	if(opt::seedLength < 0)
	opt::seedLength = 0;

	if(opt::seedLength > 0 && opt::seedStride <= 0)
	opt::seedStride = opt::seedLength;
	
	// Parse the input filenames
	opt::readsFile = argv[optind++];

	if(opt::outFile.empty())
	{
		std::string prefix = stripFilename(opt::readsFile);
		if(!opt::targetFile.empty())
		{
			prefix.append(1,'.');
			prefix.append(stripFilename(opt::targetFile));
		}
		opt::outFile = prefix + ASQG_EXT + GZIP_EXT;
	}
}
