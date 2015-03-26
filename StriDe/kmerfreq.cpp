//
// kmerfreq - get sequences kmer frequency
//
#include <iostream>
#include <fstream>
#include "SGACommon.h"
#include "Util.h"
#include "kmerfreq.h"
#include "SuffixArray.h"
#include "BWT.h"
#include "Timer.h"
#include "BWTAlgorithms.h"
#include <iomanip>
//
// Getopt
//
#define SUBPROGRAM "kmerfreq"

static const char *KMERFREQ_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"\n";

static const char *KMERFREQ_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... BWTFILE\n"
"Get sequences kmer frequency\n"
"  -k, --kmer-length                    kmer size\n"
"  -v, --verbose                        display verbose output\n"
"      --help                           display this help and exit\n";

namespace opt
{
    static unsigned int verbose;
    static std::string bwtFile;
	static size_t kmerLength = 31;
    static int sampleRate = 256;
}

static const char* shortopts = "k:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "verbose",     no_argument,       NULL, 'v' },
	{ "kmer-length",           no_argument,       NULL, 'k' },
    { "help",        no_argument,       NULL, OPT_HELP },
    { "version",     no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

int kmerfreqMain(int argc, char** argv)
{
    Timer t("sga kmerfreq");
	parseKMERFREQOptions(argc, argv);
	BWT* pBWT = new BWT(opt::bwtFile, opt::sampleRate);
    //pBWT->printInfo();
	
	std::string queryString ;
	
	//std::cout <<  opt::kmerLength << std::endl;
	
	
	while (std::cin >> queryString) 
	{
		//std::cout << queryString << std::endl;
		size_t seqLength = queryString.length();
		size_t numKmer = seqLength-opt::kmerLength+1 ;
		std::vector<std::string> kmers (numKmer);
		std::vector<size_t> kmerFreqs_same (numKmer);
		std::vector<size_t> kmerFreqs_revc (numKmer);
		
		#pragma omp parallel for
		for (size_t i = 0 ; i < numKmer  ; i++)
		{
			kmers[i] = queryString.substr (i,opt::kmerLength);
			kmerFreqs_same[i] = BWTAlgorithms::countSequenceOccurrencesSingleStrand(kmers[i], pBWT) ;
			kmerFreqs_revc[i] = BWTAlgorithms::countSequenceOccurrencesSingleStrand(reverseComplement(kmers[i]), pBWT) ;
		}	
		
		for (size_t i = 0 ; i < numKmer  ; i++)
		{
			std::cout << i << "\t" <<  kmers[i] << " " << kmerFreqs_same[i] <<  "\t" 
			          << reverseComplement(kmers[i]) <<  " " <<kmerFreqs_revc[i] << std::endl;
		
		}
		
		
		
	}
	delete pBWT;

    return 0;
}

// 
// Handle command line arguments
//
void parseKMERFREQOptions(int argc, char** argv)
{
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) 
        {
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
			case 'k': arg >> opt::kmerLength; break;
            case OPT_HELP:
                std::cout << KMERFREQ_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << KMERFREQ_VERSION_MESSAGE;
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

    if (die) 
    {
        std::cout << "\n" << KMERFREQ_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

	opt::bwtFile = argv[optind++];
}
