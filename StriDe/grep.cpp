//
// grep - grep reads by substring
//
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include "SGACommon.h"
#include "Util.h"
#include "grep.h"
#include "SuffixArray.h"
#include "BWT.h"
#include "Timer.h"
#include "BWTAlgorithms.h"
#include <iomanip>
//
// Getopt
//
#define SUBPROGRAM "grep"

static const char *GREP_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"\n";

static const char *GREP_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... READSFILE\n"
"Get sequences kmer frequency\n"
"  -v, --verbose                        display verbose output\n"
"      --help                           display this help and exit\n";

namespace opt
{
    static unsigned int verbose;
	static std::string readsFile;
}

static const char* shortopts = "v";

enum { OPT_HELP = 1, OPT_VERSION, };

static const struct option longopts[] = {
    { "verbose",     no_argument,       NULL, 'v' },
    { "help",        no_argument,       NULL, OPT_HELP },
    { "version",     no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};



const char * setYellowStart = "\033[33m";
const char * setYellowEnd = "\033[0m";




int grepMain(int argc, char** argv)
{
    Timer t("sga grep");
	parseGREPOptions(argc, argv);
	
	std::string  prefix = stripFilename(opt::readsFile);
	Timer* pLTimer = new Timer("Load Time");
	BWT* pBWT = new BWT(prefix + BWT_EXT, BWT::DEFAULT_SAMPLE_RATE_SMALL);
	SampledSuffixArray* pSSA = new SampledSuffixArray(prefix + SAI_EXT, SSA_FT_SAI);
	ReadTable* pRT = new ReadTable(opt::readsFile);
	delete pLTimer;
	
	std::vector<SeqItem> result;

	std::string queryString ;
	while (std::cin >> queryString) 
	{
		std::cout << "--" << std::endl;
		BWTInterval interval = BWTAlgorithms::findInterval(pBWT, queryString);
		if(interval.isValid()) 
		{
			for(int64_t idx = interval.lower; idx <= interval.upper; ++idx)
			{

				SeqItem  seqItem = pRT->getRead(pSSA->calcSA(idx,pBWT).getID());
				
				result.push_back(seqItem);
				
				std::cout << seqItem.id  <<"\n" ;
				//std::cout << seqItem.seq.toString() << std::endl;
				
				std::string read = seqItem.seq.toString();
				size_t found = read.find(queryString);
				
				std::cout << read.substr(0,found);
				std::cout << setYellowStart <<read.substr(found,queryString.length()) <<setYellowEnd;
				std::cout << read.substr(found+queryString.length(),read.length()-found-queryString.length()) <<"\n";
				
				
   
			}
		}
		std::cout << "--" << std::endl;
	}
	
	std::vector<SeqItem>::iterator  it = result.begin();
	
	for (;it!=result.end();)
	{
		std::vector<SeqItem>::iterator found= result.begin()  ;
		for (; found!=it ; found++ )
		{
			if (found->id == it->id)  break;
		}
		
		if (found != it )  it = result.erase(it);
		else it ++ ;
	}
	
	
	
	for (it = result.begin();it!=result.end();it++)
	{
		
		std::cout << ">" <<it->id  << "\n" << it->seq.toString() << std::endl; ;
		
	
	}

	
	/*
	pRT->indexReadsByID();
	SeqItem  si = pRT->getRead("HWI-ST180_0168:6:24:8510:155485#ACCAACT/1");
	std::cout << setYellowStart << si.id  <<"\n" << si.seq.toString() <<setYellowEnd<< std::endl;
	*/
	
	
	delete pBWT;
	delete pSSA;
	delete pRT;
	
    return 0;
}

// 
// Handle command line arguments
//
void parseGREPOptions(int argc, char** argv)
{
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) 
        {
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
            case OPT_HELP:
                std::cout << GREP_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << GREP_VERSION_MESSAGE;
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
        std::cout << "\n" << GREP_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

	opt::readsFile = argv[optind++];
}
