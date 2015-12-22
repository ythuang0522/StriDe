//-----------------------------------------------
// Copyright 2014 National Chung Cheng University
// Written by Yao-Ting Huang
// Released under the GPL
//-----------------------------------------------
//
// all: perform all assembly tasks in one click
//
#include <string>
#include <iostream>
#include <fstream>
#include "Timer.h"
#include "Util.h"

#include "strideall.h"
#include "index.h"
#include "overlap.h"
#include "assemble.h"
#include "preprocess.h"
#include "correct.h"
#include "filter.h"
#include "fm-merge.h"
#include "FMIndexWalk.h"


//
// Getopt
//
#define SUBPROGRAM "all"

static const char *STRIDE_VERSION_MESSAGE =
"StriDe Assembler  " PACKAGE_VERSION "\n"
"Written by Yao-Ting Huang.\n"
"\n"
"Copyright 2014 National Chung Cheng University\n";

static const char *STRIDE_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... READFILE (format controlled by -p) ...\n"
"Perform error correction, long-read generation, overlap, and assembly in one command.\n"
"\n"
"Mandatory arguments:\n"
"      -r, --read-length=LEN            median read length (if there are multiple libraries, set to the max read length)\n"
"      -i, --insert-size=LEN            median insert size (if there are multiple libraries, set to the max insert size)\n"
"\nOptional arguments:\n"
"      -t, --thread=N                   number of threads (default: 16)\n"
"      -p, --pe-mode=INT                1 - paired reads are separated with the first read in READS1 and the second\n"
"                                           read in READS2 of the same library (default) \n"
"                                       2 - paired reads are interleaved within a single file of the same library.\n"
"      -k, --kmer-size=N                length of kmer (default: 31)\n"
"      -c, --kmer-threshold=N           kmer frequency cutoff (default: 3)\n"
"      -m, --min-overlap=LEN            minimum reliable overlap length (default: read length * 0.8)\n"
"      --help                           display this help and exit\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
    static unsigned int verbose;

	static size_t readLength = 0 ;
	static size_t insertSize=0;
	static size_t numthread=16;
    static unsigned int peMode = 1;
	static size_t kmerLength = 31;
	static size_t kmerThreshold = 3;			//3 is good for 100x coverage 
	static size_t minOverlap = 0;				//0.8 is good for 100x coverage 

}


static const char* shortopts = "k:c:r:m:i:t:v";

enum { OPT_HELP = 1, OPT_VERSION};

static const struct option longopts[] = {
	{ "verbose",               no_argument,       NULL, 'v' },
	{ "read-length",           required_argument, NULL, 'r' },
	{ "insert-size",           required_argument, NULL, 'i' },
	{ "thread",           required_argument, NULL, 't' },
	{ "pe-mode",                required_argument, NULL, 'p' },
	{ "kmer-length",           required_argument, NULL, 'k' },
	{ "kmer-threshold",        required_argument, NULL, 'c' },
	{ "min-overlap",      required_argument, NULL, 'm' },
	{ "help",                  no_argument,       NULL, OPT_HELP },
	{ "version",               no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};


int StrideMain(int argc, char** argv)
{
    Timer* pTimer = new Timer("Stride all-in-one");
    parseStrideOptions(argc, argv);

	//$Stride preprocess --discard-quality -p 1 insert_180_1.fastq insert_180_2.fastq -o reads.fa
	std::vector <std::string> vec;
	vec.push_back("preprocess");
	vec.push_back("--discard-quality");
	vec.push_back("-p");
	vec.push_back("1");
	vec.push_back("-o");
	vec.push_back("reads.fa");

	while(optind < argc)
    {
		int numFiles = argc - optind;
		if(opt::peMode==1 && numFiles % 2 == 1)
		{
				std::cerr << "Error: An even number of files must be given for pe-mode 1\n";
				exit(EXIT_FAILURE);
		}

		if(opt::peMode == 1)
        {
			std::string filename1 = argv[optind++];
			std::string filename2 = argv[optind++];
			vec.push_back(filename1);
			vec.push_back(filename2);
		}
		else
		{
			// Read PE reads from an interleaved single file
            std::string filename = argv[optind++];
			vec.push_back(filename);
		}
	}

	char **arr = new char*[vec.size()];
	for(size_t i = 0; i < vec.size(); i++){
		arr[i] = new char[vec[i].size() + 1];
		strcpy(arr[i], vec[i].c_str());
	}
    preprocessMain(vec.size(), arr);
	free(arr);

	std::cout << "\n\n\t [ Stage I: Error correction ] \n\n";

	//$Stride index -a ropebwt -t $thread reads.fa
	std::vector <std::string> vec1;
	vec1.push_back("index");
	vec1.push_back("-t");
	std::stringstream ssnumthread;
	ssnumthread << opt::numthread;
	vec1.push_back(ssnumthread.str());
	vec1.push_back("reads.fa");

	char ** arr1 = new char*[vec1.size()];
	for(size_t i = 0; i < vec1.size(); i++){
		arr1[i] = new char[vec1[i].size() + 1];
		strcpy(arr1[i], vec1[i].c_str());
	}
    indexMain(vec1.size(), arr1);
	free(arr1);

	//$Stride correct -a overlap -r 11 -t $thread -k $kmersize -x $kmerthreshold reads.fa -o READ.ECOLr.fasta
	std::vector <std::string> vec2;
	vec2.push_back("correct");
	vec2.push_back("-a");
	vec2.push_back("overlap");
	vec2.push_back("-r");
	vec2.push_back("1");
	vec2.push_back("-t");
	vec2.push_back(ssnumthread.str());
	vec2.push_back("-k");
	std::stringstream sskmerLength;
	sskmerLength << opt::kmerLength;
	vec2.push_back(sskmerLength.str());
	vec2.push_back("-x");
	std::stringstream sskmerThreshold;
	sskmerThreshold << opt::kmerThreshold;
	vec2.push_back(sskmerThreshold.str());
	vec2.push_back("-o");
	vec2.push_back("READ.ECOLr.fasta");
	vec2.push_back("reads.fa");

	char ** arr2 = new char*[vec2.size()];
	for(size_t i = 0; i < vec2.size(); i++){
		arr2[i] = new char[vec2[i].size() + 1];
		strcpy(arr2[i], vec2[i].c_str());
	}
    correctMain(vec2.size(), arr2);
	free(arr2);

	std::cout << "\n\n\t [ Stage II: merge paired-end reads into long reads and kmerize error-prone reads ] \n\n";

	//$Stride index -a ropebwt -t $thread READ.ECOLr.fasta
	std::vector <std::string> vec3;
	vec3.push_back("index");
	vec3.push_back("-t");
	vec3.push_back(ssnumthread.str());
	vec3.push_back("READ.ECOLr.fasta");
	char ** arr3 = new char*[vec3.size()];
	for(size_t i = 0; i < vec3.size(); i++){
		arr3[i] = new char[vec3[i].size() + 1];
		strcpy(arr3[i], vec3[i].c_str());
	}
    indexMain(vec3.size(), arr3);
	free(arr3);

	//Stride fmwalk -m 80 -t $thread -L $MaxBFSLeaves -I $BFSSearchDepth -x $kmerthreshold -k $kmersize READ.ECOLr.fasta
	std::vector <std::string> vec4;
	vec4.push_back("fmwalk");
	vec4.push_back("-m");
	std::stringstream ssminOverlapLength;
	ssminOverlapLength << opt::minOverlap;
	vec4.push_back(ssminOverlapLength.str());
	vec4.push_back("-t");
	vec4.push_back(ssnumthread.str());
	vec4.push_back("-L");
	vec4.push_back("64");
	vec4.push_back("-I");
	std::stringstream ssMaxInsertSize;
	ssMaxInsertSize << opt::insertSize*2;
	vec4.push_back(ssMaxInsertSize.str());	
	vec4.push_back("-k");
	vec4.push_back(sskmerLength.str());
	vec4.push_back("-x");
	vec4.push_back(sskmerThreshold.str());
	vec4.push_back("READ.ECOLr.fasta");

	char ** arr4 = new char*[vec4.size()];
	for(size_t i = 0; i < vec4.size(); i++){
		arr4[i] = new char[vec4[i].size() + 1];
		strcpy(arr4[i], vec4[i].c_str());
	}	
    FMindexWalkMain(vec4.size(), arr4);
	free(arr4);
	
	//cat READ.ECOLr.merge.fa READ.ECOLr.kmerized.fa >merged.fa
	std::ifstream if_merge("READ.ECOLr.merge.fa", std::ios_base::binary);
	std::ifstream if_kmerize("READ.ECOLr.kmerized.fa", std::ios_base::binary);
	std::ofstream of_merged("merged.fa", std::ios_base::binary);
	if( if_merge.peek() != std::ifstream::traits_type::eof() )
		of_merged << if_merge.rdbuf();
	if( if_kmerize.peek() != std::ifstream::traits_type::eof() )
		of_merged << if_kmerize.rdbuf();
	if_merge.close();
	if_kmerize.close();
	of_merged.close();
	unlink("READ.ECOLr.kmerized.fa");
	//unlink("READ.ECOLr.merge.fa");
	
	std::cout << "\n\n\t [ Stage III:  Filter redundant reads] \n\n";
	//$Stride index -t $thread merged.fa
	std::vector <std::string> vec5;
	vec5.push_back("index");
	vec5.push_back("-t");
	vec5.push_back(ssnumthread.str());
	vec5.push_back("merged.fa");

	char ** arr5 = new char*[vec5.size()];
	for(size_t i = 0; i < vec5.size(); i++){
		arr5[i] = new char[vec5[i].size() + 1];
		strcpy(arr5[i], vec5[i].c_str());
	}
    indexMain(vec5.size(), arr5);
	free(arr5);
	
	//Stride filter -t $thread --rebuild-BWT --no-kmer-check merged.fa
	std::vector <std::string> vec6;
	vec6.push_back("filter");
	vec6.push_back("-t");
	vec6.push_back(ssnumthread.str());
	//vec6.push_back("--rebuild-BWT");
	//vec6.push_back("--no-kmer-check");
	vec6.push_back("merged.fa");

	char** arr6 = new char*[vec6.size()];
	for(size_t i = 0; i < vec6.size(); i++){
		arr6[i] = new char[vec6[i].size() + 1];
		strcpy(arr6[i], vec6[i].c_str());
	}
    filterMain(vec6.size(), arr6);
	free(arr6);

	std::cout << "\n\n\t [ Stage IV:  Compute overlap ] \n\n";
	//$Stride overlap -m $kmeroverlap -t $thread merged.filter.pass.fa
	std::stringstream sskmeroverlap;
	sskmeroverlap << opt::kmerLength-1;
	std::vector <std::string> vec7;
	vec7.push_back("overlap");
	vec7.push_back("-t");
	vec7.push_back(ssnumthread.str());
	vec7.push_back("-m");
	vec7.push_back(sskmeroverlap.str());
	vec7.push_back("merged.filter.pass.fa");

	char ** arr7 = new char*[vec7.size()];
	for(size_t i = 0; i < vec7.size(); i++){
		arr7[i] = new char[vec7[i].size() + 1];
		strcpy(arr7[i], vec7[i].c_str());
	}
    overlapMain(vec7.size(), arr7);
	free(arr7);
	unlink("merged.discard.fa");

	std::cout << "\n\n\t [ Stage V:  String Graph Assembly] \n\n";
	//$Stride assemble -k $kmersize -t $kmerthreshold -p READ.ECOLr -R $origianlReadLength -c $reliableOverlapLength merged.filter.pass.asqg.gz
	std::vector <std::string> vec8;
	vec8.push_back("assemble");
	vec8.push_back("-k");
	vec8.push_back(sskmerLength.str());
	vec8.push_back("-t");
	vec8.push_back(sskmerThreshold.str());
	vec8.push_back("-p");
	vec8.push_back("READ.ECOLr");
	std::stringstream ssinsertSize;
	ssinsertSize << opt::insertSize;
	vec8.push_back("--insert-size="+ssinsertSize.str());
	vec8.push_back("-r");
	std::stringstream ssreadLength;
	ssreadLength << opt::readLength;
	vec8.push_back(ssreadLength.str());
	vec8.push_back("-c");
	vec8.push_back(ssminOverlapLength.str());
	vec8.push_back("merged.filter.pass.asqg.gz");

	char** arr8 = new char*[vec8.size()];
	for(size_t i = 0; i < vec8.size(); i++){
		arr8[i] = new char[vec8[i].size() + 1];
		strcpy(arr8[i], vec8[i].c_str());
	}
	assembleMain(vec8.size(), arr8);
	free(arr8);
	
    delete pTimer;
    return 0;
}


// 
// Handle command line arguments
//
void parseStrideOptions(int argc, char** argv)
{
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) 
        {
		case 'k': arg >> opt::kmerLength; break;
		case 'c': arg >> opt::kmerThreshold; break;
		case 't': arg >> opt::numthread; break;
        case 'p': arg >> opt::peMode; break;
		case 'm': arg >> opt::minOverlap; break;
		case 'r': arg >> opt::readLength; break;
        case 'i': arg >> opt::insertSize; break;
		case 'v': opt::verbose++; break;
		case '?': die = true; break;
		case OPT_HELP:
			std::cout << STRIDE_USAGE_MESSAGE;
			exit(EXIT_SUCCESS);
		case OPT_VERSION:
			std::cout << STRIDE_VERSION_MESSAGE;
			exit(EXIT_SUCCESS);
        }
    }

    if (argc - optind < 2) 
    {
        std::cerr << SUBPROGRAM ": missing arguments\n";
        die = true;
    } 

    if (die) 
    {
        std::cout << "\n" << STRIDE_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

    // Parse the vertex id
    if(opt::readLength==0)
    {
        std::cerr << SUBPROGRAM ": missing maximum read length\n";
        exit(EXIT_FAILURE);
    }

	 if(opt::insertSize==0)
    {
        std::cerr << SUBPROGRAM ": missing maximum insert size\n";
        exit(EXIT_FAILURE);
    }
	
	if(opt::minOverlap==0)
    {
        opt::minOverlap=opt::readLength*0.8;
    }
}
