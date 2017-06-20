///----------------------------------------------
// Copyright 2015 National Chung Cheng University
// Written by Yao-Ting Huang & Ping-Yeh Chen
// Released under the GPL
//-----------------------------------------------

#include "HashtableSearch.h"

#include "BWTAlgorithms.h"
#include "stdaln.h"

using namespace std;


//
// Class: HashtableSearch
HashtableSearch::HashtableSearch(
					const BWT* pBWT,
                    const BWT* pRBWT
                   
) :
                    m_pBWT(pBWT), m_pRBWT(pRBWT),
					m_expectedLength(0), m_currentLength(0)
					
{
	// initialize google dense hashmap
	kmerHash.set_empty_key("");
	kmerHash.resize(6000);
}





//
HashtableSearch::~HashtableSearch()
{

   
		
	kmerHashiter hashiter = kmerHash.begin();
  
	for(; hashiter != kmerHash.end(); ++hashiter)
	{
		if(hashiter->second!=NULL)
		{
			delete hashiter->second;
			hashiter->second = NULL;
		}
	}

}


struct temphashtable
{
	
	string currentKmer;
	size_t seedStr_length;
	size_t currentLength;
	size_t smallKmerSize;
	size_t maxLength;
	size_t expectedLength;
};
// LF-mapping of each SA index in the interval independently using loop instead of BFS tree expansion
// Contaminated reads often are simple repeats C* or T* with large freq
// Give up this read if the freq is way too large
size_t HashtableSearch::addHashBySingleSeed(std::string& seedStr, size_t largeKmerSize, size_t smallKmerSize, size_t maxLength, bool skipRepeat, int expectedLength , bool istarget)
{
	// PacBio errors create repeat-like seeds with large interval
	// limit the upper bound of interval for speedup
	const int64_t maxIntervalSize = 30;
	
	// LF-mapping of each fwd index using largeKmerSize
	// assert(seedStr.length() >= largeKmerSize);
	std::string initKmer = seedStr.substr(seedStr.length() - largeKmerSize);
	BWTInterval fwdInterval=BWTAlgorithms::findInterval(m_pRBWT, reverse(initKmer));
	BWTInterval rvcInterval=BWTAlgorithms::findInterval(m_pBWT, reverseComplement(initKmer));
    
	size_t kmerFreq = 0;
	kmerFreq += fwdInterval.isValid()?fwdInterval.size():0;
	kmerFreq += fwdInterval.isValid()?rvcInterval.size():0;

	std::cout << initKmer << "\t" << smallKmerSize << "\t" << reverseComplement(initKmer) << "\t" << fwdInterval.size() << "\t" << rvcInterval.size() << "\t" << kmerFreq << " uaua\n";

	
	// skip repeat only in the 1st round
	if(skipRepeat && kmerFreq > 128) return kmerFreq;
	
	// extend each SA index and collect kmers of smallKmerSize along the extension
	for(int64_t fwdRootIndex = fwdInterval.lower; 
		fwdInterval.isValid() && fwdRootIndex <= fwdInterval.upper &&  fwdRootIndex - fwdInterval.lower < maxIntervalSize; 
		fwdRootIndex++)
	{
		// extract small hash Kmer
		std::string currentFwdKmer = seedStr.substr(seedStr.length() - smallKmerSize);

		// Bug fix: the first and last kmers in the seeds must be added.
		insertKmerToHash(currentFwdKmer, seedStr.length(), seedStr.length(), smallKmerSize, maxLength, expectedLength);		
		// std::cout << seedStr.length() << ":" << currentFwdKmer << "\n";
		// extend the fwdIndex via LF mapping
		int64_t fwdIndex = fwdRootIndex;
       std::string Fwdread = initKmer;
 
		for(int64_t currentLength = (int64_t)seedStr.length()+1; currentLength <= (int64_t)maxLength; currentLength++)
		{
			char b = m_pRBWT->getChar(fwdIndex);
			if(b == '$') break;
			
			currentFwdKmer = currentFwdKmer.substr(1) + b;
            Fwdread = Fwdread +b ;
			// std::cout << currentLength << ":" << currentFwdKmer << "\n";

			insertKmerToHash(currentFwdKmer, seedStr.length(), currentLength, smallKmerSize, maxLength, expectedLength);

			// LF mapping
            fwdIndex = m_pRBWT->getPC(b) + m_pRBWT->getOcc(b, fwdIndex - 1);			
		}
       if (istarget)
        std::cout<<  reverseComplement (Fwdread)   <<"read\n";
        else
            std::cout<<  Fwdread   <<"read\n";
	}

	
	// LF-mapping of each rvc index	
	for(int64_t rvcRootIndex=rvcInterval.lower; 
		rvcRootIndex <= rvcInterval.upper && rvcInterval.isValid() && rvcRootIndex - rvcInterval.lower < maxIntervalSize; 
		rvcRootIndex++)
	{
		std::string currentRvcKmer = reverseComplement( seedStr.substr(seedStr.length() - smallKmerSize) );
		insertKmerToHash(currentRvcKmer, seedStr.length(), seedStr.length(), smallKmerSize, maxLength, expectedLength);
		// std::cout << seedStr.length() << ":" << currentRvcKmer << "\n";
		int64_t rvcIndex = rvcRootIndex;
        std::string Rvcread = reverseComplement(initKmer);
		for(int64_t currentLength = (int64_t)seedStr.length()+1; currentLength <= (int64_t)maxLength; currentLength++)
		{
			char b = m_pBWT->getChar(rvcIndex);
			if(b == '$') break;
			
			currentRvcKmer = b + currentRvcKmer.substr(0, smallKmerSize-1);
            Rvcread = b + Rvcread;
			insertKmerToHash(currentRvcKmer, seedStr.length(), currentLength, smallKmerSize, maxLength, expectedLength);
			// std::cout << currentLength << ":" << currentRvcKmer << "\n";

			// LF mapping
            rvcIndex = m_pBWT->getPC(b) + m_pBWT->getOcc(b, rvcIndex - 1);
		}
        
        if (istarget)
            std::cout<<  Rvcread   <<" read\n";
        else 
            std::cout<<   reverseComplement(Rvcread)   <<" read\n";
	}
	
	// std::cout << kmerHash.size() << "\n";
	return kmerFreq;
}

void HashtableSearch::insertKmerToHash(std::string& insertedKmer, size_t seedStrLen, size_t currentLength, size_t smallKmerSize, size_t maxLength, int expectedLength)
{
	kmerHashiter hashiter = kmerHash.find(insertedKmer);
	
	if(hashiter == kmerHash.end())
	{
		KmerFeatures2 *newEntry;
		// source to target
		if(expectedLength<0)
			newEntry = new KmerFeatures2(currentLength - seedStrLen, maxLength);
		// target to source
		else
			newEntry = new KmerFeatures2(expectedLength - currentLength + smallKmerSize, maxLength);
			
		kmerHash.insert(std::make_pair(insertedKmer, newEntry));
	}
	else
	{
		if(expectedLength<0)
			hashiter->second->add(currentLength - seedStrLen);
		else
			hashiter->second->add(expectedLength - currentLength + smallKmerSize);
	}
}


//return hash kmer freqs
size_t HashtableSearch::hashkmerfreqs(std::string fwdkmer,size_t kemrposition)
{
    kmerHashiter iter1 = kmerHash.find(fwdkmer);	
    
    std::string rvckmer = reverseComplement(fwdkmer);
    kmerHashiter iter2 = kmerHash.find(rvckmer);
   
    size_t hashkmerfreqs;
    hashkmerfreqs = iter1==kmerHash.end()? 0 : iter1->second->getSumOfFreq(kemrposition);
    // hashkmerfreqs = iter1==kmerHash.end()? 0 : iter1->second->getTotalFreq();
    hashkmerfreqs += iter2==kmerHash.end()? 0 : iter2->second->getSumOfFreq(kemrposition);
    // hashkmerfreqs += iter2==kmerHash.end()? 0 : iter2->second->getTotalFreq();

    return hashkmerfreqs;
    
}
