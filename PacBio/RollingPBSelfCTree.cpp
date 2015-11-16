///----------------------------------------------
// Copyright 2015 National Chung Cheng University
// Written by Yao-Ting Huang & Ping-Yeh Chen
// Released under the GPL
//-----------------------------------------------

#include "RollingPBSelfCTree.h"
#include "BWTAlgorithms.h"
#include "cyclichash.h"

using namespace std;

//
// Class: RollingPBSelfCTree
RollingPBSelfCTree::RollingPBSelfCTree(
					const BWT* pBWT,
                    const BWT* pRBWT,
					size_t smallKmerSize,
					size_t min_SA_threshold,
					int maxLeavesAllowed) :
                    m_pBWT(pBWT), m_pRBWT(pRBWT), m_smallKmerSize(smallKmerSize), 
					m_min_SA_threshold(min_SA_threshold), m_maxLeavesAllowed(maxLeavesAllowed),
					m_expectedLength(0), m_currentLength(0), m_pRootNode(NULL)
{	
	// L=smallKmerLength*2 bit hash values for all possible kmers of length smallKmerLength:
	m_rollinghasher = new CyclicHash<size_t, unsigned char>(m_smallKmerSize, 63);
	assert(m_smallKmerSize<=32);
	
	// google dense hashmap requiring set empty key
	// size_t n = -1;
	// With the default value of -1, n becomes the largest value representable for size_t 
	// this special value won't be reported by rollinghasher iff L < 64
	// kmerHash.set_empty_key(n);
	
	// Estimate of largest possible hashtable size for reduce collison
	// kmerHash.resize(6000);

}

//
RollingPBSelfCTree::~RollingPBSelfCTree()
{
    // Recursively destroy the tree
	if(m_pRootNode!=NULL)
		delete m_pRootNode;
	
	kmerHashiter hashiter = kmerHash.begin();
	for(; hashiter != kmerHash.end(); ++hashiter)
	{
		delete hashiter->second;
		hashiter->second = NULL;
	}

	if(m_rollinghasher!=NULL)
		delete m_rollinghasher;
}

// create root node using src
void RollingPBSelfCTree::initializeSearchTree(std::string src, size_t hashKmerSize)
{
	m_leaves.clear();

	// Create the root node containing the seed string
    RollingNode *pRootNode = new RollingNode(&src, NULL);	
	// store initial str of root
    pRootNode->computeInitial(src);
	
	// extract suffix of length minOverlap
    std::string beginningkmer=src.substr(src.length() - hashKmerSize);
	// std::cout << beginningkmer << "\n";

	// Compute initial SA intervals in both forward and reverse BWT
    pRootNode->currFwdKmer = new CircularKmer(beginningkmer);
    pRootNode->currFwdHashValue = m_rollinghasher->hashString(beginningkmer);
	
	std::string rvckmer = reverseComplement(beginningkmer);
	pRootNode->currRvcKmer = new CircularKmer(rvckmer);
	pRootNode->currRvcHashValue = m_rollinghasher->hashString(rvckmer);
	
	m_leaves.push_back(pRootNode);
	
	m_pRootNode = pRootNode;	
	m_seedLength = src.length();
	m_currentLength = src.length();
	
	// std::cout << beginningkmer << "\t" << pRootNode->fwdInterval.size() << "\t" << pRootNode->rvcInterval.size() << "\n";
	
}

void RollingPBSelfCTree::initializeTerminalHashValues(std::string dest, size_t hashKmerSize)
{
	//ending kmer is a prefix of dest
    //initialize the ending SA intervals with kmer length=m_minOverlap	
	std::string endingkmer = dest.substr(0, hashKmerSize);
	m_fwdTerminatedValue = m_rollinghasher->hashString(endingkmer);
	endingkmer = reverseComplement(endingkmer);
	m_rvcTerminatedValue = m_rollinghasher->hashString(endingkmer); 
}

// find feasible extension using kmers collected by addHashFromSingleSeed or addHashFromPairedSeed
int RollingPBSelfCTree::mergeTwoSeedsUsingHash(const std::string &src, const std::string &dest, std::string &mergedseq, 
											size_t minLength, size_t maxLength, size_t expectedLength)
{	
	initializeSearchTree(src, m_smallKmerSize);
	initializeTerminalHashValues(dest, m_smallKmerSize);
	m_expectedLength = expectedLength;
	
	SAIntervalNodeResultVector results;
	while(!m_leaves.empty() && m_leaves.size() <= m_maxLeavesAllowed &&  (size_t)m_currentLength <= maxLength)
	{
		// std::cout << m_leaves.size() << "\t" << currentLength << "\t" << expectedLength << "\t" << maxLength << "\n";		
		RollingNodePtrList newLeaves;
	
		// attempt to extend one base for each leave
		attempToExtendUsingRollingHash(newLeaves, m_smallKmerSize);

		// extension succeed
		if(!newLeaves.empty())
			m_currentLength++;

		m_leaves.clear();
		m_leaves = newLeaves;
		
		// see if terminating string is reached
		if( (size_t) m_currentLength >= minLength)
			isTerminated(results);

		// if(src.length() == 2216 || src.length() == 2190 || src.length() == 1112){
			// std::cout << src << "\n";
			// getchar();
			// printLeaves(m_smallKmerSize);
		// }
	}
	
	string ans;
	// find the path with minimum difference or maximum kmer coverage
	if(results.size()>0)
	{
		size_t maxKmerCoverage = 0;
		int lengthDiff = 100000;

		for (size_t i = 0 ; i < results.size() ;i++)
		{
			std::string tmpseq;
			// bug fix: dest may be shorter than m_smallKmerSize
			if(dest.length() > m_smallKmerSize)
				tmpseq=results[i].thread + dest.substr(m_smallKmerSize);
			else
				tmpseq=results[i].thread;
				
			//if(debug == true) cout << ">" << i << endl << tmpseq << endl;
			
			size_t cov = results[i].SAICoverage;
			int ansDiff = tmpseq.length() - m_expectedLength;
			
			if(abs(ansDiff) <= lengthDiff && cov > maxKmerCoverage)
			{
				maxKmerCoverage=cov;
				ans=tmpseq;
				lengthDiff = abs(ansDiff);
			}
		}
	}
			
	if( ans.length()>0 )
	{
		mergedseq = ans;
		return 1;
    }
	
    // Did not reach the terminal kmer
    if(m_leaves.empty() && m_currentLength >= (int)(expectedLength-m_seedLength)/2+m_seedLength)
        return -1;	// high error
    else if( (size_t) m_currentLength > maxLength)
        return -2;	// exceed search depth
    else if(m_leaves.size() > m_maxLeavesAllowed)
        return -3;	// too much repeats
	else if(m_leaves.empty() && m_currentLength < (int)(expectedLength-m_seedLength)/2+m_seedLength)
		return -4;
	else
		return -5;
}


// LF-mapping of each SA index in the interval independently using loop instead of BFS tree expansion
// Contaminated reads often are simple repeats C* or T* with large freq
// Give up this read if the freq is way too large
bool RollingPBSelfCTree::addHashBySingleSeed(std::string& seedStr, size_t largeKmerSize, size_t smallKmerSize, size_t maxLength, int expectedLength)
{
	// PacBio errors create repeat-like seeds with large interval
	// limit the upper bound of interval for speedup
	const int64_t maxIntervalSize = 25;
	
	// LF-mapping of each fwd index
	std::string initKmer = seedStr.substr(seedStr.length() - largeKmerSize);
    BWTInterval fwdInterval=BWTAlgorithms::findInterval(m_pRBWT, reverse(initKmer));
	BWTInterval rvcInterval=BWTAlgorithms::findInterval(m_pBWT, reverseComplement(initKmer));

	std::cout << largeKmerSize << "\t" << fwdInterval.size() << "\t" << rvcInterval.size() <<"\n";
	
	for(int64_t fwdRootIndex = fwdInterval.lower; 
		fwdInterval.isValid() && fwdRootIndex <= fwdInterval.upper &&  fwdRootIndex - fwdInterval.lower < maxIntervalSize; 
		fwdRootIndex++)
	{
		// extract small hash Kmer
		// std::string currentFwdKmer = seedStr.substr(seedStr.length() - smallKmerSize);
		std::string initKmer = seedStr.substr(seedStr.length() - smallKmerSize);
		CircularKmer currentFwdKmer(initKmer);
		
		// m_rollinghasher can be only one instance in each tree due to mapping by random number generator
		m_rollinghasher->hashString(initKmer);
		
		int64_t fwdIndex = fwdRootIndex;
		for(int64_t currentLength = (int64_t)seedStr.length(); currentLength <= (int64_t)maxLength; currentLength++)
		{
			char b = m_pRBWT->getChar(fwdIndex);
			if(b == '$') break;
			
			// remove first char, append b to last char
			m_rollinghasher->update(currentFwdKmer.getFirstChar(), b);
			currentFwdKmer.shiftRight(b);
			
			// kmerHashiter hashiter = kmerHash.find(currentFwdKmer);
			kmerHashiter hashiter = kmerHash.find(m_rollinghasher->hashvalue);
			
			if(hashiter == kmerHash.end())
			{
				KmerFeatures *newEntry;
				if(expectedLength<0)
					newEntry = new KmerFeatures(currentLength-seedStr.length(), maxLength);
				else
					newEntry = new KmerFeatures(expectedLength - currentLength + smallKmerSize, maxLength);
					
				// kmerHash.insert(std::make_pair<std::string, KmerFeatures*>(currentFwdKmer, newEntry));
				kmerHash.insert(std::make_pair<size_t, KmerFeatures*>(m_rollinghasher->hashvalue, newEntry));
			}
			else
			{
				if(expectedLength<0)
					hashiter->second->add(currentLength-seedStr.length());
				else
					hashiter->second->add(expectedLength - currentLength + smallKmerSize);
			}
			// LF mapping
            fwdIndex = m_pRBWT->getPC(b) + m_pRBWT->getOcc(b, fwdIndex - 1);			
		}
	}

	std::cout << kmerHash.size() << "\n";

	// LF-mapping of each rvc index	
	for(int64_t rvcRootIndex=rvcInterval.lower; 
		rvcRootIndex <= rvcInterval.upper && rvcInterval.isValid() && rvcRootIndex - rvcInterval.lower < maxIntervalSize; 
		rvcRootIndex++)
	{
		std::string initKmer = reverseComplement( seedStr.substr(seedStr.length() - smallKmerSize) );
		CircularKmer currentRvcKmer(initKmer);
		
		m_rollinghasher->hashString(initKmer);

		int64_t rvcIndex = rvcRootIndex;
		
		for(int64_t currentLength = (int64_t)seedStr.length(); currentLength <= (int64_t)maxLength; currentLength++)
		{
			char b = m_pBWT->getChar(rvcIndex);
			if(b == '$') break;
			
			// remove last char, prepend b to first char
			m_rollinghasher->reverse_update(b, currentRvcKmer.getLastChar());
			currentRvcKmer.shiftLeft(b);

			kmerHashiter hashiter = kmerHash.find(m_rollinghasher->hashvalue);
			if(hashiter == kmerHash.end())
			{				
				KmerFeatures *newEntry;
				if(expectedLength<0)
					newEntry = new KmerFeatures(currentLength - seedStr.length(), maxLength);
				else
					newEntry = new KmerFeatures(expectedLength - currentLength + smallKmerSize, maxLength);
					
				// kmerHash.insert(std::make_pair<std::string, KmerFeatures*>(currentRvcKmer, newEntry));
				kmerHash.insert(std::make_pair<size_t, KmerFeatures*>(m_rollinghasher->hashvalue, newEntry));
			}
			else
			{
				if(expectedLength<0)
					hashiter->second->add(currentLength-seedStr.length());
				else
					hashiter->second->add(expectedLength - currentLength + smallKmerSize);
			}
			// LF mapping
            rvcIndex = m_pBWT->getPC(b) + m_pBWT->getOcc(b, rvcIndex - 1);
		}
	}
	
	std::cout << kmerHash.size() << "\n";
	return true;
}

// void RollingPBSelfCTree::printLeaves(size_t hashKmerSize)
// {
	// std::cout << m_leaves.size() << ":" << m_currentLength << "\n";
	//Bamboo.PB.45195.gap.fa
	//cout << "CGCCAGCTGAGCTGGCGGTGTGAAATCAGGCAGTCAGGCGGCTCGCGTCTTGCGCGATAACCAGTTCTTCGTTGGTTGGGATAACCACCGCAGGACGGGTACCTTCTTTGTTGATGAAACCAGATTTGCCGAAACGTGCAGCCAGGTTGCGTTCATGATCAACTTCAAAGCCCAGCACGCCCAGTTTGCCCAGAGACAGTTCACGA\n";

	//PB248_9895	
	// cout << "GCGAGACGGTATTACCCGGCCCCTGGTCGCGCGGCAGGTTATGAATATTCTGTTCATGCAGGGAAAAACTCCCCGCCAGTGTAGCGATTTCACGCTCAGCAACATGGCGCGGCACACCAGCTAATAGAACTTCTCCACGCATCTGCACAATGTTCCCGCGCTCGCCAAGTTGCAAGGTGTTAAACGATGCCACGGGCGAGACTTCCGTTGCCAC\n";
	// for(STNodePtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
	// {
		// std::string STNodeStr = (*iter)->getFullString();
		// std::string fwdrepeatunit = (*iter)->getSuffix(hashKmerSize);
		// kmerHashiter hashiter = kmerHash.find(fwdrepeatunit);
		
		// if(hashiter!=kmerHash.end() 
			// // && kmerHash[fwdrepeatunit]->getTotalFreq() > 0 
			// // && kmerHash[fwdrepeatunit]->getSumOfFreq(m_currentLength-m_seedLength)>0 
		  // )
		// {
			// std::cout << STNodeStr.substr(m_seedLength-hashKmerSize) << " " << fwdrepeatunit 
			// << ":" << kmerHash[fwdrepeatunit]->getTotalFreq() 
			// << ":" << kmerHash[fwdrepeatunit]->getSumOfFreq(m_currentLength-m_seedLength) 
			// // << ":" << kmerHash[fwdrepeatunit]->getTotalSum()/kmerHash[fwdrepeatunit]->getTotalFreq() + src.length()
			// // << ":" << kmerHash[fwdrepeatunit]->getSumOfPos(m_currentLength-m_seedLength)/ kmerHash[fwdrepeatunit]->getSumOfFreq(m_currentLength-m_seedLength)+m_seedLength
			// ;
		// }
		
		// std::string rvcrepeatunit = reverseComplement(fwdrepeatunit);
		// kmerHashiter hashiter2 = kmerHash.find(rvcrepeatunit);
		
		
		// if(hashiter2!=kmerHash.end() 
			// // && kmerHash[rvcrepeatunit]->getTotalFreq() > 0 
			// // && kmerHash[rvcrepeatunit]->getSumOfFreq(m_currentLength-m_seedLength)>0 
			// )
		// {
			// std::cout << "-"  << " " << rvcrepeatunit
			// << ":" << kmerHash[rvcrepeatunit]->getTotalFreq() 
			// << ":" << kmerHash[rvcrepeatunit]->getSumOfFreq(m_currentLength-m_seedLength) 
			// // << ":" << kmerHash[rvcrepeatunit]->getTotalSum()/kmerHash[rvcrepeatunit]->getTotalFreq() + m_seedLength
			// // << ":" << kmerHash[rvcrepeatunit]->getSumOfPos(m_currentLength-m_seedLength)/ kmerHash[rvcrepeatunit]->getSumOfFreq(m_currentLength-m_seedLength)+m_seedLength
			// ;
		// }
		
		// if(hashiter!=kmerHash.end() || hashiter2!=kmerHash.end())
		// {
			// std::cout << "--" << (double)(*iter)->getKmerCount()/m_currentLength;
			// std::cout << "\n"; 
		// }
		
		// // assert(hashiter!=kmerHash.end() || hashiter2!=kmerHash.end());
	// }
// }

// Print the string represented by every node
void RollingPBSelfCTree::printAll()
{
    std::cout << "Print all: \n";
    m_pRootNode->printAllStrings("");
}

// Extend each leaf node
void RollingPBSelfCTree::attempToExtendUsingRollingHash(RollingNodePtrList &newLeaves, size_t hashKmerSize)
{
    for(RollingNodePtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
    {
		for(int i = 1; i < BWT_ALPHABET::size; ++i) //i=A,C,G,T
		{
			size_t kmerFreq = 0;

			// compute new forward hash value using rollinghasher
			char b = BWT_ALPHABET::getChar(i);
			size_t newFwdHashValue;
			
			if((*iter)->currFwdHashValue == (size_t) -1)
			{
				std::string newKmer = (*iter)->getSuffix(hashKmerSize-1);
				newKmer = newKmer + b;
				newFwdHashValue = m_rollinghasher->hashString(newKmer);				
			}
			else 
				newFwdHashValue = m_rollinghasher->update((*iter)->currFwdHashValue, 
													 (*iter)->currFwdKmer->getFirstChar(),
													 b);
													 
			// query frequency using dense hashmap
			kmerHashiter fwditer = kmerHash.find(newFwdHashValue);
			
			// bubble removal by greedily keeping the path with larger avg frequency
			double currAvgFreq = (double)(*iter)->getKmerCount()/m_currentLength;					
			if( fwditer!=kmerHash.end() && currAvgFreq < fwditer->second->getMaxAvgFreq())
				continue;
			else if(fwditer!=kmerHash.end() && currAvgFreq > fwditer->second->getMaxAvgFreq())
				fwditer->second->setMaxAvgFreq( currAvgFreq );
				
			if(fwditer!=kmerHash.end())
				kmerFreq += fwditer->second->getSumOfFreq(m_currentLength - m_seedLength);
			else
				newFwdHashValue = (size_t)-1;

			// compute new rvc hash value using rollinghasher
			char rvcb = BWT_ALPHABET::getChar(5-i);
			size_t newRvcHashValue;

			if((*iter)->currRvcHashValue == (size_t) -1)
			{
				std::string newKmer = (*iter)->getSuffix(hashKmerSize-1);
				newKmer = reverseComplement(newKmer);
				newKmer = b + newKmer;
				newRvcHashValue = m_rollinghasher->hashString(newKmer);				
			}
			else 
				newRvcHashValue= m_rollinghasher->reverse_update((*iter)->currRvcHashValue, 
															rvcb,
														   (*iter)->currRvcKmer->getLastChar());
																	
			kmerHashiter rvciter = kmerHash.find(newRvcHashValue);
			if(rvciter!=kmerHash.end())
				kmerFreq += rvciter->second->getSumOfFreq(m_currentLength - m_seedLength);
			else
				newRvcHashValue = -1;
				
			if(kmerFreq > m_min_SA_threshold)
			{
				std::string bstr;
				bstr.push_back(b);
				RollingNode* pChildNode = (*iter)->createChild(bstr);
				pChildNode->currFwdHashValue = newFwdHashValue;
				pChildNode->currRvcHashValue = newRvcHashValue;				
				pChildNode->currFwdKmer->shiftRight(b);
				pChildNode->currRvcKmer->shiftLeft(rvcb);
								
				// inherit kmer freq from parents
				pChildNode->addKmerCount((*iter)->getKmerCount());
				pChildNode->addKmerCount(kmerFreq);
				newLeaves.push_back(pChildNode);				
			}
		}
	}
}


// Check for leaves whose extension has terminated. If the leaf has
// terminated, the walked string and coverage is pushed to the result vector

// Check for leaves whose extension has terminated. If the leaf has
// terminated, the walked string and coverage is pushed to the result vector
bool RollingPBSelfCTree::isTerminated(SAIntervalNodeResultVector& results)
{
	bool found = false;

    for(RollingNodePtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
    {
        size_t currFwdHashvalue=(*iter)->currFwdHashValue;
        size_t currRvcHashvalue=(*iter)->currRvcHashValue;


		//The current SA interval stands for a string >= terminating kmer
		//If terminating kmer is a substr, the current SA interval is a sub-interval of the terminating interval
        bool isFwdTerminated= currFwdHashvalue == m_fwdTerminatedValue;
        bool isRvcTerminated= currRvcHashvalue == m_rvcTerminatedValue;

        if(isFwdTerminated || isRvcTerminated)
        {
            std::string STNodeStr = (*iter)->getFullString();
			// std::cout << "Terminated: " << STNodeStr << "\n";
			
            SAIntervalNodeResult STresult;
            STresult.thread=STNodeStr;
			STresult.SAICoverage=(*iter)->getKmerCount();

            //compute the merged pos right next to the kmer on 2nd read.
            //STresult.index = m_minOverlap ;
            results.push_back(STresult);
            found =  true;
        }
    }

    return found;
}
