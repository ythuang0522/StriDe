///----------------------------------------------
// Copyright 2015 National Chung Cheng University
// Written by Yao-Ting Huang & Ping-Yeh Chen
// Released under the GPL
//-----------------------------------------------

#include "SAIPBSelfCTree.h"
#include "BWTAlgorithms.h"
using namespace std;

//
// Class: SAIPBSelfCorrectTree
SAIPBSelfCorrectTree::SAIPBSelfCorrectTree(
					const BWT* pBWT,
                    const BWT* pRBWT,
					size_t min_SA_threshold,
					int maxLeavesAllowed) :
                    m_pBWT(pBWT), m_pRBWT(pRBWT), m_min_SA_threshold(min_SA_threshold), m_maxLeavesAllowed(maxLeavesAllowed),
					m_expectedLength(0), m_currentLength(0), m_pRootNode(NULL), debug(false)
{

	//if(targets[0].first == 169) debug = true;
}

//
SAIPBSelfCorrectTree::~SAIPBSelfCorrectTree()
{
    // Recursively destroy the tree
	if(m_pRootNode!=NULL)
		delete m_pRootNode;
		
	SparseHashMap<std::string, KmerFeatures*> :: iterator hashiter = kmerHash.begin();
	for(; hashiter != kmerHash.end(); ++hashiter)
	{
		delete hashiter->second;
		hashiter->second = NULL;
	}

}

// create root node using src
void SAIPBSelfCorrectTree::initializeSearchTree(std::string src, size_t hashKmerSize)
{
	m_leaves.clear();

	// Create the root node containing the seed string
    SAIntervalNode *pRootNode = new SAIntervalNode(&src, NULL);	
	// store initial str of root
    pRootNode->computeInitial(src);
	
	// extract suffix of length minOverlap
    std::string beginningkmer=src.substr(src.length() - hashKmerSize);
	// std::cout << beginningkmer << "\n";

	// Compute initial SA intervals in both forward and reverse BWT
    pRootNode->fwdInterval=BWTAlgorithms::findInterval(m_pRBWT, reverse(beginningkmer));
    pRootNode->rvcInterval=BWTAlgorithms::findInterval(m_pBWT, reverseComplement(beginningkmer));
	m_leaves.push_back(pRootNode);
	
	// for cleaning memory
	m_pRootNode = pRootNode;	
	m_seedLength = src.length();
	m_currentLength = src.length();
	
	// debug
	// std::cout << beginningkmer << "\t" << pRootNode->fwdInterval.size() << "\t" << pRootNode->rvcInterval.size() << "\n";
	
}

void SAIPBSelfCorrectTree::initializeTerminalIntervals(std::string dest, size_t hashKmerSize)
{
	//ending kmer is a prefix of dest
    //initialize the ending SA intervals with kmer length=m_minOverlap	
	std::string endingkmer = dest.substr(0, hashKmerSize);
	m_fwdTerminatedInterval = BWTAlgorithms::findInterval( m_pRBWT, reverse(endingkmer));
	m_rvcTerminatedInterval = BWTAlgorithms::findInterval( m_pBWT, reverseComplement(endingkmer));
}

// find feasible extension using kmers collected by addHashFromSingleSeed or addHashFromPairedSeed
int SAIPBSelfCorrectTree::mergeTwoSeedsUsingHash(const std::string &src, const std::string &dest, std::string &mergedseq, 
											size_t hashKmerSize, size_t maxLeaves, size_t minLength, size_t maxLength, size_t expectedLength)
{	
	initializeSearchTree(src, hashKmerSize);
	initializeTerminalIntervals(dest, hashKmerSize);
	m_expectedLength = expectedLength;
	
	SAIntervalNodeResultVector results;
	while(!m_leaves.empty() && m_leaves.size() <= maxLeaves &&  (size_t)m_currentLength <= maxLength)
	{
		// std::cout << m_leaves.size() << "\t" << currentLength << "\t" << expectedLength << "\t" << maxLength << "\n";
		refineSAInterval(hashKmerSize-1);
		
		STNodePtrList newLeaves;
	
		// attempt to extend one base for each leave
		// ACGT-extend the leaf nodes via updating existing SA interval
		attempToExtendUsingHash(newLeaves, hashKmerSize);

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
			// printLeaves(hashKmerSize);
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
			// bug fix: dest may be shorter than hashKmerSize
			if(dest.length() > hashKmerSize)
				tmpseq=results[i].thread + dest.substr(hashKmerSize);
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
    else if(m_leaves.size() > maxLeaves)
        return -3;	// too much repeats
	else if(m_leaves.empty() && m_currentLength < (int)(expectedLength-m_seedLength)/2+m_seedLength)
		return -4;
	else
		return -5;
}


// LF-mapping of each SA index in the interval independently using loop instead of BFS tree expansion
// Contaminated reads often are simple repeats C* or T* with large freq
// Give up this read if the freq is way too large
bool SAIPBSelfCorrectTree::addHashFromSingleSeedUsingLFMapping(std::string& seedStr, size_t largeKmerSize, size_t smallKmerSize, size_t maxLength, int expectedLength)
{
	// LF-mapping of each fwd index
	std::string initKmer = seedStr.substr(seedStr.length() - largeKmerSize);
    BWTInterval fwdInterval=BWTAlgorithms::findInterval(m_pRBWT, reverse(initKmer));
	BWTInterval rvcInterval=BWTAlgorithms::findInterval(m_pBWT, reverseComplement(initKmer));

	// std::cout << largeKmerSize << "\t" << fwdInterval.size() << "\t" << rvcInterval.size() <<"\n";
	
	// Contamination leads to large freq
	if( (fwdInterval.isValid() && fwdInterval.size() >= (int)m_maxLeavesAllowed) || 
		(rvcInterval.isValid() && rvcInterval.size() >= (int)m_maxLeavesAllowed) ) return false;

	for(int64_t fwdRootIndex = fwdInterval.lower; fwdRootIndex <= fwdInterval.upper && fwdInterval.isValid(); fwdRootIndex++)
	{
		// extract small hash Kmer
		std::string currentFwdKmer = seedStr.substr(seedStr.length() - smallKmerSize);
		int64_t fwdIndex = fwdRootIndex;
		
		for(int64_t currentLength = (int64_t)seedStr.length(); currentLength <= (int64_t)maxLength; currentLength++)
		{
			char b = m_pRBWT->getChar(fwdIndex);
			if(b == '$') break;
			
			currentFwdKmer = currentFwdKmer.substr(1) + b;

			kmerHashiter hashiter = kmerHash.find(currentFwdKmer);
			if(hashiter == kmerHash.end())
			{
				KmerFeatures *newEntry = new KmerFeatures(maxLength);
				kmerHash.insert(std::make_pair<std::string, KmerFeatures*>(currentFwdKmer, newEntry));
			}
			

			// if(kmerHash[currentFwdKmer].second>0 && 
				// (currentLength-kmerHash[currentFwdKmer].second/kmerHash[currentFwdKmer].first)>64)
				// std::cout << currentLength-kmerHash[currentFwdKmer].second/kmerHash[currentFwdKmer].first << " haha\n";
			
			// std::cout << currentLength+1 - seedStr.length() << "\n";

			hashiter = kmerHash.find(currentFwdKmer);
			if(expectedLength<0)
				hashiter->second->add(currentLength-seedStr.length());
			else
				hashiter->second->add(expectedLength - currentLength + smallKmerSize);

			// LF mapping
            fwdIndex = m_pRBWT->getPC(b) + m_pRBWT->getOcc(b, fwdIndex - 1);			
		}
	}

	// LF-mapping of each rvc index	
	for(int64_t rvcRootIndex=rvcInterval.lower; rvcRootIndex <= rvcInterval.upper && rvcInterval.isValid(); rvcRootIndex++)
	{
		std::string currentRvcKmer = reverseComplement( seedStr.substr(seedStr.length() - smallKmerSize) );
		int64_t rvcIndex = rvcRootIndex;
		
		for(int64_t currentLength = (int64_t)seedStr.length(); currentLength <= (int64_t)maxLength; currentLength++)
		{
			char b = m_pBWT->getChar(rvcIndex);
			if(b == '$') break;
			
			currentRvcKmer = b + currentRvcKmer.substr(0, smallKmerSize-1);
			// std::cout << smallKmerSize << ":" << currentRvcKmer << "\n";

			kmerHashiter hashiter = kmerHash.find(currentRvcKmer);
			if(hashiter == kmerHash.end())
			{
				KmerFeatures *newEntry = new KmerFeatures(maxLength);
				kmerHash.insert(std::make_pair<std::string, KmerFeatures*>(currentRvcKmer, newEntry));
			}
			 
			// kmerHash[currentRvcKmer].first++;
			
			// if(expectedLength < 0)
				// kmerHash[currentRvcKmer].second += currentLength;
			// else
				// kmerHash[currentRvcKmer].second += expectedLength-currentLength+smallKmerSize;

			hashiter = kmerHash.find(currentRvcKmer);
			if(expectedLength<0)
				hashiter->second->add(currentLength-seedStr.length());
			else
				hashiter->second->add(expectedLength - currentLength + smallKmerSize);

			// LF mapping
            rvcIndex = m_pBWT->getPC(b) + m_pBWT->getOcc(b, rvcIndex - 1);
		}
	}
	
	// std::cout << kmerHash.size() << "\n";
	return true;
}


void SAIPBSelfCorrectTree::printLeaves(size_t hashKmerSize)
{
	std::cout << m_leaves.size() << ":" << m_currentLength << "\n";
	//Bamboo.PB.45195.gap.fa
	//cout << "CGCCAGCTGAGCTGGCGGTGTGAAATCAGGCAGTCAGGCGGCTCGCGTCTTGCGCGATAACCAGTTCTTCGTTGGTTGGGATAACCACCGCAGGACGGGTACCTTCTTTGTTGATGAAACCAGATTTGCCGAAACGTGCAGCCAGGTTGCGTTCATGATCAACTTCAAAGCCCAGCACGCCCAGTTTGCCCAGAGACAGTTCACGA\n";

	//PB248_9895	
	// cout << "GCGAGACGGTATTACCCGGCCCCTGGTCGCGCGGCAGGTTATGAATATTCTGTTCATGCAGGGAAAAACTCCCCGCCAGTGTAGCGATTTCACGCTCAGCAACATGGCGCGGCACACCAGCTAATAGAACTTCTCCACGCATCTGCACAATGTTCCCGCGCTCGCCAAGTTGCAAGGTGTTAAACGATGCCACGGGCGAGACTTCCGTTGCCAC\n";
	for(STNodePtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
	{
		std::string STNodeStr = (*iter)->getFullString();
		std::string fwdrepeatunit = (*iter)->getSuffix(hashKmerSize);
		kmerHashiter hashiter = kmerHash.find(fwdrepeatunit);
		
		if(hashiter!=kmerHash.end() 
			// && kmerHash[fwdrepeatunit]->getTotalFreq() > 0 
			// && kmerHash[fwdrepeatunit]->getSumOfFreq(m_currentLength-m_seedLength)>0 
		  )
		{
			std::cout << STNodeStr.substr(m_seedLength-hashKmerSize) << " " << fwdrepeatunit 
			<< ":" << kmerHash[fwdrepeatunit]->getTotalFreq() 
			<< ":" << kmerHash[fwdrepeatunit]->getSumOfFreq(m_currentLength-m_seedLength) 
			// << ":" << kmerHash[fwdrepeatunit]->getTotalSum()/kmerHash[fwdrepeatunit]->getTotalFreq() + src.length()
			// << ":" << kmerHash[fwdrepeatunit]->getSumOfPos(m_currentLength-m_seedLength)/ kmerHash[fwdrepeatunit]->getSumOfFreq(m_currentLength-m_seedLength)+m_seedLength
			;
		}
		
		std::string rvcrepeatunit = reverseComplement(fwdrepeatunit);
		kmerHashiter hashiter2 = kmerHash.find(rvcrepeatunit);
		
		
		if(hashiter2!=kmerHash.end() 
			// && kmerHash[rvcrepeatunit]->getTotalFreq() > 0 
			// && kmerHash[rvcrepeatunit]->getSumOfFreq(m_currentLength-m_seedLength)>0 
			)
		{
			std::cout << "-"  << " " << rvcrepeatunit
			<< ":" << kmerHash[rvcrepeatunit]->getTotalFreq() 
			<< ":" << kmerHash[rvcrepeatunit]->getSumOfFreq(m_currentLength-m_seedLength) 
			// << ":" << kmerHash[rvcrepeatunit]->getTotalSum()/kmerHash[rvcrepeatunit]->getTotalFreq() + m_seedLength
			// << ":" << kmerHash[rvcrepeatunit]->getSumOfPos(m_currentLength-m_seedLength)/ kmerHash[rvcrepeatunit]->getSumOfFreq(m_currentLength-m_seedLength)+m_seedLength
			;
		}
		
		if(hashiter!=kmerHash.end() || hashiter2!=kmerHash.end())
		{
			std::cout << "--" << (double)(*iter)->getKmerCount()/m_currentLength;
			std::cout << "\n"; 
		}
		
		// assert(hashiter!=kmerHash.end() || hashiter2!=kmerHash.end());
	}
}

// Print the string represented by every node
void SAIPBSelfCorrectTree::printAll()
{
    std::cout << "Print all: \n";
    m_pRootNode->printAllStrings("");
}

// Extend each leaf node
void SAIPBSelfCorrectTree::attempToExtendUsingHash(STNodePtrList &newLeaves, size_t hashKmerSize)
{
    for(STNodePtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
    {
		std::vector< std::pair<std::string, BWTIntervalPair> > extensions;
        extensions = getFMIndexRightExtensions(*iter, 1);

        // Either extend the current node or branch it
        // If no extension, do nothing and this node
        // is no longer considered a leaf
        if(extensions.size() == 1)
        {
            // Single extension, do not branch
            (*iter)->extend(extensions.front().first);
            (*iter)->fwdInterval=extensions.front().second.interval[0];
            (*iter)->rvcInterval=extensions.front().second.interval[1];
			
			std::string fwdkmer = (*iter)->getSuffix(hashKmerSize);
			double currAvgFreq = (double)(*iter)->getKmerCount()/m_currentLength;
			
			size_t kmerfreqs = 0;
			if( isExtensionValid(fwdkmer, currAvgFreq, kmerfreqs) )
			{
				(*iter)->addKmerCount(kmerfreqs);
				newLeaves.push_back(*iter);
			}
        }
        else if(extensions.size() > 1)
        {
            // Create extensions for each new base
            for(size_t i = 0; i < extensions.size(); ++i)
            {
                SAIntervalNode* pChildNode = (*iter)->createChild(extensions[i].first);
                pChildNode->fwdInterval=extensions[i].second.interval[0];
                pChildNode->rvcInterval=extensions[i].second.interval[1];
			
				std::string fwdkmer = pChildNode->getSuffix(hashKmerSize);
				double currAvgFreq = (double)(*iter)->getKmerCount()/m_currentLength;
				
				size_t kmerfreqs=0;
				if( isExtensionValid(fwdkmer, currAvgFreq, kmerfreqs) )
				{
					// inherit kmer freq from parents
					pChildNode->addKmerCount((*iter)->getKmerCount());
					pChildNode->addKmerCount(kmerfreqs);
					newLeaves.push_back(pChildNode);
				}
            }
        }
    }
}

// check if new extension satisfies frequency, position, ...
bool SAIPBSelfCorrectTree::isExtensionValid(std::string fwdkmer, double& currAvgFreq, size_t& kmerFreq)
{
		kmerHashiter iter1 = kmerHash.find(fwdkmer);			
		
		// bubble removal by removing kmer path with avg kmer freq less than previous one
		if( iter1!=kmerHash.end() && currAvgFreq < iter1->second->getMaxAvgFreq() ) 
			return false;
		
		if( iter1!=kmerHash.end() && currAvgFreq > iter1->second->getMaxAvgFreq() )
			iter1->second->setMaxAvgFreq( currAvgFreq );

		// do it again for reverse complement kmer
		std::string rvckmer = reverseComplement(fwdkmer);
		kmerHashiter iter2 = kmerHash.find(rvckmer);
		
		// bubble removal for rvc may not be consistent		
		// if( iter2!=kmerHash.end() && currAvgFreq > iter2->second->getMaxAvgFreq() )
			// iter2->second->setMaxAvgFreq( currAvgFreq );
			
		// if( iter2!=kmerHash.end() && currAvgFreq < iter2->second->getMaxAvgFreq() &&
			// iter1!=kmerHash.end() && currAvgFreq < iter1->second->getMaxAvgFreq()) 
			// return false;


		// Restricted to local kmer frequency
		kmerFreq = iter1==kmerHash.end()? 0 : iter1->second->getSumOfFreq(m_currentLength - m_seedLength);						
		kmerFreq += iter2==kmerHash.end()? 0 : iter2->second->getSumOfFreq(m_currentLength - m_seedLength);
					
		// return sum of frequency
		// kmerFreq = iter1==kmerHash.end()? 0 : iter1->second->getTotalFreq();						
		// kmerFreq += iter2==kmerHash.end()? 0 : iter2->second->getTotalFreq();
		
		if(kmerFreq >= m_min_SA_threshold)
			return true;
		else
			return false;
}

// Refine SA intervals of each leave with a new kmer
void SAIPBSelfCorrectTree::refineSAInterval(size_t newKmer)
{
    for(STNodePtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
    {
        // reset the SA intervals using original m_minOverlap
        std::string pkmer = (*iter)->getSuffix(newKmer);
        (*iter)->fwdInterval=BWTAlgorithms::findInterval(m_pRBWT, reverse(pkmer));
        (*iter)->rvcInterval=BWTAlgorithms::findInterval(m_pBWT, reverseComplement(pkmer));
    }
}

//  IntervalSizeCutoff min freq at fwd and rvc bwt
std::vector<std::pair<std::string, BWTIntervalPair> > SAIPBSelfCorrectTree::getFMIndexRightExtensions(SAIntervalNode* pNode, const size_t IntervalSizeCutoff)
{
    std::vector<std::pair<std::string, BWTIntervalPair> > out;

    for(int i = 1; i < BWT_ALPHABET::size; ++i) //i=A,C,G,T
    {
        char b = BWT_ALPHABET::getChar(i);

        //update forward Interval using extension b
        BWTInterval fwdProbe=pNode->fwdInterval;
        if(fwdProbe.isValid())
            BWTAlgorithms::updateInterval(fwdProbe,b,m_pRBWT);

        //update reverse complement Interval using extension rcb
        BWTInterval rvcProbe=pNode->rvcInterval;
		char rcb=BWT_ALPHABET::getChar(5-i);
        if(rvcProbe.isValid())
            BWTAlgorithms::updateInterval(rvcProbe,rcb,m_pBWT);

        size_t bcount = 0;
        if(fwdProbe.isValid())
            bcount += fwdProbe.size();
        if(rvcProbe.isValid())
            bcount += rvcProbe.size();

        if(bcount >= IntervalSizeCutoff)
        {
			// std::cout << m_currentKmerSize << ":" << bcount <<"\n";
            // extend to b
            std::string tmp;
            tmp.append(1,b);
            BWTIntervalPair bip;
            bip.interval[0]=fwdProbe;
            bip.interval[1]=rvcProbe;
            out.push_back(std::make_pair(tmp, bip));
        }
    }// end of ACGT

    return out;
}


bool SAIPBSelfCorrectTree::isTerminated(SAIntervalNodeResultVector& results, END target)
{
	bool found = false;
	STNodePtrList newLeaves;
    for(STNodePtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
    {
        BWTInterval currfwd=(*iter)->fwdInterval;
        BWTInterval currrvc=(*iter)->rvcInterval;

        assert(currfwd.isValid() || currrvc.isValid());

		//The current SA interval stands for a string >= terminating kmer
		//If terminating kmer is a substr, the current SA interval is a sub-interval of the terminating interval
        bool isFwdTerminated=currfwd.isValid() && currfwd.lower >= target.fwdTerminatedInterval.lower
                            && currfwd.upper <= target.fwdTerminatedInterval.upper;
        bool isRvcTerminated=currrvc.isValid() && currrvc.lower >= target.rvcTerminatedInterval.lower
                            && currrvc.upper <= target.rvcTerminatedInterval.upper;

        if(isFwdTerminated || isRvcTerminated)
        {
            std::string STNodeStr = (*iter)->getFullString();
            SAIntervalNodeResult STresult;
            STresult.thread = STNodeStr;
			STresult.SAIntervalSize = 0;
			if((*iter)->fwdInterval.isValid())
				STresult.SAIntervalSize += (*iter)->fwdInterval.size();
			if((*iter)->rvcInterval.isValid())
				STresult.SAIntervalSize += (*iter)->rvcInterval.size();
			assert(STresult.SAIntervalSize > 0);
            results.push_back(STresult);
            found =  true;
        }
		else
			newLeaves.push_back(*iter);
    }
	m_leaves.clear();
	m_leaves = newLeaves;
    return found;
}

// Check for leaves whose extension has terminated. If the leaf has
// terminated, the walked string and coverage is pushed to the result vector
bool SAIPBSelfCorrectTree::isTerminated(SAIntervalNodeResultVector& results)
{
	bool found = false;

    for(STNodePtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
    {
        BWTInterval currfwd=(*iter)->fwdInterval;
        BWTInterval currrvc=(*iter)->rvcInterval;

        assert(currfwd.isValid() || currrvc.isValid());

		//The current SA interval stands for a string >= terminating kmer
		//If terminating kmer is a substr, the current SA interval is a sub-interval of the terminating interval
        bool isFwdTerminated=currfwd.isValid() && currfwd.lower >= m_fwdTerminatedInterval.lower
                            && currfwd.upper <= m_fwdTerminatedInterval.upper;
        bool isRvcTerminated=currrvc.isValid() && currrvc.lower >= m_rvcTerminatedInterval.lower
                            && currrvc.upper <= m_rvcTerminatedInterval.upper;

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

// int SAIPBSelfCorrectTree::addHashFromPairedSeed(std::string seedStr, std::vector<std::pair<int, std::string> >  &targets, size_t hashKmerSize)
// {
	// // cut target str into shorter End str
	// std::vector<END> ENDs;
	// for(size_t i = 0 ; i < targets.size() ; i++)
	// {
		// END tmpEND;
		// int dist_between_targets = targets[i].first;
		// tmpEND.endingKmer = targets[i].second.substr(0, hashKmerSize);
		// tmpEND.maxLength = (1.4*(dist_between_targets+40)) + tmpEND.endingKmer.length()+seedStr.length();
		// tmpEND.minLength = (0.6*(dist_between_targets-40)) + tmpEND.endingKmer.length()+seedStr.length();
		// tmpEND.fwdTerminatedInterval=BWTAlgorithms::findInterval(m_pRBWT, reverse(tmpEND.endingKmer));
		// tmpEND.rvcTerminatedInterval=BWTAlgorithms::findInterval(m_pBWT, reverseComplement(tmpEND.endingKmer));
		// ENDs.push_back(tmpEND);
	// }

	// initializeSearchTree(seedStr, hashKmerSize);
	// m_fwdTerminatedInterval=ENDs[0].fwdTerminatedInterval;
    // m_rvcTerminatedInterval=ENDs[0].rvcTerminatedInterval;
	
	// int currentLength = seedStr.length();
	
	// SAIntervalNodeResultVector results;
	// while(!m_leaves.empty() && currentLength <= ENDs[ENDs.size()-1].maxLength)
    // {
		// // std::cout << m_leaves.size() << "\t" << ENDs[ENDs.size()-1].maxLength << "\n";
		// STNodePtrList newLeaves;
		// for(STNodePtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
		// {
			// std::vector< std::pair<std::string, BWTIntervalPair> > extensions;
			// extensions = getFMIndexRightExtensions(*iter, 1);

			// if(extensions.size() == 1)
			// {
				// (*iter)->extend(extensions.front().first);
				// (*iter)->fwdInterval=extensions.front().second.interval[0];
				// (*iter)->rvcInterval=extensions.front().second.interval[1];
				// newLeaves.push_back(*iter);
			// }
			// else if(extensions.size() > 1)
			// {
				// for(size_t i = 0; i < extensions.size(); ++i)
				// {
					// SAIntervalNode* pChildNode = (*iter)->createChild(extensions[i].first);
					// pChildNode->fwdInterval=extensions[i].second.interval[0];
					// pChildNode->rvcInterval=extensions[i].second.interval[1];
					// newLeaves.push_back(pChildNode);
				// }
			// }
		// }
	
		// if(!newLeaves.empty())
			// currentLength++;
		// m_leaves.clear();
		// m_leaves = newLeaves;
		
		// for(size_t i = 0 ; i < ENDs.size() ; i++)
			// if(currentLength >= ENDs[i].minLength && currentLength <= ENDs[i].maxLength)
				// isTerminated(results, ENDs[i]);
	// }
	
	// delete m_pRootNode;
	
	// int PBReads = 0;
	// if(results.size() > 0)
	// {
		// for(size_t i = 0 ; i < results.size() ; i++)
		// {
			// std::string tmpseq;
			// tmpseq = results[i].thread;
			// PBReads += results[i].SAIntervalSize;
			// //if(debug == true) std::cout << ">" << i << "_" << tmpseq.substr(seedStr.length()-hashKmerSize).length() << "\n" << tmpseq.substr(seedStr.length()-hashKmerSize) << std::endl;
			// for(size_t j = seedStr.length()-hashKmerSize+1 ; j <= tmpseq.length() - hashKmerSize && j <= ENDs[0].maxLength - hashKmerSize; j++)
			// {
				// string kmer = tmpseq.substr(j, hashKmerSize);
				// kmerHash[kmer].first+=results[i].SAIntervalSize;
			// }
		// }
		// cout << "PBReads: " << PBReads << "." << endl;
		// return 1;
	// }
	// assert(false);
	// return -1;
// }