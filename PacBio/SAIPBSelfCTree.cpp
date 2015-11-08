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
					size_t min_SA_threshold) :
                    m_pBWT(pBWT), m_pRBWT(pRBWT), m_min_SA_threshold(min_SA_threshold),
					m_PBLen(0), m_pRootNode(NULL), debug(false)
{

	//if(targets[0].first == 169) debug = true;
}

//
SAIPBSelfCorrectTree::~SAIPBSelfCorrectTree()
{
    // Recursively destroy the tree
	if(m_pRootNode!=NULL)
		delete m_pRootNode;
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
	m_PBLen = expectedLength;
	size_t currentLength = src.length();
	
	SAIntervalNodeResultVector results;
	while(!m_leaves.empty() && m_leaves.size() <= maxLeaves && currentLength <= maxLength)
	{
		// std::cout << m_leaves.size() << "\t" << currentLength << "\t" << expectedLength << "\t" << maxLength << "\n";
		refineSAInterval(hashKmerSize-1);
		
		STNodePtrList newLeaves;
	
		// attempt to extend one base for each leave
		// ACGT-extend the leaf nodes via updating existing SA interval
		attempToExtendUsingHash(newLeaves, hashKmerSize);

		// extension succeed
		if(!newLeaves.empty())
			currentLength++;

		m_leaves.clear();
		m_leaves = newLeaves;
		
		// see if terminating string is reached
		if(currentLength >= minLength)
		{
			isTerminated(results);
		}

		/*if(debug == true)*/
		// {
			// cout << "CGCCAGCTGAGCTGGCGGTGTGAAATCAGGCAGTCAGGCGGCTCGCGTCTTGCGCGATAACCAGTTCTTCGTTGGTTGGGATAACCACCGCAGGACGGGTACCTTCTTTGTTGATGAAACCAGATTTGCCGAAACGTGCAGCCAGGTTGCGTTCATGATCAACTTCAAAGCCCAGCACGCCCAGTTTGCCCAGAGACAGTTCACGA\n";
			// for(STNodePtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
			// {
				// std::string STNodeStr = (*iter)->getFullString();
				// // char DNA[4] = {'A','T','C','G'};
				// // for(int i = 0; i < 4 ; i++)
				// // {
					// // std::string fwdrepeatunit = (*iter)->getSuffix(hashKmerSize).substr(0,8) + DNA[i];
					// std::string fwdrepeatunit = (*iter)->getSuffix(hashKmerSize);
					// if(m_leaves.size() <= 50)
						// cout << STNodeStr << " " << fwdrepeatunit << ":" << kmerHash[fwdrepeatunit].first << endl;
				// // }
			// }
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
			int ansDiff = tmpseq.length() - m_PBLen;
			
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
    if(m_leaves.empty())
        return -1;	// high error
    else if(currentLength > maxLength)
        return -2;	// exceed search depth
    else if(m_leaves.size() > maxLeaves)
        return -3;	// too much repeats
	else
		return -4;
}

void SAIPBSelfCorrectTree::addHashFromSingleSeedUsingFMExtension(std::string& seedStr, size_t hashKmerSize, size_t maxLength)
{
	initializeSearchTree(seedStr, seedStr.length());

	size_t currentLength = seedStr.length();	
	while(!m_leaves.empty() && currentLength <= maxLength)
    {
		// std::cout << m_leaves.size() << "\t" << ENDs[ENDs.size()-1].maxLength << "\n";
		STNodePtrList newLeaves;
		for(STNodePtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
		{
			std::vector< std::pair<std::string, BWTIntervalPair> > extensions;
			
			// Simulate LF-mapping of each SA index by setting min frequency of 1
			extensions = getFMIndexRightExtensions(*iter, 1);

			if(extensions.size() == 1)
			{
				(*iter)->extend(extensions.front().first);
				(*iter)->fwdInterval=extensions.front().second.interval[0];
				(*iter)->rvcInterval=extensions.front().second.interval[1];
				newLeaves.push_back(*iter);

				// insert into kmer hash
				size_t kmerFreq = (*iter)->fwdInterval.isValid()?(*iter)->fwdInterval.size():0;
				kmerFreq = (*iter)->rvcInterval.isValid()?kmerFreq+(*iter)->rvcInterval.size():kmerFreq;
				assert(kmerFreq>0);
				std::string kmer = (*iter)->getSuffix(hashKmerSize);
				
				// std::cout << kmer << "\t" << kmerHash[kmer].first << "\n";
				SparseHashMap<std::string, std::pair<long long int, long long int> > :: iterator hashiter = kmerHash.find(kmer);
				if(hashiter == kmerHash.end())
				{
				  kmerHash[kmer].first = 0;
				}
				
				// hashmap[] default value of int is zero
				kmerHash[kmer].first+=kmerFreq;

			}
			else if(extensions.size() > 1)
			{
				for(size_t i = 0; i < extensions.size(); ++i)
				{
					SAIntervalNode* pChildNode = (*iter)->createChild(extensions[i].first);
					pChildNode->fwdInterval=extensions[i].second.interval[0];
					pChildNode->rvcInterval=extensions[i].second.interval[1];
					newLeaves.push_back(pChildNode);

					// insert into kmer hash
					size_t kmerFreq = (*iter)->fwdInterval.isValid()?(*iter)->fwdInterval.size():0;
					kmerFreq = (*iter)->rvcInterval.isValid()?kmerFreq+(*iter)->rvcInterval.size():kmerFreq;
					assert(kmerFreq>0);
					std::string kmer = pChildNode->getSuffix(hashKmerSize);

					// std::cout << kmer << "\t" << kmerHash[kmer].first << "\n";
					SparseHashMap<std::string, std::pair<long long int, long long int> > :: iterator hashiter = kmerHash.find(kmer);
					if(hashiter == kmerHash.end())
					{
					  kmerHash[kmer].first = 0;
					}
					kmerHash[kmer].first+=kmerFreq;
				}
			}
		}
	
		if(!newLeaves.empty())
			currentLength++;
			
		m_leaves.clear();
		m_leaves = newLeaves;
	}
	
	// m_pRootNode was pointed to root of the tree
	delete m_pRootNode;
	// cout << "CGCCAGCTGAGCTGGCGGTGTGAAATCAGGCAGTCAGGCGGCTCGCGTCTTGCGCGATAACCAGTTCTTCGTTGGTTGGGATAACCACCGCAGGACGGGTACCTTCTTTGTTGATGAAACCAGATTTGCCGAAACGTGCAGCCAGGTTGCGTTCATGATCAACTTCAAAGCCCAGCACGCCCAGTTTGCCCAGAGACAGTTCACGA\n";
	// printAll();
	// getchar();
	
}

// LF-mapping of each SA index in the interval independently using loop instead of BFS tree expansion
// Contaminated reads often are simple repeats C* or T* with large freq
// Give up this read if the freq is way too large
bool SAIPBSelfCorrectTree::addHashFromSingleSeedUsingLFMapping(std::string& seedStr, size_t hashKmerSize, size_t maxLength, int64_t contaminatedCutoff)
{
	// LF-mapping of each fwd index
    BWTInterval fwdInterval=BWTAlgorithms::findInterval(m_pRBWT, reverse(seedStr));
	BWTInterval rvcInterval=BWTAlgorithms::findInterval(m_pBWT, reverseComplement(seedStr));

	// std::cout << fwdInterval.isValid() << "\t" << fwdInterval.size() << "\t" << rvcInterval.size() <<"\n";

	// Contamination leads to large freq
	if(fwdInterval.size() >= contaminatedCutoff || rvcInterval.size() >= contaminatedCutoff) return false;

	for(int64_t fwdRootIndex = fwdInterval.lower; fwdRootIndex <= fwdInterval.upper && fwdInterval.isValid(); fwdRootIndex++)
	{
		// extract hashKmer
		std::string currentFwdKmer=seedStr.substr(seedStr.length() - hashKmerSize);
		int64_t fwdIndex = fwdRootIndex;
		
		for(size_t currentLength = seedStr.length(); currentLength <= maxLength; currentLength++)
		{
			char b = m_pRBWT->getChar(fwdIndex);
			if(b == '$') break;
			
			currentFwdKmer = currentFwdKmer.substr(1) + b;
			// std::cout << currentFwdKmer << "\n";

			// SparseHashMap<std::string, std::pair<long long int, long long int> > :: iterator hashiter = kmerHash.find(currentFwdKmer);
			// if(hashiter == kmerHash.end())
			  // kmerHash[currentFwdKmer].first = 0;

			kmerHash[currentFwdKmer].first++;

			// LF mapping
            fwdIndex = m_pRBWT->getPC(b) + m_pRBWT->getOcc(b, fwdIndex - 1);			
		}
	}

	// LF-mapping of each rvc index	
	for(int64_t rvcRootIndex=rvcInterval.lower; rvcRootIndex <= rvcInterval.upper && rvcInterval.isValid(); rvcRootIndex++)
	{
		std::string currentRvcKmer=reverseComplement(seedStr.substr(seedStr.length() - hashKmerSize));
		int64_t rvcIndex = rvcRootIndex;
		
		for(size_t currentLength = seedStr.length(); currentLength <= maxLength; currentLength++)
		{
			char b = m_pBWT->getChar(rvcIndex);
			if(b == '$') break;
			
			currentRvcKmer = b + currentRvcKmer.substr(0, hashKmerSize-1);
			// std::cout << hashKmerSize << ":" << currentRvcKmer << "\n";

			// SparseHashMap<std::string, std::pair<long long int, long long int> > :: iterator hashiter = kmerHash.find(currentRvcKmer);
			// if(hashiter == kmerHash.end())
			  // kmerHash[currentRvcKmer].first = 0;

			kmerHash[currentRvcKmer].first++;
			
			// LF mapping
            rvcIndex = m_pBWT->getPC(b) + m_pBWT->getOcc(b, rvcIndex - 1);
		}
	}
	
	return true;
}

int SAIPBSelfCorrectTree::addHashFromPairedSeed(std::string seedStr, std::vector<std::pair<int, std::string> >  &targets, size_t hashKmerSize)
{
	// cut target str into shorter End str
	std::vector<END> ENDs;
	for(size_t i = 0 ; i < targets.size() ; i++)
	{
		END tmpEND;
		int dist_between_targets = targets[i].first;
		tmpEND.endingKmer = targets[i].second.substr(0, hashKmerSize);
		tmpEND.maxLength = (1.4*(dist_between_targets+40)) + tmpEND.endingKmer.length()+seedStr.length();
		tmpEND.minLength = (0.6*(dist_between_targets-40)) + tmpEND.endingKmer.length()+seedStr.length();
		tmpEND.fwdTerminatedInterval=BWTAlgorithms::findInterval(m_pRBWT, reverse(tmpEND.endingKmer));
		tmpEND.rvcTerminatedInterval=BWTAlgorithms::findInterval(m_pBWT, reverseComplement(tmpEND.endingKmer));
		ENDs.push_back(tmpEND);
	}

	initializeSearchTree(seedStr, hashKmerSize);
	m_fwdTerminatedInterval=ENDs[0].fwdTerminatedInterval;
    m_rvcTerminatedInterval=ENDs[0].rvcTerminatedInterval;
	
	int currentLength = seedStr.length();
	
	SAIntervalNodeResultVector results;
	while(!m_leaves.empty() && currentLength <= ENDs[ENDs.size()-1].maxLength)
    {
		// std::cout << m_leaves.size() << "\t" << ENDs[ENDs.size()-1].maxLength << "\n";
		STNodePtrList newLeaves;
		for(STNodePtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
		{
			std::vector< std::pair<std::string, BWTIntervalPair> > extensions;
			extensions = getFMIndexRightExtensions(*iter, 1);

			if(extensions.size() == 1)
			{
				(*iter)->extend(extensions.front().first);
				(*iter)->fwdInterval=extensions.front().second.interval[0];
				(*iter)->rvcInterval=extensions.front().second.interval[1];
				newLeaves.push_back(*iter);
			}
			else if(extensions.size() > 1)
			{
				for(size_t i = 0; i < extensions.size(); ++i)
				{
					SAIntervalNode* pChildNode = (*iter)->createChild(extensions[i].first);
					pChildNode->fwdInterval=extensions[i].second.interval[0];
					pChildNode->rvcInterval=extensions[i].second.interval[1];
					newLeaves.push_back(pChildNode);
				}
			}
		}
	
		if(!newLeaves.empty())
			currentLength++;
		m_leaves.clear();
		m_leaves = newLeaves;
		
		for(size_t i = 0 ; i < ENDs.size() ; i++)
			if(currentLength >= ENDs[i].minLength && currentLength <= ENDs[i].maxLength)
				isTerminated(results, ENDs[i]);
	}
	
	delete m_pRootNode;
	
	int PBReads = 0;
	if(results.size() > 0)
	{
		for(size_t i = 0 ; i < results.size() ; i++)
		{
			std::string tmpseq;
			tmpseq = results[i].thread;
			PBReads += results[i].SAIntervalSize;
			//if(debug == true) std::cout << ">" << i << "_" << tmpseq.substr(seedStr.length()-hashKmerSize).length() << "\n" << tmpseq.substr(seedStr.length()-hashKmerSize) << std::endl;
			for(size_t j = seedStr.length()-hashKmerSize+1 ; j <= tmpseq.length() - hashKmerSize && j <= ENDs[0].maxLength - hashKmerSize; j++)
			{
				string kmer = tmpseq.substr(j, hashKmerSize);
				kmerHash[kmer].first+=results[i].SAIntervalSize;
			}
		}
		cout << "PBReads: " << PBReads << "." << endl;
		return 1;
	}
	assert(false);
	return -1;
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
			std::string rvckmer = reverseComplement(fwdkmer);
			
			// SparseHashMap<std::string, std::pair<long long int, long long int> > :: iterator kmeriter = kmerHash.find(kmer);
			// if(kmeriter == kmerHash.end())
			  // continue;

			size_t kmerfreqs = kmerHash[fwdkmer].first + kmerHash[rvckmer].first;
			// std::cout << kmerHash[fwdkmer].first << "\t" << kmerHash[rvckmer].first << "\n";
			// getchar();
			if(kmerfreqs >= m_min_SA_threshold)
			{
				(*iter)->addKmerCount(kmerfreqs);
				newLeaves.push_back(*iter);
			}
        }
        else if(extensions.size() > 1)
        {
            // Branch
            for(size_t i = 0; i < extensions.size(); ++i)
            {
                SAIntervalNode* pChildNode = (*iter)->createChild(extensions[i].first);
                pChildNode->fwdInterval=extensions[i].second.interval[0];
                pChildNode->rvcInterval=extensions[i].second.interval[1];
			
				std::string fwdkmer = (*iter)->getSuffix(hashKmerSize);
				std::string rvckmer = reverseComplement(fwdkmer);
				
				// SparseHashMap<std::string, std::pair<long long int, long long int> > :: iterator kmeriter = kmerHash.find(kmer);
				// if(kmeriter == kmerHash.end())
				  // continue;

				size_t kmerfreqs = kmerHash[fwdkmer].first + kmerHash[rvckmer].first;

				// std::cout << m_min_SA_threshold << ": " << kmerHash[fwdkmer].first << "\t" << kmerHash[rvckmer].first << "\n";
				// getchar();
				
				if(kmerfreqs >= m_min_SA_threshold)
				{
					pChildNode->addKmerCount((*iter)->getKmerCount());
					pChildNode->addKmerCount(kmerfreqs);
					newLeaves.push_back(pChildNode);
				}
            }
        }
    }
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