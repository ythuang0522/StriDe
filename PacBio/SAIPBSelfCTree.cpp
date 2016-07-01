///----------------------------------------------
// Copyright 2015 National Chung Cheng University
// Written by Yao-Ting Huang & Ping-Yeh Chen
// Released under the GPL
//-----------------------------------------------

#include "SAIPBSelfCTree.h"
#include "BWTAlgorithms.h"
#include "ssw_cpp.h"
#include "stdaln.h"
using namespace std;

//
// Class: SAIPBSelfCorrectTree
SAIPBSelfCorrectTree::SAIPBSelfCorrectTree(
					const BWT* pBWT,
                    			const BWT* pRBWT,
					std::string rawSeq,
					size_t min_SA_threshold,
					int maxLeavesAllowed) :
                    m_pBWT(pBWT), m_pRBWT(pRBWT), m_rawSeq(rawSeq), m_min_SA_threshold(min_SA_threshold), m_maxLeavesAllowed(maxLeavesAllowed),
					m_expectedLength(0), m_currentLength(0), m_pRootNode(NULL), m_isSourceRepeat(false), m_isTargetRepeat(false)
{
	// initialize google dense hashmap
	kmerHash.set_empty_key("");
	kmerHash.resize(6000);
}

//
SAIPBSelfCorrectTree::~SAIPBSelfCorrectTree()
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
	size_t maxUsedLeaves=0;
	
	SAIntervalNodeResultVector results;
	while(!m_leaves.empty() && m_leaves.size() <= maxLeaves &&  (size_t)m_currentLength <= maxLength)
	{
		// std::cout << m_leaves.size() << "\t" << currentLength << "\t" << expectedLength << "\t" << maxLength << "\n";
		refineSAInterval(hashKmerSize-1);
		
		STNodePtrList newLeaves;
	
		// attempt to extend one base for each leave
		// ACGT-extend the leaf nodes via updating existing SA interval
		attempToExtendUsingHash(newLeaves, hashKmerSize);

		if(m_leaves.size()> maxUsedLeaves) maxUsedLeaves=m_leaves.size();
		
		// extension succeed
		if(!newLeaves.empty())
			m_currentLength++;

		m_leaves.clear();
		m_leaves = newLeaves;
		
		// see if terminating string is reached
		if( (size_t) m_currentLength >= minLength)
			isTerminated(results);

		// if(src.length()==198)
			printLeaves(hashKmerSize);
	}
	
	if(results.size()>0)
	{
		double maxKmerCoverage = 0;
		int minLengthDiff = 100000;
		double maxMatchPercent = -100;
		for (size_t i = 0 ; i < results.size() ;i++)
		{
			std::string tmpseq;
			// bug fix: dest may be shorter than hashKmerSize
			if(dest.length() > hashKmerSize)
				tmpseq=results[i].thread + dest.substr(hashKmerSize);
			else
				tmpseq=results[i].thread;

			int currLengthDiff = std::abs((int)tmpseq.length()-(int)expectedLength);

			/**** Math bug: the division requires denominator larger than numerator for accurate average frequency
					numerator denominator division
					100	1000	0.1
					99	999	0.099099099
					100	10	10
					99	9	11
			***/			
			double avgCov = (double)results[i].SAICoverage /(tmpseq.length()+1000000);
		
			// if(src.length() == 198 )
			// {
				// std::cout << ">" << currLengthDiff << "\t" << avgCov <<  "\t" << tmpseq.length() << "\t" << expectedLength << "\n";
				// std::cout << ">" << src.length() << "\n" << tmpseq.substr(src.length()) <<  "\n";
			// }
			
			// PB135123_7431.fa
			// Rule out apparent errors due to tandem repeats with large length difference
			bool isLengthDiffBetter = currLengthDiff < minLengthDiff && std::abs(currLengthDiff - minLengthDiff) >3;
			
			// But the length diff is not informative for small indels, where kmer coverage is more informative
			bool isKmerCoverageBetter = std::abs(currLengthDiff - minLengthDiff) <=3 && (maxKmerCoverage<avgCov);
			
			// sequency identify used for repeat extensions
			if(results.size()>1)
			{
				/*
				// Declares a default Aligner
  				StripedSmithWaterman::Aligner aligner;
  				// Declares a default filter
  				StripedSmithWaterman::Filter filter(0,0,0,300);
  				// Declares an alignment that stores the result
  				StripedSmithWaterman::Alignment alignment;
  				// Aligns the query to the ref
  				aligner.Align(tmpseq.c_str(), m_rawSeq.c_str(), m_rawSeq.size(), filter, &alignment);

				//std::cout << alignment.sw_score << "\n";
				*/

				// We want to compute the total matches and percent
				int matchLen = 0;
				AlnAln *aln_global;
				aln_global = aln_stdaln(m_rawSeq.c_str(), tmpseq.c_str(), &aln_param_pacbio, 1, 1);
			
				// Calculate the alignment patterns
				for(int i = 0 ; aln_global->outm[i] != '\0' ; i++)
					if(aln_global->outm[i] == '|')
						matchLen++;
				aln_free_AlnAln(aln_global);
				double matchPercent = (double)matchLen / m_rawSeq.length();
			
			/*
			if(m_debugMode)
			{
				std::cout << ">pathBetweenSrcTarget:" << i+1 << ",len:" << pathBetweenSrcTarget.length() <<  ",identity:" << matchPercent << "\n";
				std::cout << pathBetweenSrcTarget << "\n";
			}
			*/
			
				bool isPercentMatchBetter = (maxMatchPercent < matchPercent);
				if(isPercentMatchBetter)
				{
					maxMatchPercent = matchPercent;
					mergedseq = tmpseq;
				}

/*
				if((alignment.sw_score > bestAlignmentScore) || 
				   (alignment.sw_score == bestAlignmentScore && (isLengthDiffBetter || isKmerCoverageBetter)))
				{
					bestAlignmentScore = alignment.sw_score;
					minLengthDiff = currLengthDiff;
                                	maxKmerCoverage=avgCov;
					mergedseq = tmpseq;
				}
*/
			}
			else if(isLengthDiffBetter || isKmerCoverageBetter)
			{
				minLengthDiff = currLengthDiff;
				maxKmerCoverage=avgCov;
				mergedseq=tmpseq;
			}
		}
		
		// if( (results.size()>=6 && maxUsedLeaves >=16) && mergedseq.length()-src.length()>100 
			// && hashKmerSize <=13 && ( (m_isSourceRepeat && m_sourceKmerSize>=17) || (m_isTargetRepeat && m_targetKmerSize>=17) ) )
		// {
			// std::cout << ">Selected: " << results.size() << "\t" << maxKmerCoverage << "\t" << minLengthDiff << "\t" << maxUsedLeaves <<"\n";
			// std::cout << mergedseq.substr(src.length()-hashKmerSize) << "\n";
			// std::cout << mergedseq.length()-src.length() << "\n";
		// }
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
size_t SAIPBSelfCorrectTree::addHashBySingleSeed(std::string& seedStr, size_t largeKmerSize, size_t smallKmerSize, size_t maxLength, bool skipRepeat, int expectedLength)
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

	std::cout << initKmer << "\t" << smallKmerSize << "\t" << reverseComplement(initKmer) << "\t" << fwdInterval.size() << "\t" << rvcInterval.size() << "\t" << kmerFreq << "\n";

	
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
		
		// extend the fwdIndex via LF mapping
		int64_t fwdIndex = fwdRootIndex;
		for(int64_t currentLength = (int64_t)seedStr.length(); currentLength <= (int64_t)maxLength; currentLength++)
		{
			char b = m_pRBWT->getChar(fwdIndex);
			if(b == '$') break;
			
			currentFwdKmer = currentFwdKmer.substr(1) + b;
			// std::cout << currentLength << ":" << currentFwdKmer << "\n";

			insertKmerToHash(currentFwdKmer, seedStr.length(), currentLength, smallKmerSize, maxLength, expectedLength);

			// LF mapping
            fwdIndex = m_pRBWT->getPC(b) + m_pRBWT->getOcc(b, fwdIndex - 1);			
		}
	}

	
	// LF-mapping of each rvc index	
	for(int64_t rvcRootIndex=rvcInterval.lower; 
		rvcRootIndex <= rvcInterval.upper && rvcInterval.isValid() && rvcRootIndex - rvcInterval.lower < maxIntervalSize; 
		rvcRootIndex++)
	{
		std::string currentRvcKmer = reverseComplement( seedStr.substr(seedStr.length() - smallKmerSize) );
		insertKmerToHash(currentRvcKmer, seedStr.length(), seedStr.length(), smallKmerSize, maxLength, expectedLength);
		
		int64_t rvcIndex = rvcRootIndex;
		
		for(int64_t currentLength = (int64_t)seedStr.length(); currentLength <= (int64_t)maxLength; currentLength++)
		{
			char b = m_pBWT->getChar(rvcIndex);
			if(b == '$') break;
			
			currentRvcKmer = b + currentRvcKmer.substr(0, smallKmerSize-1);
			insertKmerToHash(currentRvcKmer, seedStr.length(), currentLength, smallKmerSize, maxLength, expectedLength);
			// std::cout << currentLength << ":" << currentRvcKmer << "\n";

			// LF mapping
            rvcIndex = m_pBWT->getPC(b) + m_pBWT->getOcc(b, rvcIndex - 1);
		}
	}
	
	std::cout << kmerHash.size() << "\n";
	return kmerFreq;
}

void SAIPBSelfCorrectTree::insertKmerToHash(std::string& insertedKmer, size_t seedStrLen, size_t currentLength, size_t smallKmerSize, size_t maxLength, int expectedLength)
{
	kmerHashiter hashiter = kmerHash.find(insertedKmer);
	
	if(hashiter == kmerHash.end())
	{
		KmerFeatures *newEntry;
		// source to target
		if(expectedLength<0)
			newEntry = new KmerFeatures(currentLength - seedStrLen, maxLength);
		// target to source
		else
			newEntry = new KmerFeatures(expectedLength - currentLength + smallKmerSize, maxLength);
			
		// kmerHash.insert(std::make_pair<std::string, KmerFeatures*>(insertedKmer, newEntry));
		// C++ 11
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

void SAIPBSelfCorrectTree::printLeaves(size_t hashKmerSize)
{
	std::cout << m_leaves.size() << ":" << m_currentLength << "\n";
	cout << "CTCAAGATATTTCTTCCATCATGCAAAAAAAATTTGCAGTGCATGATGTTAATCATAAATGTCGGTGTCATCATGCGCTACGCTCTATGGCTCCCTGACGTTTTTTTAGC\n";
	for(STNodePtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
	{
		std::string STNodeStr = (*iter)->getFullString();
		std::string fwdKmer = (*iter)->getSuffix(hashKmerSize);
		kmerHashiter hashiter = kmerHash.find(fwdKmer);
		
		if(hashiter!=kmerHash.end() 
			// && kmerHash[fwdKmer]->getTotalFreq() > 0 
			// && kmerHash[fwdKmer]->getSumOfFreq(m_currentLength-m_seedLength)>0 
		  )
		{
			std::cout << STNodeStr.substr(m_seedLength-hashKmerSize) << " " << fwdKmer 
			<< ":" << kmerHash[fwdKmer]->getTotalFreq() 
			// local kmer frequency at m_currentLength-m_seedLength
			<< ":" << kmerHash[fwdKmer]->getSumOfFreq(m_currentLength-m_seedLength) 
			// << ":" << kmerHash[fwdKmer]->getTotalSum()/kmerHash[fwdKmer]->getTotalFreq() + src.length()
			// << ":" << kmerHash[fwdKmer]->getSumOfPos(m_currentLength-m_seedLength)/ kmerHash[fwdKmer]->getSumOfFreq(m_currentLength-m_seedLength)+m_seedLength
			;
		}
		
		std::string rvcKmer = reverseComplement(fwdKmer);
		kmerHashiter hashiter2 = kmerHash.find(rvcKmer);
		
		
		if(hashiter2!=kmerHash.end() 
			// && kmerHash[rvcKmer]->getTotalFreq() > 0 
			// && kmerHash[rvcKmer]->getSumOfFreq(m_currentLength-m_seedLength)>0 
			)
		{
			std::cout << "-"  << " " << rvcKmer
			<< ":" << kmerHash[rvcKmer]->getTotalFreq() 
			<< ":" << kmerHash[rvcKmer]->getSumOfFreq(m_currentLength-m_seedLength) 
			// << ":" << kmerHash[rvcKmer]->getTotalSum()/kmerHash[rvcKmer]->getTotalFreq() + m_seedLength
			// << ":" << kmerHash[rvcKmer]->getSumOfPos(m_currentLength-m_seedLength)/ kmerHash[rvcKmer]->getSumOfFreq(m_currentLength-m_seedLength)+m_seedLength
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
				
				/**** Math bug: the division requires denominator larger than numerator for accurate average frequency
					numerator denominator division
					100	1000	0.1
					99	999	0.099099099
					100	10	10
					99	9	11
				***/
				// double currAvgFreq = (double)(*iter)->getKmerCount()/m_currentLength;
				double currAvgFreq = (double)(*iter)->getKmerCount()/ (m_currentLength+1000000);
				
				size_t kmerfreqs=0;
				if( isExtensionValid(fwdkmer, currAvgFreq, kmerfreqs) )
				{
					// inherit kmer freq from parents
					// pChildNode->addKmerCount((*iter)->getKmerCount());
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
		// if( iter1!=kmerHash.end() && currAvgFreq < iter1->second->getMaxAvgFreq() ) 
		
		// Bubble removal preoduces false removal, perform only when leaves are getting larger than 8
		if( iter1!=kmerHash.end() && m_leaves.size()>8 && currAvgFreq < iter1->second->getMaxAvgFreq() ) 
			return false;
		
		if( iter1!=kmerHash.end() && currAvgFreq > iter1->second->getMaxAvgFreq() )
			iter1->second->setMaxAvgFreq( currAvgFreq );

		// // do it again for reverse complement kmer
		std::string rvckmer = reverseComplement(fwdkmer);
		kmerHashiter iter2 = kmerHash.find(rvckmer);
		
		// Restricted to local kmer frequency
		kmerFreq = iter1==kmerHash.end()? 0 : iter1->second->getSumOfFreq(m_currentLength - m_seedLength);						
		kmerFreq += iter2==kmerHash.end()? 0 : iter2->second->getSumOfFreq(m_currentLength - m_seedLength);
					
		// std::cout << fwdkmer <<": " << kmerFreq << "\n";
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


// Check for leaves whose extension has terminated. If the leaf has
// terminated, the walked string and coverage is pushed to the result vector
bool SAIPBSelfCorrectTree::isTerminated(SAIntervalNodeResultVector& results)
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

// bool SAIPBSelfCorrectTree::isTerminated(SAIntervalNodeResultVector& results, END target)
// {
	// bool found = false;
	// STNodePtrList newLeaves;
    // for(STNodePtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
    // {
        // BWTInterval currfwd=(*iter)->fwdInterval;
        // BWTInterval currrvc=(*iter)->rvcInterval;

        // assert(currfwd.isValid() || currrvc.isValid());

		// //The current SA interval stands for a string >= terminating kmer
		// //If terminating kmer is a substr, the current SA interval is a sub-interval of the terminating interval
        // bool isFwdTerminated=currfwd.isValid() && currfwd.lower >= target.fwdTerminatedInterval.lower
                            // && currfwd.upper <= target.fwdTerminatedInterval.upper;
        // bool isRvcTerminated=currrvc.isValid() && currrvc.lower >= target.rvcTerminatedInterval.lower
                            // && currrvc.upper <= target.rvcTerminatedInterval.upper;

        // if(isFwdTerminated || isRvcTerminated)
        // {
            // std::string STNodeStr = (*iter)->getFullString();
            // SAIntervalNodeResult STresult;
            // STresult.thread = STNodeStr;
			// STresult.SAIntervalSize = 0;
			// if((*iter)->fwdInterval.isValid())
				// STresult.SAIntervalSize += (*iter)->fwdInterval.size();
			// if((*iter)->rvcInterval.isValid())
				// STresult.SAIntervalSize += (*iter)->rvcInterval.size();
			// assert(STresult.SAIntervalSize > 0);
            // results.push_back(STresult);
            // found =  true;
			
        // }
		// // else
			// // newLeaves.push_back(*iter);
    // }

	// // m_leaves.clear();
	// // m_leaves = newLeaves;
    // return found;
// }
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
