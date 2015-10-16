///----------------------------------------------
// Copyright 2015 National Chung Cheng University
// Written by Yao-Ting Huang & Ping-Yeh Chen
// Released under the GPL
//-----------------------------------------------

#include "SAIPBSelfCTree.h"
#include "BWTAlgorithms.h"
using namespace std;

//
// Class: SAIntervalTreeForPBGap
SAIntervalTreeForPBGap::SAIntervalTreeForPBGap(const std::string* pQuery,
                       size_t minOverlap,
                       std::vector<std::pair<int, std::string> >targets,
                       size_t MaxLeaves,
                       const BWT* pBWT,
                       const BWT* pRBWT,
                       size_t SA_threshold) :
					   m_pQuery(pQuery), m_minOverlap(minOverlap),
                       m_MaxLeaves(MaxLeaves), m_pBWT(pBWT), m_pRBWT(pRBWT),
                       m_secondread(targets[0].second), m_min_SA_threshold(SA_threshold),
					   m_maxKmerCoverage(0), m_maxUsedLeaves(0), m_PBLen(0), debug(false)
{
	// Create the root node containing the seed string
    m_pRootNode = new SAIntervalNode(pQuery, NULL);
	// store initial str of root
    m_pRootNode->computeInitial(*pQuery);
    m_leaves.push_back(m_pRootNode);
    m_currentLength=pQuery->length();
	m_currentKmerSize=minOverlap;
	std::string beginningkmer=pQuery->substr(m_currentLength-m_minOverlap);
    m_pRootNode->fwdInterval=BWTAlgorithms::findInterval(m_pRBWT,reverse(beginningkmer));
    m_pRootNode->rvcInterval=BWTAlgorithms::findInterval(m_pBWT, reverseComplement(beginningkmer));
    
	for(size_t i = 0 ; i < targets.size() ; i++)
	{
		END tmpEND;
		tmpEND.endingKmer = targets[i].second.substr(0,m_minOverlap);
		tmpEND.maxLength = (1.4*(targets[i].first+40))+tmpEND.endingKmer.length()+m_currentLength;
		tmpEND.minLength = (0.6*(targets[i].first-40))+tmpEND.endingKmer.length()+m_currentLength;
		tmpEND.fwdTerminatedInterval=BWTAlgorithms::findInterval(m_pRBWT, reverse(tmpEND.endingKmer));
		tmpEND.rvcTerminatedInterval=BWTAlgorithms::findInterval(m_pBWT, reverseComplement(tmpEND.endingKmer));
		ENDs.push_back(tmpEND);
	}
	
	m_PBLen = targets[0].first + targets[0].second.length() + m_currentLength;
	m_MaxLength=ENDs[0].maxLength;
	m_MinLength=ENDs[0].minLength;
    m_fwdTerminatedInterval=ENDs[0].fwdTerminatedInterval;
    m_rvcTerminatedInterval=ENDs[0].rvcTerminatedInterval;
	
	cout << beginningkmer << ": " << m_pRootNode->fwdInterval.size()+m_pRootNode->rvcInterval.size() << ", " << ENDs[0].endingKmer << ": " << m_fwdTerminatedInterval.size()+m_rvcTerminatedInterval.size() << ". ";

	//if(targets[0].first == 169) debug = true;
}

//
SAIntervalTreeForPBGap::~SAIntervalTreeForPBGap()
{
    // Recursively destroy the tree
    delete m_pRootNode;
}

int SAIntervalTreeForPBGap::mergeTwoSeedsUsingHash(std::string &mergedseq)
{
	// build Hash
	assert(buildHash(ENDs) > 0);
	
	m_leaves.clear();
	m_leaves.push_back(m_pRootNode);

	SAIntervalNodeResultVector results;
	while(!m_leaves.empty() && m_leaves.size() <= m_MaxLeaves && m_currentLength <=m_MaxLength)
	{
		refineSAInterval(m_minOverlap-1);
		
		STNodePtrList newLeaves;
	
		// attempt to extend one base for each leave
		// ACGT-extend the leaf nodes via updating existing SA interval
		attempToExtendUsingHash(newLeaves);

		// extension succeed
		if(!newLeaves.empty())
			m_currentLength++;

		m_leaves.clear();
		m_leaves = newLeaves;
		
		// see if terminating string is reached
		if(m_currentLength >= m_MinLength && isTerminated(results))
			;
		
		/*if(debug == true)
		{
			cout << "CGCCAGCTGAGCTGGCGGTGTGAAATCAGGCAGTCAGGCGGCTCGCGTCTTGCGCGATAACCAGTTCTTCGTTGGTTGGGATAACCACCGCAGGACGGGTACCTTCTTTGTTGATGAAACCAGATTTGCCGAAACGTGCAGCCAGGTTGCGTTCATGATCAACTTCAAAGCCCAGCACGCCCAGTTTGCCCAGAGACAGTTCACGA\n";
			for(STNodePtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
			{
				std::string STNodeStr = (*iter)->getFullString();
				char DNA[4] = {'A','T','C','G'};
				//for(int i = 0; i < 4 ; i++)
				{
					//std::string fwdrepeatunit = (*iter)->getSuffix(m_minOverlap).substr(0,8) + DNA[i];
					std::string fwdrepeatunit = (*iter)->getSuffix(m_minOverlap);
					//if(m_leaves.size() <= 50)
						cout << STNodeStr << " " << fwdrepeatunit << ":" << kmerHash[fwdrepeatunit].first << endl;
				}
			}
		}*/
	}
	
	string ans;
	int diff = 100000;
	// find the path with minimum difference or maximum kmer coverage
	if(results.size()>0)
	{
		for (size_t i = 0 ; i < results.size() ;i++)
		{
			std::string tmpseq;
			// bug fix: m_secondread may be shorter than m_minOverlap
			if(m_secondread.length()>m_minOverlap)
				tmpseq=results[i].thread+m_secondread.substr(m_minOverlap);
			else
				tmpseq=results[i].thread;
				
			size_t cov=results[i].SAICoverage;
			int ansDiff = tmpseq.length() - m_PBLen;
			//if(debug == true) cout << ">" << i << endl << tmpseq << endl;
			if(abs(ansDiff) <= diff && cov > m_maxKmerCoverage)
			{
				m_maxKmerCoverage=cov;
				ans=tmpseq;
				diff = abs(ansDiff);
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
    else if(m_currentLength>m_MaxLength)
        return -2;	// exceed search depth
    else if(m_leaves.size() > m_MaxLeaves)
        return -3;	// too much repeats
	else
		return -4;
}

int SAIntervalTreeForPBGap::buildHash(vector<END> ENDs)
{
	int currentLength = m_currentLength;
	SAIntervalNodeResultVector results;
	
	while(!m_leaves.empty() && currentLength <= ENDs[ENDs.size()-1].maxLength)
    {
		STNodePtrList newLeaves;
		for(STNodePtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
		{
			std::vector< std::pair<std::string, BWTIntervalPair> > extensions;
			extensions = getFMIndexExtensionsForHash(*iter);

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
	
	int PBReads = 0;
	if(results.size() > 0)
	{
		for(size_t i = 0 ; i < results.size() ; i++)
		{
			std::string tmpseq;
			tmpseq = results[i].thread;
			PBReads += results[i].SAIntervalSize;
			//if(debug == true) std::cout << ">" << i << "_" << tmpseq.substr(m_pQuery->length()-m_minOverlap).length() << "\n" << tmpseq.substr(m_pQuery->length()-m_minOverlap) << std::endl;
			for(size_t j = m_pQuery->length()-m_minOverlap+1 ; j <= tmpseq.length() - m_minOverlap && j <= ENDs[0].maxLength - m_minOverlap; j++)
			{
				string kmer = tmpseq.substr(j, m_minOverlap);
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
void SAIntervalTreeForPBGap::printAll()
{
    std::cout << "Print all: \n";
    m_pRootNode->printAllStrings("");
}

// Extend each leaf node
void SAIntervalTreeForPBGap::attempToExtendUsingHash(STNodePtrList &newLeaves)
{
	//int diff = 25;
    for(STNodePtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
    {
		std::vector< std::pair<std::string, BWTIntervalPair> > extensions;
        extensions = getFMIndexExtensionsForHash(*iter);

        // Either extend the current node or branch it
        // If no extension, do nothing and this node
        // is no longer considered a leaf
        if(extensions.size() == 1)
        {
            // Single extension, do not branch
            (*iter)->extend(extensions.front().first);
            (*iter)->fwdInterval=extensions.front().second.interval[0];
            (*iter)->rvcInterval=extensions.front().second.interval[1];
			std::string kmer = (*iter)->getSuffix(m_minOverlap);
			int kmerfreqs = kmerHash[kmer].first;
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
				std::string kmer = pChildNode->getSuffix(m_minOverlap);
				int kmerfreqs = kmerHash[kmer].first;
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
void SAIntervalTreeForPBGap::refineSAInterval(size_t newKmer)
{
    for(STNodePtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
    {
        // reset the SA intervals using original m_minOverlap
        std::string pkmer = (*iter)->getSuffix(newKmer);
        (*iter)->fwdInterval=BWTAlgorithms::findInterval(m_pRBWT, reverse(pkmer));
        (*iter)->rvcInterval=BWTAlgorithms::findInterval(m_pBWT, reverseComplement(pkmer));
    }

	m_currentKmerSize=newKmer;
}

std::vector<std::pair<std::string, BWTIntervalPair> > SAIntervalTreeForPBGap::getFMIndexExtensionsForHash(SAIntervalNode* pNode)
{
    std::vector<std::pair<std::string, BWTIntervalPair> > out;
    size_t IntervalSizeCutoff=1;    //min freq at fwd and rvc bwt, >=3 is equal to >=2 kmer freq

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

bool SAIntervalTreeForPBGap::isTerminated(SAIntervalNodeResultVector& results, END target)
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
bool SAIntervalTreeForPBGap::isTerminated(SAIntervalNodeResultVector& results)
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