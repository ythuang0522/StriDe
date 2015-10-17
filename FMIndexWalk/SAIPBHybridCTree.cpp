///----------------------------------------------
// Copyright 2015 National Chung Cheng University
// Written by Yao-Ting Huang & Ping-Yeh Chen
// Released under the GPL
//-----------------------------------------------

#include "SAIPBHybridCTree.h"
#include "BWTAlgorithms.h"
using namespace std;

//
// Class: SAIntervalTree
SAIntervalPBHybridCTree::SAIntervalPBHybridCTree(const std::string* pQuery,
                               size_t minOverlap,
							   size_t maxOverlap,
                               int MaxLength,
                               size_t MaxLeaves,
                               const BWT* pBWT,
                               const BWT* pRBWT,
                               std::string secondread,
                               size_t SA_threshold,
                               bool KmerMode) :
                               m_pQuery(pQuery), m_minOverlap(minOverlap), m_maxOverlap(maxOverlap), m_MaxLength(MaxLength),
                               m_MaxLeaves(MaxLeaves), m_pBWT(pBWT), m_pRBWT(pRBWT),
                               m_secondread(secondread), m_min_SA_threshold(SA_threshold),
                               m_kmerMode(KmerMode), m_maxKmerCoverage(0), m_maxUsedLeaves(0)
{
    // Create the root node containing the seed string
    m_pRootNode = new SAIntervalNode(pQuery, NULL);
    m_pRootNode->computeInitial(*pQuery);   //store initial str of root
    m_leaves.push_back(m_pRootNode);

    m_currentLength=pQuery->length();
	//m_currentLength=m_minOverlap;
	m_currentKmerSize=minOverlap;
	
    //initialize the beginning SA intervals with kmer length=m_minOverlap
    std::string beginningkmer=pQuery->substr(m_currentLength-m_minOverlap);
    m_pRootNode->fwdInterval=BWTAlgorithms::findInterval(m_pRBWT,reverse(beginningkmer));
    m_pRootNode->rvcInterval=BWTAlgorithms::findInterval(m_pBWT, reverseComplement(beginningkmer));
	
    //initialize the ending SA intervals with kmer length=m_minOverlap
    std::string endingkmer=secondread.substr(0,m_minOverlap);
	m_MaxLength=(1.2*(MaxLength+10))+endingkmer.length()+m_currentLength;
	m_MinLength=(0.8*(MaxLength-10))+endingkmer.length()+m_currentLength;
    m_fwdTerminatedInterval=BWTAlgorithms::findInterval(m_pRBWT, reverse(endingkmer));
    m_rvcTerminatedInterval=BWTAlgorithms::findInterval(m_pBWT, reverseComplement(endingkmer));

	//std::cout << m_minOverlap << ":" << beginningkmer << ":" << endingkmer << "\n";
}

//
SAIntervalPBHybridCTree::~SAIntervalPBHybridCTree()
{
    // Recursively destroy the tree
    delete m_pRootNode;
}

//On success return the length of merged string
int SAIntervalPBHybridCTree::mergeTwoReads(std::string &mergedseq)
{
    SAIntervalNodeResultVector results;
	
	//BFS search from 1st to 2nd read via FM-index walk
    while(!m_leaves.empty() && m_leaves.size() <= m_MaxLeaves && m_currentLength <=m_MaxLength)
    {
        // ACGT-extend the leaf nodes via updating existing SA interval
        extendLeaves();
		
		//see if terminating string is reached
		if(m_currentLength >= m_MinLength && isTerminated(results))
			break;
		
		if(m_leaves.size()>m_maxUsedLeaves)	 m_maxUsedLeaves=m_leaves.size();
		// std::cout << m_currentKmerSize << ":" << m_currentLength << ":" << m_leaves.size() << "\n";	
    }
		
	//find the path with maximum kmer coverage
	if( results.size()>0 )
	{
		//if multiple paths are bubbles collapsing all together at terminal, reset m_maxUsedLeaves to indicate reliable path.
		if(results.size()==m_leaves.size() && m_maxUsedLeaves>1){
			m_maxUsedLeaves=1;
		}
		
		std::string tmpseq;
		for (size_t i = 0 ; i < results.size() ;i++){
			//bug fix: m_secondread may be shorter than m_minOverlap
			if(m_secondread.length()>m_minOverlap)
				tmpseq=results[i].thread+m_secondread.substr(m_minOverlap);
			else
				tmpseq=results[i].thread;
			
			// size_t cov = calculateKmerCoverage (tmpseq,m_minOverlap,m_pBWT);
			size_t cov=results[i].SAICoverage;
			if (  cov > m_maxKmerCoverage )
			{
				mergedseq=tmpseq;
				m_maxKmerCoverage=cov;
			}			
		}
		
		// std::cout << ">\n" << *m_pQuery << "\n>\n" << reverseComplement(m_secondread);
		// std::cout << mergedseq.length() << "\t" << m_maxKmerCoverage <<  "\t" << (double)m_maxKmerCoverage/mergedseq.length()  << "\n";
		return 1;
    }
	
    //Did not reach the terminal kmer
    if(m_leaves.empty())
        return -1;	//high error
    else if(m_currentLength>m_MaxLength)
        return -2;	//exceed search depth
    else if(m_leaves.size() > m_MaxLeaves)
        return -3;	//too much repeats
	else
		return -4;
}

// Print the string represented by every node
void SAIntervalPBHybridCTree::printAll()
{
    std::cout << "Print all: \n";
    m_pRootNode->printAllStrings("");
}

// Extend each leaf node

void SAIntervalPBHybridCTree::attempToExtend(STNodePtrList &newLeaves)
{
    for(STNodePtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
    {
        std::vector< std::pair<std::string, BWTIntervalPair> > extensions;
        extensions = getFMIndexExtensions(*iter);

        // Either extend the current node or branch it
        // If no extension, do nothing and this node
        // is no longer considered a leaf
        if(extensions.size() == 1)
        {
            // Single extension, do not branch
            (*iter)->extend(extensions.front().first);
            (*iter)->fwdInterval=extensions.front().second.interval[0];
            (*iter)->rvcInterval=extensions.front().second.interval[1];
			if((*iter)->fwdInterval.isValid()) (*iter)->addKmerCount( (*iter)->fwdInterval.size());
			if((*iter)->rvcInterval.isValid()) (*iter)->addKmerCount( (*iter)->rvcInterval.size());
			
            newLeaves.push_back(*iter);
        }
        else if(extensions.size() > 1)
        {
            // Branch
            for(size_t i = 0; i < extensions.size(); ++i)
            {
                SAIntervalNode* pChildNode = (*iter)->createChild(extensions[i].first);
                pChildNode->fwdInterval=extensions[i].second.interval[0];
                pChildNode->rvcInterval=extensions[i].second.interval[1];
				//inherit accumulated kmerCount from parent
				pChildNode->addKmerCount( (*iter)->getKmerCount() );
				if(pChildNode->fwdInterval.isValid()) pChildNode->addKmerCount( pChildNode->fwdInterval.size());
				if(pChildNode->rvcInterval.isValid()) pChildNode->addKmerCount( pChildNode->rvcInterval.size());

                newLeaves.push_back(pChildNode);
            }
        }
    }
}

void SAIntervalPBHybridCTree::extendLeaves()
{
    STNodePtrList newLeaves;
	
	//attempt to extend one base for each leave
    attempToExtend(newLeaves);
	if(m_kmerMode || m_currentKmerSize>=m_maxOverlap) 
		refineSAInterval(m_minOverlap);
	
    //shrink the SAIntervals in case overlap is larger than read length
    if(!m_kmerMode  &&  newLeaves.empty() )
    {
        refineSAInterval(m_minOverlap);
        attempToExtend(newLeaves);
    }

	//extension succeed
    if(!newLeaves.empty()){
        m_currentLength++;  
		m_currentKmerSize++;
	}

    m_leaves.clear();
    m_leaves = newLeaves;

}

// Check for leaves whose extension has terminated. If the leaf has
// terminated, the walked string and coverage is pushed to the result vector
bool SAIntervalPBHybridCTree::isTerminated(SAIntervalNodeResultVector& results)
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

std::vector<std::pair<std::string, BWTIntervalPair> > SAIntervalPBHybridCTree::getFMIndexExtensions(SAIntervalNode* pNode)
{
    std::vector<std::pair<std::string, BWTIntervalPair> > out;
    size_t IntervalSizeCutoff=m_min_SA_threshold;    //min freq at fwd and rvc bwt, >=3 is equal to >=2 kmer freq

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

size_t SAIntervalPBHybridCTree::calculateKmerCoverage (const std::string & seq , size_t kmerLength , const BWT* pBWT)
{
	if (seq.length()<kmerLength) return 0;
	size_t cov = 0 ;
	for (size_t i=0;i<=seq.length()-kmerLength;i+=kmerLength/2)
		cov += BWTAlgorithms::countSequenceOccurrences(seq.substr(i,kmerLength) , pBWT );
	return cov;
}


// Refine SA intervals of each leave with a new kmer
void SAIntervalPBHybridCTree::refineSAInterval(size_t newKmer)
{
    for(STNodePtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
    {
        // reset the SA intervals using original m_minOverlap
        std::string pkmer = (*iter)->getSuffix(newKmer);
        (*iter)->fwdInterval=BWTAlgorithms::findInterval(m_pRBWT, reverse(pkmer));
        (*iter)->rvcInterval=BWTAlgorithms::findInterval(m_pBWT, reverseComplement(pkmer));;
    }
	m_currentKmerSize=newKmer;
}