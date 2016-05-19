//----------------------------------------------
// Copyright 2016 National Chung Cheng University
// Written by Yao-Ting Huang & Ping-Yeh Chen
// Released under the GPL
//-----------------------------------------------

#include "SAIPBHybridCTree.h"
#include "BWTAlgorithms.h"
#include "stdaln.h"

using namespace std;

//
// Class: SAIntervalTree
SAIntervalPBHybridCTree::SAIntervalPBHybridCTree(FMWalkParameters parameters):
	m_pSourceSeed(&parameters.sourceSeed), 
	m_strBetweenSrcTarget(parameters.strBetweenSrcTarget),
	m_targetSeed(parameters.targetSeed),
	m_minOverlap(parameters.minOverlap), 
	m_maxOverlap(parameters.maxOverlap), 
	m_disBetweenSrcTarget(parameters.disBetweenSrcTarget),
	m_MaxLeaves(parameters.maxLeaves), 
	m_pBWT(parameters.indices.pBWT), 
	m_pRBWT(parameters.indices.pRBWT), 
	m_min_SA_threshold(parameters.SAThreshold),
	m_kmerMode(parameters.kmerMode),
	m_lowCoverageHighErrorMode(parameters.lowCoverageHighErrorMode),
	m_debugMode(parameters.debugMode),
	m_coverage(parameters.coverage)
{
	// Create the root node containing the seed string
	m_pRootNode = new SAIntervalNode(m_pSourceSeed, NULL);
	// store initial str of root
	m_pRootNode->computeInitial(*m_pSourceSeed);
	m_leaves.push_back(m_pRootNode);

	m_currentLength = m_pSourceSeed->length();
	m_currentKmerSize = m_minOverlap;
	
	// initialize the beginning SA intervals with kmer length = m_minOverlap
	std::string beginningkmer = m_pSourceSeed->substr(m_currentLength-m_minOverlap);
	m_pRootNode->fwdInterval = BWTAlgorithms::findInterval(m_pRBWT, reverse(beginningkmer));
	m_pRootNode->rvcInterval = BWTAlgorithms::findInterval(m_pBWT, reverseComplement(beginningkmer));

    // initialize the ending SA intervals with kmer length = m_minOverlap
	std::string endingkmer = m_targetSeed.substr(0, m_minOverlap);

	// PacBio reads are longer than real length due to insertions
	m_MaxLength = (1.1*(m_disBetweenSrcTarget+10))+endingkmer.length()+m_currentLength;
	m_MinLength = (0.9*(m_disBetweenSrcTarget-30))+endingkmer.length()+m_currentLength;
	//std::cout << m_MaxLength << ":\t" << m_MinLength;

	m_fwdTerminatedInterval = BWTAlgorithms::findInterval(m_pRBWT, reverse(endingkmer));
	m_rvcTerminatedInterval = BWTAlgorithms::findInterval(m_pBWT, reverseComplement(endingkmer));

	//m_expectedLength = m_pSourceSeed->length() + m_disBetweenSrcTarget + m_targetSeed.length();
	
	m_beginningIntervalSize = m_pRootNode->fwdInterval.size()+m_pRootNode->rvcInterval.size();
	m_terminatedIntervalSize = m_fwdTerminatedInterval.size()+m_rvcTerminatedInterval.size();
	
	if(m_debugMode)
	{
		int fwdKmerFreqs, rvcKmerFreqs, kmerFreqs;
		fwdKmerFreqs = BWTAlgorithms::countSequenceOccurrencesSingleStrand(beginningkmer, parameters.indices);
		rvcKmerFreqs = BWTAlgorithms::countSequenceOccurrencesSingleStrand(reverseComplement(beginningkmer), parameters.indices);
		kmerFreqs = fwdKmerFreqs+rvcKmerFreqs;
		std::cout << beginningkmer << " " << kmerFreqs << "\n";
		fwdKmerFreqs = BWTAlgorithms::countSequenceOccurrencesSingleStrand(endingkmer, parameters.indices);
		rvcKmerFreqs = BWTAlgorithms::countSequenceOccurrencesSingleStrand(reverseComplement(endingkmer), parameters.indices);
		kmerFreqs = fwdKmerFreqs+rvcKmerFreqs;
		std::cout << endingkmer << " " << kmerFreqs << "\n";
	}
	
	if(m_lowCoverageHighErrorMode == true)
		m_min_SA_threshold *= 3;
}

//
SAIntervalPBHybridCTree::~SAIntervalPBHybridCTree()
{
    // Recursively destroy the tree
    delete m_pRootNode;
}

// On success return the length of merged string
int SAIntervalPBHybridCTree::mergeTwoSeeds(FMWalkResult &FMWResult)
{
    SAIntervalNodeResultVector results;
	
	//if(m_debugMode)
	//	std::cout << m_disBetweenSrcTarget << " " << m_MinLength << "----\n";	
	
    size_t alignCount = 0;
    // BFS search from 1st to 2nd read via FM-index walk
    while(!m_leaves.empty() && m_leaves.size() <= m_MaxLeaves && m_currentLength <= m_MaxLength)
    {
		// ACGT-extend the leaf nodes via updating existing SA interval
        	extendLeaves();
		
		// see if terminating string is reached
		if(m_currentLength >= m_MinLength)
			isTerminated(results);
		
		if(m_debugMode)
		{	
			//std::cout << results.size() << "\n";
			for(STNodePtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
			{
				std::cout << "ATCAGGTAGTCGATACCGTACAGCGCAACGCGACGGTAGTCACCGATGATACGGCCACGGCCATATGCATCTGGCAGACCGGTCAGAACACCAGATTTACGGCAACGCAGGATGTCCGGAGTGTAAACGTCGAACACGCCCTGGTTGTGAGTTTTACGGTATTCAGTGAAGATTTTTTTGATCATCGGATCCAGTTCGCGGTTGTACGCTTTGCAGGAACCTTCGATCATTTTGATACCACCGAACGGGATAAG" << endl;
				std::cout << (*iter)->getSuffix(m_currentLength-m_pSourceSeed->length()) << " ";
				std::cout << m_currentLength << " " << m_currentKmerSize << " " << (*iter)->fwdInterval.size()+(*iter)->rvcInterval.size() << "\n";
			}
		}
		
		if(m_leaves.size() > m_maxUsedLeaves)
			m_maxUsedLeaves = m_leaves.size();

/*		// prune leaaves via alignment to raw read
		if(m_leaves.size()>m_MaxLeaves/2 && m_strBetweenSrcTarget.length()>10 && alignCount < 2){
			//std::cout << m_leaves.size() << ":";
			pruneLeavesByAlignment();
			alignCount++;
			//std::cout << m_leaves.size() << "\n";
		}
*/
		// std::cout << m_currentKmerSize << ":" << m_currentLength << ":" << m_leaves.size() << "\n";	
    }
	
	// reach the terminal kmer
	if(results.size() > 0)
	{
		// find the path with maximum match percent or kmer coverage
		return findTheBestPath(results, FMWResult);
	}
	
	// std::cout << m_currentLength << ":" << m_MaxLength << "\n";
/*
	if(m_disBetweenSrcTarget == 236) {
	//	std::cout << "CCCGAACCAATATCCCTTCAAAAGACGAAAAGTCACTTTTTCCAAAACTACAGTGTCCCTACAGTACCCCGATTATATCCCCCACTAACCTGAAACCATCATCTCTTCCACAAAACGAAAACTATTTTTTCCATTACTACAGTAATCCTACAGTACTCCTGCAGTACTCCTACAGTACTACAGCATCCCCTTCAGTACTCCTACTGTACCCCCCCCCCCCCCCCCGGCACTACTACAGTACCCCGACCACATCCCTCACTAACCTCAAACCAATATTTCTTTTAAAATACGAAAACCAATTTTTCCCAAACTACAGTACCCCTACAGTACTCCTACAGTATCCCTACAGTACTCCTACAGTACCCCAAACATATCCCCCTACTAACCCGAACCAATAT\n";
		std::cout << m_currentLength << ":" << m_MaxLength << "haha\n";
		printAll();
	}
*/
    // Did not reach the terminal kmer
    if(m_leaves.empty())
        return -1;	//high error
    else if(m_currentLength > m_MaxLength)
        return -2;	//exceed search depth
    else if(m_leaves.size() > m_MaxLeaves)
        return -3;	//too much repeats
	else
		return -4;
}

void SAIntervalPBHybridCTree::pruneLeavesByAlignment()
{
	int maxAlnScore = -100;
	STNodePtrList newLeaves;

	std::string rawStr = m_strBetweenSrcTarget.substr(0, m_currentLength-m_pSourceSeed->length());
	for(STNodePtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
	{
		std::string currStr = (*iter)->getSuffix(rawStr.length());
		AlnAln *aln_global = aln_stdaln(rawStr.c_str(), currStr.c_str(), &aln_param_pacbio, 1, 1);
		if( maxAlnScore <= aln_global->score )
		{
			maxAlnScore = aln_global->score;
			newLeaves.push_back(*iter);
		}	
		aln_free_AlnAln(aln_global);
	}

	m_leaves.clear();
	m_leaves = newLeaves;
}

int SAIntervalPBHybridCTree::findTheBestPath(SAIntervalNodeResultVector results, FMWalkResult &FMWResult)
{
	double maxKmerCoverage = 0;
	int maxAlgScore = -100;
	
	for (size_t i = 0 ; i < results.size() ;i++)
	{
		std::string candidateSeq;
		
		// bug fix: m_targetSeed may be shorter than m_minOverlap
		if(m_targetSeed.length() > m_minOverlap)
			candidateSeq = results[i].thread + m_targetSeed.substr(m_minOverlap);
		else
			candidateSeq = results[i].thread;		
		
		// find the path with maximum alignment score
		std::string pathBetweenSrcTarget = candidateSeq.substr(m_pSourceSeed->length()-10, candidateSeq.length()-m_pSourceSeed->length()-m_targetSeed.length()+20);
		AlnAln *aln_global;
		aln_global = aln_stdaln(m_strBetweenSrcTarget.c_str(), pathBetweenSrcTarget.c_str(), &aln_param_pacbio, 1, 1);
		
		/*
		if(m_debugMode)
		{
			std::cout << ">pathBetweenSrcTarget:" << i+1 << ",len:" << pathBetweenSrcTarget.length() 
				<<  ",identity:" << matchPercent 
				<< ",aln score:" << aln_global->score << "\n";
			std::cout << pathBetweenSrcTarget << "\n";
			printf("\n%s\n%s\n%s\n", aln_global->out1, aln_global->outm, aln_global->out2);
		}*/

		bool isAlgScoreBetter = maxAlgScore < aln_global->score;
		if(isAlgScoreBetter)
		{
			maxAlgScore = aln_global->score;
			FMWResult.alnScore = maxAlgScore;
			FMWResult.mergedSeq = candidateSeq;
		}
		
		aln_free_AlnAln(aln_global);
	}

	if(FMWResult.mergedSeq.length() != 0)
		return 1;
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
/*	if(m_kmerMode) 
		refineSAInterval(m_minOverlap);
	else if(m_lowCoverageHighErrorMode && m_currentKmerSize >= m_maxOverlap) 
		refineSAInterval(m_minOverlap);
	else */
	if(m_currentKmerSize >= m_maxOverlap)
	{
		/*
		if(m_minOverlap > 51)
			refineSAInterval(m_minOverlap);
		else if(m_beginningIntervalSize >= 80 && m_terminatedIntervalSize >= 80)
			refineSAInterval(81);
		else if(m_beginningIntervalSize >= 80 || m_terminatedIntervalSize >= 80)
			refineSAInterval(51);
		else
			refineSAInterval(m_minOverlap);
		*/
		if(m_beginningIntervalSize >= m_coverage*0.8 || m_terminatedIntervalSize >= m_coverage*0.8) // 256: 16, extension may exceed max leaves soon after
		//if(m_beginningIntervalSize >= 80 || m_terminatedIntervalSize >= 80)
			refineSAInterval(81);
		else
			refineSAInterval(m_minOverlap);
	}
	
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
        (*iter)->rvcInterval=BWTAlgorithms::findInterval(m_pBWT, reverseComplement(pkmer));

    }
	m_currentKmerSize=newKmer;
}
