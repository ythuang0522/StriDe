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
SAIntervalPBHybridCTree::SAIntervalPBHybridCTree(FMWalkParameters& parameters):
	m_pStrSourceSeed(&parameters.strSourceSeed),
	m_strTargetSeed(parameters.strTargetSeed),
	m_disBetweenSrcTarget(parameters.disBetweenSrcTarget), 
	m_rawPBStrBetweenSrcTargetWith2Minoverlap(parameters.rawPBStrBetweenSrcTargetWith2Minoverlap),
	m_minOverlap(parameters.minOverlap), 
	m_maxOverlap(parameters.maxOverlap), 
	m_MaxLeaves(parameters.maxLeaves), 
	m_pBWT(parameters.indices.pBWT), 
	m_pRBWT(parameters.indices.pRBWT), 
	m_minSAThreshold(parameters.SAThreshold),
	m_kmerMode(parameters.kmerMode),
	m_lowCoverageHighErrorMode(parameters.lowCoverageHighErrorMode),
	m_debugMode(parameters.debugMode)
{
	// Create the root node containing the seed string
	m_pRootNode = new SAIntervalNode(m_pStrSourceSeed, NULL);
	// store initial str of root
	m_pRootNode->computeInitial(*m_pStrSourceSeed);
	m_leaves.push_back(m_pRootNode);

	m_currentLength = m_pStrSourceSeed->length();
	m_currentKmerSize = m_minOverlap;
	
	// initialize the beginning SA intervals with kmer length = m_minOverlap
	std::string beginningkmer = m_pStrSourceSeed->substr(m_currentLength-m_minOverlap);
	m_pRootNode->fwdInterval = BWTAlgorithms::findInterval(m_pRBWT, reverse(beginningkmer));
	m_pRootNode->rvcInterval = BWTAlgorithms::findInterval(m_pBWT, reverseComplement(beginningkmer));

    // initialize the ending SA intervals with kmer length = m_minOverlap
	std::string endingkmer = m_strTargetSeed.substr(0, m_minOverlap);

	// PacBio reads are longer than real length due to insertions
	m_MaxLength = (1.18*(m_disBetweenSrcTarget+18))+endingkmer.length()+m_currentLength;
	m_MinLength = (0.8*(m_disBetweenSrcTarget-20))+endingkmer.length()+m_currentLength;
	if(m_debugMode)
		cout << m_disBetweenSrcTarget << "\n" << m_rawPBStrBetweenSrcTargetWith2Minoverlap << endl;
	// Errors were characterized by comparison to the known reference sequences, 
	// showing that the primary error mode was insertions at 12%, 
	// followed by deletions at 2%,  
	// and apparent mismatch errors at 1%. 
	// ref: https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-13-375
	// m_expectedLength = (0.9*m_disBetweenSrcTarget)+endingkmer.length()+m_currentLength;
	
	m_fwdTerminatedInterval = BWTAlgorithms::findInterval(m_pRBWT, reverse(endingkmer));
	m_rvcTerminatedInterval = BWTAlgorithms::findInterval(m_pBWT, reverseComplement(endingkmer));

	m_iniMinSAThreshold = parameters.SAThreshold;
	
	if(m_debugMode)
	{
		std::cout << beginningkmer << " " << m_beginningIntervalSize << "--\n";
		std::cout << endingkmer << " " << m_terminatedIntervalSize << "--\n";
	}
	
	if(m_lowCoverageHighErrorMode == true)
		m_minSAThreshold *= 3;
}

//
SAIntervalPBHybridCTree::~SAIntervalPBHybridCTree()
{
    // Recursively destroy the tree
    delete m_pRootNode;
}

// On success return the length of merged string
void SAIntervalPBHybridCTree::mergeTwoSeeds(FMWalkResult &FMWResult)
{
    SAIntervalNodeResultVector results;
	
	//if(m_debugMode)
	//	std::cout << m_disBetweenSrcTarget << " " << m_MinLength << "----\n";	
	
    // BFS search from 1st to 2nd read via FM-index walk
    while(!m_leaves.empty() && m_leaves.size() <= m_MaxLeaves && m_currentLength <= m_MaxLength)
    {
		// ACGT-extend the leaf nodes via updating existing SA interval
        	// extendLeaves();
			extendLeaves_v2();
		
		// see if terminating string is reached
		if(m_currentLength >= m_MinLength)
			isTerminated(results);
		// if(m_currentLength >= m_MinLength && isTerminated(results))
			// break;
		
		if(m_debugMode)
		{	
			std::cout << "----m_currentKmerSize: " << m_currentKmerSize
			<< ", m_leaves.size(): " << m_leaves.size()
			<< ", result: " << results.size() << "----\n";
			std::cout << "GATAGCGAACGCCCACTTTCACGCTCAAACAATAACCAAGAACCTGTGCTATGGAATTTTAAATCACTAGCAAACAATTTCAACGTTAGGTTTACGTACCATTTTCACGCCACATCGTCATCCTCTAAGATTGAGACGTATTTTCAGTTTCTAAACGATTATCTAGCGGAAAACCTATACAAGTGCATCAACATTTTTCATGATGACTGTAATGGGTTGACGAAGCCAGTTATTCATGAACAATTTA-TTAAT-TACGTCTTACAACCCATTAGGGATAAAGTAAGATCCACCCTATTTCAAAACGATTTGAAAACTTTGATCGTCCTAATTTCCCAAATCCTGGCTACAGACAAAAATTTATTGAATTCTTTTCATTACCATGGGCTAGGTTTGGTGTCGTTAATTTCCGATGAAGTATGGGAGAAATGGATCAACTATGAAGTTGAAATGGCCAATAGGCAATTCATCAATATAACTAAAAATCCGGAAGATTTCCCAAAATCTTCTCAGAATTTTGTCAAATTAATCAATAAAATTTACGATTATT" << endl;
			for(STNodePtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
			{
				std::cout << (*iter)->getSuffix(m_currentLength-m_pStrSourceSeed->length()) << " ";
				std::cout << m_currentLength << " " << (*iter)->fwdInterval.size()+(*iter)->rvcInterval.size() << " " << (*iter)->getKmerCount() << "\n";
			}
		}
		
		if(m_leaves.size() > m_maxUsedLeaves)
			m_maxUsedLeaves = m_leaves.size();

		// std::cout << m_currentKmerSize << ":" << m_currentLength << ":" << m_leaves.size() << "\n";	
    }
	// if(m_debugMode)
		// std::cout << results.size() << "\n";
	
	if(results.size() > 0)
	{
		// find the path with maximum match percent or kmer coverage
		return findTheBestPath(results, FMWResult);
	}

    // did not reach the terminal kmer
    if(m_leaves.empty())
        FMWResult.typeFMWalkResult = -1;	//high error
    else if(m_currentLength > m_MaxLength)
        FMWResult.typeFMWalkResult = -2;	//exceed search depth
    else if(m_leaves.size() > m_MaxLeaves)
        FMWResult.typeFMWalkResult = -3;	//too much repeats
	else
		FMWResult.typeFMWalkResult = -4;
	return;
}

// Print the string represented by every node
void SAIntervalPBHybridCTree::printAll()
{
    std::cout << "Print all: \n";
    m_pRootNode->printAllStrings("");
}

void SAIntervalPBHybridCTree::findTheBestPath(SAIntervalNodeResultVector results, FMWalkResult &FMWResult)
{
	int maxAlgScore = -100;
	// size_t maxKmerCoverage = 0;
	
	for (size_t i = 0 ; i < results.size() ; i++)
	{
		std::string candidateSeq = results.at(i).thread;
		std::string pathBetweenSrcTargetWith2Minoverlap = 
			candidateSeq.substr(m_pStrSourceSeed->length()-m_minOverlap);
		
		// find the path with maximum alignment score
		AlnAln *aln_global;
		aln_global = aln_stdaln(m_rawPBStrBetweenSrcTargetWith2Minoverlap.c_str(), pathBetweenSrcTargetWith2Minoverlap.c_str(), &aln_param_pacbio, 1, 1);
		
		// if(m_debugMode)
		// {
			// std::cout << ">pathBetweenSrcTargetWith2Minoverlap:" << i+1 << ",len:" << pathBetweenSrcTargetWith2Minoverlap.length() 
				// <<  ",identity:" << matchPercent 
				// << ",alnScore:" << aln_global->score << ",kmerFreq:" << results.at(i).SAICoverage <<"\n";
			// std::cout << pathBetweenSrcTargetWith2Minoverlap << "\n";
			// printf("\n%s\n%s\n%s\n", aln_global->out1, aln_global->outm, aln_global->out2);
			// cout << ">" << i << "\n" << pathBetweenSrcTargetWith2Minoverlap.c_str() << "\n";
		// }

		bool isAlgScoreBetter = maxAlgScore < aln_global->score;
		// bool isAlgScoreBetter = maxKmerCoverage < results.at(i).SAICoverage;
		if(isAlgScoreBetter)
		{
			maxAlgScore = aln_global->score;
			// maxKmerCoverage = results.at(i).SAICoverage;
			FMWResult.alnScore = maxAlgScore;
			// FMWResult.kmerFreq = maxKmerCoverage;
			FMWResult.mergedSeq = candidateSeq;
		}
		
		aln_free_AlnAln(aln_global);
	}

	if(FMWResult.mergedSeq.length() != 0)
	{	
		// if(m_debugMode)
		// {
			// std::cout << ">resultStr:YO,len:" << FMWResult.mergedSeq.length() 
				// <<  ",identity:" << matchPercent 
				// << ",alnScore:" << FMWResult.alnScore << "\n";
				// << ",kmerFreq:" << FMWResult.kmerFreq << "\n";
			// std::cout << FMWResult.mergedSeq << "\n";
			// printf("\n%s\n%s\n%s\n", aln_global->out1, aln_global->outm, aln_global->out2);
		// }
		
		// bug fix: m_targetSeed may be shorter than m_minOverlap
		if(m_strTargetSeed.length() > m_minOverlap)
			FMWResult.mergedSeq = FMWResult.mergedSeq + m_strTargetSeed.substr(m_minOverlap);
		
		FMWResult.typeFMWalkResult = 1;
		return;
	}
	
	FMWResult.typeFMWalkResult = -4;
		return;
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
			if((*iter)->fwdInterval.isValid())
				(*iter)->addKmerCount( (*iter)->fwdInterval.size());
			if((*iter)->rvcInterval.isValid())
				(*iter)->addKmerCount( (*iter)->rvcInterval.size());
			
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
				if(pChildNode->fwdInterval.isValid())
					pChildNode->addKmerCount( pChildNode->fwdInterval.size());
				if(pChildNode->rvcInterval.isValid())
					pChildNode->addKmerCount( pChildNode->rvcInterval.size());

                newLeaves.push_back(pChildNode);
            }
        }
    }
}

void SAIntervalPBHybridCTree::extendLeaves_v2()
{
	STNodePtrList newLeaves;
	
	// don't be too greedy.
	if(m_currentKmerSize > m_maxOverlap && m_leaves.size() < m_MaxLeaves*2/7)
	{
		refineSAInterval(m_maxOverlap);
		m_currentKmerSize=m_maxOverlap;
	}

	// the size of leaves is close to maxLeaves, so increase the FMW Threshold.
	if(m_leaves.size() > m_MaxLeaves/2)
		m_minSAThreshold = m_iniMinSAThreshold*3;
	else
		m_minSAThreshold = m_iniMinSAThreshold;
	
	// attempt to extend one base for each leave
	attempToExtend(newLeaves);
	
	// reduce current k-mer size until leaves appear
	while(newLeaves.empty() && m_currentKmerSize > m_minOverlap)
	{
		m_currentKmerSize--;
		refineSAInterval(m_currentKmerSize);
		attempToExtend(newLeaves);
	}

	// extension succeed
	if(!newLeaves.empty())
	{
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
    size_t IntervalSizeCutoff=m_minSAThreshold;    //min freq at fwd and rvc bwt, >=3 is equal to >=2 kmer freq

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
