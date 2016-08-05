///----------------------------------------------
// Copyright 2016 National Chung Cheng University
// Written by Yao-Ting Huang & Ping-Yeh Chen
// Released under the GPL
//-----------------------------------------------
//
// ShortReadOverlapTree - A fast overlap filter using locality-sensitive backward search.
//
//
#include "ShortReadOverlapTree.h"
#include "BWTAlgorithms.h"
#include "stdaln.h"

//
// Class: SAIOverlapTree
ShortReadOverlapTree::ShortReadOverlapTree(const std::string& sourceSeed,
				const std::string& strBetweenSrcTarget,
				const std::string& targetSeed,				
				int disBetweenSrcTarget,
				size_t minOverlap,
				size_t maxOverlap,
				const BWT* pBWT, 
				const BWT* pRBWT, 
				size_t min_SA_threshold, 				
				size_t maxIndelSize,
				double errorRate,
				size_t maxLeaves,
				size_t seedSize, 
				size_t repeatFreq):
				m_sourceSeed(sourceSeed), 
				m_strBetweenSrcTarget(strBetweenSrcTarget),
				m_targetSeed(targetSeed),				
				m_disBetweenSrcTarget(disBetweenSrcTarget),
				m_minOverlap(minOverlap), 
				m_maxOverlap(maxOverlap), 
				m_maxIndelSize(maxIndelSize), 
				m_pBWT(pBWT), 
				m_pRBWT(pRBWT),
				m_min_SA_threshold(min_SA_threshold),
				m_errorRate(errorRate),
				m_maxLeaves(maxLeaves), 
				m_seedSize(seedSize), 
				m_repeatFreq(repeatFreq)
{	
	std::string beginningkmer = m_sourceSeed.substr(m_sourceSeed.length()-m_minOverlap);
	// create one root node
	m_pRootNode = new SAIOverlapNode2(&m_sourceSeed, NULL);
	// store initial str of root
	m_pRootNode->computeInitial(m_sourceSeed);
	m_pRootNode->fwdInterval = BWTAlgorithms::findInterval(m_pRBWT, reverse(beginningkmer));
	m_pRootNode->rvcInterval = BWTAlgorithms::findInterval(m_pBWT, reverseComplement(beginningkmer));

	m_pRootNode->lastOverlapLen = m_pRootNode->currOverlapLen 
								= m_pRootNode->queryOverlapLen 
								= m_minOverlap;
	m_pRootNode->lastSeedIdx = m_pRootNode->initSeedIdx 
							= m_minOverlap - m_seedSize;
	m_pRootNode->totalSeeds = m_minOverlap - m_seedSize + 1;
	// push new node into roots and leaves vector
	m_RootNodes.push_back(m_pRootNode);
	m_leaves.push_back(m_pRootNode);
	
	// initialize the ending SA intervals with kmer length = m_minOverlap
	std::string endingkmer = m_targetSeed.substr(0, m_minOverlap);

	// PacBio reads are longer than real length due to insertions
	m_maxLength = (1.1*(m_disBetweenSrcTarget+10))+2*m_minOverlap;
	m_minLength = (0.9*(m_disBetweenSrcTarget-30))+2*m_minOverlap;
	//std::cout << m_MaxLength << ":\t" << m_MinLength;

	m_fwdTerminatedInterval = BWTAlgorithms::findInterval(m_pRBWT, reverse(endingkmer));
	m_rvcTerminatedInterval = BWTAlgorithms::findInterval(m_pBWT, reverseComplement(endingkmer));

	m_currentLength = m_currentKmerSize
					= m_minOverlap;

	m_query = beginningkmer + m_strBetweenSrcTarget + endingkmer;
	// std::cout << m_query << "\n";
	// std::cout << m_query << "\n" << beginningkmer << "\n" << endingkmer << "\n\n";
	
	// put SA intervals into m_fwdIntervals and m_rvcIntervals cache
	m_fwdIntervals.reserve(m_query.length()-m_seedSize+1);
	for(int i = 0; i <= (int)m_query.length()-(int)m_seedSize ; i++)
	{
		std::string seedStr = m_query.substr(i, m_seedSize);
		BWTInterval bi = BWTAlgorithms::findInterval( m_pRBWT, reverse(seedStr) );
		if(bi.isValid())
			m_fwdIntervals.push_back( TreeInterval<size_t>(bi.lower, bi.upper, i) );
	}
	fwdIntervalTree = IntervalTree<size_t>(m_fwdIntervals);
	m_rvcIntervals.reserve(m_query.length()-m_seedSize+1);
	for(int i = 0; i <= (int)m_query.length()-(int)m_seedSize ; i++)
	{
		std::string seedStr = m_query.substr(i, m_seedSize);
		BWTInterval bi = BWTAlgorithms::findInterval( m_pBWT, reverseComplement(seedStr) );
		if(bi.isValid())
			m_rvcIntervals.push_back( TreeInterval<size_t>(bi.lower, bi.upper, i) );
	}
	rvcIntervalTree = IntervalTree<size_t>(m_rvcIntervals);
}

//
ShortReadOverlapTree::~ShortReadOverlapTree()
{
	for (std::list<SAIOverlapNode2*>::iterator it = m_RootNodes.begin(); it != m_RootNodes.end(); ++it)
		delete *it;
	
	m_RootNodes.clear();
    // Recursively destroy the tree
    // delete m_pRootNode;
}

// comparison, not case sensitive.
static bool SeedComparator(const SAIOverlapNode2* first, const SAIOverlapNode2* second)
{
  return ( 	first->totalSeeds > second->totalSeeds );
}

//On success return the length of merged string
int ShortReadOverlapTree::extendOverlap(FMWalkResult &FMWResult)
{
	SAIntervalNodeResultVector results;
	
	//Overlap extension via FM-index walk
    while(!m_leaves.empty() && m_leaves.size() <= m_maxLeaves && m_currentLength <= m_maxLength)
    {
		// ACGT-extend the leaf nodes via updating existing SA interval
        extendLeaves();
					
		// Remove leaves without seed support within m_maxIndelSize
		PrunedBySeedSupport();
		
		// speedup by retaining the top bestN candidates after sufficient overlap length
		// This is the 3rd filter less reliable than previous ones
		const size_t bestN = 50;
		if(m_leaves.size()>= bestN)
		{
			m_leaves.sort(SeedComparator);
			SONode2PtrList::iterator iter1 = m_leaves.begin();
			SONode2PtrList::iterator iter2 = m_leaves.end();
			advance(iter1, bestN-1);
			m_leaves.erase(iter1, iter2);
		}
	
		//see if any interval reach left end $ with sufficient overlap
		if(m_currentLength >= m_minLength)
			isTerminated(results);
			
	}

	// reach the terminal kmer
	if(results.size() > 0)
	{
		// find the path with maximum match percent or kmer coverage
		return findTheBestPath(results, FMWResult);
	}
	
	// Did not reach the terminal kmer
    if(m_leaves.empty())
        return -1;	//high error
    else if(m_currentLength > m_maxLength)
        return -2;	//exceed search depth
    else if(m_leaves.size() > m_maxLeaves)
        return -3;	//too much repeats
	else
		return -4;
}

int ShortReadOverlapTree::findTheBestPath(SAIntervalNodeResultVector results, FMWalkResult &FMWResult)
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
		AlnAln *aln_global;
		aln_global = aln_stdaln(m_query.c_str(), candidateSeq.c_str(), &aln_param_pacbio, 1, 1);
		
		/*
		if(m_debugMode)
		{
			std::cout << ">pathBetweenSrcTarget:" << i+1 << ",len:" << pathBetweenSrcTarget.length() 
				<<  ",identity:" << matchPercent 
				<< ",aln score:" << aln_global->score << "\n";
			std::cout << pathBetweenSrcTarget << "\n";
			printf("\n%s\n%s\n%s\n", aln_global->out1, aln_global->outm, aln_global->out2);
		}
		*/

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
void ShortReadOverlapTree::printAll()
{
	for (std::list<SAIOverlapNode2*>::iterator it = m_RootNodes.begin(); it != m_RootNodes.end(); ++it)
		(*it)->printAllStrings("");

}

void ShortReadOverlapTree::extendLeaves()
{
    SONode2PtrList newLeaves;
	
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
		// if(m_beginningIntervalSize >= m_coverage*0.8 || m_terminatedIntervalSize >= m_coverage*0.8) // 256: 16, extension may exceed max leaves soon after
		//if(m_beginningIntervalSize >= 80 || m_terminatedIntervalSize >= 80)
			// refineSAInterval(81);
		// else
			refineSAInterval(m_minOverlap);
	}
	
    //shrink the SAIntervals in case overlap is larger than read length
    // if(!m_kmerMode  &&  newLeaves.empty() )
	if(newLeaves.empty() )
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

// Refine SA intervals of each leave with a new kmer
void ShortReadOverlapTree::refineSAInterval(size_t newKmer)
{
    for(SONode2PtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
    {
        // reset the SA intervals using original m_minOverlap
        std::string pkmer = (*iter)->getSuffix(newKmer);
        (*iter)->fwdInterval=BWTAlgorithms::findInterval(m_pRBWT, reverse(pkmer));
        (*iter)->rvcInterval=BWTAlgorithms::findInterval(m_pBWT, reverseComplement(pkmer));

    }
	m_currentKmerSize=newKmer;
}

void ShortReadOverlapTree::attempToExtend(SONode2PtrList &newLeaves)
{
    for(SONode2PtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
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
			
			// currOverlapLen/queryOverlapLen always increase wrt each extension
			// in order to know the approximate real-time matched length for terminal/containment processing
            (*iter)->currOverlapLen++;
			(*iter)->queryOverlapLen++;
			newLeaves.push_back(*iter);
        }
        else if(extensions.size() > 1)
        {
            // Branch
            for(size_t i = 0; i < extensions.size(); ++i)
            {
                SAIOverlapNode2* pChildNode = (*iter)->createChild(extensions[i].first);
                pChildNode->fwdInterval=extensions[i].second.interval[0];
                pChildNode->rvcInterval=extensions[i].second.interval[1];
				//inherit accumulated kmerCount from parent
				pChildNode->addKmerCount( (*iter)->getKmerCount() );
				if(pChildNode->fwdInterval.isValid()) pChildNode->addKmerCount( pChildNode->fwdInterval.size());
				if(pChildNode->rvcInterval.isValid()) pChildNode->addKmerCount( pChildNode->rvcInterval.size());

				pChildNode->currOverlapLen++;
				pChildNode->queryOverlapLen++;
                
				newLeaves.push_back(pChildNode);
            }
        }
    }
}

bool ShortReadOverlapTree::PrunedBySeedSupport()
{
	// the seed index in m_TerminatedIntervals for m_currentLength
	// the m_currentLength is the same for all leaves
	// which is used as the central index within the m_maxIndelSize window
	size_t currSeedIdx = m_currentLength-m_seedSize;
	
	//        ---	seed size =3. seed dist = 1;
	//         --*
	//          -*-
	//           *--
	//            ---
	// *: SNP or indel errors, ---: seed size
	// Erase the leaf if no feasible seeds are found within seedSize+maxIndelSize.
	size_t indelOffset = m_seedSize+m_maxIndelSize;
			
	// Compute the range of small and large indices for tolerating m_maxIndelSize
	size_t smallSeedIdx = currSeedIdx <= indelOffset ? 0 : currSeedIdx - indelOffset;
	size_t largeSeedIdx = (currSeedIdx+indelOffset) >= (m_query.length()-m_seedSize)?
						  (m_query.length()-m_seedSize):currSeedIdx+indelOffset;

	// check range of last seed and find new seeds for each interval
	SONode2PtrList::iterator iter = m_leaves.begin(); 
    while(iter != m_leaves.end())
    {
		double currErrorRate = computeErrorRate(*iter);

		// speedup by skipping dissimilar reads
		// This is the 2nd filter less reliable than the 1st one 
		if(currErrorRate > m_errorRate+0.04)
		{
			iter = m_leaves.erase(iter);
			continue;
		}

		if( m_currentLength - (*iter)->lastOverlapLen >= m_seedSize ||
			m_currentLength - (*iter)->lastOverlapLen <= 1)
		{
			// search for matched new seeds
			bool isNewSeedFound = isSupportedByNewSeed(*iter, smallSeedIdx, largeSeedIdx);

			// lastSeedIdxOffset records the offset between lastSeedIdx and currSeedIdx when first match is found
			if(isNewSeedFound)
				(*iter)->lastSeedIdxOffset = (int)(*iter)->lastSeedIdx - (int)currSeedIdx;
			
			// If the seed extension is stopped by SNP or indel error for the 1st time
			// increment the error number in order to distinguish two separate seeds
			// and one larger consecutive seed during error rate computation
			if(!isNewSeedFound && currSeedIdx+(*iter)->lastSeedIdxOffset == (*iter)->lastSeedIdx+1)
				(*iter)->numOfErrors++;
		}
		
		iter++;
    }
	
    return true;
}

// Identify new seeds wrt currSeedIdx
bool ShortReadOverlapTree::isSupportedByNewSeed(SAIOverlapNode2* currNode, size_t smallSeedIdx, size_t largeSeedIdx)
{
	// If there is mismatch/indel, jump to the next m_seedSize/m_seedDist, and 1 otherwise.
	size_t seedIdxOffset = currNode->lastOverlapLen < m_currentLength-m_seedSize?
							m_seedSize:m_currentLength - currNode->lastOverlapLen;

	// search for new seed starting from last matched seed or smallSeedIdx
	size_t startSeedIdx = std::max(smallSeedIdx, currNode->lastSeedIdx+seedIdxOffset);
	
	bool isNewSeedFound = false;
	BWTInterval currFwdInterval = currNode->fwdInterval;
	BWTInterval currRvcInterval = currNode->rvcInterval;

	// Binary search for new seeds using Query interval tree
	std::vector<TreeInterval<size_t> > resultsFwd, resultsRvc;
	if(currFwdInterval.isValid())
		fwdIntervalTree.findOverlapping(currFwdInterval.lower, currFwdInterval.upper, resultsFwd);
	if(currRvcInterval.isValid())
		rvcIntervalTree.findOverlapping(currRvcInterval.lower, currRvcInterval.upper, resultsRvc);
	int minIdxDiff = 10000;
	size_t currSeedIdx = m_currentLength-m_seedSize;
	for(size_t i=0 ; i<resultsFwd.size() || i<resultsRvc.size() ; i++)
	{
		if( currFwdInterval.isValid() && 
			i<resultsFwd.size() && 
			resultsFwd.at(i).value >= startSeedIdx && 
			resultsFwd.at(i).value <= largeSeedIdx )
		{
			// update currNode members
			if(std::abs(resultsFwd.at(i).value - currSeedIdx) < minIdxDiff)
			{					
				currNode->lastSeedIdx = resultsFwd.at(i).value;
				// query overlap may shift due to indels
				currNode->queryOverlapLen = resultsFwd.at(i).value+m_seedSize;
				minIdxDiff = std::abs(resultsFwd.at(i).value - currSeedIdx);
			}
			// lastOverlapLen records the overlap length of last hit
			currNode->lastOverlapLen = m_currentLength;
			// currOverlapLen is always identical to m_currentLength
			currNode->currOverlapLen = m_currentLength;
			isNewSeedFound = true;
		}
		else if( currRvcInterval.isValid() && 
			i<resultsRvc.size() && 
			resultsRvc.at(i).value >= startSeedIdx && 
			resultsRvc.at(i).value <= largeSeedIdx )
		{
			// update currNode members
			if(std::abs(resultsRvc.at(i).value - currSeedIdx) < minIdxDiff)
			{					
				currNode->lastSeedIdx = resultsRvc.at(i).value;
				// query overlap may shift due to indels
				currNode->queryOverlapLen = resultsRvc.at(i).value+m_seedSize;
				minIdxDiff = std::abs(resultsRvc.at(i).value - currSeedIdx);
			}
			// lastOverlapLen records the overlap length of last hit
			currNode->lastOverlapLen = m_currentLength;
			// currOverlapLen is always identical to m_currentLength
			currNode->currOverlapLen = m_currentLength;
			isNewSeedFound = true;
		}
	}
	
	if(isNewSeedFound)
		currNode->totalSeeds++;
	
	return isNewSeedFound;
}

double ShortReadOverlapTree::computeErrorRate(const SAIOverlapNode2* currNode)
{	
	// Compute accuracy via matched length in both query and subject
	int matchedLen = (int)currNode->totalSeeds*2;
	
	// SNP and indel over-estimate the unmatched lengths across error, ---*---
	// Restore the unmatched region via numOfErrors, which is still over-estimated
	matchedLen += (int)(currNode->numOfErrors*(m_seedSize-1)*2) ;
	
	int totalLen = (int)currNode->queryOverlapLen + currNode->currOverlapLen- (m_seedSize*2) +2;

	int unmatchedLen = totalLen - matchedLen;

	 //std::cout << currNode->queryOverlapLen  << "\t" << currNode->currOverlapLen
	//	  << "\t" << currNode->totalSeeds << "\n";

	 // std::cout << (double)unmatchedLen/totalLen << "\t" << matchedLen << "\t" << unmatchedLen 
	//		 << "\t" << totalLen << "\t"<< currNode->numOfErrors <<"\n";
	
	return unmatchedLen/(double)totalLen;
}
			
//update SA intervals of each leaf, which corresponds to one-base extension
std::vector<std::pair<std::string, BWTIntervalPair> > ShortReadOverlapTree::getFMIndexExtensions(SAIOverlapNode2* pNode)
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

// Check for leaves whose extension has terminated. If the leaf has
// terminated, the walked string and coverage is pushed to the result vector
bool ShortReadOverlapTree::isTerminated(SAIntervalNodeResultVector& results)
{
	bool found = false;

    for(SONode2PtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
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