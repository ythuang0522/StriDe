///----------------------------------------------
// Copyright 2016 National Chung Cheng University
// Written by Yao-Ting Huang & Ping-Yeh Chen
// Released under the GPL
//-----------------------------------------------
//
// LongReadSelfCorrectByOverlap - A fast overlap filter using locality-sensitive backward search.
//
//
#include "LongReadCorrectByOverlap.h"
#include "BWTAlgorithms.h"
#include "stdaln.h"
#include "LongReadOverlap.h"

//
// Class: SAIOverlapTree
LongReadSelfCorrectByOverlap::LongReadSelfCorrectByOverlap(
                const std::string& sourceSeed,
				const std::string& strBetweenSrcTarget,
				const std::string& targetSeed,				
				int disBetweenSrcTarget,
                size_t initkmersize,
                size_t maxOverlap,
                const FMextendParameters params,
                size_t min_SA_threshold, 				
				double errorRate, 
				size_t repeatFreq,
                size_t localSimilarlykmerSize
                ):
				m_sourceSeed(sourceSeed), 
				m_strBetweenSrcTarget(strBetweenSrcTarget),
				m_targetSeed(targetSeed),				
				m_disBetweenSrcTarget(disBetweenSrcTarget),
                m_initkmersize(initkmersize),
				m_minOverlap(params.minKmerLength),
                m_maxOverlap(maxOverlap),
                m_BWTindices(params.indices),
				m_pBWT(params.indices.pBWT), 
				m_pRBWT(params.indices.pRBWT),
                m_PBcoverage(params.PBcoverage),
				m_min_SA_threshold(min_SA_threshold),
				m_errorRate(errorRate),
				m_maxLeaves(params.maxLeaves), 
				m_seedSize(params.idmerLength), 
				m_repeatFreq(repeatFreq),
                m_localSimilarlykmerSize(localSimilarlykmerSize),
                m_PacBioErrorRate(params.ErrorRate),
                m_isDebug(params.Debug)
{	
	std::string beginningkmer = m_sourceSeed.substr(m_sourceSeed.length()-m_initkmersize);
    
    //if distance < 100 ,use const indel size 
    m_maxIndelSize  = m_disBetweenSrcTarget > 100 ? m_disBetweenSrcTarget * 0.2 : 20; 
    
    
    //initialRootNode
    initialRootNode(beginningkmer);
	
	// push new node into roots and leaves vector
	m_RootNodes.push_back(m_pRootNode);
	m_leaves.push_back(m_pRootNode);

    
    //frequencies of correspond k
    for(double i = m_minOverlap ; i <= m_maxOverlap+1 ; i++)
        freqsOfKmerSize.push_back(pow(1-m_PacBioErrorRate,i)*m_PBcoverage);

    
    
    if(m_isDebug)
	std::cout << "BE: " << beginningkmer << " " << m_targetSeed << " "<< m_pRootNode->fwdInterval.size() + m_pRootNode->rvcInterval.size() << "|" << disBetweenSrcTarget <<"\n"; //debugch
	
    // PacBio reads are longer than real length due to insertions
	m_maxLength = (1.2*(m_disBetweenSrcTarget+10))+2*m_initkmersize;
	m_minLength = (0.8*(m_disBetweenSrcTarget-20))+2*m_initkmersize;
    
  
    
    // initialize the ending SA intervals with kmer length = m_minOverlap
    for(int i =0 ;i <= m_targetSeed.length()-m_minOverlap; i++)
    {
        std::string endingkmer = m_targetSeed.substr(i, m_minOverlap);
        m_fwdTerminatedInterval.push_back(BWTAlgorithms::findInterval(m_pRBWT, reverse(endingkmer)));
        m_rvcTerminatedInterval.push_back(BWTAlgorithms::findInterval(m_pBWT, reverseComplement(endingkmer)));
    }
    //build overlap tree
    buildOverlapbyFMindex(beginningkmer);

}

 void LongReadSelfCorrectByOverlap::initialRootNode(std::string beginningkmer)
{	
    // create one root node
	m_pRootNode = new SAIOverlapNode3(&m_sourceSeed, NULL);
    
    // store initial str of root
	m_pRootNode->computeInitial(beginningkmer);
	m_pRootNode->fwdInterval = BWTAlgorithms::findInterval(m_pRBWT, reverse(beginningkmer));
	m_pRootNode->rvcInterval = BWTAlgorithms::findInterval(m_pBWT, reverseComplement(beginningkmer));
	m_pRootNode->lastOverlapLen = m_currentLength = m_pRootNode->currOverlapLen = m_pRootNode->queryOverlapLen = m_currentKmerSize = m_initkmersize;
	m_pRootNode->lastSeedIdx = m_pRootNode->initSeedIdx = m_initkmersize - m_seedSize;
	m_pRootNode->totalSeeds = m_initkmersize - m_seedSize + 1;
	m_pRootNode->numRedeemSeed = 0;
    m_pRootNode->LocalErrorRateRecord.push_back(0);
    m_pRootNode->GlobalErrorRateRecord.push_back(0);
    
}

void LongReadSelfCorrectByOverlap::buildOverlapbyFMindex(std::string beginningkmer)
{
   	m_query = beginningkmer + m_strBetweenSrcTarget + m_targetSeed;
	
	// put SA intervals into m_fwdIntervals and m_rvcIntervals cache 
	m_fwdIntervals.reserve(m_query.length()-m_seedSize+1);
	m_rvcIntervals.reserve(m_query.length()-m_seedSize+1);

    for(int i = 0; i <= (int)m_query.length()-(int)m_seedSize ; i++)//build overlap tree
	{
         
		std::string seedStr = m_query.substr(i, m_seedSize);
		BWTInterval bi;
		bi = BWTAlgorithms::findInterval( m_pRBWT, reverse(seedStr) );        
		if(bi.isValid())
			m_fwdIntervals.push_back( TreeInterval<size_t>(bi.lower, bi.upper, i) );
		bi = BWTAlgorithms::findInterval( m_pBWT, reverseComplement(seedStr) );       
		if(bi.isValid())
			m_rvcIntervals.push_back( TreeInterval<size_t>(bi.lower, bi.upper, i) );
	}    
	fwdIntervalTree = IntervalTree<size_t>(m_fwdIntervals);
	rvcIntervalTree = IntervalTree<size_t>(m_rvcIntervals); 
    
    
    
}
//
LongReadSelfCorrectByOverlap::~LongReadSelfCorrectByOverlap()
{
    

	for (std::list<SAIOverlapNode3*>::iterator it = m_RootNodes.begin(); it != m_RootNodes.end(); ++it)
		delete *it;
	
	m_RootNodes.clear();
    // delete hashIndex;
}

// comparison, not case sensitive.
static bool SeedComparator(const SAIOverlapNode3* first, const SAIOverlapNode3* second)
{
  return ( 	first->totalSeeds > second->totalSeeds );
}

//On success return the length of merged string
int LongReadSelfCorrectByOverlap::extendOverlap(FMWalkResult2 &FMWResult)
{
	SAIntervalNodeResultVector results;
	
	//Overlap extension via FM-index walk
    while(!m_leaves.empty() && m_leaves.size() <= m_maxLeaves && m_currentLength <= m_maxLength)
    {
        
		// ACGT-extend the leaf nodes via updating existing SA interval
        SONode3PtrList newLeaves = extendLeaves();

		
        //use overlap tree to trim branch
		PrunedBySeedSupport(newLeaves);
       
        //update leaves
        m_leaves.clear();
        m_leaves = newLeaves;

		
		// printleaves();
        if(m_isDebug)
        std::cout << "----" << std::endl;
        

		if(m_currentLength >= m_minLength)
            isTerminated(results);
        
        
	}

	// reach the terminal kmer
	if(results.size() > 0)
		// find the path with maximum match percent or kmer coverage
		return findTheBestPath(results, FMWResult);

	// Did not reach the terminal kmer
    if(m_leaves.empty())	//high error   
        return -1;
    else if(m_currentLength > m_maxLength)	//exceed search depth
        return -2;
    else if(m_leaves.size() > m_maxLeaves)	//too much repeats
        return -3;
	else
		return -4;
}


int LongReadSelfCorrectByOverlap::findTheBestPath(SAIntervalNodeResultVector results, FMWalkResult2 &FMWResult)
{
	
    double minErrorRate = 1;
    // for (size_t i = 0 ; i < results.size() ;i++)
       // minErrorRate = minErrorRate < results[i].errorRate ? minErrorRate : results[i].errorRate;
    
	for (size_t i = 0 ; i < results.size() ;i++)
	{
        std::string candidateSeq;
        if(m_targetSeed.length() > m_minOverlap)
            candidateSeq = results[i].thread ;
		else
			candidateSeq = results[i].thread;
      
        if(results[i].errorRate < minErrorRate )
            FMWResult.mergedSeq = candidateSeq;
            
        
	}
        // std::cout << FMWResult.mergedSeq << "merge seq\n"; //testcw

	if(FMWResult.mergedSeq.length() != 0)
		return 1;
	return -4;
}


SONode3PtrList LongReadSelfCorrectByOverlap::extendLeaves()
{
    SONode3PtrList newLeaves;

    if(m_currentKmerSize > m_maxOverlap) // resize when size up to upper bound
        refineSAInterval(m_maxOverlap);
    // std::cout<<   m_currentKmerSize   <<"mer\n";  
    attempToExtend(newLeaves);
    
    
    if(newLeaves.empty() ) //level 1 reduce size 
    { 
        size_t LowerBound = m_currentKmerSize - 2 >= m_minOverlap ? m_currentKmerSize - 2 : m_minOverlap;
        size_t ReduceSize = SelectFreqsOfrange(LowerBound,m_currentKmerSize,m_leaves);
        refineSAInterval(ReduceSize);
         // std::cout<<   m_currentKmerSize << " " << ReduceSize  <<"mer\n";  
        attempToExtend(newLeaves);
        
      
        if( newLeaves.empty() )//level 2 reduce threshold
        {
            
            m_min_SA_threshold--;
            attempToExtend(newLeaves);
            m_min_SA_threshold++;
        }

    }
    
    //extension succeed
	if(!newLeaves.empty())
    { 
        m_currentLength++;  
		m_currentKmerSize++;
        if( isInsufficientFreqs(newLeaves)   )// if frequency are low , relax it
        {
            
            // size_t ReduceSize = m_currentKmerSize > m_minOverlap-1 ? m_currentKmerSize-1: m_currentKmerSize;
            

            size_t LowerBound = m_currentKmerSize > m_minOverlap + 1 ? m_currentKmerSize - 2 : m_minOverlap;
            size_t ReduceSize = SelectFreqsOfrange(LowerBound,m_currentKmerSize,newLeaves);
            
            for(SONode3PtrList::iterator iter = newLeaves.begin(); iter != newLeaves.end(); ++iter)
            {   
                
                // reset the SA intervals 
                std::string pkmer = (*iter)->getSuffix(ReduceSize);
                (*iter)->fwdInterval=BWTAlgorithms::findInterval(m_pRBWT, reverse(pkmer));
                (*iter)->rvcInterval=BWTAlgorithms::findInterval(m_pBWT, reverseComplement(pkmer));

            }
            
            m_currentKmerSize = ReduceSize;

           
        }     

    }

    
       return newLeaves;
    
}
size_t LongReadSelfCorrectByOverlap::SelectFreqsOfrange(size_t LowerBound,size_t UpperBound,SONode3PtrList &newLeaves)
{
   std::vector<std::pair<std::string,BWTIntervalPair>> pkmers;
   size_t tempmaxfmfreqs = 0;
// std::cout<<      "here\n";
   for(SONode3PtrList::iterator iter = newLeaves.begin(); iter != newLeaves.end() ; ++iter)
    {
        std::string pkmer = (*iter)->getFullString().substr((*iter)->getFullString().length() -UpperBound);

        std::string startkmer = pkmer.substr(UpperBound - LowerBound); //  string of lower bound kmer size
        
        BWTInterval Fwdinterval = BWTAlgorithms::findInterval(m_pBWT, startkmer);
        BWTInterval Rvcinterval = BWTAlgorithms::findInterval(m_pRBWT, reverseComplement(reverse(startkmer)));
        BWTIntervalPair bip;
        bip.interval[0] = Fwdinterval;
        bip.interval[1] = Rvcinterval;
        

        pkmers.push_back(std::make_pair(pkmer,bip));
        // std::cout<< startkmer << "  "<< Fwdinterval.size() +  Rvcinterval.size()    <<" haha \n";
        if( Fwdinterval.size() +  Rvcinterval.size() > (int)tempmaxfmfreqs) //check interval size
            tempmaxfmfreqs = Fwdinterval.size() +  Rvcinterval.size();

    }
    // std::cout<<   (int)tempmaxfmfreqs - (int)freqsOfKmerSize.at(LowerBound - m_minOverlap)   <<"\n";
    if( (int)tempmaxfmfreqs - (int)freqsOfKmerSize.at(LowerBound - m_minOverlap) < 5 ) return LowerBound;
    
    for(size_t i=1 ; i <= UpperBound - LowerBound; i++ )
    {
       tempmaxfmfreqs = 0;
       for(size_t j = 0; j < pkmers.size(); j++)
       {
        std::string startkmer = pkmers.at(j).first.substr(UpperBound - LowerBound - i);
        BWTInterval Fwdinterval = pkmers.at(j).second.interval[0];
        BWTInterval Rvcinterval = pkmers.at(j).second.interval[1];
        
        char b = startkmer[0];
        char rcb;
        if(b == 'A' ) rcb = 'T' ;
        else if(b == 'C' ) rcb = 'G' ;
        else if(b == 'G' ) rcb = 'C' ;
        else if(b == 'T' ) rcb = 'A' ;
        BWTAlgorithms::updateInterval(Fwdinterval,b,m_pBWT);
        BWTAlgorithms::updateInterval(Rvcinterval,rcb,m_pRBWT);
        pkmers.at(j).second.interval[0] = Fwdinterval;
        pkmers.at(j).second.interval[1] = Rvcinterval;

        
        if( Fwdinterval.size() +  Rvcinterval.size()> (int)tempmaxfmfreqs)
            tempmaxfmfreqs = Fwdinterval.size() +  Rvcinterval.size();

        
       }

       
       if( (int)tempmaxfmfreqs - (int)freqsOfKmerSize.at(LowerBound - m_minOverlap + i) < 5 ) return LowerBound +i ;
       
       
       

    }
    
    
    return UpperBound;

    
}





bool LongReadSelfCorrectByOverlap::isInsufficientFreqs(SONode3PtrList &newLeaves)
{
    size_t highfreqscount = 0;
    for(SONode3PtrList::iterator iter = newLeaves.begin(); iter != newLeaves.end(); ++iter)
    {
        
       if( (*iter)->fwdInterval.size()+(*iter)->rvcInterval.size() > 3)
           highfreqscount++;

    }
    
    if(highfreqscount == 0)
        return true;
    else if (highfreqscount <= 2 && newLeaves.size()>=5 )
        return true;
    else if (highfreqscount <= 1 && newLeaves.size()>=3 )
        return true;
    return false;
}

// Refine SA intervals of each leave with a new kmer
void LongReadSelfCorrectByOverlap::refineSAInterval(size_t newKmer)
{   
       
    for(SONode3PtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
    {
        
        // reset the SA intervals using original m_minOverlap
        std::string pkmer = (*iter)->getSuffix(newKmer);
        
        (*iter)->fwdInterval=BWTAlgorithms::findInterval(m_pRBWT, reverse(pkmer));
        (*iter)->rvcInterval=BWTAlgorithms::findInterval(m_pBWT, reverseComplement(pkmer));

    }
    
	m_currentKmerSize=newKmer;
}




void LongReadSelfCorrectByOverlap::attempToExtend(SONode3PtrList &newLeaves)
{
    
    double minimumErrorRate = 1;

    std::vector<size_t> frequencies;
   
    for(SONode3PtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
    {
        if((*iter)->LocalErrorRateRecord.back() < minimumErrorRate)
            minimumErrorRate = (*iter)->LocalErrorRateRecord.back();

        
    }
    for(SONode3PtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
    {

        if((double)(*iter)->LocalErrorRateRecord.back() - (double)minimumErrorRate > 0.05 && m_currentLength > m_localSimilarlykmerSize/2 )
        {
            
           iter = m_leaves.erase(iter);
           continue;
        } 
        
        if((double)(*iter)->LocalErrorRateRecord.back() - (double)minimumErrorRate > 0.1 && m_currentLength > 15 )
        {
            
           iter = m_leaves.erase(iter);
           continue;
        } 
        
    }

    for(SONode3PtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
    {
        
        std::vector< std::pair<std::string, BWTIntervalPair> > extensions;
        char lastword = (*iter)->getFullString().back();
           
        extensions = getFMIndexExtensions(*iter);
        // std::cout<<extensions.size()<<"  size\n";
        // Either extend the current node or branch it
        // If no extension, do nothing and this node
        // is no longer considered a leaf
        
        if (extensions.size()>0) updateLeaves(newLeaves,extensions,*iter);

        else  // no extensions , give one more chance for lowest error rate path
        {
            
            if((*iter)->LocalErrorRateRecord.back() == minimumErrorRate  && m_leaves.size()>1 )
            {   

                m_min_SA_threshold--;//reduce threshold
                extensions = getFMIndexExtensions(*iter);
                
                updateLeaves(newLeaves,extensions,*iter);

                m_min_SA_threshold++;
                
            }
        }
       

    }
}

void LongReadSelfCorrectByOverlap::updateLeaves(SONode3PtrList &newLeaves,std::vector< std::pair<std::string, BWTIntervalPair> > &extensions,SAIOverlapNode3* pNode)
{
    char lastword =  pNode->getFullString().back();
    if(extensions.size() == 1)
        {
 

            // Single extension, do not branch
            pNode->extend(extensions.front().first);
            pNode->fwdInterval = extensions.front().second.interval[0];
            pNode->rvcInterval = extensions.front().second.interval[1];
			if(pNode->fwdInterval.isValid())
				pNode->addKmerCount( pNode->fwdInterval.size());
			if(pNode->rvcInterval.isValid())
				pNode->addKmerCount( pNode->rvcInterval.size());
			
			// currOverlapLen/queryOverlapLen always increase wrt each extension
			// in order to know the approximate real-time matched length for terminal/containment processing
			pNode->currOverlapLen++;
			pNode->queryOverlapLen++;
			
			newLeaves.push_back(pNode);
        }
        
    else if(extensions.size() > 1)
    {
        
        // Branch
        for(size_t i = 0; i < extensions.size(); ++i)
        {   

            SAIOverlapNode3* pChildNode = pNode->createChild(extensions[i].first);
            

            
            
            pChildNode->fwdInterval=extensions[i].second.interval[0];
            pChildNode->rvcInterval=extensions[i].second.interval[1];

            //inherit accumulated kmerCount from parent
            pChildNode->addKmerCount( pNode->getKmerCount() );
            if(pChildNode->fwdInterval.isValid())
                pChildNode->addKmerCount( pChildNode->fwdInterval.size());
            if(pChildNode->rvcInterval.isValid())
                pChildNode->addKmerCount( pChildNode->rvcInterval.size());
            pChildNode->currOverlapLen++;
            pChildNode->queryOverlapLen++;
            
            newLeaves.push_back(pChildNode);
        }
    }
    
}
bool LongReadSelfCorrectByOverlap::PrunedBySeedSupport(SONode3PtrList &newLeaves)
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
	SONode3PtrList::iterator iter = newLeaves.begin(); 
    while(iter != newLeaves.end())
    {
        bool isNewSeedFound = false;
		if( m_currentLength - (*iter)->lastOverlapLen > m_seedSize ||
			m_currentLength - (*iter)->lastOverlapLen <= 1)
		{
            size_t preSeedIdx = (*iter)->lastSeedIdx;
			// search for matched new seeds
			isNewSeedFound = isSupportedByNewSeed(*iter, smallSeedIdx, largeSeedIdx);
           
			// lastSeedIdxOffset records the offset between lastSeedIdx and currSeedIdx when first match is found
			if(isNewSeedFound)
            {
                if( currSeedIdx+(*iter)->lastSeedIdxOffset - preSeedIdx > m_seedSize )
                   (*iter)->numRedeemSeed += (m_seedSize-1)*m_PacBioErrorRate;
                (*iter)->lastSeedIdxOffset = (int)(*iter)->lastSeedIdx - (int)currSeedIdx;
            }   
            
			// If the seed extension is stopped by SNP or indel error for the 1st time
			// increment the error number in order to distinguish two separate seeds
			// and one larger consecutive seed during error rate computation
			if(!isNewSeedFound && currSeedIdx+(*iter)->lastSeedIdxOffset == (*iter)->lastSeedIdx+1)
				(*iter)->numOfErrors ++;
            
            else if((currSeedIdx+(*iter)->lastSeedIdxOffset - (*iter)->lastSeedIdx) % m_seedSize == 1 )
                (*iter)->numOfErrors ++;
			else if(!isNewSeedFound && currSeedIdx+(*iter)->lastSeedIdxOffset - (*iter)->lastSeedIdx > m_seedSize-1)
				(*iter)->numRedeemSeed += 1-m_PacBioErrorRate;
            
            // else if(!isNewSeedFound && currSeedIdx+(*iter)->lastSeedIdxOffset - (*iter)->lastSeedIdx > m_seedSize+1)
				// (*iter)->numRedeemSeed += 0.5;
		}
        
		else
			(*iter)->numRedeemSeed += 1-m_PacBioErrorRate;
		
		
        // if(m_currentLength == 200)
        // {
              
            double currErrorRate = computeErrorRate(*iter);
          
           
  


		// speedup by skipping dissimilar reads
		// This is the 2nd filter less reliable than the 1st one 
		if(currErrorRate > m_errorRate  && (*iter)->LocalErrorRateRecord.size() > m_localSimilarlykmerSize ) //testcw
		{
			iter = newLeaves.erase(iter);
           
            
			continue;
		}
        else if (currErrorRate > m_errorRate && (*iter)->LocalErrorRateRecord.size() <= m_localSimilarlykmerSize)
        {
            iter = newLeaves.erase(iter);
           
          
			continue;
        }

		
		iter++;
    }
	
    return true;
}

// Identify new seeds wrt currSeedIdx
bool LongReadSelfCorrectByOverlap::isSupportedByNewSeed(SAIOverlapNode3* currNode, size_t smallSeedIdx, size_t largeSeedIdx)
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
			if(std::abs((int)resultsFwd.at(i).value - (int)currSeedIdx) < minIdxDiff)
			{
				currNode->lastSeedIdx = resultsFwd.at(i).value;

				// query overlap may shift due to indels
				currNode->queryOverlapLen = resultsFwd.at(i).value+m_seedSize;
				minIdxDiff = std::abs((int)resultsFwd.at(i).value - (int)currSeedIdx);
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
			if(std::abs( (int)currSeedIdx - (int)resultsRvc.at(i).value ) < minIdxDiff)
			{					
				currNode->lastSeedIdx = resultsRvc.at(i).value;

				// query overlap may shift due to indels
				currNode->queryOverlapLen = resultsRvc.at(i).value+m_seedSize;
				minIdxDiff = std::abs((int)resultsRvc.at(i).value - (int)currSeedIdx);
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

double LongReadSelfCorrectByOverlap::computeErrorRate(SAIOverlapNode3* currNode)
{	
	// Compute accuracy via matched length in both query and subject
	// double matchedLen = (double)currNode->totalSeeds*2;
    // double matchedLen = (double)currNode->totalSeeds ;
	

    
    double matchedLen = (double)currNode->totalSeeds + m_seedSize-1;
    
	// SNP and indel over-estimate the unmatched lengths across error, ---*---
	// Restore the unmatched region via numOfErrors, which is still over-estimated

	matchedLen += currNode->numRedeemSeed;

	double totalLen = (double)currNode->currOverlapLen ;
    

	double unmatchedLen = totalLen - matchedLen;


    
    double  currErrorRate =  unmatchedLen/totalLen;
    currNode->GlobalErrorRateRecord.push_back(currErrorRate);
   
    if(currNode->GlobalErrorRateRecord.size() >= m_localSimilarlykmerSize)
    {
        size_t totalsize = currNode->GlobalErrorRateRecord.size();
      
        currErrorRate = ( currErrorRate*totalLen-currNode->GlobalErrorRateRecord.at( totalsize - m_localSimilarlykmerSize)*(totalLen - m_localSimilarlykmerSize) )/m_localSimilarlykmerSize;
         
    }
	currNode->LocalErrorRateRecord.push_back(currErrorRate);
	return currErrorRate;
}






std::vector<std::pair<std::string, BWTIntervalPair> > LongReadSelfCorrectByOverlap::getFMIndexExtensions(SAIOverlapNode3* pNode)
{
    std::vector<std::pair<std::string, BWTIntervalPair> > out;
    size_t IntervalSizeCutoff = m_min_SA_threshold;    //min freq at fwd and rvc bwt, >=3 is equal to >=2 kmer freq
    
    std::vector<std::pair<size_t,BWTIntervalPair>> bvector;
     
    size_t totalcount = 0;
    size_t maxfreqsofleave = 0;
    std::string currKmer = pNode->getFullString().substr( pNode->getFullString().length()-m_currentKmerSize);

    
    
    if(m_isDebug)
    
    std::cout<<pNode->getFullString()<<" || "<< pNode->LocalErrorRateRecord.back()<<"\n" << currKmer << "\t" << pNode->fwdInterval.size()+ pNode->rvcInterval.size()<<"\n";
     

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
        BWTIntervalPair bip;
        bip.interval[0]=fwdProbe;
        bip.interval[1]=rvcProbe;
        size_t bcount = 0;
        
        if(fwdProbe.isValid())
            bcount += fwdProbe.size();
        if(rvcProbe.isValid())
            bcount += rvcProbe.size();
        
       // std::cout<<    (b == lastword) <<"||" << bcount <<"||" << WeightofContinuousChar <<"   hi\n"; 
       // if (b == lastword  && WeightofContinuousChar < 1)   
       // {    std::cout<<      "ytytyt\n";
           // bcount = int(bcount+WeightofContinuousChar)  > 0 ? bcount+WeightofContinuousChar : 0 ;
           // totalcount += WeightofContinuousChar;
        // }
        
        totalcount += bcount;
        // std::string bkmer = pNode->getFullString().substr( pNode->getFullString().length()- 14 )+b;
        if(m_isDebug)
        std::cout << "bcount:" << bcount << "||" << b << "\n";//testcw
        
        bvector.push_back(std::make_pair(bcount, bip));
        if(bcount >  maxfreqsofleave)
            maxfreqsofleave = bcount;
        

    }// end of ACGT
      // std::cout << "totalcount: " << totalcount << std::endl;//testcw
     
     
    if(totalcount > 1024)  // fillter low complex repeat e.g. AAAAAAAAAAAA
        return out;
     
     

       
        
       for(int i = 1; i < BWT_ALPHABET::size; ++i)
       {
           float bratio = (float)bvector.at(i-1).first/(float)totalcount;
           size_t bdiff = std::abs((int)bvector.at(i-1).first-(int)maxfreqsofleave);
          
           char b = BWT_ALPHABET::getChar(i);
          
        if (currKmer.substr(currKmer.length()-2,1) == currKmer.substr(currKmer.length()-1,1) && currKmer.substr(currKmer.length()-3,1) ==currKmer.substr(currKmer.length()-2,1))   
        {   
            bvector.at(i-1).first = (float)bvector.at(i-1).first / (float)maxfreqsofleave >= 0.6 ? bvector.at(i-1).first : 0;
           
         }
        

       if(bratio >= 0.1 &&(bvector.at(i-1).first >= IntervalSizeCutoff || ((float)bvector.at(i-1).first / (float)maxfreqsofleave >= 0.6 && totalcount >= IntervalSizeCutoff+2)))
        {
			
            // extend to b
            
            std::string tmp;
            tmp.append(1,b);
            BWTIntervalPair bip;
            bip.interval[0] = bvector.at(i-1).second.interval[0];
            bip.interval[1] = bvector.at(i-1).second.interval[1];
            out.push_back(std::make_pair(tmp, bip));
        }
       }

    return out;
}




// Check for leaves whose extension has terminated. If the leaf has
// terminated, the walked string and coverage is pushed to the result vector
bool LongReadSelfCorrectByOverlap::isTerminated(SAIntervalNodeResultVector& results)
{
	bool found = false;
    
    for(SONode3PtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
    {
        
        
        BWTInterval currfwd=(*iter)->fwdInterval;
        BWTInterval currrvc=(*iter)->rvcInterval;

        assert(currfwd.isValid() || currrvc.isValid());
		
		//The current SA interval stands for a string >= terminating kmer
		//If terminating kmer is a substr, the current SA interval is a sub-interval of the terminating interval
         bool isFwdTerminated = false;
         bool isRvcTerminated = false;
        for(int i = (*iter)->resultindex.second > 0 ? (*iter)->resultindex.second : 0;i <= m_targetSeed.length()-(int)m_minOverlap; i++)
       {
        isFwdTerminated=currfwd.isValid() && currfwd.lower >= m_fwdTerminatedInterval.at(i).lower
                            && currfwd.upper <= m_fwdTerminatedInterval.at(i).upper;
         isRvcTerminated=currrvc.isValid() && currrvc.lower >= m_rvcTerminatedInterval.at(i).lower
                            && currrvc.upper <= m_rvcTerminatedInterval.at(i).upper;
        
       if(isFwdTerminated || isRvcTerminated)
        {
           
            std::string STNodeStr = m_targetSeed.length() > m_minOverlap?(*iter)->getFullString() + m_targetSeed.substr(i+m_minOverlap):(*iter)->getFullString();
            
            SAIntervalNodeResult STresult;
            STresult.thread=STNodeStr;
			STresult.SAICoverage=(*iter)->getKmerCount();
            STresult.errorRate = (*iter)->GlobalErrorRateRecord.back();
            //compute the merged pos right next to the kmer on 2nd read.
            //STresult.index = m_minOverlap ;
            if((*iter)->resultindex.first == -1 )
            {
                results.push_back(STresult);    
                (*iter)->resultindex = std::make_pair(results.size(),i);
            }
            else 
            {   

                results.at((*iter)->resultindex.first-1) = STresult;
                (*iter)->resultindex = std::make_pair((*iter)->resultindex.first,i);
            }
            
            found =  true;
        }
       }
       

    }

    return found;
}