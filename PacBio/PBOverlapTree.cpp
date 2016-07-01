///----------------------------------------------
// Copyright 2016 National Chung Cheng University
// Written by Yao-Ting Huang
// Released under the GPL
//-----------------------------------------------
//
// PBOverlapTree - A fast overlap filter using locality-sensitive backward search.
//
//
#include "PBOverlapTree.h"
#include "BWTAlgorithms.h"

//
// Class: SAIOverlapTree
PBOverlapTree::PBOverlapTree(const std::string& query,
			   size_t srcKmerSize,
			   size_t minOverlap,
			   size_t maxIndelSize,
			   const BWT* pBWT, const BWT* pRBWT, 
			   double errorRate,
			   size_t maxLeaves,
			   size_t seedSize, size_t seedDist, 
			   size_t repeatFreq):
                               m_Query(query), m_minOverlap(minOverlap),  
                               m_maxIndelSize(maxIndelSize), m_pBWT(pBWT), m_pRBWT(pRBWT),
								m_errorRate(errorRate),
                               m_maxLeaves(maxLeaves), 
							   m_seedSize(seedSize), m_seedDist(seedDist), 
							   m_repeatFreq(repeatFreq)
{
	std::string initSeed = m_Query.substr(m_Query.length()-srcKmerSize, srcKmerSize);
	BWTIntervalPair bip = BWTAlgorithms::findIntervalPair( m_pBWT, m_pRBWT, initSeed);

	std::cout << initSeed << ":" << srcKmerSize << ":" << bip.interval[0].size() << "\n";
			
	// create one root node
	m_pRootNode = new SAIOverlapNode(&initSeed, NULL);
	m_pRootNode->currIntervalPair = bip;
	m_pRootNode->computeInitial(reverse(initSeed));   //store initial str in the root
	m_pRootNode->lastOverlapLen = m_pRootNode->currOverlapLen = m_pRootNode->queryOverlapLen 
								= m_currentLength = srcKmerSize;

	m_pRootNode->lastSeedIdx = m_pRootNode->initSeedIdx = (srcKmerSize - m_seedSize)/m_seedDist;
	m_pRootNode->totalSeeds = srcKmerSize - m_seedSize + 1;
	
	// push new node into roots and leaves vector
	m_RootNodes.push_back(m_pRootNode);
	m_leaves.push_back(m_pRootNode);

	// initialize the seeding SA intervals with seedSize and seedDist from rightmost to leftmost
	// Note that the seed index is reverse wrt string index of query
	// size_t uniqCount = 0;

	m_TerminatedIntervals.reserve(m_Query.length()-m_seedSize+1);
	for(int i = (int)m_Query.length()-(int)m_seedSize; i >= 0 ; i -= (int)m_seedDist)
	{
		std::string seedStr = m_Query.substr(i, m_seedSize);
		// reversely put SA intervals into m_TerminatedIntervals cache
		// m_TerminatedIntervals.push_back(BWTAlgorithms::findInterval( m_pBWT, seedStr) );
		BWTInterval bi = BWTAlgorithms::findInterval( m_pBWT, seedStr);
		m_TerminatedIntervals.push_back( TreeInterval<size_t>(bi.lower, bi.upper, m_Query.length()-(int)m_seedSize-i) );
		
		// if( /*m_RootNodes.front()->lastOverlapLen == 29 &&*/ bi.size() < 10 && bi.size() > 4)
		// {
			// uniqCount++;
			// std::cout << reverseComplement(seedStr) << ":" << bi.size() <<"\n";
		// }
	}
	// std::cout << uniqCount << "\n";
	
	BWTIntervalTree = IntervalTree<size_t>(m_TerminatedIntervals);
}

//
PBOverlapTree::~PBOverlapTree()
{
	for (std::list<SAIOverlapNode*>::iterator it = m_RootNodes.begin(); it != m_RootNodes.end(); ++it)
		delete *it;
	
	m_RootNodes.clear();
    // Recursively destroy the tree
    // delete m_pRootNode;
}

// comparison, not case sensitive.
static bool SeedComparator(const SAIOverlapNode* first, const SAIOverlapNode* second)
{
  return ( 	first->totalSeeds > second->totalSeeds );
}

//On success return the length of merged string
int PBOverlapTree::extendOverlap(std::vector<std::string>& results)
{	
	//Overlap extension via FM-index walk
    while(!m_leaves.empty() && m_currentLength <= m_Query.length()+m_maxIndelSize)
    {
		// if(m_RootNodes.front()->lastOverlapLen == 31)
			// std::cout << m_currentLength << "::\t" << m_leaves.size() <<"\n";
		
		// ACGT-extend the leaf nodes via updating existing SA interval
        extendLeaves();
					
		// Remove leaves without seed support within m_maxIndelSize
		PrunedBySeedSupport();
		
		// speedup by retaining the top bestN candidates after sufficient overlap length
		// This is the 3rd filter less reliable than previous ones
		const size_t bestN = 50;
		if(m_currentLength == 1000 && m_leaves.size()>= bestN)
		{
			m_leaves.sort(SeedComparator);
			SONodePtrList::iterator iter1 = m_leaves.begin();
			SONodePtrList::iterator iter2 = m_leaves.end();
			advance(iter1, bestN-1);
			m_leaves.erase(iter1, iter2);
		}
	
		//see if any interval reach left end $ with sufficient overlap
		if(m_currentLength >= m_minOverlap)
			isTerminated(results);
			
	}
			
	if(!m_leaves.empty())
		for(SONodePtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
		{
			if(computeErrorRate(*iter) < m_errorRate)
			{
				results.push_back( reverse( (*iter)->getFullString() ) );
				std::cout << ">" << computeErrorRate(*iter) << ":" << (*iter)->totalSeeds << ":" << (*iter)->numOfErrors << "\n";
				// std::cout << reverse( (*iter)->getFullString()) << "\n";
			}
		}
	
	return 1;
}


// Print the string represented by every node
void PBOverlapTree::printAll()
{
	for (std::list<SAIOverlapNode*>::iterator it = m_RootNodes.begin(); it != m_RootNodes.end(); ++it)
		(*it)->printAllStrings("");

}


void PBOverlapTree::extendLeaves()
{	
	//attempt to extend one base for each leave
	SONodePtrList::iterator iter = m_leaves.begin();
    while(iter != m_leaves.end())
    {
        std::vector< std::pair<std::string, BWTIntervalPair> > extensions;
        extensions = getLeftFMIndexExtensions(*iter);

        // Either extend the current node or branch it
        // If no extension, do nothing and this node
        // is no longer considered a leaf
        if(extensions.size() == 1)
        {
            // Single extension, do not branch
            (*iter)->extend(extensions.front().first);
            (*iter)->currIntervalPair=extensions.front().second;
			
			// currOverlapLen/queryOverlapLen always increase wrt each extension
			// in order to know the approximate real-time matched length for terminal/containment processing
            (*iter)->currOverlapLen++;
			(*iter)->queryOverlapLen++;
			iter++;
        }
        else if(extensions.size() > 1)
        {
            // Branch for each child
            for(size_t i = 0; i < extensions.size(); ++i)
            {
                SAIOverlapNode* pChildNode = (*iter)->createChild(extensions[i].first);
                pChildNode->currIntervalPair = extensions[i].second;
				pChildNode->currOverlapLen++;
				pChildNode->queryOverlapLen++;
                
				m_leaves.insert(iter, pChildNode);
            }
			iter = m_leaves.erase(iter);
        }
		else	// no extension possible for this leaf, remove from the list
		{
			iter = m_leaves.erase(iter);
		}
    }	
	
    m_currentLength++;  
}

// Check for leaves whose extension has terminated. If the leaf has
// terminated, the walked string and coverage is pushed to the result vector
bool PBOverlapTree::PrunedBySeedSupport()
{
	// the seed index in m_TerminatedIntervals for m_currentLength
	// the m_currentLength is the same for all leaves
	// which is used as the central index within the m_maxIndelSize window
	size_t currSeedIdx = (m_currentLength-m_seedSize) / m_seedDist;
	
	//        ---	seed size =3. seed dist = 1;
	//         --*
	//          -*-
	//           *--
	//            ---
	// *: SNP or indel errors, ---: seed size
	// Erase the leaf if no feasible seeds are found within seedSize+maxIndelSize.
	size_t indelOffset = (m_seedSize+m_maxIndelSize) / m_seedDist;
			
	// Compute the range of small and large indices for tolerating m_maxIndelSize
	size_t smallSeedIdx = currSeedIdx <= indelOffset ? 0 : currSeedIdx - indelOffset;
	size_t largeSeedIdx = (currSeedIdx+indelOffset) >= (m_TerminatedIntervals.size()-1)? 
							m_TerminatedIntervals.size()-1:currSeedIdx+indelOffset;

	// check range of last seed and find new seeds for each interval
	SONodePtrList::iterator iter = m_leaves.begin(); 
    while(iter != m_leaves.end())
    {
		double currErrorRate = computeErrorRate(*iter);

		// if(m_RootNodes.front()->lastOverlapLen == 31 && m_currentLength == 100)
			// std::cout << "currSeedIdx: " << currErrorRate << "::" << smallSeedIdx << "::" 
				// << (*iter)->numOfErrors<< "::" << (*iter)->lastSeedIdx << "::" << (*iter)->totalSeeds << "\n";

		// speedup by skipping dissimilar reads after 200bp
		// This is the 2nd filter less reliable than the 1st one
		if( (m_currentLength > 200 && currErrorRate > m_errorRate+0.04)) 
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
bool PBOverlapTree::isSupportedByNewSeed(SAIOverlapNode* currNode, size_t smallSeedIdx, size_t largeSeedIdx)
{
		// If there is mismatch/indel, jump to the next m_seedSize/m_seedDist, and 1 otherwise.
		size_t seedIdxOffset = currNode->lastOverlapLen < m_currentLength-m_seedSize?
								m_seedSize/m_seedDist:
								//
								m_currentLength - currNode->lastOverlapLen;

		// search for new seed starting from last matched seed or smallSeedIdx
		size_t startSeedIdx = std::max(smallSeedIdx, currNode->lastSeedIdx+seedIdxOffset);
		
		bool isNewSeedFound = false;
		BWTIntervalPair currIntervalPair = currNode->currIntervalPair;
		
		// sequenctial search for new seeds within maxIndelSize
		/*
		// for(size_t i = currNode->lastSeedIdx+seedIdxOffset; i<=largeSeedIdx; i++)
		for(size_t i = startSeedIdx; i<=largeSeedIdx; i++)
		{
			// If seed is a substr of current extension, 
			// the current SA interval is a sub-interval of the seed interval
			isNewSeedFound = currIntervalPair.interval[0].lower >= m_TerminatedIntervals.at(i).lower
						&& currIntervalPair.interval[0].upper <= m_TerminatedIntervals.at(i).upper;

			// This extension is supported by one new seed after lastSeedIdx
			// Return the seed greedily for the first one found, which is inaccurate in tandem repeats
			if(isNewSeedFound)
			{
				std::cout << m_TerminatedIntervals.at(i).size() << "\n";
				// update currNode members
				currNode->lastSeedIdx = i;
				
				// lastOverlapLen records the overlap length of last hit
				currNode->lastOverlapLen = m_currentLength;
				
				// currOverlapLen is always identical to m_currentLength
				currNode->currOverlapLen = m_currentLength;	

				// query overlap may shift due to indels
				currNode->queryOverlapLen = i*m_seedDist+m_seedSize;
				currNode->totalSeeds++;
						
				break;
			}
		}*/

		// Binary search for new seeds using Query interval tree
		std::vector<TreeInterval<size_t> > results;
		BWTIntervalTree.findOverlapping(currIntervalPair.interval[0].lower, currIntervalPair.interval[0].upper, results);
		int minIdxDiff = 10000;
		size_t currSeedIdx = m_currentLength-m_seedSize;

		for(size_t i=0; i< results.size(); i++)
			if( results.at(i).value >= startSeedIdx && results.at(i).value <= largeSeedIdx)
			{
				// if(results.size() > 1) std::cout << results.size() << "\n";
				// update currNode members
				if(std::abs(results.at(i).value - currSeedIdx) < minIdxDiff)
				{					
					currNode->lastSeedIdx = results.at(i).value;
					// query overlap may shift due to indels
					currNode->queryOverlapLen = results.at(i).value+m_seedSize;
					minIdxDiff = std::abs(results.at(i).value - currSeedIdx);
				}

				// lastOverlapLen records the overlap length of last hit
				currNode->lastOverlapLen = m_currentLength;
				
				// currOverlapLen is always identical to m_currentLength
				currNode->currOverlapLen = m_currentLength;	

				currNode->totalSeeds++;
				
				isNewSeedFound = true;
			}
			
		return isNewSeedFound;
}

double PBOverlapTree::computeErrorRate(const SAIOverlapNode* currNode)
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
std::vector<std::pair<std::string, BWTIntervalPair> > PBOverlapTree::getLeftFMIndexExtensions(SAIOverlapNode* pNode)
{
	// if(m_RootNodes.front()->lastOverlapLen == 31 )
		// std::cout << m_currentLength << ">>" << m_leaves.size() << ">>" << reverseComplement(reverse(pNode->getFullString())) << "\n";

	std::vector<std::pair<std::string, BWTIntervalPair> > out;			

	// Only one extension is possble, do direct backward search
	if(pNode->currIntervalPair.interval[0].size()==1)
	{
		int64_t currBWTIndex = pNode->currIntervalPair.interval[0].lower;
		char b = m_pBWT->getChar(currBWTIndex);
		if(b == '$') return out;
		
		//update IntervalPair using extension b
		BWTIntervalPair probe=pNode->currIntervalPair;
		BWTAlgorithms::updateBothL(probe, b, m_pBWT);
		std::string tmp;
		tmp.append(1,b);
		out.push_back(std::make_pair(tmp, probe));

	}
	else
	// Multiple extensions are possible, try FM-index extensions
	{
		for(int i = 1; i < BWT_ALPHABET::size; ++i) //i=A,C,G,T
		{
			char b = BWT_ALPHABET::getChar(i);

			//update IntervalPair using extension b
			BWTIntervalPair probe=pNode->currIntervalPair;
			BWTAlgorithms::updateBothL(probe, b, m_pBWT);
				
			//min freq at fwd and rvc bwt
			if(probe.isValid())
			{
				// if(m_RootNodes.front()->lastOverlapLen == 31 )
					// std::cout << b << ":" << probe.interval[0].size() << "\n";
				
				if(m_currentLength <= m_pRootNode->lastOverlapLen + 8 && probe.interval[0].size() <= 8) 
					continue;
				
				std::string tmp;
				tmp.append(1,b);
				out.push_back(std::make_pair(tmp, probe));
			}
		}// end of ACGT
	}

	return out;
}

// Check if the leaves reach $
bool PBOverlapTree::isTerminated(std::vector<std::string>& results)
{
	bool found=false;
    for(SONodePtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
    {
		// Calculate which of the prefixes that match w[i, l] are terminal
		// These are the proper prefixes (they are the start of a read)
		BWTIntervalPair probe = (*iter)->currIntervalPair;
		BWTAlgorithms::updateBothL(probe, '$', m_pBWT);
		double errorRate = computeErrorRate(*iter);
		
		// The SA interval reach $
		if(probe.isValid() && (*iter)->queryOverlapLen >= m_minOverlap 
			&& (*iter)->queryOverlapLen < m_Query.length())
		{			
			if(errorRate >= m_errorRate) continue;

			std::cout << ">>" << errorRate << ":" << (*iter)->totalSeeds << ":" << m_currentLength << "\n" ;
			// std::cout << reverse( (*iter)->getFullString() ) << "\n";
						
			results.push_back( reverse( (*iter)->getFullString() ));
			found = true;
		}
	}
	
	return found;

}

