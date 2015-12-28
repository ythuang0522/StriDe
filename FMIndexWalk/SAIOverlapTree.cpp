///----------------------------------------------
// Copyright 2015 National Chung Cheng University
// Written by Yao-Ting Huang
// Released under the GPL
//-----------------------------------------------
//
// SAIOverlapTree - Iteratively construct a
// search tree representing an FM-index walk through an assembly graph
// matching a query sequence.
//
// The assembly graph is abstractly represented as
// an FM-index.
//
#include "SAIOverlapTree.h"
#include "BWTAlgorithms.h"

//
// Class: SAIOverlapTree
SAIOverlapTree::SAIOverlapTree(const std::string& query,
			   size_t minOverlap,
			   size_t maxIndelSize,
			   const BWT* pBWT, const BWT* pRBWT, 
			   AlignFlags af,
			   size_t maxLeaves,
			   size_t seedSize, size_t seedDist, 
			   size_t repeatFreq):
                               m_Query(query), m_minOverlap(minOverlap),  
                               m_maxIndelSize(maxIndelSize), m_pBWT(pBWT), m_pRBWT(pRBWT), m_af(af),
                               m_maxLeaves(maxLeaves), 
							   m_seedSize(seedSize), m_seedDist(seedDist), 
							   m_repeatFreq(repeatFreq)
{
    // Create the root node containing the rightmost seed string
	// std::string initStr = m_Query.substr(m_Query.length()-m_seedSize, m_seedSize);
	// BWTIntervalPair bip;

	// the initial seed may contain errors already, which is in real data
	// greedily find the next best seed until m_seedSize+m_maxIndelSize
	for(size_t initSeedOffset = 0; initSeedOffset<m_seedSize+m_maxIndelSize
		&& initSeedOffset+m_seedSize<=m_Query.length()	//bug fix: some short reads should be skipped
		; initSeedOffset++)
	{
		std::string initSeed = m_Query.substr(m_Query.length()-m_seedSize-initSeedOffset, m_seedSize);
		BWTIntervalPair bip = BWTAlgorithms::findIntervalPair( m_pBWT, m_pRBWT, initSeed);

		// feasible initial seeds are found if it's not a repeat seed
		if(bip.isValid() && bip.interval[0].size() < 256)
		{
			// std::cout << bip.interval[0].size() << "\n";
				
			// create one root node
			SAIOverlapNode* m_pRootNode = new SAIOverlapNode(&initSeed, NULL);
			m_pRootNode->currIntervalPair = bip;
			m_pRootNode->computeInitial(reverse(initSeed));   //store initial str of root
			m_pRootNode->lastOverlapLen = m_pRootNode->currOverlapLen = m_pRootNode->queryOverlapLen 
						= m_currentLength = m_seedSize+initSeedOffset;

			m_pRootNode->lastSeedIdx = m_pRootNode->initSeedIdx = initSeedOffset;
			m_pRootNode->totalSeeds = 1;
			
			// push new node into roots and leaves vector
			m_RootNodes.push_back(m_pRootNode);
			m_leaves.push_back(m_pRootNode);

			// initialize the seeding SA intervals with seedSize and seedDist from rightmost to leftmost
			for(int i = (int)m_Query.length()-m_seedSize; i >= 0 ; i -= (int)m_seedDist)
			{
				std::string seedStr = m_Query.substr(i, m_seedSize);
				m_TerminatedIntervals.push_back(BWTAlgorithms::findInterval( m_pBWT, seedStr) );		
			}
			
			break;
		}
	}
}

//
SAIOverlapTree::~SAIOverlapTree()
{
	for (std::list<SAIOverlapNode*>::iterator it = m_RootNodes.begin(); it != m_RootNodes.end(); ++it)
		delete *it;
	
	m_RootNodes.clear();
    // Recursively destroy the tree
    // delete m_pRootNode;
}

//On success return the length of merged string
int SAIOverlapTree::extendOverlapOneBase(std::vector<OverlapBlock>& results)
{	
	//Overlap extension via FM-index walk
    if(!m_leaves.empty() && m_leaves.size() <= m_maxLeaves && 
		m_currentLength <= m_Query.length()+m_maxIndelSize)
    {
		// if(m_currentLength==3101){
			// std::cout << m_Query.at(m_Query.length()-m_currentLength) << ":" << m_currentLength << 
				// ":" << m_leaves.size() << ":" << m_leaves.front()->currIntervalPair.interval[0]<<"\n";	

			// printAll();
			// getchar();
		// }

		// ACGT-extend the leaf nodes via updating existing SA interval
        extendLeaves();
		
		// Initial seeding region may still contain errors
		// Add seeds at m_currentLength
		if(m_currentLength == m_seedSize*2)
			addNewRootNodes();
			
		// Remove leaves without seed support within m_maxIndelSize
		PrunedBySeedSupport();

		
		//see if any interval reach left end $ with sufficient overlap
		if(m_currentLength >= m_minOverlap)
			isTerminated(results);
			
	}
			
    //Did not reach the terminal kmer
    if(m_leaves.empty())
        return -1;	//high error
    else if(m_currentLength>m_Query.length()+m_maxIndelSize)
        return -2;	//exceed search depth
    else if(m_leaves.size() > m_maxLeaves)
        return -3;	//too much repeats
	else
		return 1;
}


// Print the string represented by every node
void SAIOverlapTree::printAll()
{
	for (std::list<SAIOverlapNode*>::iterator it = m_RootNodes.begin(); it != m_RootNodes.end(); ++it)
		(*it)->printAllStrings("");
    // std::cout << "Print all: \n";
    // m_pRootNode->printAllStrings("");
}

// Extend each leaf node

void SAIOverlapTree::attempToExtend(SONodePtrList &newLeaves)
{
    for(SONodePtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
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
			
            newLeaves.push_back(*iter);
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
                
				newLeaves.push_back(pChildNode);
            }
        }
    }	
}

void SAIOverlapTree::extendLeaves()
{
    SONodePtrList newLeaves;
	
	//attempt to extend one base for each leave
    attempToExtend(newLeaves);
		
    m_currentLength++;  

    m_leaves.clear();
    m_leaves = newLeaves;
}

void SAIOverlapTree::addNewRootNodes()
{
	std::string	initStr = m_Query.substr(m_Query.length()-m_currentLength, m_seedSize);
	BWTIntervalPair bip = BWTAlgorithms::findIntervalPair( m_pBWT, m_pRBWT, initStr);
	
	if(bip.isValid() && (size_t)bip.interval[0].size() < m_repeatFreq)
	{
		// std::cout << initStr << bip.interval[0].size() << "\n";
		
	    SAIOverlapNode* m_pRootNode = new SAIOverlapNode(&initStr, NULL);

		m_pRootNode->currIntervalPair = bip;
		m_pRootNode->computeInitial(reverse(initStr));   //store initial str of root
		m_pRootNode->initSeedIdx = (m_currentLength-m_seedSize)/m_seedDist;
		m_pRootNode->lastSeedIdx = m_pRootNode->initSeedIdx-1;
		m_pRootNode->lastOverlapLen = m_pRootNode->currOverlapLen = m_pRootNode->queryOverlapLen = m_currentLength;
		m_pRootNode->totalSeeds = 1;
			
		// the redundant root nodes may lead to false substring for reads containing tandem repeats AAAAA...
		m_RootNodes.push_back(m_pRootNode);
		m_leaves.push_back(m_pRootNode);

	}

}

// Check for leaves whose extension has terminated. If the leaf has
// terminated, the walked string and coverage is pushed to the result vector
bool SAIOverlapTree::PrunedBySeedSupport()
{
	SONodePtrList newLeaves;
	// the seed index in m_TerminatedIntervals for m_currentLength
	// the m_currentLength is the same for all intervals
	// which is used as the central index within the m_maxIndelSize window
	//           ---
	//             ---	seedDist < seed Size
	// SSSDDDSSSDDDSSS
	//        |||||| -> seedDist > seed Size
	//  ||||||-> next seed idx coverage, (9-3)/(3+3)=1
	size_t currSeedIdx = (m_currentLength-m_seedSize) / m_seedDist;
	
	//        ---	seed size =3. seed dist = 1;
	//         --*
	//          -*-
	//           *--
	//            ---
	// *: SNP or indel errors, ---: seed size
	// Worst case: no feasible seeds are found for seedSize+maxIndelSize.
	size_t indelOffset = (m_seedSize+m_maxIndelSize) / m_seedDist;
	
	// If requiring no errors, indelOffset should be set to zero.
		
	// Compute the range of small and large indices for tolerating m_maxIndelSize
	size_t smallSeedIdx = currSeedIdx <= indelOffset ? 0 : currSeedIdx - indelOffset;
	size_t largeSeedIdx = (currSeedIdx+indelOffset) >= (m_TerminatedIntervals.size()-1)? 
							m_TerminatedIntervals.size()-1:currSeedIdx+indelOffset;

	// check range of last seed and find new seeds for each interval
    for(SONodePtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
    {	
		// std::cout << (*iter)->getFullString() << "\n";
		// if(m_currentLength>=270)
			// std::cout << "currSeedIdx: " << currSeedIdx << "::" << smallSeedIdx << "::" 
					// << (*iter)->initSeedIdx << "::" << (*iter)->lastSeedIdx << "::" << (*iter)->totalSeeds << "\n";

		// Require at least one seed within smallSeedIdx~largeSeedIdx i.e., within m_maxIndelSize
		bool isLastSeedWithinRange = (*iter)->lastSeedIdx >= smallSeedIdx && 
									 (*iter)->lastSeedIdx <= largeSeedIdx;

		// search for matched new seeds
		bool isNewSeedFound = isSupportedByNewSeed(*iter, largeSeedIdx);

		// the seed is overlapped with one SNP or indel error for the 1st time
		// while previous seed is a match, increment the error number
		if(isNewSeedFound)
			(*iter)->lastSeedIdxOffset = (int)(*iter)->lastSeedIdx - (int)currSeedIdx;

		if(!isNewSeedFound && currSeedIdx+(*iter)->lastSeedIdxOffset == (*iter)->lastSeedIdx+1)
			(*iter)->numOfErrors++;

		// give up if new seeds do not satisfy indel constraints
		// this is not necessary if isSupportedByNewSeed is implemented 100% accurately
		// if(!isOverlapAcceptable(*iter)) continue;

		// The current extension must be supported by sharing at least one seed with m_maxIndelSize		
		if( (isLastSeedWithinRange || isNewSeedFound))
			newLeaves.push_back(*iter);

    }
	
	// update m_leaves with newLeaves;
	bool updated = false;
	if(m_leaves.size() != newLeaves.size())
	{
		m_leaves = newLeaves;
		updated=true;
	}

    return updated;
}

// Identify new seeds wrt currSeedIdx
bool SAIOverlapTree::isSupportedByNewSeed(SAIOverlapNode* currNode, size_t largeSeedIdx)
{
		/* The starting idx must reflect the indel size
			TATATATATATATAGGAAGGG	seed size =3, dist =1
			currSeedIdx: 18::12::24::14::15
			17:     4       15			
			TATATATATATATAGGGGGGG
			currSeedIdx: 18::12::24::18::19
			21:     0       19
			
			// last GGG matched to wrong GGGs, requiring large seeds and feasible index
			TATATATATATATAGGAAGGG
			TATATATATATATAGG----GGGGG
		*/
		// size_t currOverlapDiff = m_currentLength - currNode->lastOverlapLen;
		// If there is mismatch/indel, jump to the next m_seedSize, and 1 otherwise.
		size_t seedIdxOffset = currNode->lastOverlapLen<m_currentLength-m_seedSize?
								m_seedSize/m_seedDist:
								m_currentLength-currNode->lastOverlapLen-1;

		bool isNewSeedFound = false;
		BWTIntervalPair currIntervalPair = currNode->currIntervalPair;
		
		// search for new seeds within maxIndelSize
		for(size_t i = currNode->lastSeedIdx+seedIdxOffset; i<=largeSeedIdx; i++)
		{
			// If seed is a substr of current extension, 
			// the current SA interval is a sub-interval of the seed interval
			isNewSeedFound = currIntervalPair.interval[0].lower >= m_TerminatedIntervals.at(i).lower
						&& currIntervalPair.interval[0].upper <= m_TerminatedIntervals.at(i).upper;

			// isNewSeedFound = Interval::isIntersecting(currIntervalPair.interval[0].lower, currIntervalPair.interval[0].upper, 
                                    // m_TerminatedIntervals.at(i).lower, m_TerminatedIntervals.at(i).upper);
			// This extension is supported by one new seed after lastSeedIdx
			// Return the seed greedily for the first one found, which is inaccurate in tandem repeats
			if(isNewSeedFound)
			{
				// update currNode members
				currNode->lastSeedIdx = i;
				currNode->lastOverlapLen = m_currentLength;
				currNode->currOverlapLen = m_currentLength;				
				currNode->queryOverlapLen = i*m_seedDist+m_seedSize;
				currNode->totalSeeds++;
						
				break;
			}
		}

		return isNewSeedFound;
}

bool SAIOverlapTree::isOverlapAcceptable(SAIOverlapNode* currNode)
{
	// Compute the overlap length wrt this interval
	// The overlap should be adjusted wrt the last query seed pos
	// size_t lastQuerySeedPos = currNode->lastSeedIdx*m_seedDist + m_seedSize;
	// The offset should be computed according to last subject seed pos
	// size_t lastQuerySeedPosOffset = m_currentLength - currNode->lastOverlapLen;
	// size_t queryOverlapLen = lastQuerySeedPos+lastQuerySeedPosOffset;
	
	// Estimate total errors by overlap length diff between query and subject
	size_t overlapDiff = std::abs((int)currNode->queryOverlapLen-(int)m_currentLength);

	// std::cout << queryOverlapLen << ":\t" << overlapDiff << "\t" << currNode->totalSeeds <<"\n";
	
	if(overlapDiff > m_maxIndelSize) return false;
	else return true;
}

double SAIOverlapTree::computeErrorRate(SAIOverlapNode* currNode)
{	
	// Compute accuracy via matched length in both query and subject
	// SNP and indel will over-estimate the unmatched regions across error, --*--
	int matchedLen = (int)currNode->totalSeeds*2;
	// matchedLen += (int)(currNode->numOfErrors*(m_seedSize-1)*2) ;
	int totalLen = (int)currNode->queryOverlapLen + currNode->currOverlapLen- (m_seedSize*2) +2;
	int unmatchedLen = totalLen - matchedLen;

	// if(unmatchedLen<0 || totalLen<0)
	// if(m_currentLength>=16570)
	// {
		// std::cout << currNode->queryOverlapLen  << "\t" << currNode->currOverlapLen
			 // << "\t" << currNode->totalSeeds << "\n";
	
		// std::cout << (double)unmatchedLen/totalLen << "\t" << matchedLen << "\t" << unmatchedLen 
				// << "\t" << totalLen << "\t"<< currNode->numOfErrors <<"\n";
	// }
	
	return unmatchedLen/(double)totalLen;
}
			
//update SA intervals of each leaf, which corresponds to one-base extension
std::vector<std::pair<std::string, BWTIntervalPair> > SAIOverlapTree::getLeftFMIndexExtensions(SAIOverlapNode* pNode)
{
    std::vector<std::pair<std::string, BWTIntervalPair> > out;

    for(int i = 1; i < BWT_ALPHABET::size; ++i) //i=A,C,G,T
    {
        char b = BWT_ALPHABET::getChar(i);

        //update IntervalPair using extension b
        BWTIntervalPair probe=pNode->currIntervalPair;
        BWTAlgorithms::updateBothL(probe, b, m_pBWT);
			
		//min freq at fwd and rvc bwt
        if(probe.isValid())
        {
            std::string tmp;
            tmp.append(1,b);
            BWTIntervalPair bip=probe;
            out.push_back(std::make_pair(tmp, bip));
        }
    }// end of ACGT

    return out;
}

//update SA intervals of each leaf, which corresponds to one-base extension
std::vector<std::pair<std::string, BWTIntervalPair> > SAIOverlapTree::getRightFMIndexExtensions(SAIOverlapNode* pNode)
{
    std::vector<std::pair<std::string, BWTIntervalPair> > out;

    for(int i = 1; i < BWT_ALPHABET::size; ++i) //i=A,C,G,T
    {
        char b = BWT_ALPHABET::getChar(i);

        //update IntervalPair using extension b
        BWTIntervalPair probe=pNode->currIntervalPair;
        BWTAlgorithms::updateBothR(probe, b, m_pRBWT);
			
		//min freq at fwd and rvc bwt
        if(probe.isValid())
        {
            std::string tmp;
            tmp.append(1,b);
            BWTIntervalPair bip=probe;
            out.push_back(std::make_pair(tmp, bip));
        }
    }// end of ACGT

    return out;
}

// Check if the leaves reach $
bool SAIOverlapTree::isTerminated(std::vector<OverlapBlock>& results)
{
	bool found=false;
	SONodePtrList newLeaves;
    for(SONodePtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
    {
		// Calculate which of the prefixes that match w[i, l] are terminal
		// These are the proper prefixes (they are the start of a read)
		BWTIntervalPair probe = (*iter)->currIntervalPair;
		BWTAlgorithms::updateBothL(probe, '$', m_pBWT);
		
		// The SA interval reach $
		if(probe.isValid() && (*iter)->queryOverlapLen >= m_minOverlap 
			&& (*iter)->queryOverlapLen < m_Query.length())
		{
			// substring lead to false positive overlap
			// ---------	query
			//    ----		substr
			//    ------	substr containment
			//    --------	normal overlap
			
			std::list<BWTIntervalPair> substrTerminatedReads;
			std::list<BWTIntervalPair> normalTerminatedReads = collectToRightExtreme(probe, (*iter)->initSeedIdx, substrTerminatedReads);
			
			double errorRate = computeErrorRate(*iter);
			// std::cout << errorRate << "\n";
			
			if(errorRate >= 0.3) continue;

			size_t totalErrors = errorRate * m_Query.length()*2;

			// w: AACCCTTTTGG
			//       CCTT--GG
			size_t insertionSize = (*iter)->queryOverlapLen>=m_currentLength?
									(*iter)->queryOverlapLen-m_currentLength:0;
			// w: AACCCTT--GG
			//       CCTTTTGG
			size_t deletionSize = (*iter)->queryOverlapLen<m_currentLength?
									m_currentLength-(*iter)->queryOverlapLen:0;

			// normal overlap
			for(std::list<BWTIntervalPair>::iterator it=normalTerminatedReads.begin(); it!=normalTerminatedReads.end(); it++)
			{	
				// std::cout << (*iter)->getFullString();
				// if(m_currentLength>=3071)
				// {
					// std::cout << "Terminated: " << (*iter)->queryOverlapLen << "\t" << totalErrors << "\t" << insertionSize 
					// << "\t" << deletionSize << "\t" << (*iter)->initSeedIdx << "\t"<<(*it).interval[0]<< "\t"<< errorRate << "\n";
					// getchar();
				// }
				
				OverlapBlock ob(*it, (*iter)->currIntervalPair, (*iter)->queryOverlapLen,
									totalErrors, insertionSize, deletionSize, m_af);
								
				results.push_back( ob );
									
				found = true;
			}
			
			// substring
			for(std::list<BWTIntervalPair>::iterator it=substrTerminatedReads.begin(); it!=substrTerminatedReads.end(); it++)
			{
				// std::cout << "Terminated substring: " << (*iter)->queryOverlapLen << "\t" << totalErrors  << "\t" << insertionSize 
					// << "\t" << deletionSize << "\t" << (*iter)->initSeedIdx << "\t"<<(*it).interval[0]<<"\n";

				OverlapBlock ob(*it, (*iter)->currIntervalPair, (*iter)->queryOverlapLen,
									totalErrors, insertionSize, deletionSize, m_af);
								
				ob.isTargetSubstring = true;					
				results.push_back( ob );
				found = true;
			}
		}
	}
	
	return found;

}

// return true if it's a substring
bool SAIOverlapTree::terminateContainedBlocks(std::vector<OverlapBlock>& results)
{
	SONodePtrList newLeaves;
    for(SONodePtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
    {
		if((*iter)->queryOverlapLen < m_Query.length())
		{	
			newLeaves.push_back(*iter);
			continue;
		}

		double errorRate = computeErrorRate(*iter);
		if( errorRate<0.3)
		{
			BWTIntervalPair ranges = (*iter)->currIntervalPair;
			AlphaCount64 left_ext = BWTAlgorithms::getExtCount(ranges.interval[0], m_pBWT);
			AlphaCount64 right_ext = BWTAlgorithms::getExtCount(ranges.interval[1], m_pRBWT);

			size_t totalErrors = errorRate * m_Query.length()*2;
			size_t insertionSize = (*iter)->queryOverlapLen>=m_currentLength?(*iter)->queryOverlapLen-m_currentLength:0;
			size_t deletionSize = (*iter)->queryOverlapLen<m_currentLength?m_currentLength-(*iter)->queryOverlapLen:0;

			// left substring
			if(left_ext.hasDNAChar())
			{
				// collect left extreme reads within m_maxIndelSize
				//    -----		-> query
				//  ------		case 1, remove it by extendToRightExtreme
				//  --------	case 2, substring
				std::list<BWTIntervalPair> rightTerminatedReads = extendToRightExtreme((*iter)->currIntervalPair, (*iter)->initSeedIdx);
				std::list<BWTIntervalPair> leftTerminatedReads, bothTerimnatedReads;
				bool isLeftSubstring = false;
				for(std::list<BWTIntervalPair>::iterator it=rightTerminatedReads.begin(); it!=rightTerminatedReads.end(); it++)
				{
					// may need to consider deletion size
					leftTerminatedReads = extendToLeftExtreme((*iter)->currIntervalPair, m_maxIndelSize, isLeftSubstring);
					
					//substring still has DNA char after m_maxIndelSize
					if(isLeftSubstring)
					{
						// std::cout << "Absolutely left substring:\t" << (*iter)->queryOverlapLen << ":\t" << totalErrors << "\t" << insertionSize 
						// << "\t" << deletionSize<< "\t"<< ranges.interval[0] << "\t" 
						// << (*iter)->initSeedIdx << "\t" << m_leaves.size() << "\n";
						
						return true;
					}
					else if(!leftTerminatedReads.empty())
						bothTerimnatedReads.splice(bothTerimnatedReads.end(), leftTerminatedReads);
				}

				// set overlap length to query length +1 for indicating substring
				for(std::list<BWTIntervalPair>::iterator it=bothTerimnatedReads.begin(); it!=bothTerimnatedReads.end(); it++)
				{
					// push into containment list first as this may not be the optimal alignments
					// subMaximal filter will remove these reads
					// std::cout << "has left substring:\t" << (*iter)->queryOverlapLen << ":\t" << totalErrors << "\t" << insertionSize 
					// << "\t" << deletionSize<< "\t"<< (*it).interval[0] << "\t" 
					// << (*iter)->initSeedIdx << "\t" << m_leaves.size() << "\n";

					results.push_back( OverlapBlock(*it, ranges, m_Query.length()+1,
									totalErrors, insertionSize, deletionSize, m_af) ) ;
				}
			}
			else if(right_ext.hasDNAChar())
			{
				//   -----	query
				//   ---	case 1
				//   -----  case 2
				//   ------ case 3
				
				//   ---	case 1: remove
				std::list<BWTIntervalPair> containments = extendToRightExtreme((*iter)->currIntervalPair, (*iter)->initSeedIdx);
				for(std::list<BWTIntervalPair>::iterator it1=containments.begin(); it1!=containments.end(); it1++)
				{
					BWTIntervalPair probe1 = ranges;
					BWTIntervalPair probe2 = ranges;
					BWTAlgorithms::updateBothL(probe1, '$', m_pBWT);
					BWTAlgorithms::updateBothR(probe2, '$', m_pRBWT);
					if(probe1.isValid() && probe2.isValid())
					{
						//   -----  case 2: containment
						// std::cout << "has right containment:\t"<< (*iter)->queryOverlapLen << ":\t" << totalErrors << "\t" << insertionSize 
						// << "\t" << deletionSize<< "\t"<< ranges.interval[0].size() << "\t" 
						// << (*iter)->initSeedIdx << "\t" << m_leaves.size() << "\n";

						results.push_back( OverlapBlock(probe1, ranges, m_Query.length(),
											totalErrors, insertionSize, deletionSize, m_af)) ;
						
					}
					else
					{
						assert(probe1.isValid());
						//   ------ case 3: check if right substring
						std::list<BWTIntervalPair> rightTerminals = extendToRightExtreme(*it1, 1);
						for(std::list<BWTIntervalPair>::iterator it2=rightTerminals.begin(); it2!=rightTerminals.end(); it2++)
						{
							// terminate the contained block and add it to the contained list
							// std::cout << "has right substring:\t"<< (*iter)->queryOverlapLen << ":\t" << totalErrors << "\t" << insertionSize 
							// << "\t" << deletionSize<< "\t"<< probe1.interval[0] << "\t" 
							// << (*iter)->initSeedIdx << "\t" << m_leaves.size() << "\n";
							
							// set overlap length to query length +1 for indicating substring
							results.push_back( OverlapBlock(probe1, ranges, m_Query.length()+1,
											totalErrors, insertionSize, deletionSize, m_af)) ;
						}

					}
				}
			}
			else
			{
				BWTIntervalPair probe = ranges;
				BWTAlgorithms::updateBothL(probe, '$', m_pBWT);

				//   -----	query
				//   ---	case 1, caused by initSeedIdx!=0
				//   -----  case 2
				if(probe.isValid() && (*iter)->initSeedIdx==0)
				{
					// terminate the contained block and add it to the contained list
					// BWTAlgorithms::updateBothR(probe, '$', m_pRBWT);
					// assert(probe.isValid());

					// std::cout << "containment:\t" << (*iter)->queryOverlapLen << ":\t" << totalErrors << "\t" << insertionSize 
							// << "\t" << deletionSize<< "\t"<< probe.interval[0] << "\t" 
							// << (*iter)->initSeedIdx << "\t" << m_leaves.size() << "\n";
					results.push_back( OverlapBlock(probe, ranges, m_Query.length(),
										totalErrors, insertionSize, deletionSize, m_af)) ;
					// getchar();
				}
			}
		}// end of overlap length > m_Query
	}// end of for each leave
	
	m_leaves.clear();
	m_leaves=newLeaves;
	return false;
}

// extend to left or right extreme within maxIndelSize
std::list<BWTIntervalPair> SAIOverlapTree::extendToLeftExtreme(BWTIntervalPair& currIntervalPair, int length, bool& isLeftSubstring)
{
	// assert(length > 0);
	std::list<BWTIntervalPair> currbips, results;
	currbips.push_back(currIntervalPair);
	for(int idx=0; idx<length; idx++)
	{
		std::list<BWTIntervalPair> newcurrbips;

		for(std::list<BWTIntervalPair>::iterator iter=currbips.begin(); iter!=currbips.end(); iter++)
		{

			for(int j = 1; j < BWT_ALPHABET::size; ++j) //i=A,C,G,T
			{
				char b = BWT_ALPHABET::getChar(j);

				//update IntervalPair using extension b
				BWTIntervalPair probe = *iter;
				BWTAlgorithms::updateBothL(probe, b, m_pBWT);
					
				if(probe.isValid())
					newcurrbips.push_back(probe);
			}// end of ACGT
		}

		// Extension fails
		if(newcurrbips.empty()) return results;
		
		// check for terminating intervals
		for(std::list<BWTIntervalPair>::iterator iter2=newcurrbips.begin(); iter2!=newcurrbips.end(); iter2++)
		{
			BWTIntervalPair probe = *iter2;
			BWTAlgorithms::updateBothL(probe, '$', m_pBWT);
					
			if(probe.isValid())
				results.push_back(probe);
		}

		currbips.clear();
		currbips=newcurrbips;
		newcurrbips.clear();
	}
	
	for(std::list<BWTIntervalPair>::iterator iter=currbips.begin(); iter!=currbips.end(); iter++)
	{
		AlphaCount64 left_ext = BWTAlgorithms::getExtCount((*iter).interval[0], m_pBWT);		
		if(left_ext.hasDNAChar())
			isLeftSubstring=true;
	}
		
	return results;
}

std::list<BWTIntervalPair> SAIOverlapTree::extendToRightExtreme(BWTIntervalPair& currIntervalPair, int length)
{
	// assert(length > 0);
	std::list<BWTIntervalPair> currbips;
	currbips.push_back(currIntervalPair);

	// Boundary case
	if(length == 0)
		return currbips;

	// Extend until length is reached
	for(int idx=0; idx<length; idx++)
	{
		std::list<BWTIntervalPair> newcurrbips;

		for(std::list<BWTIntervalPair>::iterator iter=currbips.begin(); iter!=currbips.end(); iter++)
		{
			for(int j = 1; j < BWT_ALPHABET::size; ++j) //i=A,C,G,T
			{
				char b = BWT_ALPHABET::getChar(j);

				//update IntervalPair using extension b
				BWTIntervalPair probe = *iter;
				BWTAlgorithms::updateBothR(probe, b, m_pRBWT);
					
				if(probe.isValid())
					newcurrbips.push_back(probe);
			}// end of ACGT
		}

		// Extension fails
		if(newcurrbips.empty()) return newcurrbips;
		
		currbips.clear();
		currbips=newcurrbips;
		newcurrbips.clear();
	}

	return currbips;
}

// ---------	query
//    ----		substr
//    ------	substr containment
//    --------	normal overlap
std::list<BWTIntervalPair> SAIOverlapTree::collectToRightExtreme(BWTIntervalPair& currIntervalPair, int length, std::list<BWTIntervalPair>& terminatedReads)
{
	std::list<BWTIntervalPair> currbips;
	currbips.push_back(currIntervalPair);

	BWTIntervalPair probe = currIntervalPair;
	BWTAlgorithms::updateBothR(probe, '$', m_pRBWT);
	if(probe.isValid())
		terminatedReads.push_back(probe);

	// Boundary case
	if(length == 0)
		return currbips;
	
	// Extend until length is reached
	for(int idx=0; idx<length; idx++)
	{
		std::list<BWTIntervalPair> newcurrbips;

		for(std::list<BWTIntervalPair>::iterator iter=currbips.begin(); iter!=currbips.end(); iter++)
		{
			for(int j = 1; j < BWT_ALPHABET::size; ++j) //i=A,C,G,T
			{
				char b = BWT_ALPHABET::getChar(j);

				//update IntervalPair using extension b
				BWTIntervalPair probe = *iter;
				BWTAlgorithms::updateBothR(probe, b, m_pRBWT);
					
				if(probe.isValid())
					newcurrbips.push_back(probe);
			}// end of ACGT
		}

		// Extension fails
		if(newcurrbips.empty()) return newcurrbips;
		
		// Collect terminated reads
		for(std::list<BWTIntervalPair>::iterator iter=currbips.begin(); iter!=currbips.end(); iter++)
		{
			BWTIntervalPair probe = *iter;
			BWTAlgorithms::updateBothR(probe, '$', m_pRBWT);
			if(probe.isValid())
				terminatedReads.push_back(probe);
		}

		currbips.clear();
		currbips=newcurrbips;
		newcurrbips.clear();
	}

	return currbips;
}