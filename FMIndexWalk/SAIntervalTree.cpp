///----------------------------------------------
// Copyright 2014 National Chung Cheng University
// Written by Yao-Ting Huang
// Released under the GPL
//-----------------------------------------------
//
// SAIntervalTree - Iteratively construct a
// string representing a walk through an assembly graph
// matching a query sequence.
//
// The assembly graph is abstractly represented as
// an FM-index.
//
#include "SAIntervalTree.h"
#include "BWTAlgorithms.h"

//
// SAIntervalNode
//
SAIntervalNode::SAIntervalNode(const std::string* pQuery, SAIntervalNode* parent) : 
									   m_kmerCount(0), m_pQuery(pQuery),m_pParent(parent)
{

}

// Destructor, recurisvely delete the children of the node
SAIntervalNode::~SAIntervalNode()
{
    // Delete children
    for(STNodePtrList::iterator iter = m_children.begin(); iter != m_children.end(); ++iter)
        delete *iter;

}

// Return a suffix of length l of the path from the root to this node
std::string SAIntervalNode::getSuffix(size_t l) const
{
    size_t n = m_label.size();
    if(l <= n)
    {
        return m_label.substr(n - l, l);
    }
    else
    {
        assert(m_pParent != NULL);
        return m_pParent->getSuffix(l - n) + m_label;
    }
}

// Return the full string of the path from the root to this node
std::string SAIntervalNode::getFullString() const
{
    if(m_pParent == NULL)
        return m_label;
    else
        return m_pParent->getFullString() + m_label;
}

// Create a new child node with the given label. Returns a pointer to the new node.
SAIntervalNode* SAIntervalNode::createChild(const std::string& label)
{
    SAIntervalNode* pAdded = new SAIntervalNode(m_pQuery, this);
    m_children.push_back(pAdded);

    //assert(!m_alignmentColumns.empty());
    //pAdded->computeExtendedAlignment(label, m_alignmentColumns.back());
    pAdded->extend(label);

    return pAdded;
}

// Extend the label of this node
void SAIntervalNode::extend(const std::string& ext)
{
    assert(!ext.empty());
    //assert(!m_alignmentColumns.empty());
    m_label.append(ext);
}


void SAIntervalNode::computeInitial(const std::string& initialLabel)
{
    m_label = initialLabel;

}


// Print the string(s) represented by this node and its children
void SAIntervalNode::printAllStrings(const std::string& parent) const
{
    if(m_children.empty())
    {
        std::cout << ">\n" << parent + m_label << "\n";
    }
    else
    {
        for(STNodePtrList::const_iterator iter = m_children.begin(); iter != m_children.end(); ++iter)
            (*iter)->printAllStrings(parent + m_label);
    }
}

//
// Class: SAIntervalTree
SAIntervalTree::SAIntervalTree(const std::string* pQuery,
                               size_t minOverlap,
							   size_t maxOverlap,
                               size_t MaxLength,
                               size_t MaxLeaves,
                               BWTIndexSet indices,
                               std::string secondread,
                               size_t SA_threshold,
                               bool KmerMode) :
                               m_pQuery(pQuery), m_minOverlap(minOverlap), m_maxOverlap(maxOverlap), m_MaxLength(MaxLength),
                               m_MaxLeaves(MaxLeaves), m_indices(indices), 
                               m_secondread(secondread), m_min_SA_threshold(SA_threshold),
                                m_kmerMode(KmerMode), m_maxKmerCoverage(0), m_maxUsedLeaves(0), m_isBubbleCollapsed(false)
{
    // Create the root node containing the seed string
    m_pRootNode = new SAIntervalNode(pQuery, NULL);
    m_pRootNode->computeInitial(*pQuery);   //store initial str of root
    m_leaves.push_back(m_pRootNode);

    m_currentLength=pQuery->length();
	m_currentKmerSize=m_minOverlap;

	//beginning kmer is a suffix of first read
    //initialize the beginning kmer SA intervals with kmer length=m_minOverlap
    std::string beginningkmer=pQuery->substr(m_currentLength-m_minOverlap);
    m_pRootNode->fwdInterval=BWTAlgorithms::findInterval( m_indices.pRBWT, reverse(beginningkmer));
    m_pRootNode->rvcInterval=BWTAlgorithms::findInterval( m_indices.pBWT, reverseComplement(beginningkmer));

	//ending kmer is a prefix of second read
    //initialize the ending SA intervals with kmer length=m_minOverlap
    std::string endingkmer=secondread.substr(0,m_minOverlap);
    m_fwdTerminatedInterval=BWTAlgorithms::findInterval( m_indices.pRBWT, reverse(endingkmer));
    m_rvcTerminatedInterval=BWTAlgorithms::findInterval( m_indices.pBWT, reverseComplement(endingkmer));

	//std::cout << m_minOverlap << ":" << beginningkmer << ":" << endingkmer << "\n";
	//getchar();
}

//
SAIntervalTree::~SAIntervalTree()
{
    // Recursively destroy the tree
    delete m_pRootNode;
}

//On success return the length of merged string
int SAIntervalTree::mergeTwoReads(std::string &mergedseq)
{
    SAIntervalNodeResultVector results;
	
    if( isTwoReadsOverlap(mergedseq))
		return 1;

	//BFS search from 1st to 2nd read via FM-index walk
    while(!m_leaves.empty() && m_leaves.size() <= m_MaxLeaves && m_currentLength <=m_MaxLength)
    {
        // ACGT-extend the leaf nodes via updating existing SA interval
        extendLeaves();
				
		if(m_leaves.size()>m_maxUsedLeaves)	 m_maxUsedLeaves=m_leaves.size();
		// std::cout << m_currentKmerSize << ":" << m_currentLength << ":" << m_leaves.size() << "\n";	

		//see if terminating string is reached
		if(isTerminated(results))
			break;		
    }
		
	//find the path with maximum kmer coverage
	if( results.size()>0 )
	{
		//if multiple paths are bubbles collapsing all together at terminal
		if(results.size()==m_leaves.size()) 
		{
			// std::cout << m_maxUsedLeaves << "\t" << results.size() <<"\n" << results[0].thread << "\n";
			m_isBubbleCollapsed=true;
		}
		std::string tmpseq;
		for (size_t i = 0 ; i < results.size() ;i++){
			//bug fix: m_secondread may be shorter than m_minOverlap
			if(m_secondread.length()>m_minOverlap)
				tmpseq=results[i].thread+m_secondread.substr(m_minOverlap);
			else
				tmpseq=results[i].thread;
			
			size_t cov = calculateKmerCoverage (tmpseq, m_minOverlap, m_indices.pBWT);
			// size_t cov=results[i].SAICoverage;
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

	// if(m_leaves.size() > 512){
		// std::cout << m_leaves.size() << "\t" << m_currentLength<< "\t" << m_currentLength-m_pQuery->length()<< "\n";
		// printAll();
		// getchar();
	// }
	
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
void SAIntervalTree::printAll()
{
    std::cout << "Print all: \n";
    m_pRootNode->printAllStrings("");
}

// Extend each leaf node

void SAIntervalTree::attempToExtend(STNodePtrList &newLeaves)
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

void SAIntervalTree::extendLeaves()
{
    STNodePtrList newLeaves;
	
	//attempt to extend one base for each leave
    attempToExtend(newLeaves);
	
    //shrink the SAIntervals in case overlap is larger than read length, which lead to empty newLeaves
    if(!m_kmerMode  &&  newLeaves.empty() )
    {
        refineSAInterval(m_minOverlap);
        attempToExtend(newLeaves);
    }	
	
	//extension succeed
    if(!newLeaves.empty()){
		m_currentKmerSize++;
        m_currentLength++;  
	}

    m_leaves.clear();
    m_leaves = newLeaves;

	if(!m_leaves.empty() && (m_kmerMode || m_currentKmerSize >= m_maxOverlap) )
		refineSAInterval(m_minOverlap);

}

// Check for leaves whose extension has terminated. If the leaf has
// terminated, the walked string and coverage is pushed to the result vector
bool SAIntervalTree::isTerminated(SAIntervalNodeResultVector& results)
{
	bool found = false;

    for(STNodePtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
    {
        BWTInterval currfwd=(*iter)->fwdInterval;
        BWTInterval currrvc=(*iter)->rvcInterval;

        // assert(currfwd.isValid() || currrvc.isValid());

		//If terminating kmer is a substr, the current SA interval is a sub-interval of the terminating interval
        bool isFwdTerminated = currfwd.isValid() && currfwd.lower >= m_fwdTerminatedInterval.lower
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
            results.push_back(STresult);
            found =  true;
        }
    }

    return found;
}

bool SAIntervalTree::isTwoReadsOverlap(std::string & mergedseq)
{
    //case 1: 1st read sense overlap to 2nd read at exact m_minOverlap bases
    if(BWTInterval::equal(m_pRootNode->fwdInterval, m_fwdTerminatedInterval))
    {
        mergedseq= (*m_pQuery)+m_secondread.substr(m_minOverlap);
        return true;
    }

    //case 2: 1st read sense overlap 2nd read
    std::string secondLeftKmer=m_secondread.substr(0,m_minOverlap);
	//assume overlap can't exceed 100 bp
    size_t pos=m_pQuery->find(secondLeftKmer, m_pQuery->length()>=200?m_pQuery->length()-200:0);
    if(pos!=std::string::npos)	
    {
		//make sure entire suffix of 1st read after pos matches the prefix of 2nd read
		if( m_pQuery->substr(pos) == m_secondread.substr(0, m_pQuery->length()-pos) )
		{
			mergedseq=m_pQuery->substr(0,pos)+m_secondread;
			return true;
		}
    }

    //case 3: 1st read antisense overlap with 2nd read, or 1st read is substr of 2nd read
	//This is rare case and we don't do this in m_kmerMode during island joint
	if(m_kmerMode) return false;
    std::string firstLeftKmer=m_pQuery->substr(0,m_minOverlap);
    pos=m_secondread.find(firstLeftKmer);
	//assume antisense overlap can't exceed 50bp due to rare cases
    if(pos!=std::string::npos && pos <=50)
    {
        //make sure entire suffix of 2nd read after pos matches the prefix of 1st read
		if( m_secondread.substr(pos) ==  m_pQuery->substr(0, m_secondread.length()-pos))
		{
			//return overlapped portion
			mergedseq=m_secondread.substr(pos);
			return true;
		}
    }

    return false;

}

//update SA intervals of each leaf, which corresponds to one-base extension
std::vector<std::pair<std::string, BWTIntervalPair> > SAIntervalTree::getFMIndexExtensions(SAIntervalNode* pNode)
{
    std::vector<std::pair<std::string, BWTIntervalPair> > out;

    for(int i = 1; i < BWT_ALPHABET::size; ++i) //i=A,C,G,T
    {
        char b = BWT_ALPHABET::getChar(i);

        //update forward Interval using extension b
        BWTInterval fwdProbe=pNode->fwdInterval;
        if(fwdProbe.isValid())
            BWTAlgorithms::updateInterval(fwdProbe, b, m_indices.pRBWT);

        //update reverse complement Interval using extension rcb
        BWTInterval rvcProbe=pNode->rvcInterval;
		char rcb=BWT_ALPHABET::getChar(5-i); //T,G,C,A
        if(rvcProbe.isValid())
            BWTAlgorithms::updateInterval(rvcProbe, rcb, m_indices.pBWT);

        size_t bcount = 0;
        if(fwdProbe.isValid())
            bcount += fwdProbe.size();
        if(rvcProbe.isValid())
            bcount += rvcProbe.size();
			
		//min freq at fwd and rvc bwt
        if(bcount >= m_min_SA_threshold)
        {
			// if(bcount>50)
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

size_t SAIntervalTree::calculateKmerCoverage (const std::string & seq , size_t kmerLength , const BWT* pBWT)
{
	if (seq.length() < kmerLength) return 0;

	size_t cov = 0 ;
	for (size_t i=0; i<=seq.length()-kmerLength;i+=kmerLength/2)
		cov += BWTAlgorithms::countSequenceOccurrences(seq.substr(i,kmerLength) , pBWT );
	
	return cov;
}

// replace each kmer with highest one at each locus
bool SAIntervalTree::replaceLowFreqKmer (std::string & seq , size_t kmerLength)
{
	bool changed = false;
	
	for (size_t i=0; i <=seq.length()-kmerLength; i++)
	{
		//Forward kmer should be computed reversely using pRBWT for forward extension
		BWTInterval fwdProbe=BWTAlgorithms::findInterval(m_indices.pRBWT, reverse(seq.substr(i, kmerLength-1)));
		BWTInterval rvcProbe=BWTAlgorithms::findInterval(m_indices.pBWT, reverseComplement(seq.substr(i, kmerLength-1)));
		
		size_t maxcov=0;
		for(int j = 1; j < BWT_ALPHABET::size; ++j) //j=A,C,G,T
		{
			char b = BWT_ALPHABET::getChar(j);

			//update forward Interval using extension b
			if(fwdProbe.isValid())
				BWTAlgorithms::updateInterval(fwdProbe, b, m_indices.pRBWT);

			//update reverse complement Interval using extension rcb
			char rcb=BWT_ALPHABET::getChar(5-i); //T,G,C,A
			if(rvcProbe.isValid())
				BWTAlgorithms::updateInterval(rvcProbe, rcb, m_indices.pBWT);

			size_t bcount = 0;
			if(fwdProbe.isValid())
				bcount += fwdProbe.size();
			if(rvcProbe.isValid())
				bcount += rvcProbe.size();

			if(bcount > maxcov) {
				maxcov = bcount;
				seq.replace(i+kmerLength-1, 1, 1, b);
				changed = true;
			}
		}
	}
	
	return changed;
}

// Refine SA intervals of each leave with a new kmer
void SAIntervalTree::refineSAInterval(size_t newKmerSize)
{
	assert(m_currentLength >= newKmerSize);

    for(STNodePtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
    {
        // reset the SA intervals using original m_minOverlap
        std::string pkmer = (*iter)->getSuffix(newKmerSize);
		(*iter)->fwdInterval=BWTAlgorithms::findInterval(m_indices.pRBWT, reverse(pkmer));
		(*iter)->rvcInterval=BWTAlgorithms::findInterval(m_indices.pBWT, reverseComplement(pkmer));
    }

	m_currentKmerSize=newKmerSize;
}

/***Dead code***/

// Remove leaves with two or more same kmers
void SAIntervalTree::removeLeavesByRepeatKmer()
{
    STNodePtrList newLeaves;

    for(STNodePtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
    {
        /*
        GAGGCAGTTGAGGCAGTTGAGGCAGTTGAGGCAGTTGAGGCAGTTGAGGCAGTTGAGGCAGTTGAGGCAGTTGAGGCAGTTGAGGCAGTTGAGGCAGTTG
        TGGATTCCAGATTGTTCGAGGAGAATTTGGTGGAGCTACGCGGGATCGAACCGCGGACCTCTTGCATGCCATGCAAGCGCTCTCCCAGCTGAGCTATAACCCCTTGGATTCCAGATTGTTCGAGGAGAATTTGGTGGAGCTACGCGGGATCGAACCGCGGACCTCTTGCATGCCATGCAAGCGCTCTCCCAGCTGAGCTATAACC
        GAGAGGGACTCGAACCCTCACACCCGGGGGGCACTAACACCTGAAGCTAGCGCGTCTACCAATTCCGCCACCTTCGCACATCGGGTTATCAGTCTGGATTTACATGCTGTCTGATAAAAGCATGGTGCGAAGAGAGGGACTCGAACCCTCACACCCGGGGGGCACTAACACCTGAAGCTAGCGCGTCTACCAATTCCGCCACCTTCGCACATCGGGTTATCAGTCTGGATTT
        GCATATCCATCCCACCAGCACATCGACCTATCGACTCATCAGTTCATCAGTTCATCAGTTCATCAGTTCATCAGTTCATCAGTTCATCAGTTCATCAGTTCATCAGTTCATCAGTTCATCAGTTCATCAGTTCATCAGTTCATC
        CATCGGCGTCAGCCTGCTGGGCTTCACCCATCAGGGCAACAAGTGGCTGTGGCAGCAGGCCAGGGCCGCTCTTCCCTCCCTCAAGGGGGAGCTGGTGGCGGGGGGGGGGGGGGGGGGGGGGGGGGCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
        */
        std::string STNodeStr = (*iter)->getFullString();
        std::string fwdrepeatunit = STNodeStr.substr(STNodeStr.size()-m_minOverlap);
        std::string revrepeatunit = reverseComplement(fwdrepeatunit);
        size_t index1=STNodeStr.find(fwdrepeatunit);
        size_t index2=STNodeStr.find(revrepeatunit);

        if(index1 == (STNodeStr.size()- m_minOverlap) && index2 == std::string::npos)
        {
            newLeaves.push_back(*iter);
        }
    }

    m_leaves=newLeaves;
}
