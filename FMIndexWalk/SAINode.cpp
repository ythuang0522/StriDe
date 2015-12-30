///----------------------------------------------
// Copyright 2015 National Chung Cheng University
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
#include "SAINode.h"
#include "BWTAlgorithms.h"

//
// SAIntervalNode
//
SAINode::SAINode(const std::string* pQuery, SAINode* parent) : 
									   m_kmerCount(0), m_pQuery(pQuery),m_pParent(parent)
{

}

// Destructor, recurisvely delete the children of the node
SAINode::~SAINode()
{
    // Delete children
    for(SAINodePtrList::iterator iter = m_children.begin(); iter != m_children.end(); ++iter)
        delete *iter;

}

// Return a suffix of length l of the path from the root to this node
std::string SAINode::getSuffix(size_t l) const
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
std::string SAINode::getFullString() const
{
    if(m_pParent == NULL)
        return m_label;
    else
        return m_pParent->getFullString() + m_label;
}

// Create a new child node with the given label. Returns a pointer to the new node.
SAINode* SAINode::createChild(const std::string& label)
{
    SAINode* pAdded = new SAINode(m_pQuery, this);
    m_children.push_back(pAdded);
    pAdded->extend(label);

    return pAdded;
}

// Extend the label of this node
void SAINode::extend(const std::string& ext)
{
    assert(!ext.empty());
    m_label.append(ext);
}


void SAINode::computeInitial(const std::string& initialLabel)
{
    m_label = initialLabel;

}


// Print the string(s) represented by this node and its children
void SAINode::printAllStrings(const std::string& parent) const
{
    if(m_children.empty())
    {
        std::cout << ">\n" << parent + m_label << "\n";
    }
    else
    {
        for(SAINodePtrList::const_iterator iter = m_children.begin(); iter != m_children.end(); ++iter)
            (*iter)->printAllStrings(parent + m_label);
    }
}

// Create a new child node with the given label. Returns a pointer to the new node.
SAIntervalNode* SAIntervalNode::createChild(const std::string& label)
{
    SAIntervalNode* pAdded = new SAIntervalNode(m_pQuery, this);
    m_children.push_back(pAdded);
    pAdded->extend(label);

	// copy each member value to child
	pAdded->addKmerCount(this->getKmerCount());
	
    return pAdded;
}

// Create a new child node with the given label. Returns a pointer to the new node.
SAIOverlapNode* SAIOverlapNode::createChild(const std::string& label)
{
    SAIOverlapNode* pAdded = new SAIOverlapNode(m_pQuery, this);
    pAdded->extend(label);
	
	// still lack of a copy constructor
	pAdded->lastSeedIdx = this->lastSeedIdx;
	pAdded->lastOverlapLen = this->lastOverlapLen;
	pAdded->totalSeeds = this->totalSeeds;
	pAdded->currOverlapLen = this->currOverlapLen;
	pAdded->queryOverlapLen = this->queryOverlapLen;
	pAdded->numOfErrors = this->numOfErrors;
	pAdded->lastSeedIdxOffset = this->lastSeedIdxOffset;
	pAdded->initSeedIdx = this->initSeedIdx;
	
	m_children.push_back(pAdded);
    
    return pAdded;
}
