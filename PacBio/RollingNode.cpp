///----------------------------------------------
// Copyright 2015 National Chung Cheng University
// Written by Yao-Ting Huang
// Released under the GPL
//-----------------------------------------------
//
// RollingNode - Iteratively construct a
// string representing a walk through an assembly graph
// matching a query sequence.
//
// The assembly graph is abstractly represented as
// an FM-index.
//
#include "RollingNode.h"
#include "BWTAlgorithms.h"

CircularKmer::CircularKmer(std::string kmer){
	MAX_SIZE = kmer.length();	
	item = new char[MAX_SIZE];
	strncpy(item, kmer.c_str(), MAX_SIZE);
	
	head = 0;
    tail = MAX_SIZE-1;
}
 
void CircularKmer::shiftRight(char c)
{
	head = (head+1) % MAX_SIZE;
	tail = (tail+1) % MAX_SIZE;
	item[tail] = c;
	
}

void CircularKmer::shiftLeft(char c)
{
	tail = tail-1;
	if(tail < 0) tail = MAX_SIZE-1;
	
	head = head-1;
	if(head < 0) head = MAX_SIZE-1;
	item[head] = c;
	
}

// Create a new child node with the given label. Returns a pointer to the new node.
RollingNode* RollingNode::createChild(const std::string& label)
{
    RollingNode* pAdded = new RollingNode(m_pQuery, this);
    pAdded->extend(label);	
	
	pAdded->currFwdKmer = new CircularKmer(*(this->currFwdKmer));
	pAdded->currRvcKmer = new CircularKmer(*(this->currRvcKmer));
	pAdded->currFwdHashValue = this->currFwdHashValue;
	pAdded->currRvcHashValue = this->currRvcHashValue;
	
	m_children.push_back(pAdded);
    
    return pAdded;
}
