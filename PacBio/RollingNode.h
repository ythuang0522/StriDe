//----------------------------------------------
// Copyright 2015 National Chung Cheng University
// Written by Yao-Ting Huang
// Released under the GPL
//-----------------------------------------------

// 
// The search tree represents a traversal through implicit FM-index graph
//
#ifndef ROLLINGODE_H
#define ROLLINGODE_H

#include <list>
#include "BWT.h"
#include "SAINode.h"
#include "cyclichash.h"

// Simulate circular queue for a kmer containing head and tail chars
class CircularKmer{
    private:
        char *item;
        int head;
        int tail;
		int MAX_SIZE;
    public:
        CircularKmer(std::string kmer);
		//need copy constructor
		
		CircularKmer(const CircularKmer& c)
		{
			head = c.head;
			tail = c.tail;
			MAX_SIZE = c.MAX_SIZE;
			item = new char[MAX_SIZE];
			strncpy(c.item, item, MAX_SIZE);
		}
		~CircularKmer(){
			delete item;
		}
        void shiftRight(char c);
		void shiftLeft(char c);
        inline char getFirstChar(){return item[head];};
        inline char getLastChar(){return item[tail];};
};
 

//
// RollingNode for implementation of FM-index extension using rolling hash
//
class RollingNode : public SAINode
{
    public:
        RollingNode(const std::string* pQuery, RollingNode* parent):SAINode(pQuery,parent)
		{
			currFwdKmer=currRvcKmer=NULL;
		}
		
        ~RollingNode()
		{
			if(currFwdKmer!=NULL)
				delete currFwdKmer;
			
			if(currRvcKmer!=NULL)
				delete currRvcKmer;
		};

		// Add a child node to this node with the given label
        // Returns a pointer to the created node
        RollingNode* createChild(const std::string& label);

        CircularKmer* currFwdKmer;
		CircularKmer* currRvcKmer;
		size_t currFwdHashValue;
		size_t currRvcHashValue;
};

// leaves of SAIOverlapNode
typedef std::list<RollingNode*> RollingNodePtrList;

#endif
