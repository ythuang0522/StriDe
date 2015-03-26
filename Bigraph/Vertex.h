//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// Vertex - Generic vertex class for bigraph
//
#ifndef VERTEX_H
#define VERTEX_H

// Includes
#include <stdio.h>
#include <map>
#include <set>
#include <vector>
#include <list>
#include <ostream>
#include <iostream>
#include <iterator>
#include "GraphCommon.h"
#include "QualityVector.h"
#include "EncodedString.h"
#include "SimpleAllocator.h"
#include "EdgeDesc.h"
#include "MultiOverlap.h"

#include <omp.h>

// Forward declare
class Edge;

// Default edge sorting function, by ID
struct EdgeIDComp
{
    bool operator()(const Edge* pA, const Edge* pB);
};

// Edge sorting function, by length
struct EdgeLenComp
{
    bool operator()(const Edge* pA, const Edge* pB);
};

// Edge sorting function, by matchlength
struct EdgeMatchLenComp
{
    bool operator()(const Edge* pA, const Edge* pB);
};

// Typedefs
typedef std::map<EdgeDesc, Edge*> EdgePtrMap;
typedef std::vector<Edge*> EdgePtrVec;
typedef std::set<EdgeDesc> EdgeDescSet;
typedef std::list<Edge*> EdgePtrList;
typedef EdgePtrMap::iterator EdgePtrMapIter;
typedef EdgePtrMap::const_iterator EdgePtrMapConstIter;
typedef EdgePtrVec::iterator EdgePtrVecIter;
typedef EdgePtrVec::const_iterator EdgePtrVecConstIter;
typedef EdgePtrList::iterator EdgePtrListIter;
typedef EdgePtrList::const_iterator EdgePtrListConstIter;

class Vertex
{
    public:

        Vertex(VertexID id, const std::string& s) : m_id(id),
                                                    m_seq(s),
                                                    m_color(GC_WHITE),
                                                    m_coverage(1),
                                                    m_isContained(false),
													m_isSuperRepeat(false)
													{
														//Added by CFL, store original read length, used by EdgeRemove Visitors
														m_originLength[ED_SENSE]=s.length();
														m_originLength[ED_ANTISENSE]=s.length();
														//Added by YTH, for parallel graph visitors
														omp_init_lock (&m_lock);
													}
        ~Vertex();

        // High-level modification functions

        // Merge another vertex into this vertex, as specified by pEdge
        void merge(Edge* pEdge);

        // For merging ARBRC: R and B will be merged into RBR
		void mergeTipVertex(Edge* pEdge);

        // sort the edges by the ID of the vertex they point to
        void sortAdjListByID();

        // sort the edges by the length of the label of the edge
        void sortAdjListByLen();

		//sort the edges by the match length of the label of the edge
		void sortAdjListByMatchLen();

        // Ensure that all the edges are unique
        bool markDuplicateEdges(GraphColor dupColor);

        // Get a multioverlap object representing the overlaps for this vertex
        MultiOverlap getMultiOverlap() const;

        // Edge list operations
        void addEdge(Edge* ep);
        void removeEdge(Edge* pEdge);
        void removeEdge(const EdgeDesc& ed);
        void deleteEdge(Edge* pEdge);
        void deleteEdges();

        int sweepEdges(GraphColor c);
        bool hasEdge(Edge* pEdge) const;
        bool hasEdge(const EdgeDesc& ed) const;
        bool hasEdgeTo(const Vertex* pY) const;

        Edge* getEdge(const EdgeDesc& ed);
        EdgePtrVec findEdgesTo(VertexID id);
        EdgePtrVec getEdges(EdgeDir dir) const;
        EdgePtrVec getEdges() const;
        EdgePtrVecIter findEdge(const EdgeDesc& ed);
        EdgePtrVecConstIter findEdge(const EdgeDesc& ed) const;
        Edge* getLongestOverlapEdge(EdgeDir dir) const;

        size_t countEdges() const;
        size_t countEdges(EdgeDir dir);

        // Calculate the difference in overlap lengths between
        // the longest and second longest edge
        int getOverlapLengthDiff(EdgeDir dir) const;

        // Ensure the vertex data is sane
        void validate() const;

        // setters
        void setID(VertexID id) { m_id = id; }
        void setEdgeColors(GraphColor c);
        void setSeq(const std::string& s) { m_seq = s; }
        void setColor(GraphColor c) { m_color = c; }
        void setContained(bool c) { m_isContained = c; }
        void setSuperRepeat(bool b) { m_isSuperRepeat = b; }
		void setOriginLength (size_t l,EdgeDir dir){ m_originLength[dir] = l;}

        // getters
        VertexID getID() const { return m_id; }
        GraphColor getColor() const { return m_color; }
        const DNAEncodedString& getSeq() const { return m_seq; }
        std::string getStr() const { return m_seq.toString(); }
        size_t getSeqLen() const { return m_seq.length(); }
        size_t getMemSize() const;
        bool isContained() const { return m_isContained; }
        bool isSuperRepeat() const { return m_isSuperRepeat; }
        uint16_t getCoverage() const { return m_coverage; }
		size_t  getOriginLength(EdgeDir dir) const { return m_originLength[dir]; }


        // Memory management
        void* operator new(size_t /*size*/, SimpleAllocator<Vertex>* pAllocator)
        {
            return pAllocator->alloc();
        }

        void operator delete(void* /*target*/, size_t /*size*/)
        {
            // delete does nothing since all allocations go through the memory pool
            // belonging to the graph. The memory allocated for the vertex will be
            // cleaned up when the graph is destroyed.
        }

        // Output edges in graphviz format
        void writeEdges(std::ostream& out, int dotFlags) const;

		// print VertexInformation
		void PrintVertexInformation()
		{
			std::cout <<  "VertexID : " << m_id <<  "\tDNAEncodedString : " << m_seq << "\tGraphColor : " << (int)m_color
				<<  "\tcoverage : " << m_coverage  <<  "\tm_isContained : "  << m_isContained << "\tm_isSuperRepeat : " << m_isSuperRepeat << "\n";
		}

		//Added by YTH, for parallel graph simplify() implementation
		void setLock()
		{
			omp_set_lock (&m_lock);
		}
		
		void releaseLock()
		{
			omp_unset_lock (&m_lock);
		}
		
		bool testLock()
		{
			return omp_test_lock (&m_lock);
		}
		
		//Read IDs in SA index mapped to this vertex, used by IslandJoin Visitor
		//[0,1,2,3]=AntisenseFwd, AntisenseRvc, SenseFwd, SenseRvc of pV
		std::vector< std::vector<int64_t> > pVReadIDs;
		
    private:

        // Global new is disallowed, all allocations must go through the pool
        void* operator new(size_t size)
        {
            return malloc(size);
        }

        // Ensure all the edges in DIR are unique
        bool markDuplicateEdges(EdgeDir dir, GraphColor dupColor);

        VertexID m_id;
        EdgePtrVec m_edges;
        DNAEncodedString m_seq;
        GraphColor m_color;

        // Counter of the number of vertices that have been merged into this one
        uint16_t m_coverage;

        bool m_isContained;
        bool m_isSuperRepeat;

		//Added by CFL, store original read length, used by EdgeRemove Visitors
		size_t m_originLength[ED_COUNT];
		
		//Added by YTH, for parallel visitor implementation
		omp_lock_t m_lock;
		
};

#endif
