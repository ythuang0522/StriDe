//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// Bidirectional graph
//
#include <assert.h>
#include <ostream>
#include <fstream>
#include <iostream>
#include <map>
#include "Bigraph.h"
#include "Timer.h"
#include "ASQG.h"

//
//
//
Bigraph::Bigraph() : m_hasContainment(false), m_hasTransitive(false), m_isExactMode(false), m_minOverlap(0), m_errorRate(0.0f)
{
	// Set up the memory pools for the graph
	m_pEdgeAllocator = new SimpleAllocator<Edge>();
	m_pVertexAllocator = new SimpleAllocator<Vertex>();

	m_vertices.set_deleted_key("");
	//WARN_ONCE("HARDCODED HASH TABLE MAX SIZE: 200000000");
	//m_vertices.resize(1000000);//produce unknown bug
}

//
//
//
Bigraph::~Bigraph()
{
	VertexPtrMap::iterator iter = m_vertices.begin();
	for(; iter != m_vertices.end(); ++iter)
	{
		delete iter->second;
		iter->second = NULL;
	}

	// Clean up the memory pools
	delete m_pEdgeAllocator;
	delete m_pVertexAllocator;
}

//
// Add a vertex
//
void Bigraph::addVertex(Vertex* pVert)
{
	std::pair<VertexPtrMapIter, bool> result =
	m_vertices.insert(std::make_pair(pVert->getID(), pVert));
	if(!result.second)
	{
		std::cerr << "Error: Attempted to insert vertex into graph with a duplicate id: " <<
		pVert->getID() << "\n";
		std::cerr << "All reads must have a unique identifier\n";
		exit(1);
	}
}

//
// Remove a vertex that is not connected to any other
//
void Bigraph::removeIslandVertex(Vertex* pVertex)
{
	assert(pVertex->countEdges() == 0);

	// Remove the vertex from the collection
	VertexID id = pVertex->getID();
	delete pVertex;
	m_vertices.erase(id);
}

//
// Remove a vertex
//
void Bigraph::removeConnectedVertex(Vertex* pVertex)
{
	// Remove the edges pointing to this Vertex
	pVertex->deleteEdges();

	// Remove the vertex from the collection
	VertexID id = pVertex->getID();
	delete pVertex;
	m_vertices.erase(id);
}


//
// Check for the existance of a vertex
//
bool Bigraph::hasVertex(VertexID id)
{
	VertexPtrMapIter iter = m_vertices.find(id);
	return iter != m_vertices.end();
}

//
// Get a vertex
//
Vertex* Bigraph::getVertex(VertexID id) const
{
	VertexPtrMapConstIter iter = m_vertices.find(id);
	if(iter == m_vertices.end())
	return NULL;
	return iter->second;
}

//
// Add an edge
//
void Bigraph::addEdge(Vertex* pVertex, Edge* pEdge)
{
	assert(pEdge->getStart() == pVertex);
	pVertex->addEdge(pEdge);
}

//
// Remove an edge
//
void Bigraph::removeEdge(const EdgeDesc& ed)
{
	ed.pVertex->removeEdge(ed);
}

//
// High level merge function that does not specify an edge
//
void Bigraph::mergeVertices(VertexID id1, VertexID id2)
{
	Vertex* pVert1 = getVertex(id1);

	// Get the edges from vertex1 to vertex2
	EdgePtrVec edgesTo = pVert1->findEdgesTo(id2);

	if(edgesTo.empty())
	{
		std::cerr << "mergeVertices: vertices are not connected\n";
		return;
	}

	if(edgesTo.size() > 1)
	{
		std::cerr << "mergeVertces: cannot merge because of ambigious edges\n";
		return;
	}

	// There is a single unique edge between the vertices
	Edge* mergeEdge = *edgesTo.begin();

	// Call the real merging function
	merge(pVert1, mergeEdge);
}

//
// Merge two vertices along the specified edge
//
void Bigraph::merge(Vertex* pV1, Edge* pEdge)
{

	Vertex* pV2 = pEdge->getEnd();
	//std::cout << "Merging " << pV1->getID() << " with " << pV2->getID() << "\n";

	// Merge the data
	pV1->merge(pEdge);

	// Get the twin edge (the edge in v2 that points to v1)
	Edge* pTwin = pEdge->getTwin();

	// Ensure v2 has the twin edge
	assert(pV2->hasEdge(pTwin));

	//
	assert(pV1->hasEdge(pEdge));
	size_t transLength = pV2->getOriginLength(!pTwin->getDir());
	pV1->setOriginLength(transLength,pEdge->getDir());


	// Get the edge set opposite of the twin edge (which will be the new edges in this direction for V1)
	EdgePtrVec transEdges = pV2->getEdges(!pTwin->getDir());

	// Move the edges from pV2 to pV1
	for(EdgePtrVecIter iter = transEdges.begin(); iter != transEdges.end(); ++iter)
	{
		Edge* pTransEdge = *iter;

		// Remove the edge from V2, this does not destroy the edge
		pV2->removeEdge(pTransEdge);

		// Join pEdge to the start of transEdge
		// This updates the starting point of pTransEdge to be V1
		// This calls Edge::extend on the twin edge
		pTransEdge->join(pEdge);
		assert(pTransEdge->getDir() == pEdge->getDir());
		pV1->addEdge(pTransEdge); // add to V1

		// Notify the edges they have been updated
		pTransEdge->update();
		pTransEdge->getTwin()->update();
	}

	// Remove the edge from pV1 to pV2
	pV1->removeEdge(pEdge);
	delete pEdge;
	pEdge = 0;

	// Remove the edge from pV2 to pV1
	pV2->removeEdge(pTwin);
	delete pTwin;
	pEdge = 0;

	// Remove V2
	// It is guarenteed to not be connected
	removeIslandVertex(pV2);
	//validate();
}


Vertex* Bigraph::buildVertex (EdgePtrVec & pathEdges,std::string vertexName)
{
	if (pathEdges.empty()) return NULL;
	std::cout << vertexName << std::endl;

	Vertex * pStartVertex = pathEdges.front()->getStart();
	Vertex * pEndVertex = pathEdges.back()->getEnd();
	assert (pStartVertex!=NULL && pEndVertex!=NULL);

	for (size_t i=1;i<pathEdges.size();i++)
	assert ( pathEdges[i-1]->getEnd() == pathEdges[i]->getStart() );


	EdgeDir extDir = pathEdges.front()->getDir();
	//get long Lable
	std::string longLable ;
	EdgeComp cmpStart = EC_SAME;
	for ( size_t i=0 ; i<pathEdges.size() ; i++ )
	{
		std::string label = pathEdges[i]->getLabel();
		if (cmpStart == EC_REVERSE) label = reverseComplement(label);
		if ( pathEdges[i]->getComp() == EC_REVERSE ) cmpStart = !cmpStart;

		if (extDir == ED_SENSE) longLable = longLable + label;
		else if (extDir == ED_ANTISENSE) longLable = label + longLable ;
		else assert(false);
	}

	//new a vertex
	std::string pathSeq = pStartVertex->getStr() ;
	//false : extend sense , true : extend antisense
	bool prepend = false;
	if (extDir == ED_SENSE)  pathSeq= pathSeq+longLable;
	else // ED_ANTISENSE
	{
		pathSeq= longLable + pathSeq ;
		prepend = true;
	}
	Vertex* pNewVertex = new(getVertexAllocator()) Vertex(vertexName, pathSeq);
	addVertex(pNewVertex);

	Edge * pFinalPathEdge = pathEdges.back();
	EdgeComp endOrientation = cmpStart;

	EdgePtrVec transEdges = pEndVertex->getEdges(!pFinalPathEdge->getTwin()->getDir());

	//ANTISENSE EXTEND
	if (prepend)
	{
		//add edges to ED_SENSE
		EdgePtrVec edges = pStartVertex->getEdges(ED_SENSE);
		for(EdgePtrVecIter iter = edges.begin(); iter != edges.end(); ++iter)
		{
			Edge* pEdge = *iter;
			Vertex* end = pEdge->getEnd();
			EdgeDir dir = pEdge->getDir();
			assert (dir==ED_SENSE);
			EdgeComp comp =pEdge->getComp();
			int offset = longLable.length() ;
			SeqCoord matchCoord (pEdge->getMatchCoord().interval.start+offset
			, pEdge->getMatchCoord().interval.end+offset
			, pathSeq.length()) ;
			Edge* pNewEdge = new (getEdgeAllocator())  Edge(end,dir,comp,matchCoord);
			Edge* pTwin = pEdge->getTwin();
			Edge* pNewTwin = new (getEdgeAllocator())  Edge(pNewVertex,pTwin->getDir(),pTwin->getComp(), pTwin->getMatchCoord());

			pNewEdge->setTwin(pNewTwin);
			pNewTwin->setTwin(pNewEdge);

			addEdge(pNewVertex, pNewEdge);
			addEdge(end,pNewTwin);

		}

		//add edges to ED_ANTISENSE
		for(EdgePtrVecIter iter = transEdges.begin(); iter != transEdges.end(); ++iter)
		{
			if (endOrientation == EC_SAME)
			{
				Edge* pEdge = *iter;
				Vertex* end = pEdge->getEnd();
				EdgeDir dir = pEdge->getDir();
				assert (dir == ED_ANTISENSE);
				EdgeComp comp =pEdge->getComp();
				SeqCoord matchCoord (pEdge->getMatchCoord().interval.start
				, pEdge->getMatchCoord().interval.end
				, pathSeq.length()) ;

				Edge* pNewEdge = new (getEdgeAllocator())  Edge(end,dir,comp,matchCoord);
				Edge* pTwin = pEdge->getTwin();
				Edge* pNewTwin = new (getEdgeAllocator())  Edge(pNewVertex,pTwin->getDir(),pTwin->getComp(), pTwin->getMatchCoord());

				pNewEdge->setTwin(pNewTwin);
				pNewTwin->setTwin(pNewEdge);
				addEdge(pNewVertex, pNewEdge);
				addEdge(end,pNewTwin);

			}

			else // EC_REVERSE
			{
				Edge* pEdge = *iter;
				Vertex* end = pEdge->getEnd();
				EdgeDir dir = !pEdge->getDir();
				assert (dir == ED_ANTISENSE);
				EdgeComp comp = !pEdge->getComp();
				SeqCoord m = pEdge->getMatchCoord();
				int matchLength = m.length() ;
				SeqCoord matchCoord (0,matchLength-1,pathSeq.length());
				Edge* pNewEdge = new (getEdgeAllocator())  Edge(end,dir,comp,matchCoord);
				Edge* pTwin = pEdge->getTwin();
				Edge* pNewTwin = new (getEdgeAllocator())  Edge(pNewVertex,pTwin->getDir(),!pTwin->getComp(), pTwin->getMatchCoord());

				pNewEdge->setTwin(pNewTwin);
				pNewTwin->setTwin(pNewEdge);
				addEdge(pNewVertex, pNewEdge);
				addEdge(end,pNewTwin);

			}
		}

		pNewVertex->validate();
	}
	else //SENSE EXTEND
	{
		//add edges to ED_ANTISENSE
		EdgePtrVec edges = pStartVertex->getEdges(ED_ANTISENSE);
		for(EdgePtrVecIter iter = edges.begin(); iter != edges.end(); ++iter)
		{
			Edge* pEdge = *iter;
			Vertex* end = pEdge->getEnd();
			EdgeDir dir = pEdge->getDir();
			EdgeComp comp =pEdge->getComp();
			assert (dir==ED_ANTISENSE);
			SeqCoord matchCoord (pEdge->getMatchCoord().interval.start, pEdge->getMatchCoord().interval.end, pathSeq.length()) ;

			Edge* pNewEdge = new (getEdgeAllocator())  Edge(end,dir,comp,matchCoord);
			Edge* pTwin = pEdge->getTwin();
			Edge* pNewTwin = new (getEdgeAllocator())  Edge(pNewVertex,pTwin->getDir(),pTwin->getComp(), pTwin->getMatchCoord());

			pNewEdge->setTwin(pNewTwin);
			pNewTwin->setTwin(pNewEdge);
			addEdge(pNewVertex, pNewEdge);
			addEdge(end,pNewTwin);

		}

		//add edges to ED_SENSE
		for(EdgePtrVecIter iter = transEdges.begin(); iter != transEdges.end(); ++iter)
		{
			if (endOrientation == EC_SAME)
			{
				Edge* pEdge = *iter;
				Vertex* end = pEdge->getEnd();
				EdgeDir dir = pEdge->getDir();
				assert (dir == ED_SENSE);
				EdgeComp comp =pEdge->getComp();
				int offset = longLable.length() ;
				SeqCoord matchCoord (pEdge->getMatchCoord().interval.start+offset
				, pEdge->getMatchCoord().interval.end+offset
				, pathSeq.length()) ;
				Edge* pNewEdge = new (getEdgeAllocator())  Edge(end,dir,comp,matchCoord);
				Edge* pTwin = pEdge->getTwin();
				Edge* pNewTwin = new (getEdgeAllocator())  Edge(pNewVertex,pTwin->getDir(),pTwin->getComp(), pTwin->getMatchCoord());

				pNewEdge->setTwin(pNewTwin);
				pNewTwin->setTwin(pNewEdge);
				addEdge(pNewVertex, pNewEdge);
				addEdge(end,pNewTwin);
			}

			else  //EC_REVERSE
			{
				Edge* pEdge = *iter;
				Vertex* end = pEdge->getEnd();
				EdgeDir dir = !pEdge->getDir();
				assert (dir == ED_SENSE);
				EdgeComp comp = !pEdge->getComp();
				SeqCoord m = pEdge->getMatchCoord();
				int matchLength = m.length() ;
				SeqCoord matchCoord (pathSeq.length()-matchLength,pathSeq.length()-1,pathSeq.length());
				Edge* pNewEdge = new (getEdgeAllocator())  Edge(end,dir,comp,matchCoord);
				Edge* pTwin = pEdge->getTwin();
				Edge* pNewTwin = new (getEdgeAllocator())  Edge(pNewVertex,pTwin->getDir(),!pTwin->getComp(), pTwin->getMatchCoord());

				pNewEdge->setTwin(pNewTwin);
				pNewTwin->setTwin(pNewEdge);
				addEdge(pNewVertex, pNewEdge);
				addEdge(end,pNewTwin);
			}
		}
	}
	pNewVertex->validate();


	return pNewVertex;

}

//
int Bigraph::sweepVertices(GraphColor c)
{
	int numRemoved = 0;
	VertexPtrMapIter iter = m_vertices.begin();
	while(iter != m_vertices.end())
	{
		VertexPtrMapIter next = iter;
		++next;
		if(iter->second->getColor() == c)
		{
			removeConnectedVertex(iter->second);
			++numRemoved;
		}
		iter = next;
	}
	return numRemoved;
}


//
int Bigraph::sweepEdges(GraphColor c)
{
	int numRemoved = 0;
	for(VertexPtrMapIter iter = m_vertices.begin(); iter != m_vertices.end(); ++iter)
	numRemoved += iter->second->sweepEdges(c);
	return numRemoved;
}

//    Simplify the graph by compacting singular edges
void Bigraph::simplify()
{
	assert(!hasContainment());
	size_t mergeCount = 0 ;
	
	//Linear time implementation by YTH
	VertexPtrMapIter iter = m_vertices.begin();

	while(iter != m_vertices.end())
	{
		mergeCount += simplify(iter->second, ED_SENSE);
		mergeCount += simplify(iter->second, ED_ANTISENSE);
		++iter;
	}
	
	if(mergeCount>0)
		std::cout << "<Simplify> Merge Vertices : "  <<  mergeCount << std::endl;
}

//merge unipaths from pV to pW in EdgeDir
size_t Bigraph::simplify(Vertex* pV, EdgeDir dir)
{
	size_t mergeCount = 0 ;
	
	// Get the edges for this direction
	EdgePtrVec edges = pV->getEdges(dir);
	
	// If there is a single edge in this direction, merge the vertices
	// Don't merge singular self edges though
	// Iteratively merge simple path pV->pW->pZ->...., where each iteration merge two vertices
	while(edges.size() == 1)
	{
		assert(!edges.front()->isSelf());
		Edge* pSingle = edges.front();
		
		// Check that the twin edge dir of pW is simple as well
		Edge* pTwin = pSingle->getTwin();
		Vertex* pW = pSingle->getEnd();
		if(pW->countEdges(pTwin->getDir()) == 1)
		{
			//statistics of overlap length and ratios, debug used only as the vector will become too large
			//m_simpleLenCount[pSingle->getMatchLength()] ++ ;
			//statsOverlapRatio(pV, pSingle); 
			
			//merge pW seq into pV and move its transitive edges to pV
			merge(pV, pSingle);
			mergeCount++ ;
			
			//remove self edges produced by V->W->V=V<->V
			edges = pV->getEdges(dir);
			size_t selfEdgeCount = 0 ;
			for(EdgePtrVecIter edge_iter = edges.begin(); edge_iter != edges.end(); ++edge_iter)
			{
				Edge* pEdge = *edge_iter;
				//self edge found
				if (pEdge->isSelf())
				{
					pV->deleteEdge(pEdge->getTwin());
					pV->deleteEdge(pEdge);
					selfEdgeCount++;
				}
			}
			//retrieve edges again if updated by selfedges
			if(selfEdgeCount >0)
				edges = pV->getEdges(dir);
		}//end of if pW single edge
		else
			break;
	}//end of if pV single edge

	
	return mergeCount ;
} 
 
 
void Bigraph::statsOverlapRatio(Vertex* pV1, Edge* pEdge)
{
	Vertex* pV2 = pEdge->getEnd();

	Edge* pTwin = pEdge->getTwin();
	assert(pV2->hasEdge(pTwin));

	size_t len1= pV1->getOriginLength(pEdge->getDir());
	size_t len2= pV2->getOriginLength(pTwin->getDir());

	size_t matchLen = pEdge->getMatchLength();
	size_t largeLen = 0 ,smallLen = 0;
	if (len1 > len2)
	{
		largeLen = len1;
		smallLen = len2;
	}
	else
	{
		largeLen = len2;
		smallLen = len1;
	}

	overlapRatioInfo oli ;
	oli.overlapLength=matchLen;
	oli.originalLength=smallLen;
	oli.ratioDiff=((double)matchLen/smallLen) - ((double)matchLen/largeLen) ;
	m_OverlapRatioStats.push_back(oli);

	//m_OverlapRatio.push_back( std::make_pair (pEdge->getMatchLength(),length) );

}

//
// Rename the vertices to have a sequential idx
// starting with prefix. This requires extra memory
// as we need to keep a vector of vertex pointers to reconstruct
// the graph after the vertices are renamed. This should
// only be done after the string graph has been simplified or
// else it could require a lot of memory.
//
void Bigraph::renameVertices(const std::string& prefix)
{
	size_t currIdx = 0;
	size_t numVertices = m_vertices.size();
	std::vector<Vertex*> vertexPtrVec(numVertices, 0);

	VertexPtrMapIter iter = m_vertices.begin();
	while(iter != m_vertices.end())
	{
		std::stringstream ss;
		ss << prefix << currIdx;
		iter->second->setID(ss.str());
		vertexPtrVec[currIdx] = iter->second;
		++iter;
		++currIdx;
	}

	// Clear the old graph
	m_vertices.clear();

	// Re-add the vertices
	for(size_t i = 0; i < numVertices; ++i)
	addVertex(vertexPtrVec[i]);
}

//
// Sort the adjacency list for each vertex by match length
//
void Bigraph::sortVertexAdjListsByMatchLen()
{
	VertexPtrMapIter iter = m_vertices.begin();
	for(; iter != m_vertices.end(); ++iter)
	iter->second->sortAdjListByMatchLen();
}


//
// Sort the adjacency list for each vertex by length
//
void Bigraph::sortVertexAdjListsByLen()
{
	VertexPtrMapIter iter = m_vertices.begin();
	for(; iter != m_vertices.end(); ++iter)
	iter->second->sortAdjListByLen();
}


//
// Sort the adjacency list for each vertex by ID
//
void Bigraph::sortVertexAdjListsByID()
{
	VertexPtrMapIter iter = m_vertices.begin();
	for(; iter != m_vertices.end(); ++iter)
	iter->second->sortAdjListByID();
}



//
// Validate that the graph is sane
//
void Bigraph::validate()
{
	VertexPtrMapIter iter = m_vertices.begin();
	for(; iter != m_vertices.end(); ++iter)
	{
		iter->second->validate();
	}
}

//
// Flip a vertex
//
void Bigraph::flip(VertexID /*id*/)
{
	assert(false);
#if 0
	// TODO: update this code
	Vertex* pVertex = getVertex(id);
	EdgePtrVec edges = pVertex->getEdges();

	for(EdgePtrVecIter iter = edges.begin(); iter != edges.end(); ++iter)
	{
		// Get the old twin
		GraphEdgeType twin = iter->getTwin();

		GraphEdgeType flipped = *iter;
		flipped.flip();

		// Remove the edge from the source ver
		pVertex->removeEdge(*iter);
		pVertex->addEdge(flipped);

		// Update the partner by deleting the old twin and
		Vertex* pV2 = getVertex(twin.getStart());
		pV2->removeEdge(twin);
		pV2->addEdge(flipped.getTwin());
	}
#endif
}

//
// Get the IDs of the vertices that do not branch (both sense/antisense degree <= 1)
//
VertexIDVec Bigraph::getNonBranchingVertices() const
{
	VertexIDVec out;
	VertexPtrMapConstIter iter = m_vertices.begin();
	for(; iter != m_vertices.end(); ++iter)
	{
		int senseEdges = iter->second->countEdges(ED_SENSE);
		int antisenseEdges = iter->second->countEdges(ED_ANTISENSE);
		if(antisenseEdges <= 1 && senseEdges <= 1)
		{
			out.push_back(iter->second->getID());
		}
	}
	return out;
}


//
// Get all the paths corresponding to the linear components of the graph
// Precondition: all vertices have in/out degree at most 1 (no branching)
//
PathVector Bigraph::getLinearComponents()
{
	PathVector outPaths;
	setColors(GC_WHITE);
	VertexPtrMapIter iter = m_vertices.begin();
	for(; iter != m_vertices.end(); ++iter)
	{
		// Output the linear path containing this vertex if it hasnt been visited already
		if(iter->second->getColor() != GC_BLACK)
		{
			outPaths.push_back(constructLinearPath(iter->second->getID()));
		}
	}
	assert(checkColors(GC_BLACK));
	return outPaths;
}

//
// Return all the path of nodes that can be linearally reached from this node
// The path expands in both directions so the first node in the path is not necessarily the source
//
Path Bigraph::constructLinearPath(VertexID id)
{
	Path sensePath;
	Path antisensePath;
	followLinear(id, ED_SENSE, sensePath);
	followLinear(id, ED_ANTISENSE, antisensePath);

	// Construct the final path
	Path final = reversePath(antisensePath);
	final.insert(final.end(), sensePath.begin(), sensePath.end());
	return final;
}

//
// Recursively follow the graph in the specified direction without branching
// outPath is an out-parameter of the edges that were followed
//
void Bigraph::followLinear(VertexID id, EdgeDir dir, Path& outPath)
{
	Vertex* pVertex = getVertex(id);
	EdgePtrVec edges = pVertex->getEdges(dir);

	// Color the vertex
	pVertex->setColor(GC_BLACK);

	if(edges.size() == 1)
	{
		Edge* pSingle = edges.front();
		outPath.push_back(pSingle);
		// Correct the direction for the comp
		assert(pSingle->getDir() == dir);
		EdgeDir corrected_dir = correctDir(pSingle->getDir(), pSingle->getComp());

		// Recurse
		followLinear(pSingle->getEndID(), corrected_dir, outPath);
	}
}

//
Path Bigraph::reversePath(const Path& path)
{
	Path out;
	for(Path::const_reverse_iterator iter = path.rbegin(); iter != path.rend(); ++iter)
	out.push_back((*iter)->getTwin());
	return out;
}

//
Vertex* Bigraph::getFirstVertex() const
{
	if(m_vertices.empty())
	return NULL;
	else
	return m_vertices.begin()->second;
}

// Returns a vector of pointers to the vertices
VertexPtrVec Bigraph::getAllVertices() const
{
	VertexPtrVec out;
	VertexPtrMapConstIter iter = m_vertices.begin();
	for(; iter != m_vertices.end(); ++iter)
	out.push_back(iter->second);
	return out;
}

//
// Append vertex sequences to the vector
void Bigraph::getVertexSequences(std::vector<std::string>& outSequences) const
{
	VertexPtrMapConstIter iter = m_vertices.begin();
	for(; iter != m_vertices.end(); ++iter)
	outSequences.push_back(iter->second->getSeq().toString());
}

//
// Set all the vertices and edges in the graph to the given color
//
void Bigraph::setColors(GraphColor c)
{
	VertexPtrMapIter iter = m_vertices.begin();
	for(; iter != m_vertices.end(); ++iter)
	{
		iter->second->setColor(c);
		iter->second->setEdgeColors(c);
	}
}

void Bigraph::setColorsP(GraphColor c)
{
	VertexPtrMapIter iter = m_vertices.begin();
	#pragma omp parallel private (iter)
	{
		for(; iter != m_vertices.end(); ++iter)
		{
			#pragma omp single nowait
			{
				iter->second->setColor(c);
				iter->second->setEdgeColors(c);
			}
		}
	}
}

//
// Set all the vertices in the graph to the given color
//
void Bigraph::setVerticesColor(GraphColor c)
{
	VertexPtrMapIter iter = m_vertices.begin();
	for(; iter != m_vertices.end(); ++iter) iter->second->setColor(c);
}

//
// Set simple vertices in the graph to the given color
//
void Bigraph::setStraightVerticesColor(GraphColor straightColor,GraphColor branchColor , GraphColor deadColor )
{
	assert (straightColor!=GC_WHITE && straightColor != GC_BLACK);
	VertexPtrMapIter iter = m_vertices.begin();
	int num_straight = 0 ;
	int num_dead = 0 ;
	int num_branch = 0 ;

	for(; iter != m_vertices.end(); ++iter)
	{
		Vertex * pVertex = iter->second;
		int s_count = pVertex->countEdges(ED_SENSE);
		int as_count = pVertex->countEdges(ED_ANTISENSE);
		if(s_count == 1 && as_count == 1) { pVertex->setColor(straightColor); ++num_straight;}
		else if (s_count == 0 || as_count == 0 )  { pVertex->setColor(deadColor); ++num_dead;}
		else { pVertex->setColor(branchColor); ++num_branch;}
	}
	std::cout << "# Straight: " << num_straight << " Dead: " << num_dead << " Branch: " << num_branch << std::endl;

}




//
// Get/Set the containment flag
//
void Bigraph::setContainmentFlag(bool b)
{
	m_hasContainment = b;
}

//
bool Bigraph::hasContainment() const
{
	return m_hasContainment;
}

//
// Get/Set the exact mode flag
//
void Bigraph::setExactMode(bool b)
{
	m_isExactMode = b;
}

//
bool Bigraph::isExactMode() const
{
	return m_isExactMode;
}

//
// Get/Set the transitive flag
//
void Bigraph::setTransitiveFlag(bool b)
{
	m_hasTransitive = b;
}

//
bool Bigraph::hasTransitive() const
{
	return m_hasTransitive;
}

//
//
//
void Bigraph::setMinOverlap(int mo)
{
	m_minOverlap = mo;
}

//
int Bigraph::getMinOverlap() const
{
	return m_minOverlap;
}

//
//
//
void Bigraph::setErrorRate(double er)
{
	m_errorRate = er;
}

//
double Bigraph::getErrorRate() const
{
	return m_errorRate;
}

//
// Check if all the vertices in the graph are the given color
//
bool Bigraph::checkColors(GraphColor c)
{
	VertexPtrMapIter iter = m_vertices.begin();
	for(; iter != m_vertices.end(); ++iter)
	{
		if(iter->second->getColor() != c)
		{
			std::cerr << "Warning vertex " << iter->second->getID() << " is color " << iter->second->getColor() << " expected " << c << "\n";
			return false;
		}
	}
	return true;
}


//
//Calculate assembly contiguity statistics.
//

void Bigraph::contigStats(size_t min_contig)  const
{
	std::vector<size_t> contigLengths;

	VertexPtrMapConstIter iter = m_vertices.begin();
	for(; iter != m_vertices.end(); ++iter)
	{
		if(iter->second->getSeqLen()>=min_contig)
			contigLengths.push_back (iter->second->getSeqLen());
	}

	std::sort(contigLengths.begin(),contigLengths.end());
	std::vector<size_t>::iterator it = contigLengths.begin();

	size_t sum = 0 ;
	for ( ; it !=  contigLengths.end(); ++it  )
		sum += *it;

	size_t curr =0 ;
	size_t n20 =0 ;
	size_t n50 =0 ;
	size_t n80 =0 ;

	std::vector<size_t>::reverse_iterator rit =contigLengths.rbegin();
	for ( ; rit!= contigLengths.rend(); ++rit)
	{
		curr+= *rit;
		if ( n20 ==0 && curr >= sum*0.2)  n20=*rit;
		if ( n50 ==0 && curr >= sum*0.5)  n50=*rit;
		if ( n80 ==0 && curr >= sum*0.8)  n80=*rit;
	}

	std::cout << "<Calculate assembly contiguity statistics>" << std::endl;
	std::cout << "Sum: " << sum  << "\tNum: " << contigLengths.size()  << "\tMAX: " << contigLengths.back() << "\tmin: " <<  contigLengths.front() << std::endl;
	std::cout << "N20: " << n20 << "\tN50: " << n50 << "\tN80: " << n80 << std::endl;

}





//
// Output simple stats
//
void Bigraph::stats() const
{
	int numVerts = 0;
	int numEdges = 0;

	VertexPtrMapConstIter iter = m_vertices.begin();
	for(; iter != m_vertices.end(); ++iter)
	{
		numEdges += iter->second->countEdges();
		++numVerts;
	}

	std::cout << "Graph has " << numVerts << " vertices and " << numEdges << " edges\n";
}

//
//
//
size_t Bigraph::getNumVertices() const
{
	return m_vertices.size();
}

//
// Output mem stats
//
void Bigraph::printMemSize() const
{
	size_t numVerts = 0;
	size_t vertMem = 0;

	size_t numEdges = 0;
	size_t edgeMem = 0;

	VertexPtrMapConstIter iter = m_vertices.begin();
	for(; iter != m_vertices.end(); ++iter)
	{
		++numVerts;
		vertMem += iter->second->getMemSize();

		EdgePtrVec edges = iter->second->getEdges();
		for(EdgePtrVecIter edgeIter = edges.begin(); edgeIter != edges.end(); ++edgeIter)
		{
			++numEdges;
			edgeMem += (*edgeIter)->getMemSize();
		}
	}
	printf("num verts: %zu using %zu bytes (%.2lf per vert)\n", numVerts, vertMem, double(vertMem) / numVerts);
	printf("num edges: %zu using %zu bytes (%.2lf per edge)\n", numEdges, edgeMem, double(edgeMem) / numEdges);
	printf("total: %zu\n", edgeMem + vertMem);
}

//
// Write the graph to a dot file
//
void Bigraph::writeDot(const std::string& filename, int dotFlags) const
{
	std::ofstream out(filename.c_str());

	std::string graphType = (dotFlags & DF_UNDIRECTED) ? "graph" : "digraph";

	out << graphType << " G\n{\n";
	VertexPtrMapConstIter iter = m_vertices.begin();
	for(; iter != m_vertices.end(); ++iter)
	{
		VertexID id = iter->second->getID();

		std::stringstream ss;
		ss << iter->second->getSeqLen();
		std::string len = ss.str();
		std::string label = (dotFlags & DF_NOID) ? "" : id+":"+len;

		out << "\"" << id << "\" [ label=\"" << label << "\" ";
		if(dotFlags & DF_COLORED)
		out << " style=\"filled\" fillcolor=\"" << getColorString(iter->second->getColor()) << "\" ";
		out << "];\n";
		iter->second->writeEdges(out, dotFlags);
	}
	out << "}\n";
	out.close();
}

//
// Write the graph to an ASQG file
//
void Bigraph::writeASQG(const std::string& filename) const
{
	std::ostream* pWriter = createWriter(filename);

	// Header
	ASQG::HeaderRecord headerRecord;
	headerRecord.setOverlapTag(m_minOverlap);
	headerRecord.setErrorRateTag(m_errorRate);
	headerRecord.setTransitiveTag(m_hasTransitive);
	headerRecord.setContainmentTag(m_hasContainment);
	headerRecord.write(*pWriter);


	VertexPtrMapConstIter iter;

	// Vertices
	for(iter = m_vertices.begin(); iter != m_vertices.end(); ++iter)
	{
		ASQG::VertexRecord vertexRecord(iter->second->getID(), iter->second->getSeq().toString());
		vertexRecord.write(*pWriter);
	}

	// Edges
	for(iter = m_vertices.begin(); iter != m_vertices.end(); ++iter)
	{
		EdgePtrVec edges = iter->second->getEdges();
		for(EdgePtrVecIter edgeIter = edges.begin(); edgeIter != edges.end(); ++edgeIter)
		{
			// We write one record for every bidirectional edge so only write edges
			// that are in canonical form (where id1 < id2)
			Overlap ovr = (*edgeIter)->getOverlap();
			if(ovr.id[0] <= ovr.id[1])
			{
				// Containment edges are in both directions so only output one
				// record if it is a containment
				if(!ovr.isContainment() || ((*edgeIter)->getDir() == ED_SENSE))
				{
					ASQG::EdgeRecord edgeRecord(ovr);
					edgeRecord.write(*pWriter);
				}
			}
		}
	}
	delete pWriter;
}

//
std::string Bigraph::getColorString(GraphColor c)
{
	switch(c)
	{
	case GC_WHITE:
		return "white";
	case GC_GRAY:
		return "gray";
	case GC_BLACK:
		return "black";
	case GC_BLUE:
		return "blue";
	case GC_RED:
		return "red";
	default:
		return "black";
	}
}

//Lock two vertices of the same edge with dead lock prevention
void Bigraph::lockTwoVertices(Vertex *pV1, Vertex* pV2)
{
	if(pV1 < pV2){
		pV1->setLock();
		pV2->setLock();
	}
	else
	{
		pV2->setLock();
		pV1->setLock();
	}

}
void Bigraph::unlockTwoVertices(Vertex *pV1, Vertex* pV2)
{
	pV1->releaseLock();
	pV2->releaseLock();
}