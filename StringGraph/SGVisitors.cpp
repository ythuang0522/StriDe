//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// SGVisitors - Algorithms that visit
// each vertex in the graph and perform some
// operation
//
#include "SGVisitors.h"
#include "SGSearch.h"
#include "stdaln.h"
#include "BWTAlgorithms.h"
#include <iomanip>
#include "SAIntervalTree.h"
//
// SGFastaVisitor - output the vertices in the graph in
// fasta format
//
bool SGFastaVisitor::visit(StringGraph* /*pGraph*/, Vertex* pVertex)
{
	size_t seqLen = pVertex->getSeqLen() ;

	m_fileHandle << ">" << pVertex->getID() << " " << seqLen;

	if (pBWT!=NULL && kmerLength!= 0)
	{
		size_t cov = 0 ;
		std::string contigs=pVertex->getStr();

		#pragma omp parallel for
		for (size_t i=0 ; i <= seqLen-kmerLength ; i++)
		{
			std::string kmer = contigs.substr (i,kmerLength);

			size_t countSame = BWTAlgorithms::countSequenceOccurrencesSingleStrand (kmer,pBWT);
			size_t countRevc = BWTAlgorithms::countSequenceOccurrencesSingleStrand (reverseComplement(kmer),pBWT);
			size_t count = countSame > countRevc ? countSame : countRevc ;


			#pragma omp atomic
			cov += count;
			//cov += BWTAlgorithms::countSequenceOccurrences(kmer,pBWT);
		}

		m_fileHandle << " " << cov << " "<< (float)cov/ (seqLen-kmerLength+1) <<"\n";
	}

	else
	m_fileHandle << " " << pVertex->getCoverage() << " " <<
				pVertex->getOriginLength(ED_ANTISENSE) << " " <<
				pVertex->getOriginLength(ED_SENSE) <<"\n";

	m_fileHandle << pVertex->getSeq() << "\n";
	return false;
}


//
// SGTransRedVisitor - Perform a transitive reduction about this vertex
// This uses Myers' algorithm (2005, The fragment assembly string graph)
// Precondition: the edge list is sorted by length (ascending)
void SGTransitiveReductionVisitor::previsit(StringGraph* pGraph)
{
	// The graph must not have containments
	assert(!pGraph->hasContainment());

	// Set all the vertices in the graph to "vacant"
	pGraph->setColors(GC_WHITE);
	pGraph->sortVertexAdjListsByLen();

	marked_verts = 0;
	marked_edges = 0;
}

//
bool SGTransitiveReductionVisitor::visit(StringGraph* /*pGraph*/, Vertex* pVertex)
{
	size_t trans_count = 0;
	static const size_t FUZZ = 10; // see myers

	for(size_t idx = 0; idx < ED_COUNT; idx++)
	{
		EdgeDir dir = EDGE_DIRECTIONS[idx];
		EdgePtrVec edges = pVertex->getEdges(dir); // These edges are already sorted
		if(edges.size() == 0)
		continue;

		for(size_t i = 0; i < edges.size(); ++i)
		(edges[i])->getEnd()->setColor(GC_GRAY);

		Edge* pLongestEdge = edges.back();
		size_t longestLen = pLongestEdge->getSeqLen() + FUZZ;

		// Stage 1
		for(size_t i = 0; i < edges.size(); ++i)
		{
			Edge* pVWEdge = edges[i];
			Vertex* pWVert = pVWEdge->getEnd();

			EdgeDir transDir = !pVWEdge->getTwinDir();
			if(pWVert->getColor() == GC_GRAY)
			{
				EdgePtrVec w_edges = pWVert->getEdges(transDir);
				for(size_t j = 0; j < w_edges.size(); ++j)
				{
					Edge* pWXEdge = w_edges[j];
					size_t trans_len = pVWEdge->getSeqLen() + pWXEdge->getSeqLen();
					if(trans_len <= longestLen)
					{
						if(pWXEdge->getEnd()->getColor() == GC_GRAY)
						{
							// X is the endpoint of an edge of V, therefore it is transitive
							pWXEdge->getEnd()->setColor(GC_BLACK);
						}
					}
					else
					break;
				}
			}
		}

		// Stage 2
		for(size_t i = 0; i < edges.size(); ++i)
		{
			Edge* pVWEdge = edges[i];
			Vertex* pWVert = pVWEdge->getEnd();

			EdgeDir transDir = !pVWEdge->getTwinDir();
			EdgePtrVec w_edges = pWVert->getEdges(transDir);
			for(size_t j = 0; j < w_edges.size(); ++j)
			{
				Edge* pWXEdge = w_edges[j];
				size_t len = pWXEdge->getSeqLen();

				if(len < FUZZ || j == 0)
				{
					if(pWXEdge->getEnd()->getColor() == GC_GRAY)
					{
						// X is the endpoint of an edge of V, therefore it is transitive
						pWXEdge->getEnd()->setColor(GC_BLACK);
					}
				}
				else
				{
					break;
				}
			}
		}

		for(size_t i = 0; i < edges.size(); ++i)
		{
			if(edges[i]->getEnd()->getColor() == GC_BLACK)
			{
				// Mark the edge and its twin for removal
				if(edges[i]->getColor() != GC_BLACK || edges[i]->getTwin()->getColor() != GC_BLACK)
				{
					edges[i]->setColor(GC_BLACK);
					edges[i]->getTwin()->setColor(GC_BLACK);
					marked_edges += 2;
					trans_count++;
				}
			}
			edges[i]->getEnd()->setColor(GC_WHITE);
		}
	}

	if(trans_count > 0)
	++marked_verts;

	return false;
}

// Remove all the marked edges
void SGTransitiveReductionVisitor::postvisit(StringGraph* pGraph)
{
	//printf("TR marked %d verts and %d edges\n", marked_verts, marked_edges);

	std::cout << "Remove " << pGraph->sweepEdges(GC_BLACK)/2 << " transitive edges." << std::endl;
	pGraph->setTransitiveFlag(false);
	assert(pGraph->checkColors(GC_WHITE));
}


//
// SGContainRemoveVisitor - Removes contained
// vertices from the graph
//
void SGContainRemoveVisitor::previsit(StringGraph* pGraph)
{
	//pGraph->setColors(GC_WHITE);
	//pGraph->setVerticesColor(GC_WHITE);
	// Clear the containment flag, if any containments are added
	// during this algorithm the flag will be reset and another
	// round must be re-run
	pGraph->setContainmentFlag(false);

	//pWriter = createWriter("rm_contain.list");
}

//
bool SGContainRemoveVisitor::visit(StringGraph* pGraph, Vertex* pVertex)
{
	if(!pVertex->isContained())
	return false;
	// Add any new irreducible edges that exist when pToRemove is deleted
	// from the graph
	EdgePtrVec neighborEdges = pVertex->getEdges();

	// If the graph has been transitively reduced, we have to check all
	// the neighbors to see if any new edges need to be added. If the graph is a
	// complete overlap graph we can just remove the edges to the deletion vertex
	if(!pGraph->hasTransitive() && !pGraph->isExactMode())
	{
		// This must be done in order of edge length or some transitive edges
		// may be created
		EdgeLenComp comp;
		std::sort(neighborEdges.begin(), neighborEdges.end(), comp);

		for(size_t j = 0; j < neighborEdges.size(); ++j)
		{
			Vertex* pRemodelVert = neighborEdges[j]->getEnd();
			Edge* pRemodelEdge = neighborEdges[j]->getTwin();
			SGAlgorithms::remodelVertexForExcision(pGraph,
			pRemodelVert,
			pRemodelEdge);
		}
	}

	// Delete the edges from the graph
	for(size_t j = 0; j < neighborEdges.size(); ++j)
	{
		Vertex* pRemodelVert = neighborEdges[j]->getEnd();
		Edge* pRemodelEdge = neighborEdges[j]->getTwin();
		pRemodelVert->deleteEdge(pRemodelEdge);
		pVertex->deleteEdge(neighborEdges[j]);
	}
	pVertex->setColor(GC_BLACK);

	//*pWriter << pVertex->getID() << std::endl;


	return false;
}

void SGContainRemoveVisitor::postvisit(StringGraph* pGraph)
{
	pGraph->sweepVertices(GC_BLACK);
	//delete pWriter;
}


//
// SGTrimVisitor - Remove "dead-end" vertices from the graph
//
void SGTrimVisitor::previsit(StringGraph* /*pGraph*/)
{
	num_island = 0;
	num_terminal = 0;
	//pGraph->setColors(GC_WHITE);

	if(!m_filename.empty())
		std::ofstream m_fileHandle(m_filename.c_str());
	
}

// Mark any nodes that either dont have edges or edges in only one direction for removal
bool SGTrimVisitor::visit(StringGraph* /*pGraph*/, Vertex* pVertex)
{
	if(pVertex->countEdges() == 0)
	{
		// Is an island, remove if the sequence length is less than the threshold
		if(pVertex->getSeqLen() < m_minLength)
		{
			pVertex->setColor(GC_BLACK);
			++num_island;
			if(!m_filename.empty())
				m_fileHandle << ">" << pVertex->getID() << " " << pVertex->getSeqLen() << " island (threshold:" << m_minLength << ")\n"<< pVertex->getSeq() << "\n";
		}
	}
	else
	{
		// Check if this node is a dead-end
		for(size_t idx = 0; idx < ED_COUNT; idx++)
		{
			EdgeDir dir = EDGE_DIRECTIONS[idx];
			if(pVertex->countEdges(dir) == 0 && pVertex->getSeqLen() < m_minLength)
			{
				pVertex->setColor(GC_BLACK);
				++num_terminal;

				if(!m_filename.empty())
				{
					m_fileHandle << ">" << pVertex->getID() << " " << pVertex->getSeqLen()
					<< " terminal (threshold:" << m_minLength << ") Direction:" << dir << " ";

					EdgePtrVec edges = pVertex->getEdges(!dir);
					assert (edges.size()!=0);
					size_t i = 0;
					for (i = 0 ;i < edges.size()-1 ; i++ )
						m_fileHandle << edges[i]->getEnd()->getID() << " , ";
					m_fileHandle << edges[i]->getEnd()->getID() << "\n";

					m_fileHandle << pVertex->getSeq() << "\n";
				}

				return true;
			}
		}
	}

	return false;
}

// Remove all the marked edges
void SGTrimVisitor::postvisit(StringGraph* pGraph)
{
	pGraph->sweepVertices(GC_BLACK);
	printf("StringGraphTrim: Removed %d island and %d dead-end short vertices\n", num_island, num_terminal);
	if(!m_filename.empty())
		m_fileHandle.close();
}

//
// SGDuplicateVisitor - Detect and remove duplicate edges
//
void SGDuplicateVisitor::previsit(StringGraph* pGraph)
{
	assert(pGraph->checkColors(GC_WHITE));
	(void)pGraph;
	m_hasDuplicate = false;
}

bool SGDuplicateVisitor::visit(StringGraph* /*pGraph*/, Vertex* pVertex)
{
	m_hasDuplicate = pVertex->markDuplicateEdges(GC_RED) || m_hasDuplicate;
	return false;
}


void SGDuplicateVisitor::postvisit(StringGraph* pGraph)
{
	//assert(pGraph->checkColors(GC_WHITE));
	if(m_hasDuplicate)
	{
		int numRemoved = pGraph->sweepEdges(GC_RED);
		if(!m_bSilent)
			std::cerr << "Warning: removed " << numRemoved << " duplicate edges\n";
	}
}



//
// SGSuperRepeatVisitor
//

// Remove all edges of nodes that have been marked as super repeats
void SGSuperRepeatVisitor::previsit(StringGraph*)
{
	m_num_superrepeats = 0;
}

//
bool SGSuperRepeatVisitor::visit(StringGraph*, Vertex* pVertex)
{
	if(pVertex->isSuperRepeat())
	{
		pVertex->deleteEdges();
		m_num_superrepeats += 1;
		return true;
	}
	return false;
}

//
void SGSuperRepeatVisitor::postvisit(StringGraph*)
{
	printf("Deleted edges for %zu super repetitive vertices\n", m_num_superrepeats);
}

//
// SGSmoothingVisitor - Find branches in the graph
// which arise from variation and remove them
//
void SGSmoothingVisitor::previsit(StringGraph* pGraph)
{
	pGraph->setVerticesColor(GC_WHITE);
	m_simpleBubblesRemoved = 0;
	m_complexBubblesRemoved = 0;
}

//
bool SGSmoothingVisitor::visit(StringGraph* pGraph, Vertex* pVertex)
{
	(void)pGraph;
	if(pVertex->getColor() == GC_RED)
	return false;

	bool found = false;
	for(size_t idx = 0; idx < ED_COUNT; idx++)
	{		
		EdgeDir dir = EDGE_DIRECTIONS[idx];
		EdgePtrVec edges = pVertex->getEdges(dir);
		if(edges.size() <= 1)
			continue;

		for(size_t i = 0; i < edges.size(); ++i)
		{
			if(edges[i]->getEnd()->getColor() == GC_RED)
			return false;
		}

		// The walks grow exponentially within repeats
		const int MAX_WALKS = 240;
		// The distance is the extended seq length excluding the starting vertex,
		const int MAX_DISTANCE = 2400;
		bool bIsDegenerate=false;
		bool bFailIndelSizeCheck=false;

		SGWalkVector variantWalks;
		SGSearch::findVariantWalks(pVertex, dir, MAX_DISTANCE, MAX_WALKS, variantWalks);

		if(variantWalks.size() > 0)
		{
			size_t selectedIdx = 0;
			size_t selectedCoverage = 0;
			
			// Calculate the minimum amount overlapped on the start/end vertex.
			// This is used to properly extract the sequences from walks that represent the variation.
			int minOverlapX = std::numeric_limits<int>::max();
			int minOverlapY = std::numeric_limits<int>::max();

			// Calculate walk coverage and min overlap
			for(size_t i = 0; i < variantWalks.size(); ++i)
			{
				if(variantWalks[i].getNumEdges() <= 1){
					bIsDegenerate = true;
					break;
				}

				// Calculate the walk coverage using the internal vertices of the walk.
				// The walk with the highest coverage will be retained
				// Walk coverage is the number of contracted vertices along the walk
				// The more contracted vertices implies larger overlap between two vertices.
				size_t walkCoverage = 0;
				for(size_t j = 1; j < variantWalks[i].getNumVertices() - 1; ++j)
					walkCoverage += variantWalks[i].getVertex(j)->getCoverage();
					
				if( walkCoverage> selectedCoverage || selectedCoverage == 0)
				{
					selectedIdx = i;
					selectedCoverage = walkCoverage;
				}

				// Calculate min overlap
				Edge* pFirstEdge = variantWalks[i].getFirstEdge();
				Edge* pLastEdge = variantWalks[i].getLastEdge();

				if((int)pFirstEdge->getMatchLength() < minOverlapX)
					minOverlapX = pFirstEdge->getMatchLength();

				if((int)pLastEdge->getTwin()->getMatchLength() < minOverlapY)
					minOverlapY = pLastEdge->getTwin()->getMatchLength();
			}
			
			if(bIsDegenerate)	
				continue;
			
			// Now check the length diff of each walk against the selected walk instead of alignment
			size_t selectedWalkLength=variantWalks[selectedIdx].getStartToEndDistance();
			for(size_t i = 0; i < variantWalks.size(); ++i)
			{
				double gapDivergence = abs( (int) variantWalks[i].getStartToEndDistance() - (int)selectedWalkLength) ;
				if ( m_bIsGapPrecent && gapDivergence/selectedWalkLength < 0.05) continue;  
				if ( gapDivergence > m_maxIndelLength )
				{
					// std::cout << variantWalks.size() << " \t" << variantWalks[i].getStartVertex()->getSeqLen() << "\t" << variantWalks[i].getLastVertex()->getSeqLen() << "\t" <<
						// variantWalks[i].getStartToEndDistance() << "\t" << selectedWalkLength << "\t" << 
						// abs( (int) variantWalks[i].getStartToEndDistance() - (int)selectedWalkLength)  << "\n";
	
					bFailIndelSizeCheck=true;
					break;
				}
			}
			
			if( bFailIndelSizeCheck ) continue;

			assert(selectedIdx != (size_t)-1);
			SGWalk& selectedWalk = variantWalks[selectedIdx];
			assert(selectedWalk.isIndexed());


			// The vertex set for each walk is not necessarily disjoint,
			// the selected walk may contain vertices that are part
			// of other paths, which should be preserved.
			for(size_t i = 0; i < variantWalks.size(); ++i)
			{
				if(i == selectedIdx)
					continue;

				SGWalk& currWalk = variantWalks[i];
				for(size_t j = 0; j < currWalk.getNumEdges() - 1; ++j)
				{
					Edge* currEdge = currWalk.getEdge(j);

					// If the vertex is also on the selected path, do not mark it
					Vertex* currVertex = currEdge->getEnd();
					if(!selectedWalk.containsVertex(currVertex->getID()))
					{
						currEdge->getEnd()->setColor(GC_RED);
						found = true;
					}
				}
			}

			if(variantWalks.size() == 2)
				m_simpleBubblesRemoved += 1;
			else
				m_complexBubblesRemoved += 1;

				++m_numRemovedTotal;
		}
	}
	return found;
}

// Remove all the marked edges
void SGSmoothingVisitor::postvisit(StringGraph* pGraph)
{
	pGraph->sweepVertices(GC_RED);
	//assert(pGraph->checkColors(GC_WHITE));

	printf("VariationSmoother: Removed %d simple and %d complex bubbles\n", m_simpleBubblesRemoved, m_complexBubblesRemoved);
}

//
// SGGraphStatsVisitor - Collect summary stasitics
// about the graph
//
void SGGraphStatsVisitor::previsit(StringGraph* /*pGraph*/)
{
	num_straight = 0 ;
	num_terminal = 0;
	num_island = 0;
	num_monobranch = 0;
	num_dibranch = 0;
	num_simple = 0;
	num_edges = 0;
	num_vertex = 0;
	sum_edgeLen = 0;
}

// Find bubbles (nodes where there is a split and then immediate rejoin) and mark them for removal
bool SGGraphStatsVisitor::visit(StringGraph* /*pGraph*/, Vertex* pVertex)
{
	int s_count = pVertex->countEdges(ED_SENSE);
	int as_count = pVertex->countEdges(ED_ANTISENSE);
	if(s_count == 0 && as_count == 0)
	{
		++num_island;
	}
	else if(s_count == 0 || as_count == 0)
	{
		++num_terminal;
	}

	if(s_count > 1 && as_count > 1)
		++num_dibranch;
	else if(s_count > 1 || as_count > 1)
		++num_monobranch;

	if(s_count == 1 || as_count == 1)
		++num_simple;

	if(s_count == 1 && as_count == 1)
		++num_straight;

	num_edges += (s_count + as_count);

	++num_vertex;

	return false;
}

//
void SGGraphStatsVisitor::postvisit(StringGraph* /*pGraph*/)
{
	printf("Vertices: %d Edges: %d Islands: %d Tips: %d Monobranch: %d Dibranch: %d Simple: %d Straight: %d\n",
	num_vertex, num_edges,num_island, num_terminal,num_monobranch, num_dibranch, num_simple,num_straight);
}






/*****************************************************************/
/*                       New Visitors   	*/
/*             written by CCU CSIE LAB 405 Tatsuki and Yao-Ting     */
/*****************************************************************/

//trim island/tip contigs terminated due to errors
bool SGFastaErosionVisitor::visit(StringGraph* /*pGraph*/, Vertex* pVertex)
{
	size_t seqLen = pVertex->getSeqLen() ;
    std::string contigs=pVertex->getStr();
    int  start = 0 , end = seqLen;

    if (pVertex->countEdges(ED_ANTISENSE)== 0)
    {
        for (int i=0 ; i <= (int)(seqLen-kmerLength) ; i++)
        {
            std::string kmer = contigs.substr (i,kmerLength);
            size_t countSame = BWTAlgorithms::countSequenceOccurrencesSingleStrand (kmer,pBWT);
            size_t countRevc = BWTAlgorithms::countSequenceOccurrencesSingleStrand (reverseComplement(kmer),pBWT);

            if ( (countSame >= threshold && countRevc >= erosion) ||  (countSame >= erosion && countRevc >= threshold) )
			// if ( (countSame >= threshold && countRevc >= threshold) )
            {
					start = i ;
					break;
            }
        }
    }

    if (pVertex->countEdges(ED_SENSE)== 0)
    {
        for (int i=(int)(seqLen-kmerLength); i >=0  ; i--)
        {
            std::string kmer = contigs.substr (i,kmerLength);
            size_t countSame = BWTAlgorithms::countSequenceOccurrencesSingleStrand (kmer,pBWT);
            size_t countRevc = BWTAlgorithms::countSequenceOccurrencesSingleStrand (reverseComplement(kmer),pBWT);
            if ( (countSame >= threshold && countRevc >= erosion) ||  (countSame >= erosion && countRevc >= threshold) )
			// if ( (countSame >= threshold) &&  (countRevc >= threshold) )
            {
                end = i+kmerLength ;
                break;
            }
        }
    }

    assert ( end >= start ) ;
    size_t length = end-start ;

    if ( length >= min_island && (pVertex->countEdges(ED_ANTISENSE)== 0 || pVertex->countEdges(ED_SENSE)== 0))
    {
			pVertex->setSeq(contigs.substr(start, length));

			//update edges coord
			EdgePtrVec edges = pVertex->getEdges(ED_SENSE) ;
			for(size_t i=0;i<edges.size();i++)
			{	
				edges[i]->updateSeqLen(length);
				edges[i]->offsetMatch(length-seqLen);//antisense may be trimmed
			}
			
			edges = pVertex->getEdges(ED_ANTISENSE) ;
			for(size_t i=0;i<edges.size();i++)
			{	//sense may be trimmed
				edges[i]->updateSeqLen(length);
			}

    }
	return false;
}


////////////////////////////////
// Remove illegal kmer edge  //
///////////////////////////////


void SGRemoveIllegalKmerEdgeVisitor::previsit(StringGraph* /*pGraph*/)
{
	//pGraph->sortVertexAdjListsByMatchLen();
	//pGraph->setColors(GC_WHITE);
	if(!m_filename.empty())
		std::ofstream m_fileHandle(m_filename.c_str());
}

bool SGRemoveIllegalKmerEdgeVisitor::visit(StringGraph* , Vertex* pVertex)
{
	bool changed=false;
	//get edges
	EdgePtrVec edges = pVertex->getEdges() ;
	std::string seq = pVertex->getStr();

	for (size_t i = 0 ; i < edges.size(); i++ )
	{
		Edge * pEdge = edges[i];
		size_t matchLen = pEdge->getMatchLength();

		//The cutoff tries to capture the lowest value of first overlap peak: kmerLength
		if( matchLen != m_kmerLength-1 ) continue;
		
		/*** Otherwise, this edge has short overlap. Proceed to check if the kmer frequency is low
		 If kmer frequency is high, this edge is due to repeat instead of high GC and should be removed.***/
		std::string Kmer;
		if(pEdge->getDir() == ED_SENSE)
			Kmer= seq.substr (seq.length()-matchLen-1, m_kmerLength);
		else
			Kmer= seq.substr (matchLen+1-m_kmerLength, m_kmerLength);

		bool IsKmerWeak =
			BWTAlgorithms::countSequenceOccurrencesSingleStrand(Kmer, pBWT) < m_threshold ||
			BWTAlgorithms::countSequenceOccurrencesSingleStrand(reverseComplement(Kmer), pBWT) < m_threshold ;

		//If this short edge was due to kmerization of low-kmer frequency, don't remove it 
		if (IsKmerWeak) continue;
		else
		{
			std::string anothorSeq = pEdge->getEnd()->getStr();
			std::string anotherKmer ;
			Edge * pTwin  = pEdge->getTwin();
			EdgeDir  anotherDir = pTwin->getDir();

			if (anotherDir == ED_SENSE )
				anotherKmer = anothorSeq.substr ( anothorSeq.length()-matchLen-1, m_kmerLength);
			else //anotherDir == ED_SENSE
				anotherKmer = anothorSeq.substr ( matchLen+1-m_kmerLength, m_kmerLength);

			bool IsAnotherKmerStrong  =
				BWTAlgorithms::countSequenceOccurrencesSingleStrand(anotherKmer, pBWT) >= m_threshold &&
				BWTAlgorithms::countSequenceOccurrencesSingleStrand(reverseComplement(anotherKmer), pBWT) >= m_threshold ;

				// //If both kmers flanking this edge is strong, this short edge is due to small overlap, remove it
			if (IsAnotherKmerStrong)
			{
				pEdge->setColor(GC_BLACK);
				pEdge->getTwin()->setColor(GC_BLACK);
				changed=true;
			}
		}
	}

	return changed;
}
void SGRemoveIllegalKmerEdgeVisitor::postvisit(StringGraph* pGraph)
{
	int num_remove = pGraph->sweepEdges(GC_BLACK);
	std::cout << "Remove " << num_remove/2 << " Edges by illegal kmer link" << std::endl;
}


//////////////////////////////
// Both Short Edges Remove  //
/////////////////////////////


void  SGBothShortEdgesRemoveVisitor::previsit(StringGraph* /*pGraph*/)
{
	std::cout << "Removing vertices with both short edges: ";
	std::cout << "VertexLength <=" << vertexLength  << " , " ;
	std::cout << "OverlapLength <=" << overlapLength ;
	if (pBWT!=NULL)
		std::cout << " ,  K: " << kmerLength  << " , T: " << threshold;
	std::cout << std::endl;

	//pGraph->setColors(GC_WHITE);
}

bool SGBothShortEdgesRemoveVisitor::visit(StringGraph* , Vertex* pVertex)
{

	if (pVertex->getSeqLen() > vertexLength || pVertex->getSeqLen() < kmerLength ||
		pVertex->countEdges(ED_ANTISENSE) ==0 ||  pVertex->countEdges(ED_SENSE) == 0 )
		return false;

	size_t maxOverlapLength [ED_COUNT] ;
	maxOverlapLength[ED_SENSE]  = pVertex->getLongestOverlapEdge(ED_SENSE)->getMatchLength();
	maxOverlapLength[ED_ANTISENSE] = pVertex->getLongestOverlapEdge(ED_ANTISENSE)->getMatchLength();

	bool changed = false ;

	if ( maxOverlapLength[ED_SENSE] <= overlapLength && maxOverlapLength[ED_ANTISENSE] <= overlapLength  )
	{
		float avgKmerFreqs = -1 ;

		if (pBWT!=NULL && kmerLength>0 && threshold>0)
		{

			std::string seq=pVertex->getStr();
			size_t numKmer = seq.length()-kmerLength+1 ;
			size_t sumKmerFreqs = 0 ;

			assert ( seq.length() >= kmerLength );

			for (size_t i = 0 ; i < numKmer  ; i++)
				sumKmerFreqs += BWTAlgorithms::countSequenceOccurrences(seq.substr (i,kmerLength), pBWT) ;
			avgKmerFreqs = (float) sumKmerFreqs/numKmer;
        }

        if(m_fileHandle!=NULL)
		m_fileHandle << ">" << pVertex->getID()  << " " << pVertex->getSeqLen()
					 << " " << maxOverlapLength[ED_SENSE]  << " " << maxOverlapLength[ED_ANTISENSE]  ;

		if ( avgKmerFreqs < 0 )
		{
			pVertex->setColor(GC_BLACK);
            if(m_fileHandle!=NULL) m_fileHandle << " BLACK";
			changed = true;
		}

		else if (avgKmerFreqs<=threshold)
		{
			pVertex->setColor(GC_BLACK);
            if(m_fileHandle!=NULL) m_fileHandle << " BLACK " << avgKmerFreqs ;
			changed = true;
        }

        if(m_fileHandle!=NULL)
            m_fileHandle << "\n" << pVertex->getSeq() << std::endl;
	}

	return changed;

}

void SGBothShortEdgesRemoveVisitor::postvisit(StringGraph* pGraph)
{
	//pGraph->sweepVertices(GC_BLACK);
	std::cout << "Remove " << pGraph->sweepVertices(GC_BLACK) << " chimera vertices"  << std::endl;
}

/////////////////////////////////////
//  Sweep Low Overlap Ratio Edge  //
////////////////////////////////////

void SGLowOverlapRatioEdgeSweepVisitor::previsit(StringGraph* pGraph)
{
	std::cout << "[ Low Overlap Ratio Edge Sweeper] :" ;
	std::cout << " min Overlap Ratio=" << m_overlapRatio;
	std::cout << ", Max Overlap Length=" << m_matchLength << std::endl;
	pGraph->setColors(GC_WHITE);
}

bool SGLowOverlapRatioEdgeSweepVisitor::visit(StringGraph* , Vertex* pVertex)
{
	bool changed = false;
	EdgePtrVec edges[ED_COUNT];
	edges[ED_SENSE] = pVertex->getEdges(ED_SENSE) ;
	edges[ED_ANTISENSE] = pVertex->getEdges(ED_ANTISENSE) ;

    if(pVertex->getSeqLen()>=m_min_vertex_size)
            return false;

	for(size_t idx = 0; idx < ED_COUNT; idx++)
	{
		EdgeDir dir = EDGE_DIRECTIONS[idx];
		size_t originLength = pVertex->getOriginLength(dir);

		for (size_t i=0;i<edges[dir].size();i++)
		{
			Edge * pEdge = edges[dir][i];
			size_t matchLen =pEdge->getMatchLength();
			
			//skip edges with match len > min matchLength
			if ( m_matchLength!=0 &&  matchLen > m_matchLength)	continue;

			Edge * pTwin = pEdge->getTwin();
			Vertex* pW = pEdge->getEnd();
			assert(pW->hasEdge(pTwin));
			size_t anotherOriginLength = pW->getOriginLength(pTwin->getDir());
			size_t minLength = originLength < anotherOriginLength ? originLength : anotherOriginLength;

			double ratio = (double) matchLen / minLength;

			if ( ratio < m_overlapRatio )
			{
				pEdge->setColor(GC_BLACK);
				pTwin->setColor(GC_BLACK);
				changed = true;

				if(m_fileHandle!=NULL)
                    m_fileHandle
                    << "Remove edge between " << pVertex->getID() <<"(" << originLength << ") and "
                    << pW->getID() <<"(" << anotherOriginLength << ") " << matchLen << "/" << minLength
                    << "(" << ratio << ")" << std::endl;
			}

		}//end of each edge

		//reset colors back to white if all edges are black
        bool isAllBlack=true;
        if(pVertex->getSeqLen()<m_min_vertex_size)
            isAllBlack=false;
        for (size_t i=0 ; i< edges[dir].size() ; i++)
        {
            if(edges[dir][i]->getColor()==GC_WHITE){
                isAllBlack=false;
                break;
            }
        }
        if(isAllBlack)
        {
            changed = false;
            //std::cout << edges[dir].size() <<"\n";
            for (size_t i=0 ; i< edges[dir].size() ; i++){
                edges[dir][i]->setColor(GC_WHITE);
                edges[dir][i]->getTwin()->setColor(GC_WHITE);
            }
        }

	}//end of dir

	return changed;
}
void SGLowOverlapRatioEdgeSweepVisitor::postvisit(StringGraph* pGraph)
{
	pGraph->sweepEdges(GC_BLACK);
	//std::cout << "Remove " << pGraph->sweepEdges(GC_BLACK)/2 << " small overlap ratio edges"  << std::endl;
}

////////////////////////////////////////////////////
// Remove Low Overlap Ratio Edge on Short Vertex  //
////////////////////////////////////////////////////

void SGLowOverlapRatioEdgeRemoveVisitor::previsit(StringGraph* pGraph)
{

	std::cout << "[ Remove Low Overlap Ratio Edge on Short Vertex ]" << std::endl;
	std::cout << "Max Vertex Length   :" << m_vetexLength << std::endl;
	std::cout << "min Ratio Difference:" << m_diffRatio << std::endl;
	pGraph->setColors(GC_WHITE);
}

//called by getDifferenceRatio
double getSingleOverlapRatio (Edge * pEdge)
{
	Vertex* pVertex = pEdge->getStart();
	assert(pVertex->hasEdge(pEdge));
	size_t matchLen = pEdge->getMatchLength();
	size_t originLen = pVertex->getOriginLength(pEdge->getDir());
	return (double) matchLen/originLen;
}

//return the overlap difference ratio between pV and pW
double getDifferenceRatio (Edge * pEdge)
{
	double r1 = getSingleOverlapRatio(pEdge);
	double r2 = getSingleOverlapRatio(pEdge->getTwin());
	double diff = r1>r2 ? r1-r2 : r2-r1;
	assert (diff >=0);
	return diff ;
}

bool SGLowOverlapRatioEdgeRemoveVisitor::visit(StringGraph* , Vertex* pVertex)
{
	//skip large or small vertices
    if(m_skipLargeVertex && pVertex->getSeqLen()>=m_vetexLength)
		return false;
	else if(!m_skipLargeVertex && pVertex->getSeqLen()<=m_vetexLength)
		return false;
		
	bool changed = false;
	EdgePtrVec edges[ED_COUNT];
	edges[ED_SENSE] = pVertex->getEdges(ED_SENSE) ;
	edges[ED_ANTISENSE] = pVertex->getEdges(ED_ANTISENSE) ;
	
	//check the overlap length diff for each edge
	for(size_t idx = 0; idx < ED_COUNT; idx++)
	{
		EdgeDir dir = EDGE_DIRECTIONS[idx];

		for (size_t i=0;i<edges[dir].size();i++)
		{
			Edge * pEdge = edges[dir][i];

			double vwDiff = getDifferenceRatio(pEdge);
			if ( vwDiff <=m_diffRatio ) continue;

			//Otherwise, vwDiff is large, continue to check any better reliable edge exist in w_edges from pW
			Vertex* pW = pEdge->getEnd();
			Edge * pTwin = pEdge->getTwin();
			assert(pW->hasEdge(pTwin));

			EdgePtrVec w_edges= pW->getEdges(pTwin->getDir()) ;

			bool existReliable =false ;
			double reliableDiff = 1 ;
			//find a better edge from pW if possible
			for (size_t i=0;i<w_edges.size();i++)
			{
				if (w_edges[i] == pTwin) continue;
				double  rDiff = getDifferenceRatio(w_edges[i]);
				if ( rDiff <= m_diffRatio )
				{
					//a better pW edge is is found, it is safe to remove this edge
					existReliable = true;
					if (rDiff<reliableDiff) reliableDiff = rDiff;
				}
			}
			if (existReliable) //reliable and no PE support, test by YT
			{
				assert (reliableDiff < 1) ;
				pEdge->setColor(GC_BLACK);
				pTwin->setColor(GC_BLACK);
				changed = true;

				if(m_fileHandle!=NULL)
                    m_fileHandle
                    << "Remove edge between " << pVertex->getID() <<" and " <<  pW->getID() << "with ratio difference: " << vwDiff  << ".\n"
                    << pW->getID() << " exist another low overlap ratio difference (" <<  reliableDiff << ")." << std::endl;
			}

		}

	}
	return changed;
}

void SGLowOverlapRatioEdgeRemoveVisitor::postvisit(StringGraph* pGraph)
{
	std::cout << "Remove " << pGraph->sweepEdges(GC_BLACK)/2 << " small overlapping "  << std::endl;
}

void SGSubGraphVisitor::previsit(StringGraph* /*pGraph*/)
{
	num_straight = 0 ;
	num_terminal = 0;
	num_island = 0;
	num_monobranch = 0;
	num_dibranch = 0;
	num_simple = 0;
	num_edges = 0;
	num_vertex = 0;
	sum_edgeLen = 0;
}

bool SGSubGraphVisitor::visit(StringGraph* /*pGraph*/, Vertex* pVertex)
{
	int s_count = pVertex->countEdges(ED_SENSE);
	int as_count = pVertex->countEdges(ED_ANTISENSE);

	if(s_count >= 1 && as_count == 0 && pVertex->getSeqLen() > 10000)
    {
        std::cout << pVertex->getID() << "\n" << s_count << ":" << as_count << "\n";
    }

	if(s_count == 0 && as_count >= 1 && pVertex->getSeqLen() > 10000)
    {
        std::cout << pVertex->getID() << "\n" << s_count << ":" << as_count << "\n";;
    }

/*
    if(s_count >1 && as_count >1)
	{
        std::cout<< ">" << pVertex->getSeqLen() << "\n" << pVertex->getStr() << "\n";
		if(s_count>1)
        {
            EdgePtrVec edges=pVertex->getEdges(ED_SENSE);
            for(size_t i=0; i<edges.size();i++)
            {
                std::cout << ">1 " << edges[i]->getMatchLength() <<"\n" << edges[i]->getEnd()->getStr() << "\n";
            }
        }
        if(as_count>1)
        {
            EdgePtrVec edges=pVertex->getEdges(ED_ANTISENSE);
            for(size_t i=0; i<edges.size();i++)
            {
                std::cout << ">2 " << edges[i]->getMatchLength() <<"\n" << edges[i]->getEnd()->getStr() << "\n";
            }
        }

		getchar();
	}
*/
/*
	if(s_count > 1 && as_count > 1)
	++num_dibranch;
	else if(s_count > 1 || as_count > 1)
	++num_monobranch;

	if(s_count == 1 || as_count == 1)
	++num_simple;

	if(s_count == 1 && as_count == 1)
    {
        std::cout<< ">1 " << pVertex->getID() << "\n" << pVertex->getStr() << "\n";
        ++num_straight;
        EdgePtrVec edges = pVertex->getEdges(ED_SENSE);
        Vertex* pWVert = edges[0]->getEnd();
        std::cout<< ">2 " << pWVert->getID() << "\n" << pWVert->getStr() << pWVert->countEdges(ED_ANTISENSE) << "\n";
        edges = pVertex->getEdges(ED_ANTISENSE);
        pWVert = edges[0]->getEnd();
        std::cout<< ">3 " << pWVert->getID() << "\n" << pWVert->getStr() << pWVert->countEdges(ED_SENSE) << "\n";
    }


	num_edges += (s_count + as_count);
	++num_vertex;

	EdgePtrVec edges = pVertex->getEdges();
	for(size_t i = 0; i < edges.size(); ++i)
	sum_edgeLen += edges[i]->getSeqLen();
*/
	return false;
}

//
void SGSubGraphVisitor::postvisit(StringGraph* pGraph)
{
	std::cout << "Remove " << pGraph->sweepVertices(GC_BLACK) << " BLACK vertices"  << std::endl;
	//printf("Vertices: %d Edges: %d Islands: %d Tips: %d Monobranch: %d Dibranch: %d Simple: %d Straight: %d\n",
//	num_vertex, num_edges,num_island, num_terminal,num_monobranch, num_dibranch, num_simple,num_straight);
}


void SGRemoveEdgeByPEVisitor::previsit(StringGraph* pGraph)
{
	std::cout << "[ SGRemoveEdgeByPEVisitor ]\t Kmer: " << m_kmerSize 
					<< "\t Insert Size: "  << m_insertSize << "\t Min PE count: " << m_minPEcount << std::endl;
    pGraph->setColors(GC_WHITE);
	pGraph->sortVertexAdjListsByMatchLen();// sort by match len in ascending order
	
    m_edgecount=0;
	
	// KmerDistribution m_kd = BWTAlgorithms::sampleKmerCounts(m_kmerSize, 100000, m_indices.pRBWT);
	// m_repeatKmerCutoff = m_kd.getMedian()*2; 

}

bool SGRemoveEdgeByPEVisitor::addReadIDsAtPos(NameSet& pVReadID, std::string& pVertexSeq, int overlapBoundaryPos)
{
	if(overlapBoundaryPos<0) overlapBoundaryPos=0;
	//get the kmer seed
	std::string seed=pVertexSeq.substr(overlapBoundaryPos, m_kmerSize);

	//add read IDs of kmer seed into hashmap
	pVReadID.addFirstReadIDs(seed);
	pVReadID.addFirstReadIDs(reverseComplement(seed));
	return true;

}

bool SGRemoveEdgeByPEVisitor::visit(StringGraph* /*pGraph*/, Vertex* pVertex)
{
    bool changed=false;

	// if(pVertex->getSeqLen() < m_insertSize) return false;
	
	for(size_t idx = 0; idx < ED_COUNT; idx++)
	{
		EdgeDir dir = EDGE_DIRECTIONS[idx];
        EdgePtrVec edges=pVertex->getEdges(dir);
		if(edges.size() < 1) continue;
		
		std::string pVertexSeq = (dir==ED_SENSE) ? pVertex->getStr() : reverseComplement(pVertex->getStr());
		NameSet pVReadID(m_indices.pBWT, m_indices.pSSA, 600);    //hashmap rapper for storing ReadIDs

		//compute the Pos left next to the matching boundary
		// int overlapBoundaryPos=pVertex->getSeqLen() - edges.back()->getMatchLength()-1;
		// addReadIDsAtPos(pVReadID, pVertexSeq, overlapBoundaryPos);

		// // move Pos left m_kmer
		// overlapBoundaryPos -= m_kmerSize/2;
		// addReadIDsAtPos(pVReadID, pVertexSeq, overlapBoundaryPos);

		// //move Pos left m_kmer*2
		// overlapBoundaryPos -= m_kmerSize/2;
		// addReadIDsAtPos(pVReadID, pVertexSeq, overlapBoundaryPos);

		// //move Pos back to middle kmer 
		// overlapBoundaryPos += m_kmerSize/2;

		// if(overlapBoundaryPos<0) overlapBoundaryPos=0;
		
		// //pVAnotherID stores the read IDs of the other end
		// std::vector<int64_t> pVAnotherID;
		// pVReadID.getAnotherReadIDs( pVAnotherID);		

		// //BFS search to insert size+variance with 256 leaves
        SGWalkVector walkVector;
        SGSearch::getTree(pVertex, dir, (size_t) m_insertSize*1.5, 128, walkVector);
		assert(walkVector.size()>=1);

		//now check each walk for existence of paired read IDs
        std::vector<NameSet> goals(walkVector.size(), NameSet(m_indices.pBWT, m_indices.pSSA, 600));
		const size_t InsertVariance=m_kmerSize/2+1;

		// scan for PE support of each edge.
        for(size_t k=0; k<edges.size(); k++)
        {
		    Edge *pEdge=edges[k];
			
			// Don't remove confident edge with sufficient overlap
		    if (pEdge->getMatchLength() >= m_insertSize*0.8)
				continue;

			std::string pVertexSeq = (dir==ED_SENSE) ? pVertex->getStr() : reverseComplement(pVertex->getStr());
			NameSet pVReadID(m_indices.pBWT, m_indices.pSSA);    //hashmap rapper for storing ReadIDs

			//compute the Pos left next to the matching boundary
			int overlapBoundaryPos=pVertex->getSeqLen() -pEdge->getMatchLength()-1;
			addReadIDsAtPos(pVReadID, pVertexSeq, overlapBoundaryPos);

			// move Pos left m_kmer
			overlapBoundaryPos -= m_kmerSize/2;
			addReadIDsAtPos(pVReadID, pVertexSeq, overlapBoundaryPos);

			//move Pos left m_kmer*2
			overlapBoundaryPos -= m_kmerSize/2;
			addReadIDsAtPos(pVReadID, pVertexSeq, overlapBoundaryPos);

			//move Pos back to middle kmer 
			overlapBoundaryPos += m_kmerSize/2;

			if(overlapBoundaryPos<0) overlapBoundaryPos=0;
			
			//pVAnotherID stores the read IDs of the other end
			std::vector<int64_t> pVAnotherID;
			pVReadID.getAnotherReadIDs( pVAnotherID);		

            //collect read IDs at possible PE ending pos
			size_t PEcount=0;
            for (size_t i=0;i<walkVector.size();i++)
            {
                if(walkVector[i].getFirstEdge()!=pEdge)
                    continue;

                //the antisense seq should be reverse complement, see the getString in SGWalk.cpp
                std::string WalkSeq = (dir==ED_SENSE) ? walkVector[i].getString(SGWT_START_TO_END, NULL)
                                                      : reverseComplement(walkVector[i].getString(SGWT_START_TO_END, NULL));

                //Find the read IDs of the other end at insert size 
				for(int targetOffset=-InsertVariance; targetOffset<=(int)InsertVariance; targetOffset+=InsertVariance)
				{
					size_t targetPos = overlapBoundaryPos+m_insertSize+targetOffset;
					if(WalkSeq.length() >= targetPos )
					{
						std::string endingkmer=WalkSeq.substr(targetPos-m_kmerSize, m_kmerSize);
						goals[i].addSecondReadIDs(endingkmer);
						goals[i].addSecondReadIDs(reverseComplement(endingkmer));
					}
					// else	// most are small repeats leading to numerous walks
						// std::cout << WalkSeq.length() << "\t" << pVertex->getSeqLen() << "\t" << edges.size() << "\t" 
								// << targetPos << "\t" << walkVector.size() << "\n";
				}
				
				//find intersection of read IDs in pVAnotherID and goals[i]
                for(size_t j=0; j<pVAnotherID.size(); j++)
                {
                   if(goals[i].exist(pVAnotherID[j]))
                        PEcount++;
					if(PEcount>=m_minPEcount) break;
                }
				
				if(PEcount>=m_minPEcount) break;
            }//end of each goal[i]

            if(PEcount < m_minPEcount)
            {
				//dangerous during multithread
                if(pEdge->getColor()==GC_WHITE)
                {
                    pEdge->setColor(GC_BLACK);
                    pEdge->getTwin()->setColor(GC_BLACK);
                    changed=true;
                    m_edgecount++;
					// std::cout << pVReadID.sizeOfFirstReadIDs() << "\t" << goals[k].sizeOfSecondReadIDs() << "\t" << pVertex->getSeqLen() 
								// << "\t" << pEdge->getMatchLength() << "\n";
                }
				// else if(pEdge->getColor()==GC_WHITE)
				// {
					// pEdge->setColor(GC_RED);
                    // pEdge->getTwin()->setColor(GC_RED);
				// }
            }

        }//end of for each edge
	}//end of for each dir

	return changed;
}
void SGRemoveEdgeByPEVisitor::postvisit(StringGraph* pGraph)
{
	std::cout << "RemoveEdgeByPE: Remove " << pGraph->sweepEdges(GC_BLACK)/2 << " edges without PE by insert size " << m_insertSize << std::endl;
	//std::cout << "MarkEdgeByPE: Mark " << m_edgecount << " PE good edges by insert size " << m_insertSize << std::endl;

}


//simply remove edges with overlap < min_overlap
void SGRemoveByOverlapLenDiffVisitor::previsit(StringGraph* pGraph)
{
	pGraph->setColors(GC_WHITE);
	pGraph->sortVertexAdjListsByMatchLen();// sort by match len in ascending order
}

bool SGRemoveByOverlapLenDiffVisitor::visit(StringGraph* , Vertex* pVertex)
{
    bool changed=false;

    if(pVertex->getSeqLen()<m_min_vertex_size)
        return false;

	EdgePtrVec edges[ED_COUNT] ;
	edges[ED_SENSE] = pVertex->getEdges(ED_SENSE) ;
	edges[ED_ANTISENSE] = pVertex->getEdges(ED_ANTISENSE) ;

	for (size_t dir=0 ; dir< ED_COUNT ; dir ++)
	{
	    if(edges[dir].size()<=1) continue;
		size_t maxlen=edges[dir].back()->getMatchLength();

		// if(m_min_overlap>0)
        if(m_min_overlap>0 && maxlen > m_min_overlap)
        {
            for (size_t i=0 ; i< edges[dir].size() ; i++)
            {
                if (edges[dir][i]->getMatchLength() < m_min_overlap)
                {
                    changed=true;
                    edges[dir][i]->setColor(GC_BLACK);
                    edges[dir][i]->getTwin()->setColor(GC_BLACK);
                }
            }
        }

		// if(m_max_overlapdiff>0)
        if(m_max_overlapdiff>0 && maxlen - edges[dir].front()->getMatchLength() >= m_max_overlapdiff)
        {
            for (size_t i=0 ; i< edges[dir].size()-1 ; i++)
            {
                size_t diff = maxlen-edges[dir][i]->getMatchLength();
                if(diff >= m_max_overlapdiff)
                {
                    changed=true;
                    edges[dir][i]->setColor(GC_BLACK);
                    edges[dir][i]->getTwin()->setColor(GC_BLACK);
                }
            }
        }

        bool isAllEdgeBlack=m_islandProtect;
        for (size_t i=0 ; i< edges[dir].size() ; i++)
        {
            if(edges[dir][i]->getColor()==GC_WHITE)
                isAllEdgeBlack=false;
        }
        if(isAllEdgeBlack)  //all edges are colored black by pVertex
        {
            for (size_t i=0 ; i< edges[dir].size() ; i++){
                edges[dir][i]->setColor(GC_WHITE);
                edges[dir][i]->getTwin()->setColor(GC_WHITE);
            }
            changed=false;
        }
	}

	return changed;
}
void SGRemoveByOverlapLenDiffVisitor::postvisit(StringGraph* pGraph)
{
	int num_remove = pGraph->sweepEdges(GC_BLACK);
	std::cout << "SGRemoveByOverlapLenDiffVisitor: Remove " << num_remove/2
	<< " Edges with min_vertex_size:min_overlap:max_diff "
	<< m_min_vertex_size << ":" << m_min_overlap <<":" << m_max_overlapdiff<<std::endl;
}



//Collect PE reads at ends of each island/tip
//Store PE read IDs into NameSet hashtable
void SGIslandCollectVisitor::previsit(StringGraph* /*pGraph*/)
{
    m_islandcount=0;
	m_kd = BWTAlgorithms::sampleKmerCounts(m_kmerSize, 100000, m_indices.pRBWT);
	m_repeatKmerCutoff = m_kd.getCutoffForProportion(0.75); 
	m_kd.computeKDAttributes();
	// m_repeatKmerCutoff =  m_kd.getMedian()*1.3;

	std::cout << "\n[ Collect paired-end reads mapped onto islands/tips ]" <<std::endl;
	std::cout << "Median kmer freq: " << m_kd.getMedian() <<  "\t Repeat kmer cutoff: " << m_repeatKmerCutoff << "\t minimum island/tip size: " << m_minIslandSize 
				  << "\t kmer size: "<< m_kmerSize << "\t insert size: "<< m_insertSize <<  "\n";
}

bool SGIslandCollectVisitor::visit(StringGraph* /*pGraph*/, Vertex* pVertex)
{
    bool changed=false;

    if( (pVertex->countEdges(ED_SENSE)==0 || pVertex->countEdges(ED_ANTISENSE)==0)  && pVertex->getSeqLen()>=m_minIslandSize)
    {
        m_islandcount++;
        std::string pVstr=pVertex->getStr();
        NameSet pVPrefixFwdID(m_indices.pBWT, m_indices.pSSA);    //hashmap rapper for storing ReadIDs
        NameSet pVPrefixRvcID(m_indices.pBWT, m_indices.pSSA);    //hashmap rapper for storing ReadIDs
        NameSet pVSuffixFwdID(m_indices.pBWT, m_indices.pSSA);    //hashmap rapper for storing ReadIDs
        NameSet pVSuffixRvcID(m_indices.pBWT, m_indices.pSSA);    //hashmap rapper for storing ReadIDs

		//collect read IDs mapped on this vertex upto insert size
        for(size_t i=0; i<m_insertSize; i+=20)
        {
			//Collect SENSE read IDs
            if(pVertex->countEdges(ED_SENSE)==0)
            {
                std::string seed = pVstr.substr(pVstr.length()-i-m_kmerSize, m_kmerSize);
				
				//skip repeat seeds
				size_t KmerFreq = BWTAlgorithms::countSequenceOccurrences( seed, m_indices.pBWT );
				if( KmerFreq < m_repeatKmerCutoff )
				{
					pVSuffixFwdID.addReadIDAndContigID(seed, m_tslv, pVertex, SenseFwd);
					pVSuffixRvcID.addReadIDAndContigID(reverseComplement(seed), m_tslv, pVertex, SenseRvc);
				}
            }

			//Collect ANTISENSE read IDs
            if(pVertex->countEdges(ED_ANTISENSE)==0)
            {
                std::string seed = pVstr.substr(i,m_kmerSize);
				
				//skip repeat seeds
				size_t KmerFreq = BWTAlgorithms::countSequenceOccurrences( seed, m_indices.pBWT );
				if( KmerFreq < m_repeatKmerCutoff )
				{
					pVPrefixFwdID.addReadIDAndContigID(seed, m_tslv, pVertex, AntisenseFwd);
					pVPrefixRvcID.addReadIDAndContigID(reverseComplement(seed), m_tslv, pVertex, AntisenseRvc);
				}
            }
        }

		(pVertex->pVReadIDs).push_back(pVPrefixFwdID.getReadIDs());
		(pVertex->pVReadIDs).push_back(pVPrefixRvcID.getReadIDs());
		(pVertex->pVReadIDs).push_back(pVSuffixFwdID.getReadIDs());
		(pVertex->pVReadIDs).push_back(pVSuffixRvcID.getReadIDs());
        changed=true;
    }

	return changed;
}
void SGIslandCollectVisitor::postvisit(StringGraph* /*pGraph*/)
{
	std::cout << "IslandCollect: Collect " << m_islandcount << " islands/tips for FM-index walk\n\n ";
}

void SGJoinIslandVisitor::previsit(StringGraph* /*pGraph*/)
{
    m_islandcount=0;	
	std::cout << "[ Perform FM-index walk between islands/tips ]" <<std::endl;
	std::cout << "Minimum PE support: " << m_minPEcount <<"\t Kmer Size: " << m_kmer<< std::endl;
}


//return SA index of the other end of a PE read
//assume PE reads are 0: R1/1, 
//								 1: R1/2, 
//								 2: R2/1, 
int64_t SGJoinIslandVisitor::getAnotherID(int64_t idx)
{
	if(idx%2==0)	//first end
		return idx+1;
	else					//second end
		return idx-1;
}

// Search for neighbor vertex pW having PE support wrt pV at islandDir
// iterate through each read ID stored in pVReadIDs[islandDir]
// return a hashmap pWIDs storing pWs with PE support
void SGJoinIslandVisitor::findNeighborWithPESupport(Vertex* pV, size_t islandDir, SparseHashMap<VertexID, size_t*, StringHasher>& pWIDs)
{
	for(size_t i=0; i<pV->pVReadIDs[islandDir].size(); i++)
	{
		//convert the read ID mapped on pV into PEID of the other end
		int64_t PEID=getAnotherID(pV->pVReadIDs[islandDir][i]);
		//retrieve the vertices pWs containing PEID into currentList 
		ThreadSafeList<std::pair<Vertex*, ReadOnContig> >  currentList= m_tslv->at(PEID);
		// std::cout << currentList.size() << "\n";
		
		//compute the PE mapping frequency for each pW in currentList
		for(std::list<std::pair<Vertex*, ReadOnContig> >::iterator slit=currentList.begin(); slit!=currentList.end(); slit++)
		{
			Vertex* pW=slit->first;
			ReadOnContig roc=slit->second;

			SparseHashMap<VertexID, size_t*, StringHasher>::iterator islandIter;
			islandIter=pWIDs.find(pW->getID());
			if(islandIter!=pWIDs.end())	//pWID already has this vertex mapped, increment the mapping frequency
			{
				islandIter->second[roc]++;
			}
			else	//first time, pWID does not have this vertex yet
			{
				size_t* readcount = new size_t[4];	//memory leak
				readcount[0]=readcount[1]=readcount[2]=readcount[3]=0;
				readcount[roc]++;
				pWIDs.insert(std::make_pair(pW->getID(), readcount));
			}
		}
	}
	
	// std::cout << pWIDs.size() << "\n";
	// if(pWIDs.size()>=3) continue;	//pV at islandDir is probably repeat and skip 
}
void SGJoinIslandVisitor::updateExtendedVertex(Vertex* pVertex, std::string& newStr, EdgeDir dir)
{
	pVertex->setSeq(newStr);
	
	// if pVertex is a tip with edges at dir, those edge lengths need to be updated for correct simplify()
	// see simplify() for detailed updating reasons
	EdgePtrVec  edges=pVertex->getEdges(dir);
	for(size_t i=0; i<edges.size(); i++)
		edges[i]->updateSeqLen(newStr.length());

}

bool SGJoinIslandVisitor::visit(StringGraph* pGraph, Vertex* pVertex)
{
    bool changed=false;

	//this visitor only processes islands or tips with size >=m_minIslandSize
    if( (pVertex->countEdges(ED_SENSE)==0 || pVertex->countEdges(ED_ANTISENSE)==0)  && pVertex->getSeqLen()>=m_minIslandSize)
    {
        Vertex *pV=pVertex;	

        //For each pV, look for PE support at islandDir
		//islandDir= 0,1,2,3=AntisenseFwd, AntisenseRvc, SenseFwd, SenseRvc of pV
        for(size_t islandDir=0; islandDir<4; islandDir++)
        {
			//Skip tips with edges already at islandDir
			if( pV->countEdges(ED_ANTISENSE)>0 && (islandDir==AntisenseFwd || islandDir==AntisenseRvc)) continue;
			if( pV->countEdges(ED_SENSE)>0 && (islandDir==SenseFwd || islandDir==SenseRvc)) continue;

			//pWIDs stores IDs of all islands/tips having PE mapped onto pV
			SparseHashMap<VertexID, size_t*, StringHasher> pWIDs;
			//Search for vertex pW having PE support wrt pV at islandDir
			findNeighborWithPESupport(pV, islandDir, pWIDs);
			
			//for each pW, perform FM-index walk between pW and pV
			SparseHashMap<VertexID, size_t*, StringHasher>::iterator islandIter;
            for(islandIter=pWIDs.begin(); islandIter!=pWIDs.end() ; islandIter++)
            {
                Vertex *pW=pGraph->getVertex(islandIter->first);
				if(pV==pW) continue;
				
				//retrieve the PE mapping frequency of pW
				size_t PreFwdCount=islandIter->second[0], PreRvcCount=islandIter->second[1], SufFwdCount=islandIter->second[2], SufRvcCount=islandIter->second[3];
				
				//Skip impossible <=minPEcount cases for speedup 
				if( islandDir==AntisenseFwd && PreFwdCount<=m_minPEcount && SufRvcCount<=m_minPEcount) continue;
				else if( islandDir==AntisenseRvc && PreRvcCount<=m_minPEcount && SufFwdCount<=m_minPEcount) continue;
				else if( islandDir==SenseFwd && SufFwdCount<=m_minPEcount && PreRvcCount<=m_minPEcount) continue;
				else if( islandDir==SenseRvc && SufRvcCount<=m_minPEcount && PreFwdCount>m_minPEcount) continue;

				// std::cout << pV->getID() << "<->" << pW->getID() <<"\t"<< PreFwdCount <<":"<< PreRvcCount <<":"<< SufFwdCount <<":"<< SufRvcCount <<"\t"<< islandDir <<"\n";
					
               	//Deadlock prevention: the vertex with smaller memory address is always locked first, which prevent mutual locking pV<->pW at two threads
				pGraph->lockTwoVertices(pV,pW);
				
				std::string pVstr=pV->getStr();
				std::string pWstr=pW->getStr();

				//Total 4*4=16 cases, 8 cases use the same operations
                //Among them, only 4 cases are feasible link between two contigs
                
				//case 1: PreFwd:PreFwd or PreRvc:PreRvc, needs reverse complement of one				
				if( ((islandDir==0 && PreFwdCount>m_minPEcount) || (islandDir==1 && PreRvcCount>m_minPEcount))
					&& pV->countEdges(ED_ANTISENSE)==0 && pW->countEdges(ED_ANTISENSE)==0 )
				{
					//Step 1: found PE links between pV and pW
					//std::cout << "case1: " << pV->getID() << ":" << pW->getID() << "\n";
					std::string pWStrRvc=reverseComplement(pWstr);
					
					//Step 2: perform FM index walk between pV and pW
					for(size_t i=0; i<m_numOfIterations; i++)	//perform 2 runs of FM index walk
					{
						std::string StartStr=pWStrRvc.substr(0,pWStrRvc.length()-i*m_kmer);
						SAIntervalTree SAITree(&StartStr, m_kmer, 100, StartStr.length()+m_SAISearchDepth, m_SAISearchLeaves, m_indices, pVstr,1,true);
						std::string mergedseq;
						SAITree.mergeTwoReads(mergedseq);
						
						if( !mergedseq.empty() )
						{
							//pV and pW are connected
							std::string pWStrNew=mergedseq.substr(0,mergedseq.length()-pVstr.length()+m_kmer);
							pW->setSeq(reverseComplement(pWStrNew));
							
							// if pW is a tip with sense edges, edge lengths need to be updated for correct simplify()
							EdgePtrVec  edges=pW->getEdges(ED_SENSE);
							for(size_t i=0; i<edges.size(); i++){
								edges[i]->updateSeqLen(pWStrNew.length());
								edges[i]->offsetMatch(pWStrNew.length()-pWstr.length());	
							}
							
							//Step 3: Create new edges between pV and pW
							SeqCoord coordV(0,m_kmer-1,pVstr.length());
							SeqCoord coordW(0,m_kmer-1,pWStrNew.length());
							Edge* edgeVW=new Edge(pW,ED_ANTISENSE,EC_REVERSE,coordV);
							Edge* edgeWV=new Edge(pV,ED_ANTISENSE,EC_REVERSE,coordW);
							edgeVW->setTwin(edgeWV);
							edgeWV->setTwin(edgeVW);
							pGraph->addEdge(pV,edgeVW);
							pGraph->addEdge(pW,edgeWV);
							
							m_islandcount++;
							break;
						} 
					}// end of multi rounds
				}//end of case 1
				
				//case4: PreFwd:SufRvc or PreRvc:SufFwd
				else if( ((islandDir==0 && SufRvcCount>m_minPEcount) || (islandDir==1 && SufFwdCount>m_minPEcount))
						&& pV->countEdges(ED_ANTISENSE)==0 && pW->countEdges(ED_SENSE)==0 )
				{
					//std::cout << "case4: " << pV->getID() << ":" << pW->getID()<< "\n";
					//Step 2: perform FM index walk between pV and pW
					for(size_t i=0; i<m_numOfIterations; i++)	//perform 2 runs of FM index walk
					{
						std::string StartStr=pWstr.substr(0,pWstr.length()-i*m_kmer);
						SAIntervalTree SAITree(&StartStr, m_kmer, 100, StartStr.length()+m_SAISearchDepth, m_SAISearchLeaves,m_indices, pVstr,1,true);
						std::string mergedseq;
						SAITree.mergeTwoReads(mergedseq);

						if( !mergedseq.empty() ){
							//extend and update pW and its edges
							std::string pWStrNew=mergedseq.substr(0,mergedseq.length()-pVstr.length()+m_kmer);
							updateExtendedVertex(pW,  pWStrNew, ED_ANTISENSE);
							
							//Step 3: Create new edges between pV and pW
							SeqCoord coordV(0,m_kmer-1,pVstr.length());
							SeqCoord coordW(pWStrNew.length()-m_kmer,pWStrNew.length()-1,pWStrNew.length());
							Edge* edgeVW=new Edge(pW,ED_ANTISENSE,EC_SAME,coordV);
							Edge* edgeWV=new Edge(pV,ED_SENSE,EC_SAME,coordW);
							edgeVW->setTwin(edgeWV);
							edgeWV->setTwin(edgeVW);
							pGraph->addEdge(pV,edgeVW);
							pGraph->addEdge(pW,edgeWV);
							
							m_islandcount++;
							break;
						}
					}
				}

				//case5: SufFwd:SufFwd or SufRvc:SufRvc
				else if( ((islandDir==2 && SufFwdCount>m_minPEcount) || (islandDir==3 && SufRvcCount>m_minPEcount))
						&& pV->countEdges(ED_SENSE)==0 && pW->countEdges(ED_SENSE)==0 )
				{
					//std::cout << "case5: " << pV->getID() << ":" << pW->getID()<< ":"<< islandDir << ":" << SufRvcCount << "\n";
					std::string pWStrRvc=reverseComplement(pWstr);
					
					//Step 2: perform FM index walk between pV and pW
					for(size_t i=0; i<m_numOfIterations; i++)	//perform 2 runs of FM index walk
					{
						std::string StartStr=pVstr.substr(0,pVstr.length()-i*m_kmer);
						SAIntervalTree SAITree(&StartStr, m_kmer, 100, StartStr.length()+m_SAISearchDepth, m_SAISearchLeaves, m_indices, pWStrRvc,1,true);
						std::string mergedseq; 
						SAITree.mergeTwoReads(mergedseq);

						if( !mergedseq.empty() )
						{
							//extend and update pV and its edges
							std::string pVStrNew=mergedseq.substr(0, mergedseq.length()-pWstr.length()+m_kmer);
							updateExtendedVertex(pV,  pVStrNew, ED_ANTISENSE);

							//Step 3: create new edges for these two vertices
							SeqCoord coordV(pVStrNew.length()-m_kmer,pVStrNew.length()-1,pVStrNew.length());
							SeqCoord coordW(pWstr.length()-m_kmer,pWstr.length()-1,pWstr.length());
							Edge* edgeVW=new Edge(pW,ED_SENSE,EC_REVERSE,coordV);
							Edge* edgeWV=new Edge(pV,ED_SENSE,EC_REVERSE,coordW);
							edgeVW->setTwin(edgeWV);
							edgeWV->setTwin(edgeVW);
							pGraph->addEdge(pV,edgeVW);
							pGraph->addEdge(pW,edgeWV);
							
							m_islandcount++;
							break;
						}
					}
				} 
				
				//case8: SufFwd:PreRvc or SufRvc:PreFwd
				else if( ((islandDir==2 && PreRvcCount>m_minPEcount) || (islandDir==3 && PreFwdCount>m_minPEcount))
						&& pV->countEdges(ED_SENSE)==0 && pW->countEdges(ED_ANTISENSE)==0 )
				{
					//std::cout << "case8: " << pV->getID() << ":" << pW->getID()<< ":"<< islandDir << ":" << PreFwdCount <<"\n";
					//Step 2: perform FM index walk between pV and pW
					for(size_t i=0; i<m_numOfIterations; i++)	//perform 2 runs of FM index walk
					{
						std::string StartStr=pVstr.substr(0,pVstr.length()-i*m_kmer);
						SAIntervalTree SAITree(&StartStr, m_kmer, 100, StartStr.length()+m_SAISearchDepth, m_SAISearchLeaves,m_indices, pWstr,1,true);
						std::string mergedseq; 
						SAITree.mergeTwoReads(mergedseq);
						
						if( !mergedseq.empty() ){
							//extend and update pV and its edges
							std::string pVStrNew=mergedseq.substr(0, mergedseq.length()-pWstr.length()+m_kmer);
							updateExtendedVertex(pV,  pVStrNew, ED_ANTISENSE);
							
							//Step 3: create new edges for these two vertices
							SeqCoord coordV(pVStrNew.length()-m_kmer,pVStrNew.length()-1,pVStrNew.length()); //dummy coord
							SeqCoord coordW(0,m_kmer-1,pWstr.length());
							Edge* edgeVW=new Edge(pW,ED_SENSE,EC_SAME,coordV);
							Edge* edgeWV=new Edge(pV,ED_ANTISENSE,EC_SAME,coordW);
							edgeVW->setTwin(edgeWV);
							edgeWV->setTwin(edgeVW);
							pGraph->addEdge(pV,edgeVW);
							pGraph->addEdge(pW,edgeWV);
							
							m_islandcount++;
							break;
						}
					}
				}
				//remember to unlock vertices flanking this edge
				pGraph->unlockTwoVertices(pV,pW);
				
				/*************** remaining cases are infeasible
                else if((islandDir==0 && PreRvcCount>0) || (islandDir==1 && PreFwdCount>0))
                    std::cout << "case2: " << pV->getID() << ":" << pW->getID() << "\n";
                else if((islandDir==0 && SufFwdCount>0) || (islandDir==1 && SufRvcCount>0))
                    std::cout << "case3: " << pV->getID() << ":" << pW->getID()<< "\n";
                //SufFwd:SufRvc or SufRvc:SufFwd
                else if((islandDir==2 && SufRvcCount>0) || (islandDir==3 && SufFwdCount>0))
                    std::cout << "case6: " << pV->getID() << ":" << pW->getID()<< "\n";
                //SufFwd:PreFwd or SufRvc:PreRvc
                else if((islandDir==2 && PreFwdCount>0) || (islandDir==3 && PreRvcCount>0))
                    std::cout << "case7: " << pV->getID() << ":" << pW->getID()<< "\n";
				****************/
            }//end of for each neighbor pW
        }

        changed=true;
    }

	return changed;
}
void SGJoinIslandVisitor::postvisit(StringGraph* pGraph)
{
	std::cout << "SGJoinIslandVisitor: Join " << m_islandcount << " islands/tips using FM-index walk " << std::endl;
	pGraph->simplify();

}



/********** NameSet Class Body****************/

void NameSet::addFirstReadIDs(std::string seed)
{
	BWTInterval interval = BWTAlgorithms::findInterval(pBWT, seed);
	// if(interval.isValid())
	if(interval.isValid() /*&& interval.size()<200*/)
	{
		for(int64_t j = interval.lower; j <= interval.upper && (j-interval.lower<m_maxIDs); ++j)
		{
		    int64_t SAindex = pSSA->calcSA(j,pBWT).getID();
			m_SAindicesSet1.insert(SAindex);
		}
	}
}

void NameSet::addSecondReadIDs(std::string seed)
{
	BWTInterval interval = BWTAlgorithms::findInterval(pBWT, seed);
	// if(interval.isValid())
	if(interval.isValid() /*&& interval.size()<200*/)
	{
		for(int64_t j = interval.lower; j <= interval.upper && (j-interval.lower<m_maxIDs); ++j)
		{
		    int64_t SAindex = pSSA->calcSA(j,pBWT).getID();
            m_SAindicesSet2.insert(SAindex);
		}
	}
}

void NameSet::addReadIDAndContigID(std::string seed, ThreadSafeListVector* tslv, Vertex* pVertex, ReadOnContig roc)
{
	BWTInterval interval = BWTAlgorithms::findInterval(pBWT, seed);
	
	//std::cout << interval.size() << "\n";
	if(interval.isValid()) 
	{
		for(int64_t j = interval.lower; j <= interval.upper && (j-interval.lower<m_maxIDs); ++j)
		{
			//add read IDs into hashset
		    int64_t SAindex = pSSA->calcSA(j, pBWT).getID();
			m_SAindicesSet1.insert( SAindex);
			
			//Also add this pVertex into list of read SAindex
			tslv->at(SAindex).setLock();
			tslv->at(SAindex).push_back( std::make_pair(pVertex,roc) );
			tslv->at(SAindex).releaseLock();
		}
	}
}

std::vector<int64_t> NameSet::getReadIDs()
{
	std::vector<int64_t> ReadIDs;
	HashSet<int64_t>::const_iterator cit = m_SAindicesSet1.begin();
	for (; cit != m_SAindicesSet1.end(); ++cit )
		ReadIDs.push_back(*cit);
		
	return ReadIDs;
}

void NameSet::getAnotherReadIDs(std::vector<int64_t>& anotherIDs)
{
	HashSet<int64_t>::const_iterator cit = m_SAindicesSet1.begin();
	for (; cit != m_SAindicesSet1.end(); ++cit )
	{
		if(*cit % 2 ==0)
			anotherIDs.push_back(*cit+1);
		else
			anotherIDs.push_back(*cit-1);
	}
	// for (size_t i=0; i<m_SAindices.size(); i++)
	// {
		// if(m_SAindices[i]%2==0)
			// anotherIDs.push_back(m_SAindices[i]+1);
		// else
			// anotherIDs.push_back(m_SAindices[i]-1);
	// }
}

bool NameSet::exist (int64_t idx)
{
	HashSet<int64_t>::const_iterator cit = m_SAindicesSet2.find (idx);
	if(cit==m_SAindicesSet2.end())
		return false;
	else
		return true;
}
