//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// SGUtils - Data structures/Functions related
// to building and manipulating string graphs
//
#include "SGUtil.h"
#include "SeqReader.h"
#include "SGAlgorithms.h"
#include "SGVisitors.h"
#include "../StriDe/SGACommon.h"

StringGraph* SGUtil::loadASQG(const StringVector & filenameList, const unsigned int minOverlap, bool allowContainments , size_t maxEdges ,GraphColor c)
{
	// Initialize graph
	StringGraph* pGraph = new StringGraph;
	int stage = 0;
	int line = 0;
	for(size_t i=0 ; i<filenameList.size() ;i++)
	{		
		
		std::istream* pReader = createReader(filenameList[i]);
		std::string recordLine;
		while(getline(*pReader, recordLine))
		{   
			ASQG::RecordType rt = ASQG::getRecordType(recordLine);
			switch(rt)
			{
			case ASQG::RT_HEADER:
				{
					if(stage != 0)
					{
						std::cerr << "Error: Unexpected header record found at line " << line << "\n";
						exit(EXIT_FAILURE);
					}

					ASQG::HeaderRecord headerRecord(recordLine);
					const SQG::IntTag& overlapTag = headerRecord.getOverlapTag();
					if(overlapTag.isInitialized())
					pGraph->setMinOverlap(overlapTag.get());
					else
					pGraph->setMinOverlap(0);

					const SQG::FloatTag& errorRateTag = headerRecord.getErrorRateTag();
					if(errorRateTag.isInitialized())
					pGraph->setErrorRate(errorRateTag.get());
					
					const SQG::IntTag& containmentTag = headerRecord.getContainmentTag();
					if(containmentTag.isInitialized())
						pGraph->setContainmentFlag(containmentTag.get());
					else
						pGraph->setContainmentFlag(false); // conservatively assume containments are present

					const SQG::IntTag& transitiveTag = headerRecord.getTransitiveTag();
					if(!transitiveTag.isInitialized())
					{
						std::cerr << "Warning: ASQG does not have transitive tag\n";
						pGraph->setTransitiveFlag(true);
					}
					else
					{
						pGraph->setTransitiveFlag(transitiveTag.get());
					}

					break;
				}
			case ASQG::RT_VERTEX:
				{
					// progress the stage if we are done the header
					if(stage == 0)
					stage = 1;

					if(stage != 1)
					{
						std::cerr << "Error: Unexpected vertex record found at line " << line << "\n";
						exit(EXIT_FAILURE);
					}

					ASQG::VertexRecord vertexRecord(recordLine);
					const SQG::IntTag& ssTag = vertexRecord.getSubstringTag();

					Vertex* pVertex = new(pGraph->getVertexAllocator()) Vertex(vertexRecord.getID(), vertexRecord.getSeq());
					if(ssTag.isInitialized() && ssTag.get() == 1)
					{
						// Vertex is a substring of some other vertex, mark it as contained
						pVertex->setContained(true);
						pGraph->setContainmentFlag(true);
					}
					pGraph->addVertex(pVertex);
					break;
				}
			case ASQG::RT_EDGE:
				{
					if(stage == 1)
					stage = 2;
					
					if(stage != 2)
					{
						std::cerr << "Error: Unexpected edge record found at line " << line << "\n";
						exit(EXIT_FAILURE);
					}

					ASQG::EdgeRecord edgeRecord(recordLine);
					const Overlap& ovr = edgeRecord.getOverlap();

	
					// Add the edge to the graph
					if(ovr.match.getMinOverlapLength() >= (int)minOverlap)
						SGAlgorithms::createEdgesFromOverlap(pGraph, ovr, allowContainments, maxEdges,c);
					break;
				}
			}
			++line;
		}
		delete pReader;
	}	
	
	// Completely delete the edges for all nodes that were marked as super-repetitive in the graph
	//SGSuperRepeatVisitor superRepeatVisitor;
	//pGraph->visitP(superRepeatVisitor);

	// Remove any duplicate edges
	SGDuplicateVisitor dupVisit;
	pGraph->visit(dupVisit);

	return pGraph;
}

StringGraph* SGUtil::loadASQGVertex(const std::string& filename, const unsigned int minOverlap, bool allowContainments, size_t maxEdges)
{
	// Initialize graph
	StringGraph* pGraph = new StringGraph;
	std::istream* pReader = createReader(filename);
	int stage = 0;
	int line = 0;
	std::string recordLine;

	while(getline(*pReader, recordLine))
	{ 	
		ASQG::RecordType rt = ASQG::getRecordType(recordLine);
		switch(rt)
		{
			case ASQG::RT_HEADER:
			{
				if(stage != 0)
				{
					std::cerr << "Error: Unexpected header record found at line " << line << "\n";
					exit(EXIT_FAILURE);
				}

				ASQG::HeaderRecord headerRecord(recordLine);
				const SQG::IntTag& overlapTag = headerRecord.getOverlapTag();
				if(overlapTag.isInitialized())
				pGraph->setMinOverlap(overlapTag.get());
				else
				pGraph->setMinOverlap(0);

				const SQG::FloatTag& errorRateTag = headerRecord.getErrorRateTag();
				if(errorRateTag.isInitialized())
				pGraph->setErrorRate(errorRateTag.get());
				
				const SQG::IntTag& containmentTag = headerRecord.getContainmentTag();
				if(containmentTag.isInitialized())
				pGraph->setContainmentFlag(containmentTag.get());
				else
				pGraph->setContainmentFlag(true); // conservatively assume containments are present

				const SQG::IntTag& transitiveTag = headerRecord.getTransitiveTag();
				if(!transitiveTag.isInitialized())
				{
					std::cerr << "Warning: ASQG does not have transitive tag\n";
					pGraph->setTransitiveFlag(true);
				}
				else
				{
					pGraph->setTransitiveFlag(transitiveTag.get());
				}

				break;
			}
			case ASQG::RT_VERTEX:
			{
				// progress the stage if we are done the header
				if(stage == 0)
				stage = 1;

				if(stage != 1)
				{
					std::cerr << "Error: Unexpected vertex record found at line " << line << "\n";
					exit(EXIT_FAILURE);
				}

				ASQG::VertexRecord vertexRecord(recordLine);
				const SQG::IntTag& ssTag = vertexRecord.getSubstringTag();

				Vertex* pVertex = new(pGraph->getVertexAllocator()) Vertex(vertexRecord.getID(), vertexRecord.getSeq());
				if(ssTag.isInitialized() && ssTag.get() == 1)
				{
					// Vertex is a substring of some other vertex, mark it as contained
					pVertex->setContained(true);
					pGraph->setContainmentFlag(true);
				}
				pGraph->addVertex(pVertex);
				break;
			}
			case ASQG::RT_EDGE:
			{
				if(stage == 1)
				stage = 2;
				
				if(stage != 2)
				{
					std::cerr << "Error: Unexpected edge record found at line " << line << "\n";
					exit(EXIT_FAILURE);
				}

				ASQG::EdgeRecord edgeRecord(recordLine);
				const Overlap& ovr = edgeRecord.getOverlap();


				// Add the edge to the graph
				if(ovr.match.getMinOverlapLength() >= (int)minOverlap)
					SGAlgorithms::createEdgesFromOverlap(pGraph, ovr, allowContainments, maxEdges);
				break;
			}
		}//end of switch
		++line;
	}//end of each line in the main asqg

	delete pReader;
	return pGraph;
}

StringGraph* SGUtil::loadASQGEdge(std::string ASQGFileName, const unsigned int minOverlap, bool allowContainments, size_t maxEdges, StringGraph* pGraph)
{
	std::string edgeFilePrefix = stripFilename(ASQGFileName);
	std::vector<std::istream*> EdgeFileVec;

	//search for the edges files named with xxxx-thread??.edges.gz, if existed
	std::stringstream ss;
	size_t fileIdx=0;
	ss << edgeFilePrefix << "-thread" << fileIdx << HITS_EXT << GZIP_EXT;
	std::string edgefile = ss.str();
	std::istream* pistream = createReader(edgefile);
	while(pistream->peek()!=std::istream::traits_type::eof())
	{
		std::cout << edgefile << std::endl;
		EdgeFileVec.push_back(pistream);
		//search for next edge file
		ss.str("");
		fileIdx++;
		ss << edgeFilePrefix << "-thread" << fileIdx << HITS_EXT << GZIP_EXT;
		edgefile = ss.str();
		pistream = createReader(edgefile);		
	}

	int line = 0;

	#pragma omp parallel for
	for(size_t i=0 ; i< EdgeFileVec.size() ;i++)
	{
		std::string recordLine;
		while(getline(*EdgeFileVec[i],recordLine))
		{
			ASQG::RecordType rt = ASQG::getRecordType(recordLine);
			if (rt != ASQG::RT_EDGE)
			{	
				std::cerr << "Error: Unexpected edge record found at line " << line << "\n";
				exit(EXIT_FAILURE);
			}
			ASQG::EdgeRecord edgeRecord(recordLine);
			const Overlap& ovr = edgeRecord.getOverlap();
			
			// Add the edge to the graph
			if(ovr.match.getMinOverlapLength() >= (int)minOverlap)
				SGAlgorithms::createEdgesFromOverlap(pGraph, ovr, allowContainments, maxEdges);
			
			++line;
		}
		delete EdgeFileVec[i];
	}

	// Completely delete the edges for all nodes that were marked as super-repetitive in the graph
	//SGSuperRepeatVisitor superRepeatVisitor;
	//pGraph->visitP(superRepeatVisitor);

	// Remove any duplicate edges
	SGDuplicateVisitor dupVisit;
	pGraph->visit(dupVisit);

	return pGraph;
}

StringGraph*  SGUtil::loadASQG_Parallel(const StringVector & filenameList, const unsigned int minOverlap, 
bool allowContainments , size_t maxEdges ,GraphColor c)
{
	//RED_EDGE for dupEdge
	assert (c != GC_RED);
	
	StringGraph* pGraph = new StringGraph;
	size_t fileNumber = filenameList.size();
	std::vector<StringVector> decodeFiles(fileNumber);

	std::cout << "#Parallel decoding" << std::endl;
	#pragma omp parallel for num_threads(fileNumber)
	for ( size_t i = 0; i < fileNumber ;i ++)
	{
		std::istream* pReader = createReader(filenameList[i]);
		std::string recordLine;
		while(getline(*pReader, recordLine)) 
			decodeFiles[i].push_back(recordLine);
		std::cout << filenameList[i] << " is decoded." << std::endl;
		delete pReader;
		
	}

	std::cout << "#Serial loading." << std::endl;
	for (size_t i=0;i<fileNumber;i++)
	{
		for (size_t j=0;j<decodeFiles[i].size();j++)
		{
			ASQG::RecordType rt = ASQG::getRecordType(decodeFiles[i][j]);
			switch(rt)
			{
			case ASQG::RT_HEADER:
				{
					ASQG::HeaderRecord headerRecord(decodeFiles[i][j]);
					const SQG::IntTag& overlapTag = headerRecord.getOverlapTag();
					if(overlapTag.isInitialized())
						pGraph->setMinOverlap(overlapTag.get());
					else
						pGraph->setMinOverlap(0);

					const SQG::FloatTag& errorRateTag = headerRecord.getErrorRateTag();
					if(errorRateTag.isInitialized())
					pGraph->setErrorRate(errorRateTag.get());
					
					const SQG::IntTag& containmentTag = headerRecord.getContainmentTag();
					if(containmentTag.isInitialized())
					pGraph->setContainmentFlag(containmentTag.get());
					else
					pGraph->setContainmentFlag(true); // conservatively assume containments are present

					const SQG::IntTag& transitiveTag = headerRecord.getTransitiveTag();
					if(!transitiveTag.isInitialized())
					{
						std::cerr << "Warning: ASQG does not have transitive tag\n";
						pGraph->setTransitiveFlag(true);
					}
					else
					{
						pGraph->setTransitiveFlag(transitiveTag.get());
					}
					break;
				}
			case ASQG::RT_VERTEX:
				{
					// progress the stage if we are done the header
					ASQG::VertexRecord vertexRecord(decodeFiles[i][j]);
					const SQG::IntTag& ssTag = vertexRecord.getSubstringTag();

					Vertex* pVertex = new(pGraph->getVertexAllocator()) Vertex(vertexRecord.getID(), vertexRecord.getSeq());
					if(ssTag.isInitialized() && ssTag.get() == 1)
					{
						// Vertex is a substring of some other vertex, mark it as contained
						pVertex->setContained(true);
						pGraph->setContainmentFlag(true);
					}
					pGraph->addVertex(pVertex);
					break;
				}
			case ASQG::RT_EDGE:
				{
					ASQG::EdgeRecord edgeRecord(decodeFiles[i][j]);
					const Overlap& ovr = edgeRecord.getOverlap();

						// Add the edge to the graph
					if(ovr.match.getMinOverlapLength() >= (int)minOverlap)
							SGAlgorithms::createEdgesFromOverlap(pGraph, ovr, allowContainments, maxEdges,c);
					break;
				}
			
			}
		
		}
		std::cout << filenameList[i] << " is finish." << std::endl;
	}
	
	

		
	// Completely delete the edges for all nodes that were marked as super-repetitive in the graph
	//SGSuperRepeatVisitor superRepeatVisitor;
	//pGraph->visitP(superRepeatVisitor);

	// Remove any duplicate edges
	SGDuplicateVisitor dupVisit;
	pGraph->visit(dupVisit);

	
	//SGGraphStatsVisitor statsVisit;
	//pGraph->visit(statsVisit);


	return pGraph;


}



StringGraph* SGUtil::loadASQG_EDGE_Parallel(const StringVector & filenameList, const unsigned int minOverlap, bool allowContainments, size_t maxEdges ,GraphColor c , StringGraph* pGraph)
{
	//RED_EDGE for dupEdge
	assert (c != GC_RED);
	
	size_t fileNumber = filenameList.size();
	std::vector<StringVector> decodeFiles(fileNumber);
	
	std::cout << "#Parallel decoding" << std::endl;
	#pragma omp parallel for num_threads(fileNumber)
	for ( size_t i = 0; i < fileNumber ;i ++)
	{
		std::istream* pReader = createReader(filenameList[i]);
		std::string recordLine;
		while(getline(*pReader, recordLine)) 
			decodeFiles[i].push_back(recordLine);
		std::cout << filenameList[i] << " is decoded." << std::endl;
		delete pReader;
		
	}
	
	std::cout << "#Serial loading." << std::endl;
	for (size_t i=0;i<fileNumber;i++)
	{
		for (size_t j=0;j<decodeFiles[i].size();j++)
		{
			ASQG::RecordType rt = ASQG::getRecordType(decodeFiles[i][j]);
			if (rt != ASQG::RT_EDGE)
			{	
				std::cerr << "Error: Unexpected edge record found at " << filenameList[i] << "\n";
				exit(EXIT_FAILURE);
			}
			ASQG::EdgeRecord edgeRecord(decodeFiles[i][j]);
			const Overlap& ovr = edgeRecord.getOverlap();
			
			// Add the edge to the graph
			if(ovr.match.getMinOverlapLength() >= (int)minOverlap)
				SGAlgorithms::createEdgesFromOverlap(pGraph, ovr, allowContainments, maxEdges,c);
				
		}
		std::cout << filenameList[i] << " is finish." << std::endl;
	}
	
	
	// Completely delete the edges for all nodes that were marked as super-repetitive in the graph
	//SGSuperRepeatVisitor superRepeatVisitor;
	//pGraph->visitP(superRepeatVisitor);

	// Remove any duplicate edges
	SGDuplicateVisitor dupVisit;
	pGraph->visit(dupVisit);

	
	//SGGraphStatsVisitor statsVisit;
	//pGraph->visit(statsVisit);

	return pGraph;
}


StringGraph* SGUtil::loadASQG(const std::string& filename, const unsigned int minOverlap, bool allowContainments, size_t maxEdges)
{
	// Initialize graph
	StringGraph* pGraph = new StringGraph;

	std::istream* pReader = createReader(filename);

	int stage = 0;
	int line = 0;
	std::string recordLine;

	while(getline(*pReader, recordLine))
	{ 	
		ASQG::RecordType rt = ASQG::getRecordType(recordLine);
		switch(rt)
		{
		case ASQG::RT_HEADER:
			{
				if(stage != 0)
				{
					std::cerr << "Error: Unexpected header record found at line " << line << "\n";
					exit(EXIT_FAILURE);
				}

				ASQG::HeaderRecord headerRecord(recordLine);
				const SQG::IntTag& overlapTag = headerRecord.getOverlapTag();
				if(overlapTag.isInitialized())
				pGraph->setMinOverlap(overlapTag.get());
				else
				pGraph->setMinOverlap(0);

				const SQG::FloatTag& errorRateTag = headerRecord.getErrorRateTag();
				if(errorRateTag.isInitialized())
				pGraph->setErrorRate(errorRateTag.get());
				
				const SQG::IntTag& containmentTag = headerRecord.getContainmentTag();
				if(containmentTag.isInitialized())
				pGraph->setContainmentFlag(containmentTag.get());
				else
				pGraph->setContainmentFlag(true); // conservatively assume containments are present

				const SQG::IntTag& transitiveTag = headerRecord.getTransitiveTag();
				if(!transitiveTag.isInitialized())
				{
					std::cerr << "Warning: ASQG does not have transitive tag\n";
					pGraph->setTransitiveFlag(true);
				}
				else
				{
					pGraph->setTransitiveFlag(transitiveTag.get());
				}

				break;
			}
		case ASQG::RT_VERTEX:
			{
				// progress the stage if we are done the header
				if(stage == 0)
				stage = 1;

				if(stage != 1)
				{
					std::cerr << "Error: Unexpected vertex record found at line " << line << "\n";
					exit(EXIT_FAILURE);
				}

				ASQG::VertexRecord vertexRecord(recordLine);
				const SQG::IntTag& ssTag = vertexRecord.getSubstringTag();

				Vertex* pVertex = new(pGraph->getVertexAllocator()) Vertex(vertexRecord.getID(), vertexRecord.getSeq());
				if(ssTag.isInitialized() && ssTag.get() == 1)
				{
					// Vertex is a substring of some other vertex, mark it as contained
					pVertex->setContained(true);
					pGraph->setContainmentFlag(true);
				}
				pGraph->addVertex(pVertex);
				break;
			}
		case ASQG::RT_EDGE:
			{
				if(stage == 1)
				stage = 2;
				
				if(stage != 2)
				{
					std::cerr << "Error: Unexpected edge record found at line " << line << "\n";
					exit(EXIT_FAILURE);
				}

				ASQG::EdgeRecord edgeRecord(recordLine);
				const Overlap& ovr = edgeRecord.getOverlap();

				// Add the edge to the graph
				if(ovr.match.getMinOverlapLength() >= (int)minOverlap)
					SGAlgorithms::createEdgesFromOverlap(pGraph, ovr, allowContainments, maxEdges);
				break;
			}
		}
		++line;
	}
	
	// Completely delete the edges for all nodes that were marked as super-repetitive in the graph
	//SGSuperRepeatVisitor superRepeatVisitor;
	//pGraph->visitP(superRepeatVisitor);

	
	// Remove any duplicate edges
	SGDuplicateVisitor dupVisit;
	pGraph->visit(dupVisit);

	// Remove identical vertices
	// This is much cheaper to do than remove via
	// SGContainRemove as no remodelling needs to occur
	/*	
	SGIdenticalRemoveVisitor irv;
	pGraph->visit(irv);
	*/
	delete pReader;
	return pGraph;
}

// Load a graph (with no edges) from a fasta file
StringGraph* SGUtil::loadFASTA(const std::string& filename)
{
	StringGraph* pGraph = new StringGraph;
	SeqReader reader(filename);
	SeqRecord record;

	while(reader.get(record))
	{
		Vertex* pVertex = new(pGraph->getVertexAllocator()) Vertex(record.id, record.seq.toString());
		pGraph->addVertex(pVertex);
	}
	return pGraph;
}
