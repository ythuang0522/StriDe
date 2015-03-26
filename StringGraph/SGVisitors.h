//-----------------------------------------------
// Originally written by Jared Simpson (js18@sanger.ac.uk)
// Old visitors were revised by YTH while adding new visitors specific for the hybrid graph
// Released under the GPL
//-----------------------------------------------
//
// SGVisitors - Algorithms that visit
// each vertex in the graph and perform some
// operation
//
#include "SGAlgorithms.h"
#include "SGUtil.h"

#include "BWTIndexSet.h"
#include "BWTAlgorithms.h"

#ifndef SGVISITORS_H
#define SGVISITORS_H


//a hashtable rapper for query and conversion of PE read IDs
//linked lists of contig IDs mapped by this read
typedef std::vector<ThreadSafeList<std::pair<Vertex*, ReadOnContig> > > ThreadSafeListVector;

class NameSet
{
public:

	NameSet(const BWT* bwt, const SampledSuffixArray* ssa, int64_t maxIDs=200)
	:pBWT(bwt),pSSA(ssa),m_maxIDs(maxIDs)
	{
	}

	//Direct SA index implementation by YTH
	void addFirstReadIDs(std::string seed);
	void addSecondReadIDs(std::string seed);
	void addReadIDAndContigID(std::string seed, ThreadSafeListVector* tslv, Vertex* pVertex, ReadOnContig roc);
	std::vector<int64_t> getReadIDs(); 
	void getAnotherReadIDs(std::vector<int64_t>& anotherIDs);
	bool exist (int64_t idx);
	size_t sizeOfFirstReadIDs(){return m_SAindicesSet1.size();};
	size_t sizeOfSecondReadIDs(){return m_SAindicesSet2.size();};

private:

	const BWT* pBWT ;
	const SampledSuffixArray* pSSA ;

	// std::vector<int64_t> m_SAindices;
	HashSet<int64_t> m_SAindicesSet1;
	HashSet<int64_t> m_SAindicesSet2;
	int64_t m_maxIDs;
};

// Visit each node, writing it to a file as a fasta record
struct SGFastaVisitor
{
    // constructor
    SGFastaVisitor(std::string filename,BWT* bwt=NULL,size_t k=0) : m_fileHandle(filename.c_str()),pBWT(bwt),kmerLength(k) {}
    ~SGFastaVisitor() { m_fileHandle.close(); }

    // functions
    void previsit(StringGraph* /*pGraph*/) {}
    bool visit(StringGraph* pGraph, Vertex* pVertex);
    void postvisit(StringGraph* /*pGraph*/) {}

    // data
    std::ofstream m_fileHandle;
	BWT* pBWT ;
	size_t kmerLength;
};

// Run the Myers transitive reduction algorithm on each node
struct SGTransitiveReductionVisitor
{
    SGTransitiveReductionVisitor() {}
    void previsit(StringGraph* pGraph);
    bool visit(StringGraph* pGraph, Vertex* pVertex);
    void postvisit(StringGraph*);

    int marked_verts;
    int marked_edges;
};

// Remove identical vertices from the graph
struct SGIdenticalRemoveVisitor
{
    SGIdenticalRemoveVisitor() {}
    void previsit(StringGraph* pGraph);
    bool visit(StringGraph* pGraph, Vertex* pVertex);
    void postvisit(StringGraph* pGraph);
    int count;
};

// Remove contained vertices from the graph
struct SGContainRemoveVisitor
{
    SGContainRemoveVisitor() {}
    void previsit(StringGraph* pGraph);
    bool visit(StringGraph* pGraph, Vertex* pVertex);
    void postvisit(StringGraph* pGraph);

	std::ostream* pWriter ;
};


// Detects and removes small "tip" vertices from the graph
// when they are less than minLength in size
struct SGTrimVisitor
{
    SGTrimVisitor(std::string filename="" , size_t minLength=400)
	: m_filename(filename), m_minLength(minLength)  {}
	~SGTrimVisitor(){
		if(!m_filename.empty())	m_fileHandle.close();
	}

    void previsit(StringGraph* pGraph);
    bool visit(StringGraph* pGraph, Vertex* pVertex);
    void postvisit(StringGraph*);

	std::string m_filename;
	std::ofstream m_fileHandle;
    size_t m_minLength;
    int num_island;
    int num_terminal;

};

// Detect and remove duplicate edges
struct SGDuplicateVisitor
{
    SGDuplicateVisitor(bool silent = false) : m_bSilent(silent) {}
    void previsit(StringGraph*);
    bool visit(StringGraph* pGraph, Vertex* pVertex);
    void postvisit(StringGraph*);

    bool m_hasDuplicate;
    bool m_bSilent;
};

// Remove the edges of super-repetitive vertices in the graph
struct SGSuperRepeatVisitor
{
    SGSuperRepeatVisitor() {}
    void previsit(StringGraph* pGraph);
    bool visit(StringGraph* pGraph, Vertex* pVertex);
    void postvisit(StringGraph*);

    size_t m_num_superrepeats;
};

// Smooth out variation in the graph
struct SGSmoothingVisitor
{
    SGSmoothingVisitor(int maxIndelLength, BWT* pBWT, bool bIsGapPrecent=false) : m_numRemovedTotal(0),
                                             m_maxIndelLength(maxIndelLength), m_pBWT(pBWT), m_bIsGapPrecent(bIsGapPrecent)
											 {}
	~SGSmoothingVisitor(){ 
	}

    void previsit(StringGraph* pGraph);
    bool visit(StringGraph* pGraph, Vertex* pVertex);
    void postvisit(StringGraph*);

    int m_simpleBubblesRemoved;
    int m_complexBubblesRemoved;
    int m_numRemovedTotal;

    int m_maxIndelLength;
	BWT* m_pBWT;
	bool m_bIsGapPrecent;

};

// Compile summary statistics for the graph
struct SGGraphStatsVisitor
{
    SGGraphStatsVisitor() {}
    void previsit(StringGraph* pGraph);
    bool visit(StringGraph* pGraph, Vertex* pVertex);
    void postvisit(StringGraph*);

	int num_straight;
    int num_terminal;
    int num_island;
    int num_monobranch;
    int num_dibranch;
    int num_simple;
    int num_edges;
    int num_vertex;
    size_t sum_edgeLen;
};


/*****************************************************************/
/*                       New Visitor   		                     */
/*             write by CCU CSIE LAB 405 YTH and Tatsuki                 */
/*****************************************************************/


// Visit each node, writing it to a file as a fasta record
struct SGFastaErosionVisitor
{
    // constructor
    SGFastaErosionVisitor(BWT* bwt,size_t k,size_t t,size_t m=500,size_t e=1)
	: pBWT(bwt),kmerLength(k),threshold(t),min_island(m), erosion(e) {}
    ~SGFastaErosionVisitor() {  }

    // functions
    void previsit(StringGraph* /*pGraph*/) {}
    bool visit(StringGraph* pGraph, Vertex* pVertex);
    void postvisit(StringGraph* /*pGraph*/) {}

    // data
	BWT* pBWT ;
	size_t kmerLength;
	size_t threshold;
	size_t min_island;
	size_t erosion;

};



////////////////////////////////
// Remove illegal kmer edge  //
///////////////////////////////


struct SGRemoveIllegalKmerEdgeVisitor
{
	SGRemoveIllegalKmerEdgeVisitor (BWT* bwt , size_t kmerLength, size_t threshold, size_t minOverlapLength, std::string filename="")
	: pBWT(bwt), m_kmerLength(kmerLength), m_threshold(threshold), m_minOverlapLength(minOverlapLength), m_filename(filename) {}
	~SGRemoveIllegalKmerEdgeVisitor(){ 
		if(!m_filename.empty())
			m_fileHandle.close(); 
	}

	void previsit(StringGraph* pGraph);
    bool visit(StringGraph* , Vertex* pVertex);
    void postvisit(StringGraph* pGraph);

	BWT* pBWT ;
	size_t m_kmerLength ;
	size_t m_threshold ;
	size_t m_minOverlapLength;
	std::string m_filename;
	std::ofstream m_fileHandle;

};



//////////////////////////////
// Both Short Edges Remove  //
/////////////////////////////


struct  SGBothShortEdgesRemoveVisitor
{
	SGBothShortEdgesRemoveVisitor (size_t vl,size_t ol,BWT* bwt=NULL,size_t kl=0,float t=0)
	: vertexLength(vl),overlapLength(ol), pBWT(bwt),kmerLength(kl),threshold(t)
	{
		if (pBWT !=NULL)
			assert (kmerLength>0);
	}

	SGBothShortEdgesRemoveVisitor (size_t vl,size_t ol,std::string filename,BWT* bwt=NULL,size_t kl=0,float t=0)
	: vertexLength(vl),overlapLength(ol), m_fileHandle(filename.c_str()),pBWT(bwt),kmerLength(kl),threshold(t)
	{
		if (pBWT !=NULL)
			assert (kmerLength>0);
	}
	~SGBothShortEdgesRemoveVisitor(){
        if(m_fileHandle!=NULL)
            m_fileHandle.close();
    }

	void previsit(StringGraph* pGraph);
    bool visit(StringGraph* , Vertex* pVertex);
    void postvisit(StringGraph* pGraph);


	size_t vertexLength;
	size_t overlapLength;
	std::ofstream m_fileHandle;
	BWT* pBWT ;
	size_t kmerLength ;
	float threshold ;
};


/////////////////////////
// Remove low overlap //
////////////////////////


struct SGRemoveByOverlapLenDiffVisitor
{
	SGRemoveByOverlapLenDiffVisitor (size_t vertex_size=3000, size_t min_overlap=80, size_t max_overlapdiff=0, bool islandProtect=true):
	    m_min_vertex_size(vertex_size), m_min_overlap(min_overlap), m_max_overlapdiff(max_overlapdiff), m_islandProtect(islandProtect)
	    {}
	void previsit(StringGraph* pGraph);
    bool visit(StringGraph* , Vertex* pVertex);
    void postvisit(StringGraph* pGraph);

    size_t m_min_vertex_size ;
	size_t m_min_overlap ;
    size_t m_max_overlapdiff ;
    bool m_islandProtect;
};


/////////////////////////////////////
//  Sweep Low Overlap Ratio Edge  //
////////////////////////////////////

struct SGLowOverlapRatioEdgeSweepVisitor
{
	SGLowOverlapRatioEdgeSweepVisitor(size_t vertex_size=0, double minOverlapRatio=0.8,size_t maxMatchLength=0)
	: m_min_vertex_size(vertex_size), m_overlapRatio(minOverlapRatio),m_matchLength(maxMatchLength){}

	SGLowOverlapRatioEdgeSweepVisitor(std::string filename,double minOverlapRatio,size_t maxMatchLength=0)
	: m_fileHandle(filename.c_str()),m_overlapRatio(minOverlapRatio),m_matchLength(maxMatchLength){}

	~SGLowOverlapRatioEdgeSweepVisitor () {
	    if(m_fileHandle!=NULL)
            m_fileHandle.close();
	}

    void previsit(StringGraph* pGraph);
    bool visit(StringGraph* , Vertex* pVertex);
    void postvisit(StringGraph* pGraph);

	std::ofstream m_fileHandle;
	size_t m_min_vertex_size;
	double m_overlapRatio;
	size_t m_matchLength ;
};


////////////////////////////////////////////////////
// Remove Low Overlap Ratio Edge on Short Vertex  //
////////////////////////////////////////////////////
struct SGLowOverlapRatioEdgeRemoveVisitor
{
    SGLowOverlapRatioEdgeRemoveVisitor(size_t maxVertexLength=0, double minDiffRatio=0.2, bool skipLargeVertex=true)
	: m_vetexLength(maxVertexLength), m_diffRatio(minDiffRatio), m_skipLargeVertex(skipLargeVertex){}

	SGLowOverlapRatioEdgeRemoveVisitor(std::string filename,size_t maxVertexLength=0,double minDiffRatio=0.2)
	: m_fileHandle(filename.c_str()),m_vetexLength(maxVertexLength),m_diffRatio(minDiffRatio){}

	~SGLowOverlapRatioEdgeRemoveVisitor () {
    if(m_fileHandle!=NULL)
        m_fileHandle.close();
    }

    void previsit(StringGraph* pGraph);
    bool visit(StringGraph* , Vertex* pVertex);
    void postvisit(StringGraph* pGraph);

	std::ofstream m_fileHandle;
	size_t m_vetexLength ;
	double m_diffRatio;
	bool m_skipLargeVertex ;
};


struct SGSubGraphVisitor
{
    SGSubGraphVisitor() {}
    void previsit(StringGraph* pGraph);
    bool visit(StringGraph* pGraph, Vertex* pVertex);
    void postvisit(StringGraph*);

	int num_straight;
    int num_terminal;
    int num_island;
    int num_monobranch;
    int num_dibranch;
    int num_simple;
    int num_edges;
    int num_vertex;
    size_t sum_edgeLen;
};

struct SGRemoveEdgeByPEVisitor
{
    SGRemoveEdgeByPEVisitor(BWTIndexSet indices, size_t insertSize, size_t kmerSize, size_t threshold)
	:m_indices(indices), m_insertSize(insertSize), m_kmerSize(kmerSize), m_minPEcount(threshold){}

	void previsit(StringGraph* pGraph);
    bool visit(StringGraph* pGraph, Vertex* pVertex);
    void postvisit(StringGraph* pGraph);

	bool addReadIDsAtPos(NameSet& pVReadID, std::string& pVertexSeq, int overlapBoundaryPos);

	BWTIndexSet m_indices;

    size_t m_insertSize;
	size_t m_kmerSize;
	size_t m_minPEcount;
	size_t m_edgecount;

	size_t m_repeatKmerCutoff;

};


//Collect PE reads at ends of each island/tip
//Store PE read IDs into NameSet hashtable
struct SGIslandCollectVisitor
{
    SGIslandCollectVisitor(ThreadSafeListVector* tslv, BWTIndexSet indices, size_t insertSize, size_t kmerSize, size_t islandSize)
	:m_tslv(tslv), m_indices(indices),m_insertSize(insertSize),m_kmerSize(kmerSize),m_minIslandSize(islandSize){}

	void previsit(StringGraph* pGraph);
    bool visit(StringGraph* pGraph, Vertex* pVertex);
    void postvisit(StringGraph* pGraph);

	ThreadSafeListVector* m_tslv;
	BWTIndexSet m_indices;

    size_t m_insertSize;
	size_t m_kmerSize;
	size_t m_minIslandSize;
	size_t m_islandcount;
    std::vector<std::pair<std::string, std::vector<NameSet> > > *m_islandReadIDs;
	KmerDistribution m_kd;
	size_t m_repeatKmerCutoff;
};

//Join islands/tips with PE link support via SAI walk
//SAI walk is revised to walk through high-error gaps
struct SGJoinIslandVisitor
{
    SGJoinIslandVisitor(size_t SAISearchDepth, size_t SAISearchLeaves, size_t kmer, size_t islandSize, ThreadSafeListVector* tslv, BWTIndexSet indices,
            size_t minPEcount=5)
	:m_SAISearchDepth(SAISearchDepth),m_SAISearchLeaves(SAISearchLeaves),
	m_kmer(kmer),m_minIslandSize(islandSize), 
	m_tslv(tslv), m_indices(indices),
	m_minPEcount(minPEcount)
	{
		m_numOfIterations=2;
	}

	void previsit(StringGraph* pGraph);
    bool visit(StringGraph* pGraph, Vertex* pVertex);
    void postvisit(StringGraph* pGraph);
	
	void findNeighborWithPESupport(Vertex* pVertex, size_t islandDir, SparseHashMap<VertexID, size_t*, StringHasher>& pWIDs);
	void updateExtendedVertex(Vertex* pVertex,  std::string& newStr, EdgeDir dir);
	int64_t getAnotherID(int64_t idx);
		
	size_t m_SAISearchDepth;
    size_t m_SAISearchLeaves;
    size_t m_kmer;	
	size_t m_minIslandSize;

	ThreadSafeListVector* m_tslv;
	BWTIndexSet m_indices;

	size_t m_islandcount;
	size_t m_minPEcount;	
	size_t m_numOfIterations;
	KmerDistribution m_kd;
};

#endif
