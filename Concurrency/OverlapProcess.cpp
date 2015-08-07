///-----------------------------------------------
// Copyright 2015 National Chung Cheng University
// Written by Yao-Ting Huang
// Revised from Jared Simpson's overlapprocess
// Implement distributed graphs for parallelization
// Released under the GPL
//-----------------------------------------------
//
// OverlapProcess - write overlap blocks into edges of graphs
//
#include "OverlapProcess.h"
#include "../SQG/ASQG.h"

//
//
//
OverlapProcess::OverlapProcess(const std::string& outFile, 
                               const OverlapAlgorithm* pOverlapper, 
                               int minOverlap) : m_pOverlapper(pOverlapper), 
                                                 m_minOverlap(minOverlap)
{
    m_pWriter = createWriter(outFile);
}

//
OverlapProcess::~OverlapProcess()
{
    delete m_pWriter;
}

//
OverlapResult OverlapProcess::process(const SequenceWorkItem& workItem)
{
	//compute overlap of workItem.read with results stored in m_blockList, where each block stores SA intervals of overlapping reads
    OverlapResult result = m_pOverlapper->overlapRead(workItem.read, m_minOverlap, &m_blockList);
	
	// Give up substring overlap
	if(result.isSubstring) 
	{	
		// std::cout << "substring after submaximal removal\n";
		return result;
	}

	//Convert list of overlap blocks into Edges
    OverlapBlockList::iterator it = m_blockList.begin();
	for(; it != m_blockList.end(); it++)
    {
        // Read one block
        OverlapBlock record = *it;
		
		// skip target substr
		if(record.isTargetSubstring) continue;
			
        // Iterate through the SA interval range and write the overlaps
        for(int64_t j = record.ranges.interval[0].lower; j <= record.ranges.interval[0].upper; ++j)
        {
            const SuffixArray* pCurrSAI = (record.flags.isTargetRev()) ? 
												m_pOverlapper->getRevSAI() : m_pOverlapper->getFwdSAI();
            const ReadInfo& queryInfo = m_pOverlapper->getQueryRIT()->getReadInfo(workItem.idx);

            int64_t saIdx = j;

            // The index of the second read is given as the position in the SuffixArray index
            const ReadInfo& targetInfo = m_pOverlapper->getTargetRIT()->getReadInfo(pCurrSAI->get(saIdx).getID());

			// std::cout << queryInfo.id << "\t" << targetInfo.id;

            // Skip self alignments and non-canonical (where the query read has a lexo. higher name)
            if(queryInfo.id != targetInfo.id)
            {    
				// std::cout << "##" << queryInfo.id << "\t" << targetInfo.id << "\n";
                Overlap o = record.toOverlap(queryInfo.id, targetInfo.id, queryInfo.length, targetInfo.length);

                // The alignment logic above produce duplicated overlaps for each edge
                // To avoid this, we skip overlaps where the id of the first coord is lexo. lower than 
                // the second or the match is a containment and the query is reversed (containments can be 
                // output up to 4 times total) or the query is super repeat or o.id[0] > o.id[1]
				if(o.match.isContainment() && record.flags.isQueryRev())
					continue;

				if(o.id[0] < o.id[1])continue;
								
				//assert(isQuerySuperRepeat || o.id[0] > o.id[1]);
				
                ASQG::EdgeRecord edgeRecord( o );
				edgeRecord.write( *m_pWriter );
            }
        }
    }
	
    m_blockList.clear();
    return result;
}

//
//
//
OverlapPostProcess::OverlapPostProcess(std::ostream* pASQGWriter, 
                                       const OverlapAlgorithm* pOverlapper) : m_pASQGWriter(pASQGWriter),
                                                                              m_pOverlapper(pOverlapper)
{

}

//
void OverlapPostProcess::process(const SequenceWorkItem& item, const OverlapResult& result)
{
    m_pOverlapper->writeResultASQG(*m_pASQGWriter, item.read, result);
}
