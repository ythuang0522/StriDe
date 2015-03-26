///-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// OverlapProcess - Wrapper for the overlap computation
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

	//Convert list of overlap blocks into Edges
    OverlapBlockList::iterator it = m_blockList.begin();
	for(; it != m_blockList.end(); it++)
    {
        // Read one block
        OverlapBlock record = *it;

        // Iterate through the SA interval range and write the overlaps
        for(int64_t j = record.ranges.interval[0].lower; j <= record.ranges.interval[0].upper; ++j)
        {
            const SuffixArray* pCurrSAI = (record.flags.isTargetRev()) ? 
												m_pOverlapper->getRevSAI() : m_pOverlapper->getFwdSAI();
            const ReadInfo& queryInfo = m_pOverlapper->getQueryRIT()->getReadInfo(workItem.idx);

            int64_t saIdx = j;

            // The index of the second read is given as the position in the SuffixArray index
            const ReadInfo& targetInfo = m_pOverlapper->getTargetRIT()->getReadInfo(pCurrSAI->get(saIdx).getID());

            // Skip self alignments and non-canonical (where the query read has a lexo. higher name)
            if(queryInfo.id != targetInfo.id)
            {    
                Overlap o = record.toOverlap(queryInfo.id, targetInfo.id, queryInfo.length, targetInfo.length);

				//Don't push edges of target of super repeat vertices and small overlap
				// bool isQuerySuperRepeat = m_pOverlapper->isSuperRepeatVertex(o.id[0]);
				// bool isTargetSuperRepeat = m_pOverlapper->isSuperRepeatVertex(o.id[1]);

                // The alignment logic above produce duplicated overlaps for each edge
                // To avoid this, we skip overlaps where the id of the first coord is lexo. lower than 
                // the second or the match is a containment and the query is reversed (containments can be 
                // output up to 4 times total) or the query is super repeat or o.id[0] > o.id[1]
				if(o.match.isContainment() && record.flags.isQueryRev())
					continue;

				if(o.id[0] < o.id[1])continue;

				// //always create edges from super repeat vertex if target is super repeat
				// if(!isQuerySuperRepeat && isTargetSuperRepeat)
					// continue;	
				
				// if(isQuerySuperRepeat && isTargetSuperRepeat && o.id[0] < o.id[1])
					// continue;

				// However, the SuperRepeat of query may still produce duplicated edges while target is super repeat but not processed yet
				// it's unknown why this lead to better results although subsequent duplicate visitor remove them 
                // if(!isQuerySuperRepeat && !isTargetSuperRepeat && o.id[0] < o.id[1])
					// continue;
								
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
