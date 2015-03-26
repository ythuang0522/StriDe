///-----------------------------------------------
///-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ErrorCorrectProcess - Wrapper to perform error correction
// for a sequence work item
//
#include "ErrorCorrectProcess.h"
#include "CorrectionThresholds.h"
#include "HashMap.h"
#include "multiple_alignment.h"
#include "KmerOverlaps.h"
#include <iomanip>
#include "FMIndexWalkProcess.h"

//#define KMER_TESTING 1



//
//
//
ErrorCorrectProcess::ErrorCorrectProcess(const ErrorCorrectParameters params) : m_params(params)
{
	m_params.depthFilter = 10000;
}

//
ErrorCorrectProcess::~ErrorCorrectProcess()
{

}

ErrorCorrectResult ErrorCorrectProcess::process(const SequenceWorkItem& workItem)
{
        ErrorCorrectResult result = correct(workItem);
        if(!result.kmerQC && !result.overlapQC && m_params.printOverlaps)
        std::cout << workItem.read.id << " failed error correction QC\n";
        return result;
}

ErrorCorrectResult ErrorCorrectProcess::correct(const SequenceWorkItem& workItem)
{
	switch(m_params.algorithm)
	{
	case ECA_HYBRID:
		{
			ErrorCorrectResult result = kmerCorrection(workItem);
			if(!result.kmerQC)
			return overlapCorrectionNew(workItem);
			else
			return result;
			break;
		}
	case ECA_KMER:
		{
			return kmerCorrection(workItem);
			break;
		}
	case ECA_OVERLAP:
		{
			return overlapCorrectionNew(workItem);
			break;
		}
	case ECA_THREAD:
		{
			//return threadingCorrection(workItem);
			break;
		}
	default:
		{
			assert(false);
		}
	}
	ErrorCorrectResult result;
	return result;
}

//
ErrorCorrectResult ErrorCorrectProcess::overlapCorrectionNew(const SequenceWorkItem& workItem)
{
	assert(m_params.indices.pBWT != NULL);
	assert(m_params.indices.pSSA != NULL);

	ErrorCorrectResult result;
	std::string current_sequence = workItem.read.seq.toString();
	std::string consensus;

	/**************************************************************************************/
	int parameterThreshold = CorrectionThresholds::Instance().getRequiredSupport(0)-1;
	if  ( parameterThreshold < 0 ) parameterThreshold=0;
		size_t threshold = (size_t)parameterThreshold ;
	/**************************************************************************************/

	int num_rounds = m_params.numOverlapRounds;
	bool isFirstRound=true;
	for(int round = 0; round < num_rounds; ++round)
	{
        KmerContext  kc (current_sequence,m_params.kmerLength,m_params.indices);
		bool allGoodKmer = true;
        int ErrorIdx=-1;
		//Locate the error index in the read via (1) kmer freq difference between two adjacent kmers and (2) kmer frequency centered at the error base
		/***   kmer freq diff between two adjacent kmers
				49:32
				50:33
				51:14   --> The end base of 51th kmer has touched the error base but kmer freq is still > kmer threshold
				52:14
				53: 9 
				============================================================================
				46:854:48
				47:836:62
				...
				60:45:52
				61:1:0		->The end base of 61th kmer touched the error
        ***/
		// std::cout << current_sequence << "\n";
        for (size_t i = 0 ; i < kc.numKmer; i++){
            // std::cout << i <<":" << kc.kmerFreqs_same[i] << ":" << kc.kmerFreqs_revc [i]  << ":" <<kc.numKmer << "\n";
            // if ( kc.kmerFreqs_same.at(i)< threshold || kc.kmerFreqs_revc.at(i) < threshold) allGoodKmer=false;
			if ( kc.kmerFreqs_same.at(i) + kc.kmerFreqs_revc.at(i) < threshold*2) allGoodKmer=false;
            
            if(i<kc.numKmer-1)
			{
				//(1) Compute kmer freq decrement between two adjacent kmers, note that kmerFreq is unsigned and required casting to int
				bool isFwdFreqLargeDiff= kc.kmerFreqs_same.at(i)>threshold?((int)kc.kmerFreqs_same.at(i)-(int)kc.kmerFreqs_same.at(i+1))/(double)kc.kmerFreqs_same.at(i)>=0.5 : false;
				bool isRvcFreqLargeDiff= kc.kmerFreqs_revc.at(i)>threshold?((int)kc.kmerFreqs_revc[i]-(int)kc.kmerFreqs_revc.at(i+1))/(double)kc.kmerFreqs_revc[i]>=0.5 : false;
				isFwdFreqLargeDiff= (int)kc.kmerFreqs_same.at(i)-(int)kc.kmerFreqs_same.at(i+1) >10 && isFwdFreqLargeDiff;
				isRvcFreqLargeDiff= (int)kc.kmerFreqs_revc.at(i)-(int)kc.kmerFreqs_revc.at(i+1) >10 && isRvcFreqLargeDiff;
				
				if ( isFwdFreqLargeDiff && isRvcFreqLargeDiff ) 
				{
					int tmpErrorIdx=i+m_params.kmerLength;
					int kmer_idx = tmpErrorIdx - m_params.kmerLength/2;
					if(kmer_idx>=(int)kc.numKmer) kmer_idx=kc.numKmer-1;
					// size_t kmer_idx_freq=kc.kmerFreqs_same.at(kmer_idx) + kc.kmerFreqs_revc.at(kmer_idx);
					
					//Flanking kmer freq of correct bases
					// size_t avgCount = (kc.kmerFreqs_same.at(i)+kc.kmerFreqs_revc.at(i))*3/2;

					//case 1: correction of sequencing error
					//if( /*isFirstRound &&*/ (kc.kmerFreqs_same.at(kmer_idx)< threshold || kc.kmerFreqs_revc.at(kmer_idx) < threshold))
					if( /*isFirstRound &&*/ (kc.kmerFreqs_same.at(kmer_idx) + kc.kmerFreqs_revc.at(kmer_idx) < threshold*2))
					{
						allGoodKmer=false;
						if(attemptKmerCorrection(tmpErrorIdx, kmer_idx, threshold, current_sequence))
							break;
						else if(!isFirstRound)
						{
							ErrorIdx=((int)i-4>=0)?i-4:0;
							break;
						}
					}
					//case 2: Convert heterozygous SNPs of diploid genomes into homozygotes
					// else if(isFirstRound && m_params.isDiploid && kmer_idx_freq < avgCount && avgCount<200)
					// {
						// if(tmpErrorIdx+1> (int)kc.numKmer-1) continue;
						// if(kc.kmerFreqs_same.at(i+1)<threshold || kc.kmerFreqs_revc.at(i+1)<threshold ) continue;
						// isFwdFreqLargeDiff= kc.kmerFreqs_same.at(tmpErrorIdx+1)>threshold?((int)kc.kmerFreqs_same.at(tmpErrorIdx+1)-(int)kc.kmerFreqs_same[tmpErrorIdx])/(double)kc.kmerFreqs_same.at(tmpErrorIdx+1)>=0.5 : false;
						// isRvcFreqLargeDiff= kc.kmerFreqs_revc.at(tmpErrorIdx+1)>threshold?((int)kc.kmerFreqs_revc.at(tmpErrorIdx+1)-(int)kc.kmerFreqs_revc[tmpErrorIdx])/(double)kc.kmerFreqs_revc.at(tmpErrorIdx+1)>=0.5 : false;
						// isFwdFreqLargeDiff= (int)kc.kmerFreqs_same.at(tmpErrorIdx+1)-(int)kc.kmerFreqs_same[tmpErrorIdx] >10 && (tmpErrorIdx+1)< (int)kc.numKmer;
						// isRvcFreqLargeDiff= (int)kc.kmerFreqs_revc.at(tmpErrorIdx+1)-(int)kc.kmerFreqs_revc[tmpErrorIdx] >10  && (tmpErrorIdx+1)<(int)kc.numKmer;
						// avgCount = (kc.kmerFreqs_same.at(tmpErrorIdx+1)+kc.kmerFreqs_revc.at(tmpErrorIdx+1))*3/2;
						// if(isFwdFreqLargeDiff && isFwdFreqLargeDiff && avgCount<200)
						// {
							// // std::cout << ">" << tmpErrorIdx  << "\n" << current_sequence << "\n";
							// // getchar();
							// if(attemptHeteroCorrection(tmpErrorIdx, kmer_idx, threshold, avgCount, current_sequence))
							// {
								// // i=tmpErrorIdx;	//jump to the right of the corrected base
								// // continue;
							// }
						// }
					// }
					//case 3: Leave the correction to overlap alignment, compute the leftmost error idx
					//else if(!isFirstRound && (kc.kmerFreqs_same.at(kmer_idx)< threshold || kc.kmerFreqs_revc.at(kmer_idx) < threshold) && ErrorIdx==-1 )
					// else if(!isFirstRound && (kc.kmerFreqs_same.at(kmer_idx) + kc.kmerFreqs_revc.at(kmer_idx) < threshold*2) && ErrorIdx==-1 )
					// {
						// allGoodKmer=false;
						// ErrorIdx=((int)i-4>=0)?i-4:0;
						// break;
					// }
				}

				//(1) Compute kmer freq increment between two adjacent kmers, note that kmerFreq is unsigned and required casting to int
				isFwdFreqLargeDiff= kc.kmerFreqs_same.at(i+1)>threshold?((int)kc.kmerFreqs_same.at(i+1)-(int)kc.kmerFreqs_same.at(i))/(double)kc.kmerFreqs_same.at(i+1)>=0.5 : false;
				isRvcFreqLargeDiff= kc.kmerFreqs_revc.at(i+1)>threshold?((int)kc.kmerFreqs_revc.at(i+1)-(int)kc.kmerFreqs_revc.at(i))/(double)kc.kmerFreqs_revc.at(i+1)>=0.5 : false;
				isFwdFreqLargeDiff= (int)kc.kmerFreqs_same.at(i+1)-(int)kc.kmerFreqs_same.at(i) >10 && isFwdFreqLargeDiff;
				isRvcFreqLargeDiff= (int)kc.kmerFreqs_revc.at(i+1)-(int)kc.kmerFreqs_revc.at(i) >10 && isRvcFreqLargeDiff;
                
				if ( isFwdFreqLargeDiff  && isRvcFreqLargeDiff) 
				{
					int tmpErrorIdx=i;
					int kmer_idx = (tmpErrorIdx >= m_params.kmerLength/2 ? tmpErrorIdx - m_params.kmerLength/2 : 0);
					// size_t kmer_idx_freq=kc.kmerFreqs_same.at(kmer_idx)+kc.kmerFreqs_revc.at(kmer_idx);
					// size_t avgCount = (kc.kmerFreqs_same.at(i+1)+kc.kmerFreqs_revc.at(i+1))*3/2;

					//case 1: Sequencing error
					// if ( /*isFirstRound && */ (kc.kmerFreqs_same.at(kmer_idx)< threshold || kc.kmerFreqs_revc.at(kmer_idx) < threshold) ){
					if ( /*isFirstRound && */ (kc.kmerFreqs_same.at(kmer_idx) + kc.kmerFreqs_revc.at(kmer_idx) < threshold*2) ){
						// ErrorIdx=i+1;
						allGoodKmer=false;

						if(attemptKmerCorrection(tmpErrorIdx, kmer_idx, threshold, current_sequence))
						{
							break;
						}else if(!isFirstRound)
						{
							ErrorIdx=i+1;
							break;
						}
					}
					//case 3: Leave the correction to overlap alignment, compute the leftmost error idx
					// else if(!isFirstRound && (kc.kmerFreqs_same.at(kmer_idx)< threshold || kc.kmerFreqs_revc.at(kmer_idx) < threshold) && ErrorIdx==-1 )
					// else if(!isFirstRound && (kc.kmerFreqs_same.at(kmer_idx)+ kc.kmerFreqs_revc.at(kmer_idx) < threshold*2) && ErrorIdx==-1 )
					// {
						// allGoodKmer=false;
						// ErrorIdx=i+1;
						// break;
					// }
				}
			}
			
			//if(ErrorIdx==-1 && (kc.kmerFreqs_same.at(i)==0 || kc.kmerFreqs_revc.at(i)==0)) ErrorIdx=((int)i-4>=0)?i-4:0;;
        }// end of for each kmer

        //no need for correction if all kmers are good or bad
        if (allGoodKmer)
        {
            result.correctSequence = current_sequence;
            result.overlapQC = true;
            return result;
        }
		
		//Reset this loop if it's first round and the last kmer has been successfully corrected 
		if(isFirstRound)
		{
			isFirstRound=false;
			round--;	//redo this loop round
			continue;
		}

		if(ErrorIdx==-1)ErrorIdx=0;	//If no kmer freq diff > 10, assume 0

		// //Finally, try slower multiple alignment correction
		MultipleAlignment multiple_alignment = KmerOverlaps::buildMultipleAlignment(current_sequence,
																				m_params.kmerLength,
																				current_sequence.length()/2, 	//m_params.minOverlap
																				m_params.minIdentity - (double) (round)*0.01,
																				threshold,
																				m_params.indices,
																				ErrorIdx,	//targetidx holds error idx
																				kc); 

		bool last_round = (round == num_rounds - 1);
		if(last_round)
			consensus = multiple_alignment.calculateBaseConsensus(kc, threshold);
		else
			current_sequence = multiple_alignment.calculateBaseConsensus(kc, threshold);
			
		// if(last_round){
			// multiple_alignment.print(200);
			// std::cout << ">" <<round <<":" << m_params.minIdentity <<"\n" << consensus << "\n";
			// getchar();
		// }
	}

	if(!consensus.empty())
	{
		result.correctSequence = consensus;
		result.overlapQC = true;
	}
	else
	{
		// Return the unmodified query sequence
		result.correctSequence = current_sequence;
		result.overlapQC = true;
	}

	return result;
}


// Correct a read with a k-mer based corrector
ErrorCorrectResult ErrorCorrectProcess::kmerCorrection(const SequenceWorkItem& workItem)
{
	assert(m_params.indices.pBWT != NULL);
	assert(m_params.indices.pCache != NULL);

	ErrorCorrectResult result;

	typedef std::map<std::string, int> KmerCountMap;
	KmerCountMap kmerCache;

	SeqRecord currRead = workItem.read;
	std::string readSequence = workItem.read.seq.toString();

#ifdef KMER_TESTING
	std::cout << "Kmer correcting read " << workItem.read.id << "\n";
#endif

	if((int)readSequence.size() < m_params.kmerLength)
	{
		// The read is shorter than the kmer length, nothing can be done
		result.correctSequence = readSequence;
		result.kmerQC = false;
		return result;
	}

	int n = readSequence.size();
	int nk = n - m_params.kmerLength + 1;

	// Are all kmers in the read well-represented?
	bool allSolid = false;
	bool done = false;
	int rounds = 0;
	int maxAttempts = m_params.numKmerRounds;

	// For each kmer, calculate the minimum phred score seen in the bases
	// of the kmer
	std::vector<int> minPhredVector(nk, 0);
	for(int i = 0; i < nk; ++i)
	{
		int end = i + m_params.kmerLength - 1;
		int minPhred = std::numeric_limits<int>::max();
		for(int j = i; j <= end; ++j)
		{
			int ps = workItem.read.getPhredScore(j);
			if(ps < minPhred)
			minPhred = ps;
		}
		minPhredVector[i] = minPhred;
	}

	while(!done && nk > 0)
	{
		// Compute the kmer counts across the read
		// and determine the positions in the read that are not covered by any solid kmers
		// These are the candidate incorrect bases
		std::vector<int> countVector(nk, 0);
		std::vector<int> solidVector(n, 0);

		for(int i = 0; i < nk; ++i)
		{
			std::string kmer = readSequence.substr(i, m_params.kmerLength);

			// First check if this kmer is in the cache
			// If its not, find its count from the fm-index and cache it
			int count = 0;
			KmerCountMap::iterator iter = kmerCache.find(kmer);

			if(iter != kmerCache.end())
			{
				count = iter->second;
			}
			else
			{
				count = BWTAlgorithms::countSequenceOccurrences(kmer, m_params.indices);
				kmerCache.insert(std::make_pair(kmer, count));
			}

			// Get the phred score for the last base of the kmer
			int phred = minPhredVector[i];
			countVector[i] = count;
			//            std::cout << i << "\t" << phred << "\t" << count << "\n";

			// Determine whether the base is solid or not based on phred scores
			int threshold = CorrectionThresholds::Instance().getRequiredSupport(phred);
			if(count >= threshold)
			{
				for(int j = i; j < i + m_params.kmerLength; ++j)
				solidVector[j] = 1;
			}
		}

		allSolid = true;
		for(int i = 0; i < n; ++i)
		{
#ifdef KMER_TESTING
			std::cout << "Position[" << i << "] = " << solidVector[i] << "\n";
#endif
			if(solidVector[i] != 1)
			allSolid = false;
		}

#ifdef KMER_TESTING
		std::cout << "Read " << workItem.read.id << (allSolid ? " is solid\n" : " has potential errors\n");
#endif

		// Stop if all kmers are well represented or we have exceeded the number of correction rounds
		if(allSolid || rounds++ > maxAttempts)
		break;

		// Attempt to correct the leftmost potentially incorrect base
		bool corrected = false;
		for(int i = 0; i < n; ++i)
		{
			if(solidVector[i] != 1)
			{
				// Attempt to correct the base using the leftmost covering kmer
				int phred = workItem.read.getPhredScore(i);
				int threshold = CorrectionThresholds::Instance().getRequiredSupport(phred);

				int left_k_idx = (i + 1 >= m_params.kmerLength ? i + 1 - m_params.kmerLength : 0);
				corrected = attemptKmerCorrection(i, left_k_idx, std::max(countVector[left_k_idx], threshold), readSequence);
				if(corrected)
				break;

				// base was not corrected, try using the rightmost covering kmer
				size_t right_k_idx = std::min(i, n - m_params.kmerLength);
				corrected = attemptKmerCorrection(i, right_k_idx, std::max(countVector[right_k_idx], threshold), readSequence);
				if(corrected)
				break;
			}
		}

		// If no base in the read was corrected, stop the correction process
		if(!corrected)
		{
			assert(!allSolid);
			done = true;
		}
	}

	if(allSolid)
	{
		result.correctSequence = readSequence;
		result.kmerQC = true;
	}
	else
	{
		result.correctSequence = workItem.read.seq.toString();
		result.kmerQC = false;
	}
	return result;
}


// Attempt to correct the base at position idx in readSequence. Returns true if a correction was made
// The correction is made only if the count of the corrected kmer is at least minCount
// And there are other alleles with kmer freq <= avgCount and >= minCount
bool ErrorCorrectProcess::attemptHeteroCorrection(size_t i, size_t k_idx, size_t minCount, size_t avgCount, std::string& readSequence)
{
	assert(i >= k_idx && i < k_idx + m_params.kmerLength);
	size_t base_idx = i - k_idx;
	char originalBase = readSequence[i];

	std::string kmer = readSequence.substr(k_idx, m_params.kmerLength);
	int bestCount = -1;
	char bestBase = '$';

	bool isAnotherAlleleExisted=false;
	
	for(int j = 0; j < DNA_ALPHABET::size; ++j)
	{
		char currBase = ALPHABET[j];
		kmer[base_idx] = currBase;
		size_t count = BWTAlgorithms::countSequenceOccurrences(kmer, m_params.indices);
		//size_t count = BWTAlgorithms::countSequenceOccurrencesSingleStrand(kmer, m_params.indices);

		//Another allele must have kmer freq < avgCount and > minCount
		// std::cout << currBase << ":" << count << ":" << avgCount << "\n";
		if(count <= avgCount && count >= minCount*2)
		{
			if( currBase != originalBase)
				isAnotherAlleleExisted=true;
				
			if( (int)count > bestCount)
			{
				bestCount = count;
				bestBase = currBase;
			}
		}
	}

	if(isAnotherAlleleExisted)
	{
		readSequence[i] = bestBase;
		return true;
	}
	return false;
}

// Attempt to correct the base at position idx in readSequence. Returns true if a correction was made
// The correction is made only if the count of the corrected kmer is at least minCount
bool ErrorCorrectProcess::attemptKmerCorrection(size_t i, size_t k_idx, size_t minCount, std::string& readSequence)
{
	size_t kmerLength=m_params.kmerLength;
	assert(i >= k_idx && i < k_idx + kmerLength);
	size_t base_idx = i - k_idx;
	char originalBase = readSequence[i];
	std::string kmer = readSequence.substr(k_idx, kmerLength);
	size_t bestCount = 0;
	char bestBase = '$';

#if KMER_TESTING
	std::cout << "i: " << i << " k-idx: " << k_idx << " " << kmer << " " << reverseComplement(kmer) << "\n";
#endif

	for(int j = 0; j < DNA_ALPHABET::size; ++j)
	{
		char currBase = ALPHABET[j];
		// if(currBase == originalBase)
			// continue;
		kmer[base_idx] = currBase;
		size_t count = BWTAlgorithms::countSequenceOccurrences(kmer, m_params.indices);
		//size_t count = BWTAlgorithms::countSequenceOccurrencesSingleStrand(kmer, m_params.indices);

// #if KMER_TESTING
		// printf("%c %c %zu\n", originalBase, currBase, count);
// #endif

		if(count >= minCount*2)
		{
			// printf("%c %c %zu\n", originalBase, currBase, count);
			// Multiple corrections exist, double kmer size
			// if(bestBase != '$' && kmerLength+m_params.kmerLength < readSequence.length())
			// {
				// kmerLength=kmerLength+m_params.kmerLength;
				// k_idx = (i >= kmerLength/2)?i - kmerLength/2:0;
				// if(k_idx>(int)readSequence.length()-kmerLength) k_idx=readSequence.length()-kmerLength;
				// base_idx = i - k_idx;
				// kmer = readSequence.substr(k_idx, kmerLength);	
				// bestCount = 0;
				// bestBase = '$';
				// j=-1;
				// continue;
			// }
			
			bestCount = count;
			bestBase = currBase;
		}
	}

	if(bestCount >= minCount*2 && bestBase!=originalBase)
	{
		assert(bestBase != '$');
		readSequence[i] = bestBase;
		return true;
	}
	return false;
}


//
//
//
ErrorCorrectPostProcess::ErrorCorrectPostProcess(std::ostream* pCorrectedWriter,
std::ostream* pDiscardWriter,
bool bCollectMetrics) :
m_pCorrectedWriter(pCorrectedWriter),
m_pDiscardWriter(pDiscardWriter),
m_bCollectMetrics(bCollectMetrics),
m_totalBases(0), m_totalErrors(0),
m_readsKept(0), m_readsDiscarded(0),
m_kmerQCPassed(0), m_overlapQCPassed(0),
m_qcFail(0)
{

}

//
ErrorCorrectPostProcess::~ErrorCorrectPostProcess()
{
	std::cout << "Reads passed kmer QC check: " << m_kmerQCPassed << "\n";
	std::cout << "Reads passed overlap QC check: " << m_overlapQCPassed << "\n";
	std::cout << "Reads failed QC: " << m_qcFail << "\n";
}

//
void ErrorCorrectPostProcess::writeMetrics(std::ostream* pWriter)
{
	m_positionMetrics.write(pWriter, "Bases corrected by position\n", "pos");
	m_originalBaseMetrics.write(pWriter, "\nOriginal base that was corrected\n", "base");
	m_precedingSeqMetrics.write(pWriter, "\nkmer preceding the corrected base\n", "kmer");
	m_qualityMetrics.write(pWriter, "\nBases corrected by quality value\n\n", "quality");

	std::cout << "ErrorCorrect -- Corrected " << m_totalErrors << " out of " << m_totalBases <<
	" bases (" << (double)m_totalErrors / m_totalBases << ")\n";
	std::cout << "Kept " << m_readsKept << " reads. Discarded " << m_readsDiscarded <<
	" reads (" << (double)m_readsDiscarded / (m_readsKept + m_readsDiscarded)<< ")\n";
}


//
void ErrorCorrectPostProcess::process(const SequenceWorkItem& item, const ErrorCorrectResult& result)
{

	// Determine if the read should be discarded
	bool readQCPass = true;
	if(result.kmerQC)
	{
		m_kmerQCPassed += 1;
	}
	else if(result.overlapQC)
	{
		m_overlapQCPassed += 1;
	}
	else
	{
		readQCPass = false;
		m_qcFail += 1;
	}


	// Collect metrics for the reads that were actually corrected
	if(m_bCollectMetrics && readQCPass)
	{
		collectMetrics(item.read.seq.toString(),
		result.correctSequence.toString(),
		item.read.qual);
	}

	SeqRecord record = item.read;
	record.seq = result.correctSequence;


	if (result.correctSequence.empty()) ;

	else if  (readQCPass || m_pDiscardWriter == NULL)
	{
		record.write(*m_pCorrectedWriter);
		++m_readsKept;
	}
	else
	{
		record.write(*m_pDiscardWriter);
		++m_readsDiscarded;
	}

}


void ErrorCorrectPostProcess::collectMetrics(const std::string& originalSeq,
const std::string& correctedSeq,
const std::string& qualityStr)
{
	size_t precedingLen = 2;
	for(size_t i = 0; i < originalSeq.length(); ++i)
	{
		char qc = !qualityStr.empty() ? qualityStr[i] : '\0';
		char ob = originalSeq[i];

		++m_totalBases;

		m_positionMetrics.incrementSample(i);

		if(!qualityStr.empty())
		m_qualityMetrics.incrementSample(qc);

		m_originalBaseMetrics.incrementSample(ob);

		std::string precedingMer;
		if(i > precedingLen)
		{
			precedingMer = originalSeq.substr(i - precedingLen, precedingLen);
			m_precedingSeqMetrics.incrementSample(precedingMer);
		}

		if(originalSeq[i] != correctedSeq[i])
		{
			m_positionMetrics.incrementError(i);
			if(!qualityStr.empty())
			m_qualityMetrics.incrementError(qc);
			m_originalBaseMetrics.incrementError(ob);

			if(!precedingMer.empty())
			{
				m_precedingSeqMetrics.incrementError(precedingMer);
			}
			++m_totalErrors;
		}
	}
}
