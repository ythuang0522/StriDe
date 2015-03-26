///-----------------------------------------------
// Copyright 2015 National Chung Cheng University
// Written by Yao-Ting Huang
// Released under the GPL
//-----------------------------------------------
//
// FMIndexWalkProcess.cpp - Implementation of FM-index walk and kmerization of PE reads
//
#include "FMIndexWalkProcess.h"
#include "CorrectionThresholds.h"
#include "HashMap.h"
#include <iomanip>
#include "SAIntervalTree.h"


//#define KMER_TESTING 1


FMIndexWalkProcess::FMIndexWalkProcess(const FMIndexWalkParameters params) : m_params(params)
{
}

FMIndexWalkProcess::~FMIndexWalkProcess()
{

}


FMIndexWalkResult FMIndexWalkProcess::MergeAndKmerize(const SequenceWorkItemPair& workItemPair)
{
	FMIndexWalkResult result;

	//get parameters
	size_t kmerLength = m_params.kmerLength ;
	size_t threshold = (size_t)CorrectionThresholds::Instance().getRequiredSupport(0)-1;

	std::string seqFirst  = workItemPair.first.read.seq.toString() ;
	std::string seqSecond = workItemPair.second.read.seq.toString();

	//Trim head and tail from both ends if there is low-frequency kmer
	seqFirst = trimRead(seqFirst, m_params.kmerLength, threshold, m_params.indices);
	seqSecond = trimRead(seqSecond, m_params.kmerLength, threshold, m_params.indices);
	
	std::string firstKRstr = seqFirst.substr(0, m_params.minOverlap);
	std::string secondKRstr  = seqSecond.substr(0, m_params.minOverlap);
	if( isSuitableForFMWalk(firstKRstr, secondKRstr) )
    {	
		//maxOverlap is limited to 90% of read length which aims to prevent over-greedy search
		size_t maxOverlap = m_params.maxOverlap!=-1?m_params.maxOverlap:
											((workItemPair.first.read.seq.length()+workItemPair.second.read.seq.length())/2)*0.95;

		std::string mergedseq1, mergedseq2;
		//Walk from the 1st end to 2nd end											
        SAIntervalTree SAITree1(&firstKRstr, m_params.minOverlap, maxOverlap, m_params.maxInsertSize, m_params.maxLeaves,
                                            m_params.indices, reverseComplement(secondKRstr));
        SAITree1.mergeTwoReads(mergedseq1);

		//Walk from the 2nd end to 1st end using the other strand
		SAIntervalTree SAITree2(&secondKRstr, m_params.minOverlap, maxOverlap, m_params.maxInsertSize, m_params.maxLeaves,
											m_params.indices, reverseComplement(firstKRstr));

		SAITree2.mergeTwoReads(mergedseq2);
		
		//Only one successful walk from first end, requiring maxUsedLeaves <=1 in order to avoid walking over chimera PE read.
		if(!mergedseq1.empty() && mergedseq2.empty() && SAITree1.getMaxUsedLeaves()<=1 && SAITree2.getMaxUsedLeaves()<=1)
		{
			// std::cout << ">" << SAITree.getKmerCoverage()<< "\n" << mergedseq << "\n" ;
			// std::cout << SAITree1.getMaxUsedLeaves() << "\t" << SAITree1.isBubbleCollapsed() << "\t" << SAITree2.getMaxUsedLeaves() << "\n";
			// getchar();
			result.merge = true ;
			result.correctSequence = mergedseq1 ;
			return result;
		//Only one successful walk from second end
		}else if( mergedseq1.empty() && !mergedseq2.empty() && SAITree2.getMaxUsedLeaves()<=1 && SAITree1.getMaxUsedLeaves() <=1)
		{
			result.merge = true ;
			result.correctSequence = mergedseq2 ;
			return result;
		}
		else if( !mergedseq1.empty() && !mergedseq2.empty() && ( mergedseq1 == reverseComplement(mergedseq2) )  )
		{
			// std::cout << ">" << SAITree.getKmerCoverage()<< "\n" << mergedseq << "\n>" << SAITree2.getKmerCoverage()<< "\n" << mergedseq2 << "\n";
			// std::cout << SAITree.getMaxUsedLeaves() << "\t" << SAITree.isBubbleCollapsed() << "\t" << SAITree2.getMaxUsedLeaves() << "\t" << SAITree2.isBubbleCollapsed()<<"\n";
			// getchar();
			result.merge = true ;
			result.correctSequence = (SAITree1.getKmerCoverage()>SAITree2.getKmerCoverage())? mergedseq1:mergedseq2 ;
			return result;
		}
		// else if( !mergedseq1.empty() && !mergedseq2.empty() && (std::abs( (int)mergedseq1.length() - (int)mergedseq2.length()) <= 3))
		// {
			// 30~40% are chimera and often mixed with repeats leading to many leaves when walking from both ends
			// std::cout << ">" << SAITree.getKmerCoverage()<< "\n" << mergedseq << "\n>" << SAITree2.getKmerCoverage()<< "\n" << mergedseq2 << "\n";
			// std::cout << SAITree.getMaxUsedLeaves() << "\t" << SAITree.isBubbleCollapsed() << "\t" << SAITree2.getMaxUsedLeaves() << "\t" << SAITree2.isBubbleCollapsed()<<"\n";
		// }
    }
	
	
    /** Case 3: kmerize the remaining reads **/
	//Compute kmer freq of each kmer
	KmerContext seqFirstKC(seqFirst, kmerLength, m_params.indices);
	KmerContext seqSecondKC(seqSecond, kmerLength, m_params.indices);

	std::vector<std::string> firstKR ;
	std::vector<std::string> secondKR ;
	int firstMainIdx=-1, secondMainIdx=-1;

	if(seqFirst.length()>=(size_t) kmerLength) 
		firstMainIdx = splitRead( seqFirstKC, firstKR, threshold, m_params.indices);
	if(seqSecond.length()>=(size_t) kmerLength)
		secondMainIdx = splitRead( seqSecondKC, secondKR, threshold, m_params.indices);
    // /***trim and kmerize reads***/

    /*** write kmernized results***/
	if (!firstKR.empty()) result.kmerize =true ;
	if (!secondKR.empty()) result.kmerize2 =true ;
	
	for (int i = 0 ; i<(int)firstKR.size() ; i ++ )
	{
		std::string kmerRead = firstKR.at(i);
		float GCratio =0 ;
		if ( isLowComplexity (kmerRead,GCratio) ) continue;
		if ( maxCon(kmerRead)*3 > kmerRead.length() ) continue;
		if (i==firstMainIdx)  result.correctSequence = kmerRead ;
		else
			result.kmerizedReads.push_back(kmerRead);
	}

	for (int i = 0 ; i<(int)secondKR.size() ; i ++ )
	{
		std::string kmerRead = secondKR.at(i);
		float GCratio =0 ;
		if ( isLowComplexity (kmerRead,GCratio) ) continue;
		if ( maxCon(kmerRead)*3 > kmerRead.length() ) continue;
		if (i==secondMainIdx)  result.correctSequence2 = kmerRead ;
		else
			result.kmerizedReads2.push_back(kmerRead);
	}
	
	return result;
	
}


FMIndexWalkResult FMIndexWalkProcess::MergePairedReads(const SequenceWorkItemPair& workItemPair)
{
	FMIndexWalkResult result;

	//get parameters
	// size_t kmerLength = m_params.kmerLength ;
	size_t threshold = (size_t)CorrectionThresholds::Instance().getRequiredSupport(0)-1;

	std::string seqFirstOriginal  = workItemPair.first.read.seq.toString() ;
	std::string seqSecondOriginal = workItemPair.second.read.seq.toString();

	//Trim head and tail of both ends containing errors
	std::string seqFirst = trimRead(seqFirstOriginal, m_params.kmerLength, threshold,m_params.indices);
	std::string seqSecond = trimRead(seqSecondOriginal, m_params.kmerLength, threshold,m_params.indices);

	if(isSuitableForFMWalk(seqFirst, seqSecond))
    {
		//extract prefix of seqFirst
		std::string firstKRstr = seqFirst.substr(0, m_params.minOverlap);	
		//extract suffix of seqSecond
		std::string secondKRstr = seqSecond.substr(0, m_params.minOverlap);
	
		//Walk from the 1st end to 2nd end
		size_t maxOverlap = m_params.maxOverlap!=-1?m_params.maxOverlap:
											((workItemPair.first.read.seq.length()+workItemPair.second.read.seq.length())/2)*0.9;
											
        SAIntervalTree SAITree(&firstKRstr, m_params.minOverlap, maxOverlap, m_params.maxInsertSize, m_params.maxLeaves,
                                            m_params.indices, reverseComplement(secondKRstr), threshold);
        std::string mergedseq;
        SAITree.mergeTwoReads(mergedseq);

		//Walk from the 2nd end to 1st end 
		SAIntervalTree SAITree2(&secondKRstr, m_params.minOverlap, maxOverlap, m_params.maxInsertSize, m_params.maxLeaves,
											m_params.indices, reverseComplement(firstKRstr), threshold);
		std::string mergedseq2;
		SAITree2.mergeTwoReads(mergedseq2);
			
		//Unipath from 1st end but no path from 2nd end
		if(!mergedseq.empty() && mergedseq2.empty() )
		{
			if(SAITree.getMaxUsedLeaves()>1 || SAITree2.getMaxUsedLeaves()>1){
				// std::cout << ">First"  << "\n" << mergedseq <<  "\n";
			}else if(mergedseq.length()>maxOverlap){
				//uniq path
				result.merge = true ;
				result.correctSequence = mergedseq ;
				return result;
			}
		//Unipath from 2nd end but no path from 1st end
		}
		else if(mergedseq.empty() && !mergedseq2.empty() )
		{
			if(SAITree.getMaxUsedLeaves()>1 || SAITree2.getMaxUsedLeaves()>1){
				// std::cout << ">Second"  << "\n" << mergedseq2 <<  "\n";
			}else if(mergedseq2.length()>maxOverlap){
				result.merge = true ;
				result.correctSequence = mergedseq2 ;
				return result;
			}
		}
		else if( !mergedseq.empty() && !mergedseq2.empty() && (mergedseq.length()==mergedseq2.length())  )	//mergedseq.length()-mergedseq2.length()<=3
		{
			// if(SAITree.getMaxUsedLeaves()>100){
				// std::cout << ">" << mergedseq.length() << "\n" << mergedseq << "\n";
				// getchar();
			// }
			if(mergedseq.length()>maxOverlap){
				result.merge = true ;
				result.correctSequence = (SAITree.getKmerCoverage()>SAITree2.getKmerCoverage())? mergedseq:mergedseq2 ;
				return result;
			}
		}
		else if(!mergedseq.empty() && !mergedseq2.empty() )
		{
			// std::cout << "Complex walks\t" << mergedseq.length() << "\t" <<mergedseq2.length() << "\n" ;
			// std::cout << ">" << SAITree.getMaxUsedLeaves()<< "\n" << mergedseq << "\n";
			// std::cout << ">" << SAITree2.getMaxUsedLeaves()  << "\n" << mergedseq2 <<  "\n";
			// getchar();
		}
		else	//failure FMwalk
		{
			// std::cout << mergedFlag1 << ":" << mergedFlag2 <<"\n";
		}
    }//end of length > min overlap
	// else
		// std::cout << ">First\n" << seqFirstOriginal << "\n>Second\n" << seqSecondOriginal<< "\n";
	
	return result;
}

//
FMIndexWalkResult FMIndexWalkProcess::process(const SequenceWorkItem& workItem)
{
	FMIndexWalkResult result = correct(workItem);
	return result;
}

FMIndexWalkResult FMIndexWalkProcess::correct(const SequenceWorkItem& /*workItem*/)
{
	switch(m_params.algorithm)
	{
	// case FMW_KMERIZE:
		// {
			// return kmerizeLowKmerReadCorrection(workItem);
			// break;
		// }
	default:
		{
			assert(false);
		}
	}
	FMIndexWalkResult result;
	return result;
}


//check necessary conditions for FM-index walk
bool FMIndexWalkProcess::isSuitableForFMWalk(std::string& seqFirst, std::string& seqSecond)
{
	//check minimum read length
	bool isBothReadLengthSufficient = seqFirst.length()>= (size_t) m_params.minOverlap && seqSecond.length() >= (size_t)m_params.minOverlap ;
	if(!isBothReadLengthSufficient) return false;

	//estimate repeat kmer
	// const size_t RepeatKmerFreq = m_params.kd.getRepeatKmerCutoff();
	const size_t RepeatKmerFreq = m_params.kd.getMedian()*1.3; 
	size_t KmerFreq1 = BWTAlgorithms::countSequenceOccurrences( seqFirst, m_params.indices.pBWT );
	bool isFirstReadUnique =  KmerFreq1 < RepeatKmerFreq;

	size_t KmerFreq2 = BWTAlgorithms::countSequenceOccurrences( seqSecond, m_params.indices.pBWT );
	bool isSecondReadUnique =  KmerFreq2 < RepeatKmerFreq;
	
	 if( (isFirstReadUnique && isSecondReadUnique) ) return true;
	 
	 //
	 // if(KmerFreq1 < m_params.kd.getMedian()*1.8 && KmerFreq2 < m_params.kd.getMedian()*1.8) return 2;
	 
	 return false;
}

// return complexity of seq, default: 0.9
bool  FMIndexWalkProcess::isLowComplexity (std::string seq , float & GCratio)
{
	size_t seqLen = seq.length();
	size_t countG =0 ;
	size_t countC =0 ;
	size_t countT =0 ;
	size_t countA =0 ;

	for (size_t i=0; i<seqLen; i++)
	{
		switch(seq[i]){
			case 'A': countA ++ ;break;
			case 'T': countT ++ ;break;
			case 'C': countC ++ ;break;
			case 'G': countG ++ ;break;
			default:  assert(false);
		}
	}

	GCratio = (float)(countG+countC)/seqLen ;

	if (  ((float) countA/seqLen >=0.9 ) || ((float) countT/seqLen >=0.9 )
	   || ((float) countC/seqLen >=0.9 ) || ((float) countG/seqLen >=0.9 ) )
	   return true;

	return false;

}

//compute the maximum length of consecutive letter AATTTTTTTTTTTCCC
size_t FMIndexWalkProcess::maxCon (std::string s)
{
	size_t c = 1 ;
	size_t max =1 ;
	for (size_t i =1 ;i<s.length();i++)
	{
		if (s[i] =='N') continue;

		if (s[i]!=s[i-1])
		{
			if (c>max) max = c ;
			c=1;
		}
		else
		{
			c++;
			if ( i == s.length() -1 )
			{
				if (c>max) max = c ;
			}
			else if ( s[i]!=s[i+1])
			{
				if (c>max) max = c ;
			}
		}
	}

	return max;

}

//search for a strong interval having high-frequent kmers at both strands 
bool FMIndexWalkProcess::isIntervalExistStrongKmer (std::pair<size_t,size_t> interval,std::vector<size_t> & countQualified)
{
	for(size_t i = interval.first ; i<=interval.second ; i++ )
	{
		//find a high-frequent kmer 
		if (countQualified.at(i)== 2 ) return true;
	}
	return false;
}

//Determine reliability between two intervals defined by existence of strong kmer at one strand
bool FMIndexWalkProcess::isPathReliable(std::pair<size_t,size_t> intervalX, std::pair<size_t,size_t> intervalY,std::vector<size_t> & countQualified)
{
	//Two adjacent strong intervals are assumed to be reliable
	if (intervalX.second+1==intervalY.first) return true;

	size_t start = intervalX.second + 1 ;
	size_t end = intervalY.first-1;
	assert (start<=end);

	//if the path exists strong kmer at one strand, it is reliable.
	//Otherwise, 
	for (size_t i =start  ; i <=end ;i++)
		if (countQualified[i]==0) return false;

	return true;
}

//Merge back two splitted intervals if there exists strong kmer link at one strand
bool FMIndexWalkProcess::isIntervalMerge (std::vector< std::pair<size_t,size_t> > & intervals , std::vector<size_t> & countQualified )
{
	std::vector<bool> stongInterval(intervals.size());
	size_t count = 0 ;
	for (size_t i =0 ;i <intervals.size();i++)
	{
		stongInterval.at(i)= isIntervalExistStrongKmer (intervals[i],countQualified);
		if (stongInterval[i]) count ++ ;
	}
	if (count <2 ) return false;
	else	//This read exists two strong kmers in two intervals, may be over splitted
	{
		int s = -1 ;
		int e = -1 ;
		for (int i=0 ;i<(int)intervals.size();i++)
		{
			//Find the first strong interval
			if (stongInterval.at(i)&& s<0)
			{
				s=i;
				continue;
			}

			//Find the second strong interval
			if (stongInterval[i])
			{
				e=i;
				//check if there exist strong kmers at one strand between them
				if (isPathReliable(intervals[s],intervals[e],countQualified))
				{
					//Extend 1st interval end to 2nd end
					intervals[s].second = intervals[e].second;
					//Erase the intervals after 1st interval (s+1) till 2nd interval e
					intervals.erase(intervals.begin()+s+1,intervals.begin()+e+1);
					return true;
				}
				s=e;
			}
		}
	}
	return false;

}


//Kmerize the read into subreads at potential error bases
int FMIndexWalkProcess::splitRead (KmerContext& seq, std::vector<std::string> & kmerReads, size_t threshold, BWTIndexSet & index)
{
	if (seq.empty()) return -1 ;

	std::vector<size_t> countQualified (seq.numKmer,0) ;
	for (size_t i=0 ;i<seq.numKmer;i++)
	{
		if (seq.kmerFreqs_same.at(i)>= threshold) countQualified[i]++;
		if (seq.kmerFreqs_revc.at(i)>= threshold) countQualified[i]++;
	}

	//Split the reads into intervals
	std::vector< std::pair<size_t,size_t> > intervals ;
	size_t start = 0 ;
	size_t end = seq.numKmer-1 ;
	for (size_t p = 1; p< seq.numKmer ; p++)
	{
	    //don't split if both strands have kmer feq >= threshold
		if (countQualified[p-1]==2 && countQualified[p]==2) continue;
				
		//kmerize read at pos p if the path is not simple
		if ( !isSimple(seq.kmers[p-1], seq.kmers[p], index, 1) )
		{
			intervals.push_back(std::make_pair (start,p-1));
			start =p;
		}
	}
	
	intervals.push_back(std::make_pair (start,end));
	
	//merge two strong intervals with one strong kmer link in between
	// while (isIntervalMerge(intervals,countQualified));

	size_t maxIntervalNum = 0;
	int mainSeedIdx = -1 ;
	for (size_t i=0;i<intervals.size();i++)
	{
		if (isIntervalExistStrongKmer(intervals[i],countQualified))
		{
			size_t intervalNum = intervals[i].second - intervals[i].first ;
			if (maxIntervalNum < intervalNum)
			{
				maxIntervalNum = intervalNum;
				mainSeedIdx = i;
			}
		}

		std::string k=seq.readSeq.substr(intervals[i].first, intervals[i].second-intervals[i].first+seq.kmerLength);
		kmerReads.push_back(k);
	}

	return mainSeedIdx;
}


// bool FMIndexWalkProcess::existStrongLink (std::string Lkmer,std::string Rkmer,BWTIndexSet & index,size_t threshold)
// {
	// assert (Lkmer.substr(1,Lkmer.length()-1) == Rkmer.substr(0,Lkmer.length()-1) );

	// size_t LkmerSameFreq = BWTAlgorithms::countSequenceOccurrencesSingleStrand(Lkmer, index);
	// size_t LkmerRevcFreq = BWTAlgorithms::countSequenceOccurrencesSingleStrand(reverseComplement(Lkmer), index);

	// size_t RkmerSameFreq = BWTAlgorithms::countSequenceOccurrencesSingleStrand(Rkmer, index);
	// size_t RkmerRevcFreq = BWTAlgorithms::countSequenceOccurrencesSingleStrand(reverseComplement(Rkmer), index);

	// bool LKmerStrong = (LkmerSameFreq>=threshold && LkmerRevcFreq>=threshold) ;
	// bool RKmerStrong = (RkmerSameFreq>=threshold && RkmerRevcFreq>=threshold) ;

	// assert (  LKmerStrong != RKmerStrong  )  ;


	// if ( LKmerStrong )	return  existNextStrongKmer (Lkmer,NK_END,index,threshold) && (numNextKmer(Rkmer,NK_START,index) == 1);
	// //RKmerStrong
	// else 	return  existNextStrongKmer (Rkmer,NK_START,index,threshold) && (numNextKmer(Lkmer,NK_END,index) == 1 );

// }

//return true or false if there is a strong/better kmer 
bool FMIndexWalkProcess::existNextStrongKmer(std::string kmer , NextKmerDir dir ,BWTIndexSet & index,size_t threshold)
{
	char nBases[4] = {'A','T','C','G'} ;
	int kmerLength = kmer.length() ;

	for (size_t i = 0 ; i < 4 ; i++)
	{
		std::string nextKmer;
		if (dir == NK_START) nextKmer = nBases[i] + kmer.substr(0,kmerLength-1);
		else if  (dir == NK_END) nextKmer = kmer.substr(1,kmerLength-1) + nBases[i];
		else assert(false);

		size_t nextKmerSameFreq = BWTAlgorithms::countSequenceOccurrencesSingleStrand(nextKmer, index);
		size_t nextKmerRevcFreq = BWTAlgorithms::countSequenceOccurrencesSingleStrand(reverseComplement(nextKmer), index);
		if (nextKmerSameFreq>=threshold && nextKmerRevcFreq>=threshold) return true ;

	}
	return false;
}


std::string FMIndexWalkProcess::trimRead( std::string readSeq, size_t kmerLength, size_t /*threshold*/, BWTIndexSet & index)
{
	int head=0, tail=readSeq.length()-kmerLength;
	// if (( kc.kmerFreqs_same[head] < threshold || kc.kmerFreqs_revc [head] < threshold ) && ( numNextKmer(kc.kmers[head],NK_START,index) == 0 ))
	
	//Dead end
	if (numNextKmer( readSeq.substr(head, kmerLength), NK_START, index, 1) == 0  )
	{
		//search for the first pos with >=2 branching kmer path, note that the min kmer threshold should be 1 as error kmer has low freq
		for (head = head+1 ; head <= tail ; head++ )
			if (numNextKmer(readSeq.substr(head,kmerLength), NK_START, index, 1) >= 2) break ;
	}

	// if (( kc.kmerFreqs_same[tail] < threshold || kc.kmerFreqs_revc [tail] < threshold )  && ( numNextKmer(kc.kmers[tail],NK_END,index) == 0 ) )
	if ( numNextKmer( readSeq.substr(tail, kmerLength), NK_END, index, 1) == 0  )
	{
		for (tail = tail-1 ; tail >= head  ; tail-- )
			if (numNextKmer( readSeq.substr(tail, kmerLength), NK_END,index, 1) >= 2) break ;
	}

	//all kemrs are dirty , return empty
	if (head > tail)
		return "";
	else
		return readSeq.substr(head,tail-head+kmerLength);

}



size_t FMIndexWalkProcess::numNextKmer(std::string kmer , NextKmerDir dir ,BWTIndexSet & index, size_t threshold)
{
	size_t num = 0 ;
	char nBases[4] = {'A','T','C','G'} ;
	int kmerLength = kmer.length() ;

	for (size_t i = 0 ; i < 4 ; i++)
	{
		std::string next_mer;
		if (dir == NK_START) next_mer = nBases[i] + kmer.substr(0,kmerLength-1);
		else if  (dir == NK_END) next_mer = kmer.substr(1,kmerLength-1) + nBases[i];

		if ( BWTAlgorithms::countSequenceOccurrences(next_mer,index) >= threshold) num++;
	}
	return num ;
}


bool FMIndexWalkProcess::isSimple (std::string Lkmer, std::string Rkmer, BWTIndexSet & index, size_t threshold)
{
	size_t LKmerPathCount = numNextKmer(Lkmer, NK_END, index, threshold);
	size_t RKmerPathCount = numNextKmer(Rkmer, NK_START, index, threshold);
	if ( LKmerPathCount == 1 &&   RKmerPathCount == 1 )  
		return true;
	else 
		return false;
	
	// if( (LKmerPathCount > 1 && LCount == 0) || (RKmerPathCount>1 && RCount == 0) )
		// return false;
	// else if( LKmerPathCount > 1 && (LCount ==1 && existNextStrongKmer (Lkmer,NK_END, index, threshold)) )
		// return false;
	// else if( RKmerPathCount > 1 && RCount ==1 && existNextStrongKmer (Rkmer,NK_START, index, threshold))
		// return false;
	// else 
		// return true;
}



//
//
//
FMIndexWalkPostProcess::FMIndexWalkPostProcess(std::ostream* pCorrectedWriter,
																				std::ostream* pDiscardWriter,
																				const FMIndexWalkParameters params) :
																													m_pCorrectedWriter(pCorrectedWriter),
																													m_pDiscardWriter(pDiscardWriter),
																													m_params(params),
																													m_kmerizePassed(0),
																													m_mergePassed(0),
																													m_qcFail(0)
																				{
																					//m_ptmpWriter = createWriter("NoPESupport.fa");
																				}

//
FMIndexWalkPostProcess::~FMIndexWalkPostProcess()
{
	std::cout << "Reads are kmerized: " << m_kmerizePassed << "\n";
	std::cout << "Reads are merged : "<< m_mergePassed << "\n";
	std::cout << "Reads failed to kmerize or merge: " << m_qcFail << "\n";
}


//
void FMIndexWalkPostProcess::process(const SequenceWorkItem& item, const FMIndexWalkResult& result)
{
	// Determine if the read should be discarded
	bool readQCPass = true;
	if (result.kmerize)
	{
		m_kmerizePassed += 1;
	}
	else
	{
		readQCPass = false;
		m_qcFail += 1;
	}

	SeqRecord record = item.read;
	record.seq = result.correctSequence;

	if (result.correctSequence.empty() && result.kmerizedReads.empty()) ;

	else if (result.kmerize)
	{
		if (!result.correctSequence.empty())
			record.write(*m_pCorrectedWriter);


		for (size_t i=0 ; i< result.kmerizedReads.size() ; i++)
		{
			record.seq = result.kmerizedReads[i];
			record.writeFasta(*m_pDiscardWriter,i);
		}

	}
	else if  (readQCPass || m_pDiscardWriter == NULL)
	{
		record.write(*m_pCorrectedWriter);
	}
	else
	{
		record.write(*m_pDiscardWriter);
	}

}

// Writting results of FMW_HYBRID and FMW_MERGE
void FMIndexWalkPostProcess::process(const SequenceWorkItemPair& itemPair, const FMIndexWalkResult& result)
{
	if (result.merge)
        m_mergePassed += 1;
	else if (  (m_params.algorithm == FMW_HYBRID)  && (result.kmerize ||  result.kmerize2) )
	{
		if (result.kmerize) m_kmerizePassed += 1;
		else m_qcFail += 1;
		if (result.kmerize2) m_kmerizePassed += 1;
		else m_qcFail += 1;
	}
	else
        m_qcFail += 2;

	SeqRecord firstRecord  = itemPair.first.read;
	SeqRecord secondRecord  = itemPair.second.read;

	if (result.merge)
	{
		SeqItem mergeRecord ;
		mergeRecord.id = firstRecord.id.substr (0, firstRecord.id.find('/') ) ;
		mergeRecord.seq = result.correctSequence;
		mergeRecord.write(*m_pCorrectedWriter);
	}
	else if(m_params.algorithm == FMW_HYBRID )
	{
		if (!result.correctSequence.empty())
		{
			firstRecord.seq = result.correctSequence;
			firstRecord.write(*m_pDiscardWriter);
		}
		for (size_t i=0 ; i< result.kmerizedReads.size() ; i++)
		{
			firstRecord.seq = result.kmerizedReads[i];
			firstRecord.writeFasta(*m_pDiscardWriter,i);
		}

		if (!result.correctSequence2.empty())
		{
			secondRecord.seq = result.correctSequence2;
			secondRecord.write(*m_pDiscardWriter);
		}
		for (size_t i=0 ; i< result.kmerizedReads2.size() ; i++)
		{
			secondRecord.seq = result.kmerizedReads2[i];
			secondRecord.writeFasta(*m_pDiscardWriter,i);
		}
	}
}
