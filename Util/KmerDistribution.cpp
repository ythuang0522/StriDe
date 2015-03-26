//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// KmerDistribution - Histogram of kmer frequencies
//
#include "KmerDistribution.h"
#include <assert.h>
#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <limits>
#include <algorithm>
#include <vector>       // std::vector
#include <math.h>

KmerDistribution::KmerDistribution()
{

}

void KmerDistribution::add(int kcount)
{
    m_data[kcount]++;
	m_rawdata.push_back(kcount);
}

double KmerDistribution::getCumulativeProportionLEQ(int n) const
{
    std::vector<int> countVector = toCountVector(1000);
    int64_t sum = 0;
    for(size_t i = 0; i < countVector.size(); ++i)
    {
        sum += countVector[i];
    }

    std::vector<double> cumulativeVector(countVector.size());
    int64_t runningSum = 0;
    for(size_t i = 0; i < cumulativeVector.size(); ++i)
    {
        runningSum += countVector[i];
        cumulativeVector[i] = (double)runningSum / sum;
        //std::cout << "CUMUL " << i << " " << cumulativeVector[i] << "\n";
    }

    assert(n < (int)cumulativeVector.size());
    return cumulativeVector[n];
}

size_t KmerDistribution::getCutoffForProportion(double p) const
{
    int MAX_COUNT = 1000;
    std::vector<int> countVector = toCountVector(MAX_COUNT);
    int64_t sum = 0;
    for(size_t i = 0; i < countVector.size(); ++i)
        sum += countVector[i];

    double c = 0.0f;
    for(size_t i = 0; i < countVector.size(); ++i)
    {
		// std::cout << i << ":" << countVector[i] << "\n";
		// getchar();
        c += (double)countVector[i] / sum;
        if(c > p)
            return i;
    }

    return MAX_COUNT;
}

//
int KmerDistribution::findFirstLocalMinimum() const
{
    std::vector<int> countVector = toCountVector(1000);
    if(countVector.empty())
        return -1;

    std::cout << "CV: " << countVector.size() << "\n";
    int prevCount = countVector[1];
    double threshold = 0.75;
    for(size_t i = 2; i < countVector.size(); ++i)
    {
        int currCount = countVector[i];
        double ratio = (double)currCount / prevCount;
        std::cout << i << " " << currCount << " " << ratio << "\n";
        if(ratio > threshold)
            return i;
        prevCount = currCount;
    }
    return -1;
}

int KmerDistribution::getCensoredMode(size_t n) const
{
    std::vector<int> countVector = toCountVector(1000);
    if(countVector.size() < n)
        return -1;
    int modeIdx = -1;
    int modeCount = -1;

    for(; n < countVector.size(); ++n)
    {
        if(countVector[n] > modeCount)
        {
            modeCount = countVector[n];
            modeIdx = n;
        }
    }
    return modeIdx;
}

// Find the boundary of the kmers that are likely erroneous
// We do this by finding the value between 1 and the trusted mode
// that contributes the fewest 
int KmerDistribution::findErrorBoundary() const
{
    int mode = getCensoredMode(5);
    if(mode == -1)
        return -1;

    std::cerr << "Trusted kmer mode: " << mode  << "\n";
    std::vector<int> countVector = toCountVector(1000);
    if(countVector.empty())
        return -1;

    int runningSum = 0;
    double minContrib = std::numeric_limits<double>::max();
    int idx = -1;
    for(int i = 1; i < mode; ++i)
    {
        runningSum += countVector[i];
        double v = (double)countVector[i] / runningSum;
        if(v < minContrib)
        {
            minContrib = v;
            idx = i;
        }
    }
    return idx;
}

// Find the boundary of the kmers that are likely erroneous
// We do this by finding the value between 1 and the trusted mode
// that contributes the fewest 
int KmerDistribution::findErrorBoundaryByRatio(double ratio) const
{
    int mode = getCensoredMode(5);
    if(mode == -1)
        return -1;

    std::cerr << "Trusted kmer mode: " << mode  << "\n";
    std::vector<int> countVector = toCountVector(1000);
    if(countVector.empty())
        return -1;

    for(int i = 1; i < mode - 1; ++i)
    {
        int currCount = countVector[i];
        int nextCount  = countVector[i+1];
        double cr = (double)currCount / nextCount;
        if(cr < ratio)
            return i;
    }
    return -1;
}

// 
std::vector<int> KmerDistribution::toCountVector(int max) const
{
    std::vector<int> out;
    if(m_data.empty())
        return out;

    int min = 0;

    for(int i = min; i <= max; ++i)
    {
        std::map<int,int>::const_iterator iter = m_data.find(i);
        int v = (iter != m_data.end()) ? iter->second : 0;
        out.push_back(v);
    }
    return out;
}

size_t KmerDistribution::getTotalKmers() const
{
    size_t sum = 0;
    std::map<int,int>::const_iterator iter = m_data.begin();
    for(; iter != m_data.end(); ++iter)
        sum += iter->second;
    return sum;
}

size_t KmerDistribution::getNumberWithCount(size_t c) const
{
    std::map<int,int>::const_iterator iter = m_data.find(c);
    if(iter != m_data.end())
        return iter->second;
    else
        return 0;
}

// for compatibility with old code
void KmerDistribution::print(int max) const
{
    print(stdout, max);
}

void KmerDistribution::print(FILE* fp, int max) const
{
    fprintf(fp, "Kmer coverage histogram\n");
    fprintf(fp, "cov\tcount\n");

    int maxCount = 0;
    std::map<int,int>::const_iterator iter = m_data.begin();
    for(; iter != m_data.end(); ++iter)
    {
        if(iter->first <= max)
            fprintf(fp, "%d\t%d\n", iter->first, iter->second);
        else
            maxCount += iter->second;
    }
    fprintf(fp, ">%d\t%d\n", max, maxCount);

}

//compute median and std
void KmerDistribution::computeKDAttributes()
{
	std::sort ( m_rawdata.begin(), m_rawdata.end()) ;
    
	//compute median
	m_median = m_rawdata.at(m_rawdata.size()/2);
	
	//compute standard deviation
	int dev = 0;
	double sdev=0;
	for(size_t i=0; i < m_rawdata.size(); ++i)
	{
		dev = (m_rawdata.at(i) - m_median)*(m_rawdata.at(i) - m_median);
        sdev = sdev + dev;
	}
    double var = sdev / (m_rawdata.size() - 1);
    m_std = sqrt(var);

	// double freq95 = getCutoffForProportion(0.95);
	// m_repeatKmerCutoff = m_std > m_median*2? m_median*1.5: m_median*1.3;
	m_repeatKmerCutoff = m_median*1.3;
	// m_repeatKmerCutoff = getCutoffForProportion(0.8);
	// m_repeatKmerCutoff = (double) m_median*(0.39+0.53* (freq95/(double)m_median));
	
}
