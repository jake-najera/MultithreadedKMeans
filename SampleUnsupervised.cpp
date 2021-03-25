/*-------------------------------------------------------
File Name: SampleUnsupervised.cpp

Purpose: Implementation of ParallelClustering for unsupervised learning

Author: Jake Najera

Creation date: 3/15/2020
---------------------------------------------------------*/
#include <iostream>
#include "SampleUnsupervised.h"
#include <limits.h>

using std::cout;
using std::endl;

// ------------------------------------------------------------------------------------

SampleUnsupervised::SampleUnsupervised()
{
	mpFeatures = new std::vector<double>();
	mpID = new std::string();
}

// ------------------------------------------------------------------------------------

SampleUnsupervised::~SampleUnsupervised()
{
	delete mpFeatures;
	delete mpID;
}

// ------------------------------------------------------------------------------------

SampleUnsupervised::SampleUnsupervised(const SampleUnsupervised& Rhs)
{
	mpFeatures = new std::vector<double>();
	mpID = new std::string();
	*mpFeatures = *Rhs.mpFeatures;
	*mpID = *Rhs.mpID;
}

// ------------------------------------------------------------------------------------

SampleUnsupervised::SampleUnsupervised(const std::vector<double>& Features, const std::string& ID)
{
	mpFeatures = new std::vector<double>();
	mpID = new std::string();
	*mpFeatures = Features;
	*mpID = ID;
}

// ------------------------------------------------------------------------------------

SampleUnsupervised& SampleUnsupervised::operator=(const SampleUnsupervised& Rhs)
{
	if (this != &Rhs)
	{
		*mpFeatures = *Rhs.mpFeatures;
		*mpID = *Rhs.mpID;
	}

	return *this;
}

// ------------------------------------------------------------------------------------

SampleSetUnsupervised::SampleSetUnsupervised()
{
	//mIsSvdDirty = true;
	mpSamples = new std::vector<SampleUnsupervised>();
	//mpSingularValueDecomposition_S = new std::vector<double>();
}

// ------------------------------------------------------------------------------------

SampleSetUnsupervised::~SampleSetUnsupervised()
{
	delete mpSamples;
	//delete mpSingularValueDecomposition_S;
}

// ------------------------------------------------------------------------------------

SampleSetUnsupervised::SampleSetUnsupervised(const SampleSetUnsupervised& Rhs)
{
	mpSamples = new std::vector<SampleUnsupervised>();
	//mpSingularValueDecomposition_S = new std::vector<double>();

	*mpSamples = *Rhs.mpSamples;
	//*mpSingularValueDecomposition_S = *Rhs.mpSingularValueDecomposition_S;
	//mIsSvdDirty = Rhs.mIsSvdDirty;
	//mSingularValueDecomposition_U = Rhs.mSingularValueDecomposition_U;
	//mSingularValueDecomposition_V = Rhs.mSingularValueDecomposition_V;
}

// ------------------------------------------------------------------------------------

SampleSetUnsupervised& SampleSetUnsupervised::operator=(const SampleSetUnsupervised& Rhs)
{
	if (this != &Rhs)
	{
		*mpSamples = *Rhs.mpSamples;
		//*mpSingularValueDecomposition_S = *Rhs.mpSingularValueDecomposition_S;
		//mIsSvdDirty = Rhs.mIsSvdDirty;
		//mSingularValueDecomposition_U = Rhs.mSingularValueDecomposition_U;
		//mSingularValueDecomposition_V = Rhs.mSingularValueDecomposition_V;
	}

	return *this;
}

// ------------------------------------------------------------------------------------

bool SampleSetUnsupervised::AddSample(const SampleUnsupervised& s)
{
	//mIsSvdDirty = true;
	if (true == mpSamples->empty())
	{
		mpSamples->push_back(s);
		return true;
	}

	// All samples have to have the same number of input values
	if (s.mpFeatures->size() == mpSamples->at(0).mpFeatures->size())
	{
		mpSamples->push_back(s);
		return true;
	}

	return false;
}

// ------------------------------------------------------------------------------------

void SampleSetUnsupervised::FeatureScaleAndMeanNormalize(const std::vector<Common::RangeAndAverage>& RangesAndAverages)
{
	//mIsSvdDirty = true;
	size_t inputSize = 0;
	if (mpSamples->empty() == false)
		inputSize = mpSamples->at(0).mpFeatures->size();

	for (size_t inputIdx = 0; inputIdx < inputSize; ++inputIdx)
	{
		double range = RangesAndAverages[inputIdx].mRange;
		double average = RangesAndAverages[inputIdx].mAverage;
		if (0.0 == range) continue;

		for (auto& s : *mpSamples)
			s.mpFeatures->at(inputIdx) = (s.mpFeatures->at(inputIdx) - average) / range;
	}
}

// ------------------------------------------------------------------------------------

void SampleSetUnsupervised::MinimumsAndMaximums(SampleUnsupervised& Min, SampleUnsupervised& Max) const
{
	Min.mpFeatures->clear();		Max.mpFeatures->clear();

	if (mpSamples->empty() == false)
	{
		size_t numberOfFeatures = mpSamples->at(0).mpFeatures->size();
		Min.mpFeatures->resize(numberOfFeatures, std::numeric_limits<double>::max());
		Max.mpFeatures->resize(numberOfFeatures, std::numeric_limits<double>::min());

		for (auto& sample : *mpSamples)
			for (size_t idx = 0; idx < numberOfFeatures; ++idx)
			{
				double value = sample.mpFeatures->at(idx);
				if (value < Min.mpFeatures->at(idx)) Min.mpFeatures->at(idx) = value;
				if (value > Max.mpFeatures->at(idx)) Max.mpFeatures->at(idx) = value;
			}
	}
}

// ------------------------------------------------------------------------------------

std::vector<double> SampleSetUnsupervised::RetainedVariancePerDimensionCount()
{
	return std::vector<double>();
}

// ------------------------------------------------------------------------------------

std::vector<Common::RangeAndAverage> SampleSetUnsupervised::ComputeRangesAndAverages(const std::vector<SampleSetUnsupervised>& TrainingSetsUnsupervised)
{
	std::vector<Common::RangeAndAverage> result;

	// Make sure all the samples of all training sets have the same number of features
	// Get the 1st one that has at least 1 sample in it
	size_t numberOfFeatures = 0;
	for (auto& ts : TrainingSetsUnsupervised)
		if (ts.Samples().size() != 0)
		{
			numberOfFeatures = ts.Samples()[0].mpFeatures->size();
			break;
		}

	// Now loop through all of them to check that they all have the same number of features
	for (auto& ts : TrainingSetsUnsupervised)
		for (auto& s : ts.Samples())
			if (s.mpFeatures->size() != numberOfFeatures)
				return result; // Empty vector

	// At this point, all samples in all training sets have the same number of features
	for (size_t featureIdx = 0; featureIdx < numberOfFeatures; ++featureIdx)
	{
		double total = 0.0, minValue = std::numeric_limits<double>::max(), maxValue = std::numeric_limits<double>::min();
		unsigned int totalNumberOfSamples = 0;
		for (auto& ts : TrainingSetsUnsupervised)
		{
			for (auto& s : ts.Samples())
			{
				double value = s.mpFeatures->at(featureIdx);

				if (value < minValue) minValue = value;
				if (value > maxValue) maxValue = value;
				total += value;
				++totalNumberOfSamples;
			}
		}

		result.push_back(Common::RangeAndAverage(maxValue - minValue, total / (double)(totalNumberOfSamples)));
	}

	return result;
}

// ------------------------------------------------------------------------------------

std::vector<double> SampleSetUnsupervised::Averages(const std::vector<const SampleUnsupervised*>& Samples)
{
	std::vector<double> averages;
	size_t numberOfSamples = Samples.size();
	if (0 == numberOfSamples) return averages;
	size_t numberOfFeatures = Samples[0]->mpFeatures->size();
	if (0 == numberOfFeatures) return averages;

	averages.resize(numberOfFeatures, 0.0);
	for (auto pSample : Samples)
		for (size_t idx = 0; idx < numberOfFeatures; ++idx)
			averages[idx] += pSample->mpFeatures->at(idx);

	for (size_t idx = 0; idx < numberOfFeatures; ++idx)
		averages[idx] /= (double)numberOfSamples;

	return averages;
}

// ------------------------------------------------------------------------------------