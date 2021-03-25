#pragma once

#include <vector>
#include "Common.h" //needed for feature scaling, range and abverage fxn


class SampleUnsupervised
{
public:
	/**
	*\SampleUnsupervised constructor
	*/
	SampleUnsupervised();

	/**
	*\SampleUnsupervised destructor
	*/
	~SampleUnsupervised();

	/**
	*\SampleUnsupervised copy constructor
	*\param: Rhs is a const reference to SampleUnsupervised
	*/
	SampleUnsupervised(const SampleUnsupervised& Rhs);

	/**
	*\SampleUnsupervised constructor, using features and an optional ID
	*\param: Features is a reference to constant vector of doubles, which will be used to initialize this newly created sample.
	*\param: ID is an optional argument which can be used to set a specific ID for this sample. Useful to keep track of specific samples.
	*/
	SampleUnsupervised(const std::vector<double>& Features, const std::string& ID = "");

	/**
	*\SampleUnsupervised = operator
	*\param: Rhs is a const reference to SampleUnsupervised, which "this" will be equal to.
	*/
	SampleUnsupervised& operator=(const SampleUnsupervised& Rhs);


public:
	std::vector<double>* mpFeatures;///< Pointer to a vector of features/coordinates
	std::string* mpID;///< std::string variable used to keep track of specific samples if needed. Default value is ""
};


// ------------------------------------------------------------------------------------
class SampleSetUnsupervised
{
public:
	/**
	*SampleSetUnsupervised Constructor
	*/
	SampleSetUnsupervised();

	/**
	*SampleSetUnsupervised destructor
	*/
	~SampleSetUnsupervised();

	/**
	*SampleSetUnsupervised copy constructor
	*\param: Rhs is a const reference to SampleSetUnsupervised
	*/
	SampleSetUnsupervised(const SampleSetUnsupervised& Rhs);

	/**
	*SampleSetUnsupervised assignment operator
	*\param: Rhs is a const reference to SampleSetUnsupervised
	*/
	SampleSetUnsupervised& operator=(const SampleSetUnsupervised& Rhs);

	/**
	*\fn: AddSample
	*\brief: Call this method to add sample to the training set.
	*\param: s is a const reference to SampleUnsupervised.
	*\return: Returns false if sample doesn't have the same number of features as previously added samples.
	*/
	bool AddSample(const SampleUnsupervised& s);

	/**
	*\fn: FeatureScaleAndMeanNormalize
	*\brief: Call this function every time when you have to add a new sample to the training set.
	*\param: Pass the range and average of the sample to apply feature scaling
	*/
	void FeatureScaleAndMeanNormalize(const std::vector<Common::RangeAndAverage>& RangesAndAverages);

	/**
	*\fn: Samples
	*\brief: Call this function to get a reference to this set's samples.
	*\param: A constant reference to a vector of SampleUnsupervised object, found in this set.
	*/
	const std::vector<SampleUnsupervised>& Samples() const { return *mpSamples; }

	/**
	*\fn: Samples
	*\brief: Computes the minimum and maximum values per coordinate among all samples in this set.
	*\param: "Min" is used to store the minimum value of all coordinates. Treat as a returned value.
	*\param: "Min" is used to store the maximum value of all coordinates. Treat as a returned value.
	*/
	void MinimumsAndMaximums(SampleUnsupervised& Min, SampleUnsupervised& Max) const;

	/**
	*\fn: RetainedVariancePerDimensionCount
	*\brief: The retained variance value depends on the # of coordinates when apply PCA using the "ReduceDimensionsUsingPrincipleComponentAnalysis" function.
	*		This functions computes the retained variance per dimensions count. If the samples have 4 coordinates, then it will return the retained variances for
	*		0, 1, 2, 3 & 4 dimensions.
	*\return: Returns a vector of the retained varainces, each for a different dimension count.
	*/
	std::vector<double> RetainedVariancePerDimensionCount();

	/**
	*\fn: ComputeRangesAndAverages
	*\brief: Call this method to compute range and average samples in multiple sets.
	*\param: TrainingSetsSupervised is a const reference of a vector of SampleSetUnsupervised.
	*\return: Returns a vector of RangeAndAverage.
	*/
	static std::vector<Common::RangeAndAverage> ComputeRangesAndAverages(const std::vector<SampleSetUnsupervised>& TrainingSetsUnsupervised);

	// 
	/**
	*\fn: Averages
	*\brief: Compute averages of the samples, where each value is the average of a specific cooridnate among all samples.
	*\param: The vector of the samples for which you want to compute the averages.
	*\return: A vector containing the average per coordinate among all samples.
	*/
	static std::vector<double> Averages(const std::vector<const SampleUnsupervised*>& Samples);

private:
	std::vector<SampleUnsupervised>* mpSamples;
};

// ------------------------------------------------------------------------------------