#pragma once

#include "Common.h"
#include "SampleUnsupervised.h"
#include "K_meansClustering.h"
#include <unordered_map>
#include <thread>
#include <mutex>
#include <atomic>
#include <future>
#include <condition_variable>


// ------------------------------------------------------------------------------------

struct ThreadData
{
	std::pair<size_t, size_t> sampleRange;
	std::pair<size_t, size_t> centroidRange;
};


class  ParallelClustering
{
public:
	/**
	*\ParalellClustering constructor.
	*\param: Unsupervised data set, which will be used to adjust the centroids' position in order to reduce the cost.
	*\param: The required number of centroids.
	*/
	ParallelClustering(const SampleSetUnsupervised& TS, unsigned int NumberOfCentroids, unsigned int NumberOfThreads);

	/**
	*\ParallelClustering destructor.
	*/
	~ParallelClustering();

	/**
	*\fn: GetCentroids.
	*\brief: Get a constant pointer to this object's vector of centroids.
	*\return: A constant vector of centroids. Can be used to track the centroids' positions.
	*/
	std::vector<Centroid>* GetCentroids() const { return mpCentroids; }

	/**
	*\fn: ClosestCentroidToSample_ByValue.
	*\brief: Get the closest centroid to a test sample.
	*\param: The sample to find its closest centroid.
	*\return: A constant pointer to the closest centroid to "Sample"
	*/
	const Centroid* ClosestCentroidToSample_ByValue(const SampleUnsupervised& Sample) const;

	/**
	*\fn: ClosestCentroidToSample_ByIndex.
	*\brief: Get the closest centroid to a test sample.
	*\param: The sample to find its closest centroid.
	*\return: The index of the closest centroid to "Sample", which can be used to index the vector returned by "GetCentroids"
	*/
	unsigned int ClosestCentroidToSample_ByIndex(const SampleUnsupervised& Sample) const;

	/**
	*\fn: Iterate.
	*\brief: Call this function to perform a single iteration. With every iteration step, the cost should decrease.
	*/
	void Execute(int iterations);

	/**
	*\fn: Cost.
	*\brief: Get the current cost, based on the centroids' current positions.
	*\return: The current cost, based on the samples and the current centroids' positions.
	*/
	double Cost();

	/**
	*\fn: PrintInfo.
	*\brief: Use this function to print total cost, centroids and elements in each cluster.
	*/
	void PrintInfo();

	/**
	*\fn: CentroidsNum.
	*\brief: Get the number of centroids.
	*\return: returns the total number of centroids/clusters in this object.
	*/
	unsigned int CentroidsNum() const { return (unsigned int)mpCentroids->size(); }

	void ThreadWork(const ThreadData* data);
	void Barrier();

	void AssignSamplesToTheirClosestCentroid(unsigned int begin, unsigned int end);
	void UpdateCentroidCoordinates(unsigned int begin, unsigned int end);

private:
	size_t ClosestCentroidToSample_Internal(const SampleUnsupervised& Sample) const;
	// Not allowed
	ParallelClustering(const ParallelClustering& Rhs) = delete;
	ParallelClustering& operator=(const ParallelClustering& Rhs);

private:
	SampleSetUnsupervised mTrainingSet;
	std::vector<Centroid>* mpCentroids;
	//std::unordered_map<const SampleUnsupervised*, const Centroid*>* mpCentroidPerSample;
	std::mutex CPSMutex;
	std::mutex SPCMutex;
	unsigned int numOfThreads;
	unsigned int numOfSamples;
	unsigned int iterations;
};

