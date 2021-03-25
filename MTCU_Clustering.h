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
class Semaphore
{
public:
	Semaphore(int threads) : m_allottedThreads(threads) {};
	void wait()
	{
		std::unique_lock<std::mutex> lock(m_mutex);
		m_allottedThreads--;


		if (m_allottedThreads < 0)
		{
			do {
				m_cond.wait(lock);
			} while (m_wakeups < 1);
			m_wakeups--;
		}

		lock.unlock();
	};

	void signal()
	{
		m_mutex.lock();
		m_allottedThreads++;

		if (m_allottedThreads <= 0)
		{
			m_wakeups++;
			m_cond.notify_one();
		}

		m_mutex.unlock();
	};

private:
	int m_allottedThreads = 0;
	int m_wakeups = 0;
	std::mutex m_mutex;
	std::condition_variable m_cond;
};

class MTCU_Clustering
{
public:
	/**
	*\MTCU_Clustering constructor.
	*\param: Unsupervised data set, which will be used to adjust the centroids' position in order to reduce the cost.
	*\param: The required number of centroids.
	*/
	MTCU_Clustering(const SampleSetUnsupervised& TS, unsigned int NumberOfCentroids);

	/**
	*\MTCU_Clustering destructor.
	*/
	~MTCU_Clustering();

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
	*\brief: Call this function to perform iterations. With every iteration step, the cost should decrease.
	*/
	void Iterate(int interations);

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
	unsigned int CentroidsNum() const { return (unsigned int)mpCentroids->size(); };

	void AssignSamplesToTheirClosestCentroid();
  void MultiUpdateCentroids(int centroidIndex);
	void ThreadWork(int index);
	void Barrier();

private:
	size_t ClosestCentroidToSample_Internal(const SampleUnsupervised& Sample) const;
	// Not allowed
	MTCU_Clustering(const MTCU_Clustering& Rhs) = delete;
	MTCU_Clustering& operator=(const MTCU_Clustering& Rhs);

private:
	SampleSetUnsupervised mTrainingSet;
	std::vector<Centroid>* mpCentroids;
	std::unordered_map<const Centroid*, std::vector<const SampleUnsupervised*>>* mpSamplesPerCentroid;
	std::unordered_map<const SampleUnsupervised*, const Centroid*>* mpCentroidPerSample;
	//bool updateFlag = false;
	//bool assignFlag = false;

};