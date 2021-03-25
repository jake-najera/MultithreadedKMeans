/*-------------------------------------------------------
File Name: ParallelClustering.cpp

Purpose: Implementation of ParallelClustering for unsupervised learning

Author: Jake Najera

Creation date: 3/17/2020
---------------------------------------------------------*/
#include <iostream>
#include <limits>
#include "ParallelClustering.h"
#include "K_meansClustering.h"
#include "Common.h"

using std::cout;
using std::endl;

static int iter = 0;
static std::mutex barrierMutex;
static std::condition_variable cv;
static unsigned int barrierThreads = 0;
static std::vector<std::vector<const SampleUnsupervised* >> mpSamplesPerCentroid;
static bool go = false;
static bool done = false;
static bool allAssigned = false;
// ------------------------------------------------------------------------------------
ParallelClustering::ParallelClustering(const SampleSetUnsupervised& TS, unsigned int NumberOfCentroids, unsigned int NumberOfThreads)
{
	mpCentroids = new std::vector<Centroid>();
	//mpSamplesPerCentroid = new std::vector<std::vector<const SampleUnsupervised* >>(NumberOfCentroids);
	//mpCentroidPerSample = new std::unordered_map<const SampleUnsupervised*, const Centroid*>();

	mTrainingSet = TS;
	size_t numberOfFeatures = 0;
	if (TS.Samples().empty() == false)
		numberOfFeatures = TS.Samples()[0].mpFeatures->size();

	numOfThreads = NumberOfThreads;
	numOfSamples = TS.Samples().size();
	iterations = 0;

	Centroid c((unsigned int)numberOfFeatures);
	mpCentroids->resize(NumberOfCentroids, c);
	mpSamplesPerCentroid.resize(NumberOfCentroids);
	//mpSamplesPerCentroid->resize(NumberOfCentroids);


	// Random positions for centroids
	// Find min & max of each coordinates
	if (0 != NumberOfCentroids && 0 != numberOfFeatures)
	{
		SampleUnsupervised minimums, maximums;
		TS.MinimumsAndMaximums(minimums, maximums);

		// Random for now
		std::vector<double> coordinates;
		coordinates.resize(numberOfFeatures, 0.0);
		for (auto& centroid : *mpCentroids)
		{
			for (size_t coordIdx = 0; coordIdx < numberOfFeatures; ++coordIdx)
			{
				double range = maximums.mpFeatures->at(coordIdx) - minimums.mpFeatures->at(coordIdx);
				double percentageOfRange = Common::RandomDouble();
				coordinates[coordIdx] = minimums.mpFeatures->at(coordIdx) + percentageOfRange * range;
			}
			centroid.SetCoordinates(coordinates);
		}
	}
}

// ------------------------------------------------------------------------------------

ParallelClustering::~ParallelClustering()
{
	delete mpCentroids;
	//delete mpSamplesPerCentroid;
	//delete mpCentroidPerSample;
}

// ------------------------------------------------------------------------------------

const Centroid* ParallelClustering::ClosestCentroidToSample_ByValue(const SampleUnsupervised& Sample) const
{
	size_t indexOfChosen = ClosestCentroidToSample_Internal(Sample);
	if (SIZE_MAX != indexOfChosen)
		return &mpCentroids->at(indexOfChosen);

	return nullptr;
}

// ------------------------------------------------------------------------------------

unsigned int ParallelClustering::ClosestCentroidToSample_ByIndex(const SampleUnsupervised& Sample) const
{
	return (unsigned int)ClosestCentroidToSample_Internal(Sample);
}

// ------------------------------------------------------------------------------------

void ParallelClustering::Execute(int N)
{
	iterations = N;

	size_t numOfSamples = mTrainingSet.Samples().size();
	size_t numOfCentroids = mpCentroids->size();

	std::vector<std::thread> threadPool;
	std::vector<ThreadData*> dataPool;
	dataPool.resize(numOfThreads);



	//Precalculate ranges for both samples and centroids!!!
	for (size_t i = 0; i < numOfThreads; i++)
	{
		ThreadData* data = new ThreadData;

		size_t sampleBegin = (i * numOfSamples) / numOfThreads;
		size_t sampleEnd = (((i + 1) * numOfSamples) / numOfThreads);
		size_t centroidBegin = (i * numOfCentroids) / numOfThreads;
		size_t centroidEnd = (((i + 1) * numOfCentroids) / numOfThreads);

		data->sampleRange = std::make_pair(sampleBegin, sampleEnd);
		data->centroidRange = std::make_pair(centroidBegin, centroidEnd);
		dataPool[i] = data;

		threadPool.push_back(std::thread(&ParallelClustering::ThreadWork, (this), dataPool[i]));
	}

	for (int i = 0; i < numOfThreads; i++)
	{
		threadPool[i].join();
	}

	for (int i = 0; i < numOfThreads; i++)
	{
		delete dataPool[i];
	}
	
}

// ------------------------------------------------------------------------------------

void ParallelClustering::ThreadWork(const ThreadData* info)
{
	size_t sampleBegin = info->sampleRange.first;
	size_t sampleEnd = info->sampleRange.second;
	size_t centroidBegin = info->centroidRange.first;
	size_t centroidEnd = info->centroidRange.second;

	//for N iterations
	while (!done)
	{
		AssignSamplesToTheirClosestCentroid(sampleBegin, sampleEnd); 

	  //Barrier
		Barrier();

    UpdateCentroidCoordinates(centroidBegin, centroidEnd);

		//Barrier
		Barrier();
	}
	
	

	//Think of GOL approach to threading!
}

// ------------------------------------------------------------------------------------

void ParallelClustering::AssignSamplesToTheirClosestCentroid(unsigned int begin, unsigned int end)
{

	size_t numOfSamples = mTrainingSet.Samples().size();
	size_t numOfCentroids = mpCentroids->size();

	std::vector < std::vector<const SampleUnsupervised*> > assignments; //Index is directly related to the closest centroid 
	assignments.resize(numOfCentroids);

	for (int i = begin; i < end; i++)
	{
		const SampleUnsupervised& sample = mTrainingSet.Samples().at(i); //Get i-th sample
		double closest = std::numeric_limits<double>::max();
		int closestCentroidIndex = -1;

		//Loop through centroids to find shortest distance
		for (int j = 0; j < numOfCentroids; j++)
		{
			const Centroid& centroid = mpCentroids->at(j);
			double distance = Common::SquaredDistance(*sample.mpFeatures, centroid.Coordinates());
			if (distance < closest)
			{
				closest = distance;
				closestCentroidIndex = j;
			}
		}

		const Centroid& c = mpCentroids->at(closestCentroidIndex);
		assignments[closestCentroidIndex].push_back(&sample);
	}

	for (int i = 0; i < numOfCentroids; i++)
	{
		const Centroid& c = mpCentroids->at(i);
		SPCMutex.lock();
		mpSamplesPerCentroid[i].insert(mpSamplesPerCentroid[i].end(), assignments[i].begin(), assignments[i].end());
		SPCMutex.unlock();
	}
}

// ------------------------------------------------------------------------------------

void ParallelClustering::UpdateCentroidCoordinates(unsigned int begin, unsigned int end)
{
	size_t numOfCentroids = mpCentroids->size();
	size_t numOfCoords = mpCentroids->at(0).Coordinates().size();
	std::vector<double> averagesOfCoords;
	averagesOfCoords.resize(numOfCoords);

	// Per centroid:
	//  • Compute the average position of the samples that were assigned to it
	//  • Set the result as the centroid’s position
	for (int i = begin; i < end; i++)
	{
		Centroid& centroid = mpCentroids->at(i);
		double average = 0;

		if (!mpSamplesPerCentroid[i].empty())
		{
			for (auto& sample : mpSamplesPerCentroid[i])
			{
				for (int k = 0; k < numOfCoords; k++)
				{
					averagesOfCoords[k] += sample->mpFeatures->at(k);
				}
			}

			for (int k = 0; k < numOfCoords; k++)
			{
				averagesOfCoords[k] = averagesOfCoords[k] / mpSamplesPerCentroid[i].size();
			}

			centroid.SetCoordinates(averagesOfCoords);
		}
	}
}

// ------------------------------------------------------------------------------------

void ParallelClustering::Barrier()
{
	std::unique_lock<std::mutex> lock(barrierMutex);
	barrierThreads++;
	if (barrierThreads == numOfThreads)
	{
		barrierThreads = 0;
		if (!allAssigned)
		{
			allAssigned = true;
		}
		else
		{
			iter++;
			if (iterations == iter)
			{
				done = true;
			}
			else
			{
				allAssigned = false;
				//mpCentroidPerSample->clear();

				for (std::vector<const SampleUnsupervised* >& samplesToCentroid : mpSamplesPerCentroid)
				{
					samplesToCentroid.clear();
				}
			}
		}
		
		cv.notify_all();
	}
	else
	{
		cv.wait(lock);
	}
}

// ------------------------------------------------------------------------------------

double ParallelClustering::Cost()
{
	//Distortion fxn
	size_t numOfSamples = mTrainingSet.Samples().size();
	double cost = 0.0;

	for (int i = 0; i < numOfSamples; i++)
	{
		const SampleUnsupervised& sampleI = mTrainingSet.Samples().at(i);
		unsigned centroidIndex = ClosestCentroidToSample_ByIndex(sampleI);
		Centroid& assignedCentroid = mpCentroids->at(centroidIndex);

		//Calculate distance between sample and centroid
		cost += Common::SquaredDistance(*sampleI.mpFeatures, assignedCentroid.Coordinates());
	}

	cost = cost / numOfSamples;

	return cost;
}

// ------------------------------------------------------------------------------------

void ParallelClustering::PrintInfo()
{
	cout << "Cost = " << Cost() << endl;

	for (int i = 0; i < mpCentroids->size(); i++)
	{
		bool first = true;
		const Centroid& centroid = mpCentroids->at(i);

		cout << "Centroid at [";
		for (auto coord : centroid.Coordinates())
			if (true == first)
			{
				first = false;
				cout << coord;
			}
			else
				cout << " ; " << coord;
		cout << "]" << endl;

		for (auto& sample : mpSamplesPerCentroid[i])
		{
			first = true;
			cout << "[";
			for (auto feature : *sample->mpFeatures)
			{
				if (true == first)
				{
					first = false;
					cout << feature;
				}
				else
					cout << " ; " << feature;
			}
			//cout << "]\t" << Common::SquaredDistance(sample->mFeatures, centroid.Coordinates()) << endl;
			cout << "]\t" << endl;
		}

		cout << endl << endl;
	} 
}

// ------------------------------------------------------------------------------------

size_t ParallelClustering::ClosestCentroidToSample_Internal(const SampleUnsupervised& Sample) const
{
	size_t indexOfChosen = SIZE_MAX;
	double smallestSquaredDistance = std::numeric_limits<double>::max();

	size_t index, num = mpCentroids->size();
	for (index = 0; index < num; ++index)
	{
		double sqDist = Common::SquaredDistance(*Sample.mpFeatures, mpCentroids->at(index).Coordinates());
		if (sqDist < smallestSquaredDistance)
		{
			smallestSquaredDistance = sqDist;
			indexOfChosen = index;
		}
	}

	return indexOfChosen;
}

// ------------------------------------------------------------------------------------

ParallelClustering& ParallelClustering::operator=(const ParallelClustering& Rhs)
{
	return *this;
}

// ------------------------------------------------------------------------------------
