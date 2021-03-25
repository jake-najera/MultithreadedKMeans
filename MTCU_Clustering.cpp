/*-------------------------------------------------------
File Name: MTCU_Clustering.cpp

Purpose: Implementation of MTCU_Clustering

Author: Jake Najera

Creation date: 3/4/2020
---------------------------------------------------------*/
#include <iostream>
#include <limits>
#include "MTCU_Clustering.h"
#include "Common.h"

using std::cout;
using std::endl;

std::mutex barrierMutex;
std::mutex updateMutex;
std::mutex assignMutex;
std::condition_variable assignCV;
std::condition_variable cv;
std::condition_variable updateCV;
unsigned int barrierThreads = 0;
int iterations = 0;
int numberOfIterations = 0;
bool done = false;


MTCU_Clustering::MTCU_Clustering(const SampleSetUnsupervised& TS, unsigned int NumberOfCentroids)
{
	mpCentroids = new std::vector<Centroid>();
	mpSamplesPerCentroid = new std::unordered_map<const Centroid*, std::vector<const SampleUnsupervised*>>();
	mpCentroidPerSample = new std::unordered_map<const SampleUnsupervised*, const Centroid*>();

	mTrainingSet = TS;
	size_t numberOfFeatures = 0;
	if (TS.Samples().empty() == false)
		numberOfFeatures = TS.Samples()[0].mpFeatures->size();

	Centroid c((unsigned int)numberOfFeatures);
	mpCentroids->resize(NumberOfCentroids, c);

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

	//AssignSamplesToTheirClosestCentroid();
}

// ------------------------------------------------------------------------------------

MTCU_Clustering::~MTCU_Clustering()
{
	delete mpCentroids;
	delete mpSamplesPerCentroid;
	delete mpCentroidPerSample;
}

// ------------------------------------------------------------------------------------

const Centroid* MTCU_Clustering::ClosestCentroidToSample_ByValue(const SampleUnsupervised& Sample) const
{
	size_t indexOfChosen = ClosestCentroidToSample_Internal(Sample);
	if (SIZE_MAX != indexOfChosen)
		return &mpCentroids->at(indexOfChosen);

	return nullptr;
}

// ------------------------------------------------------------------------------------

unsigned int MTCU_Clustering::ClosestCentroidToSample_ByIndex(const SampleUnsupervised& Sample) const
{
	return (unsigned int)ClosestCentroidToSample_Internal(Sample);
}

// ------------------------------------------------------------------------------------

/*
    * Starting at parent thread, then branching to N threads, N = number of clusters
    * Updating centroid coordinates using the threads
    * Meet afterward
    
    											Parent
    												|
    											/ | \  Child threads
    											\ | /  updating centroids
    												|
    											Parent
    
    * Approach is reliant on parent performance for speedup, since parent handles
    	step one of algorithm; assigning samples to centroids
*/
void MTCU_Clustering::Iterate(int iterations)
{
	numberOfIterations = iterations;
	size_t numOfSamples = mTrainingSet.Samples().size();
	size_t numOfCentroids = mpCentroids->size();
	std::vector <std::thread> threadPool;
	
	AssignSamplesToTheirClosestCentroid();

	//while (!done)
	//{
	//	MultiUpdateCentroids(0);
	//	Barrier();
	//}

	for (int i = 0; i < numOfCentroids; i++)
	{
		threadPool.push_back(std::thread(&MTCU_Clustering::ThreadWork, (this), i));
	}

	for (int i = 0; i < numOfCentroids; i++)
	{
		threadPool[i].join();
	}
	//cout << "done" << endl;
}

void MTCU_Clustering::ThreadWork(int index)
{
	while (!done)
	{
		MultiUpdateCentroids(index);
		Barrier();
	}
}

void MTCU_Clustering::Barrier()
{
	std::unique_lock<std::mutex> lock(barrierMutex);
	barrierThreads++;
	if (barrierThreads == mpCentroids->size())
	{
		barrierThreads = 0;
		iterations++;
		if (iterations == numberOfIterations)
		{
			done = true;
		}
		AssignSamplesToTheirClosestCentroid();
		cv.notify_all();
	}
	else
	{
		cv.wait(lock);
	}
}

// ------------------------------------------------------------------------------------

double MTCU_Clustering::Cost()
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

void MTCU_Clustering::PrintInfo()
{
	cout << "Cost = " << Cost() << endl;
	for (auto& centroid : *mpCentroids)
	{
		bool first = true;
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

		for (auto& sample : mpSamplesPerCentroid->at(&centroid))
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

size_t MTCU_Clustering::ClosestCentroidToSample_Internal(const SampleUnsupervised& Sample) const
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

void MTCU_Clustering::AssignSamplesToTheirClosestCentroid()
{
	size_t numOfSamples = mTrainingSet.Samples().size();
	size_t numOfCentroids = mpCentroids->size();

	mpCentroidPerSample->clear();
	mpSamplesPerCentroid->clear();

	std::vector < std::vector<const SampleUnsupervised*> > assignments; //Index is directly related to the closest centroid 
	assignments.resize(numOfCentroids);

	for (int i = 0; i < numOfSamples; i++)
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
		//mpCentroidPerSample->insert(std::pair<const SampleUnsupervised*, const Centroid*>(&sample, &c));

	}

	for (int i = 0; i < numOfCentroids; i++)
	{
		const Centroid& c = mpCentroids->at(i);
		mpSamplesPerCentroid->insert(std::pair<const Centroid*, std::vector<const SampleUnsupervised*> >(&c, (assignments[i])));
	}
}

// ------------------------------------------------------------------------------------

void MTCU_Clustering::MultiUpdateCentroids(int centroidIndex)
{
	size_t numOfCoords = mpCentroids->at(0).Coordinates().size();
	std::vector<double> averagesOfCoords;
	averagesOfCoords.resize(numOfCoords);

	// Per centroid:
	//  • Compute the average position of the samples that were assigned to it
	//  • Set the result as the centroid’s position
	Centroid& centroid = mpCentroids->at(centroidIndex);
	double average = 0;

	if (!mpSamplesPerCentroid->at(&centroid).empty())
	{
		for (auto& sample : mpSamplesPerCentroid->at(&centroid))
		{
			for (int k = 0; k < numOfCoords; k++)
			{
				averagesOfCoords[k] += sample->mpFeatures->at(k);
			}
		}

		for (int k = 0; k < numOfCoords; k++)
		{
			averagesOfCoords[k] = averagesOfCoords[k] / mpSamplesPerCentroid->at(&centroid).size();
		}

		centroid.SetCoordinates(averagesOfCoords);
	}
}

// ------------------------------------------------------------------------------------

MTCU_Clustering& MTCU_Clustering::operator=(const MTCU_Clustering& Rhs)
{
	return *this;
}

// ------------------------------------------------------------------------------------
