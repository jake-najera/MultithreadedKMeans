/*-------------------------------------------------------
   Filename:           driver.cpp

	 Author:             Jake Najera

	 Date:               2/29/2020
---------------------------------------------------------*/
#include <iostream>

#include "K_meansClustering.h"
#include "MTCU_Clustering.h"
#include "ParallelClustering.h"

#include <time.h> /*  time  */
#include <cstdio> /* sscanf */

#pragma warning(disable:4996)

#define SEED 20192020
#define ITERATIONS 10000

#define SINGLE 0
#define MTCU 1 // Multithreaded Centroid Updating
#define MTWP 2 // Multithreaded Work Partitioning


using std::cout;
using std::endl;

static int numOfClusters = 0;
static int numOfThreads = 0;

int main0(const SampleSetUnsupervised& SSU);
int main1(const SampleSetUnsupervised& SSU);
int main2(const SampleSetUnsupervised& SSU);

int main(int argc, char** argv)
{
	// Five Arguments
	// What flavor of threading? 0-2
	//   If 2: take in how many threads to use
	// Number of Clusters K
	// Number of Samples N
	// Minimum of Random Values Min
	// Maximum of Random Values Max

	
	//Read in cmd args
	int multithread = 0;
	std::sscanf(argv[1], "%i", &multithread);

	if (multithread == SINGLE || multithread == MTCU)
	{
		if (argc != 6)
		{
			std::cout << "Need exactly 5 arguments, read driver.cpp for specifics" << std::endl;
			return 1;
		}
	}
	else
	{
		if (multithread == MTWP)
		{
			if (argc != 7)
			{
				std::cout << "Need exactly 6 arguments, read driver.cpp for specifics" << std::endl;
				return 1;
			}
		}
	}
	

	int numOfSamples = 0;
	int minimum = 0;
	int maximum = 0;

	if (multithread == MTWP)
	{
		std::sscanf(argv[2], "%i", &numOfThreads);
		std::sscanf(argv[3], "%i", &numOfClusters);
		std::sscanf(argv[4], "%i", &numOfSamples);
		std::sscanf(argv[5], "%d", &minimum);
		std::sscanf(argv[6], "%d", &maximum);
	}
	else // MTCU
	{
		std::sscanf(argv[2], "%i", &numOfClusters);
		std::sscanf(argv[3], "%i", &numOfSamples);
		std::sscanf(argv[4], "%d", &minimum);
		std::sscanf(argv[5], "%d", &maximum);
	}

	std::default_random_engine generator(SEED); //Just for data collection purposes keeping this the same
  std::uniform_real_distribution<double> distribution(minimum, maximum);

	// Create a SampleSetUnsupervised object, and add random 2D samples
	SampleSetUnsupervised sampleSet;

	for (int i = 0; i < numOfSamples; i++)
	{
		SampleUnsupervised sample;
		
		sample.mpFeatures->push_back(distribution(generator));
		distribution.reset();
		sample.mpFeatures->push_back(distribution(generator));
		distribution.reset();

		sampleSet.AddSample(sample);
	}

	// Feature processing is done once here, before passing the samples to be clustered
	auto rangeAndAverages = SampleSetUnsupervised::ComputeRangesAndAverages({ sampleSet });
	sampleSet.FeatureScaleAndMeanNormalize(rangeAndAverages);

	if (multithread == SINGLE) 
	{
		return main0(sampleSet);
	}
	else if (multithread == MTCU) 
	{
		return main1(sampleSet);
	}
	else if (multithread == MTWP) 
	{
		return main2(sampleSet);
	}
}

// Single threaded approach
int main0(const SampleSetUnsupervised& SSU)
{
	K_meansClustering* kmeansCluster = new K_meansClustering(SSU, numOfClusters);
	cout << "Start Single-Threaded" << endl;
	for (unsigned int it = 0; it < ITERATIONS; ++it)
	{
		kmeansCluster->Iterate();
	}
	cout << "End Single-Threaded" << endl;
	//kmeansCluster->PrintInfo();

	delete kmeansCluster;
	return 0;
}

// Updating centroid coordinates on threads
// Suffers from parent thread bottlenecking!!!
// Multi-threaded step two of the algorithm
/*
		Assignment
  		  | 
  		/ | \
  	 |  |  |  Updating Centroid Coords.
  		\ | /
  		  |
	  Assignment
*/
int main1(const SampleSetUnsupervised& SSU)
{

	MTCU_Clustering* mtcCluster = new MTCU_Clustering(SSU, numOfClusters);
	cout << "Start Multi-Threaded Centroid Updating" << endl;
	mtcCluster->Iterate(ITERATIONS); //Impactful on LARGE sample size (20000+), thread per centroid
	cout << "End Multi-Threaded Centroid Updating" << endl;
	//mtcCluster->PrintInfo();

	delete mtcCluster;
	return 0;
}

// Paralellism by divying both steps in algorithm to workers
// Using ranges approach
// Multi-thread both step one and two of the algorithm
/*
    Delegates slices of samples array to T worker threads, to assign to centroids
    Barrier
    Delegates slices of centroids array to T worker threads, to update coordinates
    Barrier
    
    | | |
    V V V
    -----   Done assigning
    | | |
    V V V
    -----   Done updating centroid coordinates
*/
int main2(const SampleSetUnsupervised& SSU)
{
	ParallelClustering* paraCluster = new ParallelClustering(SSU, numOfClusters, numOfThreads);
	//The lower the cost, the better the clustering
	cout << "Start Paralellism with Workers" << endl;
	paraCluster->Execute(ITERATIONS);
	cout << "End Paralellism with Workers" << endl;
	//paraCluster->PrintInfo(); for fun

	delete paraCluster;
	return 0;
}

