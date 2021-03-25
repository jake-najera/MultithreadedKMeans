/*-------------------------------------------------------
File Name: K_meansClustering.cpp

Purpose: Implementation of K_means Clustering for unsupervised learning

Author: Jake Najera

Creation date: 2/29/2020
---------------------------------------------------------*/
#include <iostream>
#include <limits>
#include "K_meansClustering.h"
#include "Common.h"

using std::cout;
using std::endl;

// ------------------------------------------------------------------------------------

K_meansClustering::K_meansClustering(const SampleSetUnsupervised& TS, unsigned int NumberOfCentroids)
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

K_meansClustering::~K_meansClustering()
{
	delete mpCentroids;
	delete mpSamplesPerCentroid;
	delete mpCentroidPerSample;
}

// ------------------------------------------------------------------------------------

const Centroid* K_meansClustering::ClosestCentroidToSample_ByValue(const SampleUnsupervised& Sample) const
{
	size_t indexOfChosen = ClosestCentroidToSample_Internal(Sample);
	if (SIZE_MAX != indexOfChosen)
		return &mpCentroids->at(indexOfChosen);

	return nullptr;
}

// ------------------------------------------------------------------------------------

unsigned int K_meansClustering::ClosestCentroidToSample_ByIndex(const SampleUnsupervised& Sample) const
{
	return (unsigned int)ClosestCentroidToSample_Internal(Sample);
}

// ------------------------------------------------------------------------------------
void K_meansClustering::Iterate()
{
	AssignSamplesToTheirClosestCentroid();
	UpdateCentroidCoordinates();
}

// ------------------------------------------------------------------------------------

double K_meansClustering::Cost()
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

void K_meansClustering::PrintInfo()
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

size_t K_meansClustering::ClosestCentroidToSample_Internal(const SampleUnsupervised& Sample) const
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

void K_meansClustering::AssignSamplesToTheirClosestCentroid()
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
		mpCentroidPerSample->insert(std::pair<const SampleUnsupervised*, const Centroid*>(&sample, &c));

	}

	for (int i = 0; i < numOfCentroids; i++)
	{
		const Centroid& c = mpCentroids->at(i);
		mpSamplesPerCentroid->insert(std::pair<const Centroid*, std::vector<const SampleUnsupervised*> >(&c, (assignments[i])));
	}
}

// ------------------------------------------------------------------------------------

void K_meansClustering::UpdateCentroidCoordinates()
{
	size_t numOfCentroids = mpCentroids->size();
	size_t numOfCoords = mpCentroids->at(0).Coordinates().size();
	std::vector<double> averagesOfCoords;
	averagesOfCoords.resize(numOfCoords);

	// Per centroid:
	//  • Compute the average position of the samples that were assigned to it
	//  • Set the result as the centroid’s position
	for (int i = 0; i < numOfCentroids; i++)
	{
		Centroid& centroid = mpCentroids->at(i);
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
}

// ------------------------------------------------------------------------------------
K_meansClustering& K_meansClustering::operator=(const K_meansClustering& Rhs)
{
	return *this;
}

// ------------------------------------------------------------------------------------
