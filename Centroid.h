#pragma once

#include <vector>

class Centroid
{
public:
	/**
	*\Centroid constructor
	*\param: Number if needed coordinates, which should match the number of coordinates of each sample in the K or C means clustering object.
	*/
	Centroid(unsigned int NumberOfCoordinates);

	/**
	*\Centroid destructor
	*/
	~Centroid();

	/**
	*\Centroid copy constructor
	*\param: The object to copy from
	*/
	Centroid(const Centroid& Rhs);

	/**
	*\Centroid = operation.
	*\brief: Copies a centroid's data to another one.
	*\param: The centroid to copy from.
	*\return: The resultant centroid.
	*/
	Centroid& operator=(const Centroid& Rhs);


	/**
	*\fn: Coordinates
	*\brief: Get the current coordinates of this centroid object.
	*\return: Returns a copy of this centroid's current coordinates.
	*/
	const std::vector<double>& Coordinates() const { return *mpCoordinates; }

	/**
	*\fn: SetCoordinates
	*\brief: Set the coordinates of this centroid object. Used by the K/C means algorithms to adjust the centroids' position in order to reduce cost.
	*/
	void SetCoordinates(const std::vector<double>& NewCoordinates);

private:
	std::vector<double>* mpCoordinates;
};
