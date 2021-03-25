#include "Centroid.h"
// ------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------

Centroid::Centroid(unsigned int NumberOfCoordinates)
{
	mpCoordinates = new std::vector<double>();
	mpCoordinates->resize(NumberOfCoordinates, 0.0);
}

// ------------------------------------------------------------------------------------

Centroid::~Centroid()
{
	delete mpCoordinates;
}

// ------------------------------------------------------------------------------------

Centroid::Centroid(const Centroid& Rhs)
{
	mpCoordinates = new std::vector<double>();
	*mpCoordinates = *Rhs.mpCoordinates;
}

// ------------------------------------------------------------------------------------

Centroid& Centroid::operator=(const Centroid& Rhs)
{
	if (this != &Rhs)
		*mpCoordinates = *Rhs.mpCoordinates;

	return *this;
}

// ------------------------------------------------------------------------------------

void Centroid::SetCoordinates(const std::vector<double>& NewCoordinates)
{
	if (mpCoordinates->size() == NewCoordinates.size())
		*mpCoordinates = NewCoordinates;
}