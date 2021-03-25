#pragma once

#include <vector>
#include <random>

// ------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------
/**
\namespace: Common
*/
namespace Common
{
	class RangeAndAverage
	{
	public:
		/**
		*\RangeAndAverage constructor
		*\param: Range is the difference between the maximum and minimum values.
		*\param: Average is the mean of the values.
		*/
		RangeAndAverage(double Range, double Average)
		{
			mRange = Range;
			mAverage = Average;
		}

		double mRange;
		double mAverage;
	};

	// ------------------------------------------------------------------------------------
	// ------------------------------------------------------------------------------------
	// ------------------------------------------------------------------------------------

	/**
	*\fn: DotProduct
	*\brief: Calculates the dot product of two vectors. DotProduct = x1*y1 + x2*y2 + ... + xN*yN
	*\param: Provide two vectors of doubles as input.
	*\return: The dot product result
	*/
	double DotProduct(const std::vector<double> &V1, const std::vector<double> &V2);

	/**
	*\fn: Distance
	*\brief: Calculates the distance between two vectors. Distance = Square Root[(x1-y1)^2 + (x2-y2)^2 + ... + (xN-yN)^2]
	*\param: Provide two vectors of doubles as input.
	*\return: The distance value
	*/
	double Distance(const std::vector<double> &V1, const std::vector<double> &V2);

	/**
	*\fn: SquaredDistance
	*\brief: Calculates the squared distance between two vectors. SquaredDistance = (x1-y1)^2 + (x2-y2)^2 + ... + (xN-yN)^2
	*\param: Provide two vectors of template T as input.
	*\return: The squared distance value.
	*/
	template <class T>
	T SquaredDistance(const std::vector<T> &V1, const std::vector<T> &V2)
	{
		size_t size1 = V1.size();
		if (size1 != V2.size() || 0 == size1) return 0;

		T result = 0;
		for (size_t idx = 0; idx < size1; ++idx)
		{
			T diff = V1[idx] - V2[idx];
			result += diff * diff;
		}
		return result;
	}

	/**
	*\fn: Length
	*\brief: Calculates the length of the vector. Length = Square root[(x1)^2 + (x2)^2 + ... + (xN)^2]
	*\param: A vector of doubles.
	*\return: The length value of the vector.
	*/
	double Length(const std::vector<double> &V);

	/**
	*\fn: SquaredLength
	*\brief: Calculates the square length of the vector. SquaredLength = (x1)^2 + (x2)^2 + ... + (xN)^2
	*\param: A vector of doubles.
	*\return: The squared length value of the vector.
	*/
	double SquaredLength(const std::vector<double> &V);

	/**
	*\fn: VectorSum
	*\brief: Calculate the sum of the elements in a vector.
	*\param: A vector of doubles.
	*\return: The sum of the vector's elements
	*/
	double VectorSum(const std::vector<double> &V);

	/**
	*\fn: VectorPlusVector
	*\brief: Calculates the element wise addition of two vectors. V = [x1+y1, x2+y2,....,xN+yN]
	*\param: 2 vectors of doubles, to be added.
	*\return: A vector of doubles.
	*/
	std::vector<double> VectorPlusVector(const std::vector<double> &V1, const std::vector<double> &V2);


	/**
	*\fn: VectorAdd
	*\brief: In-place element wise addition of 2 vectors. V = [x1+y1, x2+y2,....,xN+yN]
	*\param: 2 vectors of doubles, to be added. Result will be stored in V1.
	*/
  void VectorAdd(std::vector<double> &V1, const std::vector<double> &V2);


	/**
	*\fn: VectorMinusVector
	*\brief: Calculates the difference between two vectors. V = [x1-y1, x2-y2,....,xN-yN]
	*\param: 2 vectors of doubles, to be subtracted.
	*\return: A vector of doubles.
	*/
	std::vector<double> VectorMinusVector(const std::vector<double> &V1, const std::vector<double> &V2);

	/**
	*\fn: VectorScale
	*\brief: Scales a vector using a scalar value.
	*\param: A vector of doubles and the scalar value.
	*\return: The scaled vector of doubles.
	*/
	std::vector<double> VectorScale(const std::vector<double> &V1, double Scale);


	/**
	*\fn: VectorScale2
	*\brief: In-place vector scale.
	*\param: A vector of doubles and the scalar value. Result will be stored in V1.
	*/
	void VectorScale2(std::vector<double> &V1, double Scale);


	/**
	*\fn: RandomDouble
	*\brief: Generates a random double value in the range of [0...1]
	*\return: A double value between [0.0, 1.0]
	*/
	double RandomDouble();

	/**
	*\fn: RandomVector
	*\brief: This function will fill up the vector with random values ranging from 0.0 to MaxValue.
	*\param: A vector of doubles,
	*\param: The maximum value that can be generated.
	*/
	void RandomVector(std::vector<double> &V, double MaximumValue);	// [0.0 ; MaximumValue]

	/**
	*\fn: ComputeGaussianParameters
	*\brief: Calculates the Mean and Variance of a vector.
	*\param: A vector of doubles containing the values.
	*\param: Mean: Reference to a double. The mean value will be stored here.
	*\param: Variance: Reference to a double. The variance value will be stored here.
	*/
	void ComputeGaussianParameters(const std::vector<double> &Values, double &Mean, double &Variance);

	/**
	*\fn: ProbabilityDensity
	*\brief: Calculates the probability of a value.
	*\param: Value: The value to be tested
	*\param: Mean: The mean of the Gaussian distribution
	*\param: StandardDeviation: The standard deviation of the Gaussian distribution
	*\return: The probability density of Value in the distribution defined by Mean and Standard Deviation
	*/
	double ProbabilityDensity(double Value, double Mean, double StandardDeviation);


}