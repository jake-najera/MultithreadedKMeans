#include "Common.h"
#include <time.h>
#include <iostream>

#define SQRT_OF_2_PI 2.506628274631000



// ------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------

namespace Common
{
	double DotProduct(const std::vector<double> &V1, const std::vector<double> &V2)
	{
		size_t size1 = V1.size();
		if (V2.size() != size1) return 0.0;

		double result = 0.0;
		for (size_t idx = 0; idx < size1; ++idx)
			result += V1[idx] * V2[idx];
		return result;
	}

	// ------------------------------------------------------------------------------------
	
	double Distance(const std::vector<double> &V1, const std::vector<double> &V2)
	{
		return sqrt(SquaredDistance(V1, V2));
	}

	// ------------------------------------------------------------------------------------
	
	double Length(const std::vector<double> &V)
	{
		return sqrt(SquaredLength(V));
	}

	// ------------------------------------------------------------------------------------

	double SquaredLength(const std::vector<double> &V)
	{
		double result = 0;
		for (auto value : V)
			result += value * value;
		return result;
	}

	// ------------------------------------------------------------------------------------
	
	double VectorSum(const std::vector<double> &V)
	{
		double sum = 0;
		for (auto v : V)
			sum += v;
		return sum;
	}
	
	// ------------------------------------------------------------------------------------

	std::vector<double> VectorPlusVector(const std::vector<double> &V1, const std::vector<double> &V2)
	{
		std::vector<double> result;
		size_t size1 = V1.size();
		if (V2.size() != size1)
			return result;

		result.resize(size1);
		for (size_t idx = 0; idx < size1; ++idx)
			result[idx] = V1[idx] + V2[idx];
		return result;
	}

	// ------------------------------------------------------------------------------------

	void VectorAdd(std::vector<double> &V1, const std::vector<double> &V2)
	{
		size_t size1 = V1.size();
		if (V2.size() != size1)
			return;

		for (size_t idx = 0; idx < size1; ++idx)
			V1[idx] += V2[idx];
	}

	// ------------------------------------------------------------------------------------

	std::vector<double> VectorMinusVector(const std::vector<double> &V1, const std::vector<double> &V2)
	{
		std::vector<double> result;
		size_t size1 = V1.size();
		if (V2.size() != size1)
			return result;

		result.resize(size1);
		for (size_t idx = 0; idx < size1; ++idx)
			result[idx] = V1[idx] - V2[idx];
		return result;
	}

	// ------------------------------------------------------------------------------------

	std::vector<double> VectorScale(const std::vector<double> &V1, double Scale)
	{
		std::vector<double> result = V1;

		for (auto &v : result)
			v *= Scale;

		return result;
	}

	// ------------------------------------------------------------------------------------

	void VectorScale2(std::vector<double> &V1, double Scale)
	{
		for (auto &v : V1)
			v *= Scale;
	}

	// ------------------------------------------------------------------------------------

	double RandomDouble()
	{
		return (double)(rand() / (float)(RAND_MAX));		// [0.0 ; 1.0] range
	}

	// ------------------------------------------------------------------------------------

	void RandomVector(std::vector<double> &V, double MaximumValue)	// [0.0 ; MaximumValue]
	{
		for (auto &v : V)
			v = RandomDouble() * MaximumValue;
	}

	// ------------------------------------------------------------------------------------

	// [STUDENT]
	void ComputeGaussianParameters(const std::vector<double> &Values, double &Mean, double &Variance)
	{
	}

	// ------------------------------------------------------------------------------------

	// [STUDENT]
	double ProbabilityDensity(double Value, double Mean, double StandardDeviation)
	{
		return 0.0;
	}

	// ------------------------------------------------------------------------------------

}