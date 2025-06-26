#pragma once
#include <math.h>
#include <algorithm>

#define my_PI 3.141592653589793238463

/********************math*******************/
template<typename T, typename Container>
inline T dot(Container& a, Container& b, int d)
{

	T re = 0;
	for (int i = 0; i < d; ++i)
	{
		re += a[i] * b[i];
	}
	return re;
}

template<typename T>
inline T cross2D(T a[], T b[])
{
	return a[0] * b[1] - a[1] * b[0];
}

template <class T>
inline void cross3D(const T* e1, const T* e2, T* pn) 
{
	pn[0] = e1[1] * e2[2] - e2[1] * e1[2];
	pn[1] = e1[2] * e2[0] - e2[2] * e1[0];
	pn[2] = e1[0] * e2[1] - e2[0] * e1[1];
}

template<typename T, typename Container>
inline T length(Container& a, int d) 
{
	T re = 0;
	for (int i = 0; i < d; ++i) 
	{
		re += a[i] * a[i];
	}
	return sqrt(re);
}

template<typename T, typename Container>
inline void normalize(Container& a, int d) 
{
	T len = length<T>(a, d);
	if (len == 0) 
	{
		for (int i = 0; i < d; ++i) 
		{
			a[i] = 0;
		}
	}
	else {
		for (int i = 0; i < d; ++i) 
		{
			a[i] = a[i] / len;
		}
	}
}

template <typename T, typename Container>
inline void standardize(Container& vec)
{
	if (vec.empty())
		return; 

	T minVal = *min_element(vec.begin(), vec.end());
	T maxVal = *max_element(vec.begin(), vec.end());

	if (minVal == maxVal)
		return;

	for (size_t i = 0; i < vec.size(); ++i)
	{
		vec[i] = (vec[i] - minVal) / (maxVal - minVal);
	}
}



template <typename T>
T optimized_pow(T base, int exp) 
{
	if (exp < 0)
		return T(1) / optimized_pow(base, -exp);
	T result = T(1);
	while (exp > 0) {
		if (exp % 2 == 1)
			result *= base;
		base *= base;
		exp /= 2;
	}
	return result;
}

