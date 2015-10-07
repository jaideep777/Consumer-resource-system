#ifndef MATH_UTILS_H_
#define MATH_UTILS_H_

#include <fstream>
#include <iostream>
#include <map>
#include <cmath>
#include <cstdlib>
#include <algorithm>
using namespace std;

const float pi = 3.14159265358;

#ifdef __CUDACC__
#define DEVICE_NAME __device__ __host__
#else
#define DEVICE_NAME 
#endif

// =============================================================================
// 		Simple indexing utils
// 		ADDED by : JAIDEEP
// 		6 May 2013
// =============================================================================

inline DEVICE_NAME int ix2(int ix, int iy, int nx){
	return iy*nx + ix;
}

inline DEVICE_NAME int ix3(int ix, int iy, int iz, int nx, int ny){
	return iz*nx*ny + iy*nx + ix;
}

// =============================================================================
// 		Simple math functions
// 		ADDED by : JAIDEEP
// 		18 Nov 2014
// =============================================================================
inline int discretize_index(float f, int n, float minf, float maxf){
	if (fabs(minf - maxf) <1e-6) return 0;
	else return int((f-minf)/(maxf-minf)*(n-1));
}


inline DEVICE_NAME int pos2cell(float x, float dL){
	return int((x-dL/2)/dL+0.5); //+1
}

inline DEVICE_NAME int2 pos2cell(float2 v, float dL){
	int cellx = pos2cell(v.x, dL);
	int celly = pos2cell(v.y, dL);
	return make_int2(cellx, celly);
}


inline DEVICE_NAME float cell2pos(int i, float dL){
	return i*dL + dL/2;
}

inline DEVICE_NAME float2 cell2pos(int2 v, float dL){
	float x = cell2pos(v.x, dL);
	float y = cell2pos(v.y, dL);
	return make_float2(x, y);
}



// =============================================================================
// 		Simple array operations
// 		ADDED by : JAIDEEP
// 		4 Mar 2013
// =============================================================================

template <class T>
T arraySum(T *v, int n){
	T sum=0;
	for (int i=0; i<n; ++i) sum += v[i];
	return sum;
}

template <class T>
T arrayMin(T *v, int n){
	T amin=v[0];
	for (int i=1; i<n; ++i) amin = min(amin, v[i]);
	return amin;
}

template <class T>
T arrayMax(T *v, int n){
	T amax=v[0];
	for (int i=1; i<n; ++i) amax = max(amax, v[i]);
	return amax;
}

// =============================================================================
// 		Random numbers from the default C++ generator
// 		ADDED by : JAIDEEP
// 		6 May 2013
// =============================================================================

//inline float runif_cpp(float min = 0, float max=1){
//	float r = float(rand())/RAND_MAX; 
//	return min + (max-min)*r;
//}

//inline float rnorm_cpp(float mu = 0, float sd = 1){
//	float u = runif_cpp(), v = runif_cpp();		// uniform rn's [0,1] for box-muller
//	float x = sqrt(-2.0*log(u)) * cos(2*pi*v);
//	return mu + sd*x;
//}


// =============================================================================
// 		Map operations
// 		ADDED by : JAIDEEP
// 		21 Dec 2013
// =============================================================================

template <class T>
T mapSum(map <int, T> &m){
	T sum = 0;
	for (typename map<int,T>::iterator it= m.begin(); it != m.end(); ++it){
		sum += it->second;
	}
	return sum;
}

template <class T>
float mapAvg(map <int, T> &m){
	return mapSum(m)/float(m.size());
}



#endif

