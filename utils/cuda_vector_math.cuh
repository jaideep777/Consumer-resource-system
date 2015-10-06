/*
 * Copyright 1993-2012 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 *
 */

/*
    This file implements common mathematical operations on vector types
    (float3, float4 etc.) since these are not provided as standard by CUDA.

    The syntax is modelled on the Cg standard library.

    This is part of the CUTIL library and is not supported by NVIDIA.

    Thanks to Linh Hah for additions and fixes.
*/

#ifndef CUTIL_MATH_H
#define CUTIL_MATH_H

#include <cuda_runtime.h>
#include <cstdlib>

typedef unsigned int uint;
typedef unsigned short ushort;

#ifndef __CUDACC__
#include <math.h>

////////////////////////////////////////////////////////////////////////////////
// host implementations of CUDA functions
////////////////////////////////////////////////////////////////////////////////

inline float fminf(float a, float b){
  return a < b ? a : b;
}

inline float fmaxf(float a, float b){
  return a > b ? a : b;
}

inline int max(int a, int b){
  return a > b ? a : b;
}

inline int min(int a, int b){
  return a < b ? a : b;
}

inline float rsqrtf(float x){
    return 1.0f / sqrtf(x);	// JAIDEEP: prevent division by 0
}
#endif

////////////////////////////////////////////////////////////////////////////////
// constructors
////////////////////////////////////////////////////////////////////////////////

inline __host__ __device__ float2 make_float2(float s){
    return make_float2(s, s);
}
inline __host__ __device__ float2 make_float2(float3 a){
    return make_float2(a.x, a.y);
}
inline __host__ __device__ float2 make_float2(int2 a){
    return make_float2(float(a.x), float(a.y));
}
inline __host__ __device__ float2 make_float2(uint2 a){
    return make_float2(float(a.x), float(a.y));
}

inline __host__ __device__ int2 make_int2(int s){
    return make_int2(s, s);
}
inline __host__ __device__ int2 make_int2(int3 a){
    return make_int2(a.x, a.y);
}
inline __host__ __device__ int2 make_int2(uint2 a){
    return make_int2(int(a.x), int(a.y));
}
inline __host__ __device__ int2 make_int2(float2 a){
    return make_int2(int(a.x), int(a.y));
}

inline __host__ __device__ uint2 make_uint2(uint s){
    return make_uint2(s, s);
}
inline __host__ __device__ uint2 make_uint2(uint3 a){
    return make_uint2(a.x, a.y);
}
inline __host__ __device__ uint2 make_uint2(int2 a){
    return make_uint2(uint(a.x), uint(a.y));
}

inline __host__ __device__ float3 make_float3(float s){
    return make_float3(s, s, s);
}
inline __host__ __device__ float3 make_float3(float2 a){
    return make_float3(a.x, a.y, 0.0f);
}
inline __host__ __device__ float3 make_float3(float2 a, float s){
    return make_float3(a.x, a.y, s);
}
inline __host__ __device__ float3 make_float3(float4 a){
    return make_float3(a.x, a.y, a.z);
}
inline __host__ __device__ float3 make_float3(int3 a){
    return make_float3(float(a.x), float(a.y), float(a.z));
}
inline __host__ __device__ float3 make_float3(uint3 a){
    return make_float3(float(a.x), float(a.y), float(a.z));
}

inline __host__ __device__ int3 make_int3(int s){
    return make_int3(s, s, s);
}
inline __host__ __device__ int3 make_int3(int2 a){
    return make_int3(a.x, a.y, 0);
}
inline __host__ __device__ int3 make_int3(int2 a, int s){
    return make_int3(a.x, a.y, s);
}
inline __host__ __device__ int3 make_int3(uint3 a){
    return make_int3(int(a.x), int(a.y), int(a.z));
}
inline __host__ __device__ int3 make_int3(float3 a){
    return make_int3(int(a.x), int(a.y), int(a.z));
}

inline __host__ __device__ uint3 make_uint3(uint s){
    return make_uint3(s, s, s);
}
inline __host__ __device__ uint3 make_uint3(uint2 a){
    return make_uint3(a.x, a.y, 0);
}
inline __host__ __device__ uint3 make_uint3(uint2 a, uint s){
    return make_uint3(a.x, a.y, s);
}
inline __host__ __device__ uint3 make_uint3(uint4 a){
    return make_uint3(a.x, a.y, a.z);
}
inline __host__ __device__ uint3 make_uint3(int3 a){
    return make_uint3(uint(a.x), uint(a.y), uint(a.z));
}

inline __host__ __device__ float4 make_float4(float s){
    return make_float4(s, s, s, s);
}
inline __host__ __device__ float4 make_float4(float3 a){
    return make_float4(a.x, a.y, a.z, 0.0f);
}
inline __host__ __device__ float4 make_float4(float3 a, float w){
    return make_float4(a.x, a.y, a.z, w);
}
inline __host__ __device__ float4 make_float4(int4 a){
    return make_float4(float(a.x), float(a.y), float(a.z), float(a.w));
}
inline __host__ __device__ float4 make_float4(uint4 a){
    return make_float4(float(a.x), float(a.y), float(a.z), float(a.w));
}

inline __host__ __device__ int4 make_int4(int s){
    return make_int4(s, s, s, s);
}
inline __host__ __device__ int4 make_int4(int3 a){
    return make_int4(a.x, a.y, a.z, 0);
}
inline __host__ __device__ int4 make_int4(int3 a, int w){
    return make_int4(a.x, a.y, a.z, w);
}
inline __host__ __device__ int4 make_int4(uint4 a){
    return make_int4(int(a.x), int(a.y), int(a.z), int(a.w));
}
inline __host__ __device__ int4 make_int4(float4 a){
    return make_int4(int(a.x), int(a.y), int(a.z), int(a.w));
}


inline __host__ __device__ uint4 make_uint4(uint s){
    return make_uint4(s, s, s, s);
}
inline __host__ __device__ uint4 make_uint4(uint3 a){
    return make_uint4(a.x, a.y, a.z, 0);
}
inline __host__ __device__ uint4 make_uint4(uint3 a, uint w){
    return make_uint4(a.x, a.y, a.z, w);
}
inline __host__ __device__ uint4 make_uint4(int4 a){
    return make_uint4(uint(a.x), uint(a.y), uint(a.z), uint(a.w));
}

////////////////////////////////////////////////////////////////////////////////
// negate
////////////////////////////////////////////////////////////////////////////////

inline __host__ __device__ float2 operator-(float2 &a){
    return make_float2(-a.x, -a.y);
}
inline __host__ __device__ int2 operator-(int2 &a){
    return make_int2(-a.x, -a.y);
}
inline __host__ __device__ float3 operator-(float3 &a){
    return make_float3(-a.x, -a.y, -a.z);
}
inline __host__ __device__ int3 operator-(int3 &a){
    return make_int3(-a.x, -a.y, -a.z);
}
inline __host__ __device__ float4 operator-(float4 &a){
    return make_float4(-a.x, -a.y, -a.z, -a.w);
}
inline __host__ __device__ int4 operator-(int4 &a){
    return make_int4(-a.x, -a.y, -a.z, -a.w);
}

////////////////////////////////////////////////////////////////////////////////
// addition
////////////////////////////////////////////////////////////////////////////////

inline __host__ __device__ float2 operator+(float2 a, float2 b){
    return make_float2(a.x + b.x, a.y + b.y);
}
inline __host__ __device__ void operator+=(float2 &a, float2 b){
    a.x += b.x; a.y += b.y;
}
inline __host__ __device__ float2 operator+(float2 a, float b){
    return make_float2(a.x + b, a.y + b);
}
inline __host__ __device__ float2 operator+(float b, float2 a){
    return make_float2(a.x + b, a.y + b);
}
inline __host__ __device__ void operator+=(float2 &a, float b){
    a.x += b; a.y += b;
}

inline __host__ __device__ int2 operator+(int2 a, int2 b){
    return make_int2(a.x + b.x, a.y + b.y);
}
inline __host__ __device__ void operator+=(int2 &a, int2 b){
    a.x += b.x; a.y += b.y;
}
inline __host__ __device__ int2 operator+(int2 a, int b){
    return make_int2(a.x + b, a.y + b);
}
inline __host__ __device__ int2 operator+(int b, int2 a){
    return make_int2(a.x + b, a.y + b);
}
inline __host__ __device__ void operator+=(int2 &a, int b){
    a.x += b; a.y += b;
}

inline __host__ __device__ uint2 operator+(uint2 a, uint2 b){
    return make_uint2(a.x + b.x, a.y + b.y);
}
inline __host__ __device__ void operator+=(uint2 &a, uint2 b){
    a.x += b.x; a.y += b.y;
}
inline __host__ __device__ uint2 operator+(uint2 a, uint b){
    return make_uint2(a.x + b, a.y + b);
}
inline __host__ __device__ uint2 operator+(uint b, uint2 a){
    return make_uint2(a.x + b, a.y + b);
}
inline __host__ __device__ void operator+=(uint2 &a, uint b){
    a.x += b; a.y += b;
}


inline __host__ __device__ float3 operator+(float3 a, float3 b){
    return make_float3(a.x + b.x, a.y + b.y, a.z + b.z);
}
inline __host__ __device__ void operator+=(float3 &a, float3 b){
    a.x += b.x; a.y += b.y; a.z += b.z;
}
inline __host__ __device__ float3 operator+(float3 a, float b){
    return make_float3(a.x + b, a.y + b, a.z + b);
}
inline __host__ __device__ void operator+=(float3 &a, float b){
    a.x += b; a.y += b; a.z += b;
}

inline __host__ __device__ int3 operator+(int3 a, int3 b){
    return make_int3(a.x + b.x, a.y + b.y, a.z + b.z);
}
inline __host__ __device__ void operator+=(int3 &a, int3 b){
    a.x += b.x; a.y += b.y; a.z += b.z;
}
inline __host__ __device__ int3 operator+(int3 a, int b){
    return make_int3(a.x + b, a.y + b, a.z + b);
}
inline __host__ __device__ void operator+=(int3 &a, int b){
    a.x += b; a.y += b; a.z += b;
}

inline __host__ __device__ uint3 operator+(uint3 a, uint3 b){
    return make_uint3(a.x + b.x, a.y + b.y, a.z + b.z);
}
inline __host__ __device__ void operator+=(uint3 &a, uint3 b){
    a.x += b.x; a.y += b.y; a.z += b.z;
}
inline __host__ __device__ uint3 operator+(uint3 a, uint b){
    return make_uint3(a.x + b, a.y + b, a.z + b);
}
inline __host__ __device__ void operator+=(uint3 &a, uint b){
    a.x += b; a.y += b; a.z += b;
}

inline __host__ __device__ int3 operator+(int b, int3 a){
    return make_int3(a.x + b, a.y + b, a.z + b);
}
inline __host__ __device__ uint3 operator+(uint b, uint3 a){
    return make_uint3(a.x + b, a.y + b, a.z + b);
}
inline __host__ __device__ float3 operator+(float b, float3 a){
    return make_float3(a.x + b, a.y + b, a.z + b);
}

inline __host__ __device__ float4 operator+(float4 a, float4 b){
    return make_float4(a.x + b.x, a.y + b.y, a.z + b.z,  a.w + b.w);
}
inline __host__ __device__ void operator+=(float4 &a, float4 b){
    a.x += b.x; a.y += b.y; a.z += b.z; a.w += b.w;
}
inline __host__ __device__ float4 operator+(float4 a, float b){
    return make_float4(a.x + b, a.y + b, a.z + b, a.w + b);
}
inline __host__ __device__ float4 operator+(float b, float4 a){
    return make_float4(a.x + b, a.y + b, a.z + b, a.w + b);
}
inline __host__ __device__ void operator+=(float4 &a, float b){
    a.x += b; a.y += b; a.z += b; a.w += b;
}

inline __host__ __device__ int4 operator+(int4 a, int4 b){
    return make_int4(a.x + b.x, a.y + b.y, a.z + b.z,  a.w + b.w);
}
inline __host__ __device__ void operator+=(int4 &a, int4 b){
    a.x += b.x; a.y += b.y; a.z += b.z; a.w += b.w;
}
inline __host__ __device__ int4 operator+(int4 a, int b){
    return make_int4(a.x + b, a.y + b, a.z + b,  a.w + b);
}
inline __host__ __device__ int4 operator+(int b, int4 a){
    return make_int4(a.x + b, a.y + b, a.z + b,  a.w + b);
}
inline __host__ __device__ void operator+=(int4 &a, int b){
    a.x += b; a.y += b; a.z += b; a.w += b;
}

inline __host__ __device__ uint4 operator+(uint4 a, uint4 b){
    return make_uint4(a.x + b.x, a.y + b.y, a.z + b.z,  a.w + b.w);
}
inline __host__ __device__ void operator+=(uint4 &a, uint4 b){
    a.x += b.x; a.y += b.y; a.z += b.z; a.w += b.w;
}
inline __host__ __device__ uint4 operator+(uint4 a, uint b){
    return make_uint4(a.x + b, a.y + b, a.z + b,  a.w + b);
}
inline __host__ __device__ uint4 operator+(uint b, uint4 a){
    return make_uint4(a.x + b, a.y + b, a.z + b,  a.w + b);
}
inline __host__ __device__ void operator+=(uint4 &a, uint b){
    a.x += b; a.y += b; a.z += b; a.w += b;
}

////////////////////////////////////////////////////////////////////////////////
// subtract
////////////////////////////////////////////////////////////////////////////////

inline __host__ __device__ float2 operator-(float2 a, float2 b){
    return make_float2(a.x - b.x, a.y - b.y);
}
inline __host__ __device__ void operator-=(float2 &a, float2 b){
    a.x -= b.x; a.y -= b.y;
}
inline __host__ __device__ float2 operator-(float2 a, float b){
    return make_float2(a.x - b, a.y - b);
}
inline __host__ __device__ float2 operator-(float b, float2 a){
    return make_float2(b - a.x, b - a.y);
}
inline __host__ __device__ void operator-=(float2 &a, float b){
    a.x -= b; a.y -= b;
}

inline __host__ __device__ int2 operator-(int2 a, int2 b){
    return make_int2(a.x - b.x, a.y - b.y);
}
inline __host__ __device__ void operator-=(int2 &a, int2 b){
    a.x -= b.x; a.y -= b.y;
}
inline __host__ __device__ int2 operator-(int2 a, int b){
    return make_int2(a.x - b, a.y - b);
}
inline __host__ __device__ int2 operator-(int b, int2 a){
    return make_int2(b - a.x, b - a.y);
}
inline __host__ __device__ void operator-=(int2 &a, int b){
    a.x -= b; a.y -= b;
}

inline __host__ __device__ uint2 operator-(uint2 a, uint2 b){
    return make_uint2(a.x - b.x, a.y - b.y);
}
inline __host__ __device__ void operator-=(uint2 &a, uint2 b){
    a.x -= b.x; a.y -= b.y;
}
inline __host__ __device__ uint2 operator-(uint2 a, uint b){
    return make_uint2(a.x - b, a.y - b);
}
inline __host__ __device__ uint2 operator-(uint b, uint2 a){
    return make_uint2(b - a.x, b - a.y);
}
inline __host__ __device__ void operator-=(uint2 &a, uint b){
    a.x -= b; a.y -= b;
}

inline __host__ __device__ float3 operator-(float3 a, float3 b){
    return make_float3(a.x - b.x, a.y - b.y, a.z - b.z);
}
inline __host__ __device__ void operator-=(float3 &a, float3 b){
    a.x -= b.x; a.y -= b.y; a.z -= b.z;
}
inline __host__ __device__ float3 operator-(float3 a, float b){
    return make_float3(a.x - b, a.y - b, a.z - b);
}
inline __host__ __device__ float3 operator-(float b, float3 a){
    return make_float3(b - a.x, b - a.y, b - a.z);
}
inline __host__ __device__ void operator-=(float3 &a, float b){
    a.x -= b; a.y -= b; a.z -= b;
}

inline __host__ __device__ int3 operator-(int3 a, int3 b){
    return make_int3(a.x - b.x, a.y - b.y, a.z - b.z);
}
inline __host__ __device__ void operator-=(int3 &a, int3 b){
    a.x -= b.x; a.y -= b.y; a.z -= b.z;
}
inline __host__ __device__ int3 operator-(int3 a, int b){
    return make_int3(a.x - b, a.y - b, a.z - b);
}
inline __host__ __device__ int3 operator-(int b, int3 a){
    return make_int3(b - a.x, b - a.y, b - a.z);
}
inline __host__ __device__ void operator-=(int3 &a, int b){
    a.x -= b; a.y -= b; a.z -= b;
}

inline __host__ __device__ uint3 operator-(uint3 a, uint3 b){
    return make_uint3(a.x - b.x, a.y - b.y, a.z - b.z);
}
inline __host__ __device__ void operator-=(uint3 &a, uint3 b){
    a.x -= b.x; a.y -= b.y; a.z -= b.z;
}
inline __host__ __device__ uint3 operator-(uint3 a, uint b){
    return make_uint3(a.x - b, a.y - b, a.z - b);
}
inline __host__ __device__ uint3 operator-(uint b, uint3 a){
    return make_uint3(b - a.x, b - a.y, b - a.z);
}
inline __host__ __device__ void operator-=(uint3 &a, uint b){
    a.x -= b; a.y -= b; a.z -= b;
}

inline __host__ __device__ float4 operator-(float4 a, float4 b){
    return make_float4(a.x - b.x, a.y - b.y, a.z - b.z,  a.w - b.w);
}
inline __host__ __device__ void operator-=(float4 &a, float4 b){
    a.x -= b.x; a.y -= b.y; a.z -= b.z; a.w -= b.w;
}
inline __host__ __device__ float4 operator-(float4 a, float b){
    return make_float4(a.x - b, a.y - b, a.z - b,  a.w - b);
}
inline __host__ __device__ void operator-=(float4 &a, float b){
    a.x -= b; a.y -= b; a.z -= b; a.w -= b;
}

inline __host__ __device__ int4 operator-(int4 a, int4 b){
    return make_int4(a.x - b.x, a.y - b.y, a.z - b.z,  a.w - b.w);
}
inline __host__ __device__ void operator-=(int4 &a, int4 b){
    a.x -= b.x; a.y -= b.y; a.z -= b.z; a.w -= b.w;
}
inline __host__ __device__ int4 operator-(int4 a, int b){
    return make_int4(a.x - b, a.y - b, a.z - b,  a.w - b);
}
inline __host__ __device__ int4 operator-(int b, int4 a){
    return make_int4(b - a.x, b - a.y, b - a.z, b - a.w);
}
inline __host__ __device__ void operator-=(int4 &a, int b){
    a.x -= b; a.y -= b; a.z -= b; a.w -= b;
}

inline __host__ __device__ uint4 operator-(uint4 a, uint4 b){
    return make_uint4(a.x - b.x, a.y - b.y, a.z - b.z,  a.w - b.w);
}
inline __host__ __device__ void operator-=(uint4 &a, uint4 b){
    a.x -= b.x; a.y -= b.y; a.z -= b.z; a.w -= b.w;
}
inline __host__ __device__ uint4 operator-(uint4 a, uint b){
    return make_uint4(a.x - b, a.y - b, a.z - b,  a.w - b);
}
inline __host__ __device__ uint4 operator-(uint b, uint4 a){
    return make_uint4(b - a.x, b - a.y, b - a.z, b - a.w);
}
inline __host__ __device__ void operator-=(uint4 &a, uint b){
    a.x -= b; a.y -= b; a.z -= b; a.w -= b;
}

////////////////////////////////////////////////////////////////////////////////
// multiply
////////////////////////////////////////////////////////////////////////////////

inline __host__ __device__ float2 operator*(float2 a, float2 b){
    return make_float2(a.x * b.x, a.y * b.y);
}
inline __host__ __device__ void operator*=(float2 &a, float2 b){
    a.x *= b.x; a.y *= b.y;
}
inline __host__ __device__ float2 operator*(float2 a, float b){
    return make_float2(a.x * b, a.y * b);
}
inline __host__ __device__ float2 operator*(float b, float2 a){
    return make_float2(b * a.x, b * a.y);
}
inline __host__ __device__ void operator*=(float2 &a, float b){
    a.x *= b; a.y *= b;
}

inline __host__ __device__ int2 operator*(int2 a, int2 b){
    return make_int2(a.x * b.x, a.y * b.y);
}
inline __host__ __device__ void operator*=(int2 &a, int2 b){
    a.x *= b.x; a.y *= b.y;
}
inline __host__ __device__ int2 operator*(int2 a, int b){
    return make_int2(a.x * b, a.y * b);
}
inline __host__ __device__ int2 operator*(int b, int2 a){
    return make_int2(b * a.x, b * a.y);
}
inline __host__ __device__ void operator*=(int2 &a, int b){
    a.x *= b; a.y *= b;
}

inline __host__ __device__ uint2 operator*(uint2 a, uint2 b){
    return make_uint2(a.x * b.x, a.y * b.y);
}
inline __host__ __device__ void operator*=(uint2 &a, uint2 b){
    a.x *= b.x; a.y *= b.y;
}
inline __host__ __device__ uint2 operator*(uint2 a, uint b){
    return make_uint2(a.x * b, a.y * b);
}
inline __host__ __device__ uint2 operator*(uint b, uint2 a){
    return make_uint2(b * a.x, b * a.y);
}
inline __host__ __device__ void operator*=(uint2 &a, uint b){
    a.x *= b; a.y *= b;
}

inline __host__ __device__ float3 operator*(float3 a, float3 b){
    return make_float3(a.x * b.x, a.y * b.y, a.z * b.z);
}
inline __host__ __device__ void operator*=(float3 &a, float3 b){
    a.x *= b.x; a.y *= b.y; a.z *= b.z;
}
inline __host__ __device__ float3 operator*(float3 a, float b){
    return make_float3(a.x * b, a.y * b, a.z * b);
}
inline __host__ __device__ float3 operator*(float b, float3 a){
    return make_float3(b * a.x, b * a.y, b * a.z);
}
inline __host__ __device__ void operator*=(float3 &a, float b){
    a.x *= b; a.y *= b; a.z *= b;
}

inline __host__ __device__ int3 operator*(int3 a, int3 b){
    return make_int3(a.x * b.x, a.y * b.y, a.z * b.z);
}
inline __host__ __device__ void operator*=(int3 &a, int3 b){
    a.x *= b.x; a.y *= b.y; a.z *= b.z;
}
inline __host__ __device__ int3 operator*(int3 a, int b){
    return make_int3(a.x * b, a.y * b, a.z * b);
}
inline __host__ __device__ int3 operator*(int b, int3 a){
    return make_int3(b * a.x, b * a.y, b * a.z);
}
inline __host__ __device__ void operator*=(int3 &a, int b){
    a.x *= b; a.y *= b; a.z *= b;
}

inline __host__ __device__ uint3 operator*(uint3 a, uint3 b){
    return make_uint3(a.x * b.x, a.y * b.y, a.z * b.z);
}
inline __host__ __device__ void operator*=(uint3 &a, uint3 b){
    a.x *= b.x; a.y *= b.y; a.z *= b.z;
}
inline __host__ __device__ uint3 operator*(uint3 a, uint b){
    return make_uint3(a.x * b, a.y * b, a.z * b);
}
inline __host__ __device__ uint3 operator*(uint b, uint3 a){
    return make_uint3(b * a.x, b * a.y, b * a.z);
}
inline __host__ __device__ void operator*=(uint3 &a, uint b){
    a.x *= b; a.y *= b; a.z *= b;
}

inline __host__ __device__ float4 operator*(float4 a, float4 b){
    return make_float4(a.x * b.x, a.y * b.y, a.z * b.z,  a.w * b.w);
}
inline __host__ __device__ void operator*=(float4 &a, float4 b){
    a.x *= b.x; a.y *= b.y; a.z *= b.z; a.w *= b.w;
}
inline __host__ __device__ float4 operator*(float4 a, float b){
    return make_float4(a.x * b, a.y * b, a.z * b,  a.w * b);
}
inline __host__ __device__ float4 operator*(float b, float4 a){
    return make_float4(b * a.x, b * a.y, b * a.z, b * a.w);
}
inline __host__ __device__ void operator*=(float4 &a, float b){
    a.x *= b; a.y *= b; a.z *= b; a.w *= b;
}

inline __host__ __device__ int4 operator*(int4 a, int4 b){
    return make_int4(a.x * b.x, a.y * b.y, a.z * b.z,  a.w * b.w);
}
inline __host__ __device__ void operator*=(int4 &a, int4 b){
    a.x *= b.x; a.y *= b.y; a.z *= b.z; a.w *= b.w;
}
inline __host__ __device__ int4 operator*(int4 a, int b){
    return make_int4(a.x * b, a.y * b, a.z * b,  a.w * b);
}
inline __host__ __device__ int4 operator*(int b, int4 a){
    return make_int4(b * a.x, b * a.y, b * a.z, b * a.w);
}
inline __host__ __device__ void operator*=(int4 &a, int b){
    a.x *= b; a.y *= b; a.z *= b; a.w *= b;
}

inline __host__ __device__ uint4 operator*(uint4 a, uint4 b){
    return make_uint4(a.x * b.x, a.y * b.y, a.z * b.z,  a.w * b.w);
}
inline __host__ __device__ void operator*=(uint4 &a, uint4 b){
    a.x *= b.x; a.y *= b.y; a.z *= b.z; a.w *= b.w;
}
inline __host__ __device__ uint4 operator*(uint4 a, uint b){
    return make_uint4(a.x * b, a.y * b, a.z * b,  a.w * b);
}
inline __host__ __device__ uint4 operator*(uint b, uint4 a){
    return make_uint4(b * a.x, b * a.y, b * a.z, b * a.w);
}
inline __host__ __device__ void operator*=(uint4 &a, uint b){
    a.x *= b; a.y *= b; a.z *= b; a.w *= b;
}

////////////////////////////////////////////////////////////////////////////////
// divide
////////////////////////////////////////////////////////////////////////////////

inline __host__ __device__ float2 operator/(float2 a, float2 b){
    return make_float2(a.x / b.x, a.y / b.y);
}
inline __host__ __device__ void operator/=(float2 &a, float2 b){
    a.x /= b.x; a.y /= b.y;
}
inline __host__ __device__ float2 operator/(float2 a, float b){
    return make_float2(a.x / b, a.y / b);
}
inline __host__ __device__ void operator/=(float2 &a, float b){
    a.x /= b; a.y /= b;
}
inline __host__ __device__ float2 operator/(float b, float2 a){
    return make_float2(b / a.x, b / a.y);
}

inline __host__ __device__ float3 operator/(float3 a, float3 b){
    return make_float3(a.x / b.x, a.y / b.y, a.z / b.z);
}
inline __host__ __device__ void operator/=(float3 &a, float3 b){
    a.x /= b.x; a.y /= b.y; a.z /= b.z;
}
inline __host__ __device__ float3 operator/(float3 a, float b){
    return make_float3(a.x / b, a.y / b, a.z / b);
}
inline __host__ __device__ void operator/=(float3 &a, float b){
    a.x /= b; a.y /= b; a.z /= b;
}
inline __host__ __device__ float3 operator/(float b, float3 a){
    return make_float3(b / a.x, b / a.y, b / a.z);
}

inline __host__ __device__ float4 operator/(float4 a, float4 b){
    return make_float4(a.x / b.x, a.y / b.y, a.z / b.z,  a.w / b.w);
}
inline __host__ __device__ void operator/=(float4 &a, float4 b){
    a.x /= b.x; a.y /= b.y; a.z /= b.z; a.w /= b.w;
}
inline __host__ __device__ float4 operator/(float4 a, float b){
    return make_float4(a.x / b, a.y / b, a.z / b,  a.w / b);
}
inline __host__ __device__ void operator/=(float4 &a, float b){
    a.x /= b; a.y /= b; a.z /= b; a.w /= b;
}
inline __host__ __device__ float4 operator/(float b, float4 a){
    return make_float4(b / a.x, b / a.y, b / a.z, b / a.w);
}

////////////////////////////////////////////////////////////////////////////////
// min
////////////////////////////////////////////////////////////////////////////////

inline  __host__ __device__ float2 fminf(float2 a, float2 b){
	return make_float2(fminf(a.x,b.x), fminf(a.y,b.y));
}
inline __host__ __device__ float3 fminf(float3 a, float3 b){
	return make_float3(fminf(a.x,b.x), fminf(a.y,b.y), fminf(a.z,b.z));
}
inline  __host__ __device__ float4 fminf(float4 a, float4 b){
	return make_float4(fminf(a.x,b.x), fminf(a.y,b.y), fminf(a.z,b.z), fminf(a.w,b.w));
}

inline __host__ __device__ int2 min(int2 a, int2 b){
    return make_int2(min(a.x,b.x), min(a.y,b.y));
}
inline __host__ __device__ int3 min(int3 a, int3 b){
    return make_int3(min(a.x,b.x), min(a.y,b.y), min(a.z,b.z));
}
inline __host__ __device__ int4 min(int4 a, int4 b){
    return make_int4(min(a.x,b.x), min(a.y,b.y), min(a.z,b.z), min(a.w,b.w));
}

inline __host__ __device__ uint2 min(uint2 a, uint2 b){
    return make_uint2(min(a.x,b.x), min(a.y,b.y));
}
inline __host__ __device__ uint3 min(uint3 a, uint3 b){
    return make_uint3(min(a.x,b.x), min(a.y,b.y), min(a.z,b.z));
}
inline __host__ __device__ uint4 min(uint4 a, uint4 b){
    return make_uint4(min(a.x,b.x), min(a.y,b.y), min(a.z,b.z), min(a.w,b.w));
}

////////////////////////////////////////////////////////////////////////////////
// max
////////////////////////////////////////////////////////////////////////////////

inline __host__ __device__ float2 fmaxf(float2 a, float2 b){
	return make_float2(fmaxf(a.x,b.x), fmaxf(a.y,b.y));
}
inline __host__ __device__ float3 fmaxf(float3 a, float3 b){
	return make_float3(fmaxf(a.x,b.x), fmaxf(a.y,b.y), fmaxf(a.z,b.z));
}
inline __host__ __device__ float4 fmaxf(float4 a, float4 b){
	return make_float4(fmaxf(a.x,b.x), fmaxf(a.y,b.y), fmaxf(a.z,b.z), fmaxf(a.w,b.w));
}

inline __host__ __device__ int2 max(int2 a, int2 b){
    return make_int2(max(a.x,b.x), max(a.y,b.y));
}
inline __host__ __device__ int3 max(int3 a, int3 b){
    return make_int3(max(a.x,b.x), max(a.y,b.y), max(a.z,b.z));
}
inline __host__ __device__ int4 max(int4 a, int4 b){
    return make_int4(max(a.x,b.x), max(a.y,b.y), max(a.z,b.z), max(a.w,b.w));
}

inline __host__ __device__ uint2 max(uint2 a, uint2 b){
    return make_uint2(max(a.x,b.x), max(a.y,b.y));
}
inline __host__ __device__ uint3 max(uint3 a, uint3 b){
    return make_uint3(max(a.x,b.x), max(a.y,b.y), max(a.z,b.z));
}
inline __host__ __device__ uint4 max(uint4 a, uint4 b){
    return make_uint4(max(a.x,b.x), max(a.y,b.y), max(a.z,b.z), max(a.w,b.w));
}

////////////////////////////////////////////////////////////////////////////////
// lerp
// - linear interpolation between a and b, based on value t in [0, 1] range
////////////////////////////////////////////////////////////////////////////////

inline __device__ __host__ float lerp(float a, float b, float t){
    return a + t*(b-a);
}
inline __device__ __host__ float2 lerp(float2 a, float2 b, float t){
    return a + t*(b-a);
}
inline __device__ __host__ float3 lerp(float3 a, float3 b, float t){
    return a + t*(b-a);
}
inline __device__ __host__ float4 lerp(float4 a, float4 b, float t){
    return a + t*(b-a);
}

////////////////////////////////////////////////////////////////////////////////
// clamp
// - clamp the value v to be in the range [a, b]
////////////////////////////////////////////////////////////////////////////////

inline __device__ __host__ float clamp(float f, float a, float b){
    return fmaxf(a, fminf(f, b));
}
inline __device__ __host__ int clamp(int f, int a, int b){
    return max(a, min(f, b));
}
inline __device__ __host__ uint clamp(uint f, uint a, uint b){
    return max(a, min(f, b));
}

inline __device__ __host__ float2 clamp(float2 v, float a, float b){
    return make_float2(clamp(v.x, a, b), clamp(v.y, a, b));
}
inline __device__ __host__ float2 clamp(float2 v, float2 a, float2 b){
    return make_float2(clamp(v.x, a.x, b.x), clamp(v.y, a.y, b.y));
}
inline __device__ __host__ float3 clamp(float3 v, float a, float b){
    return make_float3(clamp(v.x, a, b), clamp(v.y, a, b), clamp(v.z, a, b));
}
inline __device__ __host__ float3 clamp(float3 v, float3 a, float3 b){
    return make_float3(clamp(v.x, a.x, b.x), clamp(v.y, a.y, b.y), clamp(v.z, a.z, b.z));
}
inline __device__ __host__ float4 clamp(float4 v, float a, float b){
    return make_float4(clamp(v.x, a, b), clamp(v.y, a, b), clamp(v.z, a, b), clamp(v.w, a, b));
}
inline __device__ __host__ float4 clamp(float4 v, float4 a, float4 b){
    return make_float4(clamp(v.x, a.x, b.x), clamp(v.y, a.y, b.y), clamp(v.z, a.z, b.z), clamp(v.w, a.w, b.w));
}

inline __device__ __host__ int2 clamp(int2 v, int a, int b){
    return make_int2(clamp(v.x, a, b), clamp(v.y, a, b));
}
inline __device__ __host__ int2 clamp(int2 v, int2 a, int2 b){
    return make_int2(clamp(v.x, a.x, b.x), clamp(v.y, a.y, b.y));
}
inline __device__ __host__ int3 clamp(int3 v, int a, int b){
    return make_int3(clamp(v.x, a, b), clamp(v.y, a, b), clamp(v.z, a, b));
}
inline __device__ __host__ int3 clamp(int3 v, int3 a, int3 b){
    return make_int3(clamp(v.x, a.x, b.x), clamp(v.y, a.y, b.y), clamp(v.z, a.z, b.z));
}
inline __device__ __host__ int4 clamp(int4 v, int a, int b){
    return make_int4(clamp(v.x, a, b), clamp(v.y, a, b), clamp(v.z, a, b), clamp(v.w, a, b));
}
inline __device__ __host__ int4 clamp(int4 v, int4 a, int4 b){
    return make_int4(clamp(v.x, a.x, b.x), clamp(v.y, a.y, b.y), clamp(v.z, a.z, b.z), clamp(v.w, a.w, b.w));
}

inline __device__ __host__ uint2 clamp(uint2 v, uint a, uint b){
    return make_uint2(clamp(v.x, a, b), clamp(v.y, a, b));
}
inline __device__ __host__ uint2 clamp(uint2 v, uint2 a, uint2 b){
    return make_uint2(clamp(v.x, a.x, b.x), clamp(v.y, a.y, b.y));
}
inline __device__ __host__ uint3 clamp(uint3 v, uint a, uint b){
    return make_uint3(clamp(v.x, a, b), clamp(v.y, a, b), clamp(v.z, a, b));
}
inline __device__ __host__ uint3 clamp(uint3 v, uint3 a, uint3 b){
    return make_uint3(clamp(v.x, a.x, b.x), clamp(v.y, a.y, b.y), clamp(v.z, a.z, b.z));
}
inline __device__ __host__ uint4 clamp(uint4 v, uint a, uint b){
    return make_uint4(clamp(v.x, a, b), clamp(v.y, a, b), clamp(v.z, a, b), clamp(v.w, a, b));
}
inline __device__ __host__ uint4 clamp(uint4 v, uint4 a, uint4 b){
    return make_uint4(clamp(v.x, a.x, b.x), clamp(v.y, a.y, b.y), clamp(v.z, a.z, b.z), clamp(v.w, a.w, b.w));
}

////////////////////////////////////////////////////////////////////////////////
// dot product
////////////////////////////////////////////////////////////////////////////////

inline __host__ __device__ float dot(float2 a, float2 b){ 
    return a.x * b.x + a.y * b.y;
}
inline __host__ __device__ float dot(float3 a, float3 b){ 
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
inline __host__ __device__ float dot(float4 a, float4 b){ 
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}

inline __host__ __device__ int dot(int2 a, int2 b){ 
    return a.x * b.x + a.y * b.y;
}
inline __host__ __device__ int dot(int3 a, int3 b){ 
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
inline __host__ __device__ int dot(int4 a, int4 b){ 
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}

inline __host__ __device__ uint dot(uint2 a, uint2 b){ 
    return a.x * b.x + a.y * b.y;
}
inline __host__ __device__ uint dot(uint3 a, uint3 b){ 
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
inline __host__ __device__ uint dot(uint4 a, uint4 b){ 
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}

////////////////////////////////////////////////////////////////////////////////
// length
////////////////////////////////////////////////////////////////////////////////

inline __host__ __device__ float length(float2 v){
    return sqrtf(dot(v, v));
}
inline __host__ __device__ float length(float3 v){
    return sqrtf(dot(v, v));
}
inline __host__ __device__ float length(float4 v){
    return sqrtf(dot(v, v));
}

////////////////////////////////////////////////////////////////////////////////
// normalize
////////////////////////////////////////////////////////////////////////////////

inline __host__ __device__ float2 normalize(float2 v){
    float invLen = rsqrtf(dot(v, v));
    return v * invLen;
}
inline __host__ __device__ float3 normalize(float3 v){
    float invLen = rsqrtf(dot(v, v));
    return v * invLen;
}
inline __host__ __device__ float4 normalize(float4 v){
    float invLen = rsqrtf(dot(v, v));
    return v * invLen;
}

////////////////////////////////////////////////////////////////////////////////
// floor
////////////////////////////////////////////////////////////////////////////////

inline __host__ __device__ float2 floorf(float2 v){
    return make_float2(floorf(v.x), floorf(v.y));
}
inline __host__ __device__ float3 floorf(float3 v){
    return make_float3(floorf(v.x), floorf(v.y), floorf(v.z));
}
inline __host__ __device__ float4 floorf(float4 v){
    return make_float4(floorf(v.x), floorf(v.y), floorf(v.z), floorf(v.w));
}

////////////////////////////////////////////////////////////////////////////////
// frac - returns the fractional portion of a scalar or each vector component
////////////////////////////////////////////////////////////////////////////////

inline __host__ __device__ float fracf(float v){
    return v - floorf(v);
}
inline __host__ __device__ float2 fracf(float2 v){
    return make_float2(fracf(v.x), fracf(v.y));
}
inline __host__ __device__ float3 fracf(float3 v){
    return make_float3(fracf(v.x), fracf(v.y), fracf(v.z));
}
inline __host__ __device__ float4 fracf(float4 v){
    return make_float4(fracf(v.x), fracf(v.y), fracf(v.z), fracf(v.w));
}

////////////////////////////////////////////////////////////////////////////////
// fmod
////////////////////////////////////////////////////////////////////////////////

inline __host__ __device__ float2 fmodf(float2 a, float2 b){
    return make_float2(fmodf(a.x, b.x), fmodf(a.y, b.y));
}
inline __host__ __device__ float3 fmodf(float3 a, float3 b){
    return make_float3(fmodf(a.x, b.x), fmodf(a.y, b.y), fmodf(a.z, b.z));
}
inline __host__ __device__ float4 fmodf(float4 a, float4 b){
    return make_float4(fmodf(a.x, b.x), fmodf(a.y, b.y), fmodf(a.z, b.z), fmodf(a.w, b.w));
}

////////////////////////////////////////////////////////////////////////////////
// absolute value
////////////////////////////////////////////////////////////////////////////////

inline __host__ __device__ float2 fabs(float2 v){
	return make_float2(fabs(v.x), fabs(v.y));
}
inline __host__ __device__ float3 fabs(float3 v){
	return make_float3(fabs(v.x), fabs(v.y), fabs(v.z));
}
inline __host__ __device__ float4 fabs(float4 v){
	return make_float4(fabs(v.x), fabs(v.y), fabs(v.z), fabs(v.w));
}

inline __host__ __device__ int2 abs(int2 v){
	return make_int2(abs(v.x), abs(v.y));
}
inline __host__ __device__ int3 abs(int3 v){
	return make_int3(abs(v.x), abs(v.y), abs(v.z));
}
inline __host__ __device__ int4 abs(int4 v){
	return make_int4(abs(v.x), abs(v.y), abs(v.z), abs(v.w));
}

////////////////////////////////////////////////////////////////////////////////
// reflect
// - returns reflection of incident ray I around surface normal N
// - N should be normalized, reflected vector's length is equal to length of I
////////////////////////////////////////////////////////////////////////////////

inline __host__ __device__ float3 reflect(float3 i, float3 n){
	return i - 2.0f * n * dot(n,i);
}

////////////////////////////////////////////////////////////////////////////////
// cross product
////////////////////////////////////////////////////////////////////////////////

inline __host__ __device__ float3 cross(float3 a, float3 b){ 
    return make_float3(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x); 
}

////////////////////////////////////////////////////////////////////////////////
// smoothstep
// - returns 0 if x < a
// - returns 1 if x > b
// - otherwise returns smooth interpolation between 0 and 1 based on x
////////////////////////////////////////////////////////////////////////////////

inline __device__ __host__ float smoothstep(float a, float b, float x){
	float y = clamp((x - a) / (b - a), 0.0f, 1.0f);
	return (y*y*(3.0f - (2.0f*y)));
}
inline __device__ __host__ float2 smoothstep(float2 a, float2 b, float2 x){
	float2 y = clamp((x - a) / (b - a), 0.0f, 1.0f);
	return (y*y*(make_float2(3.0f) - (make_float2(2.0f)*y)));
}
inline __device__ __host__ float3 smoothstep(float3 a, float3 b, float3 x){
	float3 y = clamp((x - a) / (b - a), 0.0f, 1.0f);
	return (y*y*(make_float3(3.0f) - (make_float3(2.0f)*y)));
}
inline __device__ __host__ float4 smoothstep(float4 a, float4 b, float4 x){
	float4 y = clamp((x - a) / (b - a), 0.0f, 1.0f);
	return (y*y*(make_float4(3.0f) - (make_float4(2.0f)*y)));
}










// =============================================================================
// 		Some useful 2D vector utilities
// 		ADDED by : JAIDEEP
// 		19 Feb 2013
// =============================================================================
#include "simple_math.h"
#include <curand.h>
#include <cuda_runtime.h>
#include <thrust/fill.h>

// =============================================================================
// 		Periodic Distance related utilities
// 		ADDED by : JAIDEEP
// 		19 Feb 2013
// =============================================================================

// indicator variable I{x<a} = {1, x < a; 0 otherwise}
inline __device__ __host__ float indicator(float x, float a){
	return float(x < a);
}

// brings x to principal range defined by (xmin, xmax)
inline __device__ __host__ void makePeriodic(float &x, float xmin, float xmax){
	x = xmin + fmodf(x-xmin, xmax-xmin);
	x += xmax-xmin;							// this second fmodf is required because some implementations..
	x = xmin + fmodf(x-xmin, xmax-xmin);	// ..of fmodf use negative remainders in principle range 
//	x = x + (float(x < xmin) - float(x > xmax))*(xmax-xmin);
}

// brings x to principal range defined by (xmin, xmax)
inline __device__ __host__ int makePeriodicID(int ix, int n){
	ix = ix % n;
	ix += n;							// this second fmodf is required because some implementations..
	ix = ix % n;	// ..of fmodf use negative remainders in principle range 
	return ix;
//	x = x + (float(x < xmin) - float(x > xmax))*(xmax-xmin);
}

// calculate displacement x1 ---> x2 on a cylinder with circumference X
// limits (-X/2 , X/2) are inspired from principle range of theta
inline __device__ __host__ float periodicDisplacement(float x1, float x2, float X){
	float dx = (x2-x1);
	makePeriodic(dx, -X/2, X/2);
	return dx;
}

// periodic distance is minimum distance on a cylinder
inline __device__ __host__ float periodicDistance(float x1, float x2, float X){
	return fabs(periodicDisplacement(x1,x2,X));
}

// calculate displacement r1 ---> r2 on a hyper-cylinder with circumference X, Y
inline __device__ __host__ float2 periodicDisplacement(float2 r1, float2 r2, float X, float Y){
	return make_float2(  periodicDisplacement(r1.x, r2.x, X),
						 periodicDisplacement(r1.y, r2.y, Y) );
}


// =============================================================================
// 		Random numbers from CURAND generator (on host)
// 		ADDED by : JAIDEEP
// 		6 May 2013
// =============================================================================

// Implementations that use the specified CUDA generator
// --- 1-D ---
inline __host__ float runif(curandGenerator_t rgen, float rmin=0, float rmax=1){
	static float * rands = new float;	// static ensures memory allocation only once
	curandGenerateUniform(rgen, (float*)rands, 1);
	return rmin + (*rands)*(rmax-rmin);
}

inline __host__  float rnorm(curandGenerator_t rgen, float mu=0, float sd=1){
	static float * rands = new float;
	curandGenerateNormal(rgen, (float*)rands, 1, mu, sd);
	return *rands;
}

// --- 2-D ---
inline __host__  float2 rnorm2(curandGenerator_t rgen, float mu=0, float sd=1){
	static float * rands = new float[2];
	curandGenerateNormal(rgen, (float*)rands, 2, mu, sd);
	return make_float2(rands[0], rands[1]);
}

inline __host__  float2 runif2(curandGenerator_t rgen, float _norm){
	float theta = runif(rgen, 0, 2*pi);
	float x = _norm*cos(theta);
	float y = _norm*sin(theta); 
	return make_float2(x,y);
}

inline __host__  float2 runif2(curandGenerator_t rgen, float xlim, float ylim){
	float x = runif(rgen, -xlim, xlim);
	float y = runif(rgen, -ylim, ylim);
	return make_float2(x,y);
}

inline __host__  float2 runif2(curandGenerator_t rgen, float xmin, float xmax, float ymin, float ymax){
	float x = runif(rgen, xmin, xmax);
	float y = runif(rgen, ymin, ymax);
	return make_float2(x,y);
}

inline __host__  float2 rnorm2_bounded(curandGenerator_t rgen, float mu, float sd, float rmin, float rmax){
	float2 r = rnorm2(rgen, mu, sd);
	r.x = clamp(r.x, rmin, rmax);
	r.y = clamp(r.y, rmin, rmax);
	return r;
}

// Implementations that use the default C++ generator
// Added on: 9 Nov 2014 ... :/ sigh!
// --- 1-D ---
inline __host__ float runif(float rmin=0, float rmax=1){
	float r = float(rand())/RAND_MAX; 
	return rmin + (rmax-rmin)*r;
}

inline __host__  float rnorm(float mu=0, float sd=1){
	float u = runif(), v = runif();		// uniform rn's [0,1] for box-muller
	float x = sqrt(-2.0*log(u)) * cos(2*pi*v);
	return mu + sd*x;
}

// --- 2-D ---
inline __host__  float2 rnorm2(float mu=0, float sd=1){
	float u = runif(), v = runif();		// uniform rn's [0,1] for box-muller
	float r = sqrt(-2.0*log(u));
	float x = r * cos(2*pi*v);
	float y = r * sin(2*pi*v);
	return mu + sd*make_float2(x,y);
}

inline __host__  float2 runif2(float _norm){
	float theta = runif(0.f, 2*pi);
	float x = _norm*cos(theta);
	float y = _norm*sin(theta); 
	return make_float2(x,y);
}

inline __host__  float2 runif2(float xlim, float ylim){
	float x = runif(-xlim, xlim);
	float y = runif(-ylim, ylim);
	return make_float2(x,y);
}

inline __host__  float2 runif2(float xmin, float xmax, float ymin, float ymax){
	float x = runif(xmin, xmax);
	float y = runif(ymin, ymax);
	return make_float2(x,y);
}

inline __host__  float2 rnorm2_bounded(float mu, float sd, float rmin, float rmax){
	float2 r = rnorm2(mu, sd);
	r.x = clamp(r.x, rmin, rmax);
	r.y = clamp(r.y, rmin, rmax);
	return r;
}

inline __host__ float2 rotateRandom(float2 v, float theta_mu, float theta_sd, float theta_max){
	float theta = rnorm(theta_mu, theta_sd);
	theta = clamp(theta, -theta_max, theta_max);
	float cosL = cos(theta);
	float sinL = sin(theta);
	return make_float2(v.x*cosL - v.y*sinL, v.x*sinL + v.y*cosL);	
}

// =============================================================================
// 		Vector length / normalization related functions
// 		ADDED by : JAIDEEP
// 		8 May 2013
// =============================================================================

// this functions will not give divide by zero error. But if v is zero, the 
// normalized vector is meaningless. User should deal with this case.
inline __host__ __device__ float2 normalizeSafe(float2 v, float tol=1e-12){
	v = v+tol;
    float invLen = rsqrtf(dot(v, v));
    return v * invLen;
}

// this functions will not give divide by zero error. But if v close to zero, 
// zero vector is returned instead of noramlized vector
inline __host__ __device__ float2 normalizeSafeZero(float2 v, float l2Tol = 1e-12, float tol=1e-12){
	v = v+tol;
	float l2 = dot(v, v);
    float invLen = rsqrtf(l2);
    return (v * invLen)*float(l2 > l2Tol);
}


// =============================================================================
// 		Some simple wrappers for thrust functions
// =============================================================================

//template <class T>
//inline __device__ void thrust_setZero(T * v, int n){
//	thrust::fill((thrust::device_ptr <T>)v, (thrust::device_ptr <T>)v+n, (T)0);
//}

#endif



