# pragma once

// GPU Macro Definitions
#if defined(__CUDA_ARCH__)
#define YAMS_IS_GPU_CODE
#endif

#if defined(YAMS_BUILD_GPU_RENDERER) && defined(__CUDACC__)
#ifndef YAMS_NOINLINE
#define YAMS_NOINLINE __attribute__((noinline))
#endif
#define YAMS_CPU_GPU __host__ __device__
#define YAMS_GPU __device__
#if defined(YAMS_IS_GPU_CODE)
#define YAMS_CONST __device__ const
#else
#define YAMS_CONST const
#endif
#else
#define YAMS_CONST const
#define YAMS_CPU_GPU
#define YAMS_GPU
#endif

#ifdef _WIN32 && defined(YAMS_IS_GPU_CODE)
#define YAMS_CPU_GPU_LAMBDA(...) [ =, *this ] YAMS_CPU_GPU(__VA_ARGS__) mutable
#else
#define YAMS_CPU_GPU_LAMBDA(...) [=] YAMS_CPU_GPU(__VA_ARGS__)
#endif