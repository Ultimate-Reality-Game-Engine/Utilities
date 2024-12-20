#ifndef ULTREALITY_MATH_SSE2_CONFIG_H
#define ULTREALITY_MATH_SSE2_CONFIG_H

#ifndef __cplusplus
#error This library requires C++
#endif

// Library version
#define VECTOR_MATH_LIB_VERSION 100 // 1.0.0

#include <SIMDLevel.h>

#if !defined(_NO_VECTOR_OVERLOADS_) && (defined(__clang__) || defined(__GNUC__)) && !defined(_NO_INTRINSICS_)
#define _NO_VEC_OVERLOADS_
#endif

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable:4514 4820)
// C4514/4820: Off by default noise
#endif
#include <math.h>
#include <float.h>
#ifdef _MSC_VER
#pragma warning(pop)
#endif

#if !defined(_NO_INTRINSICS_)


#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4987)
// C4987: Off by default noise
#endif
#if defined(_MSC_VER) || defined(__MINGW32__)
#include <intrin.h>
#endif
#ifdef _MSC_VER
#pragma warning(pop)
#endif

#if (defined(__clang__) || defined(__GNUC__)) && (__x86_64__ || __i386__) && !defined(__MINGW32__)
#include <cpuid.h>
#endif

#if defined(_SSE2_INTRINSICS_)
#include <xmmintrin.h>
#include <emmintrin.h>

#if defined(_SSE3_INTRINSICS_)
#include <pmmintrin.h>
#endif

#if defined(_SSE4_INTRINSICS_)
#include <smmintrin.h>
#endif

#if defined(_AVX_INTRINSICS_)
#include <immintrin.h>
#endif
#endif // _SSE2_INTRINSICS_

#endif // !_NO_INTRINSICS_

#if (__cplusplus >= 201703L)
#define ALIGNED(x) alignas(x)
#define ALIGNED_STRUCT(x) struct alignas(x)
#elif defined(__GNUC__)
#define ALIGNED(x) __attribute__ ((aligned(x)))
#define ALIGNED_STRUCT(x) struct __attribute__ ((aligned(x)))
#else
#define ALIGNED(x) __declspec(align(x))
#define ALIGNED_STRUCT(x) __declspec(align(x)) struct
#endif

#if (__cplusplus >= 202002L)
#include <compare>
#endif

#if defined(_SSE2_INTRINSICS_) && !defined(_NO_INTRINSICS_)

#if defined(_NO_MOVNT_)
#define STREAM_PS(p, a) _mm_store_ps((p), (a))
#define _256_STREAM_PS(p, a) _mm256_store_ps((p), (a))
#define SFENCE()
#else
#define STREAM_PS(p, a) _mm_stream_ps((p), (a))
#define _256_STREAM_PS(p, a) _mm256_stream_ps((p), (a))
#define SFENCE() _mm_sfence()
#endif

#if defined(_FMA3_INTRINSICS_)
#define FMADD_PS(a, b, c) _mm_fmadd_ps((a), (b), (c))
#define FNMADD_PS(a, b, c) _mm_fnmadd_ps((a), (b), (c))
#else
#define FMADD_PS(a, b, c) _mm_add_ps(_mm_mul_ps((a), (b)), (c))
#define FNMADD_PS(a, b, c) _mm_sub_ps((c), _mm_mul_ps((a), (b)))
#endif

#if defined(_AVX_INTRINSICS_) && defined(_FAVOR_INTEL_)
#define PERMUTE_PS(v, c) _mm_permute_ps((v), c)
#else
#define PERMUTE_PS(v, c) _mm_shuffle_ps((v), (v), c)
#endif

#if defined(__GNUC__) && !defined(__clang__) && (__GNUC__ < 11)
#define LOADU_SI16( p ) _mm_cvtsi32_si128(*reinterpret_cast<unsigned short const*>(p))
#else
#define LOADU_SI16( p ) _mm_loadu_si16(p)
#endif

#endif // _SSE2_INTRINSICS_ && !_NO_INTRINSICS_

#include <sal.h>
#if defined(DEBUG) || defined(_DEBUG)
#include <assert.h>
#endif

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4005 4668)
// C4005/4668: Old header issue
#endif
#include <stdint.h>
#ifdef _MSC_VER
#pragma warning(pop)
#endif

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable:4068 4201 4365 4324 4820)
	 // C4068: ignore unknown pragmas
	 // C4201: nonstandard extension used : nameless struct/union
	 // C4365: Off by default noise
	 // C4324/4820: padding warnings
#endif

#ifdef _PREFAST_
#pragma prefast(push)
#pragma prefast(disable : 25000, "A_VECTOR is 16 bytes")
#endif

#include <VectorMathConstants.h>

namespace UltReality::Math
{
#if defined(_SSE2_INTRINSICS_) && !defined(_NO_INTRINSICS_)
	using VECTOR = __m128;
#else
	struct __vector4
	{
		union
		{
			float vector4_f32[4];
			uint32_t vector4_u32[4];
		};
	};

	using VECTOR = __vector4;
#endif

#if (defined(_M_IX86) || _VECTORCALL_ || __i386__) && !defined(_NO_INTRINSICS_)
    // Define alias to be used for 1st-3rd vector type arguments. Passed in register for x86, and vector call convention; by reference otherwise
    typedef const VECTOR A_VECTOR;
#else
    // Define alias to be used for 1st-3rd vector type arguments. Passed in register for x86, and vector call convention; by reference otherwise
	typedef const VECTOR& A_VECTOR;
#endif

#if _VECTORCALL_ && !defined(_NO_INTRINSICS_)
    // Define alias to be used for 4th vector type argument. Passed in register for vector call convention; by reference otherwise
	typedef const VECTOR B_VECTOR;
#else
    // Define alias to be used for 4th vector type argument. Passed in register for vector call convention; by reference otherwise
	typedef const VECTOR& B_VECTOR;
#endif

#if _VECTORCALL_ && !defined(_NO_INTRINSICS_)
    // Define alias to be used for 5th & 6th vector type arguments. Passes in register for vector call convention; by reference otherwise
	typedef VECTOR C_VECTOR;
#else
    // Define alias to be used for 5th & 6th vector type arguments. Passes in register for vector call convention; by reference otherwise
	typedef VECTOR& C_VECTOR;
#endif

	// Define alias to be used for 7th+ vector type arguments. Passed by reference
	typedef VECTOR& D_VECTOR;

	ALIGNED_STRUCT(16) VECTOR_F32
	{
		union
		{
			float f[4];
			VECTOR v;
		};

		operator VECTOR() const noexcept;

		operator const float*() const noexcept;

#if defined(_SSE2_INTRINSICS_)
		operator __m128i() const noexcept;

		operator __m128d() const noexcept;
#endif
	};

	ALIGNED_STRUCT(16) VECTOR_I32
	{
		union
		{
			int32_t i[4];
			VECTOR v;
		};

		operator VECTOR() const noexcept;

#if defined(_SSE2_INTRINSICS_)
		operator __m128i() const noexcept;

		operator __m128d() const noexcept;
#endif
	};

	ALIGNED_STRUCT(16) VECTOR_U32
	{
		union
		{
			uint32_t u[4];
			VECTOR v;
		};

		operator VECTOR() const noexcept;

#if defined(_SSE2_INTRINSICS_)
		operator __m128i() const noexcept;

		operator __m128d() const noexcept;
#endif
	};

	ALIGNED_STRUCT(16) VECTOR_U8
	{
		union
		{
			uint8_t u[16];
			VECTOR v;
		};

		operator VECTOR() const noexcept;

#if defined(_SSE2_INTRINSICS_)
		operator __m128i() const noexcept;

		operator __m128d() const noexcept;
#endif
	};

#if !defined(_NO_VECTOR_OVERLOADS_)
	VECTOR VEC_CALLCONV operator+ (A_VECTOR v) noexcept;

	VECTOR VEC_CALLCONV operator- (A_VECTOR v) noexcept;

	VECTOR& VEC_CALLCONV operator+= (VECTOR& v1, A_VECTOR v2) noexcept;

	VECTOR& VEC_CALLCONV operator-= (VECTOR& v1, A_VECTOR v2) noexcept;

	VECTOR& VEC_CALLCONV operator*= (VECTOR& v1, A_VECTOR v2) noexcept;

	VECTOR& VEC_CALLCONV operator/= (VECTOR& v1, A_VECTOR v2) noexcept;

	VECTOR& operator*= (VECTOR& v, const float s) noexcept;

	VECTOR& operator/= (VECTOR& v, const float s) noexcept;

	VECTOR VEC_CALLCONV operator+ (A_VECTOR v1, A_VECTOR v2) noexcept;

	VECTOR VEC_CALLCONV operator- (A_VECTOR v1, A_VECTOR v2) noexcept;

	VECTOR VEC_CALLCONV operator* (A_VECTOR v1, A_VECTOR v2) noexcept;

	VECTOR VEC_CALLCONV operator/ (A_VECTOR v1, A_VECTOR v2) noexcept;

	VECTOR VEC_CALLCONV operator* (A_VECTOR v, const float s) noexcept;

	VECTOR VEC_CALLCONV operator/ (A_VECTOR v, const float s) noexcept;

	VECTOR VEC_CALLCONV operator* (const float s, A_VECTOR v) noexcept;
#endif // !_NO_VECTOR_OVERLOADS

	namespace VEC
	{
		VECTOR VEC_CALLCONV ConvertIntToFloat(A_VECTOR vInt, uint32_t divExponent) noexcept;
		VECTOR VEC_CALLCONV ConvertFloatToInt(A_VECTOR vFloat, uint32_t mulExponent) noexcept;
		VECTOR VEC_CALLCONV ConvertUIntToFloat(A_VECTOR vUInt, uint32_t divExponent) noexcept;
		VECTOR VEC_CALLCONV ConvertFloatToUInt(A_VECTOR vFloat, uint32_t mulExponent) noexcept;

        VECTOR VEC_CALLCONV SetBinaryConstant(uint32_t C0, uint32_t C1, uint32_t C2, uint32_t C3) noexcept;
        VECTOR VEC_CALLCONV SplatConstant(int32_t intConstant, uint32_t divExponent) noexcept;
        VECTOR VEC_CALLCONV SplatConstantInt(int32_t intConstant) noexcept;

        VECTOR VEC_CALLCONV LoadInt(_In_ const uint32_t* pSource) noexcept;
        VECTOR VEC_CALLCONV LoadFloat(_In_ const float* pSource) noexcept;

        VECTOR VEC_CALLCONV LoadInt2(_In_reads_(2) const uint32_t* pSource) noexcept;
        VECTOR VEC_CALLCONV LoadAInt2(_In_reads_(2) const uint32_t* pSource) noexcept;

        VECTOR VEC_CALLCONV LoadInt3(_In_reads_(3) const uint32_t* pSource) noexcept;
        VECTOR VEC_CALLCONV LoadAInt3(_In_reads_(3) const uint32_t* pSource) noexcept;

        VECTOR VEC_CALLCONV LoadInt4(_In_reads_(4) const uint32_t* pSource) noexcept;
        VECTOR VEC_CALLCONV LoadAInt4(_In_reads_(4) const uint32_t* pSource) noexcept;

        void VEC_CALLCONV StoreInt(_Out_ uint32_t* pDestination, _In_ A_VECTOR v) noexcept;
        void VEC_CALLCONV StoreFloat(_Out_ float* pDestination, _In_ A_VECTOR v) noexcept;

        void VEC_CALLCONV StoreInt2(_Out_writes_(2) uint32_t* pDestination, _In_ A_VECTOR v) noexcept;
        void VEC_CALLCONV StoreAInt2(_Out_writes_(2) uint32_t* pDestination, _In_ A_VECTOR v) noexcept;

        void VEC_CALLCONV StoreInt3(_Out_writes_(3) uint32_t* pDestination, _In_ A_VECTOR v) noexcept;
        void VEC_CALLCONV StoreAInt3(_Out_writes_(3) uint32_t* pDestination, _In_ A_VECTOR v) noexcept;

        void VEC_CALLCONV StoreInt4(_Out_writes_(4) uint32_t* pDestination, _In_ A_VECTOR v) noexcept;
        void VEC_CALLCONV StoreAInt4(_Out_writes_(4) uint32_t* pDestination, _In_ A_VECTOR v) noexcept;

        /// <summary>
        /// Returns a vector of 0.0f, 0.0f, 0.0f, 0.0f
        /// </summary>
        /// <returns></returns>
        VECTOR VEC_CALLCONV Zero() noexcept;
        /// <summary>
        /// Returns a vector initialized with the four float values passed
        /// </summary>
        /// <param name="x">X component</param>
        /// <param name="y">Y component</param>
        /// <param name="z">Z component</param>
        /// <param name="w">W component</param>
        /// <returns></returns>
        VECTOR VEC_CALLCONV Set(float x, float y, float z, float w) noexcept;
        /// <summary>
        /// Returns a vector initialized with the four integer values passed
        /// </summary>
        /// <param name="x">X component</param>
        /// <param name="y">Y component</param>
        /// <param name="z">Z component</param>
        /// <param name="w">W component</param>
        /// <returns></returns>
        VECTOR VEC_CALLCONV SetInt(uint32_t x, uint32_t y, uint32_t z, uint32_t w) noexcept;
        /// <summary>
        /// Returns a vector with all elements replicating the float value passed
        /// </summary>
        /// <param name="value">Float to replicate</param>
        /// <returns></returns>
        VECTOR VEC_CALLCONV Replicate(float value) noexcept;
        /// <summary>
        /// Returns a vector with all elements replicating the float passed by pointer
        /// </summary>
        /// <param name="pValue">Pointer to the float to replicate</param>
        /// <returns></returns>
        VECTOR VEC_CALLCONV ReplicatePtr(_In_ const float* pValue) noexcept;
        /// <summary>
        /// Returns a vector with all elements replicating the value passed
        /// </summary>
        /// <param name="value">Unsigned int to replicate</param>
        /// <returns></returns>
        VECTOR VEC_CALLCONV ReplicateInt(uint32_t value) noexcept;
        /// <summary>
        /// Returns a vector with all elements replicating the value passed by pointer
        /// </summary>
        /// <param name="pValue">Pointer to unsinged int to replicate</param>
        /// <returns></returns>
        VECTOR VEC_CALLCONV ReplicateIntPtr(_In_ const uint32_t pValue) noexcept;
        /// <summary>
        /// Returns a vector with bits true (true mask)
        /// </summary>
        /// <returns></returns>
        VECTOR VEC_CALLCONV TrueInt() noexcept;
        /// <summary>
        /// Returns a vector with all bits clear (false mask)
        /// </summary>
        /// <returns></returns>
        VECTOR VEC_CALLCONV FalseInt() noexcept;
        /// <summary>
        /// Returns a vector with all elements replicating the x component of the provided vector
        /// </summary>
        /// <param name="v">Vector to copy the x component from</param>
        /// <returns></returns>
        VECTOR VEC_CALLCONV SplatX(A_VECTOR v) noexcept;
        /// <summary>
        /// Returns a vector with all elements replicating the y component of the provided vector
        /// </summary>
        /// <param name="v">Vector to copy the y component from</param>
        /// <returns></returns>
        VECTOR VEC_CALLCONV SplatY(A_VECTOR v) noexcept;
        /// <summary>
        /// Returns a vector with all elements replicating the z component of the provided vector
        /// </summary>
        /// <param name="v">Vector to copy the z component from</param>
        /// <returns></returns>
        VECTOR VEC_CALLCONV SplatZ(A_VECTOR v) noexcept;
        /// <summary>
        /// Returns a vector with all elements replicating the w component of the provided vector
        /// </summary>
        /// <param name="v">Vector to copy the w component from</param>
        /// <returns></returns>
        VECTOR VEC_CALLCONV SplatW(A_VECTOR v) noexcept;
        /// <summary>
        /// Returns a vector of 1.0f, 1.0f, 1.0f, 1.0f
        /// </summary>
        /// <returns></returns>
        VECTOR VEC_CALLCONV SplatOne() noexcept;
        /// <summary>
        /// Returns a vector of INF, INF, INF, INF
        /// </summary>
        /// <returns></returns>
        VECTOR VEC_CALLCONV SplatInfinity() noexcept;
        /// <summary>
        /// Returns a vector of Q_NAN, Q_NAN, Q_NAN, Q_NAN
        /// </summary>
        /// <returns></returns>
        VECTOR VEC_CALLCONV SplatQNaN() noexcept;
        /// <summary>
        /// Returns a vector of 1.192092896e-7f, 1.192092896e-7f, 1.192092896e-7f, 1.192092896e-7f
        /// </summary>
        /// <returns></returns>
        VECTOR VEC_CALLCONV SplatEpsilon() noexcept;
        /// <summary>
        /// Returns a vector of -0.0f, -0.0f, -0.0f -0.0f
        /// </summary>
        /// <returns></returns>
        VECTOR VEC_CALLCONV SplatSignMask() noexcept;

        /// <summary>
        /// <para>Returns the float value at the given index in the provided vector</para>
        /// <para>**It is not recommended to use this function as it is less performant that its specified component counterparts**</para>
        /// </summary>
        /// <param name="v">Vector to get element value from</param>
        /// <param name="index">Index of element to get, 0 for X</param>
        /// <returns></returns>
        float VEC_CALLCONV GetByIndex(A_VECTOR v, size_t index) noexcept;
        /// <summary>
        /// Return the X component of the provided vector in a FPU register
        /// </summary>
        /// <param name="v">Vector to get X component from</param>
        /// <returns></returns>
        float VEC_CALLCONV GetX(A_VECTOR v) noexcept;
        /// <summary>
        /// Return the Y component of the provided vector in a FPU register
        /// </summary>
        /// <param name="v">Vector to get Y component from</param>
        /// <returns></returns>
        float VEC_CALLCONV GetY(A_VECTOR v) noexcept;
        /// <summary>
        /// Return the Z component of the provided vector in a FPU register
        /// </summary>
        /// <param name="v">Vector to get Z component from</param>
        /// <returns></returns>
        float VEC_CALLCONV GetZ(A_VECTOR v) noexcept;
        /// <summary>
        /// Return the W component of the provided vector in a FPU register
        /// </summary>
        /// <param name="v">Vector to get W component from</param>
        /// <returns></returns>
        float VEC_CALLCONV GetW(A_VECTOR v) noexcept;

        /// <summary>
        /// <para>Stores the float value at the given index in the provided vector into the memory address pointed to by the float pointer</para>
        /// <para>**It is not recommended to use this function as it is less performant that its specified component counterparts**</para>
        /// </summary>
        /// <param name="f">Float pointer to store component at</param>
        /// <param name="v">Vector to get component from</param>
        /// <param name="index">Index of component in vector to get, 0 for X</param>
        /// <returns></returns>
        void VEC_CALLCONV GetByIndexPtr(_Out_ float* f, _In_ A_VECTOR v, _In_ size_t index) noexcept;
        /// <summary>
        /// Store the X component of the provided vector into the address pointed to by the float pointer
        /// </summary>
        /// <param name="x">Float pointer to store X at</param>
        /// <param name="v">Vector to get X from</param>
        /// <returns></returns>
        void VEC_CALLCONV GetXPtr(_Out_ float* x, _In_ A_VECTOR v) noexcept;
        /// <summary>
        /// Store the Y component of the provided vector into the address pointed to by the float pointer
        /// </summary>
        /// <param name="y">Float pointer to store Y at</param>
        /// <param name="v">Vector to get Y from</param>
        /// <returns></returns>
        void VEC_CALLCONV GetYPtr(_Out_ float* y, _In_ A_VECTOR v) noexcept;
        /// <summary>
        /// Store the Z component of the provided vector into the address pointed to by the float pointer
        /// </summary>
        /// <param name="z">Float pointer to store Z at</param>
        /// <param name="v">Vector to get Z from</param>
        /// <returns></returns>
        void VEC_CALLCONV GetZPtr(_Out_ float* z, _In_ A_VECTOR v) noexcept;
        /// <summary>
        /// Store the W component of the provided vector into the address pointed to by the float pointer
        /// </summary>
        /// <param name="w">Float pointer to store W at</param>
        /// <param name="v">Vector to get W from</param>
        /// <returns></returns>
        void VEC_CALLCONV GetWPtr(_Out_ float* w, _In_ A_VECTOR v) noexcept;

        /// <summary>
        /// <para>Returns an integer value from the provided vector at the index specified</para>
        /// <para>**It is not recommended to use this function as it is less performant that its specified component counterparts**</para>
        /// </summary>
        /// <param name="v">VectorThe vector to get the element from</param>
        /// <param name="index">The index of the element to get, 0 for X</param>
        /// <returns></returns>
        uint32_t VEC_CALLCONV GetIntByIndex(A_VECTOR v, size_t index) noexcept;
        /// <summary>
        /// Return the X component of the provided vector in an integer register
        /// </summary>
        /// <param name="v">Vector to get X component from</param>
        /// <returns></returns>
        uint32_t VEC_CALLCONV GetIntX(A_VECTOR v) noexcept;
        /// <summary>
        /// Return the Y component of the provided vector in an integer register
        /// </summary>
        /// <param name="v">Vector to get Y component from</param>
        /// <returns></returns>
        uint32_t VEC_CALLCONV GetIntY(A_VECTOR v) noexcept;
        /// <summary>
        /// Return the Z component of the provided vector in an integer register
        /// </summary>
        /// <param name="v">Vector to get Z component from</param>
        /// <returns></returns>
        uint32_t VEC_CALLCONV GetIntZ(A_VECTOR v) noexcept;
        /// <summary>
        /// Return the W component of the provided vector in an integer register
        /// </summary>
        /// <param name="v">Vector to get W component from</param>
        /// <returns></returns>
        uint32_t VEC_CALLCONV GetIntW(A_VECTOR v) noexcept;

        /// <summary>
        /// Store the element of the provided vector, indexed by index into a 32 bit integer location in memory
        /// </summary>
        /// <param name="i">Integer pointer to store element at</param>
        /// <param name="v">Vector to extract element from</param>
        /// <param name="index">Index of vector to extract, 0 for X</param>
        /// <returns></returns>
        void VEC_CALLCONV GetIntByIndexPtr(_Out_ uint32_t* i, _In_ A_VECTOR v, _In_ size_t index) noexcept;
        /// <summary>
        /// Store the X component of the provided vector into a 32 pit integer location in memory
        /// </summary>
        /// <param name="x">Integer pointer to store X into</param>
        /// <param name="v">Vector to extract X from</param>
        /// <returns></returns>
        void VEC_CALLCONV GetIntXPtr(_Out_ uint32_t* x, _In_ A_VECTOR v) noexcept;
        /// <summary>
        /// Store the Y component of the provided vector into a 32 pit integer location in memory
        /// </summary>
        /// <param name="y">Integer pointer to store Y into</param>
        /// <param name="v">Vector to extract Y from</param>
        /// <returns></returns>
        void VEC_CALLCONV GetIntYPtr(_Out_ uint32_t* y, _In_ A_VECTOR v) noexcept;
        /// <summary>
        /// Store the Z component of the provided vector into a 32 pit integer location in memory
        /// </summary>
        /// <param name="z">Integer pointer to store Z into</param>
        /// <param name="v">Vector to extract Z from</param>
        /// <returns></returns>
        void VEC_CALLCONV GetIntZPtr(_Out_ uint32_t* z, _In_ A_VECTOR v) noexcept;
        /// <summary>
        /// Store the W component of the provided vector into a 32 pit integer location in memory
        /// </summary>
        /// <param name="w">Integer pointer to store W into</param>
        /// <param name="v">Vector to extract W from</param>
        /// <returns></returns>
        void VEC_CALLCONV GetIntWPtr(_Out_ uint32_t* w, _In_ A_VECTOR v) noexcept;

        /// <summary>
        /// Set a single indexed floating point component
        /// </summary>
        /// <param name="v">Vector to start from, returned vector will copy this vector in all but the specified indicies</param>
        /// <param name="f">Value to set specified index with</param>
        /// <param name="index">Index to set</param>
        /// <returns></returns>
        VECTOR VEC_CALLCONV SetByIndex(A_VECTOR v, float f, size_t index) noexcept;
        /// <summary>
        /// Set the X component of a vector with the passed float value
        /// </summary>
        /// <param name="v">Vector to start from, returned vector will copy this vector in all the but the X component</param>
        /// <param name="x">Value to set X component to</param>
        /// <returns></returns>
        VECTOR VEC_CALLCONV SetX(A_VECTOR v, float x) noexcept;
        /// <summary>
        /// Set the Y component of a vector with the passed float value
        /// </summary>
        /// <param name="v">Vector to start from, returned vector will copy this vector in all the but the Y component</param>
        /// <param name="y">Value to set Y component to</param>
        /// <returns></returns>
        VECTOR VEC_CALLCONV SetY(A_VECTOR v, float y) noexcept;
        /// <summary>
        /// Set the Z component of a vector with the passed float value
        /// </summary>
        /// <param name="v">Vector to start from, returned vector will copy this vector in all the but the Z component</param>
        /// <param name="z">Value to set Z component to</param>
        /// <returns></returns>
        VECTOR VEC_CALLCONV SetZ(A_VECTOR v, float z) noexcept;
        /// <summary>
        /// Set the W component of a vector with the passed float value
        /// </summary>
        /// <param name="v">Vector to start from, returned vector will copy this vector in all the but the W component</param>
        /// <param name="w">Value to set W component to</param>
        /// <returns></returns>
        VECTOR VEC_CALLCONV SetW(A_VECTOR v, float w) noexcept;

        /// <summary>
        /// Sets a single indexed component of a vector to a floating point value passed by pointer
        /// </summary>
        /// <param name="v">Vector to start from, returned vector will copy this vector in all the but the indexed component</param>
        /// <param name="f">Pointer to value to set indexed component to</param>
        /// <param name="index">Index of component to set, 0 for X</param>
        /// <returns></returns>
        VECTOR VEC_CALLCONV SetByIndexPtr(_In_ A_VECTOR v, _In_ const float* f, _In_ size_t index) noexcept;
        /// <summary>
        /// Sets the X component of a vector to a floating point value passed by pointer
        /// </summary>
        /// <param name="v">Vector to start from, returned vector will copy this vector in all the but the X component</param>
        /// <param name="x">Pointer to value to set the X component to</param>
        /// <returns></returns>
        VECTOR VEC_CALLCONV SetXPtr(_In_ A_VECTOR v, _In_ const float* x) noexcept;
        /// <summary>
        /// Sets the Y component of a vector to a floating point value passed by pointer
        /// </summary>
        /// <param name="v">Vector to start from, returned vector will copy this vector in all the but the Y component</param>
        /// <param name="y">Pointer to value to set the Y component to</param>
        /// <returns></returns>
        VECTOR VEC_CALLCONV SetYPtr(_In_ A_VECTOR v, _In_ const float* y) noexcept;
        /// <summary>
        /// Sets the Z component of a vector to a floating point value passed by pointer
        /// </summary>
        /// <param name="v">Vector to start from, returned vector will copy this vector in all the but the Z component</param>
        /// <param name="z">Pointer to value to set the Z component to</param>
        /// <returns></returns>
        VECTOR VEC_CALLCONV SetZPtr(_In_ A_VECTOR v, _In_ const float* z) noexcept;
        /// <summary>
        /// Sets the W component of a vector to a floating point value passed by pointer
        /// </summary>
        /// <param name="v">Vector to start from, returned vector will copy this vector in all the but the W component</param>
        /// <param name="w">Pointer to value to set the W component to</param>
        /// <returns></returns>
        VECTOR VEC_CALLCONV SetWPtr(_In_ A_VECTOR v, _In_ const float* w) noexcept;

        /// <summary>
        /// Set a single indexed integer component
        /// </summary>
        /// <param name="v">Vector to start from, returned vector will copy this vector in all but the specified indicies</param>
        /// <param name="i">Value to set specified index with</param>
        /// <param name="index">Index to set</param>
        /// <returns></returns>
        VECTOR VEC_CALLCONV SetIntByIndex(A_VECTOR v, uint32_t i, size_t index) noexcept;
        /// <summary>
        /// Set the X component of a vector with the passed int value
        /// </summary>
        /// <param name="v">Vector to start from, returned vector will copy this vector in all the but the X component</param>
        /// <param name="x">Value to set X component to</param>
        /// <returns></returns>
        VECTOR VEC_CALLCONV SetIntX(A_VECTOR v, uint32_t x) noexcept;
        /// <summary>
        /// Set the Y component of a vector with the passed int value
        /// </summary>
        /// <param name="v">Vector to start from, returned vector will copy this vector in all the but the Y component</param>
        /// <param name="y">Value to set Y component to</param>
        /// <returns></returns>
        VECTOR VEC_CALLCONV SetIntY(A_VECTOR v, uint32_t y) noexcept;
        /// <summary>
        /// Set the Z component of a vector with the passed int value
        /// </summary>
        /// <param name="v">Vector to start from, returned vector will copy this vector in all the but the Z component</param>
        /// <param name="z">Value to set Z component to</param>
        /// <returns></returns>
        VECTOR VEC_CALLCONV SetIntZ(A_VECTOR v, uint32_t z) noexcept;
        /// <summary>
        /// Set the W component of a vector with the passed int value
        /// </summary>
        /// <param name="v">Vector to start from, returned vector will copy this vector in all the but the W component</param>
        /// <param name="w">Value to set W component to</param>
        /// <returns></returns>
        VECTOR VEC_CALLCONV SetIntW(A_VECTOR v, uint32_t w) noexcept;

        /// <summary>
        /// Set a single indexed integer component to the value passed by pointer
        /// </summary>
        /// <param name="v">Vector to start from, returned vector will copy this vector in all but the specified indicies</param>
        /// <param name="i">Pointer to value to set specified index with</param>
        /// <param name="index">Index to set</param>
        /// <returns></returns>
        VECTOR VEC_CALLCONV SetIntByIndexPtr(_In_ A_VECTOR v, _In_ const uint32_t* i, _In_ size_t index) noexcept;
        /// <summary>
        /// Set the X component of a vector with the value passed by pointer
        /// </summary>
        /// <param name="v">Vector to start from, returned vector will copy this vector in all the but the X component</param>
        /// <param name="x">Pointer to value to set X component to</param>
        /// <returns></returns>
        VECTOR VEC_CALLCONV SetIntXPtr(_In_ A_VECTOR v, _In_ const uint32_t* x) noexcept;
        /// <summary>
        /// Set the Y component of a vector with the value passed by pointer
        /// </summary>
        /// <param name="v">Vector to start from, returned vector will copy this vector in all the but the Y component</param>
        /// <param name="y">Pointer to value to set Y component to</param>
        /// <returns></returns>
        VECTOR VEC_CALLCONV SetIntYPtr(_In_ A_VECTOR v, _In_ const uint32_t* y) noexcept;
        /// <summary>
        /// Set the Z component of a vector with the value passed by pointer
        /// </summary>
        /// <param name="v">Vector to start from, returned vector will copy this vector in all the but the Z component</param>
        /// <param name="z">Pointer to value to set Z component to</param>
        /// <returns></returns>
        VECTOR VEC_CALLCONV SetIntZPtr(_In_ A_VECTOR v, _In_ const uint32_t* z) noexcept;
        /// <summary>
        /// Set the W component of a vector with the value passed by pointer
        /// </summary>
        /// <param name="v">Vector to start from, returned vector will copy this vector in all the but the W component</param>
        /// <param name="w">Pointer to value to set W component to</param>
        /// <returns></returns>
        VECTOR VEC_CALLCONV SetIntWPtr(_In_ A_VECTOR v, _In_ const uint32_t* w) noexcept;

        VECTOR VEC_CALLCONV Swizzle(A_VECTOR v, uint32_t E0, uint32_t E1, uint32_t E2, uint32_t E3) noexcept;
        VECTOR VEC_CALLCONV Permute(A_VECTOR V1, A_VECTOR V2, uint32_t permuteX, uint32_t permuteY, uint32_t permuteZ, uint32_t permuteW) noexcept;
        /// <summary>
        /// <para>Creates a control vector to be used in <seealso cref="Select"/></para>
        /// <para>A value of zero for an index causes the corresponding component from the first vector to be selected, whereas a one causes the component from the second vector to be selected</para>
        /// </summary>
        /// <param name="vectorIndex0">Controls selection for the first component of the vectors involved in the selection operation</param>
        /// <param name="vectorIndex1">Controls selection for the second component of the vectors involved in the selection operation</param>
        /// <param name="vectorIndex2">Controls selection for the third component of the vectors involved in the selection operation</param>
        /// <param name="vectorIndex3">Controls selection for the fourth component of the vectors involved in the selection operation</param>
        /// <returns>Control vector to be passed to <seealso cref="Select"/></returns>
        VECTOR VEC_CALLCONV SelectControl(uint32_t vectorIndex0, uint32_t vectorIndex1, uint32_t vectorIndex2, uint32_t vectorIndex3) noexcept;
        VECTOR VEC_CALLCONV Select(A_VECTOR V1, A_VECTOR V2, A_VECTOR control) noexcept;
        VECTOR VEC_CALLCONV MergeXY(A_VECTOR V1, A_VECTOR V2) noexcept;
        VECTOR VEC_CALLCONV MergeZW(A_VECTOR V1, A_VECTOR V2) noexcept;

        VECTOR VEC_CALLCONV ShiftLeft(A_VECTOR V1, A_VECTOR V2, uint32_t elements) noexcept;
        VECTOR VEC_CALLCONV RotateLeft(A_VECTOR v, uint32_t elements) noexcept;
        VECTOR VEC_CALLCONV RotateRight(A_VECTOR v, uint32_t elements) noexcept;
        VECTOR VEC_CALLCONV Insert(A_VECTOR vDestination, A_VECTOR vSource, uint32_t VSLeftRotateElements, uint32_t select0, uint32_t select1, uint32_t select2, uint32_t select3) noexcept;

        VECTOR VEC_CALLCONV Equal(A_VECTOR V1, A_VECTOR V2) noexcept;
        VECTOR VEC_CALLCONV EqualR(_Out_ uint32_t* pCR, _In_ A_VECTOR V1, _In_ A_VECTOR V2) noexcept;
        /// <summary>
        /// <para>Treating the components of the vectors as integers, compare them for equality</para>
        /// <para>This is useful for comparing control vectors and result vectors returned from other comparison operators</para>
        /// </summary>
        /// <param name="V1">First integer vector to compare</param>
        /// <param name="V2">Second integer vector to compare</param>
        /// <returns>
        /// <para> whose components indicate which element pairs of the argument vectors are equal or not-equal</para>
        /// <para>Equal component pairs indicated with value of 0xFFFFFFFF, not-equal with 0x00000000</para>
        /// </returns>
        VECTOR VEC_CALLCONV EqualInt(A_VECTOR V1, A_VECTOR V2) noexcept;
        VECTOR VEC_CALLCONV EqualIntR(_Out_ uint32_t* pCR, _In_ A_VECTOR V1, _In_ A_VECTOR V2) noexcept;
        VECTOR VEC_CALLCONV NearEqual(A_VECTOR V1, A_VECTOR V2, A_VECTOR epsilon) noexcept;
        VECTOR VEC_CALLCONV NotEqual(A_VECTOR V1, A_VECTOR V2) noexcept;
        VECTOR VEC_CALLCONV NotEqualInt(A_VECTOR V1, A_VECTOR V2) noexcept;
        VECTOR VEC_CALLCONV Greater(A_VECTOR V1, A_VECTOR V2) noexcept;
        VECTOR VEC_CALLCONV GreaterR(_Out_ uint32_t* pCR, _In_ A_VECTOR V1, _In_ A_VECTOR V2) noexcept;
        VECTOR VEC_CALLCONV GreaterOrEqual(A_VECTOR V1, A_VECTOR V2) noexcept;
        VECTOR VEC_CALLCONV GreaterOrEqualR(_Out_ uint32_t* pCR, _In_ A_VECTOR V1, _In_ A_VECTOR V2) noexcept;
        VECTOR VEC_CALLCONV Less(A_VECTOR V1, A_VECTOR V2) noexcept;
        VECTOR VEC_CALLCONV LessOrEqual(A_VECTOR V1, A_VECTOR V2) noexcept;
        VECTOR VEC_CALLCONV InBounds(A_VECTOR v, A_VECTOR bounds) noexcept;
        VECTOR VEC_CALLCONV InBoundsR(_Out_ uint32_t* pCR, _In_ A_VECTOR v, _In_ A_VECTOR bounds) noexcept;

        VECTOR VEC_CALLCONV IsNaN(A_VECTOR v) noexcept;
        VECTOR VEC_CALLCONV IsInfinite(A_VECTOR v) noexcept;

        VECTOR VEC_CALLCONV Min(A_VECTOR V1, A_VECTOR V2) noexcept;
        VECTOR VEC_CALLCONV Max(A_VECTOR V1, A_VECTOR V2) noexcept;
        VECTOR VEC_CALLCONV Round(A_VECTOR v) noexcept;
        VECTOR VEC_CALLCONV Truncate(A_VECTOR v) noexcept;
        VECTOR VEC_CALLCONV Floor(A_VECTOR v) noexcept;
        VECTOR VEC_CALLCONV Ceiling(A_VECTOR v) noexcept;
        VECTOR VEC_CALLCONV Clamp(A_VECTOR v, A_VECTOR min, A_VECTOR max) noexcept;
        VECTOR VEC_CALLCONV Saturate(A_VECTOR v) noexcept;

        VECTOR VEC_CALLCONV AndInt(A_VECTOR V1, A_VECTOR V2) noexcept;
        VECTOR VEC_CALLCONV AndCInt(A_VECTOR V1, A_VECTOR V2) noexcept;
        VECTOR VEC_CALLCONV OrInt(A_VECTOR V1, A_VECTOR V2) noexcept;
        VECTOR VEC_CALLCONV NorInt(A_VECTOR V1, A_VECTOR V2) noexcept;
        VECTOR VEC_CALLCONV XorInt(A_VECTOR V1, A_VECTOR V2) noexcept;

        VECTOR VEC_CALLCONV Negate(A_VECTOR v) noexcept;
        VECTOR VEC_CALLCONV Add(A_VECTOR V1, A_VECTOR V2) noexcept;
        VECTOR VEC_CALLCONV Sum(A_VECTOR v) noexcept;
        VECTOR VEC_CALLCONV AddAngles(A_VECTOR V1, A_VECTOR V2) noexcept;
        VECTOR VEC_CALLCONV Subtract(A_VECTOR V1, A_VECTOR V2) noexcept;
        VECTOR VEC_CALLCONV SubtractAngles(A_VECTOR V1, A_VECTOR V2) noexcept;
        VECTOR VEC_CALLCONV Multiply(A_VECTOR V1, A_VECTOR V2) noexcept;
        /// <summary>
        /// Multiply the first two vectors and add the third to that product
        /// </summary>
        /// <param name="V1">First of two vectors to be multiplied</param>
        /// <param name="V2">Second of two vectors to be multiplied</param>
        /// <param name="V3">Third vector to add to the product</param>
        /// <returns>Vector calculated from expression V1 * V2 + V3</returns>
        VECTOR VEC_CALLCONV MultiplyAdd(A_VECTOR V1, A_VECTOR V2, A_VECTOR V3) noexcept;
        VECTOR VEC_CALLCONV Divide(A_VECTOR V1, A_VECTOR V2) noexcept;
        /// <summary>
        /// Multiply the first two vectors and subtract that product from the third
        /// </summary>
        /// <param name="V1">First of two vectors to be multiplied</param>
        /// <param name="V2">Second of two vectors to be multiplied</param>
        /// <param name="V3">Third vector to subtract the product of the first two from</param>
        /// <returns>Vector calculated from expression V3 - V1 * V2</returns>
        VECTOR VEC_CALLCONV NegativeMultiplySubtract(A_VECTOR V1, A_VECTOR V2, A_VECTOR V3) noexcept;
        VECTOR VEC_CALLCONV Scale(A_VECTOR v, float scaleFactor) noexcept;
        VECTOR VEC_CALLCONV ReciprocalEst(A_VECTOR v) noexcept;
        VECTOR VEC_CALLCONV Reciprocal(A_VECTOR v) noexcept;
        VECTOR VEC_CALLCONV SqrtEst(A_VECTOR v) noexcept;
        VECTOR VEC_CALLCONV Sqrt(A_VECTOR v) noexcept;
        VECTOR VEC_CALLCONV ReciprocalSrtEst(A_VECTOR v) noexcept;
        VECTOR VEC_CALLCONV ReciprocalSqrt(A_VECTOR v) noexcept;
        VECTOR VEC_CALLCONV Exp2(A_VECTOR v) noexcept;
        VECTOR VEC_CALLCONV Exp10(A_VECTOR v) noexcept;
        VECTOR VEC_CALLCONV ExpE(A_VECTOR v) noexcept;
        VECTOR VEC_CALLCONV Exp(A_VECTOR v) noexcept;
        VECTOR VEC_CALLCONV Log2(A_VECTOR v) noexcept;
        VECTOR VEC_CALLCONV Log10(A_VECTOR v) noexcept;
        VECTOR VEC_CALLCONV LogE(A_VECTOR v) noexcept;
        VECTOR VEC_CALLCONV Log(A_VECTOR v) noexcept;
        VECTOR VEC_CALLCONV Pow(A_VECTOR V1, A_VECTOR V2) noexcept;
        VECTOR VEC_CALLCONV Abs(A_VECTOR v) noexcept;
        VECTOR VEC_CALLCONV Mod(A_VECTOR V1, A_VECTOR V2) noexcept;
        VECTOR VEC_CALLCONV ModAngles(A_VECTOR angles) noexcept;
        VECTOR VEC_CALLCONV Sine(A_VECTOR v) noexcept;
        VECTOR VEC_CALLCONV SineEst(A_VECTOR v) noexcept;
        VECTOR VEC_CALLCONV Cos(A_VECTOR v) noexcept;
        VECTOR VEC_CALLCONV CosEst(A_VECTOR v) noexcept;
        void VEC_CALLCONV SineCos(_Out_ VECTOR* pSine, _Out_ VECTOR* pCos, _In_ A_VECTOR v) noexcept;
        void VEC_CALLCONV SineCosEst(_Out_ VECTOR* pSine, _Out_ VECTOR* pCos, _In_ A_VECTOR v) noexcept;
        VECTOR VEC_CALLCONV Tan(A_VECTOR v) noexcept;
        VECTOR VEC_CALLCONV TanEst(A_VECTOR v) noexcept;
        VECTOR VEC_CALLCONV SineH(A_VECTOR v) noexcept;
        VECTOR VEC_CALLCONV CosH(A_VECTOR v) noexcept;
        VECTOR VEC_CALLCONV TanH(A_VECTOR v) noexcept;
        VECTOR VEC_CALLCONV ASine(A_VECTOR v) noexcept;
        VECTOR VEC_CALLCONV ASineEst(A_VECTOR v) noexcept;
        VECTOR VEC_CALLCONV ACos(A_VECTOR v) noexcept;
        VECTOR VEC_CALLCONV ACosEst(A_VECTOR v) noexcept;
        VECTOR VEC_CALLCONV ATan(A_VECTOR v) noexcept;
        VECTOR VEC_CALLCONV ATanEst(A_VECTOR v) noexcept;
        VECTOR VEC_CALLCONV ATan2(A_VECTOR y, A_VECTOR x) noexcept;
        VECTOR VEC_CALLCONV ATan2Est(A_VECTOR y, A_VECTOR x) noexcept;
        VECTOR VEC_CALLCONV Lerp(A_VECTOR V1, A_VECTOR V2, float t) noexcept;
        VECTOR VEC_CALLCONV LerpV(A_VECTOR V1, A_VECTOR V2, A_VECTOR VT) noexcept;
        VECTOR VEC_CALLCONV Hermite(A_VECTOR postion1, A_VECTOR tangent1, A_VECTOR postion2, B_VECTOR tangent2, float t) noexcept;
        VECTOR VEC_CALLCONV HermiteV(A_VECTOR postion1, A_VECTOR tangent1, A_VECTOR position2, B_VECTOR tangent2, C_VECTOR VT) noexcept;
        VECTOR VEC_CALLCONV CatmullRom(A_VECTOR position1, A_VECTOR position2, A_VECTOR position3, B_VECTOR postion4, float t) noexcept;
        VECTOR VEC_CALLCONV CatmullRomV(A_VECTOR postion1, A_VECTOR position2, A_VECTOR position3, B_VECTOR position4, C_VECTOR VT) noexcept;
        VECTOR VEC_CALLCONV BaryCentric(A_VECTOR postion1, A_VECTOR postion2, A_VECTOR postion3, float f, float g) noexcept;
        VECTOR VEC_CALLCONV BaryCentricV(A_VECTOR position1, A_VECTOR postion2, A_VECTOR postion3, B_VECTOR VF, C_VECTOR VG) noexcept;
	}

	// The purpose of the following global constants is to prevent redundant
	// reloading of the constants when they are referenced by more than one
	// separate inline math routine called within the same function. Declaring
	// a constant locally within a routine is sufficient to prevent redundant
	// reloads of that constant when that single routine is called multiple
	// times in a function, but if the constant is used (and declared) in a
	// separate math routine it would be reloaded.

#ifndef VEC_GLOBCONST
#if defined(__GNUC__) && !defined(__MINGW32__)
#define VEC_GLOBCONST extern const __attribute__((weak))
#else
#define VEC_GLOBCONST extern const __declspec(selectany)
#endif
#endif

    VEC_GLOBCONST VECTOR_F32 g_SinCoefficients0 = { { { -0.16666667f, +0.0083333310f, -0.00019840874f, +2.7525562e-06f } } };
    VEC_GLOBCONST VECTOR_F32 g_SinCoefficients1 = { { { -2.3889859e-08f, -0.16665852f /*Est1*/, +0.0083139502f /*Est2*/, -0.00018524670f /*Est3*/ } } };
    VEC_GLOBCONST VECTOR_F32 g_CosCoefficients0 = { { { -0.5f, +0.041666638f, -0.0013888378f, +2.4760495e-05f } } };
    VEC_GLOBCONST VECTOR_F32 g_CosCoefficients1 = { { { -2.6051615e-07f, -0.49992746f /*Est1*/, +0.041493919f /*Est2*/, -0.0012712436f /*Est3*/ } } };
    VEC_GLOBCONST VECTOR_F32 g_TanCoefficients0 = { { { 1.0f, 0.333333333f, 0.133333333f, 5.396825397e-2f } } };
    VEC_GLOBCONST VECTOR_F32 g_TanCoefficients1 = { { { 2.186948854e-2f, 8.863235530e-3f, 3.592128167e-3f, 1.455834485e-3f } } };
    VEC_GLOBCONST VECTOR_F32 g_TanCoefficients2 = { { { 5.900274264e-4f, 2.391290764e-4f, 9.691537707e-5f, 3.927832950e-5f } } };
    VEC_GLOBCONST VECTOR_F32 g_ArcCoefficients0 = { { { +1.5707963050f, -0.2145988016f, +0.0889789874f, -0.0501743046f } } };
    VEC_GLOBCONST VECTOR_F32 g_ArcCoefficients1 = { { { +0.0308918810f, -0.0170881256f, +0.0066700901f, -0.0012624911f } } };
    VEC_GLOBCONST VECTOR_F32 g_ATanCoefficients0 = { { { -0.3333314528f, +0.1999355085f, -0.1420889944f, +0.1065626393f } } };
    VEC_GLOBCONST VECTOR_F32 g_ATanCoefficients1 = { { { -0.0752896400f, +0.0429096138f, -0.0161657367f, +0.0028662257f } } };
    VEC_GLOBCONST VECTOR_F32 g_ATanEstCoefficients0 = { { { +0.999866f, +0.999866f, +0.999866f, +0.999866f } } };
    VEC_GLOBCONST VECTOR_F32 g_ATanEstCoefficients1 = { { { -0.3302995f, +0.180141f, -0.085133f, +0.0208351f } } };
    VEC_GLOBCONST VECTOR_F32 g_TanEstCoefficients = { { { 2.484f, -1.954923183e-1f, 2.467401101f, _1OVERPI } } };
    VEC_GLOBCONST VECTOR_F32 g_ArcEstCoefficients = { { { +1.5707288f, -0.2121144f, +0.0742610f, -0.0187293f } } };
    VEC_GLOBCONST VECTOR_F32 g_PiConstants0 = { { { _PI, _2PI, _1OVERPI, _1OVER2PI } } };
    VEC_GLOBCONST VECTOR_F32 g_IdentityR0 = { { { 1.0f, 0.0f, 0.0f, 0.0f } } };
    VEC_GLOBCONST VECTOR_F32 g_IdentityR1 = { { { 0.0f, 1.0f, 0.0f, 0.0f } } };
    VEC_GLOBCONST VECTOR_F32 g_IdentityR2 = { { { 0.0f, 0.0f, 1.0f, 0.0f } } };
    VEC_GLOBCONST VECTOR_F32 g_IdentityR3 = { { { 0.0f, 0.0f, 0.0f, 1.0f } } };
    VEC_GLOBCONST VECTOR_F32 g_NegIdentityR0 = { { { -1.0f, 0.0f, 0.0f, 0.0f } } };
    VEC_GLOBCONST VECTOR_F32 g_NegIdentityR1 = { { { 0.0f, -1.0f, 0.0f, 0.0f } } };
    VEC_GLOBCONST VECTOR_F32 g_NegIdentityR2 = { { { 0.0f, 0.0f, -1.0f, 0.0f } } };
    VEC_GLOBCONST VECTOR_F32 g_NegIdentityR3 = { { { 0.0f, 0.0f, 0.0f, -1.0f } } };
    VEC_GLOBCONST VECTOR_U32 g_NegativeZero = { { { 0x80000000, 0x80000000, 0x80000000, 0x80000000 } } };
    VEC_GLOBCONST VECTOR_U32 g_Negate3 = { { { 0x80000000, 0x80000000, 0x80000000, 0x00000000 } } };
    VEC_GLOBCONST VECTOR_U32 g_MaskXY = { { { 0xFFFFFFFF, 0xFFFFFFFF, 0x00000000, 0x00000000 } } };
    VEC_GLOBCONST VECTOR_U32 g_Mask3 = { { { 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0x00000000 } } };
    VEC_GLOBCONST VECTOR_U32 g_MaskX = { { { 0xFFFFFFFF, 0x00000000, 0x00000000, 0x00000000 } } };
    VEC_GLOBCONST VECTOR_U32 g_MaskY = { { { 0x00000000, 0xFFFFFFFF, 0x00000000, 0x00000000 } } };
    VEC_GLOBCONST VECTOR_U32 g_MaskZ = { { { 0x00000000, 0x00000000, 0xFFFFFFFF, 0x00000000 } } };
    VEC_GLOBCONST VECTOR_U32 g_MaskW = { { { 0x00000000, 0x00000000, 0x00000000, 0xFFFFFFFF } } };
    VEC_GLOBCONST VECTOR_F32 g_One = { { { 1.0f, 1.0f, 1.0f, 1.0f } } };
    VEC_GLOBCONST VECTOR_F32 g_One3 = { { { 1.0f, 1.0f, 1.0f, 0.0f } } };
    VEC_GLOBCONST VECTOR_F32 g_Zero = { { { 0.0f, 0.0f, 0.0f, 0.0f } } };
    VEC_GLOBCONST VECTOR_F32 g_Two = { { { 2.f, 2.f, 2.f, 2.f } } };
    VEC_GLOBCONST VECTOR_F32 g_Four = { { { 4.f, 4.f, 4.f, 4.f } } };
    VEC_GLOBCONST VECTOR_F32 g_Six = { { { 6.f, 6.f, 6.f, 6.f } } };
    VEC_GLOBCONST VECTOR_F32 g_NegativeOne = { { { -1.0f, -1.0f, -1.0f, -1.0f } } };
    VEC_GLOBCONST VECTOR_F32 g_OneHalf = { { { 0.5f, 0.5f, 0.5f, 0.5f } } };
    VEC_GLOBCONST VECTOR_F32 g_NegativeOneHalf = { { { -0.5f, -0.5f, -0.5f, -0.5f } } };
    VEC_GLOBCONST VECTOR_F32 g_NegativeTwoPi = { { { -_2PI, -_2PI, -_2PI, -_2PI } } };
    VEC_GLOBCONST VECTOR_F32 g_NegativePi = { { { -_PI, -_PI, -_PI, -_PI } } };
    VEC_GLOBCONST VECTOR_F32 g_HalfPi = { { { _PIOVER2, _PIOVER2, _PIOVER2, _PIOVER2 } } };
    VEC_GLOBCONST VECTOR_F32 g_Pi = { { { _PI, _PI, _PI, _PI } } };
    VEC_GLOBCONST VECTOR_F32 g_ReciprocalPi = { { { _1OVERPI, _1OVERPI, _1OVERPI, _1OVERPI } } };
    VEC_GLOBCONST VECTOR_F32 g_TwoPi = { { { _2PI, _2PI, _2PI, _2PI } } };
    VEC_GLOBCONST VECTOR_F32 g_ReciprocalTwoPi = { { { _1OVER2PI, _1OVER2PI, _1OVER2PI, _1OVER2PI } } };
    VEC_GLOBCONST VECTOR_F32 g_Epsilon = { { { 1.192092896e-7f, 1.192092896e-7f, 1.192092896e-7f, 1.192092896e-7f } } };
    VEC_GLOBCONST VECTOR_I32 g_Infinity = { { { 0x7F800000, 0x7F800000, 0x7F800000, 0x7F800000 } } };
    VEC_GLOBCONST VECTOR_I32 g_QNaN = { { { 0x7FC00000, 0x7FC00000, 0x7FC00000, 0x7FC00000 } } };
    VEC_GLOBCONST VECTOR_I32 g_QNaNTest = { { { 0x007FFFFF, 0x007FFFFF, 0x007FFFFF, 0x007FFFFF } } };
    VEC_GLOBCONST VECTOR_I32 g_AbsMask = { { { 0x7FFFFFFF, 0x7FFFFFFF, 0x7FFFFFFF, 0x7FFFFFFF } } };
    VEC_GLOBCONST VECTOR_I32 g_FltMin = { { { 0x00800000, 0x00800000, 0x00800000, 0x00800000 } } };
    VEC_GLOBCONST VECTOR_I32 g_FltMax = { { { 0x7F7FFFFF, 0x7F7FFFFF, 0x7F7FFFFF, 0x7F7FFFFF } } };
    VEC_GLOBCONST VECTOR_U32 g_NegOneMask = { { { 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF } } };
    VEC_GLOBCONST VECTOR_U32 g_MaskA8R8G8B8 = { { { 0x00FF0000, 0x0000FF00, 0x000000FF, 0xFF000000 } } };
    VEC_GLOBCONST VECTOR_U32 g_FlipA8R8G8B8 = { { { 0x00000000, 0x00000000, 0x00000000, 0x80000000 } } };
    VEC_GLOBCONST VECTOR_F32 g_FixAA8R8G8B8 = { { { 0.0f, 0.0f, 0.0f, float(0x80000000U) } } };
    VEC_GLOBCONST VECTOR_F32 g_NormalizeA8R8G8B8 = { { { 1.0f / (255.0f * float(0x10000)), 1.0f / (255.0f * float(0x100)), 1.0f / 255.0f, 1.0f / (255.0f * float(0x1000000)) } } };
    VEC_GLOBCONST VECTOR_U32 g_MaskA2B10G10R10 = { { { 0x000003FF, 0x000FFC00, 0x3FF00000, 0xC0000000 } } };
    VEC_GLOBCONST VECTOR_U32 g_FlipA2B10G10R10 = { { { 0x00000200, 0x00080000, 0x20000000, 0x80000000 } } };
    VEC_GLOBCONST VECTOR_F32 g_FixAA2B10G10R10 = { { { -512.0f, -512.0f * float(0x400), -512.0f * float(0x100000), float(0x80000000U) } } };
    VEC_GLOBCONST VECTOR_F32 g_NormalizeA2B10G10R10 = { { { 1.0f / 511.0f, 1.0f / (511.0f * float(0x400)), 1.0f / (511.0f * float(0x100000)), 1.0f / (3.0f * float(0x40000000)) } } };
    VEC_GLOBCONST VECTOR_U32 g_MaskX16Y16 = { { { 0x0000FFFF, 0xFFFF0000, 0x00000000, 0x00000000 } } };
    VEC_GLOBCONST VECTOR_I32 g_FlipX16Y16 = { { { 0x00008000, 0x00000000, 0x00000000, 0x00000000 } } };
    VEC_GLOBCONST VECTOR_F32 g_FixX16Y16 = { { { -32768.0f, 0.0f, 0.0f, 0.0f } } };
    VEC_GLOBCONST VECTOR_F32 g_NormalizeX16Y16 = { { { 1.0f / 32767.0f, 1.0f / (32767.0f * 65536.0f), 0.0f, 0.0f } } };
    VEC_GLOBCONST VECTOR_U32 g_MaskX16Y16Z16W16 = { { { 0x0000FFFF, 0x0000FFFF, 0xFFFF0000, 0xFFFF0000 } } };
    VEC_GLOBCONST VECTOR_I32 g_FlipX16Y16Z16W16 = { { { 0x00008000, 0x00008000, 0x00000000, 0x00000000 } } };
    VEC_GLOBCONST VECTOR_F32 g_FixX16Y16Z16W16 = { { { -32768.0f, -32768.0f, 0.0f, 0.0f } } };
    VEC_GLOBCONST VECTOR_F32 g_NormalizeX16Y16Z16W16 = { { { 1.0f / 32767.0f, 1.0f / 32767.0f, 1.0f / (32767.0f * 65536.0f), 1.0f / (32767.0f * 65536.0f) } } };
    VEC_GLOBCONST VECTOR_F32 g_NoFraction = { { { 8388608.0f, 8388608.0f, 8388608.0f, 8388608.0f } } };
    VEC_GLOBCONST VECTOR_I32 g_MaskByte = { { { 0x000000FF, 0x000000FF, 0x000000FF, 0x000000FF } } };
    VEC_GLOBCONST VECTOR_F32 g_NegateX = { { { -1.0f, 1.0f, 1.0f, 1.0f } } };
    VEC_GLOBCONST VECTOR_F32 g_NegateY = { { { 1.0f, -1.0f, 1.0f, 1.0f } } };
    VEC_GLOBCONST VECTOR_F32 g_NegateZ = { { { 1.0f, 1.0f, -1.0f, 1.0f } } };
    VEC_GLOBCONST VECTOR_F32 g_NegateW = { { { 1.0f, 1.0f, 1.0f, -1.0f } } };
    VEC_GLOBCONST VECTOR_U32 g_Select0101 = { { { SELECT_NONE, SELECT_ALL, SELECT_NONE, SELECT_ALL } } };
    VEC_GLOBCONST VECTOR_U32 g_Select1010 = { { { SELECT_ALL, SELECT_NONE, SELECT_ALL, SELECT_NONE } } };
    VEC_GLOBCONST VECTOR_I32 g_OneHalfMinusEpsilon = { { { 0x3EFFFFFD, 0x3EFFFFFD, 0x3EFFFFFD, 0x3EFFFFFD } } };
    VEC_GLOBCONST VECTOR_U32 g_Select1000 = { { { SELECT_ALL, SELECT_NONE, SELECT_NONE, SELECT_NONE } } };
    VEC_GLOBCONST VECTOR_U32 g_Select1100 = { { { SELECT_ALL, SELECT_ALL, SELECT_NONE, SELECT_NONE } } };
    VEC_GLOBCONST VECTOR_U32 g_Select1110 = { { { SELECT_ALL, SELECT_ALL, SELECT_ALL, SELECT_NONE } } };
    VEC_GLOBCONST VECTOR_U32 g_Select1011 = { { { SELECT_ALL, SELECT_NONE, SELECT_ALL, SELECT_ALL } } };
    VEC_GLOBCONST VECTOR_F32 g_FixupY16 = { { { 1.0f, 1.0f / 65536.0f, 0.0f, 0.0f } } };
    VEC_GLOBCONST VECTOR_F32 g_FixupY16W16 = { { { 1.0f, 1.0f, 1.0f / 65536.0f, 1.0f / 65536.0f } } };
    VEC_GLOBCONST VECTOR_U32 g_FlipY = { { { 0, 0x80000000, 0, 0 } } };
    VEC_GLOBCONST VECTOR_U32 g_FlipZ = { { { 0, 0, 0x80000000, 0 } } };
    VEC_GLOBCONST VECTOR_U32 g_FlipW = { { { 0, 0, 0, 0x80000000 } } };
    VEC_GLOBCONST VECTOR_U32 g_FlipYZ = { { { 0, 0x80000000, 0x80000000, 0 } } };
    VEC_GLOBCONST VECTOR_U32 g_FlipZW = { { { 0, 0, 0x80000000, 0x80000000 } } };
    VEC_GLOBCONST VECTOR_U32 g_FlipYW = { { { 0, 0x80000000, 0, 0x80000000 } } };
    VEC_GLOBCONST VECTOR_I32 g_MaskDec4 = { { { 0x3FF, 0x3FF << 10, 0x3FF << 20, static_cast<int>(0xC0000000) } } };
    VEC_GLOBCONST VECTOR_I32 g_XorDec4 = { { { 0x200, 0x200 << 10, 0x200 << 20, 0 } } };
    VEC_GLOBCONST VECTOR_F32 g_AddUDec4 = { { { 0, 0, 0, 32768.0f * 65536.0f } } };
    VEC_GLOBCONST VECTOR_F32 g_AddDec4 = { { { -512.0f, -512.0f * 1024.0f, -512.0f * 1024.0f * 1024.0f, 0 } } };
    VEC_GLOBCONST VECTOR_F32 g_MulDec4 = { { { 1.0f, 1.0f / 1024.0f, 1.0f / (1024.0f * 1024.0f), 1.0f / (1024.0f * 1024.0f * 1024.0f) } } };
    VEC_GLOBCONST VECTOR_U32 g_MaskByte4 = { { { 0xFF, 0xFF00, 0xFF0000, 0xFF000000 } } };
    VEC_GLOBCONST VECTOR_I32 g_XorByte4 = { { { 0x80, 0x8000, 0x800000, 0x00000000 } } };
    VEC_GLOBCONST VECTOR_F32 g_AddByte4 = { { { -128.0f, -128.0f * 256.0f, -128.0f * 65536.0f, 0 } } };
    VEC_GLOBCONST VECTOR_F32 g_FixUnsigned = { { { 32768.0f * 65536.0f, 32768.0f * 65536.0f, 32768.0f * 65536.0f, 32768.0f * 65536.0f } } };
    VEC_GLOBCONST VECTOR_F32 g_MaxInt = { { { 65536.0f * 32768.0f - 128.0f, 65536.0f * 32768.0f - 128.0f, 65536.0f * 32768.0f - 128.0f, 65536.0f * 32768.0f - 128.0f } } };
    VEC_GLOBCONST VECTOR_F32 g_MaxUInt = { { { 65536.0f * 65536.0f - 256.0f, 65536.0f * 65536.0f - 256.0f, 65536.0f * 65536.0f - 256.0f, 65536.0f * 65536.0f - 256.0f } } };
    VEC_GLOBCONST VECTOR_F32 g_UnsignedFix = { { { 32768.0f * 65536.0f, 32768.0f * 65536.0f, 32768.0f * 65536.0f, 32768.0f * 65536.0f } } };
    VEC_GLOBCONST VECTOR_F32 g_srgbScale = { { { 12.92f, 12.92f, 12.92f, 1.0f } } };
    VEC_GLOBCONST VECTOR_F32 g_srgbA = { { { 0.055f, 0.055f, 0.055f, 0.0f } } };
    VEC_GLOBCONST VECTOR_F32 g_srgbA1 = { { { 1.055f, 1.055f, 1.055f, 1.0f } } };
    VEC_GLOBCONST VECTOR_I32 g_ExponentBias = { { { 127, 127, 127, 127 } } };
    VEC_GLOBCONST VECTOR_I32 g_SubnormalExponent = { { { -126, -126, -126, -126 } } };
    VEC_GLOBCONST VECTOR_I32 g_NumTrailing = { { { 23, 23, 23, 23 } } };
    VEC_GLOBCONST VECTOR_I32 g_MinNormal = { { { 0x00800000, 0x00800000, 0x00800000, 0x00800000 } } };
    VEC_GLOBCONST VECTOR_U32 g_NegInfinity = { { { 0xFF800000, 0xFF800000, 0xFF800000, 0xFF800000 } } };
    VEC_GLOBCONST VECTOR_U32 g_NegQNaN = { { { 0xFFC00000, 0xFFC00000, 0xFFC00000, 0xFFC00000 } } };
    VEC_GLOBCONST VECTOR_I32 g_Bin128 = { { { 0x43000000, 0x43000000, 0x43000000, 0x43000000 } } };
    VEC_GLOBCONST VECTOR_U32 g_BinNeg150 = { { { 0xC3160000, 0xC3160000, 0xC3160000, 0xC3160000 } } };
    VEC_GLOBCONST VECTOR_I32 g_253 = { { { 253, 253, 253, 253 } } };
    VEC_GLOBCONST VECTOR_F32 g_ExpEst1 = { { { -6.93147182e-1f, -6.93147182e-1f, -6.93147182e-1f, -6.93147182e-1f } } };
    VEC_GLOBCONST VECTOR_F32 g_ExpEst2 = { { { +2.40226462e-1f, +2.40226462e-1f, +2.40226462e-1f, +2.40226462e-1f } } };
    VEC_GLOBCONST VECTOR_F32 g_ExpEst3 = { { { -5.55036440e-2f, -5.55036440e-2f, -5.55036440e-2f, -5.55036440e-2f } } };
    VEC_GLOBCONST VECTOR_F32 g_ExpEst4 = { { { +9.61597636e-3f, +9.61597636e-3f, +9.61597636e-3f, +9.61597636e-3f } } };
    VEC_GLOBCONST VECTOR_F32 g_ExpEst5 = { { { -1.32823968e-3f, -1.32823968e-3f, -1.32823968e-3f, -1.32823968e-3f } } };
    VEC_GLOBCONST VECTOR_F32 g_ExpEst6 = { { { +1.47491097e-4f, +1.47491097e-4f, +1.47491097e-4f, +1.47491097e-4f } } };
    VEC_GLOBCONST VECTOR_F32 g_ExpEst7 = { { { -1.08635004e-5f, -1.08635004e-5f, -1.08635004e-5f, -1.08635004e-5f } } };
    VEC_GLOBCONST VECTOR_F32 g_LogEst0 = { { { +1.442693f, +1.442693f, +1.442693f, +1.442693f } } };
    VEC_GLOBCONST VECTOR_F32 g_LogEst1 = { { { -0.721242f, -0.721242f, -0.721242f, -0.721242f } } };
    VEC_GLOBCONST VECTOR_F32 g_LogEst2 = { { { +0.479384f, +0.479384f, +0.479384f, +0.479384f } } };
    VEC_GLOBCONST VECTOR_F32 g_LogEst3 = { { { -0.350295f, -0.350295f, -0.350295f, -0.350295f } } };
    VEC_GLOBCONST VECTOR_F32 g_LogEst4 = { { { +0.248590f, +0.248590f, +0.248590f, +0.248590f } } };
    VEC_GLOBCONST VECTOR_F32 g_LogEst5 = { { { -0.145700f, -0.145700f, -0.145700f, -0.145700f } } };
    VEC_GLOBCONST VECTOR_F32 g_LogEst6 = { { { +0.057148f, +0.057148f, +0.057148f, +0.057148f } } };
    VEC_GLOBCONST VECTOR_F32 g_LogEst7 = { { { -0.010578f, -0.010578f, -0.010578f, -0.010578f } } };
    VEC_GLOBCONST VECTOR_F32 g_LgE = { { { +1.442695f, +1.442695f, +1.442695f, +1.442695f } } };
    VEC_GLOBCONST VECTOR_F32 g_InvLgE = { { { +6.93147182e-1f, +6.93147182e-1f, +6.93147182e-1f, +6.93147182e-1f } } };
    VEC_GLOBCONST VECTOR_F32 g_Lg10 = { { { +3.321928f, +3.321928f, +3.321928f, +3.321928f } } };
    VEC_GLOBCONST VECTOR_F32 g_InvLg10 = { { { +3.010299956e-1f, +3.010299956e-1f, +3.010299956e-1f, +3.010299956e-1f } } };
    VEC_GLOBCONST VECTOR_F32 g_UByteMax = { { { 255.0f, 255.0f, 255.0f, 255.0f } } };
    VEC_GLOBCONST VECTOR_F32 g_ByteMin = { { { -127.0f, -127.0f, -127.0f, -127.0f } } };
    VEC_GLOBCONST VECTOR_F32 g_ByteMax = { { { 127.0f, 127.0f, 127.0f, 127.0f } } };
    VEC_GLOBCONST VECTOR_F32 g_ShortMin = { { { -32767.0f, -32767.0f, -32767.0f, -32767.0f } } };
    VEC_GLOBCONST VECTOR_F32 g_ShortMax = { { { 32767.0f, 32767.0f, 32767.0f, 32767.0f } } };
    VEC_GLOBCONST VECTOR_F32 g_UShortMax = { { { 65535.0f, 65535.0f, 65535.0f, 65535.0f } } };
}

#include <SSE2VectorConfig.inl>

#ifdef __clang__
#pragma clang diagnostic pop
#endif
#ifdef _PREFAST_
#pragma prefast(pop)
#endif
#ifdef _MSC_VER
#pragma warning(pop)
#endif

#endif // !ULTREALITY_MATH_SSE2_CONFIG_H
