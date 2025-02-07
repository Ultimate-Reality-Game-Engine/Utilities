#ifndef ULTREALITY_MATH_FLOAT3_H
#define ULTREALITY_MATH_FLOAT3_H

#include <SSE2VectorConfig.h>
#include <MATRIX.h>
#include <Float4.h>

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

namespace UltReality::Math
{
	/// <summary>
	/// 3D vector of 32-bit floating point components
	/// </summary>
	struct Float3
	{
		float x;
		float y;
		float z;

		Float3() = default;

		Float3(const Float3&) = default;
		Float3& operator=(const Float3&) = default;

		Float3(Float3&&) = default;
		Float3& operator=(Float3&&) = default;

		constexpr Float3(float _x, float _y, float _z) noexcept : x(_x), y(_y), z(_z) {}
		explicit Float3(_In_reads_(3) const float* pArray) noexcept : x(pArray[0]), y(pArray[1]), z(pArray[2]) {}

		explicit Float3(_In_ A_VECTOR v) noexcept;
		Float3& VEC_CALLCONV operator=(_In_ A_VECTOR v) noexcept;

#if (__cplusplus >= 202002L)
		bool operator==(const Float3&) const = default;
		auto operator<=>(const Float3&) const = default;
#endif

		VECTOR VEC_CALLCONV Load() noexcept;
		void VEC_CALLCONV Store(_In_ A_VECTOR v) noexcept;
	};

	/// <summary>
	/// 3D vector of 32-bit floating point components aligned on a 16 byte boundary
	/// </summary>
	ALIGNED_STRUCT(16) AFloat3 : public Float3
	{
		using Float3::Float3;

		explicit AFloat3(_In_ A_VECTOR v) noexcept;
		AFloat3& VEC_CALLCONV operator=(_In_ A_VECTOR v) noexcept;

		VECTOR VEC_CALLCONV Load() noexcept;
		void VEC_CALLCONV Store(_In_ A_VECTOR v) noexcept;
	};

	namespace VEC
	{
		VECTOR VEC_CALLCONV LoadFloat3(_In_ const Float3* pSource) noexcept;
		VECTOR VEC_CALLCONV LoadAFloat3(_In_ const AFloat3* pSource) noexcept;

		void VEC_CALLCONV StoreFloat3(_Out_ Float3* pDestination, _In_ A_VECTOR v) noexcept;
		void VEC_CALLCONV StoreAFloat3(_Out_ AFloat3* pDestination, _In_ A_VECTOR v) noexcept;
	}

	namespace VEC3
	{
		bool VEC_CALLCONV Equal(A_VECTOR V1, A_VECTOR V2) noexcept;
		uint32_t VEC_CALLCONV EqualR(A_VECTOR V1, A_VECTOR V2) noexcept;
		bool VEC_CALLCONV NearEqual(A_VECTOR V1, A_VECTOR V2, A_VECTOR epsilon) noexcept;
		bool VEC_CALLCONV NotEqual(A_VECTOR V1, A_VECTOR V2) noexcept;
		bool VEC_CALLCONV Greater(A_VECTOR V1, A_VECTOR V2) noexcept;
		uint32_t VEC_CALLCONV GreaterR(A_VECTOR V1, A_VECTOR V2) noexcept;
		bool VEC_CALLCONV GreaterOrEqual(A_VECTOR V1, A_VECTOR V2) noexcept;
		uint32_t VEC_CALLCONV GreaterOrEqualR(A_VECTOR V1, A_VECTOR V2) noexcept;
		bool VEC_CALLCONV Less(A_VECTOR V1, A_VECTOR V2) noexcept;
		bool VEC_CALLCONV LessOrEqual(A_VECTOR V1, A_VECTOR V2) noexcept;
		bool VEC_CALLCONV InBounds(A_VECTOR v, A_VECTOR bounds) noexcept;

		bool VEC_CALLCONV IsNaN(A_VECTOR v) noexcept;
		bool VEC_CALLCONV IsInfinite(A_VECTOR v) noexcept;

		VECTOR VEC_CALLCONV Dot(A_VECTOR V1, A_VECTOR V2) noexcept;
		VECTOR VEC_CALLCONV Cross(A_VECTOR V1, A_VECTOR V2) noexcept;
		VECTOR VEC_CALLCONV LengthSq(A_VECTOR v) noexcept;
		VECTOR VEC_CALLCONV ReciprocalLengthEst(A_VECTOR v) noexcept;
		VECTOR VEC_CALLCONV ReciprocalLength(A_VECTOR v) noexcept;
		VECTOR VEC_CALLCONV LengthEst(A_VECTOR v) noexcept;
		VECTOR VEC_CALLCONV Length(A_VECTOR v) noexcept;
		// Uses a reciprocal estimate and returns QNaN on zero and infinite vectors
		VECTOR VEC_CALLCONV NormalizeEst(A_VECTOR v) noexcept;
		VECTOR VEC_CALLCONV Normalize(A_VECTOR v) noexcept;
		VECTOR VEC_CALLCONV ClampLength(A_VECTOR v, float lengthMin, float lengthMax) noexcept;
		VECTOR VEC_CALLCONV ClampLengthV(A_VECTOR v, A_VECTOR lengthMin, A_VECTOR lengthMax) noexcept;
		VECTOR VEC_CALLCONV Reflect(A_VECTOR incident, A_VECTOR normal) noexcept;
		VECTOR VEC_CALLCONV Refract(A_VECTOR incident, A_VECTOR normal, float refractionIndex) noexcept;
		VECTOR VEC_CALLCONV RefractV(A_VECTOR incident, A_VECTOR normal, A_VECTOR refractionIndex) noexcept;
		VECTOR VEC_CALLCONV Orthogonal(A_VECTOR v) noexcept;
		VECTOR VEC_CALLCONV AngleBetweenNormalsEst(A_VECTOR N1, A_VECTOR N2) noexcept;
		VECTOR VEC_CALLCONV AngleBetweenNormals(A_VECTOR N1, A_VECTOR N2) noexcept;
		VECTOR VEC_CALLCONV AngleBetweenVectors(A_VECTOR V1, A_VECTOR V2) noexcept;
		VECTOR VEC_CALLCONV LinePointDistance(A_VECTOR linePoint1, A_VECTOR linePoint2, A_VECTOR point) noexcept;
		void VEC_CALLCONV ComponentsFromNormal(_Out_ VECTOR* pParallel, _Out_ VECTOR* pPerpendicular, _In_ A_VECTOR normal) noexcept;
		// Transform a vector using a rotation expressed as a unit quaternion
		VECTOR VEC_CALLCONV Rotate(A_VECTOR v, A_VECTOR rotationQuaternion) noexcept;
		// Transform a vector using the inverse of a rotation expressed as a unity quaternion
		VECTOR VEC_CALLCONV InverseRotate(A_VECTOR v, A_VECTOR rotationQuaternion) noexcept;
		VECTOR VEC_CALLCONV Transform(A_VECTOR v, A_MATRIX m) noexcept;
		Float4* VEC_CALLCONV TransformStream(_Out_writes_bytes_(sizeof(Float4) + outputStride * (vectorCount - 1)) Float4* pOutputStream, 
			_In_ size_t outputStride, 
			_In_reads_bytes_(sizeof(Float3) + inputStride * (vectorCount - 1)) const Float3* pInputStream, 
			_In_ size_t inputStride, _In_ size_t vectorCount, _In_ A_MATRIX m) noexcept;
		VECTOR VEC_CALLCONV TransformCoord(A_VECTOR v, A_MATRIX m) noexcept;
		Float3* VEC_CALLCONV TransformCoordStream(_Out_writes_bytes_(sizeof(Float3) + outputStride * (vectorCount - 1)) Float3* pOutputStream, 
			_In_ size_t outputStride, 
			_In_reads_bytes_(sizeof(Float3) + inputStride * (vectorCount - 1)) const Float3* pInputStream, 
			_In_ size_t inputStride, _In_ size_t vectorCount, _In_ A_MATRIX m) noexcept;
		VECTOR VEC_CALLCONV TransformNormal(A_VECTOR v, A_MATRIX m) noexcept;
		Float3* VEC_CALLCONV TransformNormalStream(_Out_writes_bytes_(sizeof(Float3) + outputStride * (vectorCount - 1)) Float3* pOutputStream, 
			_In_ size_t outputStride, 
			_In_reads_bytes_(sizeof(Float3) + inputStride * (vectorCount - 1)) const Float3* pInputStream, 
			_In_ size_t inputStride, _In_ size_t vectorCount, _In_ A_MATRIX m) noexcept;
		VECTOR VEC_CALLCONV Project(A_VECTOR v, 
			float viewportX, float viewportY, 
			float viewportWidth, float viewportHeight, 
			float viewportMinZ, float viewportMaxZ, 
			A_MATRIX projection, B_MATRIX view, B_MATRIX world);
		Float3* VEC_CALLCONV ProjectStream(_Out_writes_bytes_(sizeof(Float3) + outputStride * (vectorCount - 1)) Float3* pOutputStream, 
			_In_ size_t outputStride, 
			_In_reads_bytes_(sizeof(Float3) + inputStride * (vectorCount - 1)) const Float3* pInputStream, 
			_In_ size_t inputStride, _In_ size_t vectorCount, 
			_In_ float viewportX, _In_ float viewportY, 
			_In_ float viewportWidth, _In_ float viewportHeight, 
			_In_ float viewportMinZ, _In_ float viewportMaxZ, 
			_In_ A_MATRIX projection, _In_ B_MATRIX view, _In_ B_MATRIX world) noexcept;
		VECTOR VEC_CALLCONV Unproject(A_VECTOR v, 
			float viewportX, float viewportY, 
			float viewportWidth, float viewportHeight, 
			float viewportMinZ, float viewportMaxZ, 
			A_MATRIX projection, B_MATRIX view, B_MATRIX world);
		Float3* VEC_CALLCONV UnprojectStream(_Out_writes_bytes_(sizeof(Float3) + outputStride * (vectorCount - 1)) Float3* pOutputStream, 
			_In_ size_t outputStride, 
			_In_reads_bytes_(sizeof(Float3) + inputStride * (vectorCount - 1)) const Float3* pInputStream, 
			_In_ size_t inputStride, _In_ size_t vectorCount, 
			_In_ float viewportX, _In_ float viewportY, 
			_In_ float viewportWidth, _In_ float viewportHeight, 
			_In_ float viewportMinZ, _In_ float viewportMaxZ, 
			_In_ A_MATRIX projection, _In_ B_MATRIX view, _In_ B_MATRIX world) noexcept;
	}
}

#include <Float3.inl>

#ifdef _PREFAST_
#pragma prefast(pop)
#endif
#ifdef _MSC_VER
#pragma warning(pop)
#endif

#endif // !ULTREALITY_MATH_FLOAT3_H
