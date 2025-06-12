#ifndef ULTREALITY_MATH_VECTOR_TYPES_H
#define ULTREALITY_MATH_VECTOR_TYPES_H

#include <SIMDVector.h>
#include <SIMDVectorTemplates.h>
#include <Utility.h>
#include <Random.h>

/*
	This file contains declarations for the following types, as well as 
	namespace interpretations of Float4/VECTOR4 as Plane and Quaternion.
	There implementation details are held in the corresponding .inl files
	They are included at the end of this file
*/
namespace UltReality::Math
{
	struct MATRIX;
	
	struct Int2;
	struct Int3;
	struct Int4;
	struct AInt2;
	struct AInt3;
	struct AInt4;

	struct UInt2;
	struct UInt3;
	struct UInt4;
	struct AUInt2;
	struct AUInt3;
	struct AUInt4;

	struct Float2;
	struct Float3;
	struct Float4;
	struct AFloat2;
	struct AFloat3;
	struct AFloat4;

	struct Float3x3;
	struct Float3x4;
	struct Float4x3;
	struct Float4x4;
	struct AFloat3x3;
	struct AFloat3x4;
	struct AFloat4x3;
	struct AFloat4x4;
}

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

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wgnu-anonymous-struct"
#pragma clang diagnostic ignored "-Wnested-anon-types"
#pragma clang diagnostic ignored "-Wunknown-warning-option"
#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"
#endif

namespace UltReality::Math
{
	struct MATRIX;

	// Define alias to be used for 1st matrix type argument. Passed in registers for x86, and vector call convention; by reference otherwise
#if (defined(_M_IX86) || _VECTORCALL_ || __i386__) && !defined(_NO_INTRINSICS_)
	typedef const MATRIX A_MATRIX;
#else
	typedef const MATRIX& A_MATRIX;
#endif

	// Define alias to be used for 2nd+ matrix type arguments. Passed by reference
	typedef const MATRIX& B_MATRIX;

#if defined(_SSE2_INTRINSICS_) && !defined(_NO_INTRINSICS_)
	ALIGNED_STRUCT(16) MATRIX
#else
	struct MATRIX
#endif
	{
#if defined(_NO_INTRINSICS_)
		union
		{
			VECTOR r[4];
			struct
			{
				float _00, _01, _02, _03;
				float _10, _11, _12, _13;
				float _20, _21, _22, _23;
				float _30, _31, _32, _33;
			};
			float m[4][4];
		};
#else
		VECTOR r[4];
#endif

		MATRIX() = default;

		MATRIX(const MATRIX&) = default;

#if defined(_MSC_VER) && (_MSC_FULL_VER < 191426431)
		MATRIX& operator=(const MATRIX& m) noexcept
		{
			r[0] = m.r[0];
			r[1] = m.r[1];
			r[2] = m.r[2];
			r[3] = m.r[3];

			return *this;
		}
#else
		MATRIX& operator=(const MATRIX&) = default;

		MATRIX(MATRIX&&) = default;
		MATRIX& operator=(MATRIX&&) = default;
#endif

		constexpr MATRIX(A_VECTOR r0, A_VECTOR r1, A_VECTOR r2, B_VECTOR r3) noexcept : r{ r0, r1, r2, r3 } {}

		MATRIX(
			float m00, float m01, float m02, float m03,
			float m10, float m11, float m12, float m13,
			float m20, float m21, float m22, float m23,
			float m30, float m31, float m32, float m33
		) noexcept;

		explicit MATRIX(_In_reads_(16) const float* pArray) noexcept;

#if defined(_NO_INTRINSICS_)
		float operator() (size_t row, size_t column) const noexcept { return m[row][column]; }
		float& operator() (size_t row, size_t column) noexcept { return m[row][column]; }

		const float* operator[](size_t row) const noexcept { return m[row]; }
		float* operator[](size_t row) noexcept { return m[row]; }
#endif

		MATRIX operator+() const noexcept { return *this; }
		MATRIX operator-() const noexcept;

		MATRIX& VEC_CALLCONV operator+=(A_MATRIX m) noexcept;
		MATRIX& VEC_CALLCONV operator-=(A_MATRIX m) noexcept;
		MATRIX& VEC_CALLCONV operator*=(A_MATRIX m) noexcept;
		MATRIX& operator*=(float s) noexcept;
		MATRIX& operator/=(float s) noexcept;

		MATRIX VEC_CALLCONV operator+(A_MATRIX m) const noexcept;
		MATRIX VEC_CALLCONV operator-(A_MATRIX m) const noexcept;
		MATRIX VEC_CALLCONV operator*(A_MATRIX m) const noexcept;
		MATRIX operator*(float s) const noexcept;
		MATRIX operator/(float s) const noexcept;

		friend MATRIX VEC_CALLCONV operator*(float s, A_MATRIX m) noexcept;

		static const MATRIX Zero;
		static const MATRIX Identity;
	};

	namespace Matrix
	{
		// Return true if any entry in the matrix is NaN
		bool VEC_CALLCONV IsNaN(A_MATRIX m) noexcept;
		// Return true if any entry in the matrix is +/- INF
		bool VEC_CALLCONV IsInfinite(A_MATRIX m) noexcept;
		// Return true if the MATRIX is equal to identity
		bool VEC_CALLCONV IsIdentity(A_MATRIX m) noexcept;

		// Perform a 4x4 matrix multiply
		MATRIX VEC_CALLCONV Multiply(A_MATRIX M1, B_MATRIX M2) noexcept;
		MATRIX VEC_CALLCONV MultiplyTranspose(A_MATRIX M1, B_MATRIX M2) noexcept;
		MATRIX VEC_CALLCONV Transpose(A_MATRIX m) noexcept;
		// Returns the inverse and the determinant of a 4x4 matrix
		MATRIX VEC_CALLCONV Inverse(_Out_opt_ VECTOR* pDeterminant, _In_ A_MATRIX m) noexcept;
		MATRIX VEC_CALLCONV InverseTranspose(A_MATRIX m) noexcept;
		MATRIX VEC_CALLCONV VectorTensorProduct(A_VECTOR V1, A_VECTOR V2) noexcept;
		VECTOR VEC_CALLCONV Determinant(A_MATRIX m) noexcept;

		_Success_(return)
			bool VEC_CALLCONV Decompose(_Out_ VECTOR* outScale, _Out_ VECTOR* outRotationQuaternion, _Out_ VECTOR* outTranslation, _In_ A_MATRIX m) noexcept;

		MATRIX VEC_CALLCONV Identity() noexcept;
		MATRIX VEC_CALLCONV Set(
			float m00, float m01, float m02, float m03,
			float m10, float m11, float m12, float m13,
			float m20, float m21, float m22, float m23,
			float m30, float m31, float m32, float m33
		) noexcept;
		MATRIX VEC_CALLCONV Translation(float offsetX, float offsetY, float offsetZ) noexcept;
		MATRIX VEC_CALLCONV TranslationFromVector(A_VECTOR offset) noexcept;
		MATRIX VEC_CALLCONV Scaling(float scaleX, float scaleY, float scaleZ) noexcept;
		MATRIX VEC_CALLCONV ScalingFromVector(A_VECTOR scale) noexcept;
		MATRIX VEC_CALLCONV RotationX(float angle) noexcept;
		MATRIX VEC_CALLCONV RotationY(float angle) noexcept;
		MATRIX VEC_CALLCONV RotationZ(float angle) noexcept;

		// Rotates about y-axis (yaw), then x-axis (pitch), then z-axis (roll)
		MATRIX VEC_CALLCONV RotationPitchYawRoll(float pitch, float yaw, float roll) noexcept;

		// Rotates about y-axis (yaw) from angles.y, then x-axis (pitch) from angles.x, then z-axis (roll) from angles.z
		MATRIX VEC_CALLCONV RotationPitchYawRoll(A_VECTOR angles) noexcept;

		MATRIX VEC_CALLCONV RotationNormal(A_VECTOR normalAxis, float angle) noexcept;
		MATRIX VEC_CALLCONV RotationAxis(A_VECTOR axis, float angle) noexcept;
		MATRIX VEC_CALLCONV RotationQuaternion(A_VECTOR quaternion) noexcept;
		MATRIX VEC_CALLCONV Transformation2D(A_VECTOR scalingOrigin, float scalingOrientation, A_VECTOR scaling,
			A_VECTOR rotationOrigin, float rotation, B_VECTOR translation) noexcept;
		MATRIX VEC_CALLCONV Transformation(A_VECTOR scalingOrigin, A_VECTOR scalingOrientationQuaternion, A_VECTOR scaling,
			B_VECTOR rotationOrigin, C_VECTOR rotationQuaternion, C_VECTOR translation) noexcept;
		MATRIX VEC_CALLCONV AffineTransformation2D(A_VECTOR scaling, A_VECTOR rotationOrigin, float rotation, A_VECTOR translation) noexcept;
		MATRIX VEC_CALLCONV AffineTransformation(A_VECTOR scaling, A_VECTOR rotationOrigin, A_VECTOR rotationQuaternion, B_VECTOR translation) noexcept;
		MATRIX VEC_CALLCONV Reflect(A_VECTOR reflectionPlane) noexcept;
		MATRIX VEC_CALLCONV Shadow(A_VECTOR shadowPlane, A_VECTOR lightPosition) noexcept;

		MATRIX VEC_CALLCONV LookAtLH(A_VECTOR eyePosition, A_VECTOR focusPoint, A_VECTOR upDirection) noexcept;
		MATRIX VEC_CALLCONV LookAtRH(A_VECTOR eyePosition, A_VECTOR focusPoint, A_VECTOR upDirection) noexcept;
		MATRIX VEC_CALLCONV LookToLH(A_VECTOR eyePosition, A_VECTOR eyeDirection, A_VECTOR upDirection) noexcept;
		MATRIX VEC_CALLCONV LookToRH(A_VECTOR eyePosition, A_VECTOR eyeDirection, A_VECTOR upDirection) noexcept;
		MATRIX VEC_CALLCONV PerspectiveLH(float viewWidth, float viewHeight, float nearZ, float farZ) noexcept;
		MATRIX VEC_CALLCONV PerspectiveRH(float viewWidth, float viewHeight, float nearZ, float farZ) noexcept;
		MATRIX VEC_CALLCONV PerspectiveFovLH(float fovAngleY, float aspectRatio, float nearZ, float farZ) noexcept;
		MATRIX VEC_CALLCONV PerspectiveFovRH(float fovAngleY, float aspectRatio, float nearZ, float farZ) noexcept;
		MATRIX VEC_CALLCONV PerspectiveOffCenterLH(float viewLeft, float viewRight, float viewBottom, float viewTop, float nearZ, float farZ) noexcept;
		MATRIX VEC_CALLCONV PerspectiveOffCenterRH(float viewLeft, float viewRight, float viewBottom, float viewTop, float nearZ, float farZ) noexcept;
		MATRIX VEC_CALLCONV OrthographicLH(float viewWidth, float viewHeight, float nearZ, float farZ) noexcept;
		MATRIX VEC_CALLCONV OrthographicRH(float viewWidth, float viewHeight, float nearZ, float farZ) noexcept;
		MATRIX VEC_CALLCONV OrthographicOffCenterLH(float viewLeft, float viewRight, float viewBottom, float viewTop, float nearZ, float farZ) noexcept;
		MATRIX VEC_CALLCONV OrthographicOffCenterRH(float viewLeft, float viewRight, float viewBottom, float viewTop, float nearZ, float farZ) noexcept;
	}

	/// <summary>
	/// 2D vector of 32-bit signed integer components
	/// </summary>
	struct Int2
	{
		int32_t x;
		int32_t y;

		Int2() = default;

		Int2(const Int2&) = default;
		Int2& operator=(const Int2&) = default;

		Int2(Int2&&) = default;
		Int2& operator=(Int2&&) = default;

		constexpr Int2(int32_t _x, int32_t _y) noexcept : x(_x), y(_y) {}
		explicit Int2(_In_reads_(2) const int32_t* pArray) noexcept : x(pArray[0]), y(pArray[1]) {}

		explicit Int2(_In_ A_VECTOR v) noexcept;
		Int2& VEC_CALLCONV operator=(_In_ A_VECTOR v) noexcept;

#if (__cplusplus >= 202002L)
		bool operator==(const Int2&) const = default;
		auto operator<=>(const Int2&) const = default;
#endif
		VECTOR VEC_CALLCONV Load() noexcept;
		void VEC_CALLCONV Store(_In_ A_VECTOR v) noexcept;
	};

	/// <summary>
	/// 2D vector of 32-bit signed integer components aligned on a 16 byte boundary
	/// </summary>
	ALIGNED_STRUCT(16) AInt2 : public Int2
	{
		using Int2::Int2;

		explicit AInt2(_In_ A_VECTOR v) noexcept;
		AInt2& VEC_CALLCONV operator=(_In_ A_VECTOR v) noexcept;

		VECTOR VEC_CALLCONV Load() noexcept;
		void VEC_CALLCONV Store(_In_ A_VECTOR v) noexcept;
	};

	namespace Vector
	{
		VECTOR VEC_CALLCONV LoadInt2(_In_ const Int2* pSource) noexcept;
		VECTOR VEC_CALLCONV LoadAInt2(_In_ const AInt2* pSource) noexcept;

		void VEC_CALLCONV StoreInt2(_Out_ Int2* pDestination, _In_ A_VECTOR v) noexcept;
		void VEC_CALLCONV StoreAInt2(_Out_ AInt2* pDestination, _In_ A_VECTOR v) noexcept;
	}

	namespace Vector2
	{
		bool VEC_CALLCONV EqualInt(A_VECTOR V1, A_VECTOR V2) noexcept;
		uint32_t VEC_CALLCONV EqualIntR(A_VECTOR V1, A_VECTOR V2) noexcept;
		bool VEC_CALLCONV NotEqualInt(A_VECTOR V1, A_VECTOR V2) noexcept;
	}

	/// <summary>
	/// 3D vector of 32-bit integer components
	/// </summary>
	struct Int3
	{
		int32_t x;
		int32_t y;
		int32_t z;

		Int3() = default;

		Int3(const Int3&) = default;
		Int3& operator=(const Int3&) = default;

		Int3(Int3&&) = default;
		Int3& operator=(Int3&&) = default;

		constexpr Int3(int32_t _x, int32_t _y, int32_t _z) noexcept : x(_x), y(_y), z(_z) {}
		explicit Int3(_In_reads_(3) const int32_t* pArray) noexcept : x(pArray[0]), y(pArray[1]), z(pArray[2]) {}

		explicit Int3(_In_ A_VECTOR v) noexcept;
		Int3& VEC_CALLCONV operator=(_In_ A_VECTOR v) noexcept;

#if (__cplusplus >= 202002L)
		bool operator==(const Int3&) const = default;
		auto operator<=>(const Int3&) const = default;
#endif

		VECTOR VEC_CALLCONV Load() noexcept;
		void VEC_CALLCONV Store(_In_ A_VECTOR v) noexcept;
	};

	/// <summary>
	/// 3D vector of 32-bit integer components aligned on a 16 byte boundary
	/// </summary>
	ALIGNED_STRUCT(16) AInt3 : public Int3
	{
		using Int3::Int3;

		explicit AInt3(_In_ A_VECTOR v) noexcept;
		AInt3& VEC_CALLCONV operator=(_In_ A_VECTOR v) noexcept;

		VECTOR VEC_CALLCONV Load() noexcept;
		void VEC_CALLCONV Store(_In_ A_VECTOR v) noexcept;
	};

	namespace Vector
	{
		VECTOR VEC_CALLCONV LoadInt3(_In_ const Int3* pSource) noexcept;
		VECTOR VEC_CALLCONV LoadAInt3(_In_ const AInt3* pSource) noexcept;

		void VEC_CALLCONV StoreInt3(_Out_ Int3* pDestination, _In_ A_VECTOR v) noexcept;
		void VEC_CALLCONV StoreAInt3(_Out_ AInt3* pDestination, _In_ A_VECTOR v) noexcept;
	}

	namespace Vector3
	{
		bool VEC_CALLCONV EqualInt(A_VECTOR V1, A_VECTOR V2) noexcept;
		uint32_t VEC_CALLCONV EqualIntR(A_VECTOR V1, A_VECTOR V2) noexcept;
		bool VEC_CALLCONV NotEqualInt(A_VECTOR V1, A_VECTOR V2) noexcept;
	}

	/// <summary>
	/// 4D vector of 32-bit integer components
	/// </summary>
	struct Int4
	{
		int32_t x;
		int32_t y;
		int32_t z;
		int32_t w;

		Int4() = default;

		Int4(const Int4&) = default;
		Int4& operator=(const Int4&) = default;

		Int4(Int4&&) = default;
		Int4& operator=(Int4&&) = default;

		constexpr Int4(int32_t _x, int32_t _y, int32_t _z, int32_t _w) : x(_x), y(_y), z(_z), w(_w) {}
		explicit Int4(_In_reads_(4) const int32_t* pArray) : x(pArray[0]), y(pArray[1]), z(pArray[2]), w(pArray[3]) {}

		explicit Int4(_In_ A_VECTOR v) noexcept;
		Int4& VEC_CALLCONV operator=(_In_ A_VECTOR v) noexcept;

#if (__cplusplus >= 202002L)
		bool operator==(const Int4&) const = default;
		auto operator<=>(const Int4&) const = default;
#endif

		VECTOR VEC_CALLCONV Load() noexcept;
		void VEC_CALLCONV Store(_In_ A_VECTOR v) noexcept;
	};

	/// <summary>
	/// 4D vector of 32-bit integer components aligned on a 16 byte boundary
	/// </summary>
	ALIGNED_STRUCT(16) AInt4 : public Int4
	{
		using Int4::Int4;

		explicit AInt4(_In_ A_VECTOR v) noexcept;
		AInt4& VEC_CALLCONV operator=(_In_ A_VECTOR v) noexcept;

		VECTOR VEC_CALLCONV Load() noexcept;
		void VEC_CALLCONV Store(_In_ A_VECTOR v) noexcept;
	};

	namespace Vector
	{
		VECTOR VEC_CALLCONV LoadInt4(_In_ const Int4* pSource) noexcept;
		VECTOR VEC_CALLCONV LoadAInt4(_In_ const AInt4* pSource) noexcept;

		void VEC_CALLCONV StoreInt4(_Out_ Int4* pDestination, _In_ A_VECTOR v) noexcept;
		void VEC_CALLCONV StoreAInt4(_Out_ AInt4* pDestination, _In_ A_VECTOR v) noexcept;
	}

	namespace Vector4
	{
		bool VEC_CALLCONV EqualInt(A_VECTOR V1, A_VECTOR V2) noexcept;
		uint32_t VEC_CALLCONV EqualIntR(A_VECTOR V1, A_VECTOR V2) noexcept;
		bool VEC_CALLCONV NotEqualInt(A_VECTOR V1, A_VECTOR V2) noexcept;
	}

	/// <summary>
	/// 2D vector of 32-bit unsigned integer components
	/// </summary>
	struct UInt2
	{
		uint32_t x;
		uint32_t y;

		UInt2() = default;

		UInt2(const UInt2&) = default;
		UInt2& operator=(const UInt2&) = default;

		UInt2(UInt2&&) = default;
		UInt2& operator=(UInt2&&) = default;

		constexpr UInt2(uint32_t _x, uint32_t _y) noexcept : x(_x), y(_y) {}
		explicit UInt2(_In_reads_(2) const uint32_t* pArray) noexcept : x(pArray[0]), y(pArray[1]) {}

		explicit UInt2(_In_ A_VECTOR v) noexcept;
		UInt2& VEC_CALLCONV operator=(_In_ A_VECTOR v) noexcept;

#if (__cplusplus >= 202002L)
		bool operator==(const UInt2&) const = default;
		auto operator<=>(const UInt2&) const = default;
#endif

		VECTOR VEC_CALLCONV Load() noexcept;
		void VEC_CALLCONV Store(_In_ A_VECTOR v) noexcept;
	};

	/// <summary>
	/// 2D vector of 32-bit unsigned integer components aligned on a 16 byte boundary
	/// </summary>
	ALIGNED_STRUCT(16) AUInt2 : public UInt2
	{
		using UInt2::UInt2;

		explicit AUInt2(_In_ A_VECTOR v) noexcept;
		AUInt2& VEC_CALLCONV operator=(_In_ A_VECTOR v) noexcept;

		VECTOR VEC_CALLCONV Load() noexcept;
		void VEC_CALLCONV Store(_In_ A_VECTOR v) noexcept;
	};

	namespace Vector
	{
		VECTOR VEC_CALLCONV LoadUInt2(_In_ const UInt2* pSource) noexcept;
		VECTOR VEC_CALLCONV LoadAUInt2(_In_ const AUInt2* pSource) noexcept;

		void VEC_CALLCONV StoreUInt2(_Out_ UInt2* pDestination, _In_ A_VECTOR v) noexcept;
		void VEC_CALLCONV StoreAUInt2(_Out_ AUInt2* pDestination, _In_ A_VECTOR v) noexcept;
	}

	/// <summary>
	/// 3D vector of 32-bit unsigned integer components
	/// </summary>
	struct UInt3
	{
		uint32_t x;
		uint32_t y;
		uint32_t z;

		UInt3() = default;

		UInt3(const UInt3&) = default;
		UInt3& operator=(const UInt3&) = default;

		UInt3(UInt3&&) = default;
		UInt3& operator=(UInt3&&) = default;

		constexpr UInt3(uint32_t _x, uint32_t _y, uint32_t _z) noexcept : x(_x), y(_y), z(_z) {}
		explicit UInt3(_In_reads_(3) const uint32_t* pArray) noexcept : x(pArray[0]), y(pArray[1]), z(pArray[2]) {}

		explicit UInt3(_In_ A_VECTOR v) noexcept;
		UInt3& VEC_CALLCONV operator=(_In_ A_VECTOR v) noexcept;

#if (__cplusplus >= 202002L)
		bool operator==(const UInt3&) const = default;
		auto operator<=>(const UInt3&) const = default;
#endif

		VECTOR VEC_CALLCONV Load() noexcept;
		void VEC_CALLCONV Store(_In_ A_VECTOR v) noexcept;
	};

	/// <summary>
	/// 3D vector of 32-bit unsigned integer components aligned on a 16 byte boundary
	/// </summary>
	ALIGNED_STRUCT(16) AUInt3 : public UInt3
	{
		using UInt3::UInt3;

		explicit AUInt3(_In_ A_VECTOR v) noexcept;
		AUInt3& VEC_CALLCONV operator=(_In_ A_VECTOR v) noexcept;

		VECTOR VEC_CALLCONV Load() noexcept;
		void VEC_CALLCONV Store(_In_ A_VECTOR v) noexcept;
	};

	namespace Vector
	{
		VECTOR VEC_CALLCONV LoadUInt3(_In_ const UInt3* pSource) noexcept;
		VECTOR VEC_CALLCONV LoadAUInt3(_In_ const AUInt3* pSource) noexcept;

		void VEC_CALLCONV StoreUInt3(_Out_ UInt3* pDestination, _In_ A_VECTOR v) noexcept;
		void VEC_CALLCONV StoreAUInt3(_Out_ AUInt3* pDestination, _In_ A_VECTOR v) noexcept;
	}

	/// <summary>
	/// 4D vector of 32-bit unsigned integer components
	/// </summary>
	struct UInt4
	{
		uint32_t x;
		uint32_t y;
		uint32_t z;
		uint32_t w;

		UInt4() = default;

		UInt4(const UInt4&) = default;
		UInt4& operator=(const UInt4&) = default;

		UInt4(UInt4&&) = default;
		UInt4& operator=(UInt4&&) = default;

		constexpr UInt4(uint32_t _x, uint32_t _y, uint32_t _z, uint32_t _w) noexcept : x(_x), y(_y), z(_z), w(_w) {}
		explicit UInt4(_In_reads_(4) const uint32_t* pArray) noexcept : x(pArray[0]), y(pArray[1]), z(pArray[2]), w(pArray[3]) {}

		explicit UInt4(_In_ A_VECTOR v) noexcept;
		UInt4& VEC_CALLCONV operator=(_In_ A_VECTOR v) noexcept;

#if (__cplusplus >= 202002L)
		bool operator==(const UInt4&) const = default;
		auto operator<=>(const UInt4&) const = default;
#endif

		VECTOR VEC_CALLCONV Load() noexcept;
		void VEC_CALLCONV Store(_In_ A_VECTOR v) noexcept;
	};

	/// <summary>
	/// 4D vector of 32-bit unsigned integer components aligned on a 16 byte boundary
	/// </summary>
	ALIGNED_STRUCT(16) AUInt4 : public UInt4
	{
		using UInt4::UInt4;

		explicit AUInt4(_In_ A_VECTOR v) noexcept;
		AUInt4& VEC_CALLCONV operator=(_In_ A_VECTOR v) noexcept;

		VECTOR VEC_CALLCONV Load() noexcept;
		void VEC_CALLCONV Store(_In_ A_VECTOR v) noexcept;
	};

	namespace Vector
	{
		VECTOR VEC_CALLCONV LoadUInt4(_In_ const UInt4* pSource) noexcept;
		VECTOR VEC_CALLCONV LoadAUInt4(_In_ const AUInt4* pSource) noexcept;

		void VEC_CALLCONV StoreUInt4(_Out_ UInt4* pDestination, _In_ A_VECTOR v) noexcept;
		void VEC_CALLCONV StoreAUInt4(_Out_ AUInt4* pDestination, _In_ A_VECTOR v) noexcept;
	}

	/// <summary>
	/// 2D vector of 32-bit floating point components
	/// </summary>
	struct Float2
	{
		float x;
		float y;

		Float2() = default;

		Float2(const Float2&) = default;
		Float2& operator=(const Float2&) = default;

		Float2(Float2&&) = default;
		Float2& operator=(Float2&&) = default;

		constexpr Float2(float _x, float _y) noexcept : x(_x), y(_y) {}
		explicit Float2(_In_reads_(2) const float* pArray) noexcept : x(pArray[0]), y(pArray[1]) {}

		explicit Float2(_In_ A_VECTOR v) noexcept;
		Float2& VEC_CALLCONV operator=(_In_ A_VECTOR v) noexcept;

#if(__cplusplus >= 202002L)
		bool operator==(const Float2&) const = default;
		auto operator<=>(const Float2&) const = default;
#endif

		VECTOR VEC_CALLCONV Load() noexcept;
		void VEC_CALLCONV Store(_In_ A_VECTOR v) noexcept;
	};

	/// <summary>
	/// 2D vector of 32-bit floating point components aligned on a 16 byte boundary
	/// </summary>
	ALIGNED_STRUCT(16) AFloat2 : public Float2
	{
		using Float2::Float2;

		explicit AFloat2(_In_ A_VECTOR v) noexcept;
		AFloat2& VEC_CALLCONV operator=(_In_ A_VECTOR v) noexcept;

		VECTOR VEC_CALLCONV Load() noexcept;
		void VEC_CALLCONV Store(_In_ A_VECTOR v) noexcept;
	};

	namespace Vector
	{
		VECTOR VEC_CALLCONV LoadFloat2(_In_ const Float2* pSource) noexcept;
		VECTOR VEC_CALLCONV LoadAFloat2(_In_ const AFloat2* pSource) noexcept;

		void VEC_CALLCONV StoreFloat2(_Out_ Float2* pDestination, _In_ A_VECTOR v) noexcept;
		void VEC_CALLCONV StoreAFloat2(_Out_ AFloat2* pDestination, _In_ A_VECTOR v) noexcept;
	}

	namespace Vector2
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
		// NormalizeEst uses a reciprocal estimate and returns QNaN on zero and infinite vectors
		VECTOR VEC_CALLCONV NormalizeEst(A_VECTOR v) noexcept;
		VECTOR VEC_CALLCONV Normalize(A_VECTOR v) noexcept;
		VECTOR VEC_CALLCONV ClampLength(A_VECTOR v, float lengthMin, float lengthMax) noexcept;
		VECTOR VEC_CALLCONV ClampLengthV(A_VECTOR v, A_VECTOR lengthMin, A_VECTOR lengthMax) noexcept;
		VECTOR VEC_CALLCONV Reflect(A_VECTOR incident, A_VECTOR normal) noexcept;
		VECTOR VEC_CALLCONV Refract(A_VECTOR incident, A_VECTOR normal, float refractionIndex) noexcept;
		// Return the refraction of of a 2D vector
		VECTOR VEC_CALLCONV RefractV(A_VECTOR incident, A_VECTOR normal, A_VECTOR refractionIndex) noexcept;
		VECTOR VEC_CALLCONV Orthogonal(A_VECTOR v) noexcept;
		VECTOR VEC_CALLCONV AngleBetweenNormalsEst(A_VECTOR N1, A_VECTOR N2) noexcept;
		VECTOR VEC_CALLCONV AngleBetweenNormals(A_VECTOR N1, A_VECTOR N2) noexcept;
		VECTOR VEC_CALLCONV AngleBetweenVectors(A_VECTOR V1, A_VECTOR V2) noexcept;
		VECTOR VEC_CALLCONV LinePointDistance(A_VECTOR linePoint1, A_VECTOR linePoint2, A_VECTOR point) noexcept;
		VECTOR VEC_CALLCONV IntersectLine(A_VECTOR line1Point1, A_VECTOR line1Point2, A_VECTOR line2Point1, B_VECTOR line2Point2) noexcept;
		VECTOR VEC_CALLCONV Transform(A_VECTOR v, A_MATRIX m) noexcept;
		Float4* VEC_CALLCONV TransformStream(_Out_writes_bytes_(sizeof(Float4) + outputStride * (vectorCount - 1)) Float4* pOutputStream,
			_In_ size_t outputStride,
			_In_reads_bytes_(sizeof(Float2) + inputStride * (vectorCount - 1)) const Float2* pInputStream,
			_In_ size_t inputStride, _In_ size_t vectorCount, _In_ A_MATRIX m) noexcept;
		VECTOR VEC_CALLCONV TransformCoord(A_VECTOR v, A_VECTOR m) noexcept;
		Float2* VEC_CALLCONV TransformCoordStream(_Out_writes_bytes_(sizeof(Float2) + outputStream * (vectorCount - 1)) Float2* pOutputStream,
			_In_ size_t outputStride,
			_In_reads_bytes_(sizeof(Float2) + inputStride * (vectorCount - 1)) const Float2* pInputStream,
			_In_ size_t inputStride, _In_ size_t vectorCount, _In_ A_MATRIX m) noexcept;
		VECTOR VEC_CALLCONV TransformNormal(A_VECTOR v, A_MATRIX m) noexcept;
		Float2* VEC_CALLCONV TransformNormalStream(_Out_writes_bytes_(sizeof(Float2) + outputStream * (vectorCount - 1)) Float2* pOutputStream,
			_In_ size_t outputStream,
			_In_reads_bytes_(sizeof(Float2) + inputStride * (vectorCount - 1)) const Float2* pInputStream,
			_In_ size_t inputStride, _In_ size_t vectorCount, _In_ A_MATRIX m) noexcept;
	}

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

	namespace Vector
	{
		VECTOR VEC_CALLCONV LoadFloat3(_In_ const Float3* pSource) noexcept;
		VECTOR VEC_CALLCONV LoadAFloat3(_In_ const AFloat3* pSource) noexcept;

		void VEC_CALLCONV StoreFloat3(_Out_ Float3* pDestination, _In_ A_VECTOR v) noexcept;
		void VEC_CALLCONV StoreAFloat3(_Out_ AFloat3* pDestination, _In_ A_VECTOR v) noexcept;
	}

	namespace Vector3
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
			A_MATRIX projection, B_MATRIX view, B_MATRIX world) noexcept;
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
			A_MATRIX projection, B_MATRIX view, B_MATRIX world) noexcept;
		Float3* VEC_CALLCONV UnprojectStream(_Out_writes_bytes_(sizeof(Float3) + outputStride * (vectorCount - 1)) Float3* pOutputStream,
			_In_ size_t outputStride,
			_In_reads_bytes_(sizeof(Float3) + inputStride * (vectorCount - 1)) const Float3* pInputStream,
			_In_ size_t inputStride, _In_ size_t vectorCount,
			_In_ float viewportX, _In_ float viewportY,
			_In_ float viewportWidth, _In_ float viewportHeight,
			_In_ float viewportMinZ, _In_ float viewportMaxZ,
			_In_ A_MATRIX projection, _In_ B_MATRIX view, _In_ B_MATRIX world) noexcept;
		VECTOR VEC_CALLCONV RandUnit() noexcept;
		VECTOR VEC_CALLCONV RandHemisphereUnit(A_VECTOR n) noexcept;
	}

	/// <summary>
	/// 4D vector of 32-bit floating point components
	/// </summary>
	struct Float4
	{
		float x;
		float y;
		float z;
		float w;

		Float4() = default;

		Float4(const Float4&) = default;
		Float4& operator=(const Float4&) = default;

		Float4(Float4&&) = default;
		Float4& operator=(Float4&&) = default;

		constexpr Float4(float _x, float _y, float _z, float _w) : x(_x), y(_y), z(_z), w(_w) {}
		explicit Float4(_In_reads_(4) const float* pArray) noexcept : x(pArray[0]), y(pArray[1]), z(pArray[2]), w(pArray[3]) {}

		explicit Float4(_In_ A_VECTOR v) noexcept;
		Float4& VEC_CALLCONV operator=(_In_ A_VECTOR v) noexcept;

#if (__cplusplus >= 202002L)
		bool operator==(const Float4&) const = default;
		auto operator<=>(const Float4&) const = default;
#endif

		VECTOR VEC_CALLCONV Load() noexcept;
		void VEC_CALLCONV Store(_In_ A_VECTOR v) noexcept;
	};

	/// <summary>
	/// 4D vector of 32-bit floating point components aligned on a 16 byte boundary
	/// </summary>
	ALIGNED_STRUCT(16) AFloat4 : public Float4
	{
		using Float4::Float4;

		explicit AFloat4(_In_ A_VECTOR v) noexcept;
		AFloat4& VEC_CALLCONV operator=(_In_ A_VECTOR v) noexcept;

		VECTOR VEC_CALLCONV Load() noexcept;
		void VEC_CALLCONV Store(_In_ A_VECTOR v) noexcept;
	};

	namespace Vector
	{
		VECTOR VEC_CALLCONV LoadFloat4(_In_ const Float4* pSource) noexcept;
		VECTOR VEC_CALLCONV LoadAFloat4(_In_ const AFloat4* pSource) noexcept;

		void VEC_CALLCONV StoreFloat4(_Out_ Float4* pDestination, _In_ A_VECTOR v) noexcept;
		void VEC_CALLCONV StoreAFloat4(_Out_ AFloat4* pDestination, _In_ A_VECTOR v) noexcept;
	}

	namespace Vector4
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
		VECTOR VEC_CALLCONV Cross(A_VECTOR V1, A_VECTOR V2, A_VECTOR V3) noexcept;
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
		VECTOR VEC_CALLCONV Transform(A_VECTOR v, A_MATRIX m) noexcept;
		Float4* VEC_CALLCONV TransformStream(_Out_writes_bytes_(sizeof(Float4) + outputStride * (vectorCount - 1)) Float4* pOutputStream,
			_In_ size_t outputStride,
			_In_reads_bytes_(sizeof(Float4) inputStride* (vectorCount - 1)) const Float4* pInputStream,
			_In_ size_t inputStride, _In_ size_t vectorCount, _In_ A_MATRIX m) noexcept;
	}

	namespace Plane
	{
		bool VEC_CALLCONV Equal(A_VECTOR P1, A_VECTOR P2) noexcept;
		bool VEC_CALLCONV NearEqual(A_VECTOR P1, A_VECTOR P2, A_VECTOR epsilon) noexcept;
		bool VEC_CALLCONV NotEqual(A_VECTOR P1, A_VECTOR P2) noexcept;

		bool VEC_CALLCONV IsNaN(A_VECTOR p) noexcept;
		bool VEC_CALLCONV IsInfinite(A_VECTOR p) noexcept;

		VECTOR VEC_CALLCONV Dot(A_VECTOR p, A_VECTOR v) noexcept;
		VECTOR VEC_CALLCONV DotCoord(A_VECTOR p, A_VECTOR v) noexcept;
		VECTOR VEC_CALLCONV DotNormal(A_VECTOR p, A_VECTOR v) noexcept;
		// Uses a reciprocal estimate and return QNaN on zero and infinite vectors
		VECTOR VEC_CALLCONV NormalizeEst(A_VECTOR p) noexcept;
		VECTOR VEC_CALLCONV Normalize(A_VECTOR p) noexcept;
		VECTOR VEC_CALLCONV IntersectLine(A_VECTOR p, A_VECTOR linePoint1, A_VECTOR linePoint2) noexcept;
		void VEC_CALLCONV IntersectPlane(_Out_ VECTOR* pLinePoint1, _Out_ VECTOR* pLinePoint2, _In_ A_VECTOR P1, _In_ A_VECTOR P2) noexcept;

		// Transforms a plane given an inverse transpose matrix
		VECTOR VEC_CALLCONV Transform(A_VECTOR p, A_MATRIX ITM) noexcept;

		// Transforms an array of planes given an inverse transpose matrix
		Float4* VEC_CALLCONV TransformStream(_Out_writes_bytes_(sizeof(Float4) + outputStride * (planeCount - 1)) Float4* pOutputStream,
			_In_ size_t outputStride,
			_In_reads_bytes_(sizeof(Float4) + inputStride * (planeCount - 1)) const Float4* pInputStream,
			_In_ size_t inputStride, _In_ size_t planeCount, _In_ A_MATRIX ITM) noexcept;

		VECTOR VEC_CALLCONV FromPointNormal(A_VECTOR point, A_VECTOR normal) noexcept;
		VECTOR VEC_CALLCONV FromPoints(A_VECTOR point1, A_VECTOR point2, A_VECTOR point3) noexcept;
	} // namespace Plane

	namespace Quaternion
	{
		bool VEC_CALLCONV Equal(A_VECTOR Q1, A_VECTOR Q2) noexcept;
		bool VEC_CALLCONV NotEqual(A_VECTOR Q1, A_VECTOR Q2) noexcept;

		bool VEC_CALLCONV IsNaN(A_VECTOR q) noexcept;
		bool VEC_CALLCONV IsInfinite(A_VECTOR q) noexcept;
		bool VEC_CALLCONV IsIdentity(A_VECTOR q) noexcept;

		VECTOR VEC_CALLCONV Dot(A_VECTOR Q1, A_VECTOR Q2) noexcept;
		VECTOR VEC_CALLCONV Multiply(A_VECTOR Q1, A_VECTOR Q2) noexcept;
		VECTOR VEC_CALLCONV LengthSq(A_VECTOR q) noexcept;
		VECTOR VEC_CALLCONV ReciprocalLength(A_VECTOR q) noexcept;
		VECTOR VEC_CALLCONV Length(A_VECTOR q) noexcept;
		VECTOR VEC_CALLCONV NormalizeEst(A_VECTOR q) noexcept;
		VECTOR VEC_CALLCONV Normalize(A_VECTOR q) noexcept;
		VECTOR VEC_CALLCONV Conjugate(A_VECTOR q) noexcept;
		VECTOR VEC_CALLCONV Inverse(A_VECTOR q) noexcept;
		VECTOR VEC_CALLCONV Ln(A_VECTOR q) noexcept;
		VECTOR VEC_CALLCONV Exp(A_VECTOR q) noexcept;
		VECTOR VEC_CALLCONV Slerp(A_VECTOR Q1, A_VECTOR Q2, float t) noexcept;
		VECTOR VEC_CALLCONV SlerpV(A_VECTOR Q1, A_VECTOR Q2, A_VECTOR T) noexcept;
		VECTOR VEC_CALLCONV Squad(A_VECTOR Q1, A_VECTOR Q2, A_VECTOR Q3, B_VECTOR Q4, float t) noexcept;
		VECTOR VEC_CALLCONV SquadV(A_VECTOR Q1, A_VECTOR Q2, A_VECTOR Q3, B_VECTOR Q4, C_VECTOR T) noexcept;
		void VEC_CALLCONV SquadSetup(_Out_ VECTOR* pA, _Out_ VECTOR* pB, _Out_ VECTOR* pC,
			_In_ A_VECTOR Q1, _In_ A_VECTOR Q2, _In_ A_VECTOR Q3, _In_ B_VECTOR Q4) noexcept;
		VECTOR VEC_CALLCONV BaryCentric(A_VECTOR Q1, A_VECTOR Q2, A_VECTOR Q3, float f, float g) noexcept;
		VECTOR VEC_CALLCONV BaryCentricV(A_VECTOR Q1, A_VECTOR Q2, A_VECTOR Q3, B_VECTOR F, C_VECTOR G) noexcept;

		VECTOR VEC_CALLCONV Identity() noexcept;

		// Rotates about y-axis (yaw), then x_axis (pitch), then z-axis (roll)
		VECTOR VEC_CALLCONV RotationPitchYawRoll(float pitch, float yaw, float roll) noexcept;

		// Rotates about y-axis (yaw) from angles.y, then x-axis (pitch) from angles.x, then z-axis (roll) from angles.z
		VECTOR VEC_CALLCONV RotationPitchYawRollFromVector(A_VECTOR angles) noexcept;

		VECTOR VEC_CALLCONV RotationNormal(A_VECTOR normalAxis, float angle) noexcept;
		VECTOR VEC_CALLCONV RotationAxis(A_VECTOR axis, float angle) noexcept;
		VECTOR VEC_CALLCONV RotationMatrix(A_MATRIX m) noexcept;

		void VEC_CALLCONV ToAxisAngle(_Out_ VECTOR* pAxis, _Out_ float* pAngle, _In_ A_VECTOR q) noexcept;
	}

	/// <summary>
	/// 3x3 matrix of 32-bit floating point components
	/// </summary>
	struct Float3x3
	{
		union
		{
			struct
			{
				float _00, _01, _02;
				float _10, _11, _12;
				float _20, _21, _22;
			};
			float m[3][3];
		};

		Float3x3() = default;

		Float3x3(const Float3x3&) = default;
		Float3x3& operator=(const Float3x3&) = default;

		Float3x3(Float3x3&&) = default;
		Float3x3& operator=(Float3x3&&) = default;

		constexpr Float3x3(
			float m00, float m01, float m02,
			float m10, float m11, float m12,
			float m20, float m21, float m22
		) noexcept :
			_00(m00), _01(m01), _02(m02),
			_10(m10), _11(m11), _12(m12),
			_20(m20), _21(m21), _22(m22)
		{
		}

		explicit Float3x3(_In_reads_(9) const float* pArray) noexcept;

		explicit Float3x3(_In_ A_MATRIX m) noexcept;
		Float3x3& VEC_CALLCONV operator=(_In_ A_MATRIX m) noexcept;

		float operator() (size_t row, size_t column) const noexcept { return m[row][column]; }
		float& operator() (size_t row, size_t column) noexcept { return m[row][column]; }

		const float* operator[](size_t row) const noexcept { return m[row]; }
		float* operator[](size_t row) noexcept { return m[row]; }

#if (__cplusplus >= 202002L)
		bool operator==(const Float3x3&) const = default;
		auto operator<=>(const Float3x3&) const = default;
#endif

		MATRIX VEC_CALLCONV Load() noexcept;
		void VEC_CALLCONV Store(_In_ A_MATRIX m) noexcept;
	};

	/// <summary>
	/// 3x3 matrix of 32-bit floating point components aligned on a 16 byte boundary
	/// </summary>
	ALIGNED_STRUCT(16) AFloat3x3 : public Float3x3
	{
		using Float3x3::Float3x3;

		explicit AFloat3x3(_In_ A_MATRIX m) noexcept;
		AFloat3x3& VEC_CALLCONV operator=(_In_ A_MATRIX m) noexcept;

		MATRIX VEC_CALLCONV Load() noexcept;
		void VEC_CALLCONV Store(_In_ A_MATRIX m) noexcept;
	};

	namespace Matrix
	{
		MATRIX VEC_CALLCONV LoadFloat3x3(_In_ const Float3x3* pSource) noexcept;
		MATRIX VEC_CALLCONV LoadAFloat3x3(_In_ const AFloat3x3* pSource) noexcept;

		void VEC_CALLCONV StoreFloat3x3(_Out_ Float3x3* pDestination, _In_ A_MATRIX m) noexcept;
		void VEC_CALLCONV StoreAFloat3x3(_Out_ AFloat3x3* pDestination, _In_ A_MATRIX m) noexcept;
	}

	/// <summary>
	/// 3x4 matrix of 32-bit floating point components
	/// </summary>
	struct Float3x4
	{
		union
		{
			struct
			{
				float _00, _01, _02, _03;
				float _10, _11, _12, _13;
				float _20, _21, _22, _23;
			};
			float m[3][4];
		};

		Float3x4() = default;

		Float3x4(const Float3x4&) = default;
		Float3x4& operator=(const Float3x4&) = default;

		Float3x4(Float3x4&&) = default;
		Float3x4& operator=(Float3x4&&) = default;

		constexpr Float3x4(
			float m00, float m01, float m02, float m03,
			float m10, float m11, float m12, float m13,
			float m20, float m21, float m22, float m23
		) noexcept :
			_00(m00), _01(m01), _02(m02), _03(m03),
			_10(m10), _11(m11), _12(m12), _13(m13),
			_20(m20), _21(m21), _22(m22), _23(m23)
		{
		}

		explicit Float3x4(_In_reads_(12) const float* pArray) noexcept;

		explicit Float3x4(_In_ A_MATRIX m) noexcept;
		Float3x4& VEC_CALLCONV operator=(_In_ A_MATRIX m) noexcept;

		float operator() (size_t row, size_t column) const noexcept { return m[row][column]; }
		float& operator() (size_t row, size_t column) noexcept { return m[row][column]; }

		const float* operator[](size_t row) const noexcept { return m[row]; }
		float* operator[](size_t row) noexcept { return m[row]; }

#if (__cplusplus >= 202002L)
		bool operator==(const Float3x4&) const = default;
		auto operator<=>(const Float3x4&) const = default;
#endif

		MATRIX VEC_CALLCONV Load() noexcept;
		void VEC_CALLCONV Store(_In_ A_MATRIX m) noexcept;
	};

	/// <summary>
	/// 3x4 matrix of 32-bit floating point components aligned on a 16 byte boundary
	/// </summary>
	ALIGNED_STRUCT(16) AFloat3x4 : public Float3x4
	{
		using Float3x4::Float3x4;

		explicit AFloat3x4(_In_ A_MATRIX m) noexcept;
		AFloat3x4& VEC_CALLCONV operator=(_In_ A_MATRIX m) noexcept;

		MATRIX VEC_CALLCONV Load() noexcept;
		void VEC_CALLCONV Store(_In_ A_MATRIX m) noexcept;
	};

	namespace Matrix
	{
		MATRIX VEC_CALLCONV LoadFloat3x4(_In_ const Float3x4* pSource) noexcept;
		MATRIX VEC_CALLCONV LoadAFloat3x4(_In_ const AFloat3x4* pSource) noexcept;

		void VEC_CALLCONV StoreFloat3x4(_Out_ Float3x4* pDestination, _In_ A_MATRIX m) noexcept;
		void VEC_CALLCONV StoreAFloat3x4(_Out_ AFloat3x4* pDestination, _In_ A_MATRIX m) noexcept;
	}

	/// <summary>
	/// 4x3 matrix of 32-bit floating point components
	/// </summary>
	struct Float4x3
	{
		union
		{
			struct
			{
				float _00, _01, _02;
				float _10, _11, _12;
				float _20, _21, _22;
				float _30, _31, _32;
			};
			float m[4][3];
		};

		Float4x3() = default;

		Float4x3(const Float4x3&) = default;
		Float4x3& operator=(const Float4x3&) = default;

		Float4x3(Float4x3&&) = default;
		Float4x3& operator=(Float4x3&&) = default;

		constexpr Float4x3(
			float m00, float m01, float m02,
			float m10, float m11, float m12,
			float m20, float m21, float m22,
			float m30, float m31, float m32
		) noexcept :
			_00(m00), _01(m01), _02(m02),
			_10(m10), _11(m11), _12(m12),
			_20(m20), _21(m21), _22(m22),
			_30(m30), _31(m31), _32(m32)
		{
		}

		explicit Float4x3(_In_reads_(12) const float* pArray) noexcept;

		explicit Float4x3(_In_ A_MATRIX m) noexcept;
		Float4x3& VEC_CALLCONV operator=(_In_ A_MATRIX m) noexcept;

		float operator() (size_t row, size_t column) const noexcept { return m[row][column]; }
		float& operator() (size_t row, size_t column) noexcept { return m[row][column]; }

		const float* operator[](size_t row) const noexcept { return m[row]; }
		float* operator[](size_t row) noexcept { return m[row]; }

#if (__cplusplus >= 202002L)
		bool operator==(const Float4x3&) const = default;
		auto operator<=>(const Float4x3&) const = default;
#endif

		MATRIX VEC_CALLCONV Load() noexcept;
		void VEC_CALLCONV Store(_In_ A_MATRIX m) noexcept;
	};

	/// <summary>
	/// 4x3 matrix of 32-bit floating point components aligned on a 16 byte boundary
	/// </summary>
	ALIGNED_STRUCT(16) AFloat4x3 : public Float4x3
	{
		using Float4x3::Float4x3;

		explicit AFloat4x3(_In_ A_MATRIX m) noexcept;
		AFloat4x3& VEC_CALLCONV operator=(_In_ A_MATRIX m) noexcept;

		MATRIX VEC_CALLCONV Load() noexcept;
		void VEC_CALLCONV Store(_In_ A_MATRIX m) noexcept;
	};

	namespace Matrix
	{
		MATRIX VEC_CALLCONV LoadFloat4x3(_In_ const Float4x3* pSource) noexcept;
		MATRIX VEC_CALLCONV LoadAFloat4x3(_In_ const AFloat4x3* pSource) noexcept;

		void VEC_CALLCONV StoreFloat4x3(_Out_ Float4x3* pDestination, _In_ A_MATRIX m) noexcept;
		void VEC_CALLCONV StoreAFloat4x3(_Out_ AFloat4x3* pDestination, _In_ A_MATRIX m) noexcept;
	}

	/// <summary>
	/// 4x4 matrix of 32-bit floating point components
	/// </summary>
	struct Float4x4
	{
		union
		{
			struct
			{
				float _00, _01, _02, _03;
				float _10, _11, _12, _13;
				float _20, _21, _22, _23;
				float _30, _31, _32, _33;
			};
			float m[4][4];
		};

		Float4x4() = default;

		Float4x4(const Float4x4&) = default;
		Float4x4& operator=(const Float4x4&) = default;

		Float4x4(Float4x4&&) = default;
		Float4x4& operator=(Float4x4&&) = default;

		constexpr Float4x4(
			float m00, float m01, float m02, float m03,
			float m10, float m11, float m12, float m13,
			float m20, float m21, float m22, float m23,
			float m30, float m31, float m32, float m33
		) noexcept :
			_00(m00), _01(m01), _02(m02), _03(m03),
			_10(m10), _11(m11), _12(m12), _13(m13),
			_20(m20), _21(m21), _22(m22), _23(m23),
			_30(m30), _31(m31), _32(m32), _33(m33)
		{
		}

		explicit Float4x4(_In_reads_(16) const float* pArray) noexcept;

		explicit Float4x4(_In_ A_MATRIX m) noexcept;
		Float4x4& VEC_CALLCONV operator=(_In_ A_MATRIX m) noexcept;

		float operator() (size_t row, size_t column) const noexcept { return m[row][column]; }
		float& operator() (size_t row, size_t column) noexcept { return m[row][column]; }

		const float* operator[](size_t row) const noexcept { return m[row]; }
		float* operator[](size_t row) noexcept { return m[row]; }

#if (__cplusplus >= 202002L)
		bool operator==(const Float4x4&) const = default;
		auto operator<=>(const Float4x4&) const = default;
#endif

		MATRIX VEC_CALLCONV Load() noexcept;
		void VEC_CALLCONV Store(_In_ A_MATRIX m) noexcept;

		static const Float4x4 Identity;
	};

	/// <summary>
	/// 4x4 matrix of 32-bit floating point components aligned on a 16 byte boundary
	/// </summary>
	ALIGNED_STRUCT(16) AFloat4x4 : public Float4x4
	{
		using Float4x4::Float4x4;

		explicit AFloat4x4(_In_ A_MATRIX m) noexcept;
		AFloat4x4& VEC_CALLCONV operator=(_In_ A_MATRIX m) noexcept;

		MATRIX VEC_CALLCONV Load() noexcept;
		void VEC_CALLCONV Store(_In_ A_MATRIX m) noexcept;
	};

	namespace Matrix
	{
		MATRIX VEC_CALLCONV LoadFloat4x4(_In_ const Float4x4* pSource) noexcept;
		MATRIX VEC_CALLCONV LoadAFloat4x4(_In_ const AFloat4x4* pSource) noexcept;

		void VEC_CALLCONV StoreFloat4x4(_Out_ Float4x4* pDestination, _In_ A_MATRIX m) noexcept;
		void VEC_CALLCONV StoreAFloat4x4(_Out_ AFloat4x4* pDestination, _In_ A_MATRIX m) noexcept;
	}
}

#include <MATRIX.inl>

#include <Int2.inl>
#include <Int3.inl>
#include <Int4.inl>

#include <UInt2.inl>
#include <UInt3.inl>
#include <UInt4.inl>

#include <Float2.inl>
#include <Float3.inl>
#include <Float4.inl>
#include <Plane.inl>
#include <Quaternion.inl>

#include <Float3x3.inl>
#include <Float3x4.inl>
#include <Float4x3.inl>
#include <Float4x4.inl>

#ifdef __clang__
#pragma clang diagnostic pop
#endif
#ifdef _PREFAST_
#pragma prefast(pop)
#endif
#ifdef _MSC_VER
#pragma warning(pop)
#endif

#endif // !ULTREALITY_MATH_VECTOR_TYPES_H
