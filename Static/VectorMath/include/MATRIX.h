#ifndef ULTREALITY_MATH_SSE2_MATRIX_H
#define ULTREALITY_MATH_SSE2_MATRIX_H

#include <SIMDVectorConfig.h>
#include <Utility.h>
#include <Float3.h>
#include <Float4.h>
#include <Quaternion.h>
#include <Plane.h>

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

		const float* operator[]size_t row) const noexcept { return m[row]; }
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
}

#include <MATRIX.inl>

#endif // !ULTREALITY_MATH_SSE2_MATRIX_H
