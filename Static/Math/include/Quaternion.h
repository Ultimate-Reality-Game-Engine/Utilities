#ifndef ULTREALITY_MATH_QUATERNION_H
#define ULTREALITY_MATH_QUATERNION_H

#include <SSE2VectorConfig.h>
#include <Float3.h>
#include <Float4.h>
#include <MATRIX.h>

namespace UltReality::Math
{
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
}

#endif // !ULTREALITY_MATH_QUATERNION_H