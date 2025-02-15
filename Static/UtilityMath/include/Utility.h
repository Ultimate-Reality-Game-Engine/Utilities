#ifndef ULTREALITY_MATH_UTILITY_H
#define ULTREALITY_MATH_UTILITY_H

#include <concepts>

namespace UltReality::Math
{
	/// <summary>
	/// Returns the minimum of the provided objects for totally ordered types (e.g. built-in numeric types, or user-defined types
	/// that satisfy the requirements for <, >, etc.)
	/// </summary>
	/// <typeparam name="T">The type of objects to compare</typeparam>
	/// <param name="a">First object to compare</param>
	/// <param name="b">Second object to compare</param>
	/// <returns>The minimum of the provided objects</returns>
	template<std::totally_ordered T>
	constexpr T Min(const T& a, const T& b) noexcept
	{
		return a < b ? a : b;
	}

	/// <summary>
	/// Max for totally ordered types (e.g. built-in numeric types, or user-defined types
	/// that satisfy the requirements for <, >, etc.)
	/// </summary>
	/// <typeparam name="T">The type of objects to compare</typeparam>
	/// <param name="a">First object to compare</param>
	/// <param name="b">Second object to compare</param>
	/// <returns>The maximum of the provided objects</returns>
	template<std::totally_ordered T>
	constexpr T Max(const T& a, const T& b) noexcept
	{
		return a > b ? a : b;
	}

	/// <summary>
	/// Clamp for totally ordered types (e.g. built-in numeric types, or user-defined types
	/// that satisfy the requirements for <, >, etc.)
	/// </summary>
	/// <typeparam name="T">The type of objects to compare</typeparam>
	/// <param name="x">Object to clamp to lower and upper thresholds</param>
	/// <param name="low">Lower threshold to clamp to if object value is below</param>
	/// <param name="high">Upper threshold to clamp to if object value is above</param>
	/// <returns>low if object is below low, high if object is above high, x if object is within thresholds</returns>
	template<std::totally_ordered T>
	constexpr T Clamp(const T& x, const T& low, const T& high) noexcept
	{
		return x < low ? low : (x > high ? high : x);
	}

	// Define a concept for types that can be linearly interpolated.
	// They must support subtraction, multiplication by a float, and addition.
	template<typename T>
	concept Lerpable = requires(T a, T b, float t)
	{
		{ a + (b - a) * t } -> std::convertible_to<T>;
	};

	/// <summary>
	/// Linear interpolation from a to b weighted by t
	/// </summary>
	/// <typeparam name="T">The type of object to linearly interpolate</typeparam>
	/// <param name="a">First end of interpolation</param>
	/// <param name="b">Second end of interpolation</param>
	/// <param name="t">Weight factor for interpolation</param>
	/// <returns>Weighted interpolation between points</returns>
	template<Lerpable T>
	constexpr T Lerp(const T& a, const T& b, float t) noexcept
	{
		return a + (b - a) * t;
	}

	/// <summary>
	/// Returns the polar angle of the point (x, y) in [0, 2PI)
	/// </summary>
	/// <param name="x">X coordinate of point</param>
	/// <param name="y">Y coordinate of point</param>
	constexpr float AngleFromXY(float x, float y) noexcept;

	bool ScalarNearEqual(float S1, float S2, float epsilon) noexcept;
	// Modulo the range of the given angle such that -PI <= angle <= PI
	float ScalarModAngle(float value) noexcept;

	float ScalarSine(float value) noexcept;
	float ScalarSineEst(float value) noexcept;

	float ScalarCos(float value) noexcept;
	float ScalarCosEst(float value) noexcept;

	float ScalarTan(float value) noexcept;
	float ScalarTanEst(float value) noexcept;

	void ScalarSineCos(_Out_ float* pSine, _Out_ float* pCos, float value) noexcept;
	void ScalarSineCosEst(_Out_ float* pSine, _Out_ float* pCos, float value) noexcept;

	float ScalarASine(float value) noexcept;
	float ScalarASineEst(float value) noexcept;

	float ScalarACos(float value) noexcept;
	float ScalarACosEst(float value) noexcept;

	float ScalarATan(float value) noexcept;
	float ScalarATanEst(float value) noexcept;

	constexpr float ConvertToRadians(float degrees) noexcept;
	
	constexpr float ConvertToDegrees(float radians) noexcept;
}

#include <Utility.inl>

#endif // !ULTREALITY_MATH_UTILITY_H
