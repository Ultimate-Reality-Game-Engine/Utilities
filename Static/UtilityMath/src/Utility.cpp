#include <Utility.h>

#include <cmath>
#include <Constants.h>

namespace UltReality::Math
{
	constexpr float AngleFromXY(float x, float y) noexcept
	{
		float theta = 0.0f;

		// Quadrant I or IV
		if (x >= 0.0f)
		{
			// If x == 0, then atanf(y / x) = +pi/2 if y > 0
			//				   atanf(y / x) = -pi/2 if y < 0
			theta = atanf(y / x); // in [-pi/2, +pi/2]

			if (theta < 0.0f)
				theta += 2.0f * _PI;
		}
	}

	constexpr float ConvertToRadians(float degrees) noexcept
	{
		return degrees * (_PI / 180.0f);
	}

	constexpr float ConvertToDegrees(float radians) noexcept
	{
		return radians * (180.0f / _PI);
	}
}