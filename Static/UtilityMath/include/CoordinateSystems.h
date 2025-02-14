#ifndef ULTREALITY_MATH_COORDINATE_SYSTEMS_H
#define ULTREALITY_MATH_COORDINATE_SYSTEMS_H

#include <VectorMath.h>

namespace UltReality::Math
{
	/// <summary>
	/// Convert a spherical coordinate to a cartesian coordinate
	/// </summary>
	/// <param name="radius">The radius (distance) of the point in spherical space</param>
	/// <param name="theta">Angle from positive x-axis in XY plane to point</param>
	/// <param name="phi">Angle between +z-axis and the line between the origin and the point</param>
	/// <returns>VECTOR with x, y, and z components set to the cartesian equivalent of the spherical point</returns>
	VECTOR SphericalToCartesian(float radius, float theta, float phi);
}

#endif // !ULTREALITY_MATH_COORDINATE_SYSTEMS_H
