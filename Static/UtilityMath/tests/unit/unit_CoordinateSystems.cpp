#include <gtest/gtest.h>

#include <cmath>

#include <VectorMath.h>
#include <CoordinateSystems.h>

using namespace UltReality::Math;

VECTOR vTrue = Vector::TrueInt();

VECTOR vEPSILON = Vector::Replicate(1e-5f);

static void VEC_CALLCONV ExpectVectorNear(A_VECTOR result, A_VECTOR expected)
{
	VECTOR cmp = Vector::NearEqual(result, expected, vEPSILON);
	uint32_t CR = 0;
	Vector::EqualIntR(&CR, cmp, vTrue);

	ASSERT_TRUE(CompareAllTrue(CR));
}

// Test case: With a radius of zero, all spatial components should be zero and w should be 1.
TEST(SphericalToCartesianTests, ZeroRadius)
{
    VECTOR result = SphericalToCartesian(0.0f, 0.0f, 0.0f);
    VECTOR expected = Vector::Set(0.0f, 0.0f, 0.0f, 1.0f);
    ExpectVectorNear(result, expected);
}

// Test case: For radius=1, theta=0, phi=pi/2.
// Expected: 
//   x = 1 * sin(pi/2) * cos(0) = 1 * 1 * 1 = 1,
//   y = 1 * cos(pi/2) = 1 * 0 = 0,
//   z = 1 * sin(pi/2) * sin(0) = 1 * 1 * 0 = 0.
TEST(SphericalToCartesianTests, UnitSphereThetaZeroPhiHalfPi)
{
    float radius = 1.0f;
    float theta = 0.0f;
    float phi = static_cast<float>(_PI / 2.0);
    VECTOR result = SphericalToCartesian(radius, theta, phi);
    VECTOR expected = Vector::Set(1.0f, 0.0f, 0.0f, 1.0f);
    ExpectVectorNear(result, expected);
}

// Test case: For radius=1, theta=pi/2, phi=pi/2.
// Expected: 
//   x = 1 * sin(pi/2) * cos(pi/2) = 1 * 1 * 0 = 0,
//   y = 1 * cos(pi/2) = 0,
//   z = 1 * sin(pi/2) * sin(pi/2) = 1 * 1 * 1 = 1.
TEST(SphericalToCartesianTests, UnitSphereThetaHalfPiPhiHalfPi)
{
    float radius = 1.0f;
    float theta = static_cast<float>(_PI / 2.0);
    float phi = static_cast<float>(_PI / 2.0);
    VECTOR result = SphericalToCartesian(radius, theta, phi);
    VECTOR expected = Vector::Set(0.0f, 0.0f, 1.0f, 1.0f);
    ExpectVectorNear(result, expected);
}

// Test case: For radius=1, theta=0, phi=0.
// Expected: 
//   x = 1 * sin(0) * cos(0) = 0,
//   y = 1 * cos(0) = 1,
//   z = 1 * sin(0) * sin(0) = 0.
TEST(SphericalToCartesianTests, UnitSphereThetaZeroPhiZero)
{
    float radius = 1.0f;
    float theta = 0.0f;
    float phi = 0.0f;
    VECTOR result = SphericalToCartesian(radius, theta, phi);
    VECTOR expected = Vector::Set(0.0f, 1.0f, 0.0f, 1.0f);
    ExpectVectorNear(result, expected);
}