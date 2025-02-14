#include <gtest/gtest.h>
#include <CoordinateSystems.h>

#include <cmath>

using namespace UltReality::Math;

VECTOR vTrue = VEC::TrueInt();

VECTOR vEPSILON = = VEC::Replicate(1e-5f);

static void VEC_CALLCONV ExpectVectorNear(A_VECTOR result, A_VECTOR expected)
{
	VECTOR cmp = VEC::NearEqual(result, expected, vEPSILON);

	ASSERT_TRUE(VEC::EqualInt(cmp, vTrue));
}

TEST(SphericalToCartesianTests, ZeroRadius)
{
	VECTOR vResult = SphericalToCartesian(0.0f, 0.0f, 0.0f);
	VECTOR vExpected = VEC::Set(0.0f, 0.0f, 0.0f, 1.0f);

	ExpectVectorNear(result, vExpected);
}