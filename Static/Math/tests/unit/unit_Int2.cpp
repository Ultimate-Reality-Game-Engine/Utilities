#include <gtest/gtest.h>
#include <Int2.h>

#include <array>

using namespace UltReality::Math;

class Int2Test : public ::testing::Test
{
protected:
	Int2 int2{ 3, 4 };
	AInt2 aint2{ 5, 6 };

	// Helper function to extract vector components
	static std::array<int32_t, 4> ExtractComponents(const A_VECTOR vec)
	{
#if defined(_SSE2_INTRINSICS_)
		alignas(16) int32_t components[4];
		_mm_storeu_ps(reinterpret_cast<float*>(components), vec);
		return { components[0], components[1], components[2], components[3] };
#else
		return { vec.vector4_u32[0], vec.vector4_u32[1], vec.vector4_u32[2], vec.vector4_u32[3] };
#endif
	}
};

TEST_F(Int2Test, LoadInt2)
{
	VECTOR vec = int2.Load();
	auto comp = ExtractComponents(vec);

	EXPECT_EQ(comp[0], 3);
	EXPECT_EQ(comp[1], 4);
	EXPECT_EQ(comp[2], 0);
	EXPECT_EQ(comp[3], 0);
}