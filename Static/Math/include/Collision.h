#ifndef ULTREALITY_MATH_COLLISION_H
#define ULTREALITY_MATH_COLLISION_H

#include <SSE2VectorConfig.h>
#include <Miscellaneous.h>
#include <MATRIX.h>
#include <Float3.h>
#include <Float4.h>

namespace UltReality::Math
{
    namespace Collision
    {
        enum ContainmentType
        {
            DISJOINT = 0, 
            INTERSECTS = 1, 
            CONTAINS = 2
        };

        enum PlaneIntersectionType
        {
            FRONT = 0, 
            INTERSECTING = 1, 
            BACK = 2
        };

        struct BoundingBox;
        struct BoundingOrientedBox;
        struct BoundingFrustum;

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable:4324 4820)
    // C4324: alignment padding warnings
    // C4820: Off by default noise
#endif

        struct BoundingSphere
        {
            Float3 center; // Center of the sphere
            float radius;  // Radius of the sphere

            BoundingSphere() noexcept : center(0, 0, 0), radius(1.0f) {}

            BoundingSphere(const BoundingSphere&) = default;
            BoundingSphere& operator=(const BoundingSphere&) = default;

            BoundingSphere(BoundingSphere&&) = default;
            BoundingSphere& operator=(BoundingSphere&&) = default;

            constexpr BoundingSphere(_In_ const Float3& center, _In_ float radius) noexcept
                : center(center), radius(radius) {}

            // Methods
            // Transform a sphere by an angle preserving transform
            void VEC_CALLCONV Transform(_Out_ BoundingSphere& out, _In_ A_MATRIX m) const noexcept;
            void VEC_CALLCONV Transform(_Out_ BoundingSphere& out, _In_ float scale, _In_ A_VECTOR rotation, _In_ A_VECTOR translation) const noexcept;
            // Transform the sphere

            // Point in sphere test
            ContainmentType VEC_CALLCONV Contains(_In_ A_VECTOR point) const noexcept;
            // Triangle in sphere test
            ContainmentType VEC_CALLCONV Contains(_In_ A_VECTOR V0, _In_ A_VECTOR V1, _In_ A_VECTOR V2) const noexcept;
            // Sphere in sphere test
            ContainmentType Contains(_In_ const BoundingSphere& sh) const noexcept;
            // Axis-aligned box in sphere test
            ContainmentType Contains(_In_ const BoundingBox& box) const noexcept;
            // Oriented box in sphere test
            ContainmentType Contains(_In_ const BoundingOrientedBox& box) const noexcept;
            // Frustum in sphere test
            ContainmentType Contains(_In_ const BoundingFrustum& fr) const noexcept;

            // Sphere vs. sphere test
            bool Intersects(_In_ const BoundingSphere& sh) const noexcept;
            // Box vs. sphere test
            bool Intersects(_In_ const BoundingBox& box) const noexcept;
            bool Intersects(_In_ const BoundingOrientedBox& box) const noexcept;
            // Frustum vs. sphere test
            bool Intersects(_In_ const BoundingFrustum& fr) const noexcept;

            // Triangle vs. sphere test
            bool VEC_CALLCONV Intersects(_In_ A_VECTOR V0, _In_ A_VECTOR V1, _In_ A_VECTOR V2) const noexcept;

            // Sphere-plane intersection
            PlaneIntersectionType VEC_CALLCONV Intersects(_In_ A_VECTOR plane) const noexcept;
            // Plane-sphere test

            // Compute the intersection of a ray (origin, direction) with a sphere
            bool VEC_CALLCONV Intersects(_In_ A_VECTOR origin, _In_ A_VECTOR direction, _Out_ float& dist) const noexcept;

            // Test a sphere vs. 6 planes (typically forming a frustum)
            ContainmentType VEC_CALLCONV ContainedBy(_In_ A_VECTOR plane0, _In_ A_VECTOR plane1, _In_ A_VECTOR plane2,
                _In_ B_VECTOR plane3, _In_ C_VECTOR plane4, _In_ C_VECTOR plane5) const noexcept;

            // Static methods
            // Creates a bounding sphere that contains two other bounding spheres
            static void CreateMerged(_Out_ BoundingSphere& out, _In_ const BoundingSphere& S1, _In_ const BoundingSphere& S2) noexcept;

            // Create sphere enscribing bounding box
            static void CreateFromBoundingBox(_Out_ BoundingSphere& out, _In_ const BoundingBox& box) noexcept;
            static void CreateFromBoundingBox(_Out_ BoundingSphere& out, _In_ const BoundingOrientedBox& box) noexcept;

            //-----------------------------------------------------------------------------
            // Find the approximate smallest enclosing bounding sphere for a set of
            // points. Exact computation of the smallest enclosing bounding sphere is
            // possible but is slower and requires a more complex algorithm.
            // The algorithm is based on  Jack Ritter, "An Efficient Bounding Sphere",
            // Graphics Gems.
            //-----------------------------------------------------------------------------
            static void CreateFromPoints(_Out_ BoundingSphere& out, _In_ size_t count,
                _In_reads_bytes_(sizeof(Float3) + stride * (count - 1)) const Float3* pPoints, _In_ size_t stride) noexcept;

            // Create sphere contianing frustum
            static void CreateFromFrustum(_Out_ BoundingSphere& out, _In_ const BoundingFrustum& fr) noexcept;
        };

        struct BoundingBox
        {
            static constexpr size_t CORNER_COUNT = 8;

            Float3 center; // Center of the box.
            Float3 extents; // Distance from the center to each side.

            // Creators
            BoundingBox() noexcept : center(0, 0, 0), extents(1.f, 1.f, 1.f) {}

            BoundingBox(const BoundingBox&) = default;
            BoundingBox& operator=(const BoundingBox&) = default;

            BoundingBox(BoundingBox&&) = default;
            BoundingBox& operator=(BoundingBox&&) = default;

            constexpr BoundingBox(_In_ const Float3& center, _In_ const Float3& extents) noexcept
                : center(center), extents(extents) {}

            // Methods
            // Transform an axis aligned box by an angle preserving transform
            void VEC_CALLCONV Transform(_Out_ BoundingBox& out, _In_ A_MATRIX m) const noexcept;
            void VEC_CALLCONV Transform(_Out_ BoundingBox& out, _In_ float scale, _In_ A_VECTOR rotation, _In_ A_VECTOR translation) const noexcept;

            // Get the corner points of the box
            void GetCorners(_Out_writes_(8) Float3* corners) const noexcept;
            // Gets the 8 corners of the box

            // Point in axis-aligned box test
            ContainmentType VEC_CALLCONV Contains(_In_ A_VECTOR point) const noexcept;
            // Triangle in axis-aligned box test
            ContainmentType VEC_CALLCONV Contains(_In_ A_VECTOR V0, _In_ A_VECTOR V1, _In_ A_VECTOR V2) const noexcept;
            // Sphere in axis-aligned box test
            ContainmentType Contains(_In_ const BoundingSphere& sh) const noexcept;
            // Axis-aligned box in axis-aligned box test
            ContainmentType Contains(_In_ const BoundingBox& box) const noexcept;
            // Oriented box in axis-aligned box test
            ContainmentType Contains(_In_ const BoundingOrientedBox& box) const noexcept;
            // Frustum in axis-aligned box testS
            ContainmentType Contains(_In_ const BoundingFrustum& fr) const noexcept;

            // Sphere vs axis-aligned box test
            bool Intersects(_In_ const BoundingSphere& sh) const noexcept;
            // Axis-aligned box vs. axis-aligned box test
            bool Intersects(_In_ const BoundingBox& box) const noexcept;
            // Oriented box vs. axis-aligned box test
            bool Intersects(_In_ const BoundingOrientedBox& box) const noexcept;
            // Frustum vs. axis-aligned box test
            bool Intersects(_In_ const BoundingFrustum& fr) const noexcept;
            // Triangle vs. axis aligned box test
            bool VEC_CALLCONV Intersects(_In_ A_VECTOR V0, _In_ A_VECTOR V1, _In_ A_VECTOR V2) const noexcept;
            // Triangle-Box test

            // Plane-box test
            PlaneIntersectionType VEC_CALLCONV Intersects(_In_ A_VECTOR plane) const noexcept;

            // Compute the intersection of a ray (Origin, Direction) with an axis aligned
            // box using the slabs method
            bool VEC_CALLCONV Intersects(_In_ A_VECTOR origin, _In_ A_VECTOR direction, _Out_ float& dist) const noexcept;

            // Test an axis alinged box vs 6 planes (typically forming a frustum).
            ContainmentType VEC_CALLCONV ContainedBy(_In_ A_VECTOR plane0, _In_ A_VECTOR plane1, _In_ A_VECTOR plane2, 
                _In_ B_VECTOR plane3, _In_ C_VECTOR plane4, _In_ C_VECTOR plane5) const noexcept;

            // Static methods
            // Create axis-aligned box that contains two other bounding boxes
            static void CreateMerged(_Out_ BoundingBox& out, _In_ const BoundingBox& b1, _In_ const BoundingBox& b2) noexcept;
            // Create axis-aligned box that contains a bounding sphere
            static void CreateFromSphere(_Out_ BoundingBox& out, _In_ const BoundingSphere& sh) noexcept;
            // Create axis-aligned box from min/max points
            static void VEC_CALLCONV CreateFromPoints(_Out_ BoundingBox& out, _In_ A_VECTOR pt1, _In_ A_VECTOR pt2) noexcept;
            // Find the minimum axis aligned bounding box containing a set of points
            static void CreateFromPoints(_Out_ BoundingBox& out, _In_ size_t count,
                _In_reads_bytes_(sizeof(Float3) + stride * (count - 1)) const Float3* pPoints, _In_ size_t stride) noexcept;
        };

        struct BoundingOrientedBox
        {
            static constexpr size_t CORNER_COUNT = 8;

            Float3 center;            // Center of the box.
            Float3 extents;           // Distance from the center to each side.
            Float4 orientation;       // Unit quaternion representing rotation (box -> world).

            // Creators
            BoundingOrientedBox() noexcept : center(0, 0, 0), extents(1.f, 1.f, 1.f), orientation(0, 0, 0, 1.f) {}

            BoundingOrientedBox(const BoundingOrientedBox&) = default;
            BoundingOrientedBox& operator=(const BoundingOrientedBox&) = default;

            BoundingOrientedBox(BoundingOrientedBox&&) = default;
            BoundingOrientedBox& operator=(BoundingOrientedBox&&) = default;

            constexpr BoundingOrientedBox(_In_ const Float3& center, _In_ const Float3& extents, _In_ const Float4& orientation) noexcept
                : center(center), extents(extents), orientation(orientation) {}

            // Methods
            // Transform an oriented box by an angle preserving transform
            void VEC_CALLCONV Transform(_Out_ BoundingOrientedBox& out, _In_ A_MATRIX M) const noexcept;
            void VEC_CALLCONV Transform(_Out_ BoundingOrientedBox& out, _In_ float scale, _In_ A_VECTOR rotation, _In_ A_VECTOR translation) const noexcept;

            // Get the corner points of the box
            void GetCorners(_Out_writes_(8) Float3* corners) const noexcept;

            // Point in oriented box test
            ContainmentType VEC_CALLCONV Contains(_In_ A_VECTOR point) const noexcept;
            // Triangle in oriented bounding box
            ContainmentType VEC_CALLCONV Contains(_In_ A_VECTOR V0, _In_ A_VECTOR V1, _In_ A_VECTOR V2) const noexcept;
            // Sphere in oriented bounding box
            ContainmentType Contains(_In_ const BoundingSphere& sh) const noexcept;
            // Axis aligned box vs. oriented box. Constructs an oriented box and uses
            // the oriented box vs. oriented box test
            ContainmentType Contains(_In_ const BoundingBox& box) const noexcept;
            // Oriented bounding box in oriented bounding box
            ContainmentType Contains(_In_ const BoundingOrientedBox& box) const noexcept;
            // Frustum in oriented bounding box
            ContainmentType Contains(_In_ const BoundingFrustum& fr) const noexcept;

            // Sphere vs. oriented box test
            bool Intersects(_In_ const BoundingSphere& sh) const noexcept;
            // Axis aligned box vs. oriented box. Constructs an oriented box and uses
            // the oriented box vs. oriented box test
            bool Intersects(_In_ const BoundingBox& box) const noexcept;
            // Fast oriented box / oriented box intersection test using the separating axis
            // theorem
            bool Intersects(_In_ const BoundingOrientedBox& box) const noexcept;
            // Frustum vs. oriented box test
            bool Intersects(_In_ const BoundingFrustum& fr) const noexcept;
            // Triangle vs. oriented box test
            bool VEC_CALLCONV Intersects(_In_ A_VECTOR V0, _In_ A_VECTOR V1, _In_ A_VECTOR V2) const noexcept;

            // Plane-OrientedBox test
            PlaneIntersectionType VEC_CALLCONV Intersects(_In_ A_VECTOR plane) const noexcept;

            // Compute the intersection of a ray (Origin, Direction) with an oriented box
            // using the slabs method
            bool VEC_CALLCONV Intersects(_In_ A_VECTOR origin, _In_ A_VECTOR direction, _Out_ float& dist) const noexcept;

            // Test an oriented box vs 6 planes (typically forming a frustum)
            ContainmentType VEC_CALLCONV ContainedBy(_In_ A_VECTOR plane0, _In_ A_VECTOR plane1, _In_ A_VECTOR plane2,
                _In_ B_VECTOR plane3, _In_ C_VECTOR plane4, _In_ C_VECTOR plane5) const noexcept;

            // Static methods
            // Create oriented bounding box from axis-aligned bounding box
            static void CreateFromBoundingBox(_Out_ BoundingOrientedBox& out, _In_ const BoundingBox& box) noexcept;
            
            // Find the approximate minimum oriented bounding box containing a set of
            // points.  Exact computation of minimum oriented bounding box is possible but
            // is slower and requires a more complex algorithm.
            // The algorithm works by computing the inertia tensor of the points and then
            // using the eigenvectors of the intertia tensor as the axes of the box.
            // Computing the intertia tensor of the convex hull of the points will usually
            // result in better bounding box but the computation is more complex.
            // Exact computation of the minimum oriented bounding box is possible but the
            // best know algorithm is O(N^3) and is significanly more complex to implement
            static void CreateFromPoints(_Out_ BoundingOrientedBox& out, _In_ size_t count,
                _In_reads_bytes_(sizeof(Float3) + stride * (count - 1)) const Float3* pPoints, _In_ size_t stride) noexcept;
        };

        struct BoundingFrustum
        {
            static constexpr size_t CORNER_COUNT = 8;

            Float3 origin;            // Origin of the frustum (and projection).
            Float4 orientation;       // Quaternion representing rotation.

            float rightSlope;           // Positive X (X/Z)
            float leftSlope;            // Negative X
            float topSlope;             // Positive Y (Y/Z)
            float bottomSlope;          // Negative Y
            float near, far;            // Z of the near plane and far plane.

            // Creators
            BoundingFrustum() noexcept :
                origin(0, 0, 0), orientation(0, 0, 0, 1.f), rightSlope(1.f), leftSlope(-1.f),
                topSlope(1.f), bottomSlope(-1.f), near(0), far(1.f) {}

            BoundingFrustum(const BoundingFrustum&) = default;
            BoundingFrustum& operator=(const BoundingFrustum&) = default;

            BoundingFrustum(BoundingFrustum&&) = default;
            BoundingFrustum& operator=(BoundingFrustum&&) = default;

            constexpr BoundingFrustum(_In_ const Float3& origin, _In_ const Float4& orientation,
                _In_ float rightSlope, _In_ float leftSlope, _In_ float topSlope, _In_ float bottomSlope,
                _In_ float nearPlane, _In_ float farPlane) noexcept
                : origin(origin), orientation(orientation),
                rightSlope(rightSlope), leftSlope(leftSlope), topSlope(topSlope), bottomSlope(bottomSlope),
                near(nearPlane), far(farPlane) {}
            BoundingFrustum(_In_ B_MATRIX projection, bool rhcoords = false) noexcept;

            // Methods
            // Transform a frustum by an angle preserving transform
            void VEC_CALLCONV Transform(_Out_ BoundingFrustum& out, _In_ A_MATRIX M) const noexcept;
            void VEC_CALLCONV Transform(_Out_ BoundingFrustum& out, _In_ float scale, _In_ A_VECTOR rotation, _In_ A_VECTOR translation) const noexcept;

            // Get the corner points of the frustum
            void GetCorners(_Out_writes_(8) Float3* corners) const noexcept;

            // Point in frustum test
            ContainmentType VEC_CALLCONV Contains(_In_ A_VECTOR point) const noexcept;
            // Triangle vs frustum test
            ContainmentType VEC_CALLCONV Contains(_In_ A_VECTOR V0, _In_ A_VECTOR V1, _In_ A_VECTOR V2) const noexcept;
            ContainmentType Contains(_In_ const BoundingSphere& sp) const noexcept;
            ContainmentType Contains(_In_ const BoundingBox& box) const noexcept;
            ContainmentType Contains(_In_ const BoundingOrientedBox& box) const noexcept;
            ContainmentType Contains(_In_ const BoundingFrustum& fr) const noexcept;

            // Exact sphere vs frustum test.  The algorithm first checks the sphere against
            // the planes of the frustum, then if the plane checks were indeterminate finds
            // the nearest feature (plane, line, point) on the frustum to the center of the
            // sphere and compares the distance to the nearest feature to the radius of the
            // sphere
            bool Intersects(_In_ const BoundingSphere& sh) const noexcept;
            // Exact axis aligned box vs frustum test.  Constructs an oriented box and uses
            // the oriented box vs frustum test
            bool Intersects(_In_ const BoundingBox& box) const noexcept;
            // Exact oriented box vs frustum test
            bool Intersects(_In_ const BoundingOrientedBox& box) const noexcept;
            // Exact frustum vs frustum test
            bool Intersects(_In_ const BoundingFrustum& fr) const noexcept;
            // Triangle vs frustum test
            bool VEC_CALLCONV Intersects(_In_ A_VECTOR V0, _In_ A_VECTOR V1, _In_ A_VECTOR V2) const noexcept;

            // Plane-Frustum test
            PlaneIntersectionType VEC_CALLCONV Intersects(_In_ A_VECTOR plane) const noexcept;

            // Ray vs. frustum test
            bool VEC_CALLCONV Intersects(_In_ A_VECTOR rayOrigin, _In_ A_VECTOR direction, _Out_ float& dist) const noexcept;

            // Test a frustum vs 6 planes (typically forming another frustum)
            ContainmentType VEC_CALLCONV ContainedBy(_In_ A_VECTOR plane0, _In_ A_VECTOR plane1, _In_ A_VECTOR plane2,
                _In_ B_VECTOR plane3, _In_ C_VECTOR Plane4, _In_ C_VECTOR Plane5) const noexcept;

            // Build the 6 frustum planes from a frustum.
            //
            // The intended use for these routines is for fast culling to a view frustum.
            // When the volume being tested against a view frustum is small relative to the
            // view frustum it is usually either inside all six planes of the frustum
            // (CONTAINS) or outside one of the planes of the frustum (DISJOINT). If neither
            // of these cases is true then it may or may not be intersecting the frustum
            // (INTERSECTS)
            void GetPlanes(_Out_opt_ VECTOR* nearPlane, _Out_opt_ VECTOR* farPlane, _Out_opt_ VECTOR* rightPlane,
                _Out_opt_ VECTOR* leftPlane, _Out_opt_ VECTOR* topPlane, _Out_opt_ VECTOR* bottomPlane) const noexcept;

            // Static methods
            
            // Build a frustum from a persepective projection matrix.  The matrix may only
            // contain a projection; any rotation, translation or scale will cause the
            // constructed frustum to be incorrect
            static void VEC_CALLCONV CreateFromMatrix(_Out_ BoundingFrustum& out, _In_ A_MATRIX projection, bool rhcoords = false) noexcept;
        };

        namespace TriangleTests
        {
            bool VEC_CALLCONV Intersects(_In_ A_VECTOR origin, _In_ A_VECTOR direction, _In_ A_VECTOR V0, _In_ B_VECTOR V1, _In_ C_VECTOR V2, _Out_ float& dist) noexcept;
            // Ray-Triangle

            bool VEC_CALLCONV Intersects(_In_ A_VECTOR A0, _In_ A_VECTOR A1, _In_ A_VECTOR A2, _In_ B_VECTOR B0, _In_ C_VECTOR B1, _In_ C_VECTOR B2) noexcept;
            // Triangle-Triangle

            PlaneIntersectionType VEC_CALLCONV Intersects(_In_ A_VECTOR V0, _In_ A_VECTOR V1, _In_ A_VECTOR V2, _In_ B_VECTOR plane) noexcept;
            // Plane-Triangle

            ContainmentType VEC_CALLCONV ContainedBy(_In_ A_VECTOR V0, _In_ A_VECTOR V1, _In_ A_VECTOR V2,
                _In_ B_VECTOR plane0, _In_ C_VECTOR plane1, _In_ C_VECTOR plane2,
                _In_ D_VECTOR plane3, _In_ D_VECTOR plane4, _In_ D_VECTOR plane5) noexcept;
            // Test a triangle against six planes at once (see BoundingFrustum::GetPlanes)
        } // namespace TriangleTests
        
#ifdef _MSC_VER
#pragma warning(pop)
#endif

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4068 4365 4616 6001)
     // C4068/4616: ignore unknown pragmas
     // C4365: Off by default noise
     // C6001: False positives
#endif

#ifdef _PREFAST_
#pragma prefast(push)
#pragma prefast(disable : 25000, "A_VECTOR is 16 bytes")
#pragma prefast(disable : 26495, "Union initialization confuses /analyze")
#endif

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wfloat-equal"
#pragma clang diagnostic ignored "-Wunknown-warning-option"
#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"
#endif

#include "Collision.inl"

#ifdef __clang__
#pragma clang diagnostic pop
#endif
#ifdef _PREFAST_
#pragma prefast(pop)
#endif
#ifdef _MSC_VER
#pragma warning(pop)
#endif
    } // namespace Collision
} // namespace UltReality::Math
#endif // !ULTREALITY_MATH_COLLISION_H
