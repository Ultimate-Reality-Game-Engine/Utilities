#ifndef ULTREALITY_MATH_COLLISION_INL
#define ULTREALITY_MATH_COLLISION_INL

#if defined(__GNUC__) or defined(__clang__)
#define FORCE_INLINE inline __attribute__((always_inline))
#elif defined(_MSC_VER)
#define FORCE_INLINE __forceinline
#else
#define FORCE_INLINE inline
#endif

namespace UltReality::Math
{
    namespace Collision
    {
        VEC_GLOBCONST VECTOR_F32 g_BoxOffset[8] =
        {
            { { { -1.0f, -1.0f,  1.0f, 0.0f } } },
            { { {  1.0f, -1.0f,  1.0f, 0.0f } } },
            { { {  1.0f,  1.0f,  1.0f, 0.0f } } },
            { { { -1.0f,  1.0f,  1.0f, 0.0f } } },
            { { { -1.0f, -1.0f, -1.0f, 0.0f } } },
            { { {  1.0f, -1.0f, -1.0f, 0.0f } } },
            { { {  1.0f,  1.0f, -1.0f, 0.0f } } },
            { { { -1.0f,  1.0f, -1.0f, 0.0f } } },
        };

        VEC_GLOBCONST VECTOR_F32 g_RayEpsilon = { { { 1e-20f, 1e-20f, 1e-20f, 1e-20f } } };
        VEC_GLOBCONST VECTOR_F32 g_RayNegEpsilon = { { { -1e-20f, -1e-20f, -1e-20f, -1e-20f } } };
        VEC_GLOBCONST VECTOR_F32 g_FltMin = { { { -FLT_MAX, -FLT_MAX, -FLT_MAX, -FLT_MAX } } };
        VEC_GLOBCONST VECTOR_F32 g_FltMax = { { { FLT_MAX, FLT_MAX, FLT_MAX, FLT_MAX } } };

        namespace
        {
            //-----------------------------------------------------------------------------
            // Return true if any of the elements of a 3 vector are equal to 0xffffffff.
            // Slightly more efficient than using XMVector3EqualInt.
            //-----------------------------------------------------------------------------
            FORCE_INLINE bool Vector3AnyTrue(_In_ A_VECTOR v) noexcept
            {
                // Duplicate the fourth element from the first element.
                VECTOR C = Vector::Swizzle<SWIZZLE_X, SWIZZLE_Y, SWIZZLE_Z, SWIZZLE_X>(v);

                return CompareAnyTrue(Vector4::EqualIntR(C, Vector::TrueInt()));
            }

            //-----------------------------------------------------------------------------
            // Return true if all of the elements of a 3 vector are equal to 0xffffffff.
            // Slightly more efficient than using XMVector3EqualInt.
            //-----------------------------------------------------------------------------
            FORCE_INLINE bool Vector3AllTrue(_In_ A_VECTOR v) noexcept
            {
                // Duplicate the fourth element from the first element.
                VECTOR C = Vector::Swizzle<SWIZZLE_X, SWIZZLE_Y, SWIZZLE_Z, SWIZZLE_X>(v);

                return CompareAllTrue(Vector4::EqualIntR(C, Vector::TrueInt()));
            }

#if defined(_PREFAST_) || !defined(NDEBUG)

            VEC_GLOBCONST VECTOR_F32 g_UnitVectorEpsilon = { { { 1.0e-4f, 1.0e-4f, 1.0e-4f, 1.0e-4f } } };
            VEC_GLOBCONST VECTOR_F32 g_UnitQuaternionEpsilon = { { { 1.0e-4f, 1.0e-4f, 1.0e-4f, 1.0e-4f } } };
            VEC_GLOBCONST VECTOR_F32 g_UnitPlaneEpsilon = { { { 1.0e-4f, 1.0e-4f, 1.0e-4f, 1.0e-4f } } };

            //-----------------------------------------------------------------------------
            // Return true if the vector is a unit vector (length == 1).
            //-----------------------------------------------------------------------------
            FORCE_INLINE bool Vector3IsUnit(_In_ A_VECTOR v) noexcept
            {
                VECTOR Difference = Vector::Subtract(Vector3::Length(v), Vector::SplatOne());
                return Vector4::Less(Vector::Abs(Difference), g_UnitVectorEpsilon);
            }

            //-----------------------------------------------------------------------------
            // Return true if the quaterion is a unit quaternion.
            //-----------------------------------------------------------------------------
            FORCE_INLINE bool QuaternionIsUnit(_In_ A_VECTOR q) noexcept
            {
                VECTOR Difference = Vector::Subtract(Vector4::Length(q), Vector::SplatOne());
                return Vector4::Less(Vector::Abs(Difference), g_UnitQuaternionEpsilon);
            }

            //-----------------------------------------------------------------------------
            // Return true if the plane is a unit plane.
            //-----------------------------------------------------------------------------
            FORCE_INLINE bool PlaneIsUnit(_In_ A_VECTOR plane) noexcept
            {
                VECTOR Difference = Vector::Subtract(Vector3::Length(plane), Vector::SplatOne());
                return Vector4::Less(Vector::Abs(Difference), g_UnitPlaneEpsilon);
            }

#endif // _PREFAST_ || !NDEBUG

            //-----------------------------------------------------------------------------
            FORCE_INLINE VECTOR PlaneTransform(_In_ A_VECTOR plane, _In_ A_VECTOR rotation, _In_ A_VECTOR translation) noexcept
            {
                VECTOR vNormal = Vector3::Rotate(plane, rotation);
                VECTOR vD = Vector::Subtract(Vector::SplatW(plane), Vector3::Dot(vNormal, translation));

                return Vector::Insert<0, 0, 0, 0, 1>(vNormal, vD);
            }

            //-----------------------------------------------------------------------------
            // Return the point on the line segement (S1, S2) nearest the point P.
            //-----------------------------------------------------------------------------
            FORCE_INLINE VECTOR PointOnLineSegmentNearestPoint(_In_ A_VECTOR S1, _In_ A_VECTOR S2, _In_ A_VECTOR p) noexcept
            {
                VECTOR Dir = Vector::Subtract(S2, S1);
                VECTOR Projection = Vector::Subtract(Vector3::Dot(p, Dir), Vector3::Dot(S1, Dir));
                VECTOR LengthSq = Vector3::Dot(Dir, Dir);

                VECTOR t = Vector::Multiply(Projection, Vector::Reciprocal(LengthSq));
                VECTOR Point = Vector::MultiplyAdd(t, Dir, S1);

                // t < 0
                VECTOR SelectS1 = Vector::Less(Projection, Vector::Zero());
                Point = Vector::Select(Point, S1, SelectS1);

                // t > 1
                VECTOR SelectS2 = Vector::Greater(Projection, LengthSq);
                Point = Vector::Select(Point, S2, SelectS2);

                return Point;
            }

            //-----------------------------------------------------------------------------
            // Test if the point (P) on the plane of the triangle is inside the triangle
            // (V0, V1, V2).
            //-----------------------------------------------------------------------------
            FORCE_INLINE VECTOR VEC_CALLCONV PointOnPlaneInsideTriangle(_In_ A_VECTOR p, _In_ A_VECTOR V0, _In_ A_VECTOR V1, _In_ B_VECTOR V2) noexcept
            {
                // Compute the triangle normal.
                VECTOR N = Vector3::Cross(Vector::Subtract(V2, V0), Vector::Subtract(V1, V0));

                // Compute the cross products of the vector from the base of each edge to
                // the point with each edge vector.
                VECTOR C0 = Vector3::Cross(Vector::Subtract(p, V0), Vector::Subtract(V1, V0));
                VECTOR C1 = Vector3::Cross(Vector::Subtract(p, V1), Vector::Subtract(V2, V1));
                VECTOR C2 = Vector3::Cross(Vector::Subtract(p, V2), Vector::Subtract(V0, V2));

                // If the cross product points in the same direction as the normal the the
                // point is inside the edge (it is zero if is on the edge).
                VECTOR zero = Vector::Zero();
                VECTOR Inside0 = Vector::GreaterOrEqual(Vector3::Dot(C0, N), zero);
                VECTOR Inside1 = Vector::GreaterOrEqual(Vector3::Dot(C1, N), zero);
                VECTOR Inside2 = Vector::GreaterOrEqual(Vector3::Dot(C2, N), zero);

                // If the point inside all of the edges it is inside.
                return Vector::AndInt(Vector::AndInt(Inside0, Inside1), Inside2);
            }

            //-----------------------------------------------------------------------------
            FORCE_INLINE bool SolveCubic(_In_ float e, _In_ float f, _In_ float g, _Out_ float* t, _Out_ float* u, _Out_ float* v) noexcept
            {
                float p, q, h, rc, d, theta, costh3, sinth3;

                p = f - e * e / 3.0f;
                q = g - e * f / 3.0f + e * e * e * 2.0f / 27.0f;
                h = q * q / 4.0f + p * p * p / 27.0f;

                if (h > 0)
                {
                    *t = *u = *v = 0.f;
                    return false; // only one real root
                }

                if ((h == 0) && (q == 0)) // all the same root
                {
                    *t = -e / 3;
                    *u = -e / 3;
                    *v = -e / 3;

                    return true;
                }

                d = sqrtf(q * q / 4.0f - h);
                if (d < 0)
                    rc = -powf(-d, 1.0f / 3.0f);
                else
                    rc = powf(d, 1.0f / 3.0f);

                theta = ScalarACos(-q / (2.0f * d));
                costh3 = ScalarCos(theta / 3.0f);
                sinth3 = sqrtf(3.0f) * ScalarSine(theta / 3.0f);
                *t = 2.0f * rc * costh3 - e / 3.0f;
                *u = -rc * (costh3 + sinth3) - e / 3.0f;
                *v = -rc * (costh3 - sinth3) - e / 3.0f;

                return true;
            }

            //-----------------------------------------------------------------------------
            FORCE_INLINE VECTOR CalculateEigenVector(_In_ float m11, _In_ float m12, _In_ float m13,
                _In_ float m22, _In_ float m23, _In_ float m33, _In_ float e) noexcept
            {
                float fTmp[3];
                fTmp[0] = m12 * m23 - m13 * (m22 - e);
                fTmp[1] = m13 * m12 - m23 * (m11 - e);
                fTmp[2] = (m11 - e) * (m22 - e) - m12 * m12;

                VECTOR vTmp = Vector::LoadFloat3(reinterpret_cast<const Float3*>(fTmp));

                if (Vector3::Equal(vTmp, Vector::Zero())) // planar or linear
                {
                    float f1, f2, f3;

                    // we only have one equation - find a valid one
                    if ((m11 - e != 0) || (m12 != 0) || (m13 != 0))
                    {
                        f1 = m11 - e; f2 = m12; f3 = m13;
                    }
                    else if ((m12 != 0) || (m22 - e != 0) || (m23 != 0))
                    {
                        f1 = m12; f2 = m22 - e; f3 = m23;
                    }
                    else if ((m13 != 0) || (m23 != 0) || (m33 - e != 0))
                    {
                        f1 = m13; f2 = m23; f3 = m33 - e;
                    }
                    else
                    {
                        // error, we'll just make something up - we have NO context
                        f1 = 1.0f; f2 = 0.0f; f3 = 0.0f;
                    }

                    if (f1 == 0)
                        vTmp = Vector::SetX(vTmp, 0.0f);
                    else
                        vTmp = Vector::SetX(vTmp, 1.0f);

                    if (f2 == 0)
                        vTmp = Vector::SetY(vTmp, 0.0f);
                    else
                        vTmp = Vector::SetY(vTmp, 1.0f);

                    if (f3 == 0)
                    {
                        vTmp = Vector::SetZ(vTmp, 0.0f);
                        // recalculate y to make equation work
                        if (m12 != 0)
                            vTmp = Vector::SetY(vTmp, -f1 / f2);
                    }
                    else
                    {
                        vTmp = Vector::SetZ(vTmp, (f2 - f1) / f3);
                    }
                }

                if (Vector::GetX(Vector3::LengthSq(vTmp)) > 1e-5f)
                {
                    return Vector3::Normalize(vTmp);
                }
                else
                {
                    // Multiply by a value large enough to make the vector non-zero.
                    vTmp = Vector::Scale(vTmp, 1e5f);
                    return Vector3::Normalize(vTmp);
                }
            }

            //-----------------------------------------------------------------------------
            FORCE_INLINE bool CalculateEigenVectors(_In_ float m11, _In_ float m12, _In_ float m13,
                _In_ float m22, _In_ float m23, _In_ float m33,
                _In_ float e1, _In_ float e2, _In_ float e3,
                _Out_ VECTOR* pV1, _Out_ VECTOR* pV2, _Out_ VECTOR* pV3) noexcept
            {
                *pV1 = CalculateEigenVector(m11, m12, m13, m22, m23, m33, e1);
                *pV2 = CalculateEigenVector(m11, m12, m13, m22, m23, m33, e2);
                *pV3 = CalculateEigenVector(m11, m12, m13, m22, m23, m33, e3);

                bool v1z = false;
                bool v2z = false;
                bool v3z = false;

                VECTOR zero = Vector::Zero();

                if (Vector3::Equal(*pV1, zero))
                    v1z = true;

                if (Vector3::Equal(*pV2, zero))
                    v2z = true;

                if (Vector3::Equal(*pV3, zero))
                    v3z = true;

                bool e12 = (fabsf(Vector::GetX(Vector3::Dot(*pV1, *pV2))) > 0.1f); // check for non-orthogonal vectors
                bool e13 = (fabsf(Vector::GetX(Vector3::Dot(*pV1, *pV3))) > 0.1f);
                bool e23 = (fabsf(Vector::GetX(Vector3::Dot(*pV2, *pV3))) > 0.1f);

                if ((v1z && v2z && v3z) || (e12 && e13 && e23) ||
                    (e12 && v3z) || (e13 && v2z) || (e23 && v1z)) // all eigenvectors are 0- any basis set
                {
                    *pV1 = g_IdentityR0.v;
                    *pV2 = g_IdentityR1.v;
                    *pV3 = g_IdentityR2.v;
                    return true;
                }

                if (v1z && v2z)
                {
                    VECTOR vTmp = Vector3::Cross(g_IdentityR1, *pV3);
                    if (Vector::GetX(Vector3::LengthSq(vTmp)) < 1e-5f)
                    {
                        vTmp = Vector3::Cross(g_IdentityR0, *pV3);
                    }
                    *pV1 = Vector3::Normalize(vTmp);
                    *pV2 = Vector3::Cross(*pV3, *pV1);
                    return true;
                }

                if (v3z && v1z)
                {
                    VECTOR vTmp = Vector3::Cross(g_IdentityR1, *pV2);
                    if (Vector::GetX(Vector3::LengthSq(vTmp)) < 1e-5f)
                    {
                        vTmp = Vector3::Cross(g_IdentityR0, *pV2);
                    }
                    *pV3 = Vector3::Normalize(vTmp);
                    *pV1 = Vector3::Cross(*pV2, *pV3);
                    return true;
                }

                if (v2z && v3z)
                {
                    VECTOR vTmp = Vector3::Cross(g_IdentityR1, *pV1);
                    if (Vector::GetX(Vector3::LengthSq(vTmp)) < 1e-5f)
                    {
                        vTmp = Vector3::Cross(g_IdentityR0, *pV1);
                    }
                    *pV2 = Vector3::Normalize(vTmp);
                    *pV3 = Vector3::Cross(*pV1, *pV2);
                    return true;
                }

                if ((v1z) || e12)
                {
                    *pV1 = Vector3::Cross(*pV2, *pV3);
                    return true;
                }

                if ((v2z) || e23)
                {
                    *pV2 = Vector3::Cross(*pV3, *pV1);
                    return true;
                }

                if ((v3z) || e13)
                {
                    *pV3 = Vector3::Cross(*pV1, *pV2);
                    return true;
                }

                return true;
            }

            //-----------------------------------------------------------------------------
            FORCE_INLINE bool CalculateEigenVectorsFromCovarianceMatrix(_In_ float Cxx, _In_ float Cyy, _In_ float Czz,
                _In_ float Cxy, _In_ float Cxz, _In_ float Cyz,
                _Out_ VECTOR* pV1, _Out_ VECTOR* pV2, _Out_ VECTOR* pV3) noexcept
            {
                // Calculate the eigenvalues by solving a cubic equation.
                float e = -(Cxx + Cyy + Czz);
                float f = Cxx * Cyy + Cyy * Czz + Czz * Cxx - Cxy * Cxy - Cxz * Cxz - Cyz * Cyz;
                float g = Cxy * Cxy * Czz + Cxz * Cxz * Cyy + Cyz * Cyz * Cxx - Cxy * Cyz * Cxz * 2.0f - Cxx * Cyy * Czz;

                float ev1, ev2, ev3;
                if (!SolveCubic(e, f, g, &ev1, &ev2, &ev3))
                {
                    // set them to arbitrary orthonormal basis set
                    *pV1 = g_IdentityR0.v;
                    *pV2 = g_IdentityR1.v;
                    *pV3 = g_IdentityR2.v;
                    return false;
                }

                return CalculateEigenVectors(Cxx, Cxy, Cxz, Cyy, Cyz, Czz, ev1, ev2, ev3, pV1, pV2, pV3);
            }

            //-----------------------------------------------------------------------------
            FORCE_INLINE void VEC_CALLCONV FastIntersectTrianglePlane(
                A_VECTOR V0, A_VECTOR V1, A_VECTOR V2, 
                B_VECTOR plane, 
                VECTOR& outside, VECTOR& inside) noexcept
            {
                // Plane0
                VECTOR Dist0 = Vector4::Dot(V0, plane);
                VECTOR Dist1 = Vector4::Dot(V1, plane);
                VECTOR Dist2 = Vector4::Dot(V2, plane);

                VECTOR MinDist = Vector::Min(Dist0, Dist1);
                MinDist = Vector::Min(MinDist, Dist2);

                VECTOR MaxDist = Vector::Max(Dist0, Dist1);
                MaxDist = Vector::Max(MaxDist, Dist2);

                VECTOR zero = Vector::Zero();

                // Outside the plane?
                outside = Vector::Greater(MinDist, zero);

                // Fully inside the plane?
                inside = Vector::Less(MaxDist, zero);
            }

            //-----------------------------------------------------------------------------
            FORCE_INLINE void FastIntersectSpherePlane(
                _In_ A_VECTOR center, _In_ A_VECTOR radius, _In_ A_VECTOR plane,
                _Out_ VECTOR& outside, _Out_ VECTOR& inside) noexcept
            {
                VECTOR Dist = Vector4::Dot(center, plane);

                // Outside the plane?
                outside = Vector::Greater(Dist, radius);

                // Fully inside the plane?
                inside = Vector::Less(Dist, Vector::Negate(radius));
            }

            //-----------------------------------------------------------------------------
            FORCE_INLINE void FastIntersectAxisAlignedBoxPlane(
                _In_ A_VECTOR center, _In_ A_VECTOR extents, _In_ A_VECTOR plane,
                _Out_ VECTOR& outside, _Out_ VECTOR& inside) noexcept
            {
                // Compute the distance to the center of the box.
                VECTOR Dist = Vector4::Dot(center, plane);

                // Project the axes of the box onto the normal of the plane.  Half the
                // length of the projection (sometime called the "radius") is equal to
                // h(u) * abs(n dot b(u))) + h(v) * abs(n dot b(v)) + h(w) * abs(n dot b(w))
                // where h(i) are extents of the box, n is the plane normal, and b(i) are the
                // axes of the box. In this case b(i) = [(1,0,0), (0,1,0), (0,0,1)].
                VECTOR Radius = Vector3::Dot(extents, Vector::Abs(plane));

                // Outside the plane?
                outside = Vector::Greater(Dist, Radius);

                // Fully inside the plane?
                inside = Vector::Less(Dist, Vector::Negate(Radius));
            }

            //-----------------------------------------------------------------------------
            FORCE_INLINE void VEC_CALLCONV FastIntersectOrientedBoxPlane(
                _In_ A_VECTOR center, _In_ A_VECTOR extents, _In_ A_VECTOR axis0, 
                _In_ B_VECTOR axis1, 
                _In_ C_VECTOR axis2, _In_ C_VECTOR plane,
                _Out_ VECTOR& outside, _Out_ VECTOR& inside) noexcept
            {
                // Compute the distance to the center of the box.
                VECTOR Dist = Vector4::Dot(center, plane);

                // Project the axes of the box onto the normal of the plane.  Half the
                // length of the projection (sometime called the "radius") is equal to
                // h(u) * abs(n dot b(u))) + h(v) * abs(n dot b(v)) + h(w) * abs(n dot b(w))
                // where h(i) are extents of the box, n is the plane normal, and b(i) are the
                // axes of the box.
                VECTOR Radius = Vector3::Dot(plane, axis0);
                Radius = Vector::Insert<0, 0, 1, 0, 0>(Radius, Vector3::Dot(plane, axis1));
                Radius = Vector::Insert<0, 0, 0, 1, 0>(Radius, Vector3::Dot(plane, axis2));
                Radius = Vector3::Dot(extents, Vector::Abs(Radius));

                // Outside the plane?
                outside = Vector::Greater(Dist, Radius);

                // Fully inside the plane?
                inside = Vector::Less(Dist, Vector::Negate(Radius));
            }

            //-----------------------------------------------------------------------------
            FORCE_INLINE void VEC_CALLCONV FastIntersectFrustumPlane(
                _In_ A_VECTOR point0, _In_ A_VECTOR point1, _In_ A_VECTOR point2, 
                _In_ B_VECTOR point3, 
                _In_ C_VECTOR point4, _In_ C_VECTOR point5, 
                _In_ D_VECTOR point6, _In_ D_VECTOR point7, _In_ D_VECTOR plane, 
                _Out_ VECTOR& outside, _Out_ VECTOR& inside) noexcept
            {
                // Find the min/max projection of the frustum onto the plane normal.
                VECTOR Min, Max, Dist;

                Min = Max = Vector3::Dot(plane, point0);

                Dist = Vector3::Dot(plane, point1);
                Min = Vector::Min(Min, Dist);
                Max = Vector::Max(Max, Dist);

                Dist = Vector3::Dot(plane, point2);
                Min = Vector::Min(Min, Dist);
                Max = Vector::Max(Max, Dist);

                Dist = Vector3::Dot(plane, point3);
                Min = Vector::Min(Min, Dist);
                Max = Vector::Max(Max, Dist);

                Dist = Vector3::Dot(plane, point4);
                Min = Vector::Min(Min, Dist);
                Max = Vector::Max(Max, Dist);

                Dist = Vector3::Dot(plane, point5);
                Min = Vector::Min(Min, Dist);
                Max = Vector::Max(Max, Dist);

                Dist = Vector3::Dot(plane, point6);
                Min = Vector::Min(Min, Dist);
                Max = Vector::Max(Max, Dist);

                Dist = Vector3::Dot(plane, point7);
                Min = Vector::Min(Min, Dist);
                Max = Vector::Max(Max, Dist);

                VECTOR PlaneDist = Vector::Negate(Vector::SplatW(plane));

                // Outside the plane?
                outside = Vector::Greater(Min, PlaneDist);

                // Fully inside the plane?
                inside = Vector::Less(Max, PlaneDist);
            }
        }
        
        _Use_decl_annotations_
        FORCE_INLINE void VEC_CALLCONV BoundingSphere::Transform(BoundingSphere& out, A_MATRIX m) const noexcept
        {
            // Load the center of the sphere.
            VECTOR vCenter = Vector::LoadFloat3(&center);

            // Transform the center of the sphere.
            VECTOR C = Vector3::Transform(vCenter, m);

            VECTOR dX = Vector3::Dot(m.r[0], m.r[0]);
            VECTOR dY = Vector3::Dot(m.r[1], m.r[1]);
            VECTOR dZ = Vector3::Dot(m.r[2], m.r[2]);

            VECTOR d = Vector::Max(dX, Vector::Max(dY, dZ));

            // Store the center sphere.
            Vector::StoreFloat3(&out.center, C);

            // Scale the radius of the pshere.
            float scale = sqrtf(Vector::GetX(d));
            out.radius = radius * scale;
        }

        _Use_decl_annotations_
        FORCE_INLINE void VEC_CALLCONV BoundingSphere::Transform(BoundingSphere& out, float scale, A_VECTOR rotation, A_VECTOR translation) const noexcept
        {
            // Load the center of the sphere.
            VECTOR vCenter = Vector::LoadFloat3(&Ccnter);

            // Transform the center of the sphere.
            vCenter = Vector::Add(Vector3::Rotate(Vector::Scale(vCenter, scale), rotation), translation);

            // Store the center sphere.
            Vector::StoreFloat3(&out.center, vCenter);

            // Scale the radius of the pshere.
            out.radius = radius * scale;
        }

        _Use_decl_annotations_
        FORCE_INLINE ContainmentType VEC_CALLCONV BoundingSphere::Contains(A_VECTOR point) const noexcept
        {
            VECTOR vCenter = Vector::LoadFloat3(&center);
            VECTOR vRadius = Vector::ReplicatePtr(&radius);

            VECTOR DistanceSquared = Vector3::LengthSq(Vector::Subtract(point, vCenter));
            VECTOR RadiusSquared = Vector::Multiply(vRadius, vRadius);

            return Vector3::LessOrEqual(DistanceSquared, RadiusSquared) ? CONTAINS : DISJOINT;
        }

        _Use_decl_annotations_
        FORCE_INLINE ContainmentType VEC_CALLCONV BoundingSphere::Contains(A_VECTOR V0, A_VECTOR V1, A_VECTOR V2) const noexcept
        {
            if (!Intersects(V0, V1, V2))
                return DISJOINT;

            VECTOR vCenter = Vector::LoadFloat3(&center);
            VECTOR vRadius = Vector::ReplicatePtr(&radius);
            VECTOR RadiusSquared = Vector::Multiply(vRadius, vRadius);

            VECTOR DistanceSquared = Vector3::LengthSq(Vector::Subtract(V0, vCenter));
            VECTOR Inside = Vector::LessOrEqual(DistanceSquared, RadiusSquared);

            DistanceSquared = Vector3::LengthSq(Vector::Subtract(V1, vCenter));
            Inside = Vector::AndInt(Inside, Vector::LessOrEqual(DistanceSquared, RadiusSquared));

            DistanceSquared = Vector3::LengthSq(Vector::Subtract(V2, vCenter));
            Inside = Vector::AndInt(Inside, Vector::LessOrEqual(DistanceSquared, RadiusSquared));

            return (Vector3::EqualInt(Inside, Vector::TrueInt())) ? CONTAINS : INTERSECTS;
        }

        _Use_decl_annotations_
        FORCE_INLINE ContainmentType BoundingSphere::Contains(const BoundingSphere& sh) const noexcept
        {
            VECTOR Center1 = Vector::LoadFloat3(&center);
            float r1 = radius;

            VECTOR Center2 = Vector::LoadFloat3(&sh.center);
            float r2 = sh.radius;

            VECTOR V = Vector::Subtract(Center2, Center1);

            VECTOR Dist = Vector3::Length(V);

            float d = Vector::GetX(Dist);

            return (r1 + r2 >= d) ? ((r1 - r2 >= d) ? CONTAINS : INTERSECTS) : DISJOINT;
        }

        _Use_decl_annotations_
        FORCE_INLINE ContainmentType BoundingSphere::Contains(const BoundingBox& box) const noexcept
        {
            if (!box.Intersects(*this))
                return DISJOINT;

            VECTOR vCenter = Vector::LoadFloat3(&center);
            VECTOR vRadius = Vector::ReplicatePtr(&radius);
            VECTOR RadiusSq = Vector::Multiply(vRadius, vRadius);

            VECTOR boxCenter = Vector::LoadFloat3(&box.center);
            VECTOR boxExtents = Vector::LoadFloat3(&box.extents);

            VECTOR InsideAll = Vector::TrueInt();

            VECTOR offset = Vector::Subtract(boxCenter, vCenter);

            for (size_t i = 0; i < BoundingBox::CORNER_COUNT; ++i)
            {
                VECTOR C = Vector::MultiplyAdd(boxExtents, g_BoxOffset[i], offset);
                VECTOR d = Vector3::LengthSq(C);
                InsideAll = Vector::AndInt(InsideAll, Vector::LessOrEqual(d, RadiusSq));
            }

            return (Vector3::EqualInt(InsideAll, Vector::TrueInt())) ? CONTAINS : INTERSECTS;
        }

        _Use_decl_annotations_
        FORCE_INLINE ContainmentType BoundingSphere::Contains(const BoundingOrientedBox& box) const noexcept
        {
            if (!box.Intersects(*this))
                return DISJOINT;

            VECTOR vCenter = Vector::LoadFloat3(&center);
            VECTOR vRadius = Vector::ReplicatePtr(&radius);
            VECTOR RadiusSq = Vector::Multiply(vRadius, vRadius);

            VECTOR boxCenter = Vector::LoadFloat3(&box.center);
            VECTOR boxExtents = Vector::LoadFloat3(&box.extents);
            VECTOR boxOrientation = Vector::LoadFloat4(&box.orientation);

#if defined(DEBUG) || defined(_DEBUG)
            assert(QuaternionIsUnit(boxOrientation));
#endif // DEBUG

            VECTOR InsideAll = Vector::TrueInt();

            for (size_t i = 0; i < BoundingOrientedBox::CORNER_COUNT; ++i)
            {
                VECTOR C = Vector::Add(Vector3::Rotate(Vector::Multiply(boxExtents, g_BoxOffset[i]), boxOrientation), boxCenter);
                VECTOR d = Vector3::LengthSq(Vector::Subtract(vCenter, C));
                InsideAll = Vector::AndInt(InsideAll, Vector::LessOrEqual(d, RadiusSq));
            }

            return (Vector3::EqualInt(InsideAll, Vector::TrueInt())) ? CONTAINS : INTERSECTS;
        }

        _Use_decl_annotations_
        FORCE_INLINE ContainmentType BoundingSphere::Contains(const BoundingFrustum& fr) const noexcept
        {
            if (!fr.Intersects(*this))
                return DISJOINT;

            VECTOR vCenter = Vector::LoadFloat3(&center);
            VECTOR vRadius = Vector::ReplicatePtr(&radius);
            VECTOR RadiusSq = Vector::Multiply(vRadius, vRadius);

            VECTOR vOrigin = Vector::LoadFloat3(&fr.origin);
            VECTOR vOrientation = Vector::LoadFloat4(&fr.orientation);

#if defined(DEBUG) || defined(_DEBUG)
            assert(QuaternionIsUnit(vOrientation));
#endif // DEBUG

            // Build the corners of the frustum.
            VECTOR vRightTop = Vector::Set(fr.rightSlope, fr.topSlope, 1.0f, 0.0f);
            VECTOR vRightBottom = Vector::Set(fr.rightSlope, fr.bottomSlope, 1.0f, 0.0f);
            VECTOR vLeftTop = Vector::Set(fr.leftSlope, fr.topSlope, 1.0f, 0.0f);
            VECTOR vLeftBottom = Vector::Set(fr.leftSlope, fr.bottomSlope, 1.0f, 0.0f);
            VECTOR vNear = Vector::ReplicatePtr(&fr.near);
            VECTOR vFar = Vector::ReplicatePtr(&fr.far);

            VECTOR Corners[BoundingFrustum::CORNER_COUNT];
            Corners[0] = Vector::Multiply(vRightTop, vNear);
            Corners[1] = Vector::Multiply(vRightBottom, vNear);
            Corners[2] = Vector::Multiply(vLeftTop, vNear);
            Corners[3] = Vector::Multiply(vLeftBottom, vNear);
            Corners[4] = Vector::Multiply(vRightTop, vFar);
            Corners[5] = Vector::Multiply(vRightBottom, vFar);
            Corners[6] = Vector::Multiply(vLeftTop, vFar);
            Corners[7] = Vector::Multiply(vLeftBottom, vFar);

            VECTOR InsideAll = Vector::TrueInt();
            for (size_t i = 0; i < BoundingFrustum::CORNER_COUNT; ++i)
            {
                VECTOR C = Vector::Add(Vector3::Rotate(Corners[i], vOrientation), vOrigin);
                VECTOR d = Vector3::LengthSq(Vector::Subtract(vCenter, C));
                InsideAll = Vector::AndInt(InsideAll, Vector::LessOrEqual(d, RadiusSq));
            }

            return (Vector3::EqualInt(InsideAll, Vector::TrueInt())) ? CONTAINS : INTERSECTS;
        }

        _Use_decl_annotations_
        FORCE_INLINE bool BoundingSphere::Intersects(const BoundingSphere& sh) const noexcept
        {
            // Load A.
            VECTOR vCenterA = Vector::LoadFloat3(&center);
            VECTOR vRadiusA = Vector::ReplicatePtr(&radius);

            // Load B.
            VECTOR vCenterB = Vector::LoadFloat3(&sh.center);
            VECTOR vRadiusB = Vector::ReplicatePtr(&sh.radius);

            // Distance squared between centers.
            VECTOR Delta = Vector::Subtract(vCenterB, vCenterA);
            VECTOR DistanceSquared = Vector3::LengthSq(Delta);

            // Sum of the radii squared.
            VECTOR RadiusSquared = Vector::Add(vRadiusA, vRadiusB);
            RadiusSquared = Vector::Multiply(RadiusSquared, RadiusSquared);

            return Vector3::LessOrEqual(DistanceSquared, RadiusSquared);
        }

        _Use_decl_annotations_
        FORCE_INLINE bool BoundingSphere::Intersects(const BoundingBox& box) const noexcept
        {
            return box.Intersects(*this);
        }

        _Use_decl_annotations_
        FORCE_INLINE bool BoundingSphere::Intersects(const BoundingOrientedBox& box) const noexcept
        {
            return box.Intersects(*this);
        }

        _Use_decl_annotations_
        inline bool BoundingSphere::Intersects(const BoundingFrustum& fr) const noexcept
        {
            return fr.Intersects(*this);
        }

        _Use_decl_annotations_
        FORCE_INLINE bool VEC_CALLCONV BoundingSphere::Intersects(A_VECTOR V0, A_VECTOR V1, A_VECTOR V2) const noexcept
        {
            // Load the sphere.
            VECTOR vCenter = Vector::LoadFloat3(&center);
            VECTOR vRadius = Vector::ReplicatePtr(&radius);

            // Compute the plane of the triangle (has to be normalized).
            VECTOR N = Vector3::Normalize(Vector3::Cross(Vector::Subtract(V1, V0), Vector::Subtract(V2, V0)));

#if defined(DEBUG) || defined(_DEBUG)
            // Assert that the triangle is not degenerate.
            assert(!Vector3::Equal(N, Vector::Zero()));
#endif // DEBUG

            // Find the nearest feature on the triangle to the sphere.
            VECTOR Dist = Vector3::Dot(Vector::Subtract(vCenter, V0), N);

            // If the center of the sphere is farther from the plane of the triangle than
            // the radius of the sphere, then there cannot be an intersection.
            VECTOR NoIntersection = Vector::Less(Dist, Vector::Negate(vRadius));
            NoIntersection = Vector::OrInt(NoIntersection, Vector::Greater(Dist, vRadius));

            // Project the center of the sphere onto the plane of the triangle.
            VECTOR Point = Vector::NegativeMultiplySubtract(N, Dist, vCenter);

            // Is it inside all the edges? If so we intersect because the distance
            // to the plane is less than the radius.
            VECTOR Intersection = PointOnPlaneInsideTriangle(Point, V0, V1, V2);

            // Find the nearest point on each edge.
            VECTOR RadiusSq = Vector::Multiply(vRadius, vRadius);

            // Edge 0,1
            Point = PointOnLineSegmentNearestPoint(V0, V1, vCenter);

            // If the distance to the center of the sphere to the point is less than
            // the radius of the sphere then it must intersect.
            Intersection = Vector::OrInt(Intersection, Vector::LessOrEqual(Vector3::LengthSq(Vector::Subtract(vCenter, Point)), RadiusSq));

            // Edge 1,2
            Point = PointOnLineSegmentNearestPoint(V1, V2, vCenter);

            // If the distance to the center of the sphere to the point is less than
            // the radius of the sphere then it must intersect.
            Intersection = Vector::OrInt(Intersection, Vector::LessOrEqual(Vector3::LengthSq(Vector::Subtract(vCenter, Point)), RadiusSq));

            // Edge 2,0
            Point = PointOnLineSegmentNearestPoint(V2, V0, vCenter);

            // If the distance to the center of the sphere to the point is less than
            // the radius of the sphere then it must intersect.
            Intersection = Vector::OrInt(Intersection, Vector::LessOrEqual(Vector3::LengthSq(Vector::Subtract(vCenter, Point)), RadiusSq));

            return Vector4::EqualInt(Vector::AndCInt(Intersection, NoIntersection), Vector::TrueInt());
        }

        _Use_decl_annotations_
        FORCE_INLINE PlaneIntersectionType VEC_CALLCONV BoundingSphere::Intersects(A_VECTOR plane) const noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(PlaneIsUnit(plane));
#endif // DEBUG

            // Load the sphere.
            VECTOR vCenter = Vector::LoadFloat3(&center);
            VECTOR vRadius = Vector::ReplicatePtr(&radius);

            // Set w of the center to one so we can dot4 with a plane.
            vCenter = Vector::Insert<0, 0, 0, 0, 1>(vCenter, Vector::SplatOne());

            VECTOR outside, inside;
            FastIntersectSpherePlane(vCenter, vRadius, plane, outside, inside);

            // If the sphere is outside any plane it is outside.
            if (Vector4::EqualInt(outside, Vector::TrueInt()))
                return FRONT;

            // If the sphere is inside all planes it is inside.
            if (Vector4::EqualInt(inside, Vector::TrueInt()))
                return BACK;

            // The sphere is not inside all planes or outside a plane it intersects.
            return INTERSECTING;
        }

        _Use_decl_annotations_
        FORCE_INLINE bool VEC_CALLCONV BoundingSphere::Intersects(A_VECTOR origin, A_VECTOR direction, float& dist) const noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(Vector3::IsUnit(direction));
#endif // DEBUG

            VECTOR vCenter = Vector::LoadFloat3(&center);
            VECTOR vRadius = Vector::ReplicatePtr(&radius);

            // l is the vector from the ray origin to the center of the sphere.
            VECTOR l = Vector::Subtract(vCenter, origin);

            // s is the projection of the l onto the ray direction.
            VECTOR s = Vector3::Dot(l, direction);

            VECTOR l2 = Vector3::Dot(l, l);

            VECTOR r2 = Vector::Multiply(vRadius, vRadius);

            // m2 is squared distance from the center of the sphere to the projection.
            VECTOR m2 = Vector::NegativeMultiplySubtract(s, s, l2);

            VECTOR NoIntersection;

            // If the ray origin is outside the sphere and the center of the sphere is
            // behind the ray origin there is no intersection.
            NoIntersection = Vector::AndInt(Vector::Less(s, Vector::Zero()), Vector::Greater(l2, r2));

            // If the squared distance from the center of the sphere to the projection
            // is greater than the radius squared the ray will miss the sphere.
            NoIntersection = Vector::OrInt(NoIntersection, Vector::Greater(m2, r2));

            // The ray hits the sphere, compute the nearest intersection point.
            VECTOR q = Vector::Sqrt(Vector::Subtract(r2, m2));
            VECTOR t1 = Vector::Subtract(s, q);
            VECTOR t2 = Vector::Add(s, q);

            VECTOR OriginInside = Vector::LessOrEqual(l2, r2);
            VECTOR t = Vector::Select(t1, t2, OriginInside);

            if (Vector4::NotEqualInt(NoIntersection, Vector::TrueInt()))
            {
                // Store the x-component to *pDist.
                Vector::StoreFloat(&dist, t);
                return true;
            }

            dist = 0.0f;
            return false;
        }

        _Use_decl_annotations_
        FORCE_INLINE ContainmentType VEC_CALLCONV BoundingSphere::ContainedBy(
            A_VECTOR plane0, A_VECTOR plane1, A_VECTOR plane2, 
            B_VECTOR plane3, 
            C_VECTOR plane4, C_VECTOR plane5) const noexcept
        {
            // Load the sphere.
            VECTOR vCenter = Vector::LoadFloat3(&center);
            VECTOR vRadius = Vector::ReplicatePtr(&radius);

            // Set w of the center to one so we can dot4 with a plane.
            vCenter = Vector::Insert<0, 0, 0, 0, 1>(vCenter, Vector::SplatOne());

            VECTOR outside, inside;

            // Test against each plane.
            FastIntersectSpherePlane(vCenter, vRadius, plane0, outside, inside);

            VECTOR AnyOutside = outside;
            VECTOR AllInside = inside;

            FastIntersectSpherePlane(vCenter, vRadius, plane1, outside, inside);
            AnyOutside = Vector::OrInt(AnyOutside, outside);
            AllInside = Vector::AndInt(AllInside, inside);

            FastIntersectSpherePlane(vCenter, vRadius, plane2, outside, inside);
            AnyOutside = Vector::OrInt(AnyOutside, outside);
            AllInside = Vector::AndInt(AllInside, inside);

            FastIntersectSpherePlane(vCenter, vRadius, plane3, outside, inside);
            AnyOutside = Vector::OrInt(AnyOutside, outside);
            AllInside = Vector::AndInt(AllInside, inside);

            FastIntersectSpherePlane(vCenter, vRadius, plane4, outside, inside);
            AnyOutside = Vector::OrInt(AnyOutside, outside);
            AllInside = Vector::AndInt(AllInside, inside);

            FastIntersectSpherePlane(vCenter, vRadius, plane5, outside, inside);
            AnyOutside = Vector::OrInt(AnyOutside, outside);
            AllInside = Vector::AndInt(AllInside, inside);

            // If the sphere is outside any plane it is outside.
            if (Vector4::EqualInt(AnyOutside, Vector::TrueInt()))
                return DISJOINT;

            // If the sphere is inside all planes it is inside.
            if (Vector4::EqualInt(AllInside, Vector::TrueInt()))
                return CONTAINS;

            // The sphere is not inside all planes or outside a plane, it may intersect.
            return INTERSECTS;
        }

        _Use_decl_annotations_
        FORCE_INLINE void BoundingSphere::CreateMerged(BoundingSphere& out, const BoundingSphere& S1, const BoundingSphere& S2) noexcept
        {
            VECTOR Center1 = Vector::LoadFloat3(&S1.center);
            float r1 = S1.radius;

            VECTOR Center2 = Vector::LoadFloat3(&S2.center);
            float r2 = S2.radius;

            VECTOR V = Vector::Subtract(Center2, Center1);

            VECTOR Dist = Vector3::Length(V);

            float d = Vector::GetX(Dist);

            if (r1 + r2 >= d)
            {
                if (r1 - r2 >= d)
                {
                    out = S1;
                    return;
                }
                else if (r2 - r1 >= d)
                {
                    out = S2;
                    return;
                }
            }

            VECTOR N = Vector::Divide(V, Dist);

            float t1 = Min(-r1, d - r2);
            float t2 = Max(r1, d + r2);
            float t_5 = (t2 - t1) * 0.5f;

            VECTOR NCenter = Vector::Add(Center1, Vector::Multiply(N, Vector::Replicate(t_5 + t1)));

            Vector::StoreFloat3(&out.center, NCenter);
            out.radius = t_5;
        }

        _Use_decl_annotations_
        FORCE_INLINE void BoundingSphere::CreateFromBoundingBox(BoundingSphere& out, const BoundingBox& box) noexcept
        {
            out.center = box.center;
            VECTOR vExtents = Vector::LoadFloat3(&box.extents);
            out.radius = Vector::GetX(Vector3::Length(vExtents));
        }

        _Use_decl_annotations_
        FORCE_INLINE void BoundingSphere::CreateFromBoundingBox(BoundingSphere& out, const BoundingOrientedBox& box) noexcept
        {
            // Bounding box orientation is irrelevant because a sphere is rotationally invariant
            out.center = box.center;
            VECTOR vExtents = Vector::LoadFloat3(&box.extents);
            out.radius = Vector::GetX(Vector3::Length(vExtents));
        }

        _Use_decl_annotations_
        FORCE_INLINE void BoundingSphere::CreateFromPoints(BoundingSphere& out, size_t count, const Float3* pPoints, size_t stride) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(count > 0);
            assert(pPoints);
#endif // DEBUG

            // Find the points with minimum and maximum x, y, and z
            VECTOR MinX, MaxX, MinY, MaxY, MinZ, MaxZ;

            MinX = MaxX = MinY = MaxY = MinZ = MaxZ = Vector::LoadFloat3(pPoints);

            for (size_t i = 1; i < count; ++i)
            {
                VECTOR Point = Vector::LoadFloat3(reinterpret_cast<const Float3*>(reinterpret_cast<const uint8_t*>(pPoints) + i * stride));

                float px = Vector::GetX(Point);
                float py = Vector::GetY(Point);
                float pz = Vector::GetZ(Point);

                if (px < Vector::GetX(MinX))
                    MinX = Point;

                if (px > Vector::GetX(MaxX))
                    MaxX = Point;

                if (py < Vector::GetY(MinY))
                    MinY = Point;

                if (py > Vector::GetY(MaxY))
                    MaxY = Point;

                if (pz < Vector::GetZ(MinZ))
                    MinZ = Point;

                if (pz > Vector::GetZ(MaxZ))
                    MaxZ = Point;
            }

            // Use the min/max pair that are farthest apart to form the initial sphere.
            VECTOR DeltaX = Vector::Subtract(MaxX, MinX);
            VECTOR DistX = Vector3::Length(DeltaX);

            VECTOR DeltaY = Vector::Subtract(MaxY, MinY);
            VECTOR DistY = Vector3::Length(DeltaY);

            VECTOR DeltaZ = Vector::Subtract(MaxZ, MinZ);
            VECTOR DistZ = Vector3::Length(DeltaZ);

            VECTOR vCenter;
            VECTOR vRadius;

            if (Vector3::Greater(DistX, DistY))
            {
                if (Vector3::Greater(DistX, DistZ))
                {
                    // Use min/max x.
                    vCenter = Vector::Lerp(MaxX, MinX, 0.5f);
                    vRadius = Vector::Scale(DistX, 0.5f);
                }
                else
                {
                    // Use min/max z.
                    vCenter = Vector::Lerp(MaxZ, MinZ, 0.5f);
                    vRadius = Vector::Scale(DistZ, 0.5f);
                }
            }
            else // Y >= X
            {
                if (Vector3::Greater(DistY, DistZ))
                {
                    // Use min/max y.
                    vCenter = Vector::Lerp(MaxY, MinY, 0.5f);
                    vRadius = Vector::Scale(DistY, 0.5f);
                }
                else
                {
                    // Use min/max z.
                    vCenter = Vector::Lerp(MaxZ, MinZ, 0.5f);
                    vRadius = Vector::Scale(DistZ, 0.5f);
                }
            }

            // Add any points not inside the sphere.
            for (size_t i = 0; i < count; ++i)
            {
                VECTOR Point = Vector::LoadFloat3(reinterpret_cast<const Float3*>(reinterpret_cast<const uint8_t*>(pPoints) + i * stride));

                VECTOR Delta = Vector::Subtract(Point, vCenter);

                VECTOR Dist = Vector3::Length(Delta);

                if (Vector3::Greater(Dist, vRadius))
                {
                    // Adjust sphere to include the new point.
                    vRadius = Vector::Scale(Vector::Add(vRadius, Dist), 0.5f);
                    vCenter = Vector::Add(vCenter, Vector::Multiply(Vector::Subtract(Vector::Replicate(1.0f), Vector::Divide(vRadius, Dist)), Delta));
                }
            }

            Vector::StoreFloat3(&out.center, vCenter);
            Vector::StoreFloat(&out.radius, vRadius);
        }

        _Use_decl_annotations_
        FORCE_INLINE void BoundingSphere::CreateFromFrustum(BoundingSphere& out, const BoundingFrustum& fr) noexcept
        {
            Float3 Corners[BoundingFrustum::CORNER_COUNT];
            fr.GetCorners(Corners);
            CreateFromPoints(out, BoundingFrustum::CORNER_COUNT, Corners, sizeof(Float3));
        }

        _Use_decl_annotations_
        FORCE_INLINE void VEC_CALLCONV BoundingBox::Transform(BoundingBox& out, A_MATRIX m) const noexcept
        {
            // Load center and extents.
            VECTOR vCenter = Vector::LoadFloat3(&center);
            VECTOR vExtents = Vector::LoadFloat3(&extents);

            // Compute and transform the corners and find new min/max bounds.
            VECTOR Corner = Vector::MultiplyAdd(vExtents, g_BoxOffset[0], vCenter);
            Corner = Vector3::Transform(Corner, M);

            VECTOR Min, Max;
            Min = Max = Corner;

            for (size_t i = 1; i < CORNER_COUNT; ++i)
            {
                Corner = Vector::MultiplyAdd(vExtents, g_BoxOffset[i], vCenter);
                Corner = Vector3::Transform(Corner, M);

                Min = Vector::Min(Min, Corner);
                Max = Vector::Max(Max, Corner);
            }

            // Store center and extents.
            Vector::StoreFloat3(&out.center, Vector::Scale(Vector::Add(Min, Max), 0.5f));
            Vector::StoreFloat3(&out.extents, Vector::Scale(Vector::Subtract(Max, Min), 0.5f));
        }

        _Use_decl_annotations_
        FORCE_INLINE void VEC_CALLCONV BoundingBox::Transform(BoundingBox& out, float scale, A_VECTOR rotation, A_VECTOR translation) const noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(QuaternionIsUnit(rotation));
#endif // DEBUG

            // Load center and extents.
            VECTOR vCenter = Vector::LoadFloat3(&Ccnter);
            VECTOR vExtents = Vector::LoadFloat3(&extents);

            VECTOR VectorScale = Vector::Replicate(scale);

            // Compute and transform the corners and find new min/max bounds.
            VECTOR Corner = Vector::MultiplyAdd(vExtents, g_BoxOffset[0], vCenter);
            Corner = Vector::Add(Vector3::Rotate(Vector::Multiply(Corner, VectorScale), rotation), translation);

            VECTOR Min, Max;
            Min = Max = Corner;

            for (size_t i = 1; i < CORNER_COUNT; ++i)
            {
                Corner = Vector::MultiplyAdd(vExtents, g_BoxOffset[i], vCenter);
                Corner = Vector::Add(Vector3::Rotate(Vector::Multiply(Corner, VectorScale), rotation), translation);

                Min = Vector::Min(Min, Corner);
                Max = Vector::Max(Max, Corner);
            }

            // Store center and extents.
            Vector::StoreFloat3(&out.center, Vector::Scale(Vector::Add(Min, Max), 0.5f));
            Vector::StoreFloat3(&out.extents, Vector::Scale(Vector::Subtract(Max, Min), 0.5f));
        }

        _Use_decl_annotations_
        FORCE_INLINE void BoundingBox::GetCorners(Float3* corners) const noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(corners != nullptr);
#endif // DEBUG

            // Load the box
            VECTOR vCenter = Vector::LoadFloat3(&center);
            VECTOR vExtents = Vector::LoadFloat3(&extents);

            for (size_t i = 0; i < CORNER_COUNT; ++i)
            {
                VECTOR C = Vector::MultiplyAdd(vExtents, g_BoxOffset[i], vCenter);
                Vector::StoreFloat3(&corners[i], C);
            }
        }

        _Use_decl_annotations_
        FORCE_INLINE ContainmentType VEC_CALLCONV BoundingBox::Contains(A_VECTOR point) const noexcept
        {
            VECTOR vCenter = Vector::LoadFloat3(&center);
            VECTOR vExtents = Vector::LoadFloat3(&extents);

            return Vector3::InBounds(Vector::Subtract(point, vCenter), vExtents) ? CONTAINS : DISJOINT;
        }

        _Use_decl_annotations_
        FORCE_INLINE ContainmentType VEC_CALLCONV BoundingBox::Contains(A_VECTOR V0, A_VECTOR V1, A_VECTOR V2) const noexcept
        {
            if (!Intersects(V0, V1, V2))
                return DISJOINT;

            VECTOR vCenter = Vector::LoadFloat3(&center);
            VECTOR vExtents = Vector::LoadFloat3(&extents);

            VECTOR d = Vector::Abs(Vector::Subtract(V0, vCenter));
            VECTOR Inside = Vector::LessOrEqual(d, vExtents);

            d = Vector::Abs(Vector::Subtract(V1, vCenter));
            Inside = Vector::AndInt(Inside, Vector::LessOrEqual(d, vExtents));

            d = Vector::Abs(Vector::Subtract(V2, vCenter));
            Inside = Vector::AndInt(Inside, Vector::LessOrEqual(d, vExtents));

            return (Vector3::EqualInt(Inside, Vector::TrueInt())) ? CONTAINS : INTERSECTS;
        }

        _Use_decl_annotations_
        FORCE_INLINE ContainmentType BoundingBox::Contains(const BoundingSphere& sh) const noexcept
        {
            VECTOR SphereCenter = Vector::LoadFloat3(&sh.center);
            VECTOR SphereRadius = Vector::ReplicatePtr(&sh.radius);

            VECTOR BoxCenter = Vector::LoadFloat3(&center);
            VECTOR BoxExtents = Vector::LoadFloat3(&extents);

            VECTOR BoxMin = Vector::Subtract(BoxCenter, BoxExtents);
            VECTOR BoxMax = Vector::Add(BoxCenter, BoxExtents);

            // Find the distance to the nearest point on the box.
            // for each i in (x, y, z)
            // if (SphereCenter(i) < BoxMin(i)) d2 += (SphereCenter(i) - BoxMin(i)) ^ 2
            // else if (SphereCenter(i) > BoxMax(i)) d2 += (SphereCenter(i) - BoxMax(i)) ^ 2

            VECTOR d = Vector::Zero();

            // Compute d for each dimension.
            VECTOR LessThanMin = Vector::Less(SphereCenter, BoxMin);
            VECTOR GreaterThanMax = Vector::Greater(SphereCenter, BoxMax);

            VECTOR MinDelta = Vector::Subtract(SphereCenter, BoxMin);
            VECTOR MaxDelta = Vector::Subtract(SphereCenter, BoxMax);

            // Choose value for each dimension based on the comparison.
            d = Vector::Select(d, MinDelta, LessThanMin);
            d = Vector::Select(d, MaxDelta, GreaterThanMax);

            // Use a dot-product to square them and sum them together.
            VECTOR d2 = Vector3::Dot(d, d);

            if (Vector3::Greater(d2, Vector::Multiply(SphereRadius, SphereRadius)))
                return DISJOINT;

            VECTOR InsideAll = Vector::LessOrEqual(Vector::Add(BoxMin, SphereRadius), SphereCenter);
            InsideAll = Vector::AndInt(InsideAll, Vector::LessOrEqual(SphereCenter, Vector::Subtract(BoxMax, SphereRadius)));
            InsideAll = Vector::AndInt(InsideAll, Vector::Greater(Vector::Subtract(BoxMax, BoxMin), SphereRadius));

            return (Vector3::EqualInt(InsideAll, Vector::TrueInt())) ? CONTAINS : INTERSECTS;
        }

        _Use_decl_annotations_
        FORCE_INLINE ContainmentType BoundingBox::Contains(const BoundingBox& box) const noexcept
        {
            VECTOR CenterA = Vector::LoadFloat3(&center);
            VECTOR ExtentsA = Vector::LoadFloat3(&extents);

            VECTOR CenterB = Vector::LoadFloat3(&box.center);
            VECTOR ExtentsB = Vector::LoadFloat3(&box.extents);

            VECTOR MinA = Vector::Subtract(CenterA, ExtentsA);
            VECTOR MaxA = Vector::Add(CenterA, ExtentsA);

            VECTOR MinB = Vector::Subtract(CenterB, ExtentsB);
            VECTOR MaxB = Vector::Add(CenterB, ExtentsB);

            // for each i in (x, y, z) if a_min(i) > b_max(i) or b_min(i) > a_max(i) then return false
            VECTOR Disjoint = Vector::OrInt(Vector::Greater(MinA, MaxB), Vector::Greater(MinB, MaxA));

            if (Vector3AnyTrue(Disjoint))
                return DISJOINT;

            // for each i in (x, y, z) if a_min(i) <= b_min(i) and b_max(i) <= a_max(i) then A contains B
            VECTOR Inside = Vector::AndInt(Vector::LessOrEqual(MinA, MinB), Vector::LessOrEqual(MaxB, MaxA));

            return Vector3AllTrue(Inside) ? CONTAINS : INTERSECTS;
        }

        _Use_decl_annotations_
        FORCE_INLINE ContainmentType BoundingBox::Contains(const BoundingOrientedBox& box) const noexcept
        {
            if (!box.Intersects(*this))
                return DISJOINT;

            VECTOR vCenter = Vector::LoadFloat3(&center);
            VECTOR vExtents = Vector::LoadFloat3(&extents);

            // Subtract off the AABB center to remove a subtract below
            VECTOR oCenter = Vector::Subtract(Vector::LoadFloat3(&box.center), vCenter);

            VECTOR oExtents = Vector::LoadFloat3(&box.extents);
            VECTOR oOrientation = Vector::LoadFloat4(&box.orientation);

#if defined(DEBUG) || defined(_DEBUG)
            assert(QuaternionIsUnit(oOrientation));
#endif // DEBUG

            VECTOR Inside = Vector::TrueInt();

            for (size_t i = 0; i < BoundingOrientedBox::CORNER_COUNT; ++i)
            {
                VECTOR C = Vector::Add(Vector3::Rotate(Vector::Multiply(oExtents, g_BoxOffset[i]), oOrientation), oCenter);
                VECTOR d = Vector::Abs(C);
                Inside = Vector::AndInt(Inside, Vector::LessOrEqual(d, vExtents));
            }

            return (Vector3::EqualInt(Inside, Vector::TrueInt())) ? CONTAINS : INTERSECTS;
        }

        _Use_decl_annotations_
        FORCE_INLINE ContainmentType BoundingBox::Contains(const BoundingFrustum& fr) const noexcept
        {
            if (!fr.Intersects(*this))
                return DISJOINT;

            Float3 Corners[BoundingFrustum::CORNER_COUNT];
            fr.GetCorners(Corners);

            VECTOR vCenter = Vector::LoadFloat3(&center);
            VECTOR vExtents = Vector::LoadFloat3(&extents);

            VECTOR Inside = Vector::TrueInt();

            for (size_t i = 0; i < BoundingFrustum::CORNER_COUNT; ++i)
            {
                VECTOR Point = Vector::LoadFloat3(&Corners[i]);
                VECTOR d = Vector::Abs(Vector::Subtract(Point, vCenter));
                Inside = Vector::AndInt(Inside, Vector::LessOrEqual(d, vExtents));
            }

            return (Vector3::EqualInt(Inside, Vector::TrueInt())) ? CONTAINS : INTERSECTS;
        }

        _Use_decl_annotations_
        FORCE_INLINE bool BoundingBox::Intersects(const BoundingSphere& sh) const noexcept
        {
            VECTOR SphereCenter = Vector::LoadFloat3(&sh.center);
            VECTOR SphereRadius = Vector::ReplicatePtr(&sh.radius);

            VECTOR BoxCenter = Vector::LoadFloat3(&center);
            VECTOR BoxExtents = Vector::LoadFloat3(&extents);

            VECTOR BoxMin = Vector::Subtract(BoxCenter, BoxExtents);
            VECTOR BoxMax = Vector::Add(BoxCenter, BoxExtents);

            // Find the distance to the nearest point on the box.
            // for each i in (x, y, z)
            // if (SphereCenter(i) < BoxMin(i)) d2 += (SphereCenter(i) - BoxMin(i)) ^ 2
            // else if (SphereCenter(i) > BoxMax(i)) d2 += (SphereCenter(i) - BoxMax(i)) ^ 2

            VECTOR d = Vector::Zero();

            // Compute d for each dimension.
            VECTOR LessThanMin = Vector::Less(SphereCenter, BoxMin);
            VECTOR GreaterThanMax = Vector::Greater(SphereCenter, BoxMax);

            VECTOR MinDelta = Vector::Subtract(SphereCenter, BoxMin);
            VECTOR MaxDelta = Vector::Subtract(SphereCenter, BoxMax);

            // Choose value for each dimension based on the comparison.
            d = Vector::Select(d, MinDelta, LessThanMin);
            d = Vector::Select(d, MaxDelta, GreaterThanMax);

            // Use a dot-product to square them and sum them together.
            VECTOR d2 = Vector3::Dot(d, d);

            return Vector3::LessOrEqual(d2, Vector::Multiply(SphereRadius, SphereRadius));
        }

        _Use_decl_annotations_
        FORCE_INLINE bool BoundingBox::Intersects(const BoundingBox& box) const noexcept
        {
            VECTOR CenterA = Vector::LoadFloat3(&center);
            VECTOR ExtentsA = Vector::LoadFloat3(&extents);

            VECTOR CenterB = Vector::LoadFloat3(&box.center);
            VECTOR ExtentsB = Vector::LoadFloat3(&box.extents);

            VECTOR MinA = Vector::Subtract(CenterA, ExtentsA);
            VECTOR MaxA = Vector::Add(CenterA, ExtentsA);

            VECTOR MinB = Vector::Subtract(CenterB, ExtentsB);
            VECTOR MaxB = Vector::Add(CenterB, ExtentsB);

            // for each i in (x, y, z) if a_min(i) > b_max(i) or b_min(i) > a_max(i) then return false
            VECTOR Disjoint = Vector::OrInt(Vector::Greater(MinA, MaxB), Vector::Greater(MinB, MaxA));

            return !Vector3AnyTrue(Disjoint);
        }

        _Use_decl_annotations_
        FORCE_INLINE bool BoundingBox::Intersects(const BoundingOrientedBox& box) const noexcept
        {
            return box.Intersects(*this);
        }

        _Use_decl_annotations_
        FORCE_INLINE bool BoundingBox::Intersects(const BoundingFrustum& fr) const noexcept
        {
            return fr.Intersects(*this);
        }

        _Use_decl_annotations_
        FORCE_INLINE bool VEC_CALLCONV BoundingBox::Intersects(A_VECTOR V0, A_VECTOR V1, A_VECTOR V2) const noexcept
        {
            VECTOR zero = Vector::Zero();

            // Load the box.
            VECTOR vCenter = Vector::LoadFloat3(&center);
            VECTOR vExtents = Vector::LoadFloat3(&extents);

            VECTOR BoxMin = Vector::Subtract(vCenter, vExtents);
            VECTOR BoxMax = Vector::Add(vCenter, vExtents);

            // Test the axes of the box (in effect test the AAB against the minimal AAB
            // around the triangle).
            VECTOR TriMin = Vector::Min(Vector::Min(V0, V1), V2);
            VECTOR TriMax = Vector::Max(Vector::Max(V0, V1), V2);

            // for each i in (x, y, z) if a_min(i) > b_max(i) or b_min(i) > a_max(i) then disjoint
            VECTOR Disjoint = Vector::OrInt(Vector::Greater(TriMin, BoxMax), Vector::Greater(BoxMin, TriMax));
            if (Vector3AnyTrue(Disjoint))
                return false;

            // Test the plane of the triangle.
            VECTOR Normal = Vector3::Cross(Vector::Subtract(V1, V0), Vector::Subtract(V2, V0));
            VECTOR Dist = Vector3::Dot(Normal, V0);

#if defined(DEBUG) || defined(_DEBUG)
            // Assert that the triangle is not degenerate.
            assert(!Vector3::Equal(Normal, zero));
#endif // DEBUG

            // for each i in (x, y, z) if n(i) >= 0 then v_min(i)=b_min(i), v_max(i)=b_max(i)
            // else v_min(i)=b_max(i), v_max(i)=b_min(i)
            VECTOR NormalSelect = Vector::Greater(Normal, zero);
            VECTOR V_Min = Vector::Select(BoxMax, BoxMin, NormalSelect);
            VECTOR V_Max = Vector::Select(BoxMin, BoxMax, NormalSelect);

            // if n dot v_min + d > 0 || n dot v_max + d < 0 then disjoint
            VECTOR MinDist = Vector3::Dot(V_Min, Normal);
            VECTOR MaxDist = Vector3::Dot(V_Max, Normal);

            VECTOR NoIntersection = Vector::Greater(MinDist, Dist);
            NoIntersection = Vector::OrInt(NoIntersection, Vector::Less(MaxDist, Dist));

            // Move the box center to zero to simplify the following tests.
            VECTOR TV0 = Vector::Subtract(V0, vCenter);
            VECTOR TV1 = Vector::Subtract(V1, vCenter);
            VECTOR TV2 = Vector::Subtract(V2, vCenter);

            // Test the edge/edge axes (3*3).
            VECTOR e0 = Vector::Subtract(TV1, TV0);
            VECTOR e1 = Vector::Subtract(TV2, TV1);
            VECTOR e2 = Vector::Subtract(TV0, TV2);

            // Make w zero.
            e0 = Vector::Insert<0, 0, 0, 0, 1>(e0, zero);
            e1 = Vector::Insert<0, 0, 0, 0, 1>(e1, zero);
            e2 = Vector::Insert<0, 0, 0, 0, 1>(e2, zero);

            VECTOR Axis;
            VECTOR p0, p1, p2;
            VECTOR Min, Max;
            VECTOR Radius;

            // Axis == (1,0,0) x e0 = (0, -e0.z, e0.y)
            Axis = Vector::Permute<PERMUTE_0W, PERMUTE_1Z, PERMUTE_0Y, PERMUTE_0X>(e0, Vector::Negate(e0));
            p0 = Vector3::Dot(TV0, Axis);
            // p1 = XMVector3Dot( V1, Axis ); // p1 = p0;
            p2 = Vector3::Dot(TV2, Axis);
            Min = Vector::Min(p0, p2);
            Max = Vector::Max(p0, p2);
            Radius = Vector3::Dot(vExtents, Vector::Abs(Axis));
            NoIntersection = Vector::OrInt(NoIntersection, Vector::Greater(Min, Radius));
            NoIntersection = Vector::OrInt(NoIntersection, Vector::Less(Max, Vector::Negate(Radius)));

            // Axis == (1,0,0) x e1 = (0, -e1.z, e1.y)
            Axis = Vector::Permute<PERMUTE_0W, PERMUTE_1Z, PERMUTE_0Y, PERMUTE_0X>(e1, Vector::Negate(e1));
            p0 = Vector3::Dot(TV0, Axis);
            p1 = Vector3::Dot(TV1, Axis);
            // p2 = XMVector3Dot( V2, Axis ); // p2 = p1;
            Min = Vector::Min(p0, p1);
            Max = Vector::Max(p0, p1);
            Radius = Vector3::Dot(vExtents, Vector::Abs(Axis));
            NoIntersection = Vector::OrInt(NoIntersection, Vector::Greater(Min, Radius));
            NoIntersection = Vector::OrInt(NoIntersection, Vector::Less(Max, Vector::Negate(Radius)));

            // Axis == (1,0,0) x e2 = (0, -e2.z, e2.y)
            Axis = Vector::Permute<PERMUTE_0W, PERMUTE_1Z, PERMUTE_0Y, PERMUTE_0X>(e2, Vector::Negate(e2));
            p0 = Vector3::Dot(TV0, Axis);
            p1 = Vector3::Dot(TV1, Axis);
            // p2 = XMVector3Dot( V2, Axis ); // p2 = p0;
            Min = Vector::Min(p0, p1);
            Max = Vector::Max(p0, p1);
            Radius = Vector3::Dot(vExtents, Vector::Abs(Axis));
            NoIntersection = Vector::OrInt(NoIntersection, Vector::Greater(Min, Radius));
            NoIntersection = Vector::OrInt(NoIntersection, Vector::Less(Max, Vector::Negate(Radius)));

            // Axis == (0,1,0) x e0 = (e0.z, 0, -e0.x)
            Axis = Vector::Permute<PERMUTE_0Z, PERMUTE_0W, PERMUTE_1X, PERMUTE_0Y>(e0, Vector::Negate(e0));
            p0 = Vector3::Dot(TV0, Axis);
            // p1 = XMVector3Dot( V1, Axis ); // p1 = p0;
            p2 = Vector3::Dot(TV2, Axis);
            Min = Vector::Min(p0, p2);
            Max = Vector::Max(p0, p2);
            Radius = Vector3::Dot(vExtents, Vector::Abs(Axis));
            NoIntersection = Vector::OrInt(NoIntersection, Vector::Greater(Min, Radius));
            NoIntersection = Vector::OrInt(NoIntersection, Vector::Less(Max, Vector::Negate(Radius)));

            // Axis == (0,1,0) x e1 = (e1.z, 0, -e1.x)
            Axis = Vector::Permute<PERMUTE_0Z, PERMUTE_0W, PERMUTE_1X, PERMUTE_0Y>(e1, Vector::Negate(e1));
            p0 = Vector3::Dot(TV0, Axis);
            p1 = Vector3::Dot(TV1, Axis);
            // p2 = XMVector3Dot( V2, Axis ); // p2 = p1;
            Min = Vector::Min(p0, p1);
            Max = Vector::Max(p0, p1);
            Radius = Vector3::Dot(vExtents, Vector::Abs(Axis));
            NoIntersection = Vector::OrInt(NoIntersection, Vector::Greater(Min, Radius));
            NoIntersection = Vector::OrInt(NoIntersection, Vector::Less(Max, Vector::Negate(Radius)));

            // Axis == (0,0,1) x e2 = (e2.z, 0, -e2.x)
            Axis = Vector::Permute<PERMUTE_0Z, PERMUTE_0W, PERMUTE_1X, PERMUTE_0Y>(e2, Vector::Negate(e2));
            p0 = Vector3::Dot(TV0, Axis);
            p1 = Vector3::Dot(TV1, Axis);
            // p2 = XMVector3Dot( V2, Axis ); // p2 = p0;
            Min = Vector::Min(p0, p1);
            Max = Vector::Max(p0, p1);
            Radius = Vector3::Dot(vExtents, Vector::Abs(Axis));
            NoIntersection = Vector::OrInt(NoIntersection, Vector::Greater(Min, Radius));
            NoIntersection = Vector::OrInt(NoIntersection, Vector::Less(Max, Vector::Negate(Radius)));

            // Axis == (0,0,1) x e0 = (-e0.y, e0.x, 0)
            Axis = Vector::Permute<PERMUTE_1Y, PERMUTE_0X, PERMUTE_0W, PERMUTE_0Z>(e0, Vector::Negate(e0));
            p0 = Vector3::Dot(TV0, Axis);
            // p1 = XMVector3Dot( V1, Axis ); // p1 = p0;
            p2 = Vector3::Dot(TV2, Axis);
            Min = Vector::Min(p0, p2);
            Max = Vector::Max(p0, p2);
            Radius = Vector3::Dot(vExtents, Vector::Abs(Axis));
            NoIntersection = Vector::OrInt(NoIntersection, Vector::Greater(Min, Radius));
            NoIntersection = Vector::OrInt(NoIntersection, Vector::Less(Max, Vector::Negate(Radius)));

            // Axis == (0,0,1) x e1 = (-e1.y, e1.x, 0)
            Axis = Vector::Permute<PERMUTE_1Y, PERMUTE_0X, PERMUTE_0W, PERMUTE_0Z>(e1, Vector::Negate(e1));
            p0 = Vector3::Dot(TV0, Axis);
            p1 = Vector3::Dot(TV1, Axis);
            // p2 = XMVector3Dot( V2, Axis ); // p2 = p1;
            Min = Vector::Min(p0, p1);
            Max = Vector::Max(p0, p1);
            Radius = Vector3::Dot(vExtents, Vector::Abs(Axis));
            NoIntersection = Vector::OrInt(NoIntersection, Vector::Greater(Min, Radius));
            NoIntersection = Vector::OrInt(NoIntersection, Vector::Less(Max, Vector::Negate(Radius)));

            // Axis == (0,0,1) x e2 = (-e2.y, e2.x, 0)
            Axis = Vector::Permute<PERMUTE_1Y, PERMUTE_0X, PERMUTE_0W, PERMUTE_0Z>(e2, Vector::Negate(e2));
            p0 = Vector3::Dot(TV0, Axis);
            p1 = Vector3::Dot(TV1, Axis);
            // p2 = XMVector3Dot( V2, Axis ); // p2 = p0;
            Min = Vector::Min(p0, p1);
            Max = Vector::Max(p0, p1);
            Radius = Vector3::Dot(vExtents, Vector::Abs(Axis));
            NoIntersection = Vector::OrInt(NoIntersection, Vector::Greater(Min, Radius));
            NoIntersection = Vector::OrInt(NoIntersection, Vector::Less(Max, Vector::Negate(Radius)));

            return Vector4::NotEqualInt(NoIntersection, Vector::TrueInt());
        }

        _Use_decl_annotations_
        FORCE_INLINE PlaneIntersectionType VEC_CALLCONV BoundingBox::Intersects(A_VECTOR plane) const noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(PlaneIsUnit(plane));
#endif // DEBUG

            // Load the box.
            VECTOR vCenter = Vector::LoadFloat3(&center);
            VECTOR vExtents = Vector::LoadFloat3(&extents);

            // Set w of the center to one so we can dot4 with a plane.
            vCenter = Vector::Insert<0, 0, 0, 0, 1>(vCenter, Vector::SplatOne());

            VECTOR Outside, Inside;
            FastIntersectAxisAlignedBoxPlane(vCenter, vExtents, plane, Outside, Inside);

            // If the box is outside any plane it is outside.
            if (Vector4::EqualInt(Outside, Vector::TrueInt()))
                return FRONT;

            // If the box is inside all planes it is inside.
            if (Vector4::EqualInt(Inside, Vector::TrueInt()))
                return BACK;

            // The box is not inside all planes or outside a plane it intersects.
            return INTERSECTING;
        }
        
        _Use_decl_annotations_
        FORCE_INLINE bool VEC_CALLCONV BoundingBox::Intersects(A_VECTOR origin, A_VECTOR direction, float& dist) const noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(Vector3IsUnit(direction));
#endif // DEBUG

            // Load the box.
            VECTOR vCenter = Vector::LoadFloat3(&center);
            VECTOR vExtents = Vector::LoadFloat3(&extents);

            // Adjust ray origin to be relative to center of the box.
            VECTOR TOrigin = Vector::Subtract(vCenter, origin);

            // Compute the dot product againt each axis of the box.
            // Since the axii are (1,0,0), (0,1,0), (0,0,1) no computation is necessary.
            VECTOR AxisDotOrigin = TOrigin;
            VECTOR AxisDotDirection = direction;

            // if (fabs(AxisDotDirection) <= Epsilon) the ray is nearly parallel to the slab.
            VECTOR IsParallel = Vector::LessOrEqual(Vector::Abs(AxisDotDirection), g_RayEpsilon);

            // Test against all three axii simultaneously.
            VECTOR InverseAxisDotDirection = Vector::Reciprocal(AxisDotDirection);
            VECTOR t1 = Vector::Multiply(Vector::Subtract(AxisDotOrigin, vExtents), InverseAxisDotDirection);
            VECTOR t2 = Vector::Multiply(Vector::Add(AxisDotOrigin, vExtents), InverseAxisDotDirection);

            // Compute the max of min(t1,t2) and the min of max(t1,t2) ensuring we don't
            // use the results from any directions parallel to the slab.
            VECTOR t_min = Vector::Select(Vector::Min(t1, t2), g_FltMin, IsParallel);
            VECTOR t_max = Vector::Select(Vector::Max(t1, t2), g_FltMax, IsParallel);

            // t_min.x = maximum( t_min.x, t_min.y, t_min.z );
            // t_max.x = minimum( t_max.x, t_max.y, t_max.z );
            t_min = Vector::Max(t_min, Vector::SplatY(t_min));  // x = max(x,y)
            t_min = Vector::Max(t_min, Vector::SplatZ(t_min));  // x = max(max(x,y),z)
            t_max = Vector::Min(t_max, Vector::SplatY(t_max));  // x = min(x,y)
            t_max = Vector::Min(t_max, Vector::SplatZ(t_max));  // x = min(min(x,y),z)

            // if ( t_min > t_max ) return false;
            VECTOR NoIntersection = Vector::Greater(Vector::SplatX(t_min), Vector::SplatX(t_max));

            // if ( t_max < 0.0f ) return false;
            NoIntersection = Vector::OrInt(NoIntersection, Vector::Less(Vector::SplatX(t_max), Vector::Zero()));

            // if (IsParallel && (-Extents > AxisDotOrigin || Extents < AxisDotOrigin)) return false;
            VECTOR ParallelOverlap = Vector::InBounds(AxisDotOrigin, vExtents);
            NoIntersection = Vector::OrInt(NoIntersection, Vector::AndCInt(IsParallel, ParallelOverlap));

            if (!Vector3AnyTrue(NoIntersection))
            {
                // Store the x-component to *pDist
                Vector::StoreFloat(&dist, t_min);
                return true;
            }

            dist = 0.f;
            return false;
        }

        _Use_decl_annotations_
        FORCE_INLINE ContainmentType VEC_CALLCONV BoundingBox::ContainedBy(
            A_VECTOR plane0, A_VECTOR plane1, A_VECTOR plane2,
            B_VECTOR plane3,
            C_VECTOR plane4, C_VECTOR plane5) const noexcept
        {
            // Load the box.
            VECTOR vCenter = Vector::LoadFloat3(&center);
            VECTOR vExtents = Vector::LoadFloat3(&extents);

            // Set w of the center to one so we can dot4 with a plane.
            vCenter = Vector::Insert<0, 0, 0, 0, 1>(vCenter, Vector::SplatOne());

            VECTOR Outside, Inside;

            // Test against each plane.
            FastIntersectAxisAlignedBoxPlane(vCenter, vExtents, plane0, Outside, Inside);

            VECTOR AnyOutside = Outside;
            VECTOR AllInside = Inside;

            FastIntersectAxisAlignedBoxPlane(vCenter, vExtents, plane1, Outside, Inside);
            AnyOutside = Vector::OrInt(AnyOutside, Outside);
            AllInside = Vector::AndInt(AllInside, Inside);

            FastIntersectAxisAlignedBoxPlane(vCenter, vExtents, plane2, Outside, Inside);
            AnyOutside = Vector::OrInt(AnyOutside, Outside);
            AllInside = Vector::AndInt(AllInside, Inside);

            FastIntersectAxisAlignedBoxPlane(vCenter, vExtents, plane3, Outside, Inside);
            AnyOutside = Vector::OrInt(AnyOutside, Outside);
            AllInside = Vector::AndInt(AllInside, Inside);

            FastIntersectAxisAlignedBoxPlane(vCenter, vExtents, plane4, Outside, Inside);
            AnyOutside = Vector::OrInt(AnyOutside, Outside);
            AllInside = Vector::AndInt(AllInside, Inside);

            FastIntersectAxisAlignedBoxPlane(vCenter, vExtents, plane5, Outside, Inside);
            AnyOutside = Vector::OrInt(AnyOutside, Outside);
            AllInside = Vector::AndInt(AllInside, Inside);

            // If the box is outside any plane it is outside.
            if (Vector4::EqualInt(AnyOutside, Vector::TrueInt()))
                return DISJOINT;

            // If the box is inside all planes it is inside.
            if (Vector4::EqualInt(AllInside, Vector::TrueInt()))
                return CONTAINS;

            // The box is not inside all planes or outside a plane, it may intersect.
            return INTERSECTS;
        }

        _Use_decl_annotations_
        FORCE_INLINE void BoundingBox::CreateMerged(BoundingBox& out, const BoundingBox& b1, const BoundingBox& b2) noexcept
        {
            VECTOR b1Center = Vector::LoadFloat3(&b1.center);
            VECTOR b1Extents = Vector::LoadFloat3(&b1.extents);

            VECTOR b2Center = Vector::LoadFloat3(&b2.center);
            VECTOR b2Extents = Vector::LoadFloat3(&b2.extents);

            VECTOR Min = Vector::Subtract(b1Center, b1Extents);
            Min = Vector::Min(Min, Vector::Subtract(b2Center, b2Extents));

            VECTOR Max = Vector::Add(b1Center, b1Extents);
            Max = Vector::Max(Max, Vector::Add(b2Center, b2Extents));

#if defined(DEBUG) || defined(_DEBUG)
            assert(Vector3::LessOrEqual(Min, Max));
#endif // DEBUG

            Vector::StoreFloat3(&out.center, Vector::Scale(Vector::Add(Min, Max), 0.5f));
            Vector::StoreFloat3(&out.extents, Vector::Scale(Vector::Subtract(Max, Min), 0.5f));
        }

        _Use_decl_annotations_
        FORCE_INLINE void BoundingBox::CreateFromSphere(BoundingBox& out, const BoundingSphere& sh) noexcept
        {
            VECTOR spCenter = Vector::LoadFloat3(&sh.center);
            VECTOR shRadius = Vector::ReplicatePtr(&sh.radius);

            VECTOR Min = Vector::Subtract(spCenter, shRadius);
            VECTOR Max = Vector::Add(spCenter, shRadius);

#if defined(DEBUG) || defined(_DEBUG)
            assert(Vector3::LessOrEqual(Min, Max));
#endif // DEBUG

            Vector::StoreFloat3(&out.center, Vector::Scale(Vector::Add(Min, Max), 0.5f));
            Vector::StoreFloat3(&out.extents, Vector::Scale(Vector::Subtract(Max, Min), 0.5f));
        }

        _Use_decl_annotations_
        FORCE_INLINE void VEC_CALLCONV BoundingBox::CreateFromPoints(BoundingBox& out, A_VECTOR pt1, A_VECTOR pt2) noexcept
        {
            VECTOR Min = Vector::Min(pt1, pt2);
            VECTOR Max = Vector::Max(pt1, pt2);

            // Store center and extents.
            Vector::StoreFloat3(&out.center, Vector::Scale(Vector::Add(Min, Max), 0.5f));
            Vector::StoreFloat3(&out.extents, Vector::Scale(Vector::Subtract(Max, Min), 0.5f));
        }

        _Use_decl_annotations_
        FORCE_INLINE void BoundingBox::CreateFromPoints(BoundingBox& out, size_t count, const Float3* pPoints, size_t stride) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(count > 0);
            assert(pPoints);
#endif // DEBUG

            // Find the minimum and maximum x, y, and z
            VECTOR vMin, vMax;

            vMin = vMax = Vector::LoadFloat3(pPoints);

            for (size_t i = 1; i < count; ++i)
            {
                VECTOR Point = Vector::LoadFloat3(reinterpret_cast<const Float3*>(reinterpret_cast<const uint8_t*>(pPoints) + i * stride));

                vMin = Vector::Min(vMin, Point);
                vMax = Vector::Max(vMax, Point);
            }

            // Store center and extents.
            Vector::StoreFloat3(&out.center, Vector::Scale(Vector::Add(vMin, vMax), 0.5f));
            Vector::StoreFloat3(&out.extents, Vector::Scale(Vector::Subtract(vMax, vMin), 0.5f));
        }
        
        _Use_decl_annotations_
        FORCE_INLINE void VEC_CALLCONV BoundingOrientedBox::Transform(BoundingOrientedBox& out, A_MATRIX m) const noexcept
        {
            // Load the box.
            VECTOR vCenter = Vector::LoadFloat3(&center);
            VECTOR vExtents = Vector::LoadFloat3(&extents);
            VECTOR vOrientation = Vector::LoadFloat4(&orientation);

#if defined(DEBUG) || defined(_DEBUG)
            assert(QuaternionIsUnit(vOrientation));
#endif // DEBUG

            // Composite the box rotation and the transform rotation.
            MATRIX nM;
            nM.r[0] = Vector3::Normalize(m.r[0]);
            nM.r[1] = Vector3::Normalize(m.r[1]);
            nM.r[2] = Vector3::Normalize(m.r[2]);
            nM.r[3] = g_IdentityR3;
            VECTOR Rotation = Quaternion::RotationMatrix(nM);
            vOrientation = Quaternion::Multiply(vOrientation, Rotation);

            // Transform the center.
            vCenter = Vector3::Transform(vCenter, m);

            // Scale the box extents.
            VECTOR dX = Vector3::Length(m.r[0]);
            VECTOR dY = Vector3::Length(m.r[1]);
            VECTOR dZ = Vector3::Length(m.r[2]);

            VECTOR VectorScale = Vector::Select(dY, dX, g_Select1000);
            VectorScale = Vector::Select(dZ, VectorScale, g_Select1100);
            vExtents = Vector::Multiply(vExtents, VectorScale);

            // Store the box.
            Vector::StoreFloat3(&out.center, vCenter);
            Vector::StoreFloat3(&out.extents, vExtents);
            Vector::StoreFloat4(&out.orientation, vOrientation);
        }

        _Use_decl_annotations_
        FORCE_INLINE void VEC_CALLCONV BoundingOrientedBox::Transform(BoundingOrientedBox& out, float scale, A_VECTOR rotation, A_VECTOR translation) const noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(QuaternionIsUnit(rotation));
#endif // DEBUG

            // Load the box.
            VECTOR vCenter = Vector::LoadFloat3(&center);
            VECTOR vExtents = Vector::LoadFloat3(&extents);
            VECTOR vOrientation = Vector::LoadFloat4(&orientation);

#if defined(DEBUG) || defined(_DEBUG)
            assert(QuaternionIsUnit(vOrientation));
#endif // DEBUG

            // Composite the box rotation and the transform rotation.
            vOrientation = Quaternion::Multiply(vOrientation, rotation);

            // Transform the center.
            VECTOR VectorScale = Vector::Replicate(scale);
            vCenter = Vector::Add(Vector3::Rotate(Vector::Multiply(vCenter, VectorScale), rotation), translation);

            // Scale the box extents.
            vExtents = Vector::Multiply(vExtents, VectorScale);

            // Store the box.
            Vector::StoreFloat3(&out.center, vCenter);
            Vector::StoreFloat3(&out.extents, vExtents);
            Vector::StoreFloat4(&out.orientation, vOrientation);
        }
        
        _Use_decl_annotations_
        FORCE_INLINE void BoundingOrientedBox::GetCorners(Float3* corners) const noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(corners != nullptr);
#endif // DEBUG

            // Load the box
            VECTOR vCenter = Vector::LoadFloat3(&center);
            VECTOR vExtents = Vector::LoadFloat3(&extents);
            VECTOR vOrientation = Vector::LoadFloat4(&orientation);

#if defined(DEBUG) || defined(_DEBUG)
            assert(QuaternionIsUnit(vOrientation));
#endif // DEBUG

            for (size_t i = 0; i < CORNER_COUNT; ++i)
            {
                VECTOR C = Vector::Add(Vector3::Rotate(Vector::Multiply(vExtents, g_BoxOffset[i]), vOrientation), vCenter);
                Vector::StoreFloat3(&corners[i], C);
            }
        }

        _Use_decl_annotations_
        FORCE_INLINE ContainmentType VEC_CALLCONV BoundingOrientedBox::Contains(A_VECTOR point) const noexcept
        {
            VECTOR vCenter = Vector::LoadFloat3(&center);
            VECTOR vExtents = Vector::LoadFloat3(&extents);
            VECTOR vOrientation = Vector::LoadFloat4(&orientation);

            // Transform the point to be local to the box.
            VECTOR TPoint = Vector3::InverseRotate(Vector::Subtract(point, vCenter), vOrientation);

            return Vector3::InBounds(TPoint, vExtents) ? CONTAINS : DISJOINT;
        }

        _Use_decl_annotations_
        FORCE_INLINE ContainmentType VEC_CALLCONV BoundingOrientedBox::Contains(A_VECTOR V0, A_VECTOR V1, A_VECTOR V2) const noexcept
        {
            // Load the box center & orientation.
            VECTOR vCenter = Vector::LoadFloat3(&center);
            VECTOR vOrientation = Vector::LoadFloat4(&orientation);

            // Transform the triangle vertices into the space of the box.
            VECTOR TV0 = Vector3::InverseRotate(Vector::Subtract(V0, vCenter), vOrientation);
            VECTOR TV1 = Vector3::InverseRotate(Vector::Subtract(V1, vCenter), vOrientation);
            VECTOR TV2 = Vector3::InverseRotate(Vector::Subtract(V2, vCenter), vOrientation);

            BoundingBox box;
            box.center = Float3(0.0f, 0.0f, 0.0f);
            box.extents = extents;

            // Use the triangle vs axis aligned box intersection routine.
            return box.Contains(TV0, TV1, TV2);
        }

        _Use_decl_annotations_
        FORCE_INLINE ContainmentType BoundingOrientedBox::Contains(const BoundingSphere& sh) const noexcept
        {
            VECTOR SphereCenter = Vector::LoadFloat3(&sh.center);
            VECTOR SphereRadius = Vector::ReplicatePtr(&sh.radius);

            VECTOR BoxCenter = Vector::LoadFloat3(&center);
            VECTOR BoxExtents = Vector::LoadFloat3(&extents);
            VECTOR BoxOrientation = Vector::LoadFloat4(&orientation);

#if defined(DEBUG) || defined(_DEBUG)
            assert(QuaternionIsUnit(BoxOrientation));
#endif // DEBUG

            // Transform the center of the sphere to be local to the box.
            // BoxMin = -BoxExtents
            // BoxMax = +BoxExtents
            SphereCenter = Vector3::InverseRotate(Vector::Subtract(SphereCenter, BoxCenter), BoxOrientation);

            // Find the distance to the nearest point on the box.
            // for each i in (x, y, z)
            // if (SphereCenter(i) < BoxMin(i)) d2 += (SphereCenter(i) - BoxMin(i)) ^ 2
            // else if (SphereCenter(i) > BoxMax(i)) d2 += (SphereCenter(i) - BoxMax(i)) ^ 2

            VECTOR d = Vector::Zero();

            // Compute d for each dimension.
            VECTOR LessThanMin = Vector::Less(SphereCenter, Vector::Negate(BoxExtents));
            VECTOR GreaterThanMax = Vector::Greater(SphereCenter, BoxExtents);

            VECTOR MinDelta = Vector::Add(SphereCenter, BoxExtents);
            VECTOR MaxDelta = Vector::Subtract(SphereCenter, BoxExtents);

            // Choose value for each dimension based on the comparison.
            d = Vector::Select(d, MinDelta, LessThanMin);
            d = Vector::Select(d, MaxDelta, GreaterThanMax);

            // Use a dot-product to square them and sum them together.
            VECTOR d2 = Vector3::Dot(d, d);
            VECTOR SphereRadiusSq = Vector::Multiply(SphereRadius, SphereRadius);

            if (Vector4::Greater(d2, SphereRadiusSq))
                return DISJOINT;

            // See if we are completely inside the box
            VECTOR SMin = Vector::Subtract(SphereCenter, SphereRadius);
            VECTOR SMax = Vector::Add(SphereCenter, SphereRadius);

            return (Vector3::InBounds(SMin, BoxExtents) && Vector3::InBounds(SMax, BoxExtents)) ? CONTAINS : INTERSECTS;
        }

        _Use_decl_annotations_
        FORCE_INLINE ContainmentType BoundingOrientedBox::Contains(const BoundingBox& box) const noexcept
        {
            // Make the axis aligned box oriented and do an OBB vs OBB test.
            BoundingOrientedBox obox(box.center, box.extents, Float4(0.f, 0.f, 0.f, 1.f));
            return Contains(obox);
        }

        _Use_decl_annotations_
        FORCE_INLINE ContainmentType BoundingOrientedBox::Contains(const BoundingOrientedBox& box) const noexcept
        {
            if (!Intersects(box))
                return DISJOINT;

            // Load the boxes
            VECTOR aCenter = Vector::LoadFloat3(&center);
            VECTOR aExtents = Vector::LoadFloat3(&extents);
            VECTOR aOrientation = Vector::LoadFloat4(&orientation);

#if defined(DEBUG) || defined(_DEBUG)
            assert(QuaternionIsUnit(aOrientation));
#endif // DEBUG

            VECTOR bCenter = Vector::LoadFloat3(&box.center);
            VECTOR bExtents = Vector::LoadFloat3(&box.extents);
            VECTOR bOrientation = Vector::LoadFloat4(&box.orientation);

#if defined(DEBUG) || defined(_DEBUG)
            assert(QuaternionIsUnit(bOrientation));
#endif // DEBUG

            VECTOR offset = Vector::Subtract(bCenter, aCenter);

            for (size_t i = 0; i < CORNER_COUNT; ++i)
            {
                // Cb = rotate( bExtents * corneroffset[i], bOrientation ) + bcenter
                // Ca = invrotate( Cb - aCenter, aOrientation )

                VECTOR C = Vector::Add(Vector3::Rotate(Vector::Multiply(bExtents, g_BoxOffset[i]), bOrientation), offset);
                C = Vector3::InverseRotate(C, aOrientation);

                if (!Vector3::InBounds(C, aExtents))
                    return INTERSECTS;
            }

            return CONTAINS;
        }

        _Use_decl_annotations_
        FORCE_INLINE ContainmentType BoundingOrientedBox::Contains(const BoundingFrustum& fr) const noexcept
        {
            if (!fr.Intersects(*this))
                return DISJOINT;

            Float3 Corners[BoundingFrustum::CORNER_COUNT];
            fr.GetCorners(Corners);

            // Load the box
            VECTOR vCenter = Vector::LoadFloat3(&center);
            VECTOR vExtents = Vector::LoadFloat3(&extents);
            VECTOR vOrientation = Vector::LoadFloat4(&orientation);

#if defined(DEBUG) || defined(_DEBUG)
            assert(QuaternionIsUnit(vOrientation));
#endif // DEBUG

            for (size_t i = 0; i < BoundingFrustum::CORNER_COUNT; ++i)
            {
                VECTOR C = Vector3::InverseRotate(Vector::Subtract(Vector::LoadFloat3(&Corners[i]), vCenter), vOrientation);

                if (!Vector3::InBounds(C, vExtents))
                    return INTERSECTS;
            }

            return CONTAINS;
        }

        _Use_decl_annotations_
        FORCE_INLINE bool BoundingOrientedBox::Intersects(const BoundingSphere& sh) const noexcept
        {
            VECTOR SphereCenter = Vector::LoadFloat3(&sh.center);
            VECTOR SphereRadius = Vector::ReplicatePtr(&sh.radius);

            VECTOR BoxCenter = Vector::LoadFloat3(&center);
            VECTOR BoxExtents = Vector::LoadFloat3(&extents);
            VECTOR BoxOrientation = Vector::LoadFloat4(&orientation);

#if defined(DEBUG) || defined(_DEBUG)
            assert(QuaternionIsUnit(BoxOrientation));
#endif // DEBUG

            // Transform the center of the sphere to be local to the box.
            // BoxMin = -BoxExtents
            // BoxMax = +BoxExtents
            SphereCenter = Vector3::InverseRotate(Vector::Subtract(SphereCenter, BoxCenter), BoxOrientation);

            // Find the distance to the nearest point on the box.
            // for each i in (x, y, z)
            // if (SphereCenter(i) < BoxMin(i)) d2 += (SphereCenter(i) - BoxMin(i)) ^ 2
            // else if (SphereCenter(i) > BoxMax(i)) d2 += (SphereCenter(i) - BoxMax(i)) ^ 2

            VECTOR d = Vector::Zero();

            // Compute d for each dimension.
            VECTOR LessThanMin = Vector::Less(SphereCenter, Vector::Negate(BoxExtents));
            VECTOR GreaterThanMax = Vector::Greater(SphereCenter, BoxExtents);

            VECTOR MinDelta = Vector::Add(SphereCenter, BoxExtents);
            VECTOR MaxDelta = Vector::Subtract(SphereCenter, BoxExtents);

            // Choose value for each dimension based on the comparison.
            d = Vector::Select(d, MinDelta, LessThanMin);
            d = Vector::Select(d, MaxDelta, GreaterThanMax);

            // Use a dot-product to square them and sum them together.
            VECTOR d2 = Vector3::Dot(d, d);

            return Vector4::LessOrEqual(d2, Vector::Multiply(SphereRadius, SphereRadius)) ? true : false;
        }

        _Use_decl_annotations_
        FORCE_INLINE bool BoundingOrientedBox::Intersects(const BoundingBox& box) const noexcept
        {
            // Make the axis aligned box oriented and do an OBB vs OBB test.
            BoundingOrientedBox obox(box.center, box.extents, Float4(0.f, 0.f, 0.f, 1.f));
            return Intersects(obox);
        }

        _Use_decl_annotations_
        FORCE_INLINE bool BoundingOrientedBox::Intersects(const BoundingOrientedBox& box) const noexcept
        {
            // Build the 3x3 rotation matrix that defines the orientation of B relative to A.
            VECTOR A_quat = Vector::LoadFloat4(&orientation);
            VECTOR B_quat = Vector::LoadFloat4(&box.orientation);

#if defined(DEBUG) || defined(_DEBUG)
            assert(QuaternionIsUnit(A_quat));
            assert(QuaternionIsUnit(B_quat));
#endif // DEBUG

            VECTOR Q = Quaternion::Multiply(A_quat, Quaternion::Conjugate(B_quat));
            MATRIX R = Matrix::RotationQuaternion(Q);

            // Compute the translation of B relative to A.
            VECTOR A_cent = Vector::LoadFloat3(&center);
            VECTOR B_cent = Vector::LoadFloat3(&box.center);
            VECTOR t = Vector3::InverseRotate(Vector::Subtract(B_cent, A_cent), A_quat);

            //
            // h(A) = extents of A.
            // h(B) = extents of B.
            //
            // a(u) = axes of A = (1,0,0), (0,1,0), (0,0,1)
            // b(u) = axes of B relative to A = (r00,r10,r20), (r01,r11,r21), (r02,r12,r22)
            //
            // For each possible separating axis l:
            //   d(A) = sum (for i = u,v,w) h(A)(i) * abs( a(i) dot l )
            //   d(B) = sum (for i = u,v,w) h(B)(i) * abs( b(i) dot l )
            //   if abs( t dot l ) > d(A) + d(B) then disjoint
            //

            // Load extents of A and B.
            VECTOR h_A = Vector::LoadFloat3(&extents);
            VECTOR h_B = Vector::LoadFloat3(&box.extents);

            // Rows. Note R[0,1,2]X.w = 0.
            VECTOR R0X = R.r[0];
            VECTOR R1X = R.r[1];
            VECTOR R2X = R.r[2];

            R = Matrix::Transpose(R);

            // Columns. Note RX[0,1,2].w = 0.
            VECTOR RX0 = R.r[0];
            VECTOR RX1 = R.r[1];
            VECTOR RX2 = R.r[2];

            // Absolute value of rows.
            VECTOR AR0X = Vector::Abs(R0X);
            VECTOR AR1X = Vector::Abs(R1X);
            VECTOR AR2X = Vector::Abs(R2X);

            // Absolute value of columns.
            VECTOR ARX0 = Vector::Abs(RX0);
            VECTOR ARX1 = Vector::Abs(RX1);
            VECTOR ARX2 = Vector::Abs(RX2);

            // Test each of the 15 possible seperating axii.
            VECTOR d, d_A, d_B;

            // l = a(u) = (1, 0, 0)
            // t dot l = t.x
            // d(A) = h(A).x
            // d(B) = h(B) dot abs(r00, r01, r02)
            d = Vector::SplatX(t);
            d_A = Vector::SplatX(h_A);
            d_B = Vector3::Dot(h_B, AR0X);
            VECTOR NoIntersection = Vector::Greater(Vector::Abs(d), Vector::Add(d_A, d_B));

            // l = a(v) = (0, 1, 0)
            // t dot l = t.y
            // d(A) = h(A).y
            // d(B) = h(B) dot abs(r10, r11, r12)
            d = Vector::SplatY(t);
            d_A = Vector::SplatY(h_A);
            d_B = Vector3::Dot(h_B, AR1X);
            NoIntersection = Vector::OrInt(NoIntersection,
                Vector::Greater(Vector::Abs(d), Vector::Add(d_A, d_B)));

            // l = a(w) = (0, 0, 1)
            // t dot l = t.z
            // d(A) = h(A).z
            // d(B) = h(B) dot abs(r20, r21, r22)
            d = Vector::SplatZ(t);
            d_A = Vector::SplatZ(h_A);
            d_B = Vector3::Dot(h_B, AR2X);
            NoIntersection = Vector::OrInt(NoIntersection,
                Vector::Greater(Vector::Abs(d), Vector::Add(d_A, d_B)));

            // l = b(u) = (r00, r10, r20)
            // d(A) = h(A) dot abs(r00, r10, r20)
            // d(B) = h(B).x
            d = Vector3::Dot(t, RX0);
            d_A = Vector3::Dot(h_A, ARX0);
            d_B = Vector::SplatX(h_B);
            NoIntersection = Vector::OrInt(NoIntersection,
                Vector::Greater(Vector::Abs(d), Vector::Add(d_A, d_B)));

            // l = b(v) = (r01, r11, r21)
            // d(A) = h(A) dot abs(r01, r11, r21)
            // d(B) = h(B).y
            d = Vector3::Dot(t, RX1);
            d_A = Vector3::Dot(h_A, ARX1);
            d_B = Vector::SplatY(h_B);
            NoIntersection = Vector::OrInt(NoIntersection,
                Vector::Greater(Vector::Abs(d), Vector::Add(d_A, d_B)));

            // l = b(w) = (r02, r12, r22)
            // d(A) = h(A) dot abs(r02, r12, r22)
            // d(B) = h(B).z
            d = Vector3::Dot(t, RX2);
            d_A = Vector3::Dot(h_A, ARX2);
            d_B = Vector::SplatZ(h_B);
            NoIntersection = Vector::OrInt(NoIntersection,
                Vector::Greater(Vector::Abs(d), Vector::Add(d_A, d_B)));

            // l = a(u) x b(u) = (0, -r20, r10)
            // d(A) = h(A) dot abs(0, r20, r10)
            // d(B) = h(B) dot abs(0, r02, r01)
            d = Vector3::Dot(t, Vector::Permute<PERMUTE_0W, PERMUTE_1Z, PERMUTE_0Y, PERMUTE_0X>(RX0, Vector::Negate(RX0)));
            d_A = Vector3::Dot(h_A, Vector::Swizzle<SWIZZLE_W, SWIZZLE_Z, SWIZZLE_Y, SWIZZLE_X>(ARX0));
            d_B = Vector3::Dot(h_B, Vector::Swizzle<SWIZZLE_W, SWIZZLE_Z, SWIZZLE_Y, SWIZZLE_X>(AR0X));
            NoIntersection = Vector::OrInt(NoIntersection,
                Vector::Greater(Vector::Abs(d), Vector::Add(d_A, d_B)));

            // l = a(u) x b(v) = (0, -r21, r11)
            // d(A) = h(A) dot abs(0, r21, r11)
            // d(B) = h(B) dot abs(r02, 0, r00)
            d = Vector3::Dot(t, Vector::Permute<PERMUTE_0W, PERMUTE_1Z, PERMUTE_0Y, PERMUTE_0X>(RX1, Vector::Negate(RX1)));
            d_A = Vector3::Dot(h_A, Vector::Swizzle<SWIZZLE_W, SWIZZLE_Z, SWIZZLE_Y, SWIZZLE_X>(ARX1));
            d_B = Vector3::Dot(h_B, Vector::Swizzle<SWIZZLE_Z, SWIZZLE_W, SWIZZLE_X, SWIZZLE_Y>(AR0X));
            NoIntersection = Vector::OrInt(NoIntersection,
                Vector::Greater(Vector::Abs(d), Vector::Add(d_A, d_B)));

            // l = a(u) x b(w) = (0, -r22, r12)
            // d(A) = h(A) dot abs(0, r22, r12)
            // d(B) = h(B) dot abs(r01, r00, 0)
            d = Vector3::Dot(t, Vector::Permute<PERMUTE_0W, PERMUTE_1Z, PERMUTE_0Y, PERMUTE_0X>(RX2, Vector::Negate(RX2)));
            d_A = Vector3::Dot(h_A, Vector::Swizzle<SWIZZLE_W, SWIZZLE_Z, SWIZZLE_Y, SWIZZLE_X>(ARX2));
            d_B = Vector3::Dot(h_B, Vector::Swizzle<SWIZZLE_Y, SWIZZLE_X, SWIZZLE_W, SWIZZLE_Z>(AR0X));
            NoIntersection = Vector::OrInt(NoIntersection,
                Vector::Greater(Vector::Abs(d), Vector::Add(d_A, d_B)));

            // l = a(v) x b(u) = (r20, 0, -r00)
            // d(A) = h(A) dot abs(r20, 0, r00)
            // d(B) = h(B) dot abs(0, r12, r11)
            d = Vector3::Dot(t, Vector::Permute<PERMUTE_0Z, PERMUTE_0W, PERMUTE_1X, PERMUTE_0Y>(RX0, Vector::Negate(RX0)));
            d_A = Vector3::Dot(h_A, Vector::Swizzle<SWIZZLE_Z, SWIZZLE_W, SWIZZLE_X, SWIZZLE_Y>(ARX0));
            d_B = Vector3::Dot(h_B, Vector::Swizzle<SWIZZLE_W, SWIZZLE_Z, SWIZZLE_Y, SWIZZLE_X>(AR1X));
            NoIntersection = Vector::OrInt(NoIntersection,
                Vector::Greater(Vector::Abs(d), Vector::Add(d_A, d_B)));

            // l = a(v) x b(v) = (r21, 0, -r01)
            // d(A) = h(A) dot abs(r21, 0, r01)
            // d(B) = h(B) dot abs(r12, 0, r10)
            d = Vector3::Dot(t, Vector::Permute<PERMUTE_0Z, PERMUTE_0W, PERMUTE_1X, PERMUTE_0Y>(RX1, Vector::Negate(RX1)));
            d_A = Vector3::Dot(h_A, Vector::Swizzle<SWIZZLE_Z, SWIZZLE_W, SWIZZLE_X, SWIZZLE_Y>(ARX1));
            d_B = Vector3::Dot(h_B, Vector::Swizzle<SWIZZLE_Z, SWIZZLE_W, SWIZZLE_X, SWIZZLE_Y>(AR1X));
            NoIntersection = Vector::OrInt(NoIntersection,
                Vector::Greater(Vector::Abs(d), Vector::Add(d_A, d_B)));

            // l = a(v) x b(w) = (r22, 0, -r02)
            // d(A) = h(A) dot abs(r22, 0, r02)
            // d(B) = h(B) dot abs(r11, r10, 0)
            d = Vector3::Dot(t, Vector::Permute<PERMUTE_0Z, PERMUTE_0W, PERMUTE_1X, PERMUTE_0Y>(RX2, Vector::Negate(RX2)));
            d_A = Vector3::Dot(h_A, Vector::Swizzle<SWIZZLE_Z, SWIZZLE_W, SWIZZLE_X, SWIZZLE_Y>(ARX2));
            d_B = Vector3::Dot(h_B, Vector::Swizzle<SWIZZLE_Y, SWIZZLE_X, SWIZZLE_W, SWIZZLE_Z>(AR1X));
            NoIntersection = Vector::OrInt(NoIntersection,
                Vector::Greater(Vector::Abs(d), Vector::Add(d_A, d_B)));

            // l = a(w) x b(u) = (-r10, r00, 0)
            // d(A) = h(A) dot abs(r10, r00, 0)
            // d(B) = h(B) dot abs(0, r22, r21)
            d = Vector3::Dot(t, Vector::Permute<PERMUTE_1Y, PERMUTE_0X, PERMUTE_0W, PERMUTE_0Z>(RX0, Vector::Negate(RX0)));
            d_A = Vector3::Dot(h_A, Vector::Swizzle<SWIZZLE_Y, SWIZZLE_X, SWIZZLE_W, SWIZZLE_Z>(ARX0));
            d_B = Vector3::Dot(h_B, Vector::Swizzle<SWIZZLE_W, SWIZZLE_Z, SWIZZLE_Y, SWIZZLE_X>(AR2X));
            NoIntersection = Vector::OrInt(NoIntersection,
                    Vector::Greater(Vector::Abs(d), Vector::Add(d_A, d_B)));

            // l = a(w) x b(v) = (-r11, r01, 0)
            // d(A) = h(A) dot abs(r11, r01, 0)
            // d(B) = h(B) dot abs(r22, 0, r20)
            d = Vector3::Dot(t, Vector::Permute<PERMUTE_1Y, PERMUTE_0X, PERMUTE_0W, PERMUTE_0Z>(RX1, Vector::Negate(RX1)));
            d_A = Vector3::Dot(h_A, Vector::Swizzle<SWIZZLE_Y, SWIZZLE_X, SWIZZLE_W, SWIZZLE_Z>(ARX1));
            d_B = Vector3::Dot(h_B, Vector::Swizzle<SWIZZLE_Z, SWIZZLE_W, SWIZZLE_X, SWIZZLE_Y>(AR2X));
            NoIntersection = Vector::OrInt(NoIntersection,
                Vector::Greater(Vector::Abs(d), Vector::Add(d_A, d_B)));

            // l = a(w) x b(w) = (-r12, r02, 0)
            // d(A) = h(A) dot abs(r12, r02, 0)
            // d(B) = h(B) dot abs(r21, r20, 0)
            d = Vector3::Dot(t, Vector::Permute<PERMUTE_1Y, PERMUTE_0X, PERMUTE_0W, PERMUTE_0Z>(RX2, Vector::Negate(RX2)));
            d_A = Vector3::Dot(h_A, Vector::Swizzle<SWIZZLE_Y, SWIZZLE_X, SWIZZLE_W, SWIZZLE_Z>(ARX2));
            d_B = Vector3::Dot(h_B, Vector::Swizzle<SWIZZLE_Y, SWIZZLE_X, SWIZZLE_W, SWIZZLE_Z>(AR2X));
            NoIntersection = Vector::OrInt(NoIntersection,
                Vector::Greater(Vector::Abs(d), Vector::Add(d_A, d_B)));

            // No seperating axis found, boxes must intersect.
            return Vector4::NotEqualInt(NoIntersection, Vector::TrueInt()) ? true : false;
        }

        _Use_decl_annotations_
        FORCE_INLINE bool BoundingOrientedBox::Intersects(const BoundingFrustum& fr) const noexcept
        {
            return fr.Intersects(*this);
        }

        _Use_decl_annotations_
        FORCE_INLINE bool VEC_CALLCONV BoundingOrientedBox::Intersects(A_VECTOR V0, A_VECTOR V1, A_VECTOR V2) const noexcept
        {
            // Load the box center & orientation.
            VECTOR vCenter = Vector::LoadFloat3(&center);
            VECTOR vOrientation = Vector::LoadFloat4(&orientation);

            // Transform the triangle vertices into the space of the box.
            VECTOR TV0 = Vector3::InverseRotate(Vector::Subtract(V0, vCenter), vOrientation);
            VECTOR TV1 = Vector3::InverseRotate(Vector::Subtract(V1, vCenter), vOrientation);
            VECTOR TV2 = Vector3::InverseRotate(Vector::Subtract(V2, vCenter), vOrientation);

            BoundingBox box;
            box.center = Float3(0.0f, 0.0f, 0.0f);
            box.extents = extents;

            // Use the triangle vs axis aligned box intersection routine.
            return box.Intersects(TV0, TV1, TV2);
        }

        _Use_decl_annotations_
        FORCE_INLINE PlaneIntersectionType VEC_CALLCONV BoundingOrientedBox::Intersects(A_VECTOR plane) const noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(PlaneIsUnit(plane));
#endif // DEBUG

            // Load the box.
            VECTOR vCenter = Vector::LoadFloat3(&center);
            VECTOR vExtents = Vector::LoadFloat3(&extents);
            VECTOR BoxOrientation = Vector::LoadFloat4(&orientation);

#if defined(DEBUG) || defined(_DEBUG)
            assert(QuaternionIsUnit(BoxOrientation));
#endif // DEBUG

            // Set w of the center to one so we can dot4 with a plane.
            vCenter = Vector::Insert<0, 0, 0, 0, 1>(vCenter, Vector::SplatOne());

            // Build the 3x3 rotation matrix that defines the box axes.
            MATRIX R = Matrix::RotationQuaternion(BoxOrientation);

            VECTOR Outside, Inside;
            FastIntersectOrientedBoxPlane(vCenter, vExtents, R.r[0], R.r[1], R.r[2], plane, Outside, Inside);

            // If the box is outside any plane it is outside.
            if (Vector4::EqualInt(Outside, Vector::TrueInt()))
                return FRONT;

            // If the box is inside all planes it is inside.
            if (Vector4::EqualInt(Inside, Vector::TrueInt()))
                return BACK;

            // The box is not inside all planes or outside a plane it intersects.
            return INTERSECTING;
        }

        _Use_decl_annotations_
        FORCE_INLINE bool VEC_CALLCONV BoundingOrientedBox::Intersects(A_VECTOR origin, A_VECTOR direction, float& dist) const noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(Vector3IsUnit(direction));
#endif // DEBUG

            static const VECTOR_U32 SelectY = { { { SELECT_0, SELECT_1, SELECT_0, SELECT_0 } } };
            static const VECTOR_U32 SelectZ = { { { SELECT_0, SELECT_0, SELECT_1, SELECT_0 } } };

            // Load the box.
            VECTOR vCenter = Vector::LoadFloat3(&center);
            VECTOR vExtents = Vector::LoadFloat3(&extents);
            VECTOR vOrientation = Vector::LoadFloat4(&orientation);

#if defined(DEBUG) || defined(_DEBUG)
            assert(QuaternionIsUnit(vOrientation));
#endif // DEBUG

            // Get the boxes normalized side directions.
            MATRIX R = Matrix::RotationQuaternion(vOrientation);

            // Adjust ray origin to be relative to center of the box.
            VECTOR TOrigin = Vector::Subtract(vCenter, origin);

            // Compute the dot product againt each axis of the box.
            VECTOR AxisDotOrigin = Vector3::Dot(R.r[0], TOrigin);
            AxisDotOrigin = Vector::Select(AxisDotOrigin, Vector3::Dot(R.r[1], TOrigin), SelectY);
            AxisDotOrigin = Vector::Select(AxisDotOrigin, Vector3::Dot(R.r[2], TOrigin), SelectZ);

            VECTOR AxisDotDirection = VE3::Dot(R.r[0], direction);
            AxisDotDirection = Vector::Select(AxisDotDirection, Vector3::Dot(R.r[1], direction), SelectY);
            AxisDotDirection = Vector::Select(AxisDotDirection, Vector3::Dot(R.r[2], direction), SelectZ);

            // if (fabs(AxisDotDirection) <= Epsilon) the ray is nearly parallel to the slab.
            VECTOR IsParallel = Vector::LessOrEqual(Vector::Abs(AxisDotDirection), g_RayEpsilon);

            // Test against all three axes simultaneously.
            VECTOR InverseAxisDotDirection = Vector::Reciprocal(AxisDotDirection);
            VECTOR t1 = Vector::Multiply(Vector::Subtract(AxisDotOrigin, vExtents), InverseAxisDotDirection);
            VECTOR t2 = Vector::Multiply(Vector::Add(AxisDotOrigin, vExtents), InverseAxisDotDirection);

            // Compute the max of min(t1,t2) and the min of max(t1,t2) ensuring we don't
            // use the results from any directions parallel to the slab.
            VECTOR t_min = Vector::Select(Vector::Min(t1, t2), g_FltMin, IsParallel);
            VECTOR t_max = Vector::Select(Vector::Max(t1, t2), g_FltMax, IsParallel);

            // t_min.x = maximum( t_min.x, t_min.y, t_min.z );
            // t_max.x = minimum( t_max.x, t_max.y, t_max.z );
            t_min = Vector::Max(t_min, Vector::SplatY(t_min));  // x = max(x,y)
            t_min = Vector::Max(t_min, Vector::SplatZ(t_min));  // x = max(max(x,y),z)
            t_max = Vector::Min(t_max, Vector::SplatY(t_max));  // x = min(x,y)
            t_max = Vector::Min(t_max, Vector::SplatZ(t_max));  // x = min(min(x,y),z)

            // if ( t_min > t_max ) return false;
            VECTOR NoIntersection = Vector::Greater(Vector::SplatX(t_min), Vector::SplatX(t_max));

            // if ( t_max < 0.0f ) return false;
            NoIntersection = Vector::OrInt(NoIntersection, Vector::Less(Vector::SplatX(t_max), Vector::Zero()));

            // if (IsParallel && (-Extents > AxisDotOrigin || Extents < AxisDotOrigin)) return false;
            VECTOR ParallelOverlap = Vector::InBounds(AxisDotOrigin, vExtents);
            NoIntersection = Vector::OrInt(NoIntersection, Vector::AndCInt(IsParallel, ParallelOverlap));

            if (!Vector3AnyTrue(NoIntersection))
            {
                // Store the x-component to *pDist
                StoreFloat(&dist, t_min);
                return true;
            }

            dist = 0.f;
            return false;
        }

        _Use_decl_annotations_
        FORCE_INLINE ContainmentType VEC_CALLCONV BoundingOrientedBox::ContainedBy(
            A_VECTOR plane0, A_VECTOR plane1, A_VECTOR plane2,
            B_VECTOR plane3,
            C_VECTOR plane4, C_VECTOR plane5) const noexcept
        {
            // Load the box.
            VECTOR vCenter = Vector::LoadFloat3(&center);
            VECTOR vExtents = Vector::LoadFloat3(&extents);
            VECTOR BoxOrientation = Vector::LoadFloat4(&orientation);

#if defined(DEBUG) || defined(_DEBUG)
            assert(QuaternionIsUnit(BoxOrientation));
#endif // DEBUG

            // Set w of the center to one so we can dot4 with a plane.
            vCenter = Vector::Insert<0, 0, 0, 0, 1>(vCenter, Vector::SplatOne());

            // Build the 3x3 rotation matrix that defines the box axes.
            MATRIX R = Matrix::RotationQuaternion(BoxOrientation);

            VECTOR Outside, Inside;

            // Test against each plane.
            FastIntersectOrientedBoxPlane(vCenter, vExtents, R.r[0], R.r[1], R.r[2], plane0, Outside, Inside);

            VECTOR AnyOutside = Outside;
            VECTOR AllInside = Inside;

            FastIntersectOrientedBoxPlane(vCenter, vExtents, R.r[0], R.r[1], R.r[2], plane1, Outside, Inside);
            AnyOutside = Vector::OrInt(AnyOutside, Outside);
            AllInside = Vector::AndInt(AllInside, Inside);

            FastIntersectOrientedBoxPlane(vCenter, vExtents, R.r[0], R.r[1], R.r[2], plane2, Outside, Inside);
            AnyOutside = Vector::OrInt(AnyOutside, Outside);
            AllInside = Vector::AndInt(AllInside, Inside);

            FastIntersectOrientedBoxPlane(vCenter, vExtents, R.r[0], R.r[1], R.r[2], plane3, Outside, Inside);
            AnyOutside = Vector::OrInt(AnyOutside, Outside);
            AllInside = Vector::AndInt(AllInside, Inside);

            FastIntersectOrientedBoxPlane(vCenter, vExtents, R.r[0], R.r[1], R.r[2], plane4, Outside, Inside);
            AnyOutside = Vector::OrInt(AnyOutside, Outside);
            AllInside = Vector::AndInt(AllInside, Inside);

            FastIntersectOrientedBoxPlane(vCenter, vExtents, R.r[0], R.r[1], R.r[2], plane5, Outside, Inside);
            AnyOutside = Vector::OrInt(AnyOutside, Outside);
            AllInside = Vector::AndInt(AllInside, Inside);

            // If the box is outside any plane it is outside.
            if (Vector4::EqualInt(AnyOutside, Vector::TrueInt()))
                return DISJOINT;

            // If the box is inside all planes it is inside.
            if (Vector4::EqualInt(AllInside, Vector::TrueInt()))
                return CONTAINS;

            // The box is not inside all planes or outside a plane, it may intersect.
            return INTERSECTS;
        }

        _Use_decl_annotations_
        FORCE_INLINE void BoundingOrientedBox::CreateFromBoundingBox(BoundingOrientedBox& out, const BoundingBox& box) noexcept
        {
            out.center = box.center;
            out.extents = box.extents;
            out.orientation = Float4(0.f, 0.f, 0.f, 1.f);
        }

        _Use_decl_annotations_
        FORCE_INLINE void BoundingOrientedBox::CreateFromPoints(BoundingOrientedBox& out, size_t count, const Float3* pPoints, size_t stride) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(count > 0);
            assert(pPoints != nullptr);
#endif // DEBUG

            VECTOR CenterOfMass = Vector::Zero();

            // Compute the center of mass and inertia tensor of the points.
            for (size_t i = 0; i < count; ++i)
            {
                VECTOR Point = Vector::LoadFloat3(reinterpret_cast<const Float3*>(reinterpret_cast<const uint8_t*>(pPoints) + i * stride));

                CenterOfMass = Vector::Add(CenterOfMass, Point);
            }

            CenterOfMass = Vector::Multiply(CenterOfMass, Vector::Reciprocal(Vector::Replicate(float(count))));

            // Compute the inertia tensor of the points around the center of mass.
            // Using the center of mass is not strictly necessary, but will hopefully
            // improve the stability of finding the eigenvectors.
            VECTOR XX_YY_ZZ = Vector::Zero();
            VECTOR XY_XZ_YZ = Vector::Zero();

            for (size_t i = 0; i < count; ++i)
            {
                VECTOR Point = Vector::Subtract(Vector::LoadFloat3(reinterpret_cast<const Float3*>(reinterpret_cast<const uint8_t*>(pPoints) + i * stride)), CenterOfMass);

                XX_YY_ZZ = Vector::Add(XX_YY_ZZ, Vector::Multiply(Point, Point));

                VECTOR XXY = Vector::Swizzle<SWIZZLE_X, SWIZZLE_X, SWIZZLE_Y, SWIZZLE_W>(Point);
                VECTOR YZZ = Vector::Swizzle<SWIZZLE_Y, SWIZZLE_Z, SWIZZLE_Z, SWIZZLE_W>(Point);

                XY_XZ_YZ = Vector::Add(XY_XZ_YZ, Vector::Multiply(XXY, YZZ));
            }

            VECTOR v1, v2, v3;

            // Compute the eigenvectors of the inertia tensor.
            CalculateEigenVectorsFromCovarianceMatrix(Vector::GetX(XX_YY_ZZ), Vector::GetY(XX_YY_ZZ),
                Vector::GetZ(XX_YY_ZZ),
                Vector::GetX(XY_XZ_YZ), Vector::GetY(XY_XZ_YZ),
                Vector::GetZ(XY_XZ_YZ),
                &v1, &v2, &v3);

            // Put them in a matrix.
            MATRIX R;

            R.r[0] = Vector::SetW(v1, 0.f);
            R.r[1] = Vector::SetW(v2, 0.f);
            R.r[2] = Vector::SetW(v3, 0.f);
            R.r[3] = g_IdentityR3.v;

            // Multiply by -1 to convert the matrix into a right handed coordinate
            // system (Det ~= 1) in case the eigenvectors form a left handed
            // coordinate system (Det ~= -1) because XMQuaternionRotationMatrix only
            // works on right handed matrices.
            VECTOR Det = Matrix::Determinant(R);

            if (Vector4::Less(Det, Vector::Zero()))
            {
                R.r[0] = Vector::Multiply(R.r[0], g_NegativeOne.v);
                R.r[1] = Vector::Multiply(R.r[1], g_NegativeOne.v);
                R.r[2] = Vector::Multiply(R.r[2], g_NegativeOne.v);
            }

            // Get the rotation quaternion from the matrix.
            VECTOR vOrientation = Quaternion::RotationMatrix(R);

            // Make sure it is normal (in case the vectors are slightly non-orthogonal).
            vOrientation = Quaternion::Normalize(vOrientation);

            // Rebuild the rotation matrix from the quaternion.
            R = Matrix::RotationQuaternion(vOrientation);

            // Build the rotation into the rotated space.
            MATRIX InverseR = Matrix::Transpose(R);

            // Find the minimum OBB using the eigenvectors as the axes.
            VECTOR vMin, vMax;

            vMin = vMax = Vector3::TransformNormal(Vector::LoadFloat3(pPoints), InverseR);

            for (size_t i = 1; i < count; ++i)
            {
                VECTOR Point = Vector3::TransformNormal(Vector::LoadFloat3(reinterpret_cast<const Float3*>(reinterpret_cast<const uint8_t*>(pPoints) + i * stride)),
                    InverseR);

                vMin = Vector::Min(vMin, Point);
                vMax = Vector::Max(vMax, Point);
            }

            // Rotate the center into world space.
            VECTOR vCenter = Vector::Scale(Vector::Add(vMin, vMax), 0.5f);
            vCenter = Vector3::TransformNormal(vCenter, R);

            // Store center, extents, and orientation.
            Vector::StoreFloat3(&out.center, vCenter);
            Vector::StoreFloat3(&out.extents, Vector::Scale(Vector::Subtract(vMax, vMin), 0.5f));
            Vector::StoreFloat4(&out.orientation, vOrientation);
        }

        _Use_decl_annotations_
        FORCE_INLINE BoundingFrustum::BoundingFrustum(B_MATRIX projection, bool rhcoords) noexcept
        {
            CreateFromMatrix(*this, projection, rhcoords);
        }

        _Use_decl_annotations_
        FORCE_INLINE void VEC_CALLCONV BoundingFrustum::Transform(BoundingFrustum& out, A_MATRIX m) const noexcept
        {
            // Load the frustum.
            VECTOR vOrigin = Vector::LoadFloat3(&origin);
            VECTOR vOrientation = Vector::LoadFloat4(&orientation);

#if defined(DEBUG) || defined(_DEBUG)
            assert(QuaternionIsUnit(vOrientation));
#endif // DEBUG

            // Composite the frustum rotation and the transform rotation
            MATRIX nM;
            nM.r[0] = Vector3::Normalize(m.r[0]);
            nM.r[1] = Vector3::Normalize(m.r[1]);
            nM.r[2] = Vector3::Normalize(m.r[2]);
            nM.r[3] = g_IdentityR3;
            VECTOR Rotation = Quaternion::RotationMatrix(nM);
            vOrientation = Quaternion::Multiply(vOrientation, Rotation);

            // Transform the center.
            vOrigin = Vector3::Transform(vOrigin, m);

            // Store the frustum.
            Vector::StoreFloat3(&out.origin, vOrigin);
            Vector::StoreFloat4(&out.orientation, vOrientation);

            // Scale the near and far distances (the slopes remain the same).
            VECTOR dX = Vector3::Dot(m.r[0], m.r[0]);
            VECTOR dY = Vector3::Dot(m.r[1], m.r[1]);
            VECTOR dZ = Vector3::Dot(m.r[2], m.r[2]);

            VECTOR d = Vector::Max(dX, Vector::Max(dY, dZ));
            float Scale = sqrtf(Vector::GetX(d));

            out.near = near * scale;
            out.far = far * scale;

            // Copy the slopes.
            out.rightSlope = rightSlope;
            out.leftSlope = leftSlope;
            out.topSlope = topSlope;
            out.bottomSlope = bottomSlope;
        }

        _Use_decl_annotations_
        FORCE_INLINE void VEC_CALLCONV BoundingFrustum::Transform(BoundingFrustum& out, float scale, A_VECTOR rotation, A_VECTOR translation) const noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(QuaternionIsUnit(rotation));
#endif // DEBUG

            // Load the frustum.
            VECTOR vOrigin = Vector::LoadFloat3(&origin);
            VECTOR vOrientation = Vector::LoadFloat4(&orientation);

#if defined(DEBUG) || defined(_DEBUG)
            assert(QuaternionIsUnit(vOrientation));
#endif // DEBUG

            // Composite the frustum rotation and the transform rotation.
            vOrientation = Quaternion::Multiply(vOrientation, rotation);

            // Transform the origin.
            vOrigin = Vector::Add(Vector3::Rotate(Vector::Scale(vOrigin, scale), rotation), translation);

            // Store the frustum.
            Vector::StoreFloat3(&out.origin, vOrigin);
            Vector::StoreFloat4(&out.orientation, vOrientation);

            // Scale the near and far distances (the slopes remain the same).
            out.near = near * scale;
            out.far = far * scale;

            // Copy the slopes.
            out.rightSlope = rightSlope;
            out.leftSlope = leftSlope;
            out.topSlope = topSlope;
            out.bottomSlope = bottomSlope;
        }

        _Use_decl_annotations_
        FORCE_INLINE void BoundingFrustum::GetCorners(Float3* corners) const noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(corners != nullptr);
#endif // DEBUG

            // Load origin and orientation of the frustum.
            VECTOR vOrigin = Vector::LoadFloat3(&origin);
            VECTOR vOrientation = Vector::LoadFloat4(&orientation);

#if defined(DEBUG) || defined(_DEBUG)
            assert(QuaternionIsUnit(vOrientation));
#endif // DEBUG

            // Build the corners of the frustum.
            VECTOR vRightTop = Vector::Set(rightSlope, topSlope, 1.0f, 0.0f);
            VECTOR vRightBottom = Vector::Set(rightSlope, bottomSlope, 1.0f, 0.0f);
            VECTOR vLeftTop = Vector::Set(leftSlope, topSlope, 1.0f, 0.0f);
            VECTOR vLeftBottom = Vector::Set(leftSlope, bottomSlope, 1.0f, 0.0f);
            VECTOR vNear = Vector::ReplicatePtr(&near);
            VECTOR vFar = Vector::ReplicatePtr(&far);

            // Returns 8 corners position of bounding frustum.
            //     Near    Far
            //    0----1  4----5
            //    |    |  |    |
            //    |    |  |    |
            //    3----2  7----6

            VECTOR vCorners[CORNER_COUNT];
            vCorners[0] = Vector::Multiply(vLeftTop, vNear);
            vCorners[1] = Vector::Multiply(vRightTop, vNear);
            vCorners[2] = Vector::Multiply(vRightBottom, vNear);
            vCorners[3] = Vector::Multiply(vLeftBottom, vNear);
            vCorners[4] = Vector::Multiply(vLeftTop, vFar);
            vCorners[5] = Vector::Multiply(vRightTop, vFar);
            vCorners[6] = Vector::Multiply(vRightBottom, vFar);
            vCorners[7] = Vector::Multiply(vLeftBottom, vFar);

            for (size_t i = 0; i < CORNER_COUNT; ++i)
            {
                VECTOR C = Vector::Add(Vector3::Rotate(vCorners[i], vOrientation), vOrigin);
                Vector::StoreFloat3(&corners[i], C);
            }
        }

        _Use_decl_annotations_
        FORCE_INLINE ContainmentType VEC_CALLCONV BoundingFrustum::Contains(A_VECTOR point) const noexcept
        {
            // Build frustum planes.
            VECTOR Planes[6];
            Planes[0] = Vector::Set(0.0f, 0.0f, -1.0f, near);
            Planes[1] = Vector::Set(0.0f, 0.0f, 1.0f, -far);
            Planes[2] = Vector::Set(1.0f, 0.0f, -rightSlope, 0.0f);
            Planes[3] = Vector::Set(-1.0f, 0.0f, leftSlope, 0.0f);
            Planes[4] = Vector::Set(0.0f, 1.0f, -topSlope, 0.0f);
            Planes[5] = Vector::Set(0.0f, -1.0f, bottomSlope, 0.0f);

            // Load origin and orientation.
            VECTOR vOrigin = Vector::LoadFloat3(&origin);
            VECTOR vOrientation = Vector::LoadFloat4(&orientation);

#if defined(DEBUG) || defined(_DEBUG)
            assert(QuaternionIsUnit(vOrientation));
#endif // DEBUG

            // Transform point into local space of frustum.
            VECTOR TPoint = Vector3::InverseRotate(Vector::Subtract(point, vOrigin), vOrientation);

            // Set w to one.
            TPoint = Vector::Insert<0, 0, 0, 0, 1>(TPoint, Vector::SplatOne());

            VECTOR zero = Vector::Zero();
            VECTOR outside = zero;

            // Test point against each plane of the frustum.
            for (size_t i = 0; i < 6; ++i)
            {
                VECTOR dot = Vector4::Dot(TPoint, Planes[i]);
                outside = Vector::OrInt(Outside, Vector::Greater(dot, zero));
            }

            return Vector4::NotEqualInt(outside, Vector::TrueInt()) ? CONTAINS : DISJOINT;
        }

        _Use_decl_annotations_
        FORCE_INLINE ContainmentType VEC_CALLCONV BoundingFrustum::Contains(A_VECTOR V0, A_VECTOR V1, A_VECTOR V2) const noexcept
        {
            // Load origin and orientation of the frustum.
            VECTOR vOrigin = Vector::LoadFloat3(&origin);
            VECTOR vOrientation = Vector::LoadFloat4(&orientation);

            // Create 6 planes (do it inline to encourage use of registers)
            VECTOR nearPlane = Vector::Set(0.0f, 0.0f, -1.0f, near);
            nearPlane = PlaneTransform(nearPlane, vOrientation, vOrigin);
            nearPlane = PlaneNormalize(nearPlane);

            VECTOR farPlane = Vector::Set(0.0f, 0.0f, 1.0f, -far);
            farPlane = PlaneTransform(farPlane, vOrientation, vOrigin);
            farPlane = PlaneNormalize(farPlane);

            VECTOR rightPlane = Vector::Set(1.0f, 0.0f, -rightSlope, 0.0f);
            rightPlane = PlaneTransform(rightPlane, vOrientation, vOrigin);
            rightPlane = PlaneNormalize(rightPlane);

            VECTOR leftPlane = Vector::Set(-1.0f, 0.0f, leftSlope, 0.0f);
            leftPlane = PlaneTransform(leftPlane, vOrientation, vOrigin);
            leftPlane = PlaneNormalize(leftPlane);

            VECTOR topPlane = Vector::Set(0.0f, 1.0f, -topSlope, 0.0f);
            topPlane = PlaneTransform(topPlane, vOrientation, vOrigin);
            topPlane = PlaneNormalize(topPlane);

            VECTOR bottomPlane = Vector::Set(0.0f, -1.0f, bottomSlope, 0.0f);
            bottomPlane = PlaneTransform(bottomPlane, vOrientation, vOrigin);
            bottomPlane = PlaneNormalize(bottomPlane);

            return TriangleTests::ContainedBy(V0, V1, V2, nearPlane, farPlane, rightPlane, leftPlane, topPlane, bottomPlane);
        }

        _Use_decl_annotations_
        FORCE_INLINE ContainmentType BoundingFrustum::Contains(const BoundingSphere& sh) const noexcept
        {
            // Load origin and orientation of the frustum.
            VECTOR vOrigin = Vector::LoadFloat3(&origin);
            VECTOR vOrientation = Vector::LoadFloat4(&orientation);

            // Create 6 planes (do it inline to encourage use of registers)
            VECTOR nearPlane = Vector::Set(0.0f, 0.0f, -1.0f, near);
            nearPlane = PlaneTransform(nearPlane, vOrientation, vOrigin);
            nearPlane = PlaneNormalize(nearPlane);

            VECTOR farPlane = Vector::Set(0.0f, 0.0f, 1.0f, -far);
            farPlane = PlaneTransform(farPlane, vOrientation, vOrigin);
            farPlane = PlaneNormalize(farPlane);

            VECTOR rightPlane = Vector::Set(1.0f, 0.0f, -rightSlope, 0.0f);
            rightPlane = PlaneTransform(rightPlane, vOrientation, vOrigin);
            rightPlane = PlaneNormalize(rightPlane);

            VECTOR leftPlane = Vector::Set(-1.0f, 0.0f, leftSlope, 0.0f);
            leftPlane = PlaneTransform(leftPlane, vOrientation, vOrigin);
            leftPlane = PlaneNormalize(leftPlane);

            VECTOR topPlane = Vector::Set(0.0f, 1.0f, -topSlope, 0.0f);
            topPlane = PlaneTransform(topPlane, vOrientation, vOrigin);
            topPlane = PlaneNormalize(topPlane);

            VECTOR bottomPlane = Vector::Set(0.0f, -1.0f, bottomSlope, 0.0f);
            bottomPlane = PlaneTransform(bottomPlane, vOrientation, vOrigin);
            bottomPlane = PlaneNormalize(bottomPlane);

            return sh.ContainedBy(nearPlane, farPlane, rightPlane, leftPlane, topPlane, bottomPlane);
        }

        _Use_decl_annotations_
        FORCE_INLINE ContainmentType BoundingFrustum::Contains(const BoundingBox& box) const noexcept
        {
            // Load origin and orientation of the frustum.
            VECTOR vOrigin = Vector::LoadFloat3(&origin);
            VECTOR vOrientation = Vector::LoadFloat4(&orientation);

            // Create 6 planes (do it inline to encourage use of registers)
            VECTOR nearPlane = Vector::Set(0.0f, 0.0f, -1.0f, near);
            nearPlane = PlaneTransform(nearPlane, vOrientation, vOrigin);
            nearPlane = PlaneNormalize(nearPlane);

            VECTOR farPlane = Vector::Set(0.0f, 0.0f, 1.0f, -far);
            farPlane = PlaneTransform(farPlane, vOrientation, vOrigin);
            farPlane = PlaneNormalize(farPlane);

            VECTOR rightPlane = Vector::Set(1.0f, 0.0f, -rightSlope, 0.0f);
            rightPlane = PlaneTransform(rightPlane, vOrientation, vOrigin);
            rightPlane = PlaneNormalize(rightPlane);

            VECTOR leftPlane = Vector::Set(-1.0f, 0.0f, leftSlope, 0.0f);
            leftPlane = PlaneTransform(leftPlane, vOrientation, vOrigin);
            leftPlane = PlaneNormalize(leftPlane);

            VECTOR topPlane = Vector::Set(0.0f, 1.0f, -topSlope, 0.0f);
            topPlane = PlaneTransform(topPlane, vOrientation, vOrigin);
            topPlane = PlaneNormalize(topPlane);

            VECTOR bottomPlane = Vector::Set(0.0f, -1.0f, bottomSlope, 0.0f);
            bottomPlane = PlaneTransform(bottomPlane, vOrientation, vOrigin);
            bottomPlane = PlaneNormalize(bottomPlane);

            return box.ContainedBy(nearPlane, farPlane, rightPlane, leftPlane, topPlane, bottomPlane);
        }

        _Use_decl_annotations_
        FORCE_INLINE ContainmentType BoundingFrustum::Contains(const BoundingOrientedBox& box) const noexcept
        {
            // Load origin and orientation of the frustum.
            VECTOR vOrigin = Vector::LoadFloat3(&origin);
            VECTOR vOrientation = Vector::LoadFloat4(&orientation);

            // Create 6 planes (do it inline to encourage use of registers)
            VECTOR nearPlane = Vector::Set(0.0f, 0.0f, -1.0f, near);
            nearPlane = PlaneTransform(nearPlane, vOrientation, vOrigin);
            nearPlane = PlaneNormalize(nearPlane);

            VECTOR farPlane = Vector::Set(0.0f, 0.0f, 1.0f, -far);
            farPlane = PlaneTransform(farPlane, vOrientation, vOrigin);
            farPlane = PlaneNormalize(farPlane);

            VECTOR rightPlane = Vector::Set(1.0f, 0.0f, -rightSlope, 0.0f);
            rightPlane = PlaneTransform(rightPlane, vOrientation, vOrigin);
            rightPlane = PlaneNormalize(rightPlane);

            VECTOR leftPlane = Vector::Set(-1.0f, 0.0f, leftSlope, 0.0f);
            leftPlane = PlaneTransform(leftPlane, vOrientation, vOrigin);
            leftPlane = PlaneNormalize(leftPlane);

            VECTOR topPlane = Vector::Set(0.0f, 1.0f, -topSlope, 0.0f);
            topPlane = PlaneTransform(topPlane, vOrientation, vOrigin);
            topPlane = PlaneNormalize(topPlane);

            VECTOR bottomPlane = Vector::Set(0.0f, -1.0f, bottomSlope, 0.0f);
            bottomPlane = PlaneTransform(bottomPlane, vOrientation, vOrigin);
            bottomPlane = PlaneNormalize(bottomPlane);

            return box.ContainedBy(nearPlane, farPlane, rightPlane, leftPlane, topPlane, bottomPlane);
        }

        _Use_decl_annotations_
        FORCE_INLINE ContainmentType BoundingFrustum::Contains(const BoundingFrustum& fr) const noexcept
        {
            // Load origin and orientation of the frustum.
            VECTOR vOrigin = Vector::LoadFloat3(&origin);
            VECTOR vOrientation = Vector::LoadFloat4(&orientation);

            // Create 6 planes (do it inline to encourage use of registers)
            VECTOR nearPlane = Vector::Set(0.0f, 0.0f, -1.0f, near);
            nearPlane = PlaneTransform(nearPlane, vOrientation, vOrigin);
            nearPlane = PlaneNormalize(nearPlane);

            VECTOR farPlane = Vector::Set(0.0f, 0.0f, 1.0f, -far);
            farPlane = PlaneTransform(farPlane, vOrientation, vOrigin);
            farPlane = PlaneNormalize(farPlane);

            VECTOR rightPlane = Vector::Set(1.0f, 0.0f, -rightSlope, 0.0f);
            rightPlane = PlaneTransform(rightPlane, vOrientation, vOrigin);
            rightPlane = PlaneNormalize(rightPlane);

            VECTOR leftPlane = Vector::Set(-1.0f, 0.0f, leftSlope, 0.0f);
            leftPlane = PlaneTransform(leftPlane, vOrientation, vOrigin);
            leftPlane = PlaneNormalize(leftPlane);

            VECTOR topPlane = Vector::Set(0.0f, 1.0f, -topSlope, 0.0f);
            topPlane = PlaneTransform(topPlane, vOrientation, vOrigin);
            topPlane = PlaneNormalize(topPlane);

            VECTOR bottomPlane = Vector::Set(0.0f, -1.0f, bottomSlope, 0.0f);
            bottomPlane = PlaneTransform(bottomPlane, vOrientation, vOrigin);
            bottomPlane = PlaneNormalize(bottomPlane);

            return fr.ContainedBy(nearPlane, farPlane, rightPlane, leftPlane, topPlane, bottomPlane);
        }

        _Use_decl_annotations_
        FORCE_INLINE bool BoundingFrustum::Intersects(const BoundingSphere& sh) const noexcept
        {
            VECTOR zero = Vector::Zero();

            // Build the frustum planes.
            VECTOR planes[6];
            planes[0] = Vector::Set(0.0f, 0.0f, -1.0f, near);
            planes[1] = Vector::Set(0.0f, 0.0f, 1.0f, -far);
            planes[2] = Vector::Set(1.0f, 0.0f, -rightSlope, 0.0f);
            planes[3] = Vector::Set(-1.0f, 0.0f, leftSlope, 0.0f);
            planes[4] = Vector::Set(0.0f, 1.0f, -topSlope, 0.0f);
            planes[5] = Vector::Set(0.0f, -1.0f, bottomSlope, 0.0f);

            // Normalize the planes so we can compare to the sphere radius.
            planes[2] = Vector3::Normalize(planes[2]);
            planes[3] = Vector3::Normalize(planes[3]);
            planes[4] = Vector3::Normalize(planes[4]);
            planes[5] = Vector3::Normalize(planes[5]);

            // Load origin and orientation of the frustum.
            VECTOR vOrigin = Vector::LoadFloat3(&origin);
            VECTOR vOrientation = Vector::LoadFloat4(&orientation);

#if defined(DEBUG) || defined(_DEBUG)
            assert(QuaternionIsUnit(vOrientation));
#endif // DEBUG

            // Load the sphere.
            VECTOR vCenter = Vector::LoadFloat3(&sh.center);
            VECTOR vRadius = Vector::ReplicatePtr(&sh.radius);

            // Transform the center of the sphere into the local space of frustum.
            vCenter = Vector3::InverseRotate(Vector::Subtract(vCenter, vOrigin), vOrientation);

            // Set w of the center to one so we can dot4 with the plane.
            vCenter = Vector::Insert<0, 0, 0, 0, 1>(vCenter, Vector::SplatOne());

            // Check against each plane of the frustum.
            VECTOR outside = Vector::FalseInt();
            VECTOR insideAll = Vector::TrueInt();
            VECTOR centerInsideAll = Vector::TrueInt();

            VECTOR dist[6];

            for (size_t i = 0; i < 6; ++i)
            {
                dist[i] = Vector4::Dot(vCenter, planes[i]);

                // Outside the plane?
                outside = Vector::OrInt(outside, Vector::Greater(dist[i], vRadius));

                // Fully inside the plane?
                insideAll = Vector::AndInt(insideAll, Vector::LessOrEqual(dist[i], Vector::Negate(vRadius)));

                // Check if the center is inside the plane.
                centerInsideAll = Vector::AndInt(centerInsideAll, Vector::LessOrEqual(dist[i], zero));
            }

            // If the sphere is outside any of the planes it is outside.
            if (Vector4::EqualInt(outside, Vector::TrueInt()))
                return false;

            // If the sphere is inside all planes it is fully inside.
            if (Vector4::EqualInt(insideAll, Vector::TrueInt()))
                return true;

            // If the center of the sphere is inside all planes and the sphere intersects
            // one or more planes then it must intersect.
            if (Vector4::EqualInt(centerInsideAll, Vector::TrueInt()))
                return true;

            // The sphere may be outside the frustum or intersecting the frustum.
            // Find the nearest feature (face, edge, or corner) on the frustum
            // to the sphere.

            // The faces adjacent to each face are:
            static const size_t adjacent_faces[6][4] =
            {
                { 2, 3, 4, 5 },    // 0
                { 2, 3, 4, 5 },    // 1
                { 0, 1, 4, 5 },    // 2
                { 0, 1, 4, 5 },    // 3
                { 0, 1, 2, 3 },    // 4
                { 0, 1, 2, 3 }
            };  // 5

            VECTOR intersects = Vector::FalseInt();

            // Check to see if the nearest feature is one of the planes.
            for (size_t i = 0; i < 6; ++i)
            {
                // Find the nearest point on the plane to the center of the sphere.
                VECTOR point = Vector::NegativeMultiplySubtract(planes[i], dist[i], vCenter);

                // Set w of the point to one.
                point = Vector::Insert<0, 0, 0, 0, 1>(point, Vector::SplatOne());

                // If the point is inside the face (inside the adjacent planes) then
                // this plane is the nearest feature.
                VECTOR insideFace = Vector::TrueInt();

                for (size_t j = 0; j < 4; j++)
                {
                    size_t plane_index = adjacent_faces[i][j];

                    insideFace = Vector::AndInt(insideFace,
                        Vector::LessOrEqual(Vector4::Dot(point, planes[plane_index]), zero));
                }

                // Since we have already checked distance from the plane we know that the
                // sphere must intersect if this plane is the nearest feature.
                intersects = Vector::OrInt(intersects,
                    Vector::AndInt(Vector::Greater(dist[i], zero), insideFace));
            }

            if (Vector4::EqualInt(intersects, Vector::TrueInt()))
                return true;

            // Build the corners of the frustum.
            VECTOR vRightTop = Vector::Set(rightSlope, topSlope, 1.0f, 0.0f);
            VECTOR vRightBottom = Vector::Set(rightSlope, bottomSlope, 1.0f, 0.0f);
            VECTOR vLeftTop = Vector::Set(leftSlope, topSlope, 1.0f, 0.0f);
            VECTOR vLeftBottom = Vector::Set(leftSlope, bottomSlope, 1.0f, 0.0f);
            VECTOR vNear = Vector::ReplicatePtr(&near);
            VECTOR vFar = Vector::ReplicatePtr(&far);

            VECTOR corners[CORNER_COUNT];
            corners[0] = Vector::Multiply(vRightTop, vNear);
            corners[1] = Vector::Multiply(vRightBottom, vNear);
            corners[2] = Vector::Multiply(vLeftTop, vNear);
            corners[3] = Vector::Multiply(vLeftBottom, vNear);
            corners[4] = Vector::Multiply(vRightTop, vFar);
            corners[5] = Vector::Multiply(vRightBottom, vFar);
            corners[6] = Vector::Multiply(vLeftTop, vFar);
            corners[7] = Vector::Multiply(vLeftBottom, vFar);

            // The Edges are:
            static const size_t edges[12][2] =
            {
                { 0, 1 }, { 2, 3 }, { 0, 2 }, { 1, 3 },    // Near plane
                { 4, 5 }, { 6, 7 }, { 4, 6 }, { 5, 7 },    // Far plane
                { 0, 4 }, { 1, 5 }, { 2, 6 }, { 3, 7 },
            }; // Near to far

            VECTOR radiusSq = Vector::Multiply(vRadius, vRadius);

            // Check to see if the nearest feature is one of the edges (or corners).
            for (size_t i = 0; i < 12; ++i)
            {
                size_t ei0 = edges[i][0];
                size_t ei1 = edges[i][1];

                // Find the nearest point on the edge to the center of the sphere.
                // The corners of the frustum are included as the endpoints of the edges.
                VECTOR point = PointOnLineSegmentNearestPoint(corners[ei0], corners[ei1], vCenter);

                VECTOR delta = Vector::Subtract(vCenter, point);

                VECTOR distSq = Vector3::Dot(delta, delta);

                // If the distance to the center of the sphere to the point is less than
                // the radius of the sphere then it must intersect.
                intersects = Vector::OrInt(intersects, Vector::LessOrEqual(distSq, radiusSq));
            }

            if (Vector4::EqualInt(intersects, Vector::TrueInt()))
                return true;

            // The sphere must be outside the frustum.
            return false;
        }

        _Use_decl_annotations_
        FORCE_INLINE bool BoundingFrustum::Intersects(const BoundingBox& box) const noexcept
        {
            // Make the axis aligned box oriented and do an OBB vs frustum test.
            BoundingOrientedBox obox(box.center, box.extents, Float4(0.f, 0.f, 0.f, 1.f));
            return Intersects(obox);
        }

        _Use_decl_annotations_
        FORCE_INLINE bool BoundingFrustum::Intersects(const BoundingOrientedBox& box) const noexcept
        {
            static const VECTOR_U32 selectY = { { { SELECT_0, SELECT_1, SELECT_0, SELECT_0 } } };
            static const VECTOR_U32 selectZ = { { { SELECT_0, SELECT_0, SELECT_1, SELECT_0 } } };

            VECTOR zero = Vector::Zero();

            // Build the frustum planes.
            VECTOR planes[6];
            planes[0] = Vector::Set(0.0f, 0.0f, -1.0f, near);
            planes[1] = Vector::Set(0.0f, 0.0f, 1.0f, -far);
            planes[2] = Vector::Set(1.0f, 0.0f, -rightSlope, 0.0f);
            planes[3] = Vector::Set(-1.0f, 0.0f, leftSlope, 0.0f);
            planes[4] = Vector::Set(0.0f, 1.0f, -topSlope, 0.0f);
            planes[5] = Vector::Set(0.0f, -1.0f, bottomSlope, 0.0f);

            // Load origin and orientation of the frustum.
            VECTOR vOrigin = Vector::LoadFloat3(&origin);
            VECTOR frustumOrientation = Vector::LoadFloat4(&orientation);

#if defined(DEBUG) || defined(_DEBUG)
            assert(QuaternionIsUnit(frustumOrientation));
#endif // DEBUG

            // Load the box.
            VECTOR center = Vector::LoadFloat3(&box.center);
            VECTOR extents = Vector::LoadFloat3(&box.extents);
            VECTOR boxOrientation = Vector::LoadFloat4(&box.orientation);

#if defined(DEBUG) || defined(_DEBUG)
            assert(QuaternionIsUnit(boxOrientation));
#endif // DEBUG

            // Transform the oriented box into the space of the frustum in order to
            // minimize the number of transforms we have to do.
            center = Vector3::InverseRotate(Vector::Subtract(center, vOrigin), frustumOrientation);
            boxOrientation = Quaternion::Multiply(boxOrientation, Quaternion::Conjugate(frustumOrientation));

            // Set w of the center to one so we can dot4 with the plane.
            center = Vector::Insert<0, 0, 0, 0, 1>(center, Vector::SplatOne());

            // Build the 3x3 rotation matrix that defines the box axes.
            MATRIX R = Matrix::RotationQuaternion(boxOrientation);

            // Check against each plane of the frustum.
            VECTOR outside = Vector::FalseInt();
            VECTOR insideAll = Vector::TrueInt();
            VECTOR centerInsideAll = Vector::TrueInt();

            for (size_t i = 0; i < 6; ++i)
            {
                // Compute the distance to the center of the box.
                VECTOR dist = Vector4::Dot(center, planes[i]);

                // Project the axes of the box onto the normal of the plane.  Half the
                // length of the projection (sometime called the "radius") is equal to
                // h(u) * abs(n dot b(u))) + h(v) * abs(n dot b(v)) + h(w) * abs(n dot b(w))
                // where h(i) are extents of the box, n is the plane normal, and b(i) are the
                // axes of the box.
                VECTOR radius = Vector3::Dot(planes[i], R.r[0]);
                radius = Vector::Select(radius, Vector3::Dot(planes[i], R.r[1]), selectY);
                radius = Vector::Select(radius, Vector3::Dot(planes[i], R.r[2]), selectZ);
                radius = Vector3::Dot(extents, Vector::Abs(radius));

                // Outside the plane?
                outside = Vector::OrInt(outside, Vector::Greater(dist, radius));

                // Fully inside the plane?
                insideAll = Vector::AndInt(insideAll, Vector::LessOrEqual(dist, Vector::Negate(radius)));

                // Check if the center is inside the plane.
                centerInsideAll = Vector::AndInt(centerInsideAll, Vector::LessOrEqual(dist, zero));
            }

            // If the box is outside any of the planes it is outside.
            if (Vector4::EqualInt(outside, Vector::TrueInt()))
                return false;

            // If the box is inside all planes it is fully inside.
            if (Vector4::EqualInt(insideAll, Vector::TrueInt()))
                return true;

            // If the center of the box is inside all planes and the box intersects
            // one or more planes then it must intersect.
            if (Vector4::EqualInt(centerInsideAll, Vector::TrueInt()))
                return true;

            // Build the corners of the frustum.
            VECTOR vRightTop = Vector::Set(rightSlope, topSlope, 1.0f, 0.0f);
            VECTOR vRightBottom = Vector::Set(rightSlope, bottomSlope, 1.0f, 0.0f);
            VECTOR vLeftTop = Vector::Set(leftSlope, topSlope, 1.0f, 0.0f);
            VECTOR vLeftBottom = Vector::Set(leftSlope, bottomSlope, 1.0f, 0.0f);
            VECTOR vNear = Vector::ReplicatePtr(&near);
            VECTOR vFar = Vector::ReplicatePtr(&far);

            VECTOR corners[CORNER_COUNT];
            corners[0] = Vector::Multiply(vRightTop, vNear);
            corners[1] = Vector::Multiply(vRightBottom, vNear);
            corners[2] = Vector::Multiply(vLeftTop, vNear);
            corners[3] = Vector::Multiply(vLeftBottom, vNear);
            corners[4] = Vector::Multiply(vRightTop, vFar);
            corners[5] = Vector::Multiply(vRightBottom, vFar);
            corners[6] = Vector::Multiply(vLeftTop, vFar);
            corners[7] = Vector::Multiply(vLeftBottom, vFar);

            // Test against box axes (3)
            {
                // Find the min/max values of the projection of the frustum onto each axis.
                VECTOR frustumMin, frustumMax;

                frustumMin = Vector3::Dot(corners[0], R.r[0]);
                frustumMin = Vector::Select(frustumMin, Vector3::Dot(corners[0], R.r[1]), selectY);
                frustumMin = Vector::Select(frustumMin, Vector3::Dot(corners[0], R.r[2]), selectZ);
                frustumMax = frustumMin;

                for (size_t i = 1; i < BoundingOrientedBox::CORNER_COUNT; ++i)
                {
                    VECTOR temp = Vector3::Dot(corners[i], R.r[0]);
                    temp = Vector::Select(temp, Vector3::Dot(corners[i], R.r[1]), selectY);
                    temp = Vector::Select(temp, Vector3::Dot(corners[i], R.r[2]), selectZ);

                    frustumMin = Vector::Min(frustumMin, temp);
                    frustumMax = Vector::Max(frustumMax, temp);
                }

                // Project the center of the box onto the axes.
                VECTOR boxDist = Vector3::Dot(center, R.r[0]);
                boxDist = Vector::Select(boxDist, Vector3::Dot(center, R.r[1]), selectY);
                boxDist = Vector::Select(boxDist, Vector3::Dot(center, R.r[2]), selectZ);

                // The projection of the box onto the axis is just its Center and Extents.
                // if (min > box_max || max < box_min) reject;
                VECTOR result = Vector::OrInt(Vector::Greater(frustumMin, Vector::Add(boxDist, extents)),
                    Vector::Less(frustumMax, Vector::Subtract(boxDist, extents)));

                if (Vector3AnyTrue(result))
                    return false;
            }

            // Test against edge/edge axes (3*6).
            VECTOR frustumEdgeAxis[6];

            frustumEdgeAxis[0] = vRightTop;
            frustumEdgeAxis[1] = vRightBottom;
            frustumEdgeAxis[2] = vLeftTop;
            frustumEdgeAxis[3] = vLeftBottom;
            frustumEdgeAxis[4] = Vector::Subtract(vRightTop, vLeftTop);
            frustumEdgeAxis[5] = Vector::Subtract(vLeftBottom, vLeftTop);

            for (size_t i = 0; i < 3; ++i)
            {
                for (size_t j = 0; j < 6; j++)
                {
                    // Compute the axis we are going to test.
                    VECTOR axis = Vector3::Cross(R.r[i], frustumEdgeAxis[j]);

                    // Find the min/max values of the projection of the frustum onto the axis.
                    VECTOR frustumMin, frustumMax;

                    frustumMin = frustumMax = Vector3::Dot(axis, corners[0]);

                    for (size_t k = 1; k < CORNER_COUNT; k++)
                    {
                        VECTOR temp = Vector3::Dot(axis, corners[k]);
                        frustumMin = Vector::Min(frustumMin, temp);
                        frustumMax = Vector::Max(frustumMax, temp);
                    }

                    // Project the center of the box onto the axis.
                    VECTOR dist = Vector3::Dot(center, axis);

                    // Project the axes of the box onto the axis to find the "radius" of the box.
                    VECTOR radius = Vector3::Dot(axis, R.r[0]);
                    radius = Vector::Select(radius, Vector3::Dot(axis, R.r[1]), selectY);
                    radius = Vector::Select(radius, Vector3::Dot(axis, R.r[2]), selectZ);
                    radius = Vector3::Dot(extents, Vector::Abs(radius));

                    // if (center > max + radius || center < min - radius) reject;
                    outside = Vector::OrInt(outside, Vector::Greater(dist, Vector::Add(frustumMax, radius)));
                    outside = Vector::OrInt(outside, Vector::Less(dist, Vector::Subtract(frustumMin, radius)));
                }
            }

            if (Vector4::EqualInt(outside, Vector::TrueInt()))
                return false;

            // If we did not find a separating plane then the box must intersect the frustum.
            return true;
        }

        _Use_decl_annotations_
        FORCE_INLINE bool BoundingFrustum::Intersects(const BoundingFrustum& fr) const noexcept
        {
            // Load origin and orientation of frustum B.
            VECTOR originB = Vector::LoadFloat3(&origin);
            VECTOR orientationB = Vector::LoadFloat4(&orientation);

#if defined(DEBUG) || defined(_DEBUG)
            assert(QuaternionIsUnit(orientationB));
#endif // DEBUG

            // Build the planes of frustum B.
            VECTOR axisB[6];
            axisB[0] = Vector::Set(0.0f, 0.0f, -1.0f, 0.0f);
            axisB[1] = Vector::Set(0.0f, 0.0f, 1.0f, 0.0f);
            axisB[2] = Vector::Set(1.0f, 0.0f, -rightSlope, 0.0f);
            axisB[3] = Vector::Set(-1.0f, 0.0f, leftSlope, 0.0f);
            axisB[4] = Vector::Set(0.0f, 1.0f, -topSlope, 0.0f);
            axisB[5] = Vector::Set(0.0f, -1.0f, bottomSlope, 0.0f);

            VECTOR planeDistB[6];
            planeDistB[0] = Vector::Negate(Vector::ReplicatePtr(&near));
            planeDistB[1] = Vector::ReplicatePtr(&far);
            planeDistB[2] = Vector::Zero();
            planeDistB[3] = Vector::Zero();
            planeDistB[4] = Vector::Zero();
            planeDistB[5] = Vector::Zero();

            // Load origin and orientation of frustum A.
            VECTOR originA = Vector::LoadFloat3(&fr.origin);
            VECTOR orientationA = Vector::LoadFloat4(&fr.orientation);

#if defined(DEBUG) || defined(_DEBUG)
            assert(QuaternionIsUnit(orientationA));
#endif // DEBUG

            // Transform frustum A into the space of the frustum B in order to
            // minimize the number of transforms we have to do.
            originA = Vector3::InverseRotate(Vector::Subtract(originA, originB), orientationB);
            orientationA = Quaternion::Multiply(orientationA, Quaternion::Conjugate(orientationB));

            // Build the corners of frustum A (in the local space of B).
            VECTOR rightTopA = Vector::Set(fr.rightSlope, fr.topSlope, 1.0f, 0.0f);
            VECTOR rightBottomA = Vector::Set(fr.rightSlope, fr.bottomSlope, 1.0f, 0.0f);
            VECTOR leftTopA = Vector::Set(fr.leftSlope, fr.topSlope, 1.0f, 0.0f);
            VECTOR leftBottomA = Vector::Set(fr.leftSlope, fr.bottomSlope, 1.0f, 0.0f);
            VECTOR nearA = Vector::ReplicatePtr(&fr.near);
            VECTOR farA = Vector::ReplicatePtr(&fr.far);

            rightTopA = Vector3::Rotate(rightTopA, orientationA);
            rightBottomA = Vector3::Rotate(rightBottomA, orientationA);
            leftTopA = Vector3::Rotate(leftTopA, orientationA);
            leftBottomA = Vector3::Rotate(leftBottomA, orientationA);

            VECTOR cornersA[CORNER_COUNT];
            cornersA[0] = Vector::MultiplyAdd(rightTopA, nearA, originA);
            cornersA[1] = Vector::MultiplyAdd(rightBottomA, nearA, originA);
            cornersA[2] = Vector::MultiplyAdd(leftTopA, nearA, originA);
            cornersA[3] = Vector::MultiplyAdd(leftBottomA, nearA, originA);
            cornersA[4] = Vector::MultiplyAdd(rightTopA, farA, originA);
            cornersA[5] = Vector::MultiplyAdd(rightBottomA, farA, originA);
            cornersA[6] = Vector::MultiplyAdd(leftTopA, farA, originA);
            cornersA[7] = Vector::MultiplyAdd(leftBottomA, farA, originA);

            // Check frustum A against each plane of frustum B.
            VECTOR outside = Vector::FalseInt();
            VECTOR insideAll = Vector::TrueInt();

            for (size_t i = 0; i < 6; ++i)
            {
                // Find the min/max projection of the frustum onto the plane normal.
                VECTOR min, max;

                min = max = Vector3::Dot(axisB[i], cornersA[0]);

                for (size_t j = 1; j < CORNER_COUNT; j++)
                {
                    VECTOR temp = Vector3::Dot(axisB[i], cornersA[j]);
                    min = Vector::Min(min, temp);
                    max = Vector::Max(max, temp);
                }

                // Outside the plane?
                outside = Vector::OrInt(outside, Vector::Greater(min, planeDistB[i]));

                // Fully inside the plane?
                insideAll = Vector::AndInt(insideAll, Vector::LessOrEqual(max, planeDistB[i]));
            }

            // If the frustum A is outside any of the planes of frustum B it is outside.
            if (Vector4::EqualInt(outside, Vector::TrueInt()))
                return false;

            // If frustum A is inside all planes of frustum B it is fully inside.
            if (Vector4::EqualInt(insideAll, Vector::TrueInt()))
                return true;

            // Build the corners of frustum B.
            VECTOR rightTopB = Vector::Set(rightSlope, topSlope, 1.0f, 0.0f);
            VECTOR rightBottomB = Vector::Set(rightSlope, bottomSlope, 1.0f, 0.0f);
            VECTOR leftTopB = Vector::Set(leftSlope, topSlope, 1.0f, 0.0f);
            VECTOR leftBottomB = Vector::Set(leftSlope, bottomSlope, 1.0f, 0.0f);
            VECTOR nearB = Vector::ReplicatePtr(&near);
            VECTOR farB = Vector::ReplicatePtr(&far);

            VECTOR cornersB[BoundingFrustum::CORNER_COUNT];
            cornersB[0] = Vector::Multiply(rightTopB, nearB);
            cornersB[1] = Vector::Multiply(rightBottomB, nearB);
            cornersB[2] = Vector::Multiply(leftTopB, nearB);
            cornersB[3] = Vector::Multiply(leftBottomB, nearB);
            cornersB[4] = Vector::Multiply(rightTopB, farB);
            cornersB[5] = Vector::Multiply(rightBottomB, farB);
            cornersB[6] = Vector::Multiply(leftTopB, farB);
            cornersB[7] = Vector::Multiply(leftBottomB, farB);

            // Build the planes of frustum A (in the local space of B).
            VECTOR axisA[6];
            VECTOR planeDistA[6];

            axisA[0] = Vector::Set(0.0f, 0.0f, -1.0f, 0.0f);
            axisA[1] = Vector::Set(0.0f, 0.0f, 1.0f, 0.0f);
            axisA[2] = Vector::Set(1.0f, 0.0f, -fr.rightSlope, 0.0f);
            axisA[3] = Vector::Set(-1.0f, 0.0f, fr.leftSlope, 0.0f);
            axisA[4] = Vector::Set(0.0f, 1.0f, -fr.topSlope, 0.0f);
            axisA[5] = Vector::Set(0.0f, -1.0f, fr.bottomSlope, 0.0f);

            axisA[0] = Vector3::Rotate(axisA[0], orientationA);
            axisA[1] = Vector::Negate(axisA[0]);
            axisA[2] = Vector3::Rotate(axisA[2], orientationA);
            axisA[3] = Vector3::Rotate(axisA[3], orientationA);
            axisA[4] = Vector3::Rotate(axisA[4], orientationA);
            axisA[5] = Vector3::Rotate(axisA[5], orientationA);

            planeDistA[0] = Vector3::Dot(axisA[0], cornersA[0]);  // Re-use corner on near plane.
            planeDistA[1] = Vector3::Dot(axisA[1], cornersA[4]);  // Re-use corner on far plane.
            planeDistA[2] = Vector3::Dot(axisA[2], originA);
            planeDistA[3] = Vector3::Dot(axisA[3], originA);
            planeDistA[4] = Vector3::Dot(axisA[4], originA);
            planeDistA[5] = Vector3::Dot(axisA[5], originA);

            // Check each axis of frustum A for a seperating plane (5).
            for (size_t i = 0; i < 6; ++i)
            {
                // Find the minimum projection of the frustum onto the plane normal.
                VECTOR min;

                min = Vector3::Dot(axisA[i], cornersB[0]);

                for (size_t j = 1; j < CORNER_COUNT; j++)
                {
                    VECTOR temp = Vector3::Dot(axisA[i], cornersB[j]);
                    min = Vector::Min(min, temp);
                }

                // Outside the plane?
                outside = Vector::OrInt(outside, Vector::Greater(min, planeDistA[i]));
            }

            // If the frustum B is outside any of the planes of frustum A it is outside.
            if (Vector4::EqualInt(outside, Vector::TrueInt()))
                return false;

            // Check edge/edge axes (6 * 6).
            VECTOR frustumEdgeAxisA[6];
            frustumEdgeAxisA[0] = rightTopA;
            frustumEdgeAxisA[1] = rightBottomA;
            frustumEdgeAxisA[2] = leftTopA;
            frustumEdgeAxisA[3] = leftBottomA;
            frustumEdgeAxisA[4] = Vector::Subtract(rightTopA, leftTopA);
            frustumEdgeAxisA[5] = Vector::Subtract(leftBottomA, leftTopA);

            VECTOR frustumEdgeAxisB[6];
            frustumEdgeAxisB[0] = rightTopB;
            frustumEdgeAxisB[1] = rightBottomB;
            frustumEdgeAxisB[2] = leftTopB;
            frustumEdgeAxisB[3] = leftBottomB;
            frustumEdgeAxisB[4] = Vector::Subtract(rightTopB, leftTopB);
            frustumEdgeAxisB[5] = Vector::Subtract(leftBottomB, leftTopB);

            for (size_t i = 0; i < 6; ++i)
            {
                for (size_t j = 0; j < 6; j++)
                {
                    // Compute the axis we are going to test.
                    VECTOR axis = Vector3::Cross(frustumEdgeAxisA[i], frustumEdgeAxisB[j]);

                    // Find the min/max values of the projection of both frustums onto the axis.
                    VECTOR minA, maxA;
                    VECTOR minB, maxB;

                    minA = maxA = Vector3::Dot(axis, cornersA[0]);
                    minB = maxB = Vector3::Dot(axis, cornersB[0]);

                    for (size_t k = 1; k < CORNER_COUNT; k++)
                    {
                        VECTOR tempA = Vector3::Dot(axis, cornersA[k]);
                        minA = Vector::Min(minA, tempA);
                        maxA = Vector::Max(maxA, tempA);

                        VECTOR tempB = Vector3::Dot(axis, cornersB[k]);
                        minB = Vector::Min(minB, tempB);
                        maxB = VE::Max(maxB, tempB);
                    }

                    // if (MinA > MaxB || MinB > MaxA) reject
                    outside = Vector::OrInt(outside, Vector::Greater(minA, maxB));
                    outside = Vector::OrInt(outside, Vector::Greater(minB, maxA));
                }
            }

            // If there is a seperating plane, then the frustums do not intersect.
            if (Vector4::EqualInt(outside, Vector::TrueInt()))
                return false;

            // If we did not find a separating plane then the frustums intersect.
            return true;
        }

        _Use_decl_annotations_
        FORCE_INLINE bool VEC_CALLCONV BoundingFrustum::Intersects(A_VECTOR V0, A_VECTOR V1, A_VECTOR V2) const noexcept
        {
            // Build the frustum planes (NOTE: D is negated from the usual).
            VECTOR planes[6];
            planes[0] = Vector::Set(0.0f, 0.0f, -1.0f, -near);
            planes[1] = Vector::Set(0.0f, 0.0f, 1.0f, far);
            planes[2] = Vector::Set(1.0f, 0.0f, -rightSlope, 0.0f);
            planes[3] = Vector::Set(-1.0f, 0.0f, leftSlope, 0.0f);
            planes[4] = Vector::Set(0.0f, 1.0f, -topSlope, 0.0f);
            planes[5] = Vector::Set(0.0f, -1.0f, bottomSlope, 0.0f);

            // Load origin and orientation of the frustum.
            VECTOR vOrigin = Vector::LoadFloat3(&origin);
            VECTOR vOrientation = Vector::LoadFloat4(&orientation);

#if defined(DEBUG) || defined(_DEBUG)
            assert(QuaternionIsUnit(vOrientation));
#endif // DEBUG

            // Transform triangle into the local space of frustum.
            VECTOR TV0 = Vector3::InverseRotate(Vector::Subtract(V0, vOrigin), vOrientation);
            VECTOR TV1 = Vector3::InverseRotate(Vector::Subtract(V1, vOrigin), vOrientation);
            VECTOR TV2 = Vector3::InverseRotate(Vector::Subtract(V2, vOrigin), vOrientation);

            // Test each vertex of the triangle against the frustum planes.
            VECTOR outside = Vector::FalseInt();
            VECTOR insideAll = Vector::TrueInt();

            for (size_t i = 0; i < 6; ++i)
            {
                VECTOR dist0 = Vector3::Dot(TV0, planes[i]);
                VECTOR dist1 = Vector3::Dot(TV1, planes[i]);
                VECTOR dist2 = Vector3::Dot(TV2, planes[i]);

                VECTOR minDist = Vector::Min(dist0, dist1);
                minDist = Vector::Min(minDist, dist2);
                VECTOR maxDist = Vector::Max(dist0, dist1);
                maxDist = Vector::Max(maxDist, dist2);

                VECTOR planeDist = Vector::SplatW(planes[i]);

                // Outside the plane?
                outside = Vector::OrInt(outside, Vector::Greater(minDist, planeDist));

                // Fully inside the plane?
                insideAll = Vector::AndInt(insideAll, Vector::LessOrEqual(maxDist, planeDist));
            }

            // If the triangle is outside any of the planes it is outside.
            if (Vector4::EqualInt(outside, Vector::TrueInt()))
                return false;

            // If the triangle is inside all planes it is fully inside.
            if (Vector4::EqualInt(insideAll, Vector::TrueInt()))
                return true;

            // Build the corners of the frustum.
            VECTOR vRightTop = Vector::Set(rightSlope, topSlope, 1.0f, 0.0f);
            VECTOR vRightBottom = Vector::Set(rightSlope, bottomSlope, 1.0f, 0.0f);
            VECTOR vLeftTop = Vector::Set(leftSlope, topSlope, 1.0f, 0.0f);
            VECTOR vLeftBottom = Vector::Set(leftSlope, bottomSlope, 1.0f, 0.0f);
            VECTOR vNear = Vector::ReplicatePtr(&near);
            VECTOR vFar = Vector::ReplicatePtr(&far);

            VECTOR corners[CORNER_COUNT];
            corners[0] = Vector::Multiply(vRightTop, vNear);
            corners[1] = Vector::Multiply(vRightBottom, vNear);
            corners[2] = Vector::Multiply(vLeftTop, vNear);
            corners[3] = Vector::Multiply(vLeftBottom, vNear);
            corners[4] = Vector::Multiply(vRightTop, vFar);
            corners[5] = Vector::Multiply(vRightBottom, vFar);
            corners[6] = Vector::Multiply(vLeftTop, vFar);
            corners[7] = Vector::Multiply(vLeftBottom, vFar);

            // Test the plane of the triangle.
            VECTOR normal = Vector3::Cross(Vector::Subtract(V1, V0), Vector::Subtract(V2, V0));
            VECTOR dist = Vector3::Dot(normal, V0);

            VECTOR minDist, maxDist;
            minDist = maxDist = Vector3::Dot(corners[0], normal);
            for (size_t i = 1; i < CORNER_COUNT; ++i)
            {
                VECTOR temp = Vector3::Dot(corners[i], normal);
                minDist = Vector::Min(minDist, temp);
                maxDist = Vector::Max(maxDist, temp);
            }

            outside = Vector::OrInt(Vector::Greater(minDist, dist), Vector::Less(maxDist, dist));
            if (Vector4::EqualInt(outside, Vector::TrueInt()))
                return false;

            // Check the edge/edge axes (3*6).
            VECTOR triangleEdgeAxis[3];
            triangleEdgeAxis[0] = Vector::Subtract(V1, V0);
            triangleEdgeAxis[1] = Vector::Subtract(V2, V1);
            triangleEdgeAxis[2] = Vector::Subtract(V0, V2);

            VECTOR frustumEdgeAxis[6];
            frustumEdgeAxis[0] = vRightTop;
            frustumEdgeAxis[1] = vRightBottom;
            frustumEdgeAxis[2] = vLeftTop;
            frustumEdgeAxis[3] = vLeftBottom;
            frustumEdgeAxis[4] = Vector::Subtract(vRightTop, vLeftTop);
            frustumEdgeAxis[5] = Vector::Subtract(vLeftBottom, vLeftTop);

            for (size_t i = 0; i < 3; ++i)
            {
                for (size_t j = 0; j < 6; j++)
                {
                    // Compute the axis we are going to test.
                    VECTOR axis = Vector3::Cross(triangleEdgeAxis[i], frustumEdgeAxis[j]);

                    // Find the min/max of the projection of the triangle onto the axis.
                    VECTOR minA, maxA;

                    VECTOR dist0 = Vector3::Dot(V0, axis);
                    VECTOR dist1 = Vector3::Dot(V1, axis);
                    VECTOR dist2 = Vector3::Dot(V2, axis);

                    minA = Vector::Min(dist0, dist1);
                    minA = Vector::Min(minA, dist2);
                    maxA = Vector::Max(dist0, dist1);
                    maxA = Vector::Max(maxA, dist2);

                    // Find the min/max of the projection of the frustum onto the axis.
                    VECTOR minB, maxB;

                    minB = maxB = Vector3::Dot(axis, corners[0]);

                    for (size_t k = 1; k < CORNER_COUNT; k++)
                    {
                        VECTOR temp = Vector3::Dot(axis, corners[k]);
                        minB = Vector::Min(minB, temp);
                        maxB = Vector::Max(maxB, temp);
                    }

                    // if (MinA > MaxB || MinB > MaxA) reject;
                    outside = Vector::OrInt(outside, Vector::Greater(minA, maxB));
                    outside = Vector::OrInt(outside, Vector::Greater(minB, maxA));
                }
            }

            if (Vector4::EqualInt(outside, Vector::TrueInt()))
                return false;

            // If we did not find a separating plane then the triangle must intersect the frustum.
            return true;
        }

        _Use_decl_annotations_
        FORCE_INLINE PlaneIntersectionType VEC_CALLCONV BoundingFrustum::Intersects(A_VECTOR plane) const noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(PlaneIsUnit(plane));
#endif // DEBUG

            // Load origin and orientation of the frustum.
            VECTOR vOrigin = Vector::LoadFloat3(&origin);
            VECTOR vOrientation = Vector::LoadFloat4(&orientation);

#if defined(DEBUG) || defined(_DEBUG)
            assert(QuaternionIsUnit(vOrientation));
#endif // DEBUG

            // Set w of the origin to one so we can dot4 with a plane.
            vOrigin = Vector::Insert<0, 0, 0, 0, 1>(vOrigin, Vector::SplatOne());

            // Build the corners of the frustum (in world space).
            VECTOR rightTop = Vector::Set(rightSlope, topSlope, 1.0f, 0.0f);
            VECTOR rightBottom = Vector::Set(rightSlope, bottomSlope, 1.0f, 0.0f);
            VECTOR leftTop = Vector::Set(leftSlope, topSlope, 1.0f, 0.0f);
            VECTOR leftBottom = Vector::Set(leftSlope, bottomSlope, 1.0f, 0.0f);
            VECTOR vNear = Vector::ReplicatePtr(&near);
            VECTOR vFar = Vector::ReplicatePtr(&far);

            rightTop = Vector3::Rotate(rightTop, vOrientation);
            rightBottom = Vector3::Rotate(rightBottom, vOrientation);
            leftTop = Vector3::Rotate(leftTop, vOrientation);
            leftBottom = Vector3::Rotate(leftBottom, vOrientation);

            VECTOR corners0 = Vector::MultiplyAdd(rightTop, vNear, vOrigin);
            VECTOR corners1 = Vector::MultiplyAdd(rightBottom, vNear, vOrigin);
            VECTOR corners2 = Vector::MultiplyAdd(leftTop, vNear, vOrigin);
            VECTOR corners3 = Vector::MultiplyAdd(leftBottom, vNear, vOrigin);
            VECTOR corners4 = Vector::MultiplyAdd(rightTop, vFar, vOrigin);
            VECTOR corners5 = Vector::MultiplyAdd(rightBottom, vFar, vOrigin);
            VECTOR corners6 = Vector::MultiplyAdd(leftTop, vFar, vOrigin);
            VECTOR corners7 = Vector::MultiplyAdd(leftBottom, vFar, vOrigin);

            VECTOR outside, inside;
            FastIntersectFrustumPlane(corners0, corners1, corners2, corners3,
                corners4, corners5, corners6, corners7,
                plane, outside, inside);

            // If the frustum is outside any plane it is outside.
            if (Vector4::EqualInt(outside, Vector::TrueInt()))
                return FRONT;

            // If the frustum is inside all planes it is inside.
            if (Vector4::EqualInt(inside, Vector::TrueInt()))
                return BACK;

            // The frustum is not inside all planes or outside a plane it intersects.
            return INTERSECTING;
        }

        _Use_decl_annotations_
        FORCE_INLINE bool VEC_CALLCONV BoundingFrustum::Intersects(A_VECTOR rayOrigin, A_VECTOR direction, float& dist) const noexcept
        {
            // If ray starts inside the frustum, return a distance of 0 for the hit
            if (Contains(rayOrigin) == CONTAINS)
            {
                dist = 0.0f;
                return true;
            }

            // Build the frustum planes.
            VECTOR planes[6];
            planes[0] = Vector::Set(0.0f, 0.0f, -1.0f, near);
            planes[1] = Vector::Set(0.0f, 0.0f, 1.0f, -far);
            planes[2] = Vector::Set(1.0f, 0.0f, -rightSlope, 0.0f);
            planes[3] = Vector::Set(-1.0f, 0.0f, leftSlope, 0.0f);
            planes[4] = Vector::Set(0.0f, 1.0f, -topSlope, 0.0f);
            planes[5] = Vector::Set(0.0f, -1.0f, bottomSlope, 0.0f);

            // Load origin and orientation of the frustum.
            VECTOR frOrigin = Vector::LoadFloat3(&origin);
            VECTOR frOrientation = Vector::LoadFloat4(&orientation);

            // This algorithm based on "Fast Ray-Convex Polyhedron Intersectin," in James Arvo, ed., Graphics Gems II pp. 247-250
            float tnear = -FLT_MAX;
            float tfar = FLT_MAX;

            for (size_t i = 0; i < 6; ++i)
            {
                VECTOR plane = PlaneTransform(planes[i], frOrientation, frOrigin);
                plane = PlaneNormalize(plane);

                VECTOR axisDotOrigin = PlaneDotCoord(plane, rayOrigin);
                VECTOR axisDotDirection = Vector3::Dot(plane, direction);

                if (Vector3::LessOrEqual(Vector::Abs(axisDotDirection), g_RayEpsilon))
                {
                    // Ray is parallel to plane - check if ray origin is inside plane's
                    if (Vector3::Greater(axisDotOrigin, g_Zero))
                    {
                        // Ray origin is outside half-space.
                        dist = 0.f;
                        return false;
                    }
                }
                else
                {
                    // Ray not parallel - get distance to plane.
                    float vd = Vector::GetX(axisDotDirection);
                    float vn = Vector::GetX(axisDotOrigin);
                    float t = -vn / vd;
                    if (vd < 0.0f)
                    {
                        // Front face - T is a near point.
                        if (t > tfar)
                        {
                            dist = 0.f;
                            return false;
                        }
                        if (t > tnear)
                        {
                            // Hit near face.
                            tnear = t;
                        }
                    }
                    else
                    {
                        // back face - T is far point.
                        if (t < tnear)
                        {
                            dist = 0.f;
                            return false;
                        }
                        if (t < tfar)
                        {
                            // Hit far face.
                            tfar = t;
                        }
                    }
                }
            }

            // Survived all tests.
            // Note: if ray originates on polyhedron, may want to change 0.0f to some
            // epsilon to avoid intersecting the originating face.
            float distance = (tnear >= 0.0f) ? tnear : tfar;
            if (distance >= 0.0f)
            {
                dist = distance;
                return true;
            }

            dist = 0.f;
            return false;
        }

        _Use_decl_annotations_
        FORCE_INLINE ContainmentType VEC_CALLCONV BoundingFrustum::ContainedBy(
            A_VECTOR plane0, A_VECTOR plane1, A_VECTOR plane2,
            B_VECTOR plane3,
            C_VECTOR plane4, C_VECTOR plane5) const noexcept
        {
            // Load origin and orientation of the frustum.
            VECTOR vOrigin = Vector::LoadFloat3(&origin);
            VECTOR vOrientation = Vector::LoadFloat4(&orientation);

#if defined(DEBUG) || defined(_DEBUG)
            assert(QuaternionIsUnit(vOrientation));
#endif // DEBUG

            // Set w of the origin to one so we can dot4 with a plane.
            vOrigin = Vector::Insert<0, 0, 0, 0, 1>(vOrigin, Vector::SplatOne());

            // Build the corners of the frustum (in world space).
            VECTOR rightTop = Vector::Set(rightSlope, topSlope, 1.0f, 0.0f);
            VECTOR rightBottom = Vector::Set(rightSlope, bottomSlope, 1.0f, 0.0f);
            VECTOR leftTop = Vector::Set(leftSlope, topSlope, 1.0f, 0.0f);
            VECTOR leftBottom = Vector::Set(leftSlope, bottomSlope, 1.0f, 0.0f);
            VECTOR vNear = Vector::ReplicatePtr(&near);
            VECTOR vFar = Vector::ReplicatePtr(&far);

            rightTop = Vector3::Rotate(rightTop, vOrientation);
            rightBottom = Vector3::Rotate(rightBottom, vOrientation);
            leftTop = Vector3::Rotate(leftTop, vOrientation);
            leftBottom = Vector3::Rotate(leftBottom, vOrientation);

            VECTOR corners0 = Vector::MultiplyAdd(rightTop, vNear, vOrigin);
            VECTOR corners1 = Vector::MultiplyAdd(rightBottom, vNear, vOrigin);
            VECTOR corners2 = Vector::MultiplyAdd(leftTop, vNear, vOrigin);
            VECTOR corners3 = Vector::MultiplyAdd(leftBottom, vNear, vOrigin);
            VECTOR corners4 = Vector::MultiplyAdd(rightTop, vFar, vOrigin);
            VECTOR corners5 = Vector::MultiplyAdd(rightBottom, vFar, vOrigin);
            VECTOR corners6 = Vector::MultiplyAdd(leftTop, vFar, vOrigin);
            VECTOR corners7 = Vector::MultiplyAdd(leftBottom, vFar, vOrigin);

            VECTOR outside, inside;

            // Test against each plane.
            FastIntersectFrustumPlane(corners0, corners1, corners2, corners3,
                corners4, corners5, corners6, corners7,
                plane0, outside, inside);

            VECTOR anyOutside = outside;
            VECTOR allInside = inside;

            FastIntersectFrustumPlane(corners0, corners1, corners2, corners3,
                corners4, corners5, corners6, corners7,
                plane1, outside, inside);

            anyOutside = Vector::OrInt(anyOutside, outside);
            allInside = Vector::AndInt(allInside, inside);

            FastIntersectFrustumPlane(corners0, corners1, corners2, corners3,
                corners4, corners5, corners6, corners7,
                plane2, outside, inside);

            anyOutside = Vector::OrInt(anyOutside, outside);
            allInside = Vector::AndInt(allInside, inside);

            FastIntersectFrustumPlane(corners0, corners1, corners2, corners3,
                corners4, corners5, corners6, corners7,
                plane3, outside, inside);

            anyOutside = Vector::OrInt(anyOutside, outside);
            allInside = Vector::AndInt(allInside, inside);

            FastIntersectFrustumPlane(corners0, corners1, corners2, corners3,
                corners4, corners5, corners6, corners7,
                plane4, outside, inside);

            anyOutside = Vector::OrInt(anyOutside, outside);
            allInside = Vector::AndInt(allInside, inside);

            FastIntersectFrustumPlane(corners0, corners1, corners2, corners3,
                corners4, corners5, corners6, corners7,
                plane5, outside, inside);

            anyOutside = Vector::OrInt(anyOutside, outside);
            allInside = Vector::AndInt(allInside, inside);

            // If the frustum is outside any plane it is outside.
            if (Vector4::EqualInt(anyOutside, Vector::TrueInt()))
                return DISJOINT;

            // If the frustum is inside all planes it is inside.
            if (Vector4::EqualInt(allInside, Vector::TrueInt()))
                return CONTAINS;

            // The frustum is not inside all planes or outside a plane, it may intersect.
            return INTERSECTS;
        }

        _Use_decl_annotations_
        FORCE_INLINE void BoundingFrustum::GetPlanes(VECTOR* nearPlane, VECTOR* farPlane, VECTOR* rightPlane,
            VECTOR* leftPlane, VECTOR* topPlane, VECTOR* bottomPlane) const noexcept
        {
            // Load origin and orientation of the frustum.
            VECTOR vOrigin = Vector::LoadFloat3(&origin);
            VECTOR vOrientation = Vector::LoadFloat4(&orientation);

            if (nearPlane)
            {
                VECTOR vNearPlane = Vector::Set(0.0f, 0.0f, -1.0f, near);
                vNearPlane = PlaneTransform(vNearPlane, vOrientation, vOrigin);
                *nearPlane = PlaneNormalize(vNearPlane);
            }

            if (farPlane)
            {
                VECTOR vFarPlane = Vector::Set(0.0f, 0.0f, 1.0f, -far);
                vFarPlane = PlaneTransform(vFarPlane, vOrientation, vOrigin);
                *farPlane = PlaneNormalize(vFarPlane);
            }

            if (rightPlane)
            {
                VECTOR vRightPlane = Vector::Set(1.0f, 0.0f, -rightSlope, 0.0f);
                vRightPlane = PlaneTransform(vRightPlane, vOrientation, vOrigin);
                *rightPlane = PlaneNormalize(vRightPlane);
            }

            if (leftPlane)
            {
                VECTOR vLeftPlane = Vector::Set(-1.0f, 0.0f, leftSlope, 0.0f);
                vLeftPlane = PlaneTransform(vLeftPlane, vOrientation, vOrigin);
                *leftPlane = PlaneNormalize(vLeftPlane);
            }

            if (topPlane)
            {
                VECTOR vTopPlane = Vector::Set(0.0f, 1.0f, -topSlope, 0.0f);
                vTopPlane = PlaneTransform(vTopPlane, vOrientation, vOrigin);
                *topPlane = PlaneNormalize(vTopPlane);
            }

            if (bottomPlane)
            {
                VECTOR vBottomPlane = Vector::Set(0.0f, -1.0f, bottomSlope, 0.0f);
                vBottomPlane = PlaneTransform(vBottomPlane, vOrientation, vOrigin);
                *bottomPlane = PlaneNormalize(vBottomPlane);
            }
        }

        _Use_decl_annotations_
        FORCE_INLINE void VEC_CALLCONV BoundingFrustum::CreateFromMatrix(BoundingFrustum& out, A_MATRIX projection, bool rhcoords) noexcept
        {
            // Corners of the projection frustum in NDC space.
            static VECTOR_F32 NDCPoints[6] =
            {
                { { {  1.0f,  0.0f, 1.0f, 1.0f } } },   // right (at far plane)
                { { { -1.0f,  0.0f, 1.0f, 1.0f } } },   // left
                { { {  0.0f,  1.0f, 1.0f, 1.0f } } },   // top
                { { {  0.0f, -1.0f, 1.0f, 1.0f } } },   // bottom

                { { { 0.0f, 0.0f, 0.0f, 1.0f } } },     // near
                { { { 0.0f, 0.0f, 1.0f, 1.0f } } }      // far
            };

            VECTOR determinant;
            MATRIX matInverse = Matrix::Inverse(&determinant, projection);

            // Compute the frustum corners in world space.
            VECTOR points[6];

            for (size_t i = 0; i < 6; ++i)
            {
                // Transform point.
                points[i] = Vector4::Transform(NDCPoints[i], matInverse);
            }

            out.origin = Float3(0.0f, 0.0f, 0.0f);
            out.orientation = Float4(0.0f, 0.0f, 0.0f, 1.0f);

            // Compute the slopes.
            points[0] = Vector::Multiply(points[0], Vector::Reciprocal(Vector::SplatZ(points[0])));
            points[1] = Vector::Multiply(points[1], Vector::Reciprocal(Vector::SplatZ(points[1])));
            points[2] = Vector::Multiply(points[2], Vector::Reciprocal(Vector::SplatZ(points[2])));
            points[3] = Vector::Multiply(points[3], Vector::Reciprocal(Vector::SplatZ(points[3])));

            out.rightSlope = Vector::GetX(points[0]);
            out.leftSlope = Vector::GetX(points[1]);
            out.topSlope = Vector::GetY(points[2]);
            out.bottomSlope = Vector::GetY(points[3]);

            // Compute near and far.
            points[4] = Vector::Multiply(points[4], Vector::Reciprocal(Vector::SplatW(points[4])));
            points[5] = Vector::Multiply(points[5], Vector::Reciprocal(Vector::SplatW(points[5])));

            if (rhcoords)
            {
                out.near = Vector::GetZ(points[5]);
                out.far = Vector::GetZ(points[4]);
            }
            else
            {
                out.near = Vector::GetZ(points[4]);
                out.far = Vector::GetZ(points[5]);
            }
        }

        namespace TriangleTests
        {

            //-----------------------------------------------------------------------------
            // Compute the intersection of a ray (origin, direction) with a triangle
            // (V0, V1, V2).  Return true if there is an intersection and also set *pDist
            // to the distance along the ray to the intersection.
            //
            // The algorithm is based on Moller, Tomas and Trumbore, "Fast, Minimum Storage
            // Ray-Triangle Intersection", Journal of Graphics Tools, vol. 2, no. 1,
            // pp 21-28, 1997.
            //-----------------------------------------------------------------------------
            _Use_decl_annotations_
            FORCE_INLINE bool VEC_CALLCONV Intersects(
                A_VECTOR origin, A_VECTOR direction, A_VECTOR V0,
                B_VECTOR V1,
                C_VECTOR V2, float& dist) noexcept
            {
#if defined(DEBUG) || defined(_DEBUG)
                assert(Vector3IsUnit(direction));
#endif // DEBUG

                VECTOR zero = Vector::Zero();

                VECTOR e1 = Vector::Subtract(V1, V0);
                VECTOR e2 = Vector::Subtract(V2, V0);

                // p = Direction ^ e2;
                VECTOR p = Vector3::Cross(direction, e2);

                // det = e1 * p;
                VECTOR det = Vector3::Dot(e1, p);

                VECTOR u, v, t;

                if (Vector3::GreaterOrEqual(det, g_RayEpsilon))
                {
                    // Determinate is positive (front side of the triangle).
                    VECTOR s = Vector::Subtract(origin, V0);

                    // u = s * p;
                    u = Vector3::Dot(s, p);

                    VECTOR noIntersection = Vector::Less(u, zero);
                    noIntersection = Vector::OrInt(noIntersection, Vector::Greater(u, det));

                    // q = s ^ e1;
                    VECTOR q = Vector3::Cross(s, e1);

                    // v = Direction * q;
                    v = Vector3::Dot(direction, q);

                    noIntersection = Vector::OrInt(noIntersection, Vector::Less(v, zero));
                    noIntersection = Vector::OrInt(noIntersection, Vector::Greater(Vector::Add(u, v), det));

                    // t = e2 * q;
                    t = Vector3::Dot(e2, q);

                    noIntersection = Vector::OrInt(noIntersection, Vector::Less(t, zero));

                    if (Vector4::EqualInt(noIntersection, Vector::TrueInt()))
                    {
                        dist = 0.f;
                        return false;
                    }
                }
                else if (Vector3::LessOrEqual(det, g_RayNegEpsilon))
                {
                    // Determinate is negative (back side of the triangle).
                    VECTOR s = Vector::Subtract(origin, V0);

                    // u = s * p;
                    u = Vector3::Dot(s, p);

                    VECTOR noIntersection = Vector::Greater(u, zero);
                    noIntersection = Vector::OrInt(noIntersection, Vector::Less(u, det));

                    // q = s ^ e1;
                    VECTOR q = Vector3::Cross(s, e1);

                    // v = Direction * q;
                    v = Vector3::Dot(direction, q);

                    noIntersection = Vector::OrInt(noIntersection, Vector::Greater(v, zero));
                    noIntersection = Vector::OrInt(noIntersection, Vector::Less(Vector::Add(u, v), det));

                    // t = e2 * q;
                    t = Vector3::Dot(e2, q);

                    noIntersection = Vector::OrInt(noIntersection, Vector::Greater(t, zero));

                    if (Vector4::EqualInt(noIntersection, Vector::TrueInt()))
                    {
                        dist = 0.f;
                        return false;
                    }
                }
                else
                {
                    // Parallel ray.
                    dist = 0.f;
                    return false;
                }

                t = Vector::Divide(t, det);

                // (u / det) and (v / dev) are the barycentric cooridinates of the intersection.

                // Store the x-component to *pDist
                StoreFloat(&dist, t);

                return true;
            }


            //-----------------------------------------------------------------------------
            // Test if two triangles intersect.
            //
            // The final test of algorithm is based on Shen, Heng, and Tang, "A Fast
            // Triangle-Triangle Overlap Test Using Signed Distances", Journal of Graphics
            // Tools, vol. 8, no. 1, pp 17-23, 2003 and Guigue and Devillers, "Fast and
            // Robust Triangle-Triangle Overlap Test Using Orientation Predicates", Journal
            // of Graphics Tools, vol. 8, no. 1, pp 25-32, 2003.
            //
            // The final test could be considered an edge-edge separating plane test with
            // the 9 possible cases narrowed down to the only two pairs of edges that can
            // actaully result in a seperation.
            //-----------------------------------------------------------------------------
            _Use_decl_annotations_
            FORCE_INLINE bool VEC_CALLCONV Intersects(A_VECTOR A0, A_VECTOR A1, A_VECTOR A2, B_VECTOR B0, C_VECTOR B1, C_VECTOR B2) noexcept
            {
                static const VECTOR_U32 selectY = { { { SELECT_0, SELECT_1, SELECT_0, SELECT_0 } } };
                static const VECTOR_U32 selectZ = { { { SELECT_0, SELECT_0, SELECT_1, SELECT_0 } } };
                static const VECTOR_U32 select0111 = { { { SELECT_0, SELECT_1, SELECT_1, SELECT_1 } } };
                static const VECTOR_U32 select1011 = { { { SELECT_1, SELECT_0, SELECT_1, SELECT_1 } } };
                static const VECTOR_U32 select1101 = { { { SELECT_1, SELECT_1, SELECT_0, SELECT_1 } } };

                VECTOR zero = Vector::Zero();

                // Compute the normal of triangle A.
                VECTOR N1 = Vector3::Cross(Vector::Subtract(A1, A0), Vector::Subtract(A2, A0));

#if defined(DEBUG) || defined(_DEBUG)
                // Assert that the triangle is not degenerate.
                assert(!Vector3::Equal(N1, zero));
#endif // DEBUG

                // Test points of B against the plane of A.
                VECTOR BDist = Vector3::Dot(N1, Vector::Subtract(B0, A0));
                BDist = Vector::Select(BDist, Vector3::Dot(N1, Vector::Subtract(B1, A0)), selectY);
                BDist = Vector::Select(BDist, Vector3::Dot(N1, Vector::Subtract(B2, A0)), selectZ);

                // Ensure robustness with co-planar triangles by zeroing small distances.
                uint32_t BDistIsZeroCR;
                VECTOR BDistIsZero = GreaterR(&BDistIsZeroCR, g_RayEpsilon, Vector::Abs(BDist));
                BDist = Vector::Select(BDist, zero, BDistIsZero);

                uint32_t BDistIsLessCR;
                VECTOR BDistIsLess = Vector::GreaterR(&BDistIsLessCR, zero, BDist);

                uint32_t BDistIsGreaterCR;
                VECTOR BDistIsGreater = Vector::GreaterR(&BDistIsGreaterCR, BDist, zero);

                // If all the points are on the same side we don't intersect.
                if (ComparisonAllTrue(BDistIsLessCR) || ComparisonAllTrue(BDistIsGreaterCR))
                    return false;

                // Compute the normal of triangle B.
                VECTOR N2 = Vector3::Cross(Vector::Subtract(B1, B0), Vector::Subtract(B2, B0));

#if defined(DEBUG) || defined(_DEBUG)
                // Assert that the triangle is not degenerate.
                assert(!Vector3::Equal(N2, zero));
#endif // DEBUG

                // Test points of A against the plane of B.
                VECTOR ADist = Vector3::Dot(N2, Vector::Subtract(A0, B0));
                ADist = Vector::Select(ADist, Vector3::Dot(N2, Vector::Subtract(A1, B0)), selectY);
                ADist = Vector::Select(ADist, Vector3::Dot(N2, Vector::Subtract(A2, B0)), selectZ);

                // Ensure robustness with co-planar triangles by zeroing small distances.
                uint32_t ADistIsZeroCR;
                VECTOR ADistIsZero = Vector::GreaterR(&ADistIsZeroCR, g_RayEpsilon, Vector::Abs(ADist));
                ADist = Vector::Select(ADist, zero, ADistIsZero);

                uint32_t ADistIsLessCR;
                VECTOR ADistIsLess = Vector::GreaterR(&ADistIsLessCR, zero, ADist);

                uint32_t ADistIsGreaterCR;
                VECTOR ADistIsGreater = Vector::GreaterR(&ADistIsGreaterCR, ADist, zero);

                // If all the points are on the same side we don't intersect.
                if (ComparisonAllTrue(ADistIsLessCR) || ComparisonAllTrue(ADistIsGreaterCR))
                    return false;

                // Special case for co-planar triangles.
                if (ComparisonAllTrue(ADistIsZeroCR) || ComparisonAllTrue(BDistIsZeroCR))
                {
                    VECTOR axis, dist, minDist;

                    // Compute an axis perpindicular to the edge (points out).
                    axis = Vector3::Cross(N1, Vector::Subtract(A1, A0));
                    dist = Vector3::Dot(axis, A0);

                    // Test points of B against the axis.
                    minDist = Vector3::Dot(B0, axis);
                    minDist = Vector::Min(minDist, Vector3::Dot(B1, axis));
                    minDist = Vector::Min(minDist, Vector3::Dot(B2, axis));
                    if (Vector4::GreaterOrEqual(minDist, dist))
                        return false;

                    // Edge (A1, A2)
                    axis = Vector3::Cross(N1, Vector::Subtract(A2, A1));
                    dist = Vector3::Dot(axis, A1);

                    minDist = Vector3::Dot(B0, axis);
                    minDist = Vector::Min(minDist, Vector3::Dot(B1, axis));
                    minDist = Vector::Min(minDist, Vector3::Dot(B2, axis));
                    if (Vector4::GreaterOrEqual(minDist, dist))
                        return false;

                    // Edge (A2, A0)
                    axis = Vector3::Cross(N1, Vector::Subtract(A0, A2));
                    dist = Vector3::Dot(axis, A2);

                    minDist = Vector3::Dot(B0, axis);
                    minDist = Vector::Min(minDist, Vector3::Dot(B1, axis));
                    minDist = Vector::Min(minDist, Vector3::Dot(B2, axis));
                    if (Vector4::GreaterOrEqual(minDist, dist))
                        return false;

                    // Edge (B0, B1)
                    axis = Vector3::Cross(N2, Vector::Subtract(B1, B0));
                    dist = Vector3::Dot(axis, B0);

                    minDist = Vector3::Dot(A0, axis);
                    minDist = Vector::Min(minDist, Vector3::Dot(A1, axis));
                    minDist = Vector::Min(minDist, Vector3::Dot(A2, axis));
                    if (Vector4::GreaterOrEqual(minDist, dist))
                        return false;

                    // Edge (B1, B2)
                    axis = Vector3::Cross(N2, Vector::Subtract(B2, B1));
                    dist = Vector3::Dot(axis, B1);

                    minDist = Vector3::Dot(A0, axis);
                    minDist = Vector::Min(minDist, Vector3::Dot(A1, axis));
                    minDist = Vector::Min(minDist, Vector3::Dot(A2, axis));
                    if (Vector4::GreaterOrEqual(minDist, dist))
                        return false;

                    // Edge (B2,B0)
                    axis = Vector3::Cross(N2, Vector::Subtract(B0, B2));
                    dist = Vector3::Dot(axis, B2);

                    minDist = Vector3::Dot(A0, axis);
                    minDist = Vector::Min(minDist, Vector3::Dot(A1, axis));
                    minDist = Vector::Min(minDist, Vector3::Dot(A2, axis));
                    if (Vector4::GreaterOrEqual(minDist, dist))
                        return false;

                    return true;
                }

                //
                // Find the single vertex of A and B (ie the vertex on the opposite side
                // of the plane from the other two) and reorder the edges so we can compute
                // the signed edge/edge distances.
                //
                // if ( (V0 >= 0 && V1 <  0 && V2 <  0) ||
                //      (V0 >  0 && V1 <= 0 && V2 <= 0) ||
                //      (V0 <= 0 && V1 >  0 && V2 >  0) ||
                //      (V0 <  0 && V1 >= 0 && V2 >= 0) ) then V0 is singular;
                //
                // If our singular vertex is not on the positive side of the plane we reverse
                // the triangle winding so that the overlap comparisons will compare the
                // correct edges with the correct signs.
                //
                VECTOR ADistIsLessEqual = Vector::OrInt(ADistIsLess, ADistIsZero);
                VECTOR ADistIsGreaterEqual = Vector::OrInt(ADistIsGreater, ADistIsZero);

                VECTOR AA0, AA1, AA2;
                bool bPositiveA;

                if (Vector3AllTrue(Vector::Select(ADistIsGreaterEqual, ADistIsLess, select0111)) ||
                    Vector3AllTrue(Vector::Select(ADistIsGreater, ADistIsLessEqual, select0111)))
                {
                    // A0 is singular, crossing from positive to negative.
                    AA0 = A0; AA1 = A1; AA2 = A2;
                    bPositiveA = true;
                }
                else if (Vector3AllTrue(Vector::Select(ADistIsLessEqual, ADistIsGreater, select0111)) ||
                    Vector3AllTrue(Vector::Select(ADistIsLess, ADistIsGreaterEqual, select0111)))
                {
                    // A0 is singular, crossing from negative to positive.
                    AA0 = A0; AA1 = A2; AA2 = A1;
                    bPositiveA = false;
                }
                else if (Vector3AllTrue(Vector::Select(ADistIsGreaterEqual, ADistIsLess, select1011)) ||
                    Vector3AllTrue(Vector::Select(ADistIsGreater, ADistIsLessEqual, select1011)))
                {
                    // A1 is singular, crossing from positive to negative.
                    AA0 = A1; AA1 = A2; AA2 = A0;
                    bPositiveA = true;
                }
                else if (Vector3AllTrue(Vector::Select(ADistIsLessEqual, ADistIsGreater, select1011)) ||
                    Vector3AllTrue(Vector::Select(ADistIsLess, ADistIsGreaterEqual, select1011)))
                {
                    // A1 is singular, crossing from negative to positive.
                    AA0 = A1; AA1 = A0; AA2 = A2;
                    bPositiveA = false;
                }
                else if (Vector3AllTrue(Vector::Select(ADistIsGreaterEqual, ADistIsLess, select1101)) ||
                    Vector3AllTrue(Vector::Select(ADistIsGreater, ADistIsLessEqual, select1101)))
                {
                    // A2 is singular, crossing from positive to negative.
                    AA0 = A2; AA1 = A0; AA2 = A1;
                    bPositiveA = true;
                }
                else if (Vector3AllTrue(Vector::Select(ADistIsLessEqual, ADistIsGreater, select1101)) ||
                    Vector3AllTrue(Vector::Select(ADistIsLess, ADistIsGreaterEqual, select1101)))
                {
                    // A2 is singular, crossing from negative to positive.
                    AA0 = A2; AA1 = A1; AA2 = A0;
                    bPositiveA = false;
                }
                else
                {
#if defined(DEBUG) || defined(_DEBUG)
                    assert(false);
#endif // DEBUG

                    return false;
                }

                VECTOR BDistIsLessEqual = Vector::OrInt(BDistIsLess, BDistIsZero);
                VECTOR BDistIsGreaterEqual = Vector::OrInt(BDistIsGreater, BDistIsZero);

                VECTOR BB0, BB1, BB2;
                bool bPositiveB;

                if (Vector3AllTrue(Vector::Select(BDistIsGreaterEqual, BDistIsLess, select0111)) ||
                    Vector3AllTrue(Select(BDistIsGreater, BDistIsLessEqual, select0111)))
                {
                    // B0 is singular, crossing from positive to negative.
                    BB0 = B0; BB1 = B1; BB2 = B2;
                    bPositiveB = true;
                }
                else if (Vector3AllTrue(Vector::Select(BDistIsLessEqual, BDistIsGreater, select0111)) ||
                    Vector3AllTrue(Vector::Select(BDistIsLess, BDistIsGreaterEqual, select0111)))
                {
                    // B0 is singular, crossing from negative to positive.
                    BB0 = B0; BB1 = B2; BB2 = B1;
                    bPositiveB = false;
                }
                else if (Vector3AllTrue(Vector::Select(BDistIsGreaterEqual, BDistIsLess, select1011)) ||
                    Vector3AllTrue(Vector::Select(BDistIsGreater, BDistIsLessEqual, select1011)))
                {
                    // B1 is singular, crossing from positive to negative.
                    BB0 = B1; BB1 = B2; BB2 = B0;
                    bPositiveB = true;
                }
                else if (Vector3AllTrue(Vector::Select(BDistIsLessEqual, BDistIsGreater, select1011)) ||
                    Vector3AllTrue(Vector::Select(BDistIsLess, BDistIsGreaterEqual, select1011)))
                {
                    // B1 is singular, crossing from negative to positive.
                    BB0 = B1; BB1 = B0; BB2 = B2;
                    bPositiveB = false;
                }
                else if (Vector3AllTrue(Vector::Select(BDistIsGreaterEqual, BDistIsLess, select1101)) ||
                    Vector3AllTrue(Vector::Select(BDistIsGreater, BDistIsLessEqual, select1101)))
                {
                    // B2 is singular, crossing from positive to negative.
                    BB0 = B2; BB1 = B0; BB2 = B1;
                    bPositiveB = true;
                }
                else if (Vector3AllTrue(Vector::Select(BDistIsLessEqual, BDistIsGreater, select1101)) ||
                    Vector3AllTrue(Vector::Select(BDistIsLess, BDistIsGreaterEqual, select1101)))
                {
                    // B2 is singular, crossing from negative to positive.
                    BB0 = B2; BB1 = B1; BB2 = B0;
                    bPositiveB = false;
                }
                else
                {
#if defined(DEBUG) || defined(_DEBUG)
                    assert(false);
#endif // DEBUG

                    return false;
                }

                VECTOR delta0, delta1;

                // Reverse the direction of the test depending on whether the singular vertices are
                // the same sign or different signs.
                if (bPositiveA ^ bPositiveB)
                {
                    delta0 = Vector::Subtract(BB0, AA0);
                    delta1 = Vector::Subtract(AA0, BB0);
                }
                else
                {
                    delta0 = Vector::Subtract(AA0, BB0);
                    delta1 = Vector::Subtract(BB0, AA0);
                }

                // Check if the triangles overlap on the line of intersection between the
                // planes of the two triangles by finding the signed line distances.
                VECTOR dist0 = Vector3::Dot(delta0, Vector3::Cross(Vector::Subtract(BB2, BB0), Vector::Subtract(AA2, AA0)));
                if (Vector4::Greater(dist0, zero))
                    return false;

                VECTOR dist1 = Vector3::Dot(delta1, Vector3::Cross(Vector::Subtract(BB1, BB0), Vector::Subtract(AA1, AA0)));
                if (Vector4::Greater(dist1, zero))
                    return false;

                return true;
            }


            //-----------------------------------------------------------------------------
            // Ray-triangle test
            //-----------------------------------------------------------------------------
            _Use_decl_annotations_
            FORCE_INLINE PlaneIntersectionType VEC_CALLCONV Intersects(A_VECTOR V0, A_VECTOR V1, A_VECTOR V2, B_VECTOR plane) noexcept
            {
                VECTOR one = Vector::SplatOne();

#if defined(DEBUG) || defined(_DEBUG)
                assert(PlaneIsUnit(plane));
#endif // DEBUG

                // Set w of the points to one so we can dot4 with a plane.
                VECTOR TV0 = Vector::Insert<0, 0, 0, 0, 1>(V0, one);
                VECTOR TV1 = Vector::Insert<0, 0, 0, 0, 1>(V1, one);
                VECTOR TV2 = Vector::Insert<0, 0, 0, 0, 1>(V2, one);

                VECTOR outside, inside;
                FastIntersectTrianglePlane(TV0, TV1, TV2, plane, outside, inside);

                // If the triangle is outside any plane it is outside.
                if (Vector4::EqualInt(outside, Vector::TrueInt()))
                    return FRONT;

                // If the triangle is inside all planes it is inside.
                if (Vector4::EqualInt(inside, Vector::TrueInt()))
                    return BACK;

                // The triangle is not inside all planes or outside a plane it intersects.
                return INTERSECTING;
            }


            //-----------------------------------------------------------------------------
            // Test a triangle vs 6 planes (typically forming a frustum).
            //-----------------------------------------------------------------------------
            _Use_decl_annotations_
            FORCE_INLINE ContainmentType VEC_CALLCONV ContainedBy(
                A_VECTOR V0, A_VECTOR V1, A_VECTOR V2,
                B_VECTOR plane0,
                C_VECTOR plane1, C_VECTOR plane2,
                D_VECTOR plane3, D_VECTOR plane4, D_VECTOR plane5) noexcept
            {
                VECTOR one = Vector::SplatOne();

                // Set w of the points to one so we can dot4 with a plane.
                VECTOR TV0 = Vector::Insert<0, 0, 0, 0, 1>(V0, one);
                VECTOR TV1 = Vector::Insert<0, 0, 0, 0, 1>(V1, one);
                VECTOR TV2 = Vector::Insert<0, 0, 0, 0, 1>(V2, one);

                VECTOR outside, inside;

                // Test against each plane.
                FastIntersectTrianglePlane(TV0, TV1, TV2, plane0, outside, inside);

                VECTOR anyOutside = outside;
                VECTOR allInside = inside;

                FastIntersectTrianglePlane(TV0, TV1, TV2, plane1, outside, inside);
                anyOutside = Vector::OrInt(anyOutside, outside);
                allInside = Vector::AndInt(allInside, inside);

                FastIntersectTrianglePlane(TV0, TV1, TV2, plane2, outside, inside);
                anyOutside = Vector::OrInt(anyOutside, outside);
                allInside = Vector::AndInt(allInside, inside);

                FastIntersectTrianglePlane(TV0, TV1, TV2, plane3, outside, inside);
                anyOutside = Vector::OrInt(anyOutside, outside);
                allInside = Vector::AndInt(allInside, inside);

                FastIntersectTrianglePlane(TV0, TV1, TV2, plane4, outside, inside);
                anyOutside = Vector::OrInt(anyOutside, outside);
                allInside = Vector::AndInt(allInside, inside);

                FastIntersectTrianglePlane(TV0, TV1, TV2, plane5, outside, inside);
                anyOutside = Vector::OrInt(anyOutside, outside);
                allInside = Vector::AndInt(allInside, inside);

                // If the triangle is outside any plane it is outside.
                if (Vector4::EqualInt(anyOutside, Vector::TrueInt()))
                    return DISJOINT;

                // If the triangle is inside all planes it is inside.
                if (Vector4::EqualInt(allInside, Vector::TrueInt()))
                    return CONTAINS;

                // The triangle is not inside all planes or outside a plane, it may intersect.
                return INTERSECTS;
            }
        } // namespace TriangleTests
    } // namespace Collision
} // namespace UltReality::Math

#endif // !ULTREALITY_MATH_COLLISION_INL
