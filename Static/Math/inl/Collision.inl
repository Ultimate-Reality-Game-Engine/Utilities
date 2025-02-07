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
                VECTOR C = VEC::Swizzle<SWIZZLE_X, SWIZZLE_Y, SWIZZLE_Z, SWIZZLE_X>(v);

                return ComparisonAnyTrue(VEC4::EqualIntR(C, VEC::TrueInt()));
            }

            //-----------------------------------------------------------------------------
            // Return true if all of the elements of a 3 vector are equal to 0xffffffff.
            // Slightly more efficient than using XMVector3EqualInt.
            //-----------------------------------------------------------------------------
            FORCE_INLINE bool Vector3AllTrue(_In_ A_VECTOR v) noexcept
            {
                // Duplicate the fourth element from the first element.
                VECTOR C = VEC::Swizzle<SWIZZLE_X, SWIZZLE_Y, SWIZZLE_Z, SWIZZLE_X>(v);

                return ComparisonAllTrue(VE4::EqualIntR(C, VEC::TrueInt()));
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
                VECTOR Difference = VEC::Subtract(VE3::Length(v), VEC::SplatOne());
                return VEC4::Less(VEC::Abs(Difference), g_UnitVectorEpsilon);
            }

            //-----------------------------------------------------------------------------
            // Return true if the quaterion is a unit quaternion.
            //-----------------------------------------------------------------------------
            FORCE_INLINE bool QuaternionIsUnit(_In_ A_VECTOR q) noexcept
            {
                VECTOR Difference = VEC::Subtract(VEC4::Length(q), VEC::SplatOne());
                return VEC4::Less(VEC::Abs(Difference), g_UnitQuaternionEpsilon);
            }

            //-----------------------------------------------------------------------------
            // Return true if the plane is a unit plane.
            //-----------------------------------------------------------------------------
            FORCE_INLINE bool PlaneIsUnit(_In_ A_VECTOR plane) noexcept
            {
                VECTOR Difference = VEC::Subtract(VEC3::Length(plane), VEC::SplatOne());
                return VEC4::Less(VEC::Abs(Difference), g_UnitPlaneEpsilon);
            }

#endif // _PREFAST_ || !NDEBUG

            //-----------------------------------------------------------------------------
            FORCE_INLINE VECTOR PlaneTransform(_In_ A_VECTOR plane, _In_ A_VECTOR rotation, _In_ A_VECTOR translation) noexcept
            {
                VECTOR vNormal = VEC3::Rotate(plane, rotation);
                VECTOR vD = VEC::Subtract(VEC::SplatW(plane), VEC3::Dot(vNormal, translation));

                return VEC::Insert<0, 0, 0, 0, 1>(vNormal, vD);
            }

            //-----------------------------------------------------------------------------
            // Return the point on the line segement (S1, S2) nearest the point P.
            //-----------------------------------------------------------------------------
            FORCE_INLINE VECTOR PointOnLineSegmentNearestPoint(_In_ A_VECTOR S1, _In_ A_VECTOR S2, _In_ A_VECTOR p) noexcept
            {
                VECTOR Dir = VEC::Subtract(S2, S1);
                VECTOR Projection = VEC::Subtract(VEC3::Dot(p, Dir), VEC3::Dot(S1, Dir));
                VECTOR LengthSq = VEC3::Dot(Dir, Dir);

                VECTOR t = VEC::Multiply(Projection, VEC::Reciprocal(LengthSq));
                VECTOR Point = VEC::MultiplyAdd(t, Dir, S1);

                // t < 0
                VECTOR SelectS1 = VEC::Less(Projection, VEC::Zero());
                Point = VEC::Select(Point, S1, SelectS1);

                // t > 1
                VECTOR SelectS2 = VEC::Greater(Projection, LengthSq);
                Point = VEC::Select(Point, S2, SelectS2);

                return Point;
            }

            //-----------------------------------------------------------------------------
            // Test if the point (P) on the plane of the triangle is inside the triangle
            // (V0, V1, V2).
            //-----------------------------------------------------------------------------
            FORCE_INLINE VECTOR VEC_CALLCONV PointOnPlaneInsideTriangle(_In_ A_VECTOR p, _In_ A_VECTOR V0, _In_ A_VECTOR V1, _In_ B_VECTOR V2) noexcept
            {
                // Compute the triangle normal.
                VECTOR N = VE3::Cross(VEC::Subtract(V2, V0), VEC::Subtract(V1, V0));

                // Compute the cross products of the vector from the base of each edge to
                // the point with each edge vector.
                VECTOR C0 = VE3::Cross(VEC::Subtract(p, V0), VEC::Subtract(V1, V0));
                VECTOR C1 = VEC3::Cross(VEC::Subtract(p, V1), VEC::Subtract(V2, V1));
                VECTOR C2 = VEC3::Cross(VEC::Subtract(p, V2), VEC::Subtract(V0, V2));

                // If the cross product points in the same direction as the normal the the
                // point is inside the edge (it is zero if is on the edge).
                VECTOR zero = VEC::Zero();
                VECTOR Inside0 = VEC::GreaterOrEqual(VEC3::Dot(C0, N), zero);
                VECTOR Inside1 = VEC::GreaterOrEqual(VEC3::Dot(C1, N), zero);
                VECTOR Inside2 = VEC::GreaterOrEqual(VE3::Dot(C2, N), zero);

                // If the point inside all of the edges it is inside.
                return VEC::AndInt(VEC::AndInt(Inside0, Inside1), Inside2);
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
                sinth3 = sqrtf(3.0f) * ScalarSin(theta / 3.0f);
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

                VECTOR vTmp = VEC::LoadFloat3(reinterpret_cast<const Float3*>(fTmp));

                if (VEC3::Equal(vTmp, VEC::Zero())) // planar or linear
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
                        vTmp = VEC::SetX(vTmp, 0.0f);
                    else
                        vTmp = VEC::SetX(vTmp, 1.0f);

                    if (f2 == 0)
                        vTmp = VEC::SetY(vTmp, 0.0f);
                    else
                        vTmp = VEC::SetY(vTmp, 1.0f);

                    if (f3 == 0)
                    {
                        vTmp = VEC::SetZ(vTmp, 0.0f);
                        // recalculate y to make equation work
                        if (m12 != 0)
                            vTmp = VEC::SetY(vTmp, -f1 / f2);
                    }
                    else
                    {
                        vTmp = VEC::SetZ(vTmp, (f2 - f1) / f3);
                    }
                }

                if (VEC::GetX(VEC3::LengthSq(vTmp)) > 1e-5f)
                {
                    return VEC3::Normalize(vTmp);
                }
                else
                {
                    // Multiply by a value large enough to make the vector non-zero.
                    vTmp = VEC::Scale(vTmp, 1e5f);
                    return VEC3::Normalize(vTmp);
                }
            }

            //-----------------------------------------------------------------------------
            FORCE_INLINE bool CalculateEigenVectors(_In_ float m11, _In_ float m12, _In_ float m13,
                _In_ float m22, _In_ float m23, _In_ float m33,
                _In_ float e1, _In_ float e2, _In_ float e3,
                _Out_ XMVECTOR* pV1, _Out_ XMVECTOR* pV2, _Out_ XMVECTOR* pV3) noexcept
            {
                *pV1 = CalculateEigenVector(m11, m12, m13, m22, m23, m33, e1);
                *pV2 = CalculateEigenVector(m11, m12, m13, m22, m23, m33, e2);
                *pV3 = CalculateEigenVector(m11, m12, m13, m22, m23, m33, e3);

                bool v1z = false;
                bool v2z = false;
                bool v3z = false;

                VECTOR zero = VEC::Zero();

                if (VEC3::Equal(*pV1, zero))
                    v1z = true;

                if (VEC3::Equal(*pV2, zero))
                    v2z = true;

                if (VEC3::Equal(*pV3, zero))
                    v3z = true;

                bool e12 = (fabsf(VEC::GetX(VEC3::Dot(*pV1, *pV2))) > 0.1f); // check for non-orthogonal vectors
                bool e13 = (fabsf(VEC::GetX(VEC3::Dot(*pV1, *pV3))) > 0.1f);
                bool e23 = (fabsf(VEC::GetX(VEC3::Dot(*pV2, *pV3))) > 0.1f);

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
                    VECTOR vTmp = VEC3::Cross(g_IdentityR1, *pV3);
                    if (VEC::GetX(VEC3::LengthSq(vTmp)) < 1e-5f)
                    {
                        vTmp = VEC3::Cross(g_IdentityR0, *pV3);
                    }
                    *pV1 = VEC3::Normalize(vTmp);
                    *pV2 = VEC3::Cross(*pV3, *pV1);
                    return true;
                }

                if (v3z && v1z)
                {
                    VECTOR vTmp = VEC3::Cross(g_IdentityR1, *pV2);
                    if (VEC::GetX(VEC3::LengthSq(vTmp)) < 1e-5f)
                    {
                        vTmp = VEC3::Cross(g_IdentityR0, *pV2);
                    }
                    *pV3 = VEC3::Normalize(vTmp);
                    *pV1 = VEC3::Cross(*pV2, *pV3);
                    return true;
                }

                if (v2z && v3z)
                {
                    VECTOR vTmp = VEC3::Cross(g_IdentityR1, *pV1);
                    if (VEC::GetX(VEC3::LengthSq(vTmp)) < 1e-5f)
                    {
                        vTmp = VEC3::Cross(g_IdentityR0, *pV1);
                    }
                    *pV2 = VEC3::Normalize(vTmp);
                    *pV3 = VEC3::Cross(*pV1, *pV2);
                    return true;
                }

                if ((v1z) || e12)
                {
                    *pV1 = VEC3::Cross(*pV2, *pV3);
                    return true;
                }

                if ((v2z) || e23)
                {
                    *pV2 = VEC3::Cross(*pV3, *pV1);
                    return true;
                }

                if ((v3z) || e13)
                {
                    *pV3 = VEC3::Cross(*pV1, *pV2);
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
                VECTOR Dist0 = VEC4::Dot(V0, plane);
                VECTOR Dist1 = VEC4::Dot(V1, plane);
                VECTOR Dist2 = VEC4::Dot(V2, plane);

                VECTOR MinDist = VEC::Min(Dist0, Dist1);
                MinDist = VEC::Min(MinDist, Dist2);

                VECTOR MaxDist = VEC::Max(Dist0, Dist1);
                MaxDist = VEC::Max(MaxDist, Dist2);

                VECTOR zero = VEC::Zero();

                // Outside the plane?
                Outside = VEC::Greater(MinDist, zero);

                // Fully inside the plane?
                Inside = VEC::Less(MaxDist, zero);
            }

            //-----------------------------------------------------------------------------
            FORCE_INLINE void FastIntersectSpherePlane(
                _In_ A_VECTOR center, _In_ A_VECTOR radius, _In_ A_VECTOR plane,
                _Out_ VECTOR& outside, _Out_ VECTOR& inside) noexcept
            {
                VECTOR Dist = VEC4::Dot(center, plane);

                // Outside the plane?
                outside = VEC::Greater(Dist, radius);

                // Fully inside the plane?
                inside = VEC::Less(Dist, VEC::Negate(radius));
            }

            //-----------------------------------------------------------------------------
            FORCE_INLINE void FastIntersectAxisAlignedBoxPlane(
                _In_ A_VECTOR center, _In_ A_VECTOR extents, _In_ A_VECTOR plane,
                _Out_ VECTOR& outside, _Out_ VECTOR& inside) noexcept
            {
                // Compute the distance to the center of the box.
                VECTOR Dist = VEC4::Dot(center, plane);

                // Project the axes of the box onto the normal of the plane.  Half the
                // length of the projection (sometime called the "radius") is equal to
                // h(u) * abs(n dot b(u))) + h(v) * abs(n dot b(v)) + h(w) * abs(n dot b(w))
                // where h(i) are extents of the box, n is the plane normal, and b(i) are the
                // axes of the box. In this case b(i) = [(1,0,0), (0,1,0), (0,0,1)].
                VECTOR Radius = VEC3::Dot(extents, VEC::Abs(plane));

                // Outside the plane?
                outside = VEC::Greater(Dist, Radius);

                // Fully inside the plane?
                inside = VEC::Less(Dist, VEC::Negate(Radius));
            }

            //-----------------------------------------------------------------------------
            FORCE_INLINE void VEC_CALLCONV FastIntersectOrientedBoxPlane(
                _In_ A_VECTOR center, _In_ A_VECTOR extents, _In_ A_VECTOR axis0, 
                _In_ B_VECTOR axis1, 
                _In_ C_VECTOR axis2, _In_ C_VECTOR plane,
                _Out_ VECTOR& outside, _Out_ VECTOR& inside) noexcept
            {
                // Compute the distance to the center of the box.
                VECTOR Dist = VEC4::Dot(center, plane);

                // Project the axes of the box onto the normal of the plane.  Half the
                // length of the projection (sometime called the "radius") is equal to
                // h(u) * abs(n dot b(u))) + h(v) * abs(n dot b(v)) + h(w) * abs(n dot b(w))
                // where h(i) are extents of the box, n is the plane normal, and b(i) are the
                // axes of the box.
                VECTOR Radius = VEC3::Dot(plane, axis0);
                Radius = VEC::Insert<0, 0, 1, 0, 0>(Radius, VEC3::Dot(plane, axis1));
                Radius = VEC::Insert<0, 0, 0, 1, 0>(Radius, VEC3::Dot(plane, axis2));
                Radius = VEC3::Dot(extents, VEC::Abs(Radius));

                // Outside the plane?
                outside = VEC::Greater(Dist, Radius);

                // Fully inside the plane?
                inside = VEC::Less(Dist, VEC::Negate(Radius));
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

                Min = Max = VEC3::Dot(plane, point0);

                Dist = VEC3::Dot(plane, point1);
                Min = VEC::Min(Min, Dist);
                Max = VEC::Max(Max, Dist);

                Dist = VEC3::Dot(plane, point2);
                Min = VEC::Min(Min, Dist);
                Max = VEC::Max(Max, Dist);

                Dist = VEC3::Dot(plane, point3);
                Min = VEC::Min(Min, Dist);
                Max = VEC::Max(Max, Dist);

                Dist = VEC3::Dot(plane, point4);
                Min = VEC::Min(Min, Dist);
                Max = VEC::Max(Max, Dist);

                Dist = VEC3::Dot(plane, point5);
                Min = VEC::Min(Min, Dist);
                Max = VEC::Max(Max, Dist);

                Dist = VEC3::Dot(plane, point6);
                Min = VEC::Min(Min, Dist);
                Max = VEC::Max(Max, Dist);

                Dist = VEC3::Dot(plane, point7);
                Min = VEC::Min(Min, Dist);
                Max = VEC::Max(Max, Dist);

                VECTOR PlaneDist = VEC::Negate(VEC::SplatW(plane));

                // Outside the plane?
                outside = VEC::Greater(Min, PlaneDist);

                // Fully inside the plane?
                inside = VEC::Less(Max, PlaneDist);
            }
        }
        
        _Use_decl_annotations_
        FORCE_INLINE void VEC_CALLCONV BoundingSphere::Transform(BoundingSphere& out, A_MATRIX m) const noexcept
        {
            // Load the center of the sphere.
            VECTOR vCenter = VEC::LoadFloat3(&center);

            // Transform the center of the sphere.
            VECTOR C = VEC3::Transform(vCenter, m);

            VECTOR dX = VEC3::Dot(m.r[0], m.r[0]);
            VECTOR dY = VEC3::Dot(m.r[1], m.r[1]);
            VECTOR dZ = VEC3::Dot(m.r[2], m.r[2]);

            VECTOR d = VEC::Max(dX, VEC::Max(dY, dZ));

            // Store the center sphere.
            VEC::StoreFloat3(&out.center, C);

            // Scale the radius of the pshere.
            float scale = sqrtf(VEC::GetX(d));
            out.radius = radius * scale;
        }

        _Use_decl_annotations_
        FORCE_INLINE void VEC_CALLCONV BoundingSphere::Transform(BoundingSphere& out, float scale, A_VECTOR rotation, A_VECTOR translation) const noexcept
        {
            // Load the center of the sphere.
            VECTOR vCenter = VEC::LoadFloat3(&Ccnter);

            // Transform the center of the sphere.
            vCenter = VEC::Add(VEC3::Rotate(VEC::Scale(vCenter, scale), rotation), translation);

            // Store the center sphere.
            VEC::StoreFloat3(&out.center, vCenter);

            // Scale the radius of the pshere.
            out.radius = radius * scale;
        }

        _Use_decl_annotations_
        FORCE_INLINE ContainmentType VEC_CALLCONV BoundingSphere::Contains(A_VECTOR point) const noexcept
        {
            VECTOR vCenter = VEC::LoadFloat3(&center);
            VECTOR vRadius = VEC::ReplicatePtr(&radius);

            VECTOR DistanceSquared = VEC3::LengthSq(VEC::Subtract(point, vCenter));
            VECTOR RadiusSquared = VEC::Multiply(vRadius, vRadius);

            return VEC3::LessOrEqual(DistanceSquared, RadiusSquared) ? CONTAINS : DISJOINT;
        }

        _Use_decl_annotations_
        FORCE_INLINE ContainmentType VEC_CALLCONV BoundingSphere::Contains(A_VECTOR V0, A_VECTOR V1, A_VECTOR V2) const noexcept
        {
            if (!Intersects(V0, V1, V2))
                return DISJOINT;

            VECTOR vCenter = VEC::LoadFloat3(&center);
            VECTOR vRadius = VEC::ReplicatePtr(&radius);
            VECTOR RadiusSquared = VEC::Multiply(vRadius, vRadius);

            VECTOR DistanceSquared = VEC3::LengthSq(VEC::Subtract(V0, vCenter));
            VECTOR Inside = VEC::LessOrEqual(DistanceSquared, RadiusSquared);

            DistanceSquared = VEC3::LengthSq(VEC::Subtract(V1, vCenter));
            Inside = VEC::AndInt(Inside, VEC::LessOrEqual(DistanceSquared, RadiusSquared));

            DistanceSquared = VEC3::LengthSq(VEC::Subtract(V2, vCenter));
            Inside = VEC::AndInt(Inside, VEC::LessOrEqual(DistanceSquared, RadiusSquared));

            return (VEC3::EqualInt(Inside, VEC::TrueInt())) ? CONTAINS : INTERSECTS;
        }

        _Use_decl_annotations_
        FORCE_INLINE ContainmentType BoundingSphere::Contains(const BoundingSphere& sh) const noexcept
        {
            VECTOR Center1 = VEC::LoadFloat3(&center);
            float r1 = radius;

            VECTOR Center2 = VEC::LoadFloat3(&sh.center);
            float r2 = sh.radius;

            VECTOR V = VEC::Subtract(Center2, Center1);

            VECTOR Dist = VEC3::Length(V);

            float d = VEC::GetX(Dist);

            return (r1 + r2 >= d) ? ((r1 - r2 >= d) ? CONTAINS : INTERSECTS) : DISJOINT;
        }

        _Use_decl_annotations_
        FORCE_INLINE ContainmentType BoundingSphere::Contains(const BoundingBox& box) const noexcept
        {
            if (!box.Intersects(*this))
                return DISJOINT;

            VECTOR vCenter = VEC::LoadFloat3(&center);
            VECTOR vRadius = VEC::ReplicatePtr(&radius);
            VECTOR RadiusSq = VEC::Multiply(vRadius, vRadius);

            VECTOR boxCenter = VEC::LoadFloat3(&box.center);
            VECTOR boxExtents = VEC::LoadFloat3(&box.extents);

            VECTOR InsideAll = VEC::TrueInt();

            VECTOR offset = VEC::Subtract(boxCenter, vCenter);

            for (size_t i = 0; i < BoundingBox::CORNER_COUNT; ++i)
            {
                VECTOR C = VEC::MultiplyAdd(boxExtents, g_BoxOffset[i], offset);
                VECTOR d = VEC3::LengthSq(C);
                InsideAll = VEC::AndInt(InsideAll, VEC::LessOrEqual(d, RadiusSq));
            }

            return (VEC3::EqualInt(InsideAll, VEC::TrueInt())) ? CONTAINS : INTERSECTS;
        }

        _Use_decl_annotations_
        FORCE_INLINE ContainmentType BoundingSphere::Contains(const BoundingOrientedBox& box) const noexcept
        {
            if (!box.Intersects(*this))
                return DISJOINT;

            VECTOR vCenter = VEC::LoadFloat3(&center);
            VECTOR vRadius = VEC::ReplicatePtr(&radius);
            VECTOR RadiusSq = VEC::Multiply(vRadius, vRadius);

            VECTOR boxCenter = VEC::LoadFloat3(&box.center);
            VECTOR boxExtents = VEC::LoadFloat3(&box.extents);
            VECTOR boxOrientation = VEC::LoadFloat4(&box.orientation);

#if defined(DEBUG) || defined(_DEBUG)
            assert(QuaternionIsUnit(boxOrientation));
#endif // DEBUG

            VECTOR InsideAll = VEC::TrueInt();

            for (size_t i = 0; i < BoundingOrientedBox::CORNER_COUNT; ++i)
            {
                VECTOR C = VEC::Add(VEC3::Rotate(VEC::Multiply(boxExtents, g_BoxOffset[i]), boxOrientation), boxCenter);
                VECTOR d = VEC3::LengthSq(VEC::Subtract(vCenter, C));
                InsideAll = VEC::AndInt(InsideAll, VEC::LessOrEqual(d, RadiusSq));
            }

            return (VEC3::EqualInt(InsideAll, VEC::TrueInt())) ? CONTAINS : INTERSECTS;
        }

        _Use_decl_annotations_
        FORCE_INLINE ContainmentType BoundingSphere::Contains(const BoundingFrustum& fr) const noexcept
        {
            if (!fr.Intersects(*this))
                return DISJOINT;

            VECTOR vCenter = VEC::LoadFloat3(&center);
            VECTOR vRadius = VEC::ReplicatePtr(&radius);
            VECTOR RadiusSq = VEC::Multiply(vRadius, vRadius);

            VECTOR vOrigin = VEC::LoadFloat3(&fr.origin);
            VECTOR vOrientation = VEC::LoadFloat4(&fr.orientation);

#if defined(DEBUG) || defined(_DEBUG)
            assert(QuaternionIsUnit(vOrientation));
#endif // DEBUG

            // Build the corners of the frustum.
            VECTOR vRightTop = VEC::Set(fr.rightSlope, fr.topSlope, 1.0f, 0.0f);
            VECTOR vRightBottom = VEC::Set(fr.rightSlope, fr.bottomSlope, 1.0f, 0.0f);
            VECTOR vLeftTop = VEC::Set(fr.leftSlope, fr.topSlope, 1.0f, 0.0f);
            VECTOR vLeftBottom = VEC::Set(fr.leftSlope, fr.bottomSlope, 1.0f, 0.0f);
            VECTOR vNear = VEC::ReplicatePtr(&fr.near);
            VECTOR vFar = VEC::ReplicatePtr(&fr.far);

            VECTOR Corners[BoundingFrustum::CORNER_COUNT];
            Corners[0] = VEC::Multiply(vRightTop, vNear);
            Corners[1] = VEC::Multiply(vRightBottom, vNear);
            Corners[2] = VEC::Multiply(vLeftTop, vNear);
            Corners[3] = VEC::Multiply(vLeftBottom, vNear);
            Corners[4] = VEC::Multiply(vRightTop, vFar);
            Corners[5] = VEC::Multiply(vRightBottom, vFar);
            Corners[6] = VEC::Multiply(vLeftTop, vFar);
            Corners[7] = VEC::Multiply(vLeftBottom, vFar);

            VECTOR InsideAll = VEC::TrueInt();
            for (size_t i = 0; i < BoundingFrustum::CORNER_COUNT; ++i)
            {
                VECTOR C = VEC::Add(VEC3::Rotate(Corners[i], vOrientation), vOrigin);
                VECTOR d = VEC3::LengthSq(VEC::Subtract(vCenter, C));
                InsideAll = VEC::AndInt(InsideAll, VEC::LessOrEqual(d, RadiusSq));
            }

            return (VEC3::EqualInt(InsideAll, VEC::TrueInt())) ? CONTAINS : INTERSECTS;
        }

        _Use_decl_annotations_
        FORCE_INLINE bool BoundingSphere::Intersects(const BoundingSphere& sh) const noexcept
        {
            // Load A.
            VECTOR vCenterA = VEC::LoadFloat3(&center);
            VECTOR vRadiusA = VEC::ReplicatePtr(&radius);

            // Load B.
            VECTOR vCenterB = VEC::LoadFloat3(&sh.center);
            VECTOR vRadiusB = VEC::ReplicatePtr(&sh.radius);

            // Distance squared between centers.
            VECTOR Delta = VEC::Subtract(vCenterB, vCenterA);
            VECTOR DistanceSquared = VEC3::LengthSq(Delta);

            // Sum of the radii squared.
            VECTOR RadiusSquared = VEC::Add(vRadiusA, vRadiusB);
            RadiusSquared = VEC::Multiply(RadiusSquared, RadiusSquared);

            return VEC3::LessOrEqual(DistanceSquared, RadiusSquared);
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
            VECTOR vCenter = VEC::LoadFloat3(&center);
            VECTOR vRadius = VEC::ReplicatePtr(&radius);

            // Compute the plane of the triangle (has to be normalized).
            VECTOR N = VEC3::Normalize(VEC3::Cross(VEC::Subtract(V1, V0), VEC::Subtract(V2, V0)));

#if defined(DEBUG) || defined(_DEBUG)
            // Assert that the triangle is not degenerate.
            assert(!VEC3::Equal(N, VEC::Zero()));
#endif // DEBUG

            // Find the nearest feature on the triangle to the sphere.
            VECTOR Dist = VEC3::Dot(VEC::Subtract(vCenter, V0), N);

            // If the center of the sphere is farther from the plane of the triangle than
            // the radius of the sphere, then there cannot be an intersection.
            VECTOR NoIntersection = VEC::Less(Dist, VEC::Negate(vRadius));
            NoIntersection = VEC::OrInt(NoIntersection, VEC::Greater(Dist, vRadius));

            // Project the center of the sphere onto the plane of the triangle.
            VECTOR Point = VEC::NegativeMultiplySubtract(N, Dist, vCenter);

            // Is it inside all the edges? If so we intersect because the distance
            // to the plane is less than the radius.
            VECTOR Intersection = PointOnPlaneInsideTriangle(Point, V0, V1, V2);

            // Find the nearest point on each edge.
            VECTOR RadiusSq = VEC::Multiply(vRadius, vRadius);

            // Edge 0,1
            Point = PointOnLineSegmentNearestPoint(V0, V1, vCenter);

            // If the distance to the center of the sphere to the point is less than
            // the radius of the sphere then it must intersect.
            Intersection = VEC::OrInt(Intersection, VEC::LessOrEqual(VEC3::LengthSq(VEC::Subtract(vCenter, Point)), RadiusSq));

            // Edge 1,2
            Point = PointOnLineSegmentNearestPoint(V1, V2, vCenter);

            // If the distance to the center of the sphere to the point is less than
            // the radius of the sphere then it must intersect.
            Intersection = VEC::OrInt(Intersection, VEC::LessOrEqual(VEC3::LengthSq(VEC::Subtract(vCenter, Point)), RadiusSq));

            // Edge 2,0
            Point = PointOnLineSegmentNearestPoint(V2, V0, vCenter);

            // If the distance to the center of the sphere to the point is less than
            // the radius of the sphere then it must intersect.
            Intersection = VEC::OrInt(Intersection, VEC::LessOrEqual(VEC3::LengthSq(VEC::Subtract(vCenter, Point)), RadiusSq));

            return VEC4::EqualInt(VEC::AndCInt(Intersection, NoIntersection), VEC::TrueInt());
        }

        _Use_decl_annotations_
        FORCE_INLINE PlaneIntersectionType VEC_CALLCONV BoundingSphere::Intersects(A_VECTOR plane) const noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(PlaneIsUnit(plane));
#endif // DEBUG

            // Load the sphere.
            VECTOR vCenter = VEC::LoadFloat3(&center);
            VECTOR vRadius = VEC::ReplicatePtr(&radius);

            // Set w of the center to one so we can dot4 with a plane.
            vCenter = VEC::Insert<0, 0, 0, 0, 1>(vCenter, VEC::SplatOne());

            VECTOR outside, inside;
            FastIntersectSpherePlane(vCenter, vRadius, plane, outside, inside);

            // If the sphere is outside any plane it is outside.
            if (VEC4::EqualInt(outside, VEC::TrueInt()))
                return FRONT;

            // If the sphere is inside all planes it is inside.
            if (VEC4::EqualInt(inside, VEC::TrueInt()))
                return BACK;

            // The sphere is not inside all planes or outside a plane it intersects.
            return INTERSECTING;
        }

        _Use_decl_annotations_
        FORCE_INLINE bool VEC_CALLCONV BoundingSphere::Intersects(A_VECTOR origin, A_VECTOR direction, float& dist) const noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(VEC3::IsUnit(direction));
#endif // DEBUG

            VECTOR vCenter = VEC::LoadFloat3(&center);
            VECTOR vRadius = VEC::ReplicatePtr(&radius);

            // l is the vector from the ray origin to the center of the sphere.
            VECTOR l = VEC::Subtract(vCenter, origin);

            // s is the projection of the l onto the ray direction.
            VECTOR s = VEC3::Dot(l, direction);

            VECTOR l2 = VEC3::Dot(l, l);

            VECTOR r2 = VEC::Multiply(vRadius, vRadius);

            // m2 is squared distance from the center of the sphere to the projection.
            VECTOR m2 = VEC::NegativeMultiplySubtract(s, s, l2);

            VECTOR NoIntersection;

            // If the ray origin is outside the sphere and the center of the sphere is
            // behind the ray origin there is no intersection.
            NoIntersection = VEC::AndInt(VEC::Less(s, VEC::Zero()), VEC::Greater(l2, r2));

            // If the squared distance from the center of the sphere to the projection
            // is greater than the radius squared the ray will miss the sphere.
            NoIntersection = VEC::OrInt(NoIntersection, VEC::Greater(m2, r2));

            // The ray hits the sphere, compute the nearest intersection point.
            VECTOR q = VEC::Sqrt(VEC::Subtract(r2, m2));
            VECTOR t1 = VEC::Subtract(s, q);
            VECTOR t2 = VEC::Add(s, q);

            VECTOR OriginInside = VEC::LessOrEqual(l2, r2);
            VECTOR t = VEC::Select(t1, t2, OriginInside);

            if (VEC4::NotEqualInt(NoIntersection, VEC::TrueInt()))
            {
                // Store the x-component to *pDist.
                VEC::StoreFloat(&dist, t);
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
            VECTOR vCenter = VEC::LoadFloat3(&center);
            VECTOR vRadius = VEC::ReplicatePtr(&radius);

            // Set w of the center to one so we can dot4 with a plane.
            vCenter = VEC::Insert<0, 0, 0, 0, 1>(vCenter, VEC::SplatOne());

            VECTOR outside, inside;

            // Test against each plane.
            FastIntersectSpherePlane(vCenter, vRadius, plane0, outside, inside);

            VECTOR AnyOutside = outside;
            VECTOR AllInside = inside;

            FastIntersectSpherePlane(vCenter, vRadius, plane1, outside, inside);
            AnyOutside = VEC::OrInt(AnyOutside, outside);
            AllInside = VEC::AndInt(AllInside, inside);

            FastIntersectSpherePlane(vCenter, vRadius, plane2, outside, inside);
            AnyOutside = VEC::OrInt(AnyOutside, outside);
            AllInside = VEC::AndInt(AllInside, inside);

            FastIntersectSpherePlane(vCenter, vRadius, plane3, outside, inside);
            AnyOutside = VEC::OrInt(AnyOutside, outside);
            AllInside = VEC::AndInt(AllInside, inside);

            FastIntersectSpherePlane(vCenter, vRadius, plane4, outside, inside);
            AnyOutside = VEC::OrInt(AnyOutside, outside);
            AllInside = VEC::AndInt(AllInside, inside);

            FastIntersectSpherePlane(vCenter, vRadius, plane5, outside, inside);
            AnyOutside = VEC::OrInt(AnyOutside, outside);
            AllInside = VEC::AndInt(AllInside, inside);

            // If the sphere is outside any plane it is outside.
            if (VEC4::EqualInt(AnyOutside, VEC::TrueInt()))
                return DISJOINT;

            // If the sphere is inside all planes it is inside.
            if (VEC4::EqualInt(AllInside, VEC::TrueInt()))
                return CONTAINS;

            // The sphere is not inside all planes or outside a plane, it may intersect.
            return INTERSECTS;
        }

        _Use_decl_annotations_
        FORCE_INLINE void BoundingSphere::CreateMerged(BoundingSphere& out, const BoundingSphere& S1, const BoundingSphere& S2) noexcept
        {
            VECTOR Center1 = VEC::LoadFloat3(&S1.center);
            float r1 = S1.radius;

            VECTOR Center2 = VEC::LoadFloat3(&S2.center);
            float r2 = S2.radius;

            VECTOR V = VEC::Subtract(Center2, Center1);

            VECTOR Dist = VEC3::Length(V);

            float d = VEC::GetX(Dist);

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

            VECTOR N = VEC::Divide(V, Dist);

            float t1 = Min(-r1, d - r2);
            float t2 = Max(r1, d + r2);
            float t_5 = (t2 - t1) * 0.5f;

            VECTOR NCenter = VEC::Add(Center1, VEC::Multiply(N, VEC::Replicate(t_5 + t1)));

            VEC::StoreFloat3(&out.center, NCenter);
            out.radius = t_5;
        }

        _Use_decl_annotations_
        FORCE_INLINE void BoundingSphere::CreateFromBoundingBox(BoundingSphere& out, const BoundingBox& box) noexcept
        {
            out.center = box.center;
            VECTOR vExtents = VEC::LoadFloat3(&box.extents);
            out.radius = VEC::GetX(VEC3::Length(vExtents));
        }

        _Use_decl_annotations_
        FORCE_INLINE void BoundingSphere::CreateFromBoundingBox(BoundingSphere& out, const BoundingOrientedBox& box) noexcept
        {
            // Bounding box orientation is irrelevant because a sphere is rotationally invariant
            out.center = box.center;
            VECTOR vExtents = VEC::LoadFloat3(&box.extents);
            out.radius = VEC::GetX(VEC3::Length(vExtents));
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

            MinX = MaxX = MinY = MaxY = MinZ = MaxZ = VEC::LoadFloat3(pPoints);

            for (size_t i = 1; i < count; ++i)
            {
                VECTOR Point = VEC::LoadFloat3(reinterpret_cast<const Float3*>(reinterpret_cast<const uint8_t*>(pPoints) + i * stride));

                float px = VEC::GetX(Point);
                float py = VEC::GetY(Point);
                float pz = VEC::GetZ(Point);

                if (px < VEC::GetX(MinX))
                    MinX = Point;

                if (px > VEC::GetX(MaxX))
                    MaxX = Point;

                if (py < VEC::GetY(MinY))
                    MinY = Point;

                if (py > VEC::GetY(MaxY))
                    MaxY = Point;

                if (pz < VEC::GetZ(MinZ))
                    MinZ = Point;

                if (pz > VEC::GetZ(MaxZ))
                    MaxZ = Point;
            }

            // Use the min/max pair that are farthest apart to form the initial sphere.
            VECTOR DeltaX = VEC::Subtract(MaxX, MinX);
            VECTOR DistX = VEC3::Length(DeltaX);

            VECTOR DeltaY = VEC::Subtract(MaxY, MinY);
            VECTOR DistY = VEC3::Length(DeltaY);

            VECTOR DeltaZ = VEC::Subtract(MaxZ, MinZ);
            VECTOR DistZ = VEC3::Length(DeltaZ);

            VECTOR vCenter;
            VECTOR vRadius;

            if (VEC3::Greater(DistX, DistY))
            {
                if (VEC3::Greater(DistX, DistZ))
                {
                    // Use min/max x.
                    vCenter = VEC::Lerp(MaxX, MinX, 0.5f);
                    vRadius = VEC::Scale(DistX, 0.5f);
                }
                else
                {
                    // Use min/max z.
                    vCenter = VEC::Lerp(MaxZ, MinZ, 0.5f);
                    vRadius = VEC::Scale(DistZ, 0.5f);
                }
            }
            else // Y >= X
            {
                if (VEC3::Greater(DistY, DistZ))
                {
                    // Use min/max y.
                    vCenter = VEC::Lerp(MaxY, MinY, 0.5f);
                    vRadius = VEC::Scale(DistY, 0.5f);
                }
                else
                {
                    // Use min/max z.
                    vCenter = VEC::Lerp(MaxZ, MinZ, 0.5f);
                    vRadius = VEC::Scale(DistZ, 0.5f);
                }
            }

            // Add any points not inside the sphere.
            for (size_t i = 0; i < count; ++i)
            {
                VECTOR Point = VEC::LoadFloat3(reinterpret_cast<const Float3*>(reinterpret_cast<const uint8_t*>(pPoints) + i * stride));

                VECTOR Delta = VEC::Subtract(Point, vCenter);

                VECTOR Dist = VEC3::Length(Delta);

                if (VEC3::Greater(Dist, vRadius))
                {
                    // Adjust sphere to include the new point.
                    vRadius = VEC::Scale(VEC::Add(vRadius, Dist), 0.5f);
                    vCenter = VEC::Add(vCenter, VEC::Multiply(VEC::Subtract(VEC::Replicate(1.0f), VEC::Divide(vRadius, Dist)), Delta));
                }
            }

            VEC::StoreFloat3(&out.center, vCenter);
            VEC::StoreFloat(&out.radius, vRadius);
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
            VECTOR vCenter = VEC::LoadFloat3(&center);
            VECTOR vExtents = VEC::LoadFloat3(&extents);

            // Compute and transform the corners and find new min/max bounds.
            VECTOR Corner = VEC::MultiplyAdd(vExtents, g_BoxOffset[0], vCenter);
            Corner = VEC3::Transform(Corner, M);

            VECTOR Min, Max;
            Min = Max = Corner;

            for (size_t i = 1; i < CORNER_COUNT; ++i)
            {
                Corner = VEC::MultiplyAdd(vExtents, g_BoxOffset[i], vCenter);
                Corner = VEC3::Transform(Corner, M);

                Min = VEC::Min(Min, Corner);
                Max = VEC::Max(Max, Corner);
            }

            // Store center and extents.
            VEC::StoreFloat3(&out.center, VEC::Scale(VEC::Add(Min, Max), 0.5f));
            VEC::StoreFloat3(&out.extents, VEC::Scale(VEC::Subtract(Max, Min), 0.5f));
        }

        _Use_decl_annotations_
        FORCE_INLINE void VEC_CALLCONV BoundingBox::Transform(BoundingBox& out, float scale, A_VECTOR rotation, A_VECTOR translation) const noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(QuaternionIsUnit(rotation));
#endif // DEBUG

            // Load center and extents.
            VECTOR vCenter = VEC::LoadFloat3(&Ccnter);
            VECTOR vExtents = VEC::LoadFloat3(&extents);

            VECTOR VectorScale = VEC::Replicate(scale);

            // Compute and transform the corners and find new min/max bounds.
            VECTOR Corner = VEC::MultiplyAdd(vExtents, g_BoxOffset[0], vCenter);
            Corner = VEC::Add(VEC3::Rotate(VEC::Multiply(Corner, VectorScale), rotation), translation);

            VECTOR Min, Max;
            Min = Max = Corner;

            for (size_t i = 1; i < CORNER_COUNT; ++i)
            {
                Corner = VEC::MultiplyAdd(vExtents, g_BoxOffset[i], vCenter);
                Corner = VEC::Add(VEC3::Rotate(VEC::Multiply(Corner, VectorScale), rotation), translation);

                Min = VEC::Min(Min, Corner);
                Max = VEC::Max(Max, Corner);
            }

            // Store center and extents.
            VEC::StoreFloat3(&out.center, VEC::Scale(VEC::Add(Min, Max), 0.5f));
            VEC::StoreFloat3(&out.extents, VEC::Scale(VEC::Subtract(Max, Min), 0.5f));
        }

        _Use_decl_annotations_
        FORCE_INLINE void BoundingBox::GetCorners(Float3* corners) const noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(corners != nullptr);
#endif // DEBUG

            // Load the box
            VECTOR vCenter = VEC::LoadFloat3(&center);
            VECTOR vExtents = VEC::LoadFloat3(&extents);

            for (size_t i = 0; i < CORNER_COUNT; ++i)
            {
                VECTOR C = VEC::MultiplyAdd(vExtents, g_BoxOffset[i], vCenter);
                VEC::StoreFloat3(&corners[i], C);
            }
        }

        _Use_decl_annotations_
        FORCE_INLINE ContainmentType VEC_CALLCONV BoundingBox::Contains(A_VECTOR point) const noexcept
        {
            VECTOR vCenter = VEC::LoadFloat3(&center);
            VECTOR vExtents = VEC::LoadFloat3(&extents);

            return VEC3::InBounds(VEC::Subtract(point, vCenter), vExtents) ? CONTAINS : DISJOINT;
        }

        _Use_decl_annotations_
        FORCE_INLINE ContainmentType VEC_CALLCONV BoundingBox::Contains(A_VECTOR V0, A_VECTOR V1, A_VECTOR V2) const noexcept
        {
            if (!Intersects(V0, V1, V2))
                return DISJOINT;

            VECTOR vCenter = VEC::LoadFloat3(&center);
            VECTOR vExtents = VEC::LoadFloat3(&extents);

            VECTOR d = VEC::Abs(VEC::Subtract(V0, vCenter));
            VECTOR Inside = VEC::LessOrEqual(d, vExtents);

            d = VEC::Abs(VEC::Subtract(V1, vCenter));
            Inside = VEC::AndInt(Inside, VEC::LessOrEqual(d, vExtents));

            d = VEC::Abs(VEC::Subtract(V2, vCenter));
            Inside = VEC::AndInt(Inside, VEC::LessOrEqual(d, vExtents));

            return (VEC3::EqualInt(Inside, VEC::TrueInt())) ? CONTAINS : INTERSECTS;
        }

        _Use_decl_annotations_
        FORCE_INLINE ContainmentType BoundingBox::Contains(const BoundingSphere& sh) const noexcept
        {
            VECTOR SphereCenter = VEC::LoadFloat3(&sh.center);
            VECTOR SphereRadius = VEC::ReplicatePtr(&sh.radius);

            VECTOR BoxCenter = VEC::LoadFloat3(&center);
            VECTOR BoxExtents = VEC::LoadFloat3(&extents);

            VECTOR BoxMin = VEC::Subtract(BoxCenter, BoxExtents);
            VECTOR BoxMax = VEC::Add(BoxCenter, BoxExtents);

            // Find the distance to the nearest point on the box.
            // for each i in (x, y, z)
            // if (SphereCenter(i) < BoxMin(i)) d2 += (SphereCenter(i) - BoxMin(i)) ^ 2
            // else if (SphereCenter(i) > BoxMax(i)) d2 += (SphereCenter(i) - BoxMax(i)) ^ 2

            VECTOR d = VEC::Zero();

            // Compute d for each dimension.
            VECTOR LessThanMin = VEC::Less(SphereCenter, BoxMin);
            VECTOR GreaterThanMax = VEC::Greater(SphereCenter, BoxMax);

            VECTOR MinDelta = VEC::Subtract(SphereCenter, BoxMin);
            VECTOR MaxDelta = VEC::Subtract(SphereCenter, BoxMax);

            // Choose value for each dimension based on the comparison.
            d = VEC::Select(d, MinDelta, LessThanMin);
            d = VEC::Select(d, MaxDelta, GreaterThanMax);

            // Use a dot-product to square them and sum them together.
            VECTOR d2 = VEC3::Dot(d, d);

            if (VEC3::Greater(d2, VEC::Multiply(SphereRadius, SphereRadius)))
                return DISJOINT;

            VECTOR InsideAll = VEC::LessOrEqual(VEC::Add(BoxMin, SphereRadius), SphereCenter);
            InsideAll = VEC::AndInt(InsideAll, VEC::LessOrEqual(SphereCenter, VEC::Subtract(BoxMax, SphereRadius)));
            InsideAll = VEC::AndInt(InsideAll, VEC::Greater(VEC::Subtract(BoxMax, BoxMin), SphereRadius));

            return (VEC3::EqualInt(InsideAll, VEC::TrueInt())) ? CONTAINS : INTERSECTS;
        }

        _Use_decl_annotations_
        FORCE_INLINE ContainmentType BoundingBox::Contains(const BoundingBox& box) const noexcept
        {
            VECTOR CenterA = VEC::LoadFloat3(&center);
            VECTOR ExtentsA = VEC::LoadFloat3(&extents);

            VECTOR CenterB = VEC::LoadFloat3(&box.center);
            VECTOR ExtentsB = VEC::LoadFloat3(&box.extents);

            VECTOR MinA = VEC::Subtract(CenterA, ExtentsA);
            VECTOR MaxA = VEC::Add(CenterA, ExtentsA);

            VECTOR MinB = VEC::Subtract(CenterB, ExtentsB);
            VECTOR MaxB = VEC::Add(CenterB, ExtentsB);

            // for each i in (x, y, z) if a_min(i) > b_max(i) or b_min(i) > a_max(i) then return false
            VECTOR Disjoint = VEC::OrInt(VEC::Greater(MinA, MaxB), VEC::Greater(MinB, MaxA));

            if (Vector3AnyTrue(Disjoint))
                return DISJOINT;

            // for each i in (x, y, z) if a_min(i) <= b_min(i) and b_max(i) <= a_max(i) then A contains B
            VECTOR Inside = VEC::AndInt(VEC::LessOrEqual(MinA, MinB), VEC::LessOrEqual(MaxB, MaxA));

            return Vector3AllTrue(Inside) ? CONTAINS : INTERSECTS;
        }

        _Use_decl_annotations_
        FORCE_INLINE ContainmentType BoundingBox::Contains(const BoundingOrientedBox& box) const noexcept
        {
            if (!box.Intersects(*this))
                return DISJOINT;

            VECTOR vCenter = VEC::LoadFloat3(&center);
            VECTOR vExtents = VEC::LoadFloat3(&extents);

            // Subtract off the AABB center to remove a subtract below
            VECTOR oCenter = VEC::Subtract(VEC::LoadFloat3(&box.center), vCenter);

            VECTOR oExtents = VEC::LoadFloat3(&box.extents);
            VECTOR oOrientation = VEC::LoadFloat4(&box.orientation);

#if defined(DEBUG) || defined(_DEBUG)
            assert(QuaternionIsUnit(oOrientation));
#endif // DEBUG

            VECTOR Inside = VEC::TrueInt();

            for (size_t i = 0; i < BoundingOrientedBox::CORNER_COUNT; ++i)
            {
                VECTOR C = VEC::Add(VEC3::Rotate(VEC::Multiply(oExtents, g_BoxOffset[i]), oOrientation), oCenter);
                VECTOR d = VEC::Abs(C);
                Inside = VEC::AndInt(Inside, VEC::LessOrEqual(d, vExtents));
            }

            return (VEC3::EqualInt(Inside, VEC::TrueInt())) ? CONTAINS : INTERSECTS;
        }

        _Use_decl_annotations_
        FORCE_INLINE ContainmentType BoundingBox::Contains(const BoundingFrustum& fr) const noexcept
        {
            if (!fr.Intersects(*this))
                return DISJOINT;

            Float3 Corners[BoundingFrustum::CORNER_COUNT];
            fr.GetCorners(Corners);

            VECTOR vCenter = VEC::LoadFloat3(&center);
            VECTOR vExtents = VEC::LoadFloat3(&extents);

            VECTOR Inside = VEC::TrueInt();

            for (size_t i = 0; i < BoundingFrustum::CORNER_COUNT; ++i)
            {
                VECTOR Point = VEC::LoadFloat3(&Corners[i]);
                VECTOR d = VEC::Abs(VEC::Subtract(Point, vCenter));
                Inside = VEC::AndInt(Inside, VEC::LessOrEqual(d, vExtents));
            }

            return (VEC3::EqualInt(Inside, VEC::TrueInt())) ? CONTAINS : INTERSECTS;
        }

        _Use_decl_annotations_
        FORCE_INLINE bool BoundingBox::Intersects(const BoundingSphere& sh) const noexcept
        {
            VECTOR SphereCenter = VEC::LoadFloat3(&sh.center);
            VECTOR SphereRadius = VEC::ReplicatePtr(&sh.radius);

            VECTOR BoxCenter = VEC::LoadFloat3(&center);
            VECTOR BoxExtents = VEC::LoadFloat3(&extents);

            VECTOR BoxMin = VEC::Subtract(BoxCenter, BoxExtents);
            VECTOR BoxMax = VEC::Add(BoxCenter, BoxExtents);

            // Find the distance to the nearest point on the box.
            // for each i in (x, y, z)
            // if (SphereCenter(i) < BoxMin(i)) d2 += (SphereCenter(i) - BoxMin(i)) ^ 2
            // else if (SphereCenter(i) > BoxMax(i)) d2 += (SphereCenter(i) - BoxMax(i)) ^ 2

            VECTOR d = VEC::Zero();

            // Compute d for each dimension.
            VECTOR LessThanMin = VEC::Less(SphereCenter, BoxMin);
            VECTOR GreaterThanMax = VEC::Greater(SphereCenter, BoxMax);

            VECTOR MinDelta = VEC::Subtract(SphereCenter, BoxMin);
            VECTOR MaxDelta = VEC::Subtract(SphereCenter, BoxMax);

            // Choose value for each dimension based on the comparison.
            d = VEC::Select(d, MinDelta, LessThanMin);
            d = VEC::Select(d, MaxDelta, GreaterThanMax);

            // Use a dot-product to square them and sum them together.
            VECTOR d2 = VEC3::Dot(d, d);

            return VEC3::LessOrEqual(d2, VEC::Multiply(SphereRadius, SphereRadius));
        }

        _Use_decl_annotations_
        FORCE_INLINE bool BoundingBox::Intersects(const BoundingBox& box) const noexcept
        {
            VECTOR CenterA = VEC::LoadFloat3(&center);
            VECTOR ExtentsA = VEC::LoadFloat3(&extents);

            VECTOR CenterB = VEC::LoadFloat3(&box.center);
            VECTOR ExtentsB = VEC::LoadFloat3(&box.extents);

            VECTOR MinA = VEC::Subtract(CenterA, ExtentsA);
            VECTOR MaxA = VEC::Add(CenterA, ExtentsA);

            VECTOR MinB = VEC::Subtract(CenterB, ExtentsB);
            VECTOR MaxB = VEC::Add(CenterB, ExtentsB);

            // for each i in (x, y, z) if a_min(i) > b_max(i) or b_min(i) > a_max(i) then return false
            VECTOR Disjoint = VEC::OrInt(VEC::Greater(MinA, MaxB), VEC::Greater(MinB, MaxA));

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
            VECTOR zero = VEC::Zero();

            // Load the box.
            VECTOR vCenter = VEC::LoadFloat3(&center);
            VECTOR vExtents = VEC::LoadFloat3(&extents);

            VECTOR BoxMin = VEC::Subtract(vCenter, vExtents);
            VECTOR BoxMax = VEC::Add(vCenter, vExtents);

            // Test the axes of the box (in effect test the AAB against the minimal AAB
            // around the triangle).
            VECTOR TriMin = VEC::Min(VEC::Min(V0, V1), V2);
            VECTOR TriMax = VEC::Max(VEC::Max(V0, V1), V2);

            // for each i in (x, y, z) if a_min(i) > b_max(i) or b_min(i) > a_max(i) then disjoint
            VECTOR Disjoint = VEC::OrInt(VEC::Greater(TriMin, BoxMax), VEC::Greater(BoxMin, TriMax));
            if (Vector3AnyTrue(Disjoint))
                return false;

            // Test the plane of the triangle.
            VECTOR Normal = VEC3::Cross(VEC::Subtract(V1, V0), VEC::Subtract(V2, V0));
            VECTOR Dist = VEC3::Dot(Normal, V0);

#if defined(DEBUG) || defined(_DEBUG)
            // Assert that the triangle is not degenerate.
            assert(!VEC3::Equal(Normal, zero));
#endif // DEBUG

            // for each i in (x, y, z) if n(i) >= 0 then v_min(i)=b_min(i), v_max(i)=b_max(i)
            // else v_min(i)=b_max(i), v_max(i)=b_min(i)
            VECTOR NormalSelect = VEC::Greater(Normal, zero);
            VECTOR V_Min = VEC::Select(BoxMax, BoxMin, NormalSelect);
            VECTOR V_Max = VEC::Select(BoxMin, BoxMax, NormalSelect);

            // if n dot v_min + d > 0 || n dot v_max + d < 0 then disjoint
            VECTOR MinDist = VEC3::Dot(V_Min, Normal);
            VECTOR MaxDist = VEC3::Dot(V_Max, Normal);

            VECTOR NoIntersection = VEC::Greater(MinDist, Dist);
            NoIntersection = VEC::OrInt(NoIntersection, VEC::Less(MaxDist, Dist));

            // Move the box center to zero to simplify the following tests.
            VECTOR TV0 = VEC::Subtract(V0, vCenter);
            VECTOR TV1 = VEC::Subtract(V1, vCenter);
            VECTOR TV2 = VEC::Subtract(V2, vCenter);

            // Test the edge/edge axes (3*3).
            VECTOR e0 = VEC::Subtract(TV1, TV0);
            VECTOR e1 = VEC::Subtract(TV2, TV1);
            VECTOR e2 = VEC::Subtract(TV0, TV2);

            // Make w zero.
            e0 = VEC::Insert<0, 0, 0, 0, 1>(e0, zero);
            e1 = VEC::Insert<0, 0, 0, 0, 1>(e1, zero);
            e2 = VEC::Insert<0, 0, 0, 0, 1>(e2, zero);

            VECTOR Axis;
            VECTOR p0, p1, p2;
            VECTOR Min, Max;
            VECTOR Radius;

            // Axis == (1,0,0) x e0 = (0, -e0.z, e0.y)
            Axis = VEC::Permute<PERMUTE_0W, PERMUTE_1Z, PERMUTE_0Y, PERMUTE_0X>(e0, VEC::Negate(e0));
            p0 = VEC3::Dot(TV0, Axis);
            // p1 = XMVector3Dot( V1, Axis ); // p1 = p0;
            p2 = VEC3::Dot(TV2, Axis);
            Min = VEC::Min(p0, p2);
            Max = VEC::Max(p0, p2);
            Radius = VEC3::Dot(vExtents, VEC::Abs(Axis));
            NoIntersection = VEC::OrInt(NoIntersection, VEC::Greater(Min, Radius));
            NoIntersection = VEC::OrInt(NoIntersection, VEC::Less(Max, VEC::Negate(Radius)));

            // Axis == (1,0,0) x e1 = (0, -e1.z, e1.y)
            Axis = VEC::Permute<PERMUTE_0W, PERMUTE_1Z, PERMUTE_0Y, PERMUTE_0X>(e1, VEC::Negate(e1));
            p0 = VEC3::Dot(TV0, Axis);
            p1 = VEC3::Dot(TV1, Axis);
            // p2 = XMVector3Dot( V2, Axis ); // p2 = p1;
            Min = VEC::Min(p0, p1);
            Max = VEC::Max(p0, p1);
            Radius = VEC3::Dot(vExtents, VEC::Abs(Axis));
            NoIntersection = VEC::OrInt(NoIntersection, VEC::Greater(Min, Radius));
            NoIntersection = VEC::OrInt(NoIntersection, VEC::Less(Max, VEC::Negate(Radius)));

            // Axis == (1,0,0) x e2 = (0, -e2.z, e2.y)
            Axis = VEC::Permute<PERMUTE_0W, PERMUTE_1Z, PERMUTE_0Y, PERMUTE_0X>(e2, VEC::Negate(e2));
            p0 = VEC3::Dot(TV0, Axis);
            p1 = VEC3::Dot(TV1, Axis);
            // p2 = XMVector3Dot( V2, Axis ); // p2 = p0;
            Min = VEC::Min(p0, p1);
            Max = VEC::Max(p0, p1);
            Radius = VEC3::Dot(vExtents, VEC::Abs(Axis));
            NoIntersection = VEC::OrInt(NoIntersection, VEC::Greater(Min, Radius));
            NoIntersection = VEC::OrInt(NoIntersection, VEC::Less(Max, VEC::Negate(Radius)));

            // Axis == (0,1,0) x e0 = (e0.z, 0, -e0.x)
            Axis = VEC::Permute<PERMUTE_0Z, PERMUTE_0W, PERMUTE_1X, PERMUTE_0Y>(e0, VEC::Negate(e0));
            p0 = VEC3::Dot(TV0, Axis);
            // p1 = XMVector3Dot( V1, Axis ); // p1 = p0;
            p2 = VEC3::Dot(TV2, Axis);
            Min = VEC::Min(p0, p2);
            Max = VEC::Max(p0, p2);
            Radius = VEC3::Dot(vExtents, VEC::Abs(Axis));
            NoIntersection = VEC::OrInt(NoIntersection, VEC::Greater(Min, Radius));
            NoIntersection = VEC::OrInt(NoIntersection, VEC::Less(Max, VEC::Negate(Radius)));

            // Axis == (0,1,0) x e1 = (e1.z, 0, -e1.x)
            Axis = VEC::Permute<PERMUTE_0Z, PERMUTE_0W, PERMUTE_1X, PERMUTE_0Y>(e1, VEC::Negate(e1));
            p0 = VEC3::Dot(TV0, Axis);
            p1 = VEC3::Dot(TV1, Axis);
            // p2 = XMVector3Dot( V2, Axis ); // p2 = p1;
            Min = VEC::Min(p0, p1);
            Max = VEC::Max(p0, p1);
            Radius = VEC3::Dot(vExtents, VEC::Abs(Axis));
            NoIntersection = VEC::OrInt(NoIntersection, VEC::Greater(Min, Radius));
            NoIntersection = VEC::OrInt(NoIntersection, VEC::Less(Max, VEC::Negate(Radius)));

            // Axis == (0,0,1) x e2 = (e2.z, 0, -e2.x)
            Axis = VEC::Permute<PERMUTE_0Z, PERMUTE_0W, PERMUTE_1X, PERMUTE_0Y>(e2, VEC::Negate(e2));
            p0 = VEC3::Dot(TV0, Axis);
            p1 = VEC3::Dot(TV1, Axis);
            // p2 = XMVector3Dot( V2, Axis ); // p2 = p0;
            Min = VEC::Min(p0, p1);
            Max = VEC::Max(p0, p1);
            Radius = VEC3::Dot(vExtents, VEC::Abs(Axis));
            NoIntersection = VEC::OrInt(NoIntersection, VEC::Greater(Min, Radius));
            NoIntersection = VEC::OrInt(NoIntersection, VEC::Less(Max, VEC::Negate(Radius)));

            // Axis == (0,0,1) x e0 = (-e0.y, e0.x, 0)
            Axis = VEC::Permute<PERMUTE_1Y, PERMUTE_0X, PERMUTE_0W, PERMUTE_0Z>(e0, VEC::Negate(e0));
            p0 = VEC3::Dot(TV0, Axis);
            // p1 = XMVector3Dot( V1, Axis ); // p1 = p0;
            p2 = VEC3::Dot(TV2, Axis);
            Min = VEC::Min(p0, p2);
            Max = VEC::Max(p0, p2);
            Radius = VEC3::Dot(vExtents, VEC::Abs(Axis));
            NoIntersection = VEC::OrInt(NoIntersection, VEC::Greater(Min, Radius));
            NoIntersection = VEC::OrInt(NoIntersection, VEC::Less(Max, VEC::Negate(Radius)));

            // Axis == (0,0,1) x e1 = (-e1.y, e1.x, 0)
            Axis = VEC::Permute<PERMUTE_1Y, PERMUTE_0X, PERMUTE_0W, PERMUTE_0Z>(e1, VEC::Negate(e1));
            p0 = VEC3::Dot(TV0, Axis);
            p1 = VEC3::Dot(TV1, Axis);
            // p2 = XMVector3Dot( V2, Axis ); // p2 = p1;
            Min = VEC::Min(p0, p1);
            Max = VEC::Max(p0, p1);
            Radius = VEC3::Dot(vExtents, VEC::Abs(Axis));
            NoIntersection = VEC::OrInt(NoIntersection, VEC::Greater(Min, Radius));
            NoIntersection = VEC::OrInt(NoIntersection, VEC::Less(Max, VEC::Negate(Radius)));

            // Axis == (0,0,1) x e2 = (-e2.y, e2.x, 0)
            Axis = VEC::Permute<PERMUTE_1Y, PERMUTE_0X, PERMUTE_0W, PERMUTE_0Z>(e2, VEC::Negate(e2));
            p0 = VEC3::Dot(TV0, Axis);
            p1 = VEC3::Dot(TV1, Axis);
            // p2 = XMVector3Dot( V2, Axis ); // p2 = p0;
            Min = VEC::Min(p0, p1);
            Max = VEC::Max(p0, p1);
            Radius = VEC3::Dot(vExtents, VEC::Abs(Axis));
            NoIntersection = VEC::OrInt(NoIntersection, VEC::Greater(Min, Radius));
            NoIntersection = VEC::OrInt(NoIntersection, VEC::Less(Max, VEC::Negate(Radius)));

            return VEC4::NotEqualInt(NoIntersection, VEC::TrueInt());
        }

        _Use_decl_annotations_
        FORCE_INLINE PlaneIntersectionType VEC_CALLCONV BoundingBox::Intersects(A_VECTOR plane) const noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(PlaneIsUnit(plane));
#endif // DEBUG

            // Load the box.
            VECTOR vCenter = VEC::LoadFloat3(&center);
            VECTOR vExtents = VEC::LoadFloat3(&extents);

            // Set w of the center to one so we can dot4 with a plane.
            vCenter = VEC::Insert<0, 0, 0, 0, 1>(vCenter, VEC::SplatOne());

            VECTOR Outside, Inside;
            FastIntersectAxisAlignedBoxPlane(vCenter, vExtents, plane, Outside, Inside);

            // If the box is outside any plane it is outside.
            if (VEC4::EqualInt(Outside, VEC::TrueInt()))
                return FRONT;

            // If the box is inside all planes it is inside.
            if (VEC4::EqualInt(Inside, VEC::TrueInt()))
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
            VECTOR vCenter = VEC::LoadFloat3(&center);
            VECTOR vExtents = VEC::LoadFloat3(&extents);

            // Adjust ray origin to be relative to center of the box.
            VECTOR TOrigin = VEC::Subtract(vCenter, origin);

            // Compute the dot product againt each axis of the box.
            // Since the axii are (1,0,0), (0,1,0), (0,0,1) no computation is necessary.
            VECTOR AxisDotOrigin = TOrigin;
            VECTOR AxisDotDirection = direction;

            // if (fabs(AxisDotDirection) <= Epsilon) the ray is nearly parallel to the slab.
            VECTOR IsParallel = VEC::LessOrEqual(VEC::Abs(AxisDotDirection), g_RayEpsilon);

            // Test against all three axii simultaneously.
            VECTOR InverseAxisDotDirection = VEC::Reciprocal(AxisDotDirection);
            VECTOR t1 = VEC::Multiply(VEC::Subtract(AxisDotOrigin, vExtents), InverseAxisDotDirection);
            VECTOR t2 = VEC::Multiply(VEC::Add(AxisDotOrigin, vExtents), InverseAxisDotDirection);

            // Compute the max of min(t1,t2) and the min of max(t1,t2) ensuring we don't
            // use the results from any directions parallel to the slab.
            VECTOR t_min = VEC::Select(VEC::Min(t1, t2), g_FltMin, IsParallel);
            VECTOR t_max = VEC::Select(VEC::Max(t1, t2), g_FltMax, IsParallel);

            // t_min.x = maximum( t_min.x, t_min.y, t_min.z );
            // t_max.x = minimum( t_max.x, t_max.y, t_max.z );
            t_min = VEC::Max(t_min, VEC::SplatY(t_min));  // x = max(x,y)
            t_min = VEC::Max(t_min, VEC::SplatZ(t_min));  // x = max(max(x,y),z)
            t_max = VEC::Min(t_max, VEC::SplatY(t_max));  // x = min(x,y)
            t_max = VEC::Min(t_max, VEC::SplatZ(t_max));  // x = min(min(x,y),z)

            // if ( t_min > t_max ) return false;
            VECTOR NoIntersection = VEC::Greater(VEC::SplatX(t_min), VEC::SplatX(t_max));

            // if ( t_max < 0.0f ) return false;
            NoIntersection = VEC::OrInt(NoIntersection, VEC::Less(VEC::SplatX(t_max), VEC::Zero()));

            // if (IsParallel && (-Extents > AxisDotOrigin || Extents < AxisDotOrigin)) return false;
            VECTOR ParallelOverlap = VEC::InBounds(AxisDotOrigin, vExtents);
            NoIntersection = VEC::OrInt(NoIntersection, VEC::AndCInt(IsParallel, ParallelOverlap));

            if (!Vector3AnyTrue(NoIntersection))
            {
                // Store the x-component to *pDist
                VEC::StoreFloat(&dist, t_min);
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
            VECTOR vCenter = VEC::LoadFloat3(&center);
            VECTOR vExtents = VEC::LoadFloat3(&extents);

            // Set w of the center to one so we can dot4 with a plane.
            vCenter = VEC::Insert<0, 0, 0, 0, 1>(vCenter, VEC::SplatOne());

            VECTOR Outside, Inside;

            // Test against each plane.
            FastIntersectAxisAlignedBoxPlane(vCenter, vExtents, plane0, Outside, Inside);

            VECTOR AnyOutside = Outside;
            VECTOR AllInside = Inside;

            FastIntersectAxisAlignedBoxPlane(vCenter, vExtents, plane1, Outside, Inside);
            AnyOutside = VEC::OrInt(AnyOutside, Outside);
            AllInside = VEC::AndInt(AllInside, Inside);

            FastIntersectAxisAlignedBoxPlane(vCenter, vExtents, plane2, Outside, Inside);
            AnyOutside = VEC::OrInt(AnyOutside, Outside);
            AllInside = VEC::AndInt(AllInside, Inside);

            FastIntersectAxisAlignedBoxPlane(vCenter, vExtents, plane3, Outside, Inside);
            AnyOutside = VEC::OrInt(AnyOutside, Outside);
            AllInside = VEC::AndInt(AllInside, Inside);

            FastIntersectAxisAlignedBoxPlane(vCenter, vExtents, plane4, Outside, Inside);
            AnyOutside = VEC::OrInt(AnyOutside, Outside);
            AllInside = VEC::AndInt(AllInside, Inside);

            FastIntersectAxisAlignedBoxPlane(vCenter, vExtents, plane5, Outside, Inside);
            AnyOutside = VEC::OrInt(AnyOutside, Outside);
            AllInside = VEC::AndInt(AllInside, Inside);

            // If the box is outside any plane it is outside.
            if (VEC4::EqualInt(AnyOutside, VEC::TrueInt()))
                return DISJOINT;

            // If the box is inside all planes it is inside.
            if (VEC4::EqualInt(AllInside, VEC::TrueInt()))
                return CONTAINS;

            // The box is not inside all planes or outside a plane, it may intersect.
            return INTERSECTS;
        }

        _Use_decl_annotations_
        FORCE_INLINE void BoundingBox::CreateMerged(BoundingBox& out, const BoundingBox& b1, const BoundingBox& b2) noexcept
        {
            VECTOR b1Center = VEC::LoadFloat3(&b1.center);
            VECTOR b1Extents = VEC::LoadFloat3(&b1.extents);

            VECTOR b2Center = VEC::LoadFloat3(&b2.center);
            VECTOR b2Extents = VEC::LoadFloat3(&b2.extents);

            VECTOR Min = VEC::Subtract(b1Center, b1Extents);
            Min = VEC::Min(Min, VEC::Subtract(b2Center, b2Extents));

            VECTOR Max = VEC::Add(b1Center, b1Extents);
            Max = VEC::Max(Max, VEC::Add(b2Center, b2Extents));

#if defined(DEBUG) || defined(_DEBUG)
            assert(VEC3::LessOrEqual(Min, Max));
#endif // DEBUG

            VEC::StoreFloat3(&out.center, VEC::Scale(VEC::Add(Min, Max), 0.5f));
            VEC::StoreFloat3(&out.extents, VEC::Scale(VEC::Subtract(Max, Min), 0.5f));
        }

        _Use_decl_annotations_
        FORCE_INLINE void BoundingBox::CreateFromSphere(BoundingBox& out, const BoundingSphere& sh) noexcept
        {
            VECTOR spCenter = VEC::LoadFloat3(&sh.center);
            VECTOR shRadius = VEC::ReplicatePtr(&sh.radius);

            VECTOR Min = VEC::Subtract(spCenter, shRadius);
            VECTOR Max = VEC::Add(spCenter, shRadius);

#if defined(DEBUG) || defined(_DEBUG)
            assert(VEC3::LessOrEqual(Min, Max));
#endif // DEBUG

            VEC::StoreFloat3(&out.center, VEC::Scale(VEC::Add(Min, Max), 0.5f));
            VEC::StoreFloat3(&out.extents, VEC::Scale(VEC::Subtract(Max, Min), 0.5f));
        }

        _Use_decl_annotations_
        FORCE_INLINE void VEC_CALLCONV BoundingBox::CreateFromPoints(BoundingBox& out, A_VECTOR pt1, A_VECTOR pt2) noexcept
        {
            VECTOR Min = VEC::Min(pt1, pt2);
            VECTOR Max = VEC::Max(pt1, pt2);

            // Store center and extents.
            VEC::StoreFloat3(&out.center, VEC::Scale(VEC::Add(Min, Max), 0.5f));
            VEC::StoreFloat3(&out.extents, VEC::Scale(VEC::Subtract(Max, Min), 0.5f));
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

            vMin = vMax = VEC::LoadFloat3(pPoints);

            for (size_t i = 1; i < count; ++i)
            {
                VECTOR Point = VEC::LoadFloat3(reinterpret_cast<const Float3*>(reinterpret_cast<const uint8_t*>(pPoints) + i * stride));

                vMin = VEC::Min(vMin, Point);
                vMax = VEC::Max(vMax, Point);
            }

            // Store center and extents.
            VEC::StoreFloat3(&out.center, VEC::Scale(VEC::Add(vMin, vMax), 0.5f));
            VEC::StoreFloat3(&out.extents, VEC::Scale(VEC::Subtract(vMax, vMin), 0.5f));
        }
        
        _Use_decl_annotations_
        FORCE_INLINE void VEC_CALLCONV BoundingOrientedBox::Transform(BoundingOrientedBox& out, A_MATRIX m) const noexcept
        {
            // Load the box.
            VECTOR vCenter = VEC::LoadFloat3(&center);
            VECTOR vExtents = VEC::LoadFloat3(&extents);
            VECTOR vOrientation = VEC::LoadFloat4(&orientation);

#if defined(DEBUG) || defined(_DEBUG)
            assert(QuaternionIsUnit(vOrientation));
#endif // DEBUG

            // Composite the box rotation and the transform rotation.
            MATRIX nM;
            nM.r[0] = VEC3::Normalize(m.r[0]);
            nM.r[1] = VEC3::Normalize(m.r[1]);
            nM.r[2] = VEC3::Normalize(m.r[2]);
            nM.r[3] = g_IdentityR3;
            VECTOR Rotation = Quaternion::RotationMatrix(nM);
            vOrientation = Quaternion::Multiply(vOrientation, Rotation);

            // Transform the center.
            vCenter = VEC3::Transform(vCenter, m);

            // Scale the box extents.
            VECTOR dX = VEC3::Length(m.r[0]);
            VECTOR dY = VEC3::Length(m.r[1]);
            VECTOR dZ = VEC3::Length(m.r[2]);

            VECTOR VectorScale = VEC::Select(dY, dX, g_Select1000);
            VectorScale = VEC::Select(dZ, VectorScale, g_Select1100);
            vExtents = VEC::Multiply(vExtents, VectorScale);

            // Store the box.
            VEC::StoreFloat3(&out.center, vCenter);
            VEC::StoreFloat3(&out.extents, vExtents);
            VEC::StoreFloat4(&out.orientation, vOrientation);
        }

        _Use_decl_annotations_
        FORCE_INLINE void VEC_CALLCONV BoundingOrientedBox::Transform(BoundingOrientedBox& out, float scale, A_VECTOR rotation, A_VECTOR translation) const noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(QuaternionIsUnit(rotation));
#endif // DEBUG

            // Load the box.
            VECTOR vCenter = VEC::LoadFloat3(&center);
            VECTOR vExtents = VEC::LoadFloat3(&extents);
            VECTOR vOrientation = VEC::LoadFloat4(&orientation);

#if defined(DEBUG) || defined(_DEBUG)
            assert(QuaternionIsUnit(vOrientation));
#endif // DEBUG

            // Composite the box rotation and the transform rotation.
            vOrientation = Quaternion::Multiply(vOrientation, rotation);

            // Transform the center.
            VECTOR VectorScale = VEC::Replicate(scale);
            vCenter = VEC::Add(VEC3::Rotate(VEC::Multiply(vCenter, VectorScale), rotation), translation);

            // Scale the box extents.
            vExtents = VEC::Multiply(vExtents, VectorScale);

            // Store the box.
            VEC::StoreFloat3(&out.center, vCenter);
            VEC::StoreFloat3(&out.extents, vExtents);
            VEC::StoreFloat4(&out.orientation, vOrientation);
        }
        
        _Use_decl_annotations_
        FORCE_INLINE void BoundingOrientedBox::GetCorners(Float3* corners) const noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(corners != nullptr);
#endif // DEBUG

            // Load the box
            VECTOR vCenter = VEC::LoadFloat3(&center);
            VECTOR vExtents = VEC::LoadFloat3(&extents);
            VECTOR vOrientation = VEC::LoadFloat4(&orientation);

#if defined(DEBUG) || defined(_DEBUG)
            assert(QuaternionIsUnit(vOrientation));
#endif // DEBUG

            for (size_t i = 0; i < CORNER_COUNT; ++i)
            {
                VECTOR C = VEC::Add(VEC3::Rotate(VEC::Multiply(vExtents, g_BoxOffset[i]), vOrientation), vCenter);
                VEC::StoreFloat3(&corners[i], C);
            }
        }

        _Use_decl_annotations_
        FORCE_INLINE ContainmentType VEC_CALLCONV BoundingOrientedBox::Contains(A_VECTOR point) const noexcept
        {
            VECTOR vCenter = VEC::LoadFloat3(&center);
            VECTOR vExtents = VEC::LoadFloat3(&extents);
            VECTOR vOrientation = VEC::LoadFloat4(&orientation);

            // Transform the point to be local to the box.
            VECTOR TPoint = VEC3::InverseRotate(VEC::Subtract(point, vCenter), vOrientation);

            return VEC3::InBounds(TPoint, vExtents) ? CONTAINS : DISJOINT;
        }

        _Use_decl_annotations_
        FORCE_INLINE ContainmentType VEC_CALLCONV BoundingOrientedBox::Contains(A_VECTOR V0, A_VECTOR V1, A_VECTOR V2) const noexcept
        {
            // Load the box center & orientation.
            VECTOR vCenter = VEC::LoadFloat3(&center);
            VECTOR vOrientation = VEC::LoadFloat4(&orientation);

            // Transform the triangle vertices into the space of the box.
            VECTOR TV0 = VEC3::InverseRotate(VEC::Subtract(V0, vCenter), vOrientation);
            VECTOR TV1 = VEC3::InverseRotate(VEC::Subtract(V1, vCenter), vOrientation);
            VECTOR TV2 = VEC3::InverseRotate(VEC::Subtract(V2, vCenter), vOrientation);

            BoundingBox box;
            box.center = Float3(0.0f, 0.0f, 0.0f);
            box.extents = extents;

            // Use the triangle vs axis aligned box intersection routine.
            return box.Contains(TV0, TV1, TV2);
        }

        _Use_decl_annotations_
        FORCE_INLINE ContainmentType BoundingOrientedBox::Contains(const BoundingSphere& sh) const noexcept
        {
            VECTOR SphereCenter = VEC::LoadFloat3(&sh.center);
            VECTOR SphereRadius = VEC::ReplicatePtr(&sh.radius);

            VECTOR BoxCenter = VEC::LoadFloat3(&center);
            VECTOR BoxExtents = VEC::LoadFloat3(&extents);
            VECTOR BoxOrientation = VEC::LoadFloat4(&orientation);

#if defined(DEBUG) || defined(_DEBUG)
            assert(QuaternionIsUnit(BoxOrientation));
#endif // DEBUG

            // Transform the center of the sphere to be local to the box.
            // BoxMin = -BoxExtents
            // BoxMax = +BoxExtents
            SphereCenter = VEC3::InverseRotate(VEC::Subtract(SphereCenter, BoxCenter), BoxOrientation);

            // Find the distance to the nearest point on the box.
            // for each i in (x, y, z)
            // if (SphereCenter(i) < BoxMin(i)) d2 += (SphereCenter(i) - BoxMin(i)) ^ 2
            // else if (SphereCenter(i) > BoxMax(i)) d2 += (SphereCenter(i) - BoxMax(i)) ^ 2

            VECTOR d = VEC::Zero();

            // Compute d for each dimension.
            VECTOR LessThanMin = VEC::Less(SphereCenter, VEC::Negate(BoxExtents));
            VECTOR GreaterThanMax = VEC::Greater(SphereCenter, BoxExtents);

            VECTOR MinDelta = VEC::Add(SphereCenter, BoxExtents);
            VECTOR MaxDelta = VEC::Subtract(SphereCenter, BoxExtents);

            // Choose value for each dimension based on the comparison.
            d = VEC::Select(d, MinDelta, LessThanMin);
            d = VEC::Select(d, MaxDelta, GreaterThanMax);

            // Use a dot-product to square them and sum them together.
            VECTOR d2 = VEC3::Dot(d, d);
            VECTOR SphereRadiusSq = VEC::Multiply(SphereRadius, SphereRadius);

            if (VEC4::Greater(d2, SphereRadiusSq))
                return DISJOINT;

            // See if we are completely inside the box
            VECTOR SMin = VEC::Subtract(SphereCenter, SphereRadius);
            VECTOR SMax = VEC::Add(SphereCenter, SphereRadius);

            return (VEC3::InBounds(SMin, BoxExtents) && VEC3::InBounds(SMax, BoxExtents)) ? CONTAINS : INTERSECTS;
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
            VECTOR aCenter = VEC::LoadFloat3(&center);
            VECTOR aExtents = VEC::LoadFloat3(&extents);
            VECTOR aOrientation = VEC::LoadFloat4(&orientation);

#if defined(DEBUG) || defined(_DEBUG)
            assert(QuaternionIsUnit(aOrientation));
#endif // DEBUG

            VECTOR bCenter = VEC::LoadFloat3(&box.center);
            VECTOR bExtents = VEC::LoadFloat3(&box.extents);
            VECTOR bOrientation = VEC::LoadFloat4(&box.orientation);

#if defined(DEBUG) || defined(_DEBUG)
            assert(QuaternionIsUnit(bOrientation));
#endif // DEBUG

            VECTOR offset = VEC::Subtract(bCenter, aCenter);

            for (size_t i = 0; i < CORNER_COUNT; ++i)
            {
                // Cb = rotate( bExtents * corneroffset[i], bOrientation ) + bcenter
                // Ca = invrotate( Cb - aCenter, aOrientation )

                VECTOR C = VEC::Add(VEC3::Rotate(VEC::Multiply(bExtents, g_BoxOffset[i]), bOrientation), offset);
                C = VEC3::InverseRotate(C, aOrientation);

                if (!VEC3::InBounds(C, aExtents))
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
            VECTOR vCenter = VEC::LoadFloat3(&center);
            VECTOR vExtents = VEC::LoadFloat3(&extents);
            VECTOR vOrientation = VEC::LoadFloat4(&orientation);

#if defined(DEBUG) || defined(_DEBUG)
            assert(QuaternionIsUnit(vOrientation));
#endif // DEBUG

            for (size_t i = 0; i < BoundingFrustum::CORNER_COUNT; ++i)
            {
                VECTOR C = VEC3::InverseRotate(VEC::Subtract(VEC::LoadFloat3(&Corners[i]), vCenter), vOrientation);

                if (!VEC3::InBounds(C, vExtents))
                    return INTERSECTS;
            }

            return CONTAINS;
        }

        _Use_decl_annotations_
        FORCE_INLINE bool BoundingOrientedBox::Intersects(const BoundingSphere& sh) const noexcept
        {
            VECTOR SphereCenter = VEC::LoadFloat3(&sh.center);
            VECTOR SphereRadius = VEC::ReplicatePtr(&sh.radius);

            VECTOR BoxCenter = VEC::LoadFloat3(&center);
            VECTOR BoxExtents = VEC::LoadFloat3(&extents);
            VECTOR BoxOrientation = VEC::LoadFloat4(&orientation);

#if defined(DEBUG) || defined(_DEBUG)
            assert(QuaternionIsUnit(BoxOrientation));
#endif // DEBUG

            // Transform the center of the sphere to be local to the box.
            // BoxMin = -BoxExtents
            // BoxMax = +BoxExtents
            SphereCenter = VEC3::InverseRotate(VEC::Subtract(SphereCenter, BoxCenter), BoxOrientation);

            // Find the distance to the nearest point on the box.
            // for each i in (x, y, z)
            // if (SphereCenter(i) < BoxMin(i)) d2 += (SphereCenter(i) - BoxMin(i)) ^ 2
            // else if (SphereCenter(i) > BoxMax(i)) d2 += (SphereCenter(i) - BoxMax(i)) ^ 2

            VECTOR d = VEC::Zero();

            // Compute d for each dimension.
            VECTOR LessThanMin = VEC::Less(SphereCenter, VEC::Negate(BoxExtents));
            VECTOR GreaterThanMax = VEC::Greater(SphereCenter, BoxExtents);

            VECTOR MinDelta = VEC::Add(SphereCenter, BoxExtents);
            VECTOR MaxDelta = VEC::Subtract(SphereCenter, BoxExtents);

            // Choose value for each dimension based on the comparison.
            d = VEC::Select(d, MinDelta, LessThanMin);
            d = VEC::Select(d, MaxDelta, GreaterThanMax);

            // Use a dot-product to square them and sum them together.
            VECTOR d2 = VEC3::Dot(d, d);

            return VEC4::LessOrEqual(d2, VEC::Multiply(SphereRadius, SphereRadius)) ? true : false;
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
            VECTOR A_quat = VEC::LoadFloat4(&orientation);
            VECTOR B_quat = VEC::LoadFloat4(&box.orientation);

#if defined(DEBUG) || defined(_DEBUG)
            assert(QuaternionIsUnit(A_quat));
            assert(QuaternionIsUnit(B_quat));
#endif // DEBUG

            VECTOR Q = Quaternion::Multiply(A_quat, Quaternion::Conjugate(B_quat));
            MATRIX R = MAT::RotationQuaternion(Q);

            // Compute the translation of B relative to A.
            VECTOR A_cent = VEC::LoadFloat3(&center);
            VECTOR B_cent = VEC::LoadFloat3(&box.center);
            VECTOR t = VEC3::InverseRotate(VEC::Subtract(B_cent, A_cent), A_quat);

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
            VECTOR h_A = VEC::LoadFloat3(&extents);
            VECTOR h_B = VEC::LoadFloat3(&box.extents);

            // Rows. Note R[0,1,2]X.w = 0.
            VECTOR R0X = R.r[0];
            VECTOR R1X = R.r[1];
            VECTOR R2X = R.r[2];

            R = MAT::Transpose(R);

            // Columns. Note RX[0,1,2].w = 0.
            VECTOR RX0 = R.r[0];
            VECTOR RX1 = R.r[1];
            VECTOR RX2 = R.r[2];

            // Absolute value of rows.
            VECTOR AR0X = VEC::Abs(R0X);
            VECTOR AR1X = VEC::Abs(R1X);
            VECTOR AR2X = VEC::Abs(R2X);

            // Absolute value of columns.
            VECTOR ARX0 = VEC::Abs(RX0);
            VECTOR ARX1 = VEC::Abs(RX1);
            VECTOR ARX2 = VEC::Abs(RX2);

            // Test each of the 15 possible seperating axii.
            VECTOR d, d_A, d_B;

            // l = a(u) = (1, 0, 0)
            // t dot l = t.x
            // d(A) = h(A).x
            // d(B) = h(B) dot abs(r00, r01, r02)
            d = VEC::SplatX(t);
            d_A = VEC::SplatX(h_A);
            d_B = VEC3::Dot(h_B, AR0X);
            VECTOR NoIntersection = VEC::Greater(VEC::Abs(d), VEC::Add(d_A, d_B));

            // l = a(v) = (0, 1, 0)
            // t dot l = t.y
            // d(A) = h(A).y
            // d(B) = h(B) dot abs(r10, r11, r12)
            d = VEC::SplatY(t);
            d_A = VEC::SplatY(h_A);
            d_B = VEC3::Dot(h_B, AR1X);
            NoIntersection = VEC::OrInt(NoIntersection,
                VEC::Greater(VEC::Abs(d), VEC::Add(d_A, d_B)));

            // l = a(w) = (0, 0, 1)
            // t dot l = t.z
            // d(A) = h(A).z
            // d(B) = h(B) dot abs(r20, r21, r22)
            d = VEC::SplatZ(t);
            d_A = VEC::SplatZ(h_A);
            d_B = VEC3::Dot(h_B, AR2X);
            NoIntersection = VEC::OrInt(NoIntersection,
                VEC::Greater(VEC::Abs(d), VEC::Add(d_A, d_B)));

            // l = b(u) = (r00, r10, r20)
            // d(A) = h(A) dot abs(r00, r10, r20)
            // d(B) = h(B).x
            d = VEC3::Dot(t, RX0);
            d_A = VEC3::Dot(h_A, ARX0);
            d_B = VEC::SplatX(h_B);
            NoIntersection = VEC::OrInt(NoIntersection,
                VEC::Greater(VEC::Abs(d), VEC::Add(d_A, d_B)));

            // l = b(v) = (r01, r11, r21)
            // d(A) = h(A) dot abs(r01, r11, r21)
            // d(B) = h(B).y
            d = VEC3::Dot(t, RX1);
            d_A = VEC3::Dot(h_A, ARX1);
            d_B = VEC::SplatY(h_B);
            NoIntersection = VEC::OrInt(NoIntersection,
                VEC::Greater(VEC::Abs(d), VEC::Add(d_A, d_B)));

            // l = b(w) = (r02, r12, r22)
            // d(A) = h(A) dot abs(r02, r12, r22)
            // d(B) = h(B).z
            d = VEC3::Dot(t, RX2);
            d_A = VEC3::Dot(h_A, ARX2);
            d_B = VEC::SplatZ(h_B);
            NoIntersection = VEC::OrInt(NoIntersection,
                VEC::Greater(VEC::Abs(d), VEC::Add(d_A, d_B)));

            // l = a(u) x b(u) = (0, -r20, r10)
            // d(A) = h(A) dot abs(0, r20, r10)
            // d(B) = h(B) dot abs(0, r02, r01)
            d = VEC3::Dot(t, VEC::Permute<PERMUTE_0W, PERMUTE_1Z, PERMUTE_0Y, PERMUTE_0X>(RX0, VEC::Negate(RX0)));
            d_A = VEC3::Dot(h_A, VEC::Swizzle<SWIZZLE_W, SWIZZLE_Z, SWIZZLE_Y, SWIZZLE_X>(ARX0));
            d_B = VEC3::Dot(h_B, VEC::Swizzle<SWIZZLE_W, SWIZZLE_Z, SWIZZLE_Y, SWIZZLE_X>(AR0X));
            NoIntersection = VEC::OrInt(NoIntersection,
                VEC::Greater(VEC::Abs(d), VEC::Add(d_A, d_B)));

            // l = a(u) x b(v) = (0, -r21, r11)
            // d(A) = h(A) dot abs(0, r21, r11)
            // d(B) = h(B) dot abs(r02, 0, r00)
            d = VEC3::Dot(t, VEC::Permute<PERMUTE_0W, PERMUTE_1Z, PERMUTE_0Y, PERMUTE_0X>(RX1, VEC::Negate(RX1)));
            d_A = VEC3::Dot(h_A, VEC::Swizzle<SWIZZLE_W, SWIZZLE_Z, SWIZZLE_Y, SWIZZLE_X>(ARX1));
            d_B = VEC3::Dot(h_B, VEC::Swizzle<SWIZZLE_Z, SWIZZLE_W, SWIZZLE_X, SWIZZLE_Y>(AR0X));
            NoIntersection = VEC::OrInt(NoIntersection,
                VEC::Greater(VEC::Abs(d), VEC::Add(d_A, d_B)));

            // l = a(u) x b(w) = (0, -r22, r12)
            // d(A) = h(A) dot abs(0, r22, r12)
            // d(B) = h(B) dot abs(r01, r00, 0)
            d = VEC3::Dot(t, VEC::Permute<PERMUTE_0W, PERMUTE_1Z, PERMUTE_0Y, PERMUTE_0X>(RX2, VEC::Negate(RX2)));
            d_A = VEC3::Dot(h_A, VEC::Swizzle<SWIZZLE_W, SWIZZLE_Z, SWIZZLE_Y, SWIZZLE_X>(ARX2));
            d_B = VEC3::Dot(h_B, VEC::Swizzle<SWIZZLE_Y, SWIZZLE_X, SWIZZLE_W, SWIZZLE_Z>(AR0X));
            NoIntersection = VEC::OrInt(NoIntersection,
                VEC::Greater(VEC::Abs(d), VEC::Add(d_A, d_B)));

            // l = a(v) x b(u) = (r20, 0, -r00)
            // d(A) = h(A) dot abs(r20, 0, r00)
            // d(B) = h(B) dot abs(0, r12, r11)
            d = VEC3::Dot(t, VEC::Permute<PERMUTE_0Z, PERMUTE_0W, PERMUTE_1X, PERMUTE_0Y>(RX0, VEC::Negate(RX0)));
            d_A = VEC3::Dot(h_A, VEC::Swizzle<SWIZZLE_Z, SWIZZLE_W, SWIZZLE_X, SWIZZLE_Y>(ARX0));
            d_B = VEC3::Dot(h_B, VEC::Swizzle<SWIZZLE_W, SWIZZLE_Z, SWIZZLE_Y, SWIZZLE_X>(AR1X));
            NoIntersection = VEC::OrInt(NoIntersection,
                VEC::Greater(VEC::Abs(d), VEC::Add(d_A, d_B)));

            // l = a(v) x b(v) = (r21, 0, -r01)
            // d(A) = h(A) dot abs(r21, 0, r01)
            // d(B) = h(B) dot abs(r12, 0, r10)
            d = VEC3::Dot(t, VEC::Permute<PERMUTE_0Z, PERMUTE_0W, PERMUTE_1X, PERMUTE_0Y>(RX1, VEC::Negate(RX1)));
            d_A = VEC3::Dot(h_A, VEC::Swizzle<SWIZZLE_Z, SWIZZLE_W, SWIZZLE_X, SWIZZLE_Y>(ARX1));
            d_B = VEC3::Dot(h_B, VEC::Swizzle<SWIZZLE_Z, SWIZZLE_W, SWIZZLE_X, SWIZZLE_Y>(AR1X));
            NoIntersection = VEC::OrInt(NoIntersection,
                VEC::Greater(VEC::Abs(d), VEC::Add(d_A, d_B)));

            // l = a(v) x b(w) = (r22, 0, -r02)
            // d(A) = h(A) dot abs(r22, 0, r02)
            // d(B) = h(B) dot abs(r11, r10, 0)
            d = VEC3::Dot(t, VEC::Permute<PERMUTE_0Z, PERMUTE_0W, PERMUTE_1X, PERMUTE_0Y>(RX2, VEC::Negate(RX2)));
            d_A = VEC3::Dot(h_A, VEC::Swizzle<SWIZZLE_Z, SWIZZLE_W, SWIZZLE_X, SWIZZLE_Y>(ARX2));
            d_B = VEC3::Dot(h_B, VEC::Swizzle<SWIZZLE_Y, SWIZZLE_X, SWIZZLE_W, SWIZZLE_Z>(AR1X));
            NoIntersection = VEC::OrInt(NoIntersection,
                VEC::Greater(VEC::Abs(d), VEC::Add(d_A, d_B)));

            // l = a(w) x b(u) = (-r10, r00, 0)
            // d(A) = h(A) dot abs(r10, r00, 0)
            // d(B) = h(B) dot abs(0, r22, r21)
            d = VEC3::Dot(t, VEC::Permute<PERMUTE_1Y, PERMUTE_0X, PERMUTE_0W, PERMUTE_0Z>(RX0, VEC::Negate(RX0)));
            d_A = VEC3::Dot(h_A, VEC::Swizzle<SWIZZLE_Y, SWIZZLE_X, SWIZZLE_W, SWIZZLE_Z>(ARX0));
            d_B = VEC3::Dot(h_B, VEC::Swizzle<SWIZZLE_W, SWIZZLE_Z, SWIZZLE_Y, SWIZZLE_X>(AR2X));
            NoIntersection = VEC::OrInt(NoIntersection,
                    VEC::Greater(VEC::Abs(d), VEC::Add(d_A, d_B)));

            // l = a(w) x b(v) = (-r11, r01, 0)
            // d(A) = h(A) dot abs(r11, r01, 0)
            // d(B) = h(B) dot abs(r22, 0, r20)
            d = VEC3::Dot(t, VEC::Permute<PERMUTE_1Y, PERMUTE_0X, PERMUTE_0W, PERMUTE_0Z>(RX1, VEC::Negate(RX1)));
            d_A = VEC3::Dot(h_A, VEC::Swizzle<SWIZZLE_Y, SWIZZLE_X, SWIZZLE_W, SWIZZLE_Z>(ARX1));
            d_B = VEC3::Dot(h_B, VEC::Swizzle<SWIZZLE_Z, SWIZZLE_W, SWIZZLE_X, SWIZZLE_Y>(AR2X));
            NoIntersection = VEC::OrInt(NoIntersection,
                VEC::Greater(VEC::Abs(d), VEC::Add(d_A, d_B)));

            // l = a(w) x b(w) = (-r12, r02, 0)
            // d(A) = h(A) dot abs(r12, r02, 0)
            // d(B) = h(B) dot abs(r21, r20, 0)
            d = VEC3::Dot(t, VEC::Permute<PERMUTE_1Y, PERMUTE_0X, PERMUTE_0W, PERMUTE_0Z>(RX2, VEC::Negate(RX2)));
            d_A = VEC3::Dot(h_A, VEC::Swizzle<SWIZZLE_Y, SWIZZLE_X, SWIZZLE_W, SWIZZLE_Z>(ARX2));
            d_B = VEC3::Dot(h_B, VEC::Swizzle<SWIZZLE_Y, SWIZZLE_X, SWIZZLE_W, SWIZZLE_Z>(AR2X));
            NoIntersection = VEC::OrInt(NoIntersection,
                VEC::Greater(VEC::Abs(d), VEC::Add(d_A, d_B)));

            // No seperating axis found, boxes must intersect.
            return VEC4::NotEqualInt(NoIntersection, VEC::TrueInt()) ? true : false;
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
            VECTOR vCenter = VEC::LoadFloat3(&center);
            VECTOR vOrientation = VEC::LoadFloat4(&orientation);

            // Transform the triangle vertices into the space of the box.
            VECTOR TV0 = VEC3::InverseRotate(VEC::Subtract(V0, vCenter), vOrientation);
            VECTOR TV1 = VEC3::InverseRotate(VEC::Subtract(V1, vCenter), vOrientation);
            VECTOR TV2 = VEC3::InverseRotate(VEC::Subtract(V2, vCenter), vOrientation);

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
            VECTOR vCenter = VEC::LoadFloat3(&center);
            VECTOR vExtents = VEC::LoadFloat3(&extents);
            VECTOR BoxOrientation = VEC::LoadFloat4(&orientation);

#if defined(DEBUG) || defined(_DEBUG)
            assert(QuaternionIsUnit(BoxOrientation));
#endif // DEBUG

            // Set w of the center to one so we can dot4 with a plane.
            vCenter = VEC::Insert<0, 0, 0, 0, 1>(vCenter, VEC::SplatOne());

            // Build the 3x3 rotation matrix that defines the box axes.
            MATRIX R = MAT::RotationQuaternion(BoxOrientation);

            VECTOR Outside, Inside;
            FastIntersectOrientedBoxPlane(vCenter, vExtents, R.r[0], R.r[1], R.r[2], plane, Outside, Inside);

            // If the box is outside any plane it is outside.
            if (VEC4::EqualInt(Outside, VEC::TrueInt()))
                return FRONT;

            // If the box is inside all planes it is inside.
            if (VEC4::EqualInt(Inside, VEC::TrueInt()))
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
            VECTOR vCenter = VEC::LoadFloat3(&center);
            VECTOR vExtents = VEC::LoadFloat3(&extents);
            VECTOR vOrientation = VEC::LoadFloat4(&orientation);

#if defined(DEBUG) || defined(_DEBUG)
            assert(QuaternionIsUnit(vOrientation));
#endif // DEBUG

            // Get the boxes normalized side directions.
            MATRIX R = MAT::RotationQuaternion(vOrientation);

            // Adjust ray origin to be relative to center of the box.
            VECTOR TOrigin = VEC::Subtract(vCenter, origin);

            // Compute the dot product againt each axis of the box.
            VECTOR AxisDotOrigin = VEC3::Dot(R.r[0], TOrigin);
            AxisDotOrigin = VEC::Select(AxisDotOrigin, VEC3::Dot(R.r[1], TOrigin), SelectY);
            AxisDotOrigin = VEC::Select(AxisDotOrigin, VEC3::Dot(R.r[2], TOrigin), SelectZ);

            VECTOR AxisDotDirection = VE3::Dot(R.r[0], direction);
            AxisDotDirection = VEC::Select(AxisDotDirection, VEC3::Dot(R.r[1], direction), SelectY);
            AxisDotDirection = VEC::Select(AxisDotDirection, VEC3::Dot(R.r[2], direction), SelectZ);

            // if (fabs(AxisDotDirection) <= Epsilon) the ray is nearly parallel to the slab.
            VECTOR IsParallel = VEC::LessOrEqual(VEC::Abs(AxisDotDirection), g_RayEpsilon);

            // Test against all three axes simultaneously.
            VECTOR InverseAxisDotDirection = VEC::Reciprocal(AxisDotDirection);
            VECTOR t1 = VEC::Multiply(VEC::Subtract(AxisDotOrigin, vExtents), InverseAxisDotDirection);
            VECTOR t2 = VEC::Multiply(VEC::Add(AxisDotOrigin, vExtents), InverseAxisDotDirection);

            // Compute the max of min(t1,t2) and the min of max(t1,t2) ensuring we don't
            // use the results from any directions parallel to the slab.
            VECTOR t_min = VEC::Select(VEC::Min(t1, t2), g_FltMin, IsParallel);
            VECTOR t_max = VEC::Select(VEC::Max(t1, t2), g_FltMax, IsParallel);

            // t_min.x = maximum( t_min.x, t_min.y, t_min.z );
            // t_max.x = minimum( t_max.x, t_max.y, t_max.z );
            t_min = VEC::Max(t_min, VEC::SplatY(t_min));  // x = max(x,y)
            t_min = VEC::Max(t_min, VEC::SplatZ(t_min));  // x = max(max(x,y),z)
            t_max = VEC::Min(t_max, VEC::SplatY(t_max));  // x = min(x,y)
            t_max = VEC::Min(t_max, VEC::SplatZ(t_max));  // x = min(min(x,y),z)

            // if ( t_min > t_max ) return false;
            VECTOR NoIntersection = VEC::Greater(VEC::SplatX(t_min), VEC::SplatX(t_max));

            // if ( t_max < 0.0f ) return false;
            NoIntersection = VEC::OrInt(NoIntersection, VEC::Less(VEC::SplatX(t_max), VEC::Zero()));

            // if (IsParallel && (-Extents > AxisDotOrigin || Extents < AxisDotOrigin)) return false;
            VECTOR ParallelOverlap = VEC::InBounds(AxisDotOrigin, vExtents);
            NoIntersection = VEC::OrInt(NoIntersection, VEC::AndCInt(IsParallel, ParallelOverlap));

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
            VECTOR vCenter = VEC::LoadFloat3(&center);
            VECTOR vExtents = VEC::LoadFloat3(&extents);
            VECTOR BoxOrientation = VEC::LoadFloat4(&orientation);

#if defined(DEBUG) || defined(_DEBUG)
            assert(QuaternionIsUnit(BoxOrientation));
#endif // DEBUG

            // Set w of the center to one so we can dot4 with a plane.
            vCenter = VEC::Insert<0, 0, 0, 0, 1>(vCenter, VEC::SplatOne());

            // Build the 3x3 rotation matrix that defines the box axes.
            MATRIX R = MAT::RotationQuaternion(BoxOrientation);

            VECTOR Outside, Inside;

            // Test against each plane.
            FastIntersectOrientedBoxPlane(vCenter, vExtents, R.r[0], R.r[1], R.r[2], plane0, Outside, Inside);

            VECTOR AnyOutside = Outside;
            VECTOR AllInside = Inside;

            FastIntersectOrientedBoxPlane(vCenter, vExtents, R.r[0], R.r[1], R.r[2], plane1, Outside, Inside);
            AnyOutside = VEC::OrInt(AnyOutside, Outside);
            AllInside = VEC::AndInt(AllInside, Inside);

            FastIntersectOrientedBoxPlane(vCenter, vExtents, R.r[0], R.r[1], R.r[2], plane2, Outside, Inside);
            AnyOutside = VEC::OrInt(AnyOutside, Outside);
            AllInside = VEC::AndInt(AllInside, Inside);

            FastIntersectOrientedBoxPlane(vCenter, vExtents, R.r[0], R.r[1], R.r[2], plane3, Outside, Inside);
            AnyOutside = VEC::OrInt(AnyOutside, Outside);
            AllInside = VEC::AndInt(AllInside, Inside);

            FastIntersectOrientedBoxPlane(vCenter, vExtents, R.r[0], R.r[1], R.r[2], plane4, Outside, Inside);
            AnyOutside = VEC::OrInt(AnyOutside, Outside);
            AllInside = VEC::AndInt(AllInside, Inside);

            FastIntersectOrientedBoxPlane(vCenter, vExtents, R.r[0], R.r[1], R.r[2], plane5, Outside, Inside);
            AnyOutside = VEC::OrInt(AnyOutside, Outside);
            AllInside = VEC::AndInt(AllInside, Inside);

            // If the box is outside any plane it is outside.
            if (VEC4::EqualInt(AnyOutside, VEC::TrueInt()))
                return DISJOINT;

            // If the box is inside all planes it is inside.
            if (VEC4::EqualInt(AllInside, VEC::TrueInt()))
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

            VECTOR CenterOfMass = VEC::Zero();

            // Compute the center of mass and inertia tensor of the points.
            for (size_t i = 0; i < count; ++i)
            {
                VECTOR Point = VEC::LoadFloat3(reinterpret_cast<const Float3*>(reinterpret_cast<const uint8_t*>(pPoints) + i * stride));

                CenterOfMass = VEC::Add(CenterOfMass, Point);
            }

            CenterOfMass = VEC::Multiply(CenterOfMass, VEC::Reciprocal(VEC::Replicate(float(count))));

            // Compute the inertia tensor of the points around the center of mass.
            // Using the center of mass is not strictly necessary, but will hopefully
            // improve the stability of finding the eigenvectors.
            VECTOR XX_YY_ZZ = VEC::Zero();
            VECTOR XY_XZ_YZ = VEC::Zero();

            for (size_t i = 0; i < count; ++i)
            {
                VECTOR Point = VEC::Subtract(VEC::LoadFloat3(reinterpret_cast<const Float3*>(reinterpret_cast<const uint8_t*>(pPoints) + i * stride)), CenterOfMass);

                XX_YY_ZZ = VEC::Add(XX_YY_ZZ, VEC::Multiply(Point, Point));

                VECTOR XXY = VEC::Swizzle<SWIZZLE_X, SWIZZLE_X, SWIZZLE_Y, SWIZZLE_W>(Point);
                VECTOR YZZ = VEC::Swizzle<SWIZZLE_Y, SWIZZLE_Z, SWIZZLE_Z, SWIZZLE_W>(Point);

                XY_XZ_YZ = VEC::Add(XY_XZ_YZ, VEC::Multiply(XXY, YZZ));
            }

            VECTOR v1, v2, v3;

            // Compute the eigenvectors of the inertia tensor.
            CalculateEigenVectorsFromCovarianceMatrix(VEC::GetX(XX_YY_ZZ), VEC::GetY(XX_YY_ZZ),
                VEC::GetZ(XX_YY_ZZ),
                VEC::GetX(XY_XZ_YZ), VEC::GetY(XY_XZ_YZ),
                VEC::GetZ(XY_XZ_YZ),
                &v1, &v2, &v3);

            // Put them in a matrix.
            MATRIX R;

            R.r[0] = VEC::SetW(v1, 0.f);
            R.r[1] = VEC::SetW(v2, 0.f);
            R.r[2] = VEC::SetW(v3, 0.f);
            R.r[3] = g_IdentityR3.v;

            // Multiply by -1 to convert the matrix into a right handed coordinate
            // system (Det ~= 1) in case the eigenvectors form a left handed
            // coordinate system (Det ~= -1) because XMQuaternionRotationMatrix only
            // works on right handed matrices.
            VECTOR Det = MAT::Determinant(R);

            if (VEC4::Less(Det, VEC::Zero()))
            {
                R.r[0] = VEC::Multiply(R.r[0], g_NegativeOne.v);
                R.r[1] = VEC::Multiply(R.r[1], g_NegativeOne.v);
                R.r[2] = VEC::Multiply(R.r[2], g_NegativeOne.v);
            }

            // Get the rotation quaternion from the matrix.
            VECTOR vOrientation = Quaternion::RotationMatrix(R);

            // Make sure it is normal (in case the vectors are slightly non-orthogonal).
            vOrientation = Quaternion::Normalize(vOrientation);

            // Rebuild the rotation matrix from the quaternion.
            R = MAT::RotationQuaternion(vOrientation);

            // Build the rotation into the rotated space.
            MATRIX InverseR = MAT::Transpose(R);

            // Find the minimum OBB using the eigenvectors as the axes.
            VECTOR vMin, vMax;

            vMin = vMax = VEC3::TransformNormal(VEC::LoadFloat3(pPoints), InverseR);

            for (size_t i = 1; i < count; ++i)
            {
                VECTOR Point = VEC3::TransformNormal(VEC::LoadFloat3(reinterpret_cast<const Float3*>(reinterpret_cast<const uint8_t*>(pPoints) + i * stride)),
                    InverseR);

                vMin = VEC::Min(vMin, Point);
                vMax = VEC::Max(vMax, Point);
            }

            // Rotate the center into world space.
            VECTOR vCenter = VEC::Scale(VEC::Add(vMin, vMax), 0.5f);
            vCenter = VEC3::TransformNormal(vCenter, R);

            // Store center, extents, and orientation.
            VEC::StoreFloat3(&out.center, vCenter);
            VEC::StoreFloat3(&out.extents, VEC::Scale(VEC::Subtract(vMax, vMin), 0.5f));
            VEC::StoreFloat4(&out.orientation, vOrientation);
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
            VECTOR vOrigin = VEC::LoadFloat3(&origin);
            VECTOR vOrientation = VEC::LoadFloat4(&orientation);

#if defined(DEBUG) || defined(_DEBUG)
            assert(QuaternionIsUnit(vOrientation));
#endif // DEBUG

            // Composite the frustum rotation and the transform rotation
            MATRIX nM;
            nM.r[0] = VEC3::Normalize(m.r[0]);
            nM.r[1] = VEC3::Normalize(m.r[1]);
            nM.r[2] = VEC3::Normalize(m.r[2]);
            nM.r[3] = g_IdentityR3;
            VECTOR Rotation = Quaternion::RotationMatrix(nM);
            vOrientation = Quaternion::Multiply(vOrientation, Rotation);

            // Transform the center.
            vOrigin = VEC3::Transform(vOrigin, m);

            // Store the frustum.
            VEC::StoreFloat3(&out.origin, vOrigin);
            VEC::StoreFloat4(&out.orientation, vOrientation);

            // Scale the near and far distances (the slopes remain the same).
            VECTOR dX = VEC3::Dot(m.r[0], m.r[0]);
            VECTOR dY = VEC3::Dot(m.r[1], m.r[1]);
            VECTOR dZ = VEC3::Dot(m.r[2], m.r[2]);

            VECTOR d = VEC::Max(dX, VEC::Max(dY, dZ));
            float Scale = sqrtf(VEC::GetX(d));

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
            VECTOR vOrigin = VEC::LoadFloat3(&origin);
            VECTOR vOrientation = VEC::LoadFloat4(&orientation);

#if defined(DEBUG) || defined(_DEBUG)
            assert(QuaternionIsUnit(vOrientation));
#endif // DEBUG

            // Composite the frustum rotation and the transform rotation.
            vOrientation = Quaternion::Multiply(vOrientation, rotation);

            // Transform the origin.
            vOrigin = VEC::Add(VEC3::Rotate(VEC::Scale(vOrigin, scale), rotation), translation);

            // Store the frustum.
            VEC::StoreFloat3(&out.origin, vOrigin);
            VEC::StoreFloat4(&out.orientation, vOrientation);

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
            VECTOR vOrigin = VEC::LoadFloat3(&origin);
            VECTOR vOrientation = VEC::LoadFloat4(&orientation);

#if defined(DEBUG) || defined(_DEBUG)
            assert(QuaternionIsUnit(vOrientation));
#endif // DEBUG

            // Build the corners of the frustum.
            VECTOR vRightTop = VEC::Set(rightSlope, topSlope, 1.0f, 0.0f);
            VECTOR vRightBottom = VEC::Set(rightSlope, bottomSlope, 1.0f, 0.0f);
            VECTOR vLeftTop = VEC::Set(leftSlope, topSlope, 1.0f, 0.0f);
            VECTOR vLeftBottom = VEC::Set(leftSlope, bottomSlope, 1.0f, 0.0f);
            VECTOR vNear = VEC::ReplicatePtr(&near);
            VECTOR vFar = VEC::ReplicatePtr(&far);

            // Returns 8 corners position of bounding frustum.
            //     Near    Far
            //    0----1  4----5
            //    |    |  |    |
            //    |    |  |    |
            //    3----2  7----6

            VECTOR vCorners[CORNER_COUNT];
            vCorners[0] = VEC::Multiply(vLeftTop, vNear);
            vCorners[1] = VEC::Multiply(vRightTop, vNear);
            vCorners[2] = VEC::Multiply(vRightBottom, vNear);
            vCorners[3] = VEC::Multiply(vLeftBottom, vNear);
            vCorners[4] = VEC::Multiply(vLeftTop, vFar);
            vCorners[5] = VEC::Multiply(vRightTop, vFar);
            vCorners[6] = VEC::Multiply(vRightBottom, vFar);
            vCorners[7] = VEC::Multiply(vLeftBottom, vFar);

            for (size_t i = 0; i < CORNER_COUNT; ++i)
            {
                VECTOR C = VEC::Add(VEC3::Rotate(vCorners[i], vOrientation), vOrigin);
                VEC::StoreFloat3(&corners[i], C);
            }
        }

        _Use_decl_annotations_
        FORCE_INLINE ContainmentType VEC_CALLCONV BoundingFrustum::Contains(A_VECTOR point) const noexcept
        {
            // Build frustum planes.
            VECTOR Planes[6];
            Planes[0] = VEC::Set(0.0f, 0.0f, -1.0f, near);
            Planes[1] = VEC::Set(0.0f, 0.0f, 1.0f, -far);
            Planes[2] = VEC::Set(1.0f, 0.0f, -rightSlope, 0.0f);
            Planes[3] = VEC::Set(-1.0f, 0.0f, leftSlope, 0.0f);
            Planes[4] = VEC::Set(0.0f, 1.0f, -topSlope, 0.0f);
            Planes[5] = VEC::Set(0.0f, -1.0f, bottomSlope, 0.0f);

            // Load origin and orientation.
            VECTOR vOrigin = VEC::LoadFloat3(&origin);
            VECTOR vOrientation = VEC::LoadFloat4(&orientation);

#if defined(DEBUG) || defined(_DEBUG)
            assert(QuaternionIsUnit(vOrientation));
#endif // DEBUG

            // Transform point into local space of frustum.
            VECTOR TPoint = VEC3::InverseRotate(VEC::Subtract(point, vOrigin), vOrientation);

            // Set w to one.
            TPoint = VEC::Insert<0, 0, 0, 0, 1>(TPoint, VEC::SplatOne());

            VECTOR zero = VEC::Zero();
            VECTOR outside = zero;

            // Test point against each plane of the frustum.
            for (size_t i = 0; i < 6; ++i)
            {
                VECTOR dot = VEC4::Dot(TPoint, Planes[i]);
                outside = VEC::OrInt(Outside, VEC::Greater(dot, zero));
            }

            return VEC4::NotEqualInt(outside, VEC::TrueInt()) ? CONTAINS : DISJOINT;
        }

        _Use_decl_annotations_
        FORCE_INLINE ContainmentType VEC_CALLCONV BoundingFrustum::Contains(A_VECTOR V0, A_VECTOR V1, A_VECTOR V2) const noexcept
        {
            // Load origin and orientation of the frustum.
            VECTOR vOrigin = VEC::LoadFloat3(&origin);
            VECTOR vOrientation = VEC::LoadFloat4(&orientation);

            // Create 6 planes (do it inline to encourage use of registers)
            VECTOR nearPlane = VEC::Set(0.0f, 0.0f, -1.0f, near);
            nearPlane = PlaneTransform(nearPlane, vOrientation, vOrigin);
            nearPlane = PlaneNormalize(nearPlane);

            VECTOR farPlane = VEC::Set(0.0f, 0.0f, 1.0f, -far);
            farPlane = PlaneTransform(farPlane, vOrientation, vOrigin);
            farPlane = PlaneNormalize(farPlane);

            VECTOR rightPlane = VEC::Set(1.0f, 0.0f, -rightSlope, 0.0f);
            rightPlane = PlaneTransform(rightPlane, vOrientation, vOrigin);
            rightPlane = PlaneNormalize(rightPlane);

            VECTOR leftPlane = VEC::Set(-1.0f, 0.0f, leftSlope, 0.0f);
            leftPlane = PlaneTransform(leftPlane, vOrientation, vOrigin);
            leftPlane = PlaneNormalize(leftPlane);

            VECTOR topPlane = VEC::Set(0.0f, 1.0f, -topSlope, 0.0f);
            topPlane = PlaneTransform(topPlane, vOrientation, vOrigin);
            topPlane = PlaneNormalize(topPlane);

            VECTOR bottomPlane = VEC::Set(0.0f, -1.0f, bottomSlope, 0.0f);
            bottomPlane = PlaneTransform(bottomPlane, vOrientation, vOrigin);
            bottomPlane = PlaneNormalize(bottomPlane);

            return TriangleTests::ContainedBy(V0, V1, V2, nearPlane, farPlane, rightPlane, leftPlane, topPlane, bottomPlane);
        }

        _Use_decl_annotations_
        FORCE_INLINE ContainmentType BoundingFrustum::Contains(const BoundingSphere& sh) const noexcept
        {
            // Load origin and orientation of the frustum.
            VECTOR vOrigin = VEC::LoadFloat3(&origin);
            VECTOR vOrientation = VEC::LoadFloat4(&orientation);

            // Create 6 planes (do it inline to encourage use of registers)
            VECTOR nearPlane = VEC::Set(0.0f, 0.0f, -1.0f, near);
            nearPlane = PlaneTransform(nearPlane, vOrientation, vOrigin);
            nearPlane = PlaneNormalize(nearPlane);

            VECTOR farPlane = VEC::Set(0.0f, 0.0f, 1.0f, -far);
            farPlane = PlaneTransform(farPlane, vOrientation, vOrigin);
            farPlane = PlaneNormalize(farPlane);

            VECTOR rightPlane = VEC::Set(1.0f, 0.0f, -rightSlope, 0.0f);
            rightPlane = PlaneTransform(rightPlane, vOrientation, vOrigin);
            rightPlane = PlaneNormalize(rightPlane);

            VECTOR leftPlane = VEC::Set(-1.0f, 0.0f, leftSlope, 0.0f);
            leftPlane = PlaneTransform(leftPlane, vOrientation, vOrigin);
            leftPlane = PlaneNormalize(leftPlane);

            VECTOR topPlane = VEC::Set(0.0f, 1.0f, -topSlope, 0.0f);
            topPlane = PlaneTransform(topPlane, vOrientation, vOrigin);
            topPlane = PlaneNormalize(topPlane);

            VECTOR bottomPlane = VEC::Set(0.0f, -1.0f, bottomSlope, 0.0f);
            bottomPlane = PlaneTransform(bottomPlane, vOrientation, vOrigin);
            bottomPlane = PlaneNormalize(bottomPlane);

            return sh.ContainedBy(nearPlane, farPlane, rightPlane, leftPlane, topPlane, bottomPlane);
        }

        _Use_decl_annotations_
        FORCE_INLINE ContainmentType BoundingFrustum::Contains(const BoundingBox& box) const noexcept
        {
            // Load origin and orientation of the frustum.
            VECTOR vOrigin = VEC::LoadFloat3(&origin);
            VECTOR vOrientation = VEC::LoadFloat4(&orientation);

            // Create 6 planes (do it inline to encourage use of registers)
            VECTOR nearPlane = VEC::Set(0.0f, 0.0f, -1.0f, near);
            nearPlane = PlaneTransform(nearPlane, vOrientation, vOrigin);
            nearPlane = PlaneNormalize(nearPlane);

            VECTOR farPlane = VEC::Set(0.0f, 0.0f, 1.0f, -far);
            farPlane = PlaneTransform(farPlane, vOrientation, vOrigin);
            farPlane = PlaneNormalize(farPlane);

            VECTOR rightPlane = VEC::Set(1.0f, 0.0f, -rightSlope, 0.0f);
            rightPlane = PlaneTransform(rightPlane, vOrientation, vOrigin);
            rightPlane = PlaneNormalize(rightPlane);

            VECTOR leftPlane = VEC::Set(-1.0f, 0.0f, leftSlope, 0.0f);
            leftPlane = PlaneTransform(leftPlane, vOrientation, vOrigin);
            leftPlane = PlaneNormalize(leftPlane);

            VECTOR topPlane = VEC::Set(0.0f, 1.0f, -topSlope, 0.0f);
            topPlane = PlaneTransform(topPlane, vOrientation, vOrigin);
            topPlane = PlaneNormalize(topPlane);

            VECTOR bottomPlane = VEC::Set(0.0f, -1.0f, bottomSlope, 0.0f);
            bottomPlane = PlaneTransform(bottomPlane, vOrientation, vOrigin);
            bottomPlane = PlaneNormalize(bottomPlane);

            return box.ContainedBy(nearPlane, farPlane, rightPlane, leftPlane, topPlane, bottomPlane);
        }

        _Use_decl_annotations_
        FORCE_INLINE ContainmentType BoundingFrustum::Contains(const BoundingOrientedBox& box) const noexcept
        {
            // Load origin and orientation of the frustum.
            VECTOR vOrigin = VEC::LoadFloat3(&origin);
            VECTOR vOrientation = VEC::LoadFloat4(&orientation);

            // Create 6 planes (do it inline to encourage use of registers)
            VECTOR nearPlane = VEC::Set(0.0f, 0.0f, -1.0f, near);
            nearPlane = PlaneTransform(nearPlane, vOrientation, vOrigin);
            nearPlane = PlaneNormalize(nearPlane);

            VECTOR farPlane = VEC::Set(0.0f, 0.0f, 1.0f, -far);
            farPlane = PlaneTransform(farPlane, vOrientation, vOrigin);
            farPlane = PlaneNormalize(farPlane);

            VECTOR rightPlane = VEC::Set(1.0f, 0.0f, -rightSlope, 0.0f);
            rightPlane = PlaneTransform(rightPlane, vOrientation, vOrigin);
            rightPlane = PlaneNormalize(rightPlane);

            VECTOR leftPlane = VEC::Set(-1.0f, 0.0f, leftSlope, 0.0f);
            leftPlane = PlaneTransform(leftPlane, vOrientation, vOrigin);
            leftPlane = PlaneNormalize(leftPlane);

            VECTOR topPlane = VEC::Set(0.0f, 1.0f, -topSlope, 0.0f);
            topPlane = PlaneTransform(topPlane, vOrientation, vOrigin);
            topPlane = PlaneNormalize(topPlane);

            VECTOR bottomPlane = VEC::Set(0.0f, -1.0f, bottomSlope, 0.0f);
            bottomPlane = PlaneTransform(bottomPlane, vOrientation, vOrigin);
            bottomPlane = PlaneNormalize(bottomPlane);

            return box.ContainedBy(nearPlane, farPlane, rightPlane, leftPlane, topPlane, bottomPlane);
        }

        _Use_decl_annotations_
        FORCE_INLINE ContainmentType BoundingFrustum::Contains(const BoundingFrustum& fr) const noexcept
        {
            // Load origin and orientation of the frustum.
            VECTOR vOrigin = VEC::LoadFloat3(&origin);
            VECTOR vOrientation = VEC::LoadFloat4(&orientation);

            // Create 6 planes (do it inline to encourage use of registers)
            VECTOR nearPlane = VEC::Set(0.0f, 0.0f, -1.0f, near);
            nearPlane = PlaneTransform(nearPlane, vOrientation, vOrigin);
            nearPlane = PlaneNormalize(nearPlane);

            VECTOR farPlane = VEC::Set(0.0f, 0.0f, 1.0f, -far);
            farPlane = PlaneTransform(farPlane, vOrientation, vOrigin);
            farPlane = PlaneNormalize(farPlane);

            VECTOR rightPlane = VEC::Set(1.0f, 0.0f, -rightSlope, 0.0f);
            rightPlane = PlaneTransform(rightPlane, vOrientation, vOrigin);
            rightPlane = PlaneNormalize(rightPlane);

            VECTOR leftPlane = VEC::Set(-1.0f, 0.0f, leftSlope, 0.0f);
            leftPlane = PlaneTransform(leftPlane, vOrientation, vOrigin);
            leftPlane = PlaneNormalize(leftPlane);

            VECTOR topPlane = VEC::Set(0.0f, 1.0f, -topSlope, 0.0f);
            topPlane = PlaneTransform(topPlane, vOrientation, vOrigin);
            topPlane = PlaneNormalize(topPlane);

            VECTOR bottomPlane = VEC::Set(0.0f, -1.0f, bottomSlope, 0.0f);
            bottomPlane = PlaneTransform(bottomPlane, vOrientation, vOrigin);
            bottomPlane = PlaneNormalize(bottomPlane);

            return fr.ContainedBy(nearPlane, farPlane, rightPlane, leftPlane, topPlane, bottomPlane);
        }

        _Use_decl_annotations_
        FORCE_INLINE bool BoundingFrustum::Intersects(const BoundingSphere& sh) const noexcept
        {
            VECTOR zero = VEC::Zero();

            // Build the frustum planes.
            VECTOR planes[6];
            planes[0] = VEC::Set(0.0f, 0.0f, -1.0f, near);
            planes[1] = VEC::Set(0.0f, 0.0f, 1.0f, -far);
            planes[2] = VEC::Set(1.0f, 0.0f, -rightSlope, 0.0f);
            planes[3] = VEC::Set(-1.0f, 0.0f, leftSlope, 0.0f);
            planes[4] = VEC::Set(0.0f, 1.0f, -topSlope, 0.0f);
            planes[5] = VEC::Set(0.0f, -1.0f, bottomSlope, 0.0f);

            // Normalize the planes so we can compare to the sphere radius.
            planes[2] = VEC3::Normalize(planes[2]);
            planes[3] = VEC3::Normalize(planes[3]);
            planes[4] = VEC3::Normalize(planes[4]);
            planes[5] = VEC3::Normalize(planes[5]);

            // Load origin and orientation of the frustum.
            VECTOR vOrigin = VEC::LoadFloat3(&origin);
            VECTOR vOrientation = VEC::LoadFloat4(&orientation);

#if defined(DEBUG) || defined(_DEBUG)
            assert(QuaternionIsUnit(vOrientation));
#endif // DEBUG

            // Load the sphere.
            VECTOR vCenter = VEC::LoadFloat3(&sh.center);
            VECTOR vRadius = VEC::ReplicatePtr(&sh.radius);

            // Transform the center of the sphere into the local space of frustum.
            vCenter = VEC3::InverseRotate(VEC::Subtract(vCenter, vOrigin), vOrientation);

            // Set w of the center to one so we can dot4 with the plane.
            vCenter = VEC::Insert<0, 0, 0, 0, 1>(vCenter, VEC::SplatOne());

            // Check against each plane of the frustum.
            VECTOR outside = VEC::FalseInt();
            VECTOR insideAll = VEC::TrueInt();
            VECTOR centerInsideAll = VEC::TrueInt();

            VECTOR dist[6];

            for (size_t i = 0; i < 6; ++i)
            {
                dist[i] = VEC4::Dot(vCenter, planes[i]);

                // Outside the plane?
                outside = VEC::OrInt(outside, VEC::Greater(dist[i], vRadius));

                // Fully inside the plane?
                insideAll = VEC::AndInt(insideAll, VEC::LessOrEqual(dist[i], VEC::Negate(vRadius)));

                // Check if the center is inside the plane.
                centerInsideAll = VEC::AndInt(centerInsideAll, VEC::LessOrEqual(dist[i], zero));
            }

            // If the sphere is outside any of the planes it is outside.
            if (VEC4::EqualInt(outside, VEC::TrueInt()))
                return false;

            // If the sphere is inside all planes it is fully inside.
            if (VEC4::EqualInt(insideAll, VEC::TrueInt()))
                return true;

            // If the center of the sphere is inside all planes and the sphere intersects
            // one or more planes then it must intersect.
            if (VEC4::EqualInt(centerInsideAll, VEC::TrueInt()))
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

            VECTOR intersects = VEC::FalseInt();

            // Check to see if the nearest feature is one of the planes.
            for (size_t i = 0; i < 6; ++i)
            {
                // Find the nearest point on the plane to the center of the sphere.
                VECTOR point = VEC::NegativeMultiplySubtract(planes[i], dist[i], vCenter);

                // Set w of the point to one.
                point = VEC::Insert<0, 0, 0, 0, 1>(point, VEC::SplatOne());

                // If the point is inside the face (inside the adjacent planes) then
                // this plane is the nearest feature.
                VECTOR insideFace = VEC::TrueInt();

                for (size_t j = 0; j < 4; j++)
                {
                    size_t plane_index = adjacent_faces[i][j];

                    insideFace = VEC::AndInt(insideFace,
                        VEC::LessOrEqual(VEC4::Dot(point, planes[plane_index]), zero));
                }

                // Since we have already checked distance from the plane we know that the
                // sphere must intersect if this plane is the nearest feature.
                intersects = VEC::OrInt(intersects,
                    VEC::AndInt(VEC::Greater(dist[i], zero), insideFace));
            }

            if (VEC4::EqualInt(intersects, VEC::TrueInt()))
                return true;

            // Build the corners of the frustum.
            VECTOR vRightTop = VEC::Set(rightSlope, topSlope, 1.0f, 0.0f);
            VECTOR vRightBottom = VEC::Set(rightSlope, bottomSlope, 1.0f, 0.0f);
            VECTOR vLeftTop = VEC::Set(leftSlope, topSlope, 1.0f, 0.0f);
            VECTOR vLeftBottom = VEC::Set(leftSlope, bottomSlope, 1.0f, 0.0f);
            VECTOR vNear = VEC::ReplicatePtr(&near);
            VECTOR vFar = VEC::ReplicatePtr(&far);

            VECTOR corners[CORNER_COUNT];
            corners[0] = VEC::Multiply(vRightTop, vNear);
            corners[1] = VEC::Multiply(vRightBottom, vNear);
            corners[2] = VEC::Multiply(vLeftTop, vNear);
            corners[3] = VEC::Multiply(vLeftBottom, vNear);
            corners[4] = VEC::Multiply(vRightTop, vFar);
            corners[5] = VEC::Multiply(vRightBottom, vFar);
            corners[6] = VEC::Multiply(vLeftTop, vFar);
            corners[7] = VEC::Multiply(vLeftBottom, vFar);

            // The Edges are:
            static const size_t edges[12][2] =
            {
                { 0, 1 }, { 2, 3 }, { 0, 2 }, { 1, 3 },    // Near plane
                { 4, 5 }, { 6, 7 }, { 4, 6 }, { 5, 7 },    // Far plane
                { 0, 4 }, { 1, 5 }, { 2, 6 }, { 3, 7 },
            }; // Near to far

            VECTOR radiusSq = VEC::Multiply(vRadius, vRadius);

            // Check to see if the nearest feature is one of the edges (or corners).
            for (size_t i = 0; i < 12; ++i)
            {
                size_t ei0 = edges[i][0];
                size_t ei1 = edges[i][1];

                // Find the nearest point on the edge to the center of the sphere.
                // The corners of the frustum are included as the endpoints of the edges.
                VECTOR point = PointOnLineSegmentNearestPoint(corners[ei0], corners[ei1], vCenter);

                VECTOR delta = VEC::Subtract(vCenter, point);

                VECTOR distSq = VEC3::Dot(delta, delta);

                // If the distance to the center of the sphere to the point is less than
                // the radius of the sphere then it must intersect.
                intersects = VEC::OrInt(intersects, VEC::LessOrEqual(distSq, radiusSq));
            }

            if (VEC4::EqualInt(intersects, VEC::TrueInt()))
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

            VECTOR zero = VEC::Zero();

            // Build the frustum planes.
            VECTOR planes[6];
            planes[0] = VEC::Set(0.0f, 0.0f, -1.0f, near);
            planes[1] = VEC::Set(0.0f, 0.0f, 1.0f, -far);
            planes[2] = VEC::Set(1.0f, 0.0f, -rightSlope, 0.0f);
            planes[3] = VEC::Set(-1.0f, 0.0f, leftSlope, 0.0f);
            planes[4] = VEC::Set(0.0f, 1.0f, -topSlope, 0.0f);
            planes[5] = VEC::Set(0.0f, -1.0f, bottomSlope, 0.0f);

            // Load origin and orientation of the frustum.
            VECTOR vOrigin = VEC::LoadFloat3(&origin);
            VECTOR frustumOrientation = VEC::LoadFloat4(&orientation);

#if defined(DEBUG) || defined(_DEBUG)
            assert(QuaternionIsUnit(frustumOrientation));
#endif // DEBUG

            // Load the box.
            VECTOR center = VEC::LoadFloat3(&box.center);
            VECTOR extents = VEC::LoadFloat3(&box.extents);
            VECTOR boxOrientation = VEC::LoadFloat4(&box.orientation);

#if defined(DEBUG) || defined(_DEBUG)
            assert(QuaternionIsUnit(boxOrientation));
#endif // DEBUG

            // Transform the oriented box into the space of the frustum in order to
            // minimize the number of transforms we have to do.
            center = VEC3::InverseRotate(VEC::Subtract(center, vOrigin), frustumOrientation);
            boxOrientation = Quaternion::Multiply(boxOrientation, Quaternion::Conjugate(frustumOrientation));

            // Set w of the center to one so we can dot4 with the plane.
            center = VEC::Insert<0, 0, 0, 0, 1>(center, VEC::SplatOne());

            // Build the 3x3 rotation matrix that defines the box axes.
            MATRIX R = MAT::RotationQuaternion(boxOrientation);

            // Check against each plane of the frustum.
            VECTOR outside = VEC::FalseInt();
            VECTOR insideAll = VEC::TrueInt();
            VECTOR centerInsideAll = VEC::TrueInt();

            for (size_t i = 0; i < 6; ++i)
            {
                // Compute the distance to the center of the box.
                VECTOR dist = VEC4::Dot(center, planes[i]);

                // Project the axes of the box onto the normal of the plane.  Half the
                // length of the projection (sometime called the "radius") is equal to
                // h(u) * abs(n dot b(u))) + h(v) * abs(n dot b(v)) + h(w) * abs(n dot b(w))
                // where h(i) are extents of the box, n is the plane normal, and b(i) are the
                // axes of the box.
                VECTOR radius = VEC3::Dot(planes[i], R.r[0]);
                radius = VEC::Select(radius, VEC3::Dot(planes[i], R.r[1]), selectY);
                radius = VEC::Select(radius, VEC3::Dot(planes[i], R.r[2]), selectZ);
                radius = VEC3::Dot(extents, VEC::Abs(radius));

                // Outside the plane?
                outside = VEC::OrInt(outside, VEC::Greater(dist, radius));

                // Fully inside the plane?
                insideAll = VEC::AndInt(insideAll, VEC::LessOrEqual(dist, VEC::Negate(radius)));

                // Check if the center is inside the plane.
                centerInsideAll = VEC::AndInt(centerInsideAll, VEC::LessOrEqual(dist, zero));
            }

            // If the box is outside any of the planes it is outside.
            if (VEC4::EqualInt(outside, VEC::TrueInt()))
                return false;

            // If the box is inside all planes it is fully inside.
            if (VEC4::EqualInt(insideAll, VEC::TrueInt()))
                return true;

            // If the center of the box is inside all planes and the box intersects
            // one or more planes then it must intersect.
            if (VEC4::EqualInt(centerInsideAll, VEC::TrueInt()))
                return true;

            // Build the corners of the frustum.
            VECTOR vRightTop = VEC::Set(rightSlope, topSlope, 1.0f, 0.0f);
            VECTOR vRightBottom = VEC::Set(rightSlope, bottomSlope, 1.0f, 0.0f);
            VECTOR vLeftTop = VEC::Set(leftSlope, topSlope, 1.0f, 0.0f);
            VECTOR vLeftBottom = VEC::Set(leftSlope, bottomSlope, 1.0f, 0.0f);
            VECTOR vNear = VEC::ReplicatePtr(&near);
            VECTOR vFar = VEC::ReplicatePtr(&far);

            VECTOR corners[CORNER_COUNT];
            corners[0] = VEC::Multiply(vRightTop, vNear);
            corners[1] = VEC::Multiply(vRightBottom, vNear);
            corners[2] = VEC::Multiply(vLeftTop, vNear);
            corners[3] = VEC::Multiply(vLeftBottom, vNear);
            corners[4] = VEC::Multiply(vRightTop, vFar);
            corners[5] = VEC::Multiply(vRightBottom, vFar);
            corners[6] = VEC::Multiply(vLeftTop, vFar);
            corners[7] = VEC::Multiply(vLeftBottom, vFar);

            // Test against box axes (3)
            {
                // Find the min/max values of the projection of the frustum onto each axis.
                VECTOR frustumMin, frustumMax;

                frustumMin = VEC3::Dot(corners[0], R.r[0]);
                frustumMin = VEC::Select(frustumMin, VEC3::Dot(corners[0], R.r[1]), selectY);
                frustumMin = VEC::Select(frustumMin, VEC3::Dot(corners[0], R.r[2]), selectZ);
                frustumMax = frustumMin;

                for (size_t i = 1; i < BoundingOrientedBox::CORNER_COUNT; ++i)
                {
                    VECTOR temp = VEC3::Dot(corners[i], R.r[0]);
                    temp = VEC::Select(temp, VEC3::Dot(corners[i], R.r[1]), selectY);
                    temp = VEC::Select(temp, VEC3::Dot(corners[i], R.r[2]), selectZ);

                    frustumMin = VEC::Min(frustumMin, temp);
                    frustumMax = VEC::Max(frustumMax, temp);
                }

                // Project the center of the box onto the axes.
                VECTOR boxDist = VEC3::Dot(center, R.r[0]);
                boxDist = VEC::Select(boxDist, VEC3::Dot(center, R.r[1]), selectY);
                boxDist = VEC::Select(boxDist, VEC3::Dot(center, R.r[2]), selectZ);

                // The projection of the box onto the axis is just its Center and Extents.
                // if (min > box_max || max < box_min) reject;
                VECTOR result = VEC::OrInt(VEC::Greater(frustumMin, VEC::Add(boxDist, extents)),
                    VEC::Less(frustumMax, VEC::Subtract(boxDist, extents)));

                if (Vector3AnyTrue(result))
                    return false;
            }

            // Test against edge/edge axes (3*6).
            VECTOR frustumEdgeAxis[6];

            frustumEdgeAxis[0] = vRightTop;
            frustumEdgeAxis[1] = vRightBottom;
            frustumEdgeAxis[2] = vLeftTop;
            frustumEdgeAxis[3] = vLeftBottom;
            frustumEdgeAxis[4] = VEC::Subtract(vRightTop, vLeftTop);
            frustumEdgeAxis[5] = VEC::Subtract(vLeftBottom, vLeftTop);

            for (size_t i = 0; i < 3; ++i)
            {
                for (size_t j = 0; j < 6; j++)
                {
                    // Compute the axis we are going to test.
                    VECTOR axis = VEC3::Cross(R.r[i], frustumEdgeAxis[j]);

                    // Find the min/max values of the projection of the frustum onto the axis.
                    VECTOR frustumMin, frustumMax;

                    frustumMin = frustumMax = VEC3::Dot(axis, corners[0]);

                    for (size_t k = 1; k < CORNER_COUNT; k++)
                    {
                        VECTOR temp = VEC3::Dot(axis, corners[k]);
                        frustumMin = VEC::Min(frustumMin, temp);
                        frustumMax = VEC::Max(frustumMax, temp);
                    }

                    // Project the center of the box onto the axis.
                    VECTOR dist = VEC3::Dot(center, axis);

                    // Project the axes of the box onto the axis to find the "radius" of the box.
                    VECTOR radius = VEC3::Dot(axis, R.r[0]);
                    radius = VEC::Select(radius, VEC3::Dot(axis, R.r[1]), selectY);
                    radius = VEC::Select(radius, VEC3::Dot(axis, R.r[2]), selectZ);
                    radius = VEC3::Dot(extents, VEC::Abs(radius));

                    // if (center > max + radius || center < min - radius) reject;
                    outside = VEC::OrInt(outside, VEC::Greater(dist, VEC::Add(frustumMax, radius)));
                    outside = VEC::OrInt(outside, VEC::Less(dist, VEC::Subtract(frustumMin, radius)));
                }
            }

            if (VEC4::EqualInt(outside, VEC::TrueInt()))
                return false;

            // If we did not find a separating plane then the box must intersect the frustum.
            return true;
        }

        _Use_decl_annotations_
        FORCE_INLINE bool BoundingFrustum::Intersects(const BoundingFrustum& fr) const noexcept
        {
            // Load origin and orientation of frustum B.
            VECTOR originB = VEC::LoadFloat3(&origin);
            VECTOR orientationB = VEC::LoadFloat4(&orientation);

#if defined(DEBUG) || defined(_DEBUG)
            assert(QuaternionIsUnit(orientationB));
#endif // DEBUG

            // Build the planes of frustum B.
            VECTOR axisB[6];
            axisB[0] = VEC::Set(0.0f, 0.0f, -1.0f, 0.0f);
            axisB[1] = VEC::Set(0.0f, 0.0f, 1.0f, 0.0f);
            axisB[2] = VEC::Set(1.0f, 0.0f, -rightSlope, 0.0f);
            axisB[3] = VEC::Set(-1.0f, 0.0f, leftSlope, 0.0f);
            axisB[4] = VEC::Set(0.0f, 1.0f, -topSlope, 0.0f);
            axisB[5] = VEC::Set(0.0f, -1.0f, bottomSlope, 0.0f);

            VECTOR planeDistB[6];
            planeDistB[0] = VEC::Negate(VEC::ReplicatePtr(&near));
            planeDistB[1] = VEC::ReplicatePtr(&far);
            planeDistB[2] = VEC::Zero();
            planeDistB[3] = VEC::Zero();
            planeDistB[4] = VEC::Zero();
            planeDistB[5] = VEC::Zero();

            // Load origin and orientation of frustum A.
            VECTOR originA = VEC::LoadFloat3(&fr.origin);
            VECTOR orientationA = VEC::LoadFloat4(&fr.orientation);

#if defined(DEBUG) || defined(_DEBUG)
            assert(QuaternionIsUnit(orientationA));
#endif // DEBUG

            // Transform frustum A into the space of the frustum B in order to
            // minimize the number of transforms we have to do.
            originA = VEC3::InverseRotate(VEC::Subtract(originA, originB), orientationB);
            orientationA = Quaternion::Multiply(orientationA, Quaternion::Conjugate(orientationB));

            // Build the corners of frustum A (in the local space of B).
            VECTOR rightTopA = VEC::Set(fr.rightSlope, fr.topSlope, 1.0f, 0.0f);
            VECTOR rightBottomA = VEC::Set(fr.rightSlope, fr.bottomSlope, 1.0f, 0.0f);
            VECTOR leftTopA = VEC::Set(fr.leftSlope, fr.topSlope, 1.0f, 0.0f);
            VECTOR leftBottomA = VEC::Set(fr.leftSlope, fr.bottomSlope, 1.0f, 0.0f);
            VECTOR nearA = VEC::ReplicatePtr(&fr.near);
            VECTOR farA = VEC::ReplicatePtr(&fr.far);

            rightTopA = VEC3::Rotate(rightTopA, orientationA);
            rightBottomA = VEC3::Rotate(rightBottomA, orientationA);
            leftTopA = VEC3::Rotate(leftTopA, orientationA);
            leftBottomA = VEC3::Rotate(leftBottomA, orientationA);

            VECTOR cornersA[CORNER_COUNT];
            cornersA[0] = VEC::MultiplyAdd(rightTopA, nearA, originA);
            cornersA[1] = VEC::MultiplyAdd(rightBottomA, nearA, originA);
            cornersA[2] = VEC::MultiplyAdd(leftTopA, nearA, originA);
            cornersA[3] = VEC::MultiplyAdd(leftBottomA, nearA, originA);
            cornersA[4] = VEC::MultiplyAdd(rightTopA, farA, originA);
            cornersA[5] = VEC::MultiplyAdd(rightBottomA, farA, originA);
            cornersA[6] = VEC::MultiplyAdd(leftTopA, farA, originA);
            cornersA[7] = VEC::MultiplyAdd(leftBottomA, farA, originA);

            // Check frustum A against each plane of frustum B.
            VECTOR outside = VEC::FalseInt();
            VECTOR insideAll = VEC::TrueInt();

            for (size_t i = 0; i < 6; ++i)
            {
                // Find the min/max projection of the frustum onto the plane normal.
                VECTOR min, max;

                min = max = VEC3::Dot(axisB[i], cornersA[0]);

                for (size_t j = 1; j < CORNER_COUNT; j++)
                {
                    VECTOR temp = VEC3::Dot(axisB[i], cornersA[j]);
                    min = VEC::Min(min, temp);
                    max = VEC::Max(max, temp);
                }

                // Outside the plane?
                outside = VEC::OrInt(outside, VEC::Greater(min, planeDistB[i]));

                // Fully inside the plane?
                insideAll = VEC::AndInt(insideAll, VEC::LessOrEqual(max, planeDistB[i]));
            }

            // If the frustum A is outside any of the planes of frustum B it is outside.
            if (VEC4::EqualInt(outside, VEC::TrueInt()))
                return false;

            // If frustum A is inside all planes of frustum B it is fully inside.
            if (VEC4::EqualInt(insideAll, VEC::TrueInt()))
                return true;

            // Build the corners of frustum B.
            VECTOR rightTopB = VEC::Set(rightSlope, topSlope, 1.0f, 0.0f);
            VECTOR rightBottomB = VEC::Set(rightSlope, bottomSlope, 1.0f, 0.0f);
            VECTOR leftTopB = VEC::Set(leftSlope, topSlope, 1.0f, 0.0f);
            VECTOR leftBottomB = VEC::Set(leftSlope, bottomSlope, 1.0f, 0.0f);
            VECTOR nearB = VEC::ReplicatePtr(&near);
            VECTOR farB = VEC::ReplicatePtr(&far);

            VECTOR cornersB[BoundingFrustum::CORNER_COUNT];
            cornersB[0] = VEC::Multiply(rightTopB, nearB);
            cornersB[1] = VEC::Multiply(rightBottomB, nearB);
            cornersB[2] = VEC::Multiply(leftTopB, nearB);
            cornersB[3] = VEC::Multiply(leftBottomB, nearB);
            cornersB[4] = VEC::Multiply(rightTopB, farB);
            cornersB[5] = VEC::Multiply(rightBottomB, farB);
            cornersB[6] = VEC::Multiply(leftTopB, farB);
            cornersB[7] = VEC::Multiply(leftBottomB, farB);

            // Build the planes of frustum A (in the local space of B).
            VECTOR axisA[6];
            VECTOR planeDistA[6];

            axisA[0] = VEC::Set(0.0f, 0.0f, -1.0f, 0.0f);
            axisA[1] = VEC::Set(0.0f, 0.0f, 1.0f, 0.0f);
            axisA[2] = VEC::Set(1.0f, 0.0f, -fr.rightSlope, 0.0f);
            axisA[3] = VEC::Set(-1.0f, 0.0f, fr.leftSlope, 0.0f);
            axisA[4] = VEC::Set(0.0f, 1.0f, -fr.topSlope, 0.0f);
            axisA[5] = VEC::Set(0.0f, -1.0f, fr.bottomSlope, 0.0f);

            axisA[0] = VEC3::Rotate(axisA[0], orientationA);
            axisA[1] = VEC::Negate(axisA[0]);
            axisA[2] = VEC3::Rotate(axisA[2], orientationA);
            axisA[3] = VEC3::Rotate(axisA[3], orientationA);
            axisA[4] = VEC3::Rotate(axisA[4], orientationA);
            axisA[5] = VEC3::Rotate(axisA[5], orientationA);

            planeDistA[0] = VEC3::Dot(axisA[0], cornersA[0]);  // Re-use corner on near plane.
            planeDistA[1] = VEC3::Dot(axisA[1], cornersA[4]);  // Re-use corner on far plane.
            planeDistA[2] = VEC3::Dot(axisA[2], originA);
            planeDistA[3] = VEC3::Dot(axisA[3], originA);
            planeDistA[4] = VEC3::Dot(axisA[4], originA);
            planeDistA[5] = VEC3::Dot(axisA[5], originA);

            // Check each axis of frustum A for a seperating plane (5).
            for (size_t i = 0; i < 6; ++i)
            {
                // Find the minimum projection of the frustum onto the plane normal.
                VECTOR min;

                min = VEC3::Dot(axisA[i], cornersB[0]);

                for (size_t j = 1; j < CORNER_COUNT; j++)
                {
                    VECTOR temp = VEC3::Dot(axisA[i], cornersB[j]);
                    min = VEC::Min(min, temp);
                }

                // Outside the plane?
                outside = VEC::OrInt(outside, VEC::Greater(min, planeDistA[i]));
            }

            // If the frustum B is outside any of the planes of frustum A it is outside.
            if (VEC4::EqualInt(outside, VEC::TrueInt()))
                return false;

            // Check edge/edge axes (6 * 6).
            VECTOR frustumEdgeAxisA[6];
            frustumEdgeAxisA[0] = rightTopA;
            frustumEdgeAxisA[1] = rightBottomA;
            frustumEdgeAxisA[2] = leftTopA;
            frustumEdgeAxisA[3] = leftBottomA;
            frustumEdgeAxisA[4] = VEC::Subtract(rightTopA, leftTopA);
            frustumEdgeAxisA[5] = VEC::Subtract(leftBottomA, leftTopA);

            VECTOR frustumEdgeAxisB[6];
            frustumEdgeAxisB[0] = rightTopB;
            frustumEdgeAxisB[1] = rightBottomB;
            frustumEdgeAxisB[2] = leftTopB;
            frustumEdgeAxisB[3] = leftBottomB;
            frustumEdgeAxisB[4] = VEC::Subtract(rightTopB, leftTopB);
            frustumEdgeAxisB[5] = VEC::Subtract(leftBottomB, leftTopB);

            for (size_t i = 0; i < 6; ++i)
            {
                for (size_t j = 0; j < 6; j++)
                {
                    // Compute the axis we are going to test.
                    VECTOR axis = VEC3::Cross(frustumEdgeAxisA[i], frustumEdgeAxisB[j]);

                    // Find the min/max values of the projection of both frustums onto the axis.
                    VECTOR minA, maxA;
                    VECTOR minB, maxB;

                    minA = maxA = VEC3::Dot(axis, cornersA[0]);
                    minB = maxB = VEC3::Dot(axis, cornersB[0]);

                    for (size_t k = 1; k < CORNER_COUNT; k++)
                    {
                        VECTOR tempA = VEC3::Dot(axis, cornersA[k]);
                        minA = VEC::Min(minA, tempA);
                        maxA = VEC::Max(maxA, tempA);

                        VECTOR tempB = VEC3::Dot(axis, cornersB[k]);
                        minB = VEC::Min(minB, tempB);
                        maxB = VE::Max(maxB, tempB);
                    }

                    // if (MinA > MaxB || MinB > MaxA) reject
                    outside = VEC::OrInt(outside, VEC::Greater(minA, maxB));
                    outside = VEC::OrInt(outside, VEC::Greater(minB, maxA));
                }
            }

            // If there is a seperating plane, then the frustums do not intersect.
            if (VEC4::EqualInt(outside, VEC::TrueInt()))
                return false;

            // If we did not find a separating plane then the frustums intersect.
            return true;
        }

        _Use_decl_annotations_
        FORCE_INLINE bool VEC_CALLCONV BoundingFrustum::Intersects(A_VECTOR V0, A_VECTOR V1, A_VECTOR V2) const noexcept
        {
            // Build the frustum planes (NOTE: D is negated from the usual).
            VECTOR planes[6];
            planes[0] = VEC::Set(0.0f, 0.0f, -1.0f, -near);
            planes[1] = VEC::Set(0.0f, 0.0f, 1.0f, far);
            planes[2] = VEC::Set(1.0f, 0.0f, -rightSlope, 0.0f);
            planes[3] = VEC::Set(-1.0f, 0.0f, leftSlope, 0.0f);
            planes[4] = VEC::Set(0.0f, 1.0f, -topSlope, 0.0f);
            planes[5] = VEC::Set(0.0f, -1.0f, bottomSlope, 0.0f);

            // Load origin and orientation of the frustum.
            VECTOR vOrigin = VEC::LoadFloat3(&origin);
            VECTOR vOrientation = VEC::LoadFloat4(&orientation);

#if defined(DEBUG) || defined(_DEBUG)
            assert(QuaternionIsUnit(vOrientation));
#endif // DEBUG

            // Transform triangle into the local space of frustum.
            VECTOR TV0 = VEC3::InverseRotate(VEC::Subtract(V0, vOrigin), vOrientation);
            VECTOR TV1 = VEC3::InverseRotate(VEC::Subtract(V1, vOrigin), vOrientation);
            VECTOR TV2 = VEC3::InverseRotate(VEC::Subtract(V2, vOrigin), vOrientation);

            // Test each vertex of the triangle against the frustum planes.
            VECTOR outside = VEC::FalseInt();
            VECTOR insideAll = VEC::TrueInt();

            for (size_t i = 0; i < 6; ++i)
            {
                VECTOR dist0 = VEC3::Dot(TV0, planes[i]);
                VECTOR dist1 = VEC3::Dot(TV1, planes[i]);
                VECTOR dist2 = VEC3::Dot(TV2, planes[i]);

                VECTOR minDist = VEC::Min(dist0, dist1);
                minDist = VEC::Min(minDist, dist2);
                VECTOR maxDist = VEC::Max(dist0, dist1);
                maxDist = VEC::Max(maxDist, dist2);

                VECTOR planeDist = VEC::SplatW(planes[i]);

                // Outside the plane?
                outside = VEC::OrInt(outside, VEC::Greater(minDist, planeDist));

                // Fully inside the plane?
                insideAll = VEC::AndInt(insideAll, VEC::LessOrEqual(maxDist, planeDist));
            }

            // If the triangle is outside any of the planes it is outside.
            if (VEC4::EqualInt(outside, VEC::TrueInt()))
                return false;

            // If the triangle is inside all planes it is fully inside.
            if (VEC4::EqualInt(insideAll, VEC::TrueInt()))
                return true;

            // Build the corners of the frustum.
            VECTOR vRightTop = VEC::Set(rightSlope, topSlope, 1.0f, 0.0f);
            VECTOR vRightBottom = VEC::Set(rightSlope, bottomSlope, 1.0f, 0.0f);
            VECTOR vLeftTop = VEC::Set(leftSlope, topSlope, 1.0f, 0.0f);
            VECTOR vLeftBottom = VEC::Set(leftSlope, bottomSlope, 1.0f, 0.0f);
            VECTOR vNear = VEC::ReplicatePtr(&near);
            VECTOR vFar = VEC::ReplicatePtr(&far);

            VECTOR corners[CORNER_COUNT];
            corners[0] = VEC::Multiply(vRightTop, vNear);
            corners[1] = VEC::Multiply(vRightBottom, vNear);
            corners[2] = VEC::Multiply(vLeftTop, vNear);
            corners[3] = VEC::Multiply(vLeftBottom, vNear);
            corners[4] = VEC::Multiply(vRightTop, vFar);
            corners[5] = VEC::Multiply(vRightBottom, vFar);
            corners[6] = VEC::Multiply(vLeftTop, vFar);
            corners[7] = VEC::Multiply(vLeftBottom, vFar);

            // Test the plane of the triangle.
            VECTOR normal = VEC3::Cross(VEC::Subtract(V1, V0), VEC::Subtract(V2, V0));
            VECTOR dist = VEC3::Dot(normal, V0);

            VECTOR minDist, maxDist;
            minDist = maxDist = VEC3::Dot(corners[0], normal);
            for (size_t i = 1; i < CORNER_COUNT; ++i)
            {
                VECTOR temp = VEC3::Dot(corners[i], normal);
                minDist = VEC::Min(minDist, temp);
                maxDist = VEC::Max(maxDist, temp);
            }

            outside = VEC::OrInt(VEC::Greater(minDist, dist), VEC::Less(maxDist, dist));
            if (VEC4::EqualInt(outside, VEC::TrueInt()))
                return false;

            // Check the edge/edge axes (3*6).
            VECTOR triangleEdgeAxis[3];
            triangleEdgeAxis[0] = VEC::Subtract(V1, V0);
            triangleEdgeAxis[1] = VEC::Subtract(V2, V1);
            triangleEdgeAxis[2] = VEC::Subtract(V0, V2);

            VECTOR frustumEdgeAxis[6];
            frustumEdgeAxis[0] = vRightTop;
            frustumEdgeAxis[1] = vRightBottom;
            frustumEdgeAxis[2] = vLeftTop;
            frustumEdgeAxis[3] = vLeftBottom;
            frustumEdgeAxis[4] = VEC::Subtract(vRightTop, vLeftTop);
            frustumEdgeAxis[5] = VEC::Subtract(vLeftBottom, vLeftTop);

            for (size_t i = 0; i < 3; ++i)
            {
                for (size_t j = 0; j < 6; j++)
                {
                    // Compute the axis we are going to test.
                    VECTOR axis = VEC3::Cross(triangleEdgeAxis[i], frustumEdgeAxis[j]);

                    // Find the min/max of the projection of the triangle onto the axis.
                    VECTOR minA, maxA;

                    VECTOR dist0 = VEC3::Dot(V0, axis);
                    VECTOR dist1 = VEC3::Dot(V1, axis);
                    VECTOR dist2 = VEC3::Dot(V2, axis);

                    minA = VEC::Min(dist0, dist1);
                    minA = VEC::Min(minA, dist2);
                    maxA = VEC::Max(dist0, dist1);
                    maxA = VEC::Max(maxA, dist2);

                    // Find the min/max of the projection of the frustum onto the axis.
                    VECTOR minB, maxB;

                    minB = maxB = VEC3::Dot(axis, corners[0]);

                    for (size_t k = 1; k < CORNER_COUNT; k++)
                    {
                        VECTOR temp = VEC3::Dot(axis, corners[k]);
                        minB = VEC::Min(minB, temp);
                        maxB = VEC::Max(maxB, temp);
                    }

                    // if (MinA > MaxB || MinB > MaxA) reject;
                    outside = VEC::OrInt(outside, VEC::Greater(minA, maxB));
                    outside = VEC::OrInt(outside, VEC::Greater(minB, maxA));
                }
            }

            if (VEC4::EqualInt(outside, VEC::TrueInt()))
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
            VECTOR vOrigin = VEC::LoadFloat3(&origin);
            VECTOR vOrientation = VEC::LoadFloat4(&orientation);

#if defined(DEBUG) || defined(_DEBUG)
            assert(QuaternionIsUnit(vOrientation));
#endif // DEBUG

            // Set w of the origin to one so we can dot4 with a plane.
            vOrigin = VEC::Insert<0, 0, 0, 0, 1>(vOrigin, VEC::SplatOne());

            // Build the corners of the frustum (in world space).
            VECTOR rightTop = VEC::Set(rightSlope, topSlope, 1.0f, 0.0f);
            VECTOR rightBottom = VEC::Set(rightSlope, bottomSlope, 1.0f, 0.0f);
            VECTOR leftTop = VEC::Set(leftSlope, topSlope, 1.0f, 0.0f);
            VECTOR leftBottom = VEC::Set(leftSlope, bottomSlope, 1.0f, 0.0f);
            VECTOR vNear = VEC::ReplicatePtr(&near);
            VECTOR vFar = VEC::ReplicatePtr(&far);

            rightTop = VEC3::Rotate(rightTop, vOrientation);
            rightBottom = VEC3::Rotate(rightBottom, vOrientation);
            leftTop = VEC3::Rotate(leftTop, vOrientation);
            leftBottom = VEC3::Rotate(leftBottom, vOrientation);

            VECTOR corners0 = VEC::MultiplyAdd(rightTop, vNear, vOrigin);
            VECTOR corners1 = VEC::MultiplyAdd(rightBottom, vNear, vOrigin);
            VECTOR corners2 = VEC::MultiplyAdd(leftTop, vNear, vOrigin);
            VECTOR corners3 = VEC::MultiplyAdd(leftBottom, vNear, vOrigin);
            VECTOR corners4 = VEC::MultiplyAdd(rightTop, vFar, vOrigin);
            VECTOR corners5 = VEC::MultiplyAdd(rightBottom, vFar, vOrigin);
            VECTOR corners6 = VEC::MultiplyAdd(leftTop, vFar, vOrigin);
            VECTOR corners7 = VEC::MultiplyAdd(leftBottom, vFar, vOrigin);

            VECTOR outside, inside;
            FastIntersectFrustumPlane(corners0, corners1, corners2, corners3,
                corners4, corners5, corners6, corners7,
                plane, outside, inside);

            // If the frustum is outside any plane it is outside.
            if (VEC4::EqualInt(outside, VEC::TrueInt()))
                return FRONT;

            // If the frustum is inside all planes it is inside.
            if (VEC4::EqualInt(inside, VEC::TrueInt()))
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
            planes[0] = VEC::Set(0.0f, 0.0f, -1.0f, near);
            planes[1] = VEC::Set(0.0f, 0.0f, 1.0f, -far);
            planes[2] = VEC::Set(1.0f, 0.0f, -rightSlope, 0.0f);
            planes[3] = VEC::Set(-1.0f, 0.0f, leftSlope, 0.0f);
            planes[4] = VEC::Set(0.0f, 1.0f, -topSlope, 0.0f);
            planes[5] = VEC::Set(0.0f, -1.0f, bottomSlope, 0.0f);

            // Load origin and orientation of the frustum.
            VECTOR frOrigin = VEC::LoadFloat3(&origin);
            VECTOR frOrientation = VEC::LoadFloat4(&orientation);

            // This algorithm based on "Fast Ray-Convex Polyhedron Intersectin," in James Arvo, ed., Graphics Gems II pp. 247-250
            float tnear = -FLT_MAX;
            float tfar = FLT_MAX;

            for (size_t i = 0; i < 6; ++i)
            {
                VECTOR plane = PlaneTransform(planes[i], frOrientation, frOrigin);
                plane = PlaneNormalize(plane);

                VECTOR axisDotOrigin = PlaneDotCoord(plane, rayOrigin);
                VECTOR axisDotDirection = VEC3::Dot(plane, direction);

                if (VEC3::LessOrEqual(VEC::Abs(axisDotDirection), g_RayEpsilon))
                {
                    // Ray is parallel to plane - check if ray origin is inside plane's
                    if (VEC3::Greater(axisDotOrigin, g_Zero))
                    {
                        // Ray origin is outside half-space.
                        dist = 0.f;
                        return false;
                    }
                }
                else
                {
                    // Ray not parallel - get distance to plane.
                    float vd = VEC::GetX(axisDotDirection);
                    float vn = VEC::GetX(axisDotOrigin);
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
            VECTOR vOrigin = VEC::LoadFloat3(&origin);
            VECTOR vOrientation = VEC::LoadFloat4(&orientation);

#if defined(DEBUG) || defined(_DEBUG)
            assert(QuaternionIsUnit(vOrientation));
#endif // DEBUG

            // Set w of the origin to one so we can dot4 with a plane.
            vOrigin = VEC::Insert<0, 0, 0, 0, 1>(vOrigin, VEC::SplatOne());

            // Build the corners of the frustum (in world space).
            VECTOR rightTop = VEC::Set(rightSlope, topSlope, 1.0f, 0.0f);
            VECTOR rightBottom = VEC::Set(rightSlope, bottomSlope, 1.0f, 0.0f);
            VECTOR leftTop = VEC::Set(leftSlope, topSlope, 1.0f, 0.0f);
            VECTOR leftBottom = VEC::Set(leftSlope, bottomSlope, 1.0f, 0.0f);
            VECTOR vNear = VEC::ReplicatePtr(&near);
            VECTOR vFar = VEC::ReplicatePtr(&far);

            rightTop = VEC3::Rotate(rightTop, vOrientation);
            rightBottom = VEC3::Rotate(rightBottom, vOrientation);
            leftTop = VEC3::Rotate(leftTop, vOrientation);
            leftBottom = VEC3::Rotate(leftBottom, vOrientation);

            VECTOR corners0 = VEC::MultiplyAdd(rightTop, vNear, vOrigin);
            VECTOR corners1 = VEC::MultiplyAdd(rightBottom, vNear, vOrigin);
            VECTOR corners2 = VEC::MultiplyAdd(leftTop, vNear, vOrigin);
            VECTOR corners3 = VEC::MultiplyAdd(leftBottom, vNear, vOrigin);
            VECTOR corners4 = VEC::MultiplyAdd(rightTop, vFar, vOrigin);
            VECTOR corners5 = VEC::MultiplyAdd(rightBottom, vFar, vOrigin);
            VECTOR corners6 = VEC::MultiplyAdd(leftTop, vFar, vOrigin);
            VECTOR corners7 = VEC::MultiplyAdd(leftBottom, vFar, vOrigin);

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

            anyOutside = VEC::OrInt(anyOutside, outside);
            allInside = VEC::AndInt(allInside, inside);

            FastIntersectFrustumPlane(corners0, corners1, corners2, corners3,
                corners4, corners5, corners6, corners7,
                plane2, outside, inside);

            anyOutside = VEC::OrInt(anyOutside, outside);
            allInside = VEC::AndInt(allInside, inside);

            FastIntersectFrustumPlane(corners0, corners1, corners2, corners3,
                corners4, corners5, corners6, corners7,
                plane3, outside, inside);

            anyOutside = VEC::OrInt(anyOutside, outside);
            allInside = VEC::AndInt(allInside, inside);

            FastIntersectFrustumPlane(corners0, corners1, corners2, corners3,
                corners4, corners5, corners6, corners7,
                plane4, outside, inside);

            anyOutside = VEC::OrInt(anyOutside, outside);
            allInside = VEC::AndInt(allInside, inside);

            FastIntersectFrustumPlane(corners0, corners1, corners2, corners3,
                corners4, corners5, corners6, corners7,
                plane5, outside, inside);

            anyOutside = VEC::OrInt(anyOutside, outside);
            allInside = VEC::AndInt(allInside, inside);

            // If the frustum is outside any plane it is outside.
            if (VEC4::EqualInt(anyOutside, VEC::TrueInt()))
                return DISJOINT;

            // If the frustum is inside all planes it is inside.
            if (VEC4::EqualInt(allInside, VEC::TrueInt()))
                return CONTAINS;

            // The frustum is not inside all planes or outside a plane, it may intersect.
            return INTERSECTS;
        }

        _Use_decl_annotations_
        FORCE_INLINE void BoundingFrustum::GetPlanes(VECTOR* nearPlane, VECTOR* farPlane, VECTOR* rightPlane,
            VECTOR* leftPlane, VECTOR* topPlane, VECTOR* bottomPlane) const noexcept
        {
            // Load origin and orientation of the frustum.
            VECTOR vOrigin = VEC::LoadFloat3(&origin);
            VECTOR vOrientation = VEC::LoadFloat4(&orientation);

            if (nearPlane)
            {
                VECTOR vNearPlane = VEC::Set(0.0f, 0.0f, -1.0f, near);
                vNearPlane = PlaneTransform(vNearPlane, vOrientation, vOrigin);
                *nearPlane = PlaneNormalize(vNearPlane);
            }

            if (farPlane)
            {
                VECTOR vFarPlane = VEC::Set(0.0f, 0.0f, 1.0f, -far);
                vFarPlane = PlaneTransform(vFarPlane, vOrientation, vOrigin);
                *farPlane = PlaneNormalize(vFarPlane);
            }

            if (rightPlane)
            {
                VECTOR vRightPlane = VEC::Set(1.0f, 0.0f, -rightSlope, 0.0f);
                vRightPlane = PlaneTransform(vRightPlane, vOrientation, vOrigin);
                *rightPlane = PlaneNormalize(vRightPlane);
            }

            if (leftPlane)
            {
                VECTOR vLeftPlane = VEC::Set(-1.0f, 0.0f, leftSlope, 0.0f);
                vLeftPlane = PlaneTransform(vLeftPlane, vOrientation, vOrigin);
                *leftPlane = PlaneNormalize(vLeftPlane);
            }

            if (topPlane)
            {
                VECTOR vTopPlane = VEC::Set(0.0f, 1.0f, -topSlope, 0.0f);
                vTopPlane = PlaneTransform(vTopPlane, vOrientation, vOrigin);
                *topPlane = PlaneNormalize(vTopPlane);
            }

            if (bottomPlane)
            {
                VECTOR vBottomPlane = VEC::Set(0.0f, -1.0f, bottomSlope, 0.0f);
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
            MATRIX matInverse = MAT::Inverse(&determinant, projection);

            // Compute the frustum corners in world space.
            VECTOR points[6];

            for (size_t i = 0; i < 6; ++i)
            {
                // Transform point.
                points[i] = VEC4::Transform(NDCPoints[i], matInverse);
            }

            out.origin = Float3(0.0f, 0.0f, 0.0f);
            out.orientation = Float4(0.0f, 0.0f, 0.0f, 1.0f);

            // Compute the slopes.
            points[0] = VEC::Multiply(points[0], VEC::Reciprocal(VEC::SplatZ(points[0])));
            points[1] = VEC::Multiply(points[1], VEC::Reciprocal(VEC::SplatZ(points[1])));
            points[2] = VEC::Multiply(points[2], VEC::Reciprocal(VEC::SplatZ(points[2])));
            points[3] = VEC::Multiply(points[3], VEC::Reciprocal(VEC::SplatZ(points[3])));

            out.rightSlope = VEC::GetX(points[0]);
            out.leftSlope = VEC::GetX(points[1]);
            out.topSlope = VEC::GetY(points[2]);
            out.bottomSlope = VEC::GetY(points[3]);

            // Compute near and far.
            points[4] = VEC::Multiply(points[4], VEC::Reciprocal(VEC::SplatW(points[4])));
            points[5] = VEC::Multiply(points[5], VEC::Reciprocal(VEC::SplatW(points[5])));

            if (rhcoords)
            {
                out.near = VEC::GetZ(points[5]);
                out.far = VEC::GetZ(points[4]);
            }
            else
            {
                out.near = VEC::GetZ(points[4]);
                out.far = VEC::GetZ(points[5]);
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

                VECTOR zero = VEC::Zero();

                VECTOR e1 = VEC::Subtract(V1, V0);
                VECTOR e2 = VEC::Subtract(V2, V0);

                // p = Direction ^ e2;
                VECTOR p = VEC3::Cross(direction, e2);

                // det = e1 * p;
                VECTOR det = VEC3::Dot(e1, p);

                VECTOR u, v, t;

                if (VEC3::GreaterOrEqual(det, g_RayEpsilon))
                {
                    // Determinate is positive (front side of the triangle).
                    VECTOR s = VEC::Subtract(origin, V0);

                    // u = s * p;
                    u = VEC3::Dot(s, p);

                    VECTOR noIntersection = VEC::Less(u, zero);
                    noIntersection = VEC::OrInt(noIntersection, VEC::Greater(u, det));

                    // q = s ^ e1;
                    VECTOR q = VEC3::Cross(s, e1);

                    // v = Direction * q;
                    v = VEC3::Dot(direction, q);

                    noIntersection = VEC::OrInt(noIntersection, VEC::Less(v, zero));
                    noIntersection = VEC::OrInt(noIntersection, VEC::Greater(VEC::Add(u, v), det));

                    // t = e2 * q;
                    t = VEC3::Dot(e2, q);

                    noIntersection = VEC::OrInt(noIntersection, VEC::Less(t, zero));

                    if (VEC4::EqualInt(noIntersection, VEC::TrueInt()))
                    {
                        dist = 0.f;
                        return false;
                    }
                }
                else if (VEC3::LessOrEqual(det, g_RayNegEpsilon))
                {
                    // Determinate is negative (back side of the triangle).
                    VECTOR s = VEC::Subtract(origin, V0);

                    // u = s * p;
                    u = VEC3::Dot(s, p);

                    VECTOR noIntersection = VEC::Greater(u, zero);
                    noIntersection = VEC::OrInt(noIntersection, VEC::Less(u, det));

                    // q = s ^ e1;
                    VECTOR q = VEC3::Cross(s, e1);

                    // v = Direction * q;
                    v = VEC3::Dot(direction, q);

                    noIntersection = VEC::OrInt(noIntersection, VEC::Greater(v, zero));
                    noIntersection = VEC::OrInt(noIntersection, VEC::Less(VEC::Add(u, v), det));

                    // t = e2 * q;
                    t = VEC3::Dot(e2, q);

                    noIntersection = VEC::OrInt(noIntersection, VEC::Greater(t, zero));

                    if (VEC4::EqualInt(noIntersection, VEC::TrueInt()))
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

                t = VEC::Divide(t, det);

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

                VECTOR zero = VEC::Zero();

                // Compute the normal of triangle A.
                VECTOR N1 = VEC3::Cross(VEC::Subtract(A1, A0), VEC::Subtract(A2, A0));

#if defined(DEBUG) || defined(_DEBUG)
                // Assert that the triangle is not degenerate.
                assert(!VEC3::Equal(N1, zero));
#endif // DEBUG

                // Test points of B against the plane of A.
                VECTOR BDist = VEC3::Dot(N1, VEC::Subtract(B0, A0));
                BDist = VEC::Select(BDist, VEC3::Dot(N1, VEC::Subtract(B1, A0)), selectY);
                BDist = VEC::Select(BDist, VEC3::Dot(N1, VEC::Subtract(B2, A0)), selectZ);

                // Ensure robustness with co-planar triangles by zeroing small distances.
                uint32_t BDistIsZeroCR;
                VECTOR BDistIsZero = GreaterR(&BDistIsZeroCR, g_RayEpsilon, VEC::Abs(BDist));
                BDist = VEC::Select(BDist, zero, BDistIsZero);

                uint32_t BDistIsLessCR;
                VECTOR BDistIsLess = VEC::GreaterR(&BDistIsLessCR, zero, BDist);

                uint32_t BDistIsGreaterCR;
                VECTOR BDistIsGreater = VEC::GreaterR(&BDistIsGreaterCR, BDist, zero);

                // If all the points are on the same side we don't intersect.
                if (ComparisonAllTrue(BDistIsLessCR) || ComparisonAllTrue(BDistIsGreaterCR))
                    return false;

                // Compute the normal of triangle B.
                VECTOR N2 = VEC3::Cross(VEC::Subtract(B1, B0), VEC::Subtract(B2, B0));

#if defined(DEBUG) || defined(_DEBUG)
                // Assert that the triangle is not degenerate.
                assert(!VEC3::Equal(N2, zero));
#endif // DEBUG

                // Test points of A against the plane of B.
                VECTOR ADist = VEC3::Dot(N2, VEC::Subtract(A0, B0));
                ADist = VEC::Select(ADist, VEC3::Dot(N2, VEC::Subtract(A1, B0)), selectY);
                ADist = VEC::Select(ADist, VEC3::Dot(N2, VEC::Subtract(A2, B0)), selectZ);

                // Ensure robustness with co-planar triangles by zeroing small distances.
                uint32_t ADistIsZeroCR;
                VECTOR ADistIsZero = VEC::GreaterR(&ADistIsZeroCR, g_RayEpsilon, VEC::Abs(ADist));
                ADist = VEC::Select(ADist, zero, ADistIsZero);

                uint32_t ADistIsLessCR;
                VECTOR ADistIsLess = VEC::GreaterR(&ADistIsLessCR, zero, ADist);

                uint32_t ADistIsGreaterCR;
                VECTOR ADistIsGreater = VEC::GreaterR(&ADistIsGreaterCR, ADist, zero);

                // If all the points are on the same side we don't intersect.
                if (ComparisonAllTrue(ADistIsLessCR) || ComparisonAllTrue(ADistIsGreaterCR))
                    return false;

                // Special case for co-planar triangles.
                if (ComparisonAllTrue(ADistIsZeroCR) || ComparisonAllTrue(BDistIsZeroCR))
                {
                    VECTOR axis, dist, minDist;

                    // Compute an axis perpindicular to the edge (points out).
                    axis = VEC3::Cross(N1, VEC::Subtract(A1, A0));
                    dist = VEC3::Dot(axis, A0);

                    // Test points of B against the axis.
                    minDist = VEC3::Dot(B0, axis);
                    minDist = VEC::Min(minDist, VEC3::Dot(B1, axis));
                    minDist = VEC::Min(minDist, VEC3::Dot(B2, axis));
                    if (VEC4::GreaterOrEqual(minDist, dist))
                        return false;

                    // Edge (A1, A2)
                    axis = VEC3::Cross(N1, VEC::Subtract(A2, A1));
                    dist = VEC3::Dot(axis, A1);

                    minDist = VEC3::Dot(B0, axis);
                    minDist = VEC::Min(minDist, VEC3::Dot(B1, axis));
                    minDist = VEC::Min(minDist, VEC3::Dot(B2, axis));
                    if (VEC4::GreaterOrEqual(minDist, dist))
                        return false;

                    // Edge (A2, A0)
                    axis = VEC3::Cross(N1, VEC::Subtract(A0, A2));
                    dist = VEC3::Dot(axis, A2);

                    minDist = VEC3::Dot(B0, axis);
                    minDist = VEC::Min(minDist, VEC3::Dot(B1, axis));
                    minDist = VEC::Min(minDist, VEC3::Dot(B2, axis));
                    if (VEC4::GreaterOrEqual(minDist, dist))
                        return false;

                    // Edge (B0, B1)
                    axis = VEC3::Cross(N2, VEC::Subtract(B1, B0));
                    dist = VEC3::Dot(axis, B0);

                    minDist = VEC3::Dot(A0, axis);
                    minDist = VEC::Min(minDist, VEC3::Dot(A1, axis));
                    minDist = VEC::Min(minDist, VEC3::Dot(A2, axis));
                    if (VEC4::GreaterOrEqual(minDist, dist))
                        return false;

                    // Edge (B1, B2)
                    axis = VEC3::Cross(N2, VEC::Subtract(B2, B1));
                    dist = VEC3::Dot(axis, B1);

                    minDist = VEC3::Dot(A0, axis);
                    minDist = VEC::Min(minDist, VEC3::Dot(A1, axis));
                    minDist = VEC::Min(minDist, VEC3::Dot(A2, axis));
                    if (VEC4::GreaterOrEqual(minDist, dist))
                        return false;

                    // Edge (B2,B0)
                    axis = VEC3::Cross(N2, VEC::Subtract(B0, B2));
                    dist = VEC3::Dot(axis, B2);

                    minDist = VEC3::Dot(A0, axis);
                    minDist = VEC::Min(minDist, VEC3::Dot(A1, axis));
                    minDist = VEC::Min(minDist, VEC3::Dot(A2, axis));
                    if (VEC4::GreaterOrEqual(minDist, dist))
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
                VECTOR ADistIsLessEqual = VEC::OrInt(ADistIsLess, ADistIsZero);
                VECTOR ADistIsGreaterEqual = VEC::OrInt(ADistIsGreater, ADistIsZero);

                VECTOR AA0, AA1, AA2;
                bool bPositiveA;

                if (Vector3AllTrue(VEC::Select(ADistIsGreaterEqual, ADistIsLess, select0111)) ||
                    Vector3AllTrue(VEC::Select(ADistIsGreater, ADistIsLessEqual, select0111)))
                {
                    // A0 is singular, crossing from positive to negative.
                    AA0 = A0; AA1 = A1; AA2 = A2;
                    bPositiveA = true;
                }
                else if (Vector3AllTrue(VEC::Select(ADistIsLessEqual, ADistIsGreater, select0111)) ||
                    Vector3AllTrue(VEC::Select(ADistIsLess, ADistIsGreaterEqual, select0111)))
                {
                    // A0 is singular, crossing from negative to positive.
                    AA0 = A0; AA1 = A2; AA2 = A1;
                    bPositiveA = false;
                }
                else if (Vector3AllTrue(VEC::Select(ADistIsGreaterEqual, ADistIsLess, select1011)) ||
                    Vector3AllTrue(VEC::Select(ADistIsGreater, ADistIsLessEqual, select1011)))
                {
                    // A1 is singular, crossing from positive to negative.
                    AA0 = A1; AA1 = A2; AA2 = A0;
                    bPositiveA = true;
                }
                else if (Vector3AllTrue(VEC::Select(ADistIsLessEqual, ADistIsGreater, select1011)) ||
                    Vector3AllTrue(VEC::Select(ADistIsLess, ADistIsGreaterEqual, select1011)))
                {
                    // A1 is singular, crossing from negative to positive.
                    AA0 = A1; AA1 = A0; AA2 = A2;
                    bPositiveA = false;
                }
                else if (Vector3AllTrue(VEC::Select(ADistIsGreaterEqual, ADistIsLess, select1101)) ||
                    Vector3AllTrue(VEC::Select(ADistIsGreater, ADistIsLessEqual, select1101)))
                {
                    // A2 is singular, crossing from positive to negative.
                    AA0 = A2; AA1 = A0; AA2 = A1;
                    bPositiveA = true;
                }
                else if (Vector3AllTrue(VEC::Select(ADistIsLessEqual, ADistIsGreater, select1101)) ||
                    Vector3AllTrue(VEC::Select(ADistIsLess, ADistIsGreaterEqual, select1101)))
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

                VECTOR BDistIsLessEqual = VEC::OrInt(BDistIsLess, BDistIsZero);
                VECTOR BDistIsGreaterEqual = VEC::OrInt(BDistIsGreater, BDistIsZero);

                VECTOR BB0, BB1, BB2;
                bool bPositiveB;

                if (Vector3AllTrue(VEC::Select(BDistIsGreaterEqual, BDistIsLess, select0111)) ||
                    Vector3AllTrue(Select(BDistIsGreater, BDistIsLessEqual, select0111)))
                {
                    // B0 is singular, crossing from positive to negative.
                    BB0 = B0; BB1 = B1; BB2 = B2;
                    bPositiveB = true;
                }
                else if (Vector3AllTrue(VEC::Select(BDistIsLessEqual, BDistIsGreater, select0111)) ||
                    Vector3AllTrue(VEC::Select(BDistIsLess, BDistIsGreaterEqual, select0111)))
                {
                    // B0 is singular, crossing from negative to positive.
                    BB0 = B0; BB1 = B2; BB2 = B1;
                    bPositiveB = false;
                }
                else if (Vector3AllTrue(VEC::Select(BDistIsGreaterEqual, BDistIsLess, select1011)) ||
                    Vector3AllTrue(VEC::Select(BDistIsGreater, BDistIsLessEqual, select1011)))
                {
                    // B1 is singular, crossing from positive to negative.
                    BB0 = B1; BB1 = B2; BB2 = B0;
                    bPositiveB = true;
                }
                else if (Vector3AllTrue(VEC::Select(BDistIsLessEqual, BDistIsGreater, select1011)) ||
                    Vector3AllTrue(VEC::Select(BDistIsLess, BDistIsGreaterEqual, select1011)))
                {
                    // B1 is singular, crossing from negative to positive.
                    BB0 = B1; BB1 = B0; BB2 = B2;
                    bPositiveB = false;
                }
                else if (Vector3AllTrue(VEC::Select(BDistIsGreaterEqual, BDistIsLess, select1101)) ||
                    Vector3AllTrue(VEC::Select(BDistIsGreater, BDistIsLessEqual, select1101)))
                {
                    // B2 is singular, crossing from positive to negative.
                    BB0 = B2; BB1 = B0; BB2 = B1;
                    bPositiveB = true;
                }
                else if (Vector3AllTrue(VEC::Select(BDistIsLessEqual, BDistIsGreater, select1101)) ||
                    Vector3AllTrue(VEC::Select(BDistIsLess, BDistIsGreaterEqual, select1101)))
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
                    delta0 = VEC::Subtract(BB0, AA0);
                    delta1 = VEC::Subtract(AA0, BB0);
                }
                else
                {
                    delta0 = VEC::Subtract(AA0, BB0);
                    delta1 = VEC::Subtract(BB0, AA0);
                }

                // Check if the triangles overlap on the line of intersection between the
                // planes of the two triangles by finding the signed line distances.
                VECTOR dist0 = VEC3::Dot(delta0, VEC3::Cross(VEC::Subtract(BB2, BB0), VEC::Subtract(AA2, AA0)));
                if (VEC4::Greater(dist0, zero))
                    return false;

                VECTOR dist1 = VEC3::Dot(delta1, VEC3::Cross(VEC::Subtract(BB1, BB0), VEC::Subtract(AA1, AA0)));
                if (VEC4::Greater(dist1, zero))
                    return false;

                return true;
            }


            //-----------------------------------------------------------------------------
            // Ray-triangle test
            //-----------------------------------------------------------------------------
            _Use_decl_annotations_
            FORCE_INLINE PlaneIntersectionType VEC_CALLCONV Intersects(A_VECTOR V0, A_VECTOR V1, A_VECTOR V2, B_VECTOR plane) noexcept
            {
                VECTOR one = VEC::SplatOne();

#if defined(DEBUG) || defined(_DEBUG)
                assert(PlaneIsUnit(plane));
#endif // DEBUG

                // Set w of the points to one so we can dot4 with a plane.
                VECTOR TV0 = VEC::Insert<0, 0, 0, 0, 1>(V0, one);
                VECTOR TV1 = VEC::Insert<0, 0, 0, 0, 1>(V1, one);
                VECTOR TV2 = VEC::Insert<0, 0, 0, 0, 1>(V2, one);

                VECTOR outside, inside;
                FastIntersectTrianglePlane(TV0, TV1, TV2, plane, outside, inside);

                // If the triangle is outside any plane it is outside.
                if (VEC4::EqualInt(outside, VEC::TrueInt()))
                    return FRONT;

                // If the triangle is inside all planes it is inside.
                if (VEC4::EqualInt(inside, VEC::TrueInt()))
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
                VECTOR one = VEC::SplatOne();

                // Set w of the points to one so we can dot4 with a plane.
                VECTOR TV0 = VEC::Insert<0, 0, 0, 0, 1>(V0, one);
                VECTOR TV1 = VEC::Insert<0, 0, 0, 0, 1>(V1, one);
                VECTOR TV2 = VEC::Insert<0, 0, 0, 0, 1>(V2, one);

                VECTOR outside, inside;

                // Test against each plane.
                FastIntersectTrianglePlane(TV0, TV1, TV2, plane0, outside, inside);

                VECTOR anyOutside = outside;
                VECTOR allInside = inside;

                FastIntersectTrianglePlane(TV0, TV1, TV2, plane1, outside, inside);
                anyOutside = VEC::OrInt(anyOutside, outside);
                allInside = VEC::AndInt(allInside, inside);

                FastIntersectTrianglePlane(TV0, TV1, TV2, plane2, outside, inside);
                anyOutside = VEC::OrInt(anyOutside, outside);
                allInside = VEC::AndInt(allInside, inside);

                FastIntersectTrianglePlane(TV0, TV1, TV2, plane3, outside, inside);
                anyOutside = VEC::OrInt(anyOutside, outside);
                allInside = VEC::AndInt(allInside, inside);

                FastIntersectTrianglePlane(TV0, TV1, TV2, plane4, outside, inside);
                anyOutside = VEC::OrInt(anyOutside, outside);
                allInside = VEC::AndInt(allInside, inside);

                FastIntersectTrianglePlane(TV0, TV1, TV2, plane5, outside, inside);
                anyOutside = VEC::OrInt(anyOutside, outside);
                allInside = VEC::AndInt(allInside, inside);

                // If the triangle is outside any plane it is outside.
                if (VEC4::EqualInt(anyOutside, VEC::TrueInt()))
                    return DISJOINT;

                // If the triangle is inside all planes it is inside.
                if (VEC4::EqualInt(allInside, VEC::TrueInt()))
                    return CONTAINS;

                // The triangle is not inside all planes or outside a plane, it may intersect.
                return INTERSECTS;
            }
        } // namespace TriangleTests
    } // namespace Collision
} // namespace UltReality::Math

#endif // !ULTREALITY_MATH_COLLISION_INL
