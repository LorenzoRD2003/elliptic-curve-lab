use crate::elliptic_curves::{
    CurveError, ProjectivePoint, ShortWeierstrassCurve,
    traits::{CurveModel, HasProjectiveModel, ProjectiveGroupCurveModel},
};
use crate::fields::traits::Field;

#[derive(Clone, Debug)]
enum JacobianPoint<T> {
    Infinity,
    Finite { x: T, y: T, z: T },
}

impl<F: Field> ShortWeierstrassCurve<F> {
    /// Reinterprets one public homogeneous representative `(X : Y : Z)` with
    /// `x = X / Z`, `y = Y / Z` as one internal Jacobian representative
    ///
    /// `x = X_J / Z_J^2`, `y = Y_J / Z_J^3`
    ///
    /// via
    ///
    /// `X_J = X Z`, `Y_J = Y Z^2`, `Z_J = Z`.
    ///
    /// Complexity: `O(1)` field operations, with no inversion in the field.
    fn jacobian_from_projective(&self, point: &ProjectivePoint<F>) -> JacobianPoint<F::Elem> {
        match point {
            ProjectivePoint::Infinity => JacobianPoint::Infinity,
            ProjectivePoint::Finite { x, y, z } => {
                let z2 = F::square(z);
                JacobianPoint::Finite {
                    x: F::mul(x, z),
                    y: F::mul(y, &z2),
                    z: z.clone(),
                }
            }
        }
    }

    /// Returns one public homogeneous representative from one internal
    /// Jacobian point by using the same affine coordinates
    ///
    /// `x = X_J / Z_J^2`, `y = Y_J / Z_J^3`
    ///
    /// and choosing the public chart
    ///
    /// `X = X_J Z_J`, `Y = Y_J`, `Z = Z_J^3`.
    ///
    /// Complexity: `O(1)` field operations, with no inversion in the field.
    fn projective_from_jacobian(&self, point: JacobianPoint<F::Elem>) -> ProjectivePoint<F> {
        match point {
            JacobianPoint::Infinity => ProjectivePoint::Infinity,
            JacobianPoint::Finite { x, y, z } => {
                if F::is_zero(&z) {
                    return ProjectivePoint::Infinity;
                }
                let z_squared = F::square(&z);
                ProjectivePoint::new(F::mul(&x, &z), y, F::mul(&z_squared, &z))
            }
        }
    }

    /// Doubles one internal Jacobian point on
    ///
    /// `Y^2 = X^3 + a X Z^4 + b Z^6`
    ///
    /// through the standard formulas
    ///
    /// `S = 4 X Y^2`, `M = 3 X^2 + a Z^4`,
    ///
    /// `X3 = M^2 - 2S`,
    ///
    /// `Y3 = M(S - X3) - 8 Y^4`,
    ///
    /// `Z3 = 2 Y Z`.
    ///
    /// Complexity: `O(1)` field operations, with no inversion in the field.
    fn jacobian_double_unchecked(&self, point: &JacobianPoint<F::Elem>) -> JacobianPoint<F::Elem> {
        match point {
            JacobianPoint::Infinity => JacobianPoint::Infinity,
            JacobianPoint::Finite { x, y, z } => {
                if F::is_zero(y) || F::is_zero(z) {
                    return JacobianPoint::Infinity;
                }

                let two = F::from_i64(2);
                let three = F::from_i64(3);
                let eight = F::from_i64(8);

                let x_squared = F::square(x);
                let y_squared = F::square(y);
                let y_fourth = F::square(&y_squared);
                let z_squared = F::square(z);
                let z_fourth = F::square(&z_squared);
                let x_plus_y_squared = F::add(x, &y_squared);
                let s = F::mul(
                    &two,
                    &F::sub(
                        &F::sub(&F::square(&x_plus_y_squared), &x_squared),
                        &y_fourth,
                    ),
                );
                let m = F::add(&F::mul(&three, &x_squared), &F::mul(self.a(), &z_fourth));
                let t = F::sub(&F::square(&m), &F::mul(&two, &s));
                let new_y = F::sub(&F::mul(&m, &F::sub(&s, &t)), &F::mul(&eight, &y_fourth));
                let new_z = F::sub(&F::sub(&F::square(&F::add(y, z)), &y_squared), &z_squared);

                JacobianPoint::Finite {
                    x: t,
                    y: new_y,
                    z: new_z,
                }
            }
        }
    }

    /// Adds two internal Jacobian points through the standard formulas
    ///
    /// `U1 = X1 Z2^2`, `U2 = X2 Z1^2`,
    ///
    /// `S1 = Y1 Z2^3`, `S2 = Y2 Z1^3`,
    ///
    /// `H = U2 - U1`, `I = (2H)^2`, `J = H I`,
    ///
    /// `r = 2 (S2 - S1)`, `V = U1 I`,
    ///
    /// `X3 = r^2 - J - 2V`,
    ///
    /// `Y3 = r (V - X3) - 2 S1 J`,
    ///
    /// `Z3 = ((Z1 + Z2)^2 - Z1^2 - Z2^2) H`.
    ///
    /// Complexity: `O(1)` field operations, with no inversion in the field.
    fn jacobian_add_unchecked(
        &self,
        left: &JacobianPoint<F::Elem>,
        right: &JacobianPoint<F::Elem>,
    ) -> JacobianPoint<F::Elem> {
        match (left, right) {
            (JacobianPoint::Infinity, point) | (point, JacobianPoint::Infinity) => point.clone(),
            (
                JacobianPoint::Finite {
                    x: x1,
                    y: y1,
                    z: z1,
                },
                JacobianPoint::Finite {
                    x: x2,
                    y: y2,
                    z: z2,
                },
            ) => {
                if F::is_zero(z1) || F::is_zero(z2) {
                    return JacobianPoint::Infinity;
                }

                let two = F::from_i64(2);

                let z1z1 = F::square(z1);
                let z2z2 = F::square(z2);
                let u1 = F::mul(x1, &z2z2);
                let u2 = F::mul(x2, &z1z1);
                let s1 = F::mul(y1, &F::mul(z2, &z2z2));
                let s2 = F::mul(y2, &F::mul(z1, &z1z1));

                if F::eq(&u1, &u2) {
                    if F::eq(&s1, &s2) {
                        return self.jacobian_double_unchecked(left);
                    }
                    return JacobianPoint::Infinity;
                }

                let h = F::sub(&u2, &u1);
                let i = F::square(&F::mul(&two, &h));
                let j = F::mul(&h, &i);
                let r = F::mul(&two, &F::sub(&s2, &s1));
                let v = F::mul(&u1, &i);
                let x3 = F::sub(&F::sub(&F::square(&r), &j), &F::mul(&two, &v));
                let y3 = F::sub(
                    &F::mul(&r, &F::sub(&v, &x3)),
                    &F::mul(&two, &F::mul(&s1, &j)),
                );
                let z3 = F::mul(
                    &F::sub(&F::sub(&F::square(&F::add(z1, z2)), &z1z1), &z2z2),
                    &h,
                );

                JacobianPoint::Finite {
                    x: x3,
                    y: y3,
                    z: z3,
                }
            }
        }
    }

    /// Adds one internal Jacobian point to one affine point, specialized to
    /// the mixed case where the right input already satisfies `Z_2 = 1`.
    ///
    /// For left Jacobian coordinates
    ///
    /// `x_1 = X_1 / Z_1^2`, `y_1 = Y_1 / Z_1^3`
    ///
    /// and affine right input `(x_2, y_2)`, the current formula uses
    ///
    /// `U_2 = x_2 Z_1^2`, `S_2 = y_2 Z_1^3`,
    ///
    /// `H = U_2 - X_1`, `I = (2H)^2`, `J = H I`,
    ///
    /// `r = 2(S_2 - Y_1)`, `V = X_1 I`,
    ///
    /// `X_3 = r^2 - J - 2V`,
    ///
    /// `Y_3 = r(V - X_3) - 2 Y_1 J`,
    ///
    /// `Z_3 = (Z_1 + H)^2 - Z_1^2 - H^2`.
    ///
    /// Complexity: `O(1)` field operations, with no inversion in the field.
    fn jacobian_mixed_add_unchecked(
        &self,
        left: &JacobianPoint<F::Elem>,
        right_x: &F::Elem,
        right_y: &F::Elem,
    ) -> JacobianPoint<F::Elem> {
        match left {
            JacobianPoint::Infinity => JacobianPoint::Finite {
                x: right_x.clone(),
                y: right_y.clone(),
                z: F::one(),
            },
            JacobianPoint::Finite {
                x: x1,
                y: y1,
                z: z1,
            } => {
                if F::is_zero(z1) {
                    return JacobianPoint::Finite {
                        x: right_x.clone(),
                        y: right_y.clone(),
                        z: F::one(),
                    };
                }

                let two = F::from_i64(2);

                let z1z1 = F::square(z1);
                let u2 = F::mul(right_x, &z1z1);
                let z1_cubed = F::mul(z1, &z1z1);
                let s2 = F::mul(right_y, &z1_cubed);

                if F::eq(x1, &u2) {
                    if F::eq(y1, &s2) {
                        return self.jacobian_double_unchecked(left);
                    }
                    return JacobianPoint::Infinity;
                }

                let h = F::sub(&u2, x1);
                let hh = F::square(&h);
                let i = F::square(&F::mul(&two, &h));
                let j = F::mul(&h, &i);
                let r = F::mul(&two, &F::sub(&s2, y1));
                let v = F::mul(x1, &i);
                let x3 = F::sub(&F::sub(&F::square(&r), &j), &F::mul(&two, &v));
                let y3 = F::sub(
                    &F::mul(&r, &F::sub(&v, &x3)),
                    &F::mul(&two, &F::mul(y1, &j)),
                );
                let z3 = F::sub(&F::sub(&F::square(&F::add(z1, &h)), &z1z1), &hh);

                JacobianPoint::Finite {
                    x: x3,
                    y: y3,
                    z: z3,
                }
            }
        }
    }
}

impl<F: Field> ProjectiveGroupCurveModel for ShortWeierstrassCurve<F> {
    fn neg_projective(&self, point: &Self::ProjectivePoint) -> Self::ProjectivePoint {
        match point {
            ProjectivePoint::Infinity => ProjectivePoint::Infinity,
            ProjectivePoint::Finite { x, y, z } => ProjectivePoint::Finite {
                x: x.clone(),
                y: F::neg(y),
                z: z.clone(),
            },
        }
    }

    fn add_projective(
        &self,
        left: &Self::ProjectivePoint,
        right: &Self::ProjectivePoint,
    ) -> Result<Self::ProjectivePoint, CurveError> {
        if !self.is_projective_point_on_curve(left) || !self.is_projective_point_on_curve(right) {
            return Err(CurveError::PointNotOnCurve);
        }

        Ok(self.projective_from_jacobian(self.jacobian_add_unchecked(
            &self.jacobian_from_projective(left),
            &self.jacobian_from_projective(right),
        )))
    }

    fn double_projective(
        &self,
        point: &Self::ProjectivePoint,
    ) -> Result<Self::ProjectivePoint, CurveError> {
        if !self.is_projective_point_on_curve(point) {
            return Err(CurveError::PointNotOnCurve);
        }

        Ok(self.projective_from_jacobian(
            self.jacobian_double_unchecked(&self.jacobian_from_projective(point)),
        ))
    }

    fn mixed_add_projective(
        &self,
        left: &Self::ProjectivePoint,
        right: &Self::Point,
    ) -> Result<Self::ProjectivePoint, CurveError> {
        if !self.is_projective_point_on_curve(left) || !self.contains(right) {
            return Err(CurveError::PointNotOnCurve);
        }

        match right {
            crate::elliptic_curves::AffinePoint::Infinity => Ok(left.clone()),
            crate::elliptic_curves::AffinePoint::Finite {
                x: right_x,
                y: right_y,
            } => Ok(
                self.projective_from_jacobian(self.jacobian_mixed_add_unchecked(
                    &self.jacobian_from_projective(left),
                    right_x,
                    right_y,
                )),
            ),
        }
    }

    fn mul_scalar_projective(
        &self,
        point: &Self::ProjectivePoint,
        scalar: u64,
    ) -> Result<Self::ProjectivePoint, CurveError> {
        if !self.is_projective_point_on_curve(point) {
            return Err(CurveError::PointNotOnCurve);
        }

        let mut result = self.projective_identity();
        let mut base = point.clone();
        let mut k = scalar;

        while k > 0 {
            if k & 1 == 1 {
                result = self.add_projective(&result, &base)?;
            }
            k >>= 1;
            if k > 0 {
                base = self.double_projective(&base)?;
            }
        }

        Ok(result)
    }
}
