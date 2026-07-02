use crate::elliptic_curves::{
    AffinePoint, CurveError, GeneralWeierstrassCurve, ProjectivePoint,
    traits::{CurveModel, HasProjectiveModel, ProjectiveGroupCurveModel},
};
use crate::fields::traits::*;

impl<F: Field> GeneralWeierstrassCurve<F>
where
    F::Elem: Clone,
{
    /// Adds two distinct finite projective points with different affine
    /// `x`-coordinates by homogenizing the affine line-slope formulas.
    ///
    /// If we store `x = X / Z`, `y = Y / Z`, the affine slope
    /// `λ = (y2 - y1) / (x2 - x1)` becomes
    ///
    /// `λ = (Y2 Z1 - Y1 Z2) / (X2 Z1 - X1 Z2)`.
    ///
    /// We then homogenize the standard general-Weierstrass formulas
    ///
    /// `x3 = λ^2 + a1 λ - a2 - x1 - x2`
    ///
    /// `y3 = -(λ + a1) x3 - ν - a3`, with `ν = y1 - λ x1`,
    ///
    /// choosing one common denominator so the final result can stay inside the
    /// public `(X : Y : Z)` chart without any inversion.
    ///
    /// Complexity: `O(1)` field operations, with no inversion in the field.
    fn add_distinct_finite_projective_unchecked(
        &self,
        left_x: &F::Elem,
        left_y: &F::Elem,
        left_z: &F::Elem,
        right_x: &F::Elem,
        right_y: &F::Elem,
        right_z: &F::Elem,
    ) -> ProjectivePoint<F> {
        let slope_num = F::sub(&F::mul(right_y, left_z), &F::mul(left_y, right_z));
        let slope_den = F::sub(&F::mul(right_x, left_z), &F::mul(left_x, right_z));
        let z1z2 = F::mul(left_z, right_z);
        let z1_sq = F::square(left_z);
        let z2_sq = F::square(right_z);
        let den_sq = F::square(&slope_den);
        let x3_base = {
            let quadratic_part = F::add(
                &F::add(
                    &F::square(&slope_num),
                    &F::mul(self.a1(), &F::mul(&slope_num, &slope_den)),
                ),
                &F::neg(&F::mul(self.a2(), &den_sq)),
            );
            let x_sum_scaled = F::add(&F::mul(left_x, right_z), &F::mul(right_x, left_z));
            F::sub(
                &F::mul(&quadratic_part, &z1z2),
                &F::mul(&den_sq, &x_sum_scaled),
            )
        };
        let new_z = F::mul(&F::cube(&slope_den), &F::mul(&z1_sq, &z2_sq));
        let new_x = F::mul(&x3_base, &F::mul(&slope_den, &z1z2));
        let new_y = {
            let first_term = F::mul(
                &F::neg(&F::add(&slope_num, &F::mul(self.a1(), &slope_den))),
                &F::mul(&x3_base, &z1z2),
            );
            let intercept_num = F::sub(&F::mul(left_y, &slope_den), &F::mul(&slope_num, left_x));
            let second_term = F::neg(&F::mul(
                &intercept_num,
                &F::mul(&den_sq, &F::mul(left_z, &z2_sq)),
            ));
            let third_term = F::neg(&F::mul(self.a3(), &new_z));
            F::add(&F::add(&first_term, &second_term), &third_term)
        };

        ProjectivePoint::new(new_x, new_y, new_z)
    }

    /// Adds one finite projective point to one affine point, specialized to
    /// the mixed case where the right input already satisfies `Z_2 = 1`.
    ///
    /// If we store the left input as `x_1 = X_1 / Z_1`, `y_1 = Y_1 / Z_1` and
    /// the right input as the affine point `(x_2, y_2)`, then the secant slope
    ///
    /// `λ = (y_2 - y_1) / (x_2 - x_1)`
    ///
    /// becomes
    ///
    /// `λ = (y_2 Z_1 - Y_1) / (x_2 Z_1 - X_1)`.
    ///
    /// We then homogenize the same general-Weierstrass formulas as in the
    /// full projective addition case, but simplify every occurrence of the
    /// right denominator `Z_2` to `1`.
    ///
    /// Complexity: `O(1)` field operations, with no inversion in the field.
    fn mixed_add_finite_projective_unchecked(
        &self,
        left_x: &F::Elem,
        left_y: &F::Elem,
        left_z: &F::Elem,
        right_x: &F::Elem,
        right_y: &F::Elem,
    ) -> ProjectivePoint<F> {
        let slope_num = F::sub(&F::mul(right_y, left_z), left_y);
        let slope_den = F::sub(&F::mul(right_x, left_z), left_x);
        let left_z_sq = F::square(left_z);
        let den_sq = F::square(&slope_den);
        let x3_base = {
            let quadratic_part = F::add(
                &F::add(
                    &F::square(&slope_num),
                    &F::mul(self.a1(), &F::mul(&slope_num, &slope_den)),
                ),
                &F::neg(&F::mul(self.a2(), &den_sq)),
            );
            let x_sum_scaled = F::add(left_x, &F::mul(right_x, left_z));
            F::sub(
                &F::mul(&quadratic_part, left_z),
                &F::mul(&den_sq, &x_sum_scaled),
            )
        };
        let new_z = F::mul(&F::cube(&slope_den), &left_z_sq);
        let new_x = F::mul(&x3_base, &F::mul(&slope_den, left_z));
        let new_y = {
            let first_term = F::mul(
                &F::neg(&F::add(&slope_num, &F::mul(self.a1(), &slope_den))),
                &F::mul(&x3_base, left_z),
            );
            let intercept_num = F::sub(&F::mul(left_y, &slope_den), &F::mul(&slope_num, left_x));
            let second_term = F::neg(&F::mul(&intercept_num, &F::mul(&den_sq, left_z)));
            let third_term = F::neg(&F::mul(self.a3(), &new_z));
            F::add(&F::add(&first_term, &second_term), &third_term)
        };

        ProjectivePoint::new(new_x, new_y, new_z)
    }

    /// Doubles one finite projective point by homogenizing the affine tangent
    /// formulas for the general Weierstrass model.
    ///
    /// In affine coordinates the tangent slope is
    ///
    /// `λ = (3x^2 + 2 a2 x + a4 - a1 y) / (2y + a1 x + a3)`.
    ///
    /// Under the public chart `x = X/Z`, `y = Y/Z`, the same slope becomes
    /// `λ = N/(DZ)` with
    ///
    /// `N = 3 X^2 + 2 a2 X Z + a4 Z^2 - a1 Y Z`
    ///
    /// `D = 2 Y + a1 X + a3 Z`.
    ///
    /// We then homogenize the same affine formulas for `x3` and `y3`, again
    /// keeping all work inside the field operations and avoiding inversion.
    ///
    /// Complexity: `O(1)` field operations, with no inversion in the field.
    fn double_finite_projective_unchecked(
        &self,
        x: &F::Elem,
        y: &F::Elem,
        z: &F::Elem,
    ) -> ProjectivePoint<F> {
        let two = F::from_i64(2);
        let three = F::from_i64(3);

        let tangent_den = F::add(
            &F::mul(&two, y),
            &F::add(&F::mul(self.a1(), x), &F::mul(self.a3(), z)),
        );
        if F::is_zero(&tangent_den) {
            return ProjectivePoint::Infinity;
        }

        let z_sq = F::square(z);
        let tangent_num = F::sub(
            &F::add(
                &F::add(
                    &F::mul(&three, &F::square(x)),
                    &F::mul(&F::mul(&two, self.a2()), &F::mul(x, z)),
                ),
                &F::mul(self.a4(), &z_sq),
            ),
            &F::mul(self.a1(), &F::mul(y, z)),
        );
        let den_sq = F::square(&tangent_den);
        let x3_base = F::sub(
            &F::add(
                &F::add(
                    &F::square(&tangent_num),
                    &F::mul(self.a1(), &F::mul(&tangent_num, &F::mul(&tangent_den, z))),
                ),
                &F::neg(&F::mul(self.a2(), &F::mul(&den_sq, &z_sq))),
            ),
            &F::mul(&two, &F::mul(x, &F::mul(&den_sq, z))),
        );
        let new_z = F::mul(&F::cube(&tangent_den), &F::cube(z));
        let new_x = F::mul(&x3_base, &F::mul(&tangent_den, z));
        let new_y = {
            let first_term = F::mul(
                &F::neg(&F::add(
                    &tangent_num,
                    &F::mul(self.a1(), &F::mul(&tangent_den, z)),
                )),
                &x3_base,
            );
            let intercept_num = F::sub(
                &F::mul(y, &F::mul(&tangent_den, z)),
                &F::mul(&tangent_num, x),
            );
            let second_term = F::neg(&F::mul(&intercept_num, &F::mul(&den_sq, z)));
            let third_term = F::neg(&F::mul(self.a3(), &new_z));
            F::add(&F::add(&first_term, &second_term), &third_term)
        };

        ProjectivePoint::new(new_x, new_y, new_z)
    }
}

impl<F: Field> ProjectiveGroupCurveModel for GeneralWeierstrassCurve<F>
where
    F::Elem: Clone,
{
    fn neg_projective(&self, point: &Self::ProjectivePoint) -> Self::ProjectivePoint {
        match point {
            ProjectivePoint::Infinity => ProjectivePoint::Infinity,
            ProjectivePoint::Finite { x, y, z } => ProjectivePoint::Finite {
                x: x.clone(),
                y: F::sub(
                    &F::neg(y),
                    &F::add(&F::mul(self.a1(), x), &F::mul(self.a3(), z)),
                ),
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

        if left == right {
            return self.double_projective(left);
        }

        match (left, right) {
            (ProjectivePoint::Infinity, point) | (point, ProjectivePoint::Infinity) => {
                Ok(point.clone())
            }
            (
                ProjectivePoint::Finite {
                    x: left_x,
                    y: left_y,
                    z: left_z,
                },
                ProjectivePoint::Finite {
                    x: right_x,
                    y: right_y,
                    z: right_z,
                },
            ) => {
                if F::eq(&F::mul(left_x, right_z), &F::mul(right_x, left_z)) {
                    return Ok(ProjectivePoint::Infinity);
                }

                Ok(self.add_distinct_finite_projective_unchecked(
                    left_x, left_y, left_z, right_x, right_y, right_z,
                ))
            }
        }
    }

    fn double_projective(
        &self,
        point: &Self::ProjectivePoint,
    ) -> Result<Self::ProjectivePoint, CurveError> {
        if !self.is_projective_point_on_curve(point) {
            return Err(CurveError::PointNotOnCurve);
        }

        match point {
            ProjectivePoint::Infinity => Ok(ProjectivePoint::Infinity),
            ProjectivePoint::Finite { x, y, z } => {
                Ok(self.double_finite_projective_unchecked(x, y, z))
            }
        }
    }

    fn mixed_add_projective(
        &self,
        left: &Self::ProjectivePoint,
        right: &Self::Point,
    ) -> Result<Self::ProjectivePoint, CurveError> {
        if !self.is_projective_point_on_curve(left) || !self.contains(right) {
            return Err(CurveError::PointNotOnCurve);
        }

        let right_projective = ProjectivePoint::from_affine(right);
        if left == &right_projective {
            return self.double_projective(left);
        }

        match (left, right) {
            (ProjectivePoint::Infinity, point) => self.to_projective(point),
            (point, AffinePoint::Infinity) => Ok(point.clone()),
            (
                ProjectivePoint::Finite {
                    x: left_x,
                    y: left_y,
                    z: left_z,
                },
                AffinePoint::Finite {
                    x: right_x,
                    y: right_y,
                },
            ) => {
                if F::eq(left_x, &F::mul(right_x, left_z)) {
                    return Ok(ProjectivePoint::Infinity);
                }

                Ok(self.mixed_add_finite_projective_unchecked(
                    left_x, left_y, left_z, right_x, right_y,
                ))
            }
        }
    }
}
