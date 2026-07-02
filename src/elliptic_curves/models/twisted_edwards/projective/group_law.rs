use crate::elliptic_curves::{
    CurveError, TwistedEdwardsCurve,
    traits::{CurveModel, ProjectiveGroupCurveModel},
    twisted_edwards::projective::ExtendedTwistedEdwardsPoint,
};
use crate::fields::traits::*;

type ETEP<F> = ExtendedTwistedEdwardsPoint<F>;

impl<F: Field> TwistedEdwardsCurve<F>
where
    F::Elem: Clone,
{
    /// Validates one native extended-coordinate formula output before exposing
    /// it as a public projective point.
    ///
    /// Points with `Z = 0` are allowed here when they still satisfy the
    /// projective twisted-Edwards equations; only genuinely invalid outputs,
    /// such as the all-zero tuple, are rejected.
    ///
    /// Geometrically, `Z = 0` means the point has left the affine chart
    ///
    /// `x = X / Z`, `y = Y / Z`
    ///
    /// used by the public `AffinePoint<F>` representation, but it can still
    /// define a valid rational point of the projective twisted-Edwards model.
    /// That distinction matters later: the native projective formulas may stay
    /// well-defined even when a subsequent affine recovery must fail honestly
    /// with division by zero.
    fn checked_extended_formula_output(&self, point: ETEP<F>) -> Result<ETEP<F>, CurveError> {
        if self.contains_extended_point(&point) {
            Ok(point)
        } else {
            Err(CurveError::PointNotOnCurve)
        }
    }

    /// Adds two extended points with the native twisted-Edwards formulas in
    /// `(X:Y:Z:T)` coordinates.
    ///
    /// For `E_{a,d}: a x^2 + y^2 = 1 + d x^2 y^2` with
    /// `x = X/Z`, `y = Y/Z`, `T = XY/Z`, the current addition formulas use:
    ///
    /// `A = X1 X2`, `B = Y1 Y2`, `C = d T1 T2`, `D = Z1 Z2`,
    /// `E = (X1 + Y1)(X2 + Y2) - A - B`
    /// `F = D - C`, `G = D + C`, `H = B - a A`
    ///
    /// and return
    ///
    /// `X3 = E F`, `Y3 = G H`, `T3 = E H`, `Z3 = F G`.
    ///
    /// These identities homogenize the generic affine formulas while avoiding
    /// field inversion. The implementation does not claim completeness for all
    /// generic `(a, d)` inputs; it only returns outputs that still satisfy the
    /// projective model equations.
    ///
    /// In particular, this native projective addition may legitimately return
    /// a point with `Z_3 = 0`. That does not mean "addition failed"; it means
    /// the sum exists in the projective model but does not belong to the
    /// affine chart currently used by `AffinePoint<F>`.
    ///
    /// Complexity: `Θ(1)` field operations.
    fn add_projective_unchecked(
        &self,
        left: &ETEP<F>,
        right: &ETEP<F>,
    ) -> Result<ETEP<F>, CurveError> {
        let a = F::mul(left.x(), right.x());
        let b = F::mul(left.y(), right.y());
        let c = F::mul(self.d(), &F::mul(left.t(), right.t()));
        let d = F::mul(left.z(), right.z());
        let e = F::sub(
            &F::sub(
                &F::mul(&F::add(left.x(), left.y()), &F::add(right.x(), right.y())),
                &a,
            ),
            &b,
        );
        let f = F::sub(&d, &c);
        let g = F::add(&d, &c);
        let h = F::sub(&b, &F::mul(self.a(), &a));

        self.checked_extended_formula_output(ETEP::new(
            F::mul(&e, &f),
            F::mul(&g, &h),
            F::mul(&f, &g),
            F::mul(&e, &h),
        ))
    }

    /// Adds one extended point to one affine point specialized to the chart
    /// `Z_2 = 1`, `T_2 = x_2 y_2`.
    ///
    /// The current mixed-add formulas simplify the same extended
    /// twisted-Edwards pattern directly in the affine-right-input case:
    ///
    /// Start from the full addition formulas with an affine right input
    ///
    /// `(x_2, y_2) -> (X_2 : Y_2 : Z_2 : T_2) = (x_2 : y_2 : 1 : x_2 y_2)`.
    ///
    /// Substituting `Z_2 = 1` and `T_2 = t_2 = x_2 y_2` gives
    ///
    /// `A = X_1 x_2`, `B = Y_1 y_2`, `C = d T_1 t_2`, `D = Z_1`,
    ///
    /// `E = (X_1 + Y_1)(x_2 + y_2) - A - B`,
    ///
    /// `F = D - C`, `G = D + C`, `H = B - a A`,
    ///
    /// and return
    ///
    /// `X_3 = E F`, `Y_3 = G H`, `T_3 = E H`, `Z_3 = F G`,
    ///
    /// so the specialization is not a different group law, only the
    /// `Z_2 = 1` simplification of the same extended-coordinate addition.
    ///
    /// As in the full addition case, `Z_3 = 0` is geometrically allowed: it
    /// records that the mixed sum has left the affine chart rather than that
    /// the projective formula has broken down.
    ///
    /// Complexity: `Θ(1)` field operations.
    fn mixed_add_projective_unchecked(
        &self,
        left: &ETEP<F>,
        right_x: &F::Elem,
        right_y: &F::Elem,
    ) -> Result<ETEP<F>, CurveError> {
        let a = F::mul(left.x(), right_x);
        let b = F::mul(left.y(), right_y);
        let c = F::mul(self.d(), &F::mul(left.t(), &F::mul(right_x, right_y)));
        let d = left.z().clone();
        let e = F::sub(
            &F::sub(
                &F::mul(&F::add(left.x(), left.y()), &F::add(right_x, right_y)),
                &a,
            ),
            &b,
        );
        let f = F::sub(&d, &c);
        let g = F::add(&d, &c);
        let h = F::sub(&b, &F::mul(self.a(), &a));

        self.checked_extended_formula_output(ETEP::new(
            F::mul(&e, &f),
            F::mul(&g, &h),
            F::mul(&f, &g),
            F::mul(&e, &h),
        ))
    }

    /// Doubles one extended point with the native twisted-Edwards formulas in
    /// `(X:Y:Z:T)` coordinates.
    ///
    /// The current formulas use
    ///
    /// `A = X^2`, `B = Y^2`, `C = 2 Z^2`, `D = a A`, `E = (X + Y)^2 - A - B`
    /// `G = D + B`, `F = G - C`, `H = D - B`
    ///
    /// and return
    ///
    /// `X([2]P) = E F`, `Y([2]P) = G H`, `T([2]P) = E H`, `Z([2]P) = F G`.
    ///
    /// This avoids inversion and naturally allows `Z = 0` outputs when the
    /// doubled point leaves the affine chart. In that situation the doubling
    /// itself is still a valid projective operation; only a later attempt to
    /// recover affine coordinates must fail honestly.
    ///
    /// Complexity: `Θ(1)` field operations.
    fn double_projective_unchecked(&self, point: &ETEP<F>) -> Result<ETEP<F>, CurveError> {
        let a = F::square(point.x());
        let b = F::square(point.y());
        let c = F::mul(&F::from_i64(2), &F::square(point.z()));
        let d = F::mul(self.a(), &a);
        let e = F::sub(&F::sub(&F::square(&F::add(point.x(), point.y())), &a), &b);
        let g = F::add(&d, &b);
        let f = F::sub(&g, &c);
        let h = F::sub(&d, &b);

        self.checked_extended_formula_output(ETEP::new(
            F::mul(&e, &f),
            F::mul(&g, &h),
            F::mul(&f, &g),
            F::mul(&e, &h),
        ))
    }
}

impl<F: Field> ProjectiveGroupCurveModel for TwistedEdwardsCurve<F>
where
    F::Elem: Clone,
{
    fn neg_projective(&self, point: &Self::ProjectivePoint) -> Self::ProjectivePoint {
        let negated = point.neg();
        debug_assert!(
            !self.contains_extended_point(point) || self.contains_extended_point(&negated),
            "twisted-Edwards projective negation should preserve on-curve inputs"
        );
        negated
    }

    fn add_projective(
        &self,
        left: &Self::ProjectivePoint,
        right: &Self::ProjectivePoint,
    ) -> Result<Self::ProjectivePoint, CurveError> {
        if !self.contains_extended_point(left) || !self.contains_extended_point(right) {
            return Err(CurveError::PointNotOnCurve);
        }
        self.add_projective_unchecked(left, right)
    }

    fn double_projective(
        &self,
        point: &Self::ProjectivePoint,
    ) -> Result<Self::ProjectivePoint, CurveError> {
        if !self.contains_extended_point(point) {
            return Err(CurveError::PointNotOnCurve);
        }
        self.double_projective_unchecked(point)
    }

    fn mixed_add_projective(
        &self,
        left: &Self::ProjectivePoint,
        right: &Self::Point,
    ) -> Result<Self::ProjectivePoint, CurveError> {
        if !self.contains_extended_point(left) || !self.contains(right) {
            return Err(CurveError::PointNotOnCurve);
        }

        let crate::elliptic_curves::AffinePoint::Finite {
            x: right_x,
            y: right_y,
        } = right
        else {
            return Err(CurveError::PointNotOnCurve);
        };

        self.mixed_add_projective_unchecked(left, right_x, right_y)
    }
}
