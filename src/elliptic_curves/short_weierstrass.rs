use core::fmt;

use crate::fields::{EnumerableFiniteField, Field, SqrtField};

use crate::elliptic_curves::affine::AffinePoint;
use crate::elliptic_curves::error::CurveError;
use crate::elliptic_curves::invariants::HasJInvariant;
use crate::elliptic_curves::isomorphisms::{CurveIsomorphismError, ShortWeierstrassIsomorphism};
use crate::elliptic_curves::traits::{AffineCurveModel, CurveModel, GroupCurveModel, LiftXCoordinate};

/// Short-Weierstrass curve model `y^2 = x^3 + ax + b`.
///
/// This educational implementation currently supports only fields of
/// characteristic different from `2` and `3`, where the classical short form
/// and its discriminant formula behave as expected.
pub struct ShortWeierstrassCurve<F: Field> {
    a: F::Elem,
    b: F::Elem,
}

impl<F> Clone for ShortWeierstrassCurve<F>
where
    F: Field,
    F::Elem: Clone,
{
    fn clone(&self) -> Self {
        Self {
            a: self.a.clone(),
            b: self.b.clone(),
        }
    }
}

impl<F> PartialEq for ShortWeierstrassCurve<F>
where
    F: Field,
    F::Elem: PartialEq,
{
    fn eq(&self, other: &Self) -> bool {
        self.a == other.a && self.b == other.b
    }
}

impl<F> Eq for ShortWeierstrassCurve<F>
where
    F: Field,
    F::Elem: Eq,
{
}

impl<F: Field> ShortWeierstrassCurve<F> {
    /// Builds a validated short-Weierstrass curve descriptor.
    pub fn new(a: F::Elem, b: F::Elem) -> Result<Self, CurveError> {
        let characteristic = F::characteristic();
        if matches!(characteristic, 2 | 3) {
            return Err(CurveError::UnsupportedCharacteristic { characteristic });
        }

        let curve = Self { a, b };
        if F::is_zero(&curve.discriminant()) {
            return Err(CurveError::SingularCurve);
        }

        Ok(curve)
    }

    /// Returns the `a` coefficient in `x^3 + ax + b`.
    pub fn a(&self) -> &F::Elem {
        &self.a
    }

    /// Returns the `b` coefficient in `x^3 + ax + b`.
    pub fn b(&self) -> &F::Elem {
        &self.b
    }

    /// Returns the short-Weierstrass equation as plain text.
    pub fn to_equation_string(&self) -> String
    where
        F::Elem: fmt::Display,
    {
        format!("y^2 = x^3 + ({})x + ({})", self.a, self.b)
    }

    /// Returns the discriminant `Δ = -16(4a^3 + 27b^2)`.
    ///
    /// For a short-Weierstrass equation `y^2 = x^3 + ax + b` in
    /// characteristic different from `2` and `3`, this quantity detects
    /// singularity: `Δ = 0` exactly when the cubic on the right-hand side has a
    /// repeated root, so the model fails to define a smooth elliptic curve.
    pub fn discriminant(&self) -> F::Elem {
        let four = F::from_i64(4);
        let minus_sixteen = F::from_i64(-16);
        let twenty_seven = F::from_i64(27);

        let four_a_cubed = F::mul(&four, &F::cube(&self.a));
        let twenty_seven_b_squared = F::mul(&twenty_seven, &F::square(&self.b));
        let inner = F::add(&four_a_cubed, &twenty_seven_b_squared);
        F::mul(&minus_sixteen, &inner)
    }

    /// Returns the classical Weierstrass invariant `c4 = -48a`.
    ///
    /// Together with [`Self::c6`] and [`Self::discriminant`], this invariant is
    /// part of the standard package attached to a short-Weierstrass model. It
    /// satisfies the classical relation `c4^3 - c6^2 = 1728Δ`.
    pub fn c4(&self) -> F::Elem {
        F::mul(&F::from_i64(-48), &self.a)
    }

    /// Returns the classical Weierstrass invariant `c6 = -864b`.
    ///
    /// This invariant is paired with [`Self::c4`] and
    /// [`Self::discriminant`] in the standard invariant theory of short
    /// Weierstrass equations.
    pub fn c6(&self) -> F::Elem {
        F::mul(&F::from_i64(-864), &self.b)
    }

    /// Returns the `j`-invariant `j = c4^3 / Δ`.
    ///
    /// Over an algebraic closure, the `j`-invariant classifies elliptic curves
    /// up to isomorphism. For a validated [`ShortWeierstrassCurve`], the
    /// discriminant is non-zero by construction, so this quotient is always
    /// defined.
    pub fn j_invariant(&self) -> F::Elem {
        let c4_cubed = F::cube(&self.c4());
        F::div(&c4_cubed, &self.discriminant())
            .expect("validated short Weierstrass curve has non-zero discriminant")
    }

    /// Returns whether this curve and `other` have the same `j`-invariant.
    ///
    /// For short-Weierstrass curves, equality of `j` means the two curves are
    /// isomorphic over an algebraic closure of the base field.
    ///
    /// This is weaker than being isomorphic over the current base field. Over
    /// the base field itself, one must still exhibit a scaling factor
    /// `u in F^*` such that `a' = u^4 a`, `b' = u^6 b`.
    pub fn has_same_j_invariant(&self, other: &Self) -> bool {
        F::eq(&self.j_invariant(), &other.j_invariant())
    }

    /// Returns the short-Weierstrass model obtained from the scaling
    /// parameter `u`.
    ///
    /// Under the convention
    /// `\phi_u : E -> E'`, `(x, y) -> (u^2 x, u^3 y)`,
    /// a curve `E : y^2 = x^3 + ax + b` is sent to `E' : y^2 = x^3 + a'x + b'`
    /// with `a' = u^4 a`, `b' = u^6 b`.
    ///
    /// Since this change of variables is meant to define an isomorphism,
    /// `u` must be invertible in the base field.
    pub fn scaled_by(&self, u: F::Elem) -> Result<Self, CurveIsomorphismError> {
        if F::inv(&u).is_none() {
            return Err(CurveIsomorphismError::NonInvertibleScale);
        }

        let u2 = F::square(&u);
        let u4 = F::square(&u2);
        let u6 = F::mul(&u4, &u2);

        Self::new(F::mul(&u4, &self.a), F::mul(&u6, &self.b)).map_err(Into::into)
    }

    /// Returns whether `other` is exactly the short-Weierstrass model obtained
    /// by scaling this curve with the supplied parameter `u`.
    ///
    /// This is a direct coefficient check for the chosen `u`; it does not
    /// search for any scaling factor on its own.
    pub fn isomorphic_via_scale(&self, other: &Self, u: &F::Elem) -> bool {
        match self.scaled_by(u.clone()) {
            Ok(scaled_curve) => {
                F::eq(scaled_curve.a(), other.a()) && F::eq(scaled_curve.b(), other.b())
            }
            Err(_) => false,
        }
    }

    /// Returns the quadratic twist determined by the non-zero factor `d`.
    ///
    /// For `E : y^2 = x^3 + ax + b`,
    /// the quadratic twist by `d in F^*` is the short-Weierstrass model
    /// `E^(d) : y^2 = x^3 + d^2 a x + d^3 b`.
    ///
    /// This helper constructs the twisted model directly. It does not attempt
    /// to decide whether the twist is trivial or non-trivial over the current
    /// base field. In particular, if `d` is a square in `F`, then `E` and
    /// `E^(d)` are already isomorphic over `F`.
    pub fn quadratic_twist(&self, d: F::Elem) -> Result<Self, CurveIsomorphismError> {
        if F::inv(&d).is_none() {
            return Err(CurveIsomorphismError::NonInvertibleScale);
        }

        let d2 = F::square(&d);
        let d3 = F::mul(&d2, &d);

        Self::new(F::mul(&d2, &self.a), F::mul(&d3, &self.b)).map_err(Into::into)
    }

    /// Divides two field elements under a caller-provided non-zero guarantee.
    fn divide_by_nonzero(
        &self,
        numerator: &F::Elem,
        denominator: &F::Elem,
        context: &'static str,
    ) -> F::Elem {
        let inverse = F::inv(denominator).expect(context);
        F::mul(numerator, &inverse)
    }

    /// Returns the secant slope used to add two distinct affine points.
    ///
    /// For points `P = (x1, y1)` and `Q = (x2, y2)` with `x1 != x2`, the short
    /// Weierstrass addition law uses `m = (y2 - y1) / (x2 - x1)`.
    fn slope_for_addition(
        &self,
        x_left: &F::Elem,
        y_left: &F::Elem,
        x_right: &F::Elem,
        y_right: &F::Elem,
    ) -> F::Elem {
        let numerator = F::sub(y_right, y_left);
        let denominator = F::sub(x_right, x_left);
        self.divide_by_nonzero(
            &numerator,
            &denominator,
            "distinct affine x-coordinates give a non-zero denominator in a field",
        )
    }

    /// Returns the tangent slope used to double an affine point.
    ///
    /// For a finite point `P = (x, y)` with `y != 0`, the short Weierstrass
    /// doubling law uses `m = (3x^2 + a) / (2y)`.
    fn slope_for_doubling(&self, x: &F::Elem, y: &F::Elem) -> F::Elem {
        let numerator = F::add(&F::mul(&F::from_i64(3), &F::square(x)), &self.a);
        let denominator = F::mul(&F::from_i64(2), y);
        self.divide_by_nonzero(
            &numerator,
            &denominator,
            "finite point with non-zero y-coordinate has invertible tangent denominator",
        )
    }

    /// Reconstructs the third intersection point from a known slope.
    ///
    /// Given a slope `m` coming from either a secant or tangent line and a
    /// left input point `(x1, y1)`, this helper applies the standard affine
    /// formulas `x3 = m^2 - x1 - x2` and `y3 = m(x1 - x3) - y1`.
    ///
    /// The returned point is built through [`Self::unchecked_point`] because
    /// this helper is only used from internal paths where the geometric
    /// formulas already guarantee curve membership.
    fn point_from_slope(
        &self,
        slope: &F::Elem,
        x_left: &F::Elem,
        y_left: &F::Elem,
        x_right: &F::Elem,
    ) -> Result<AffinePoint<F>, CurveError> {
        let x_result = F::sub(&F::sub(&F::square(slope), x_left), x_right);
        let y_result = F::sub(&F::mul(slope, &F::sub(x_left, &x_result)), y_left);
        let point = self.unchecked_point(x_result, y_result);
        debug_assert!(self.contains(&point));
        Ok(point)
    }

    /// Adds two curve points under the assumption that both are already valid.
    ///
    /// This helper skips public input validation and is intended for internal
    /// reuse from checked entry points such as [`GroupCurveModel::add`] and the
    /// scalar-multiplication routines.
    fn add_unchecked(
        &self,
        left: &AffinePoint<F>,
        right: &AffinePoint<F>,
    ) -> Result<AffinePoint<F>, CurveError> {
        debug_assert!(self.contains(left));
        debug_assert!(self.contains(right));

        match (left, right) {
            (AffinePoint::Infinity, _) => Ok(right.clone()),
            (_, AffinePoint::Infinity) => Ok(left.clone()),
            (
                AffinePoint::Finite {
                    x: x_left,
                    y: y_left,
                },
                AffinePoint::Finite {
                    x: x_right,
                    y: y_right,
                },
            ) => {
                if self.are_inverse_points(left, right) {
                    return Ok(self.identity());
                }

                if F::eq(x_left, x_right) {
                    return self.double_unchecked(left);
                }

                let slope = self.slope_for_addition(x_left, y_left, x_right, y_right);
                self.point_from_slope(&slope, x_left, y_left, x_right)
            }
        }
    }

    /// Returns whether two finite affine points are additive inverses.
    ///
    /// In short Weierstrass form this happens exactly when the two points have
    /// the same `x`-coordinate and opposite `y`-coordinates.
    fn are_inverse_points(&self, left: &AffinePoint<F>, right: &AffinePoint<F>) -> bool {
        match (left, right) {
            (
                AffinePoint::Finite {
                    x: x_left,
                    y: y_left,
                },
                AffinePoint::Finite {
                    x: x_right,
                    y: y_right,
                },
            ) => F::eq(x_left, x_right) && F::is_zero(&F::add(y_left, y_right)),
            _ => false,
        }
    }

    /// Doubles a curve point under the assumption that it is already valid.
    ///
    /// This helper skips public input validation and exists so checked entry
    /// points can avoid repeating curve-membership work internally.
    fn double_unchecked(&self, point: &AffinePoint<F>) -> Result<AffinePoint<F>, CurveError> {
        debug_assert!(self.contains(point));

        match point {
            AffinePoint::Infinity => Ok(self.identity()),
            AffinePoint::Finite { x, y } => {
                if F::is_zero(y) {
                    return Ok(self.identity());
                }

                let slope = self.slope_for_doubling(x, y);
                self.point_from_slope(&slope, x, y, x)
            }
        }
    }

    /// Multiplies a valid curve point by a non-negative integer without
    /// repeating curve-membership checks along the way.
    fn mul_scalar_unchecked(
        &self,
        point: &AffinePoint<F>,
        scalar: u64,
    ) -> Result<AffinePoint<F>, CurveError> {
        debug_assert!(self.contains(point));

        let mut result = self.identity();
        let mut base = point.clone();
        let mut k = scalar;

        while k > 0 {
            if k & 1 == 1 {
                result = self.add_unchecked(&result, &base)?;
            }

            k >>= 1;

            if k > 0 {
                base = self.double_unchecked(&base)?;
            }
        }

        Ok(result)
    }

    /// Builds a finite affine point without checking the curve equation.
    ///
    /// This is the internal counterpart to [`AffineCurveModel::point`]. It is
    /// intended only for call sites where membership was already validated or
    /// derived from trusted formulas.
    fn unchecked_point(&self, x: F::Elem, y: F::Elem) -> AffinePoint<F> {
        AffinePoint::new(x, y)
    }

    /// Returns the cubic right-hand side `x^3 + ax + b`.
    ///
    /// This internal helper centralizes the polynomial part shared by
    /// membership checks and `x`-coordinate lifting.
    fn rhs_value(&self, x: &F::Elem) -> F::Elem {
        let x_cubed = F::cube(x);
        let ax = F::mul(&self.a, x);
        F::add(&F::add(&x_cubed, &ax), &self.b)
    }
}

impl<F: EnumerableFiniteField> ShortWeierstrassCurve<F> {
    /// Searches exhaustively for a base-field scaling isomorphism from this
    /// curve to `other`.
    ///
    /// It scans every non-zero element `u in F^*`, constructs the scaled model
    /// determined by `u`, and returns the first matching
    /// [`ShortWeierstrassIsomorphism`] when the coefficients agree exactly.
    ///
    /// This is appropriate only for small enumerable fields. It should be read
    /// as a concrete witness search over `F`, not as a scalable algorithm for
    /// large finite fields.
    pub fn find_isomorphism_to(&self, other: &Self) -> Option<ShortWeierstrassIsomorphism<F>> {
        for u in F::elements() {
            if F::is_zero(&u) {
                continue;
            }

            let scaled_curve = match self.scaled_by(u.clone()) {
                Ok(curve) => curve,
                Err(_) => continue,
            };

            if F::eq(scaled_curve.a(), other.a()) && F::eq(scaled_curve.b(), other.b()) {
                return ShortWeierstrassIsomorphism::new(
                    Self {
                        a: self.a.clone(),
                        b: self.b.clone(),
                    },
                    u,
                )
                .ok();
            }
        }

        None
    }

    /// Returns whether this curve is isomorphic to `other` over the current
    /// enumerable base field.
    ///
    /// This method uses the same exhaustive small-field witness search as
    /// [`Self::find_isomorphism_to`], so it is intended only for
    /// work over small enumerable finite fields.
    pub fn is_isomorphic_to(&self, other: &Self) -> bool {
        self.find_isomorphism_to(other).is_some()
    }

    /// Returns every short-Weierstrass scaling automorphism of this curve over
    /// the current enumerable base field. The search is exhaustive over `u in F^*`.
    ///
    /// A scaling parameter defines an automorphism exactly when it fixes the
    /// coefficients: `u^4 a = a`, `u^6 b = b`.
    ///
    /// This is suitable only for small enumerable finite fields.
    pub fn automorphisms(&self) -> Vec<ShortWeierstrassIsomorphism<F>> {
        let mut automorphisms = Vec::new();
        for u in F::elements() {
            if !F::is_zero(&u)
                && let Ok(scaled_curve) = self.scaled_by(u.clone())
                && F::eq(scaled_curve.a(), self.a())
                && F::eq(scaled_curve.b(), self.b())
                && let Ok(isomorphism) = ShortWeierstrassIsomorphism::new(
                    Self {
                        a: self.a.clone(),
                        b: self.b.clone(),
                    },
                    u,
                )
            {
                automorphisms.push(isomorphism);
            }
        }
        automorphisms
    }
}

impl<F> fmt::Display for ShortWeierstrassCurve<F>
where
    F: Field,
    F::Elem: fmt::Display,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.to_equation_string())
    }
}

impl<F: Field> fmt::Debug for ShortWeierstrassCurve<F> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("ShortWeierstrassCurve")
            .field(
                "equation",
                &format_args!("y^2 = x^3 + ({:?})x + ({:?})", self.a, self.b),
            )
            .field("a", &self.a)
            .field("b", &self.b)
            .finish()
    }
}

impl<F: Field> CurveModel for ShortWeierstrassCurve<F> {
    type Elem = F::Elem;
    type BaseField = F;
    type Point = AffinePoint<F>;

    /// Returns the affine identity point at infinity.
    fn identity(&self) -> Self::Point {
        AffinePoint::infinity()
    }

    /// Returns whether the supplied point is the distinguished identity.
    fn is_identity(&self, point: &Self::Point) -> bool {
        point.is_identity()
    }

    /// Checks the short-Weierstrass equation or accepts the point at infinity.
    fn contains(&self, point: &Self::Point) -> bool {
        match point {
            AffinePoint::Infinity => true,
            AffinePoint::Finite { x, y } => {
                let left = F::square(y);
                let right = self.rhs_value(x);
                F::eq(&left, &right)
            }
        }
    }
}

impl<F: Field> HasJInvariant for ShortWeierstrassCurve<F> {
    fn j_invariant(&self) -> Self::Elem {
        ShortWeierstrassCurve::j_invariant(self)
    }
}

impl<F: Field> AffineCurveModel for ShortWeierstrassCurve<F> {
    /// Builds a finite affine point after validating the curve equation.
    fn point(&self, x: Self::Elem, y: Self::Elem) -> Result<Self::Point, CurveError> {
        let point = self.unchecked_point(x, y);
        if self.contains(&point) {
            Ok(point)
        } else {
            Err(CurveError::PointNotOnCurve)
        }
    }
}

impl<F: SqrtField> LiftXCoordinate for ShortWeierstrassCurve<F> {
    /// Returns the cubic right-hand side used by the short-Weierstrass model.
    fn rhs(&self, x: &Self::Elem) -> Self::Elem {
        self.rhs_value(x)
    }
}

impl<F: Field> GroupCurveModel for ShortWeierstrassCurve<F> {
    /// Negates a point by reflecting it across the `x`-axis.
    fn neg(&self, point: &Self::Point) -> Self::Point {
        point.neg()
    }

    /// Adds two affine points using the classical short-Weierstrass formulas.
    ///
    /// This handles the identity, inverse-point, secant, and tangent cases in
    /// the usual geometric way.
    fn add(&self, left: &Self::Point, right: &Self::Point) -> Result<Self::Point, CurveError> {
        if !self.contains(left) || !self.contains(right) {
            return Err(CurveError::PointNotOnCurve);
        }

        self.add_unchecked(left, right)
    }

    /// Subtracts one point from another after validating both public inputs.
    fn sub(&self, left: &Self::Point, right: &Self::Point) -> Result<Self::Point, CurveError> {
        if !self.contains(left) || !self.contains(right) {
            return Err(CurveError::PointNotOnCurve);
        }

        let negated = self.neg(right);
        self.add_unchecked(left, &negated)
    }

    /// Doubles a point using the tangent formula in affine coordinates.
    fn double(&self, point: &Self::Point) -> Result<Self::Point, CurveError> {
        if !self.contains(point) {
            return Err(CurveError::PointNotOnCurve);
        }

        self.double_unchecked(point)
    }

    /// Multiplies a point by a non-negative integer after one upfront validity
    /// check.
    fn mul_scalar(&self, point: &Self::Point, scalar: u64) -> Result<Self::Point, CurveError> {
        if !self.contains(point) {
            return Err(CurveError::PointNotOnCurve);
        }

        self.mul_scalar_unchecked(point, scalar)
    }

    /// Multiplies a point by a signed integer after one upfront validity
    /// check.
    fn mul_scalar_signed(
        &self,
        point: &Self::Point,
        scalar: i64,
    ) -> Result<Self::Point, CurveError> {
        if !self.contains(point) {
            return Err(CurveError::PointNotOnCurve);
        }

        if scalar < 0 {
            let negated = self.neg(point);
            self.mul_scalar_unchecked(&negated, scalar.unsigned_abs())
        } else {
            self.mul_scalar_unchecked(point, scalar as u64)
        }
    }
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeMap;

    use num_bigint::BigInt;
    use num_rational::BigRational;
    use proptest::prelude::*;

    use crate::elliptic_curves::ShortWeierstrassCurve;
    use crate::elliptic_curves::{
        AffineCurveModel, AffinePoint, CurveError, CurveIsomorphism, CurveIsomorphismError,
        CurveModel, EnumerableCurveModel, FiniteAbelianGroupStructure, FiniteGroupCurveModel,
        GroupCurveModel, HasJInvariant, LiftXCoordinate,
    };
    use crate::fields::{EnumerableFiniteField, Field, Fp, Q, SqrtField};
    use crate::proptest_support::non_singular_short_weierstrass_curve;

    type F2 = Fp<2>;
    type F3 = Fp<3>;
    type F5 = Fp<5>;
    type F7 = Fp<7>;
    type F13 = Fp<13>;
    type F17 = Fp<17>;
    type F19 = Fp<19>;
    type F37 = Fp<37>;
    type F43 = Fp<43>;

    fn q(numerator: i64, denominator: i64) -> BigRational {
        BigRational::new(BigInt::from(numerator), BigInt::from(denominator))
    }

    fn f7_curve() -> ShortWeierstrassCurve<F7> {
        ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3)).expect("valid curve")
    }

    fn f5_noncyclic_curve() -> ShortWeierstrassCurve<F5> {
        ShortWeierstrassCurve::<F5>::new(F5::from_i64(-1), F5::zero()).expect("valid curve")
    }

    fn f43_curve() -> ShortWeierstrassCurve<F43> {
        ShortWeierstrassCurve::<F43>::new(F43::from_i64(2), F43::from_i64(3)).expect("valid curve")
    }

    fn f13_j_1728_curve() -> ShortWeierstrassCurve<F13> {
        ShortWeierstrassCurve::<F13>::new(F13::from_i64(2), F13::zero()).expect("valid curve")
    }

    fn f13_j_zero_curve() -> ShortWeierstrassCurve<F13> {
        ShortWeierstrassCurve::<F13>::new(F13::zero(), F13::from_i64(2)).expect("valid curve")
    }

    fn f13_generic_curve() -> ShortWeierstrassCurve<F13> {
        ShortWeierstrassCurve::<F13>::new(F13::from_i64(2), F13::from_i64(3)).expect("valid curve")
    }

    fn f19_curve() -> ShortWeierstrassCurve<F19> {
        ShortWeierstrassCurve::<F19>::new(F19::from_i64(2), F19::from_i64(3)).expect("valid curve")
    }

    fn f37_curve() -> ShortWeierstrassCurve<F37> {
        ShortWeierstrassCurve::<F37>::new(F37::from_i64(2), F37::from_i64(3)).expect("valid curve")
    }

    fn f7_point(x: i64, y: i64) -> AffinePoint<F7> {
        f7_curve()
            .point(F7::from_i64(x), F7::from_i64(y))
            .expect("point should lie on the curve")
    }

    fn assert_contains(curve: &ShortWeierstrassCurve<F7>, point: &AffinePoint<F7>) {
        assert!(curve.contains(point));
    }

    fn assert_group_law(
        curve: &ShortWeierstrassCurve<F7>,
        left: &AffinePoint<F7>,
        right: &AffinePoint<F7>,
        expected: &AffinePoint<F7>,
    ) {
        assert_contains(curve, left);
        assert_contains(curve, right);
        assert_contains(curve, expected);
        assert_eq!(curve.add(left, right), Ok(expected.clone()));
        assert_eq!(curve.sub(expected, right), Ok(left.clone()));
        assert_eq!(curve.sub(expected, left), Ok(right.clone()));
    }

    fn assert_add_commutative(
        curve: &ShortWeierstrassCurve<F7>,
        left: &AffinePoint<F7>,
        right: &AffinePoint<F7>,
    ) {
        let left_right = curve
            .add(left, right)
            .expect("enumerated points should add successfully");
        let right_left = curve
            .add(right, left)
            .expect("enumerated points should add successfully");

        assert_eq!(left_right, right_left);
        assert_contains(curve, &left_right);
    }

    fn assert_add_associative(
        curve: &ShortWeierstrassCurve<F7>,
        left: &AffinePoint<F7>,
        middle: &AffinePoint<F7>,
        right: &AffinePoint<F7>,
    ) {
        let left_grouped = curve
            .add(
                &curve
                    .add(left, middle)
                    .expect("enumerated points should add successfully"),
                right,
            )
            .expect("enumerated points should add successfully");
        let right_grouped = curve
            .add(
                left,
                &curve
                    .add(middle, right)
                    .expect("enumerated points should add successfully"),
            )
            .expect("enumerated points should add successfully");

        assert_eq!(left_grouped, right_grouped);
        assert_contains(curve, &left_grouped);
    }

    fn assert_identity_law(curve: &ShortWeierstrassCurve<F7>, point: &AffinePoint<F7>) {
        assert_eq!(curve.add(&curve.identity(), point), Ok(point.clone()));
        assert_eq!(curve.add(point, &curve.identity()), Ok(point.clone()));
    }

    fn assert_inverse_law(curve: &ShortWeierstrassCurve<F7>, point: &AffinePoint<F7>) {
        let inverse = curve.neg(point);

        assert_eq!(curve.add(point, &inverse), Ok(curve.identity()));
        assert_eq!(curve.add(&inverse, point), Ok(curve.identity()));
    }

    fn assert_scalar_mul_consistent(
        curve: &ShortWeierstrassCurve<F7>,
        point: &AffinePoint<F7>,
        n: u64,
        m: u64,
    ) {
        let left = curve
            .mul_scalar(point, n + m)
            .expect("scalar multiplication should succeed");
        let right = curve
            .add(
                &curve
                    .mul_scalar(point, n)
                    .expect("scalar multiplication should succeed"),
                &curve
                    .mul_scalar(point, m)
                    .expect("scalar multiplication should succeed"),
            )
            .expect("point addition should succeed");

        assert_eq!(left, right);
        assert_contains(curve, &left);
    }

    #[test]
    fn constructor_rejects_characteristics_two_and_three() {
        assert!(matches!(
            ShortWeierstrassCurve::<F2>::new(F2::zero(), F2::one()),
            Err(CurveError::UnsupportedCharacteristic { characteristic: 2 }),
        ));
        assert!(matches!(
            ShortWeierstrassCurve::<F3>::new(F3::zero(), F3::one()),
            Err(CurveError::UnsupportedCharacteristic { characteristic: 3 }),
        ));
    }

    #[test]
    fn constructor_rejects_singular_coefficients() {
        assert!(matches!(
            ShortWeierstrassCurve::<F5>::new(F5::zero(), F5::zero()),
            Err(CurveError::SingularCurve),
        ));
    }

    #[test]
    fn accessors_discriminant_and_rhs_match_the_model() {
        let curve = ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3))
            .expect("valid curve");

        assert!(F7::eq(curve.a(), &F7::from_i64(2)));
        assert!(F7::eq(curve.b(), &F7::from_i64(3)));
        assert!(F7::eq(&curve.discriminant(), &F7::from_i64(3)));
        assert!(F7::eq(&curve.c4(), &F7::from_i64(2)));
        assert!(F7::eq(&curve.c6(), &F7::from_i64(5)));
        assert!(F7::eq(&curve.j_invariant(), &F7::from_i64(5)));
        assert!(F7::eq(
            &LiftXCoordinate::rhs(&curve, &F7::from_i64(2)),
            &F7::from_i64(1)
        ));
    }

    #[test]
    fn weierstrass_invariants_satisfy_the_classical_relation() {
        let curve = ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3))
            .expect("valid curve");
        let left = F7::sub(&F7::cube(&curve.c4()), &F7::square(&curve.c6()));
        let right = F7::mul(&F7::from_i64(1728), &curve.discriminant());

        assert!(F7::eq(&left, &right));
    }

    #[test]
    fn j_invariant_matches_a_classical_exact_example_over_q() {
        let curve = ShortWeierstrassCurve::<Q>::new(q(-1, 1), q(0, 1)).expect("non-singular curve");

        assert!(Q::eq(&curve.discriminant(), &q(64, 1)));
        assert!(Q::eq(&curve.c4(), &q(48, 1)));
        assert!(Q::eq(&curve.c6(), &q(0, 1)));
        assert!(Q::eq(&curve.j_invariant(), &q(1728, 1)));
    }

    #[test]
    fn has_j_invariant_trait_matches_the_inherent_method() {
        let curve = f7_curve();

        assert!(F7::eq(
            &HasJInvariant::j_invariant(&curve),
            &ShortWeierstrassCurve::j_invariant(&curve),
        ));
    }

    #[test]
    fn scaled_by_rejects_noninvertible_scale() {
        let curve = f7_curve();

        assert!(matches!(
            curve.scaled_by(F7::zero()),
            Err(CurveIsomorphismError::NonInvertibleScale)
        ));
    }

    #[test]
    fn scaled_by_applies_the_expected_u4_and_u6_coefficients() {
        let curve = f7_curve();
        let scaled = curve
            .scaled_by(F7::from_i64(3))
            .expect("non-zero scale should define a valid scaled model");

        assert!(F7::eq(scaled.a(), &F7::from_i64(1)));
        assert!(F7::eq(scaled.b(), &F7::from_i64(3)));
    }

    #[test]
    fn isomorphic_via_scale_matches_the_scaled_curve() {
        let curve = f7_curve();
        let scaled = curve
            .scaled_by(F7::from_i64(3))
            .expect("non-zero scale should define a valid scaled model");

        assert!(curve.isomorphic_via_scale(&scaled, &F7::from_i64(3)));
        assert!(!curve.isomorphic_via_scale(&scaled, &F7::from_i64(2)));
    }

    #[test]
    fn isomorphic_via_scale_returns_false_for_noninvertible_scale() {
        let curve = f7_curve();

        assert!(!curve.isomorphic_via_scale(&curve, &F7::zero()));
    }

    fn first_nonsquare<F>() -> F::Elem
    where
        F: EnumerableFiniteField + SqrtField,
    {
        F::elements()
            .into_iter()
            .find(|value| !F::is_zero(value) && !F::has_square_root(value))
            .expect("small odd prime fields should contain non-squares")
    }

    #[test]
    fn quadratic_twist_rejects_zero_factor() {
        let curve = f19_curve();

        assert!(matches!(
            curve.quadratic_twist(F19::zero()),
            Err(CurveIsomorphismError::NonInvertibleScale)
        ));
    }

    #[test]
    fn quadratic_twist_preserves_the_j_invariant_over_f19() {
        let curve = f19_curve();
        let twist = curve
            .quadratic_twist(F19::from_i64(2))
            .expect("non-zero twist factor should produce a valid model");

        assert!(curve.has_same_j_invariant(&twist));
    }

    #[test]
    fn quadratic_twist_by_a_square_is_base_field_isomorphic_over_f19() {
        let curve = f19_curve();
        let square = F19::from_i64(4);
        let twist = curve
            .quadratic_twist(square)
            .expect("square twist factor should produce a valid model");

        assert!(F19::has_square_root(&square));
        assert!(curve.is_isomorphic_to(&twist));
    }

    #[test]
    fn quadratic_twist_by_a_nonsquare_is_not_base_field_isomorphic_in_the_sample_f19_case() {
        let curve = f19_curve();
        let nonsquare = first_nonsquare::<F19>();
        let twist = curve
            .quadratic_twist(nonsquare)
            .expect("non-zero twist factor should produce a valid model");

        assert!(!F19::has_square_root(&nonsquare));
        assert!(curve.has_same_j_invariant(&twist));
        assert!(!curve.is_isomorphic_to(&twist));
    }

    #[test]
    fn quadratic_twist_point_count_relation_holds_over_f19() {
        let curve = f19_curve();
        let nonsquare = first_nonsquare::<F19>();
        let twist = curve
            .quadratic_twist(nonsquare)
            .expect("non-zero twist factor should produce a valid model");

        assert_eq!(curve.order() + twist.order(), 2 * 19 + 2);
    }

    #[test]
    fn quadratic_twist_point_count_relation_holds_over_f37() {
        let curve = f37_curve();
        let nonsquare = first_nonsquare::<F37>();
        let twist = curve
            .quadratic_twist(nonsquare)
            .expect("non-zero twist factor should produce a valid model");

        assert_eq!(curve.order() + twist.order(), 2 * 37 + 2);
    }

    #[test]
    fn find_isomorphism_to_recovers_a_base_field_scaling_witness() {
        let curve = f7_curve();
        let other = curve
            .scaled_by(F7::from_i64(3))
            .expect("non-zero scale should define a valid scaled model");
        let point = f7_point(2, 1);

        let isomorphism = curve
            .find_isomorphism_to(&other)
            .expect("a base-field scaling witness should exist");
        let image = isomorphism
            .evaluate(&point)
            .expect("the witness isomorphism should transport domain points");

        assert!(isomorphism.codomain().contains(&image));
        assert_eq!(
            isomorphism
                .inverse()
                .expect("inverse should exist")
                .evaluate(&image)
                .expect("inverse should recover the original point"),
            point
        );
    }

    #[test]
    fn find_isomorphism_to_returns_none_for_same_j_but_no_base_field_scale() {
        let curve = f7_curve();
        let same_j_not_base_isomorphic =
            ShortWeierstrassCurve::<F7>::new(F7::from_i64(4), F7::from_i64(4))
                .expect("valid curve");

        assert!(curve.has_same_j_invariant(&same_j_not_base_isomorphic));
        assert!(
            curve
                .find_isomorphism_to(&same_j_not_base_isomorphic)
                .is_none()
        );
    }

    #[test]
    fn is_isomorphic_to_returns_true_when_a_base_field_witness_exists() {
        let curve = f7_curve();
        let other = curve
            .scaled_by(F7::from_i64(3))
            .expect("non-zero scale should define a valid scaled model");

        assert!(curve.is_isomorphic_to(&other));
    }

    #[test]
    fn is_isomorphic_to_returns_false_when_only_the_j_invariant_matches() {
        let curve = f7_curve();
        let same_j_not_base_isomorphic =
            ShortWeierstrassCurve::<F7>::new(F7::from_i64(4), F7::from_i64(4))
                .expect("valid curve");

        assert!(curve.has_same_j_invariant(&same_j_not_base_isomorphic));
        assert!(!curve.is_isomorphic_to(&same_j_not_base_isomorphic));
    }

    #[test]
    fn generic_curve_has_only_plus_minus_one_automorphisms_over_f7() {
        let curve = f7_curve();
        let automorphisms = curve.automorphisms();

        assert_eq!(automorphisms.len(), 2);
        assert!(
            automorphisms
                .iter()
                .any(|iso| F7::eq(iso.scaling_factor(), &F7::one()))
        );
        assert!(
            automorphisms
                .iter()
                .any(|iso| F7::eq(iso.scaling_factor(), &F7::from_i64(-1)))
        );
    }

    #[test]
    fn special_j_1728_family_supports_exhaustive_base_field_isomorphism_search() {
        let curve = f13_j_1728_curve();
        let scaled = curve
            .scaled_by(F13::from_i64(2))
            .expect("non-zero scale should define a valid scaled model");

        assert!(F13::is_zero(curve.b()));
        assert!(curve.has_same_j_invariant(&scaled));
        assert!(curve.is_isomorphic_to(&scaled));
        assert!(curve.find_isomorphism_to(&scaled).is_some());
    }

    #[test]
    fn special_j_1728_curve_has_four_automorphisms_over_f13() {
        let curve = f13_j_1728_curve();
        let automorphisms = curve.automorphisms();

        assert_eq!(automorphisms.len(), 4);
        assert!(
            automorphisms
                .iter()
                .any(|iso| F13::eq(iso.scaling_factor(), &F13::one()))
        );
        assert!(
            automorphisms
                .iter()
                .any(|iso| F13::eq(iso.scaling_factor(), &F13::from_i64(-1)))
        );
        assert!(
            automorphisms
                .iter()
                .any(|iso| F13::eq(iso.scaling_factor(), &F13::from_i64(5)))
        );
        assert!(
            automorphisms
                .iter()
                .any(|iso| F13::eq(iso.scaling_factor(), &F13::from_i64(-5)))
        );
    }

    #[test]
    fn special_j_zero_family_supports_exhaustive_base_field_isomorphism_search() {
        let curve = f13_j_zero_curve();
        let scaled = curve
            .scaled_by(F13::from_i64(2))
            .expect("non-zero scale should define a valid scaled model");

        assert!(F13::is_zero(curve.a()));
        assert!(curve.has_same_j_invariant(&scaled));
        assert!(curve.is_isomorphic_to(&scaled));
        assert!(curve.find_isomorphism_to(&scaled).is_some());
    }

    #[test]
    fn special_j_zero_curve_has_six_automorphisms_over_f13() {
        let curve = f13_j_zero_curve();
        let automorphisms = curve.automorphisms();

        assert_eq!(automorphisms.len(), 6);
        assert!(
            automorphisms
                .iter()
                .any(|iso| F13::eq(iso.scaling_factor(), &F13::one()))
        );
        assert!(
            automorphisms
                .iter()
                .any(|iso| F13::eq(iso.scaling_factor(), &F13::from_i64(-1)))
        );
        assert!(
            automorphisms
                .iter()
                .any(|iso| F13::eq(iso.scaling_factor(), &F13::from_i64(4)))
        );
        assert!(
            automorphisms
                .iter()
                .any(|iso| F13::eq(iso.scaling_factor(), &F13::from_i64(-4)))
        );
        assert!(
            automorphisms
                .iter()
                .any(|iso| F13::eq(iso.scaling_factor(), &F13::from_i64(3)))
        );
        assert!(
            automorphisms
                .iter()
                .any(|iso| F13::eq(iso.scaling_factor(), &F13::from_i64(-3)))
        );
    }

    #[test]
    fn generic_family_supports_exhaustive_base_field_isomorphism_search() {
        let curve = f13_generic_curve();
        let scaled = curve
            .scaled_by(F13::from_i64(2))
            .expect("non-zero scale should define a valid scaled model");

        assert!(!F13::is_zero(curve.a()));
        assert!(!F13::is_zero(curve.b()));
        assert!(curve.has_same_j_invariant(&scaled));
        assert!(curve.is_isomorphic_to(&scaled));
        assert!(curve.find_isomorphism_to(&scaled).is_some());
    }

    #[test]
    fn special_j_1728_same_j_does_not_force_base_field_isomorphism() {
        let first = f13_j_1728_curve();
        let second =
            ShortWeierstrassCurve::<F13>::new(F13::from_i64(1), F13::zero()).expect("valid curve");

        assert!(first.has_same_j_invariant(&second));
        assert!(!first.is_isomorphic_to(&second));
        assert!(first.find_isomorphism_to(&second).is_none());
    }

    #[test]
    fn special_j_zero_same_j_does_not_force_base_field_isomorphism() {
        let first = f13_j_zero_curve();
        let second =
            ShortWeierstrassCurve::<F13>::new(F13::zero(), F13::from_i64(1)).expect("valid curve");

        assert!(first.has_same_j_invariant(&second));
        assert!(!first.is_isomorphic_to(&second));
        assert!(first.find_isomorphism_to(&second).is_none());
    }

    #[test]
    fn has_same_j_invariant_detects_scaled_models() {
        let curve = f7_curve();
        let scaled = curve
            .scaled_by(F7::from_i64(3))
            .expect("non-zero scale should define a valid scaled model");

        assert!(curve.has_same_j_invariant(&scaled));
    }

    #[test]
    fn has_same_j_invariant_distinguishes_curves_with_different_j() {
        let first = f7_curve();
        let second = ShortWeierstrassCurve::<F7>::new(F7::from_i64(1), F7::from_i64(1))
            .expect("valid curve");

        assert!(!first.has_same_j_invariant(&second));
    }

    #[test]
    fn contains_accepts_affine_and_infinite_points_on_the_curve() {
        let curve = f7_curve();
        let affine_point = f7_point(2, 1);
        let infinity = AffinePoint::<F7>::infinity();

        assert_contains(&curve, &affine_point);
        assert_contains(&curve, &infinity);
    }

    #[test]
    fn contains_rejects_points_off_the_curve() {
        let curve = ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3))
            .expect("valid curve");
        let point = AffinePoint::<F7>::new(F7::from_i64(2), F7::from_i64(2));

        assert!(!curve.contains(&point));
        assert!(!curve.is_on_curve_nonzero(&point));
    }

    #[test]
    fn point_constructor_accepts_valid_affine_coordinates() {
        let curve = f7_curve();

        let point = curve
            .point(F7::from_i64(2), F7::from_i64(1))
            .expect("point should lie on the curve");

        assert!(matches!(point, AffinePoint::Finite { .. }));
    }

    #[test]
    fn point_constructor_rejects_invalid_affine_coordinates() {
        let curve = ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3))
            .expect("valid curve");

        assert!(matches!(
            curve.point(F7::from_i64(2), F7::from_i64(2)),
            Err(CurveError::PointNotOnCurve)
        ));
    }

    #[test]
    fn characteristic_zero_fields_are_allowed() {
        let curve = ShortWeierstrassCurve::<Q>::new(q(-1, 1), q(0, 1)).expect("non-singular curve");
        let point = curve
            .point(q(0, 1), q(0, 1))
            .expect("point should lie on the curve");

        assert!(curve.contains(&point));
        assert_eq!(Q::characteristic(), 0);
    }

    #[test]
    fn point_from_x_returns_one_point_when_rhs_has_a_square_root() {
        let curve = f7_curve();

        let point = curve
            .point_from_x(F7::from_i64(2))
            .expect("x = 2 should lift to a point");

        assert_contains(&curve, &point);
        assert!(matches!(point, AffinePoint::Finite { .. }));
    }

    #[test]
    fn point_from_x_returns_none_when_rhs_is_not_a_square() {
        let curve = ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3))
            .expect("valid curve");

        assert!(curve.point_from_x(F7::from_i64(0)).is_none());
    }

    #[test]
    fn points_from_x_returns_both_points_when_they_are_distinct() {
        let curve = f7_curve();

        let (left, right) = curve
            .points_from_x(F7::from_i64(2))
            .expect("x = 2 should lift to two points");

        assert_contains(&curve, &left);
        assert_contains(&curve, &right);
        assert_ne!(left, right);
    }

    #[test]
    fn points_from_x_repeats_the_point_when_the_square_root_is_zero() {
        let curve = f7_curve();

        let (left, right) = curve
            .points_from_x(F7::from_i64(6))
            .expect("x = 6 should give y = 0");

        assert_eq!(left, right);
        assert_contains(&curve, &left);
    }

    #[test]
    fn is_on_curve_nonzero_distinguishes_identity_from_finite_points() {
        let curve = ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3))
            .expect("valid curve");
        let finite_point = curve
            .point(F7::from_i64(2), F7::from_i64(1))
            .expect("point should lie on the curve");
        let identity = AffinePoint::<F7>::infinity();

        assert!(curve.contains(&identity));
        assert!(curve.is_identity(&identity));
        assert!(!curve.is_on_curve_nonzero(&identity));

        assert!(curve.contains(&finite_point));
        assert!(!curve.is_identity(&finite_point));
        assert!(curve.is_on_curve_nonzero(&finite_point));
    }

    #[test]
    fn points_from_x_works_over_q_when_an_exact_root_exists() {
        let curve = ShortWeierstrassCurve::<Q>::new(q(-1, 1), q(0, 1)).expect("non-singular curve");

        let (left, right) = curve
            .points_from_x(q(1, 1))
            .expect("x = 1 should give y = 0 in Q");

        assert_eq!(left, right);
        assert!(curve.contains(&left));
    }

    #[test]
    fn finite_point_enumeration_lists_exactly_the_non_identity_points() {
        let curve = f7_curve();
        let finite_points = curve.finite_points();

        assert_eq!(finite_points.len(), 5);
        assert!(
            finite_points
                .iter()
                .all(|point| curve.is_on_curve_nonzero(point))
        );
    }

    #[test]
    fn full_point_enumeration_includes_identity_and_order() {
        let curve = f7_curve();
        let points = curve.points();

        assert_eq!(points.len(), 6);
        assert!(curve.is_identity(points.first().expect("identity should be present")));
        assert_eq!(curve.order(), 6);
    }

    #[test]
    fn random_point_uses_the_supplied_index_sampler() {
        let curve = f7_curve();
        let expected = curve.points()[2].clone();
        let mut sampler = |upper_bound: usize| {
            assert_eq!(upper_bound, 6);
            Some(2)
        };

        let sampled = curve
            .random_point(&mut sampler)
            .expect("sampler should choose an existing point");

        assert_eq!(sampled, expected);
    }

    #[test]
    fn random_point_propagates_sampler_failure() {
        let curve = f7_curve();
        let mut sampler = |_upper_bound: usize| None;

        assert!(curve.random_point(&mut sampler).is_none());
    }

    #[test]
    fn group_negation_matches_affine_involution() {
        let curve = f7_curve();
        let point = f7_point(2, 1);

        assert_eq!(curve.neg(&point), f7_point(2, 6));
        assert_eq!(curve.neg(&curve.identity()), curve.identity());
    }

    #[test]
    fn group_add_handles_identity_and_inverse_cases() {
        let curve = f7_curve();
        let point = f7_point(2, 1);

        assert_identity_law(&curve, &point);
        assert_inverse_law(&curve, &point);
    }

    #[test]
    fn group_add_and_double_match_known_small_field_examples() {
        let curve = f7_curve();
        let p = f7_point(2, 1);
        let q = f7_point(3, 1);
        let two_p = f7_point(3, 6);
        let p_plus_q = f7_point(2, 6);
        let torsion_point = f7_point(6, 0);

        assert_eq!(curve.double(&p), Ok(two_p));
        assert_group_law(&curve, &p, &q, &p_plus_q);
        assert_eq!(curve.sub(&p, &q), Ok(torsion_point));
    }

    #[test]
    fn doubling_a_two_torsion_point_returns_the_identity() {
        let curve = f7_curve();
        let point = f7_point(6, 0);

        assert_eq!(curve.double(&point), Ok(curve.identity()));
    }

    #[test]
    fn scalar_multiplication_matches_repeated_addition_examples() {
        let curve = f7_curve();
        let point = f7_point(2, 1);
        let three_p = f7_point(6, 0);
        let minus_two_p = f7_point(3, 1);

        assert_eq!(curve.mul_scalar(&point, 0), Ok(curve.identity()));
        assert_eq!(curve.mul_scalar(&point, 1), Ok(point.clone()));
        assert_eq!(curve.mul_scalar(&point, 3), Ok(three_p));
        assert_eq!(curve.mul_scalar(&point, 6), Ok(curve.identity()));
        assert_eq!(curve.mul_scalar_signed(&point, -2), Ok(minus_two_p));
        assert_scalar_mul_consistent(&curve, &point, 2, 3);
        assert_scalar_mul_consistent(&curve, &point, 1, 5);
    }

    #[test]
    fn group_operations_reject_points_outside_the_curve() {
        let curve = ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3))
            .expect("valid curve");
        let valid = curve
            .point(F7::from_i64(2), F7::from_i64(1))
            .expect("point should lie on the curve");
        let invalid = AffinePoint::<F7>::new(F7::from_i64(2), F7::from_i64(2));

        assert_eq!(
            curve.add(&valid, &invalid),
            Err(CurveError::PointNotOnCurve)
        );
        assert_eq!(curve.double(&invalid), Err(CurveError::PointNotOnCurve));
        assert_eq!(
            curve.sub(&valid, &invalid),
            Err(CurveError::PointNotOnCurve)
        );
        assert_eq!(
            curve.mul_scalar(&invalid, 3),
            Err(CurveError::PointNotOnCurve)
        );
        assert_eq!(
            curve.mul_scalar_signed(&invalid, -3),
            Err(CurveError::PointNotOnCurve)
        );
    }

    #[test]
    fn enumerated_points_form_an_abelian_group_in_the_small_example() {
        let curve = f7_curve();
        let points = curve.points();

        for left in &points {
            for right in &points {
                assert_add_commutative(&curve, left, right);

                for third in &points {
                    assert_add_associative(&curve, left, right, third);
                }
            }
        }
    }

    #[test]
    fn torsion_helper_detects_known_orders_in_the_small_example() {
        let curve = f7_curve();
        let order_six_point = f7_point(2, 1);
        let order_two_point = f7_point(6, 0);
        let identity = curve.identity();

        assert!(curve.is_torsion_point(&order_six_point, 6));
        assert!(!curve.is_torsion_point(&order_six_point, 3));
        assert!(curve.is_torsion_point(&order_two_point, 2));
        assert!(curve.is_torsion_point(&identity, 5));
    }

    #[test]
    fn torsion_helper_rejects_zero_order_and_invalid_points() {
        let curve = f7_curve();
        let valid = f7_point(2, 1);
        let invalid = AffinePoint::<F7>::new(F7::from_i64(2), F7::from_i64(2));

        assert!(!curve.is_torsion_point(&valid, 0));
        assert!(!curve.is_torsion_point(&invalid, 6));
    }

    #[test]
    fn point_order_matches_known_small_group_examples() {
        let curve = f7_curve();
        let order_six_point = f7_point(2, 1);
        let order_two_point = f7_point(6, 0);

        assert_eq!(curve.point_order(&curve.identity()), Some(1));
        assert_eq!(curve.point_order(&order_two_point), Some(2));
        assert_eq!(curve.point_order(&order_six_point), Some(6));
    }

    #[test]
    fn point_order_returns_none_for_points_outside_the_curve() {
        let curve = ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3))
            .expect("valid curve");
        let invalid = AffinePoint::<F7>::new(F7::from_i64(2), F7::from_i64(2));

        assert_eq!(curve.point_order(&invalid), None);
    }

    #[test]
    fn point_orders_cover_the_full_small_curve_group() {
        let curve = f7_curve();
        let point_orders = curve.point_orders();

        assert_eq!(point_orders.len(), curve.order());
        assert_eq!(point_orders[0], (curve.identity(), 1));
        assert!(point_orders.contains(&(f7_point(6, 0), 2)));
        assert!(point_orders.contains(&(f7_point(2, 1), 6)));
        assert!(point_orders.contains(&(f7_point(2, 6), 6)));
    }

    #[test]
    fn points_of_order_filters_exact_orders() {
        let curve = f7_curve();

        assert_eq!(curve.points_of_order(1), vec![curve.identity()]);
        assert_eq!(curve.points_of_order(2), vec![f7_point(6, 0)]);
        assert_eq!(
            curve.points_of_order(3),
            vec![f7_point(3, 1), f7_point(3, 6)]
        );
        assert_eq!(
            curve.points_of_order(6),
            vec![f7_point(2, 1), f7_point(2, 6)]
        );
        assert!(curve.points_of_order(4).is_empty());
    }

    #[test]
    fn order_distribution_matches_the_small_cyclic_example() {
        let curve = f7_curve();
        let expected = BTreeMap::from([(1, 1), (2, 1), (3, 2), (6, 2)]);

        assert_eq!(curve.order_distribution(), expected);
    }

    #[test]
    fn exponent_generator_and_cyclicity_match_the_small_example() {
        let curve = f7_curve();
        let generator = curve.generator().expect("group should be cyclic");
        let structure = curve.group_structure();

        assert_eq!(curve.exponent(), 6);
        assert!(curve.is_cyclic());
        assert_eq!(
            structure,
            FiniteAbelianGroupStructure {
                order: 6,
                exponent: 6,
                cyclic: true,
                invariant_factors: None,
            }
        );
        assert_eq!(curve.describe_group_structure(), "Z/6Z");
        assert_eq!(curve.point_order(&generator), Some(curve.order()));
    }

    #[test]
    fn noncyclic_f5_example_has_split_two_torsion_structure() {
        let curve = f5_noncyclic_curve();
        let expected = BTreeMap::from([(1, 1), (2, 3), (4, 4)]);
        let structure = curve.group_structure();

        assert_eq!(curve.order(), 8);
        assert_eq!(curve.order_distribution(), expected);
        assert_eq!(curve.exponent(), 4);
        assert_eq!(curve.generator(), None);
        assert!(!curve.is_cyclic());
        assert_eq!(
            structure,
            FiniteAbelianGroupStructure {
                order: 8,
                exponent: 4,
                cyclic: false,
                invariant_factors: Some((2, 4)),
            }
        );
        assert_eq!(curve.describe_group_structure(), "Z/2Z x Z/4Z");
    }

    #[test]
    fn exhaustive_group_axiom_check_passes_for_a_small_f43_curve() {
        let curve = f43_curve();

        assert_eq!(curve.check_group_axioms(), Ok(()));
    }

    fn curve_and_group_data() -> impl Strategy<
        Value = (
            ShortWeierstrassCurve<F17>,
            AffinePoint<F17>,
            AffinePoint<F17>,
            u64,
            u64,
        ),
    > {
        non_singular_short_weierstrass_curve::<17>().prop_flat_map(|curve| {
            let points = curve.points();
            let len = points.len();

            (
                Just(curve.clone()),
                Just(points),
                0usize..len,
                0usize..len,
                0u64..8,
                0u64..8,
            )
                .prop_map(|(curve, points, left_index, right_index, n, m)| {
                    (
                        curve,
                        points[left_index].clone(),
                        points[right_index].clone(),
                        n,
                        m,
                    )
                })
        })
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(28))]

        #[test]
        fn property_short_weierstrass_invariants_satisfy_the_classical_relation(
            curve in non_singular_short_weierstrass_curve::<17>(),
        ) {
            let left = F17::sub(&F17::cube(&curve.c4()), &F17::square(&curve.c6()));
            let right = F17::mul(&F17::from_i64(1728), &curve.discriminant());

            prop_assert!(F17::eq(&left, &right));
        }

        #[test]
        fn property_short_weierstrass_group_law_holds_on_enumerated_points(
            (curve, left, right, n, m) in curve_and_group_data(),
        ) {
            let left_plus_right = curve.add(&left, &right).expect("enumerated points should add");
            let right_plus_left = curve.add(&right, &left).expect("enumerated points should add");
            let inverse = curve.neg(&left);
            let scalar_sum = curve.mul_scalar(&left, n + m).expect("scalar multiplication should succeed");
            let split_scalar = curve
                .add(
                    &curve.mul_scalar(&left, n).expect("scalar multiplication should succeed"),
                    &curve.mul_scalar(&left, m).expect("scalar multiplication should succeed"),
                )
                .expect("point addition should succeed");

            prop_assert_eq!(left_plus_right, right_plus_left);
            prop_assert_eq!(curve.add(&left, &inverse).expect("inverse sum should succeed"), curve.identity());
            prop_assert_eq!(scalar_sum, split_scalar);
        }
    }
}
