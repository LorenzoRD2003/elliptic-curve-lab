use crate::elliptic_curves::{
    AffinePoint, CurveError,
    traits::{AffineCurveModel, CurveModel, LiftXCoordinate},
};
use crate::fields::{Field, SqrtField};

/// Short-Weierstrass curve model `y^2 = x^3 + ax + b`.
///
/// This educational implementation currently supports only fields of
/// characteristic different from `2` and `3`, where the classical short form
/// and its discriminant formula behave as expected.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ShortWeierstrassCurve<F: Field> {
    a: F::Elem,
    b: F::Elem,
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

    fn rhs_value(&self, x: &F::Elem) -> F::Elem {
        let x_cubed = F::cube(x);
        let ax = F::mul(&self.a, x);
        F::add(&F::add(&x_cubed, &ax), &self.b)
    }
}

impl<F: Field> CurveModel for ShortWeierstrassCurve<F> {
    type Elem = F::Elem;
    type BaseField = F;
    type Point = AffinePoint<F>;

    fn identity(&self) -> Self::Point {
        AffinePoint::infinity()
    }

    fn is_identity(&self, point: &Self::Point) -> bool {
        point.is_identity()
    }

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

impl<F: Field> AffineCurveModel for ShortWeierstrassCurve<F> {
    fn point(&self, x: Self::Elem, y: Self::Elem) -> Result<Self::Point, CurveError> {
        let point = AffinePoint::new(x, y);
        if self.contains(&point) {
            Ok(point)
        } else {
            Err(CurveError::PointNotOnCurve)
        }
    }
}

impl<F: SqrtField> LiftXCoordinate for ShortWeierstrassCurve<F> {
    fn rhs(&self, x: &Self::Elem) -> Self::Elem {
        self.rhs_value(x)
    }
}

#[cfg(test)]
mod tests {
    use num_bigint::BigInt;
    use num_rational::BigRational;

    use super::ShortWeierstrassCurve;
    use crate::elliptic_curves::{
        AffineCurveModel, AffinePoint, CurveError, CurveModel, EnumerableCurveModel,
        LiftXCoordinate,
    };
    use crate::fields::{Field, Fp, Q};

    type F2 = Fp<2>;
    type F3 = Fp<3>;
    type F5 = Fp<5>;
    type F7 = Fp<7>;

    fn q(numerator: i64, denominator: i64) -> BigRational {
        BigRational::new(BigInt::from(numerator), BigInt::from(denominator))
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
    fn contains_accepts_affine_and_infinite_points_on_the_curve() {
        let curve = ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3))
            .expect("valid curve");
        let affine_point = curve
            .point(F7::from_i64(2), F7::from_i64(1))
            .expect("point should lie on the curve");
        let infinity = AffinePoint::<F7>::infinity();

        assert!(curve.contains(&affine_point));
        assert!(curve.contains(&infinity));
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
        let curve = ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3))
            .expect("valid curve");

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
        let curve = ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3))
            .expect("valid curve");

        let point = curve
            .point_from_x(F7::from_i64(2))
            .expect("x = 2 should lift to a point");

        assert!(curve.contains(&point));
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
        let curve = ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3))
            .expect("valid curve");

        let (left, right) = curve
            .points_from_x(F7::from_i64(2))
            .expect("x = 2 should lift to two points");

        assert!(curve.contains(&left));
        assert!(curve.contains(&right));
        assert_ne!(left, right);
    }

    #[test]
    fn points_from_x_repeats_the_point_when_the_square_root_is_zero() {
        let curve = ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3))
            .expect("valid curve");

        let (left, right) = curve
            .points_from_x(F7::from_i64(6))
            .expect("x = 6 should give y = 0");

        assert_eq!(left, right);
        assert!(curve.contains(&left));
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
        let curve = ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3))
            .expect("valid curve");
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
        let curve = ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3))
            .expect("valid curve");
        let points = curve.points();

        assert_eq!(points.len(), 6);
        assert!(curve.is_identity(points.first().expect("identity should be present")));
        assert_eq!(curve.order(), 6);
    }

    #[test]
    fn random_point_uses_the_supplied_index_sampler() {
        let curve = ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3))
            .expect("valid curve");
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
        let curve = ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3))
            .expect("valid curve");
        let mut sampler = |_upper_bound: usize| None;

        assert!(curve.random_point(&mut sampler).is_none());
    }
}
