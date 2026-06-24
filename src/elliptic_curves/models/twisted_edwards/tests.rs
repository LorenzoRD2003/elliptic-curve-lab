use crate::elliptic_curves::{
    CurveError, TwistedEdwardsCurve,
    traits::{AffineCurveModel, CurveModel, HasJInvariant, LiftXCoordinate, LiftedPoints},
};
use crate::fields::{Fp, traits::Field};

type F2 = Fp<2>;
type F5 = Fp<5>;
type F13 = Fp<13>;

#[test]
fn constructor_rejects_characteristic_two() {
    assert!(matches!(
        TwistedEdwardsCurve::<F2>::new(F2::one(), F2::zero()),
        Err(CurveError::UnsupportedCharacteristic { characteristic: 2 })
    ));
}

#[test]
fn constructor_rejects_zero_a_as_singular() {
    assert!(matches!(
        TwistedEdwardsCurve::<F5>::new(F5::zero(), F5::one()),
        Err(CurveError::SingularCurve)
    ));
}

#[test]
fn constructor_rejects_zero_d_as_singular() {
    assert!(matches!(
        TwistedEdwardsCurve::<F5>::new(F5::one(), F5::zero()),
        Err(CurveError::SingularCurve)
    ));
}

#[test]
fn constructor_rejects_equal_a_and_d_as_singular() {
    assert!(matches!(
        TwistedEdwardsCurve::<F5>::new(F5::one(), F5::one()),
        Err(CurveError::SingularCurve)
    ));
}

#[test]
fn valid_curve_exposes_coefficients_and_equation_string() {
    let curve = TwistedEdwardsCurve::<F5>::new(F5::one(), F5::from_i64(2))
        .expect("sample twisted-Edwards curve should be non-singular");

    assert!(F5::eq(curve.a(), &F5::one()));
    assert!(F5::eq(curve.d(), &F5::from_i64(2)));
    assert_eq!(
        curve.to_equation_string(),
        "(1 (mod 5))x^2 + y^2 = 1 + (2 (mod 5))x^2y^2"
    );
}

#[test]
fn invariants_match_known_small_example() {
    let curve = TwistedEdwardsCurve::<F5>::new(F5::one(), F5::from_i64(2))
        .expect("sample twisted-Edwards curve should be non-singular");

    assert!(F5::eq(&curve.c4(), &F5::from_i64(3)));
    assert!(F5::eq(&curve.c6(), &F5::one()));
    assert!(F5::eq(&curve.discriminant(), &F5::from_i64(2)));
    assert!(F5::eq(&curve.j_invariant(), &F5::one()));
}

#[test]
fn invariants_satisfy_weierstrass_identity() {
    let curve = TwistedEdwardsCurve::<F5>::new(F5::one(), F5::from_i64(2))
        .expect("sample twisted-Edwards curve should be non-singular");

    let left = F5::sub(&F5::cube(&curve.c4()), &F5::square(&curve.c6()));
    let right = F5::mul(&F5::from_i64(1728), &curve.discriminant());

    assert!(F5::eq(&left, &right));
}

#[test]
fn same_curve_has_same_j_invariant() {
    let left = TwistedEdwardsCurve::<F5>::new(F5::one(), F5::from_i64(2))
        .expect("sample twisted-Edwards curve should be non-singular");
    let right = TwistedEdwardsCurve::<F5>::new(F5::one(), F5::from_i64(2))
        .expect("sample twisted-Edwards curve should be non-singular");

    assert!(left.has_same_j_invariant(&right));
}

#[test]
fn curve_model_identity_is_the_finite_point_zero_one() {
    let curve = TwistedEdwardsCurve::<F5>::new(F5::one(), F5::from_i64(2))
        .expect("sample twisted-Edwards curve should be non-singular");
    let identity = curve.identity();

    assert!(curve.is_identity(&identity));
    assert!(curve.contains(&identity));
}

#[test]
fn affine_infinity_is_not_on_the_twisted_edwards_model() {
    let curve = TwistedEdwardsCurve::<F5>::new(F5::one(), F5::from_i64(2))
        .expect("sample twisted-Edwards curve should be non-singular");

    assert!(!curve.contains(&crate::elliptic_curves::AffinePoint::<F5>::infinity()));
}

#[test]
fn point_constructor_accepts_points_on_the_curve_and_rejects_points_off_it() {
    let curve = TwistedEdwardsCurve::<F5>::new(F5::one(), F5::from_i64(2))
        .expect("sample twisted-Edwards curve should be non-singular");

    assert!(curve.point(F5::zero(), F5::one()).is_ok());
    assert_eq!(
        curve.point(F5::one(), F5::one()),
        Err(CurveError::PointNotOnCurve)
    );
}

#[test]
fn has_j_invariant_trait_matches_the_inherent_method() {
    let curve = TwistedEdwardsCurve::<F5>::new(F5::one(), F5::from_i64(2))
        .expect("sample twisted-Edwards curve should be non-singular");

    assert!(F5::eq(
        &<TwistedEdwardsCurve<F5> as HasJInvariant>::j_invariant(&curve),
        &curve.j_invariant()
    ));
}

#[test]
fn lift_x_at_zero_returns_the_identity_and_its_inverse() {
    let curve = TwistedEdwardsCurve::<F5>::new(F5::one(), F5::from_i64(2))
        .expect("sample twisted-Edwards curve should be non-singular");

    match curve.lift_x(F5::zero()).expect("lifting should succeed") {
        LiftedPoints::TwoPoints(left, right) => {
            assert!(curve.is_identity(&left) || curve.is_identity(&right));
            assert!(
                left == crate::elliptic_curves::AffinePoint::new(F5::zero(), F5::one())
                    || right == crate::elliptic_curves::AffinePoint::new(F5::zero(), F5::one())
            );
            assert!(
                left == crate::elliptic_curves::AffinePoint::new(F5::zero(), F5::from_i64(-1))
                    || right
                        == crate::elliptic_curves::AffinePoint::new(F5::zero(), F5::from_i64(-1))
            );
        }
        other => panic!("expected two lifted points above x = 0, got {other:?}"),
    }
}

#[test]
fn lift_x_returns_no_point_when_the_fiber_denominator_vanishes() {
    let curve = TwistedEdwardsCurve::<F13>::new(F13::from_i64(2), F13::one())
        .expect("sample twisted-Edwards curve should be non-singular");

    assert_eq!(
        curve.lift_x(F13::one()).expect("lifting should succeed"),
        LiftedPoints::NoPoint
    );
}

#[test]
fn point_from_x_reports_none_when_no_affine_point_exists() {
    let curve = TwistedEdwardsCurve::<F13>::new(F13::from_i64(2), F13::one())
        .expect("sample twisted-Edwards curve should be non-singular");

    assert_eq!(curve.point_from_x(F13::one()), Ok(None));
}
