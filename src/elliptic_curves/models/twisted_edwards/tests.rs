use crate::elliptic_curves::{
    AffinePoint, CurveError, GeneralWeierstrassCurve, MontgomeryCurve, ShortWeierstrassCurve,
    TwistedEdwardsCurve,
    traits::{AffineCurveModel, CurveModel, HasJInvariant, LiftXCoordinate, LiftedPoints},
};
use crate::elliptic_curves::models::twisted_edwards::TwistedEdwardsBirationalMapError;
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

    assert!(!curve.contains(&AffinePoint::<F5>::infinity()));
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
                left == AffinePoint::new(F5::zero(), F5::one())
                    || right == AffinePoint::new(F5::zero(), F5::one())
            );
            assert!(
                left == AffinePoint::new(F5::zero(), F5::from_i64(-1))
                    || right == AffinePoint::new(F5::zero(), F5::from_i64(-1))
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

#[test]
fn twisted_edwards_to_montgomery_matches_known_small_example() {
    let curve = TwistedEdwardsCurve::<F5>::new(F5::one(), F5::from_i64(2))
        .expect("sample twisted-Edwards curve should be non-singular");
    let montgomery = curve.as_montgomery();

    assert!(F5::eq(montgomery.a(), &F5::from_i64(4)));
    assert!(F5::eq(montgomery.b(), &F5::one()));
}

#[test]
fn twisted_edwards_to_montgomery_roundtrip_preserves_coefficients() {
    let twisted = TwistedEdwardsCurve::<F13>::new(F13::from_i64(3), F13::from_i64(5))
        .expect("sample twisted-Edwards curve should be non-singular");
    let roundtrip = twisted.as_montgomery().as_twisted_edwards();

    assert!(F13::eq(twisted.a(), roundtrip.a()));
    assert!(F13::eq(twisted.d(), roundtrip.d()));
}

#[test]
fn montgomery_to_twisted_edwards_roundtrip_preserves_coefficients() {
    let montgomery = MontgomeryCurve::<F13>::new(F13::from_i64(3), F13::from_i64(2))
        .expect("sample Montgomery curve should be non-singular");
    let roundtrip = montgomery.as_twisted_edwards().as_montgomery();

    assert!(F13::eq(montgomery.a(), roundtrip.a()));
    assert!(F13::eq(montgomery.b(), roundtrip.b()));
}

#[test]
fn direct_and_composed_montgomery_j_invariants_agree() {
    let twisted = TwistedEdwardsCurve::<F13>::new(F13::from_i64(3), F13::from_i64(5))
        .expect("sample twisted-Edwards curve should be non-singular");
    let montgomery = twisted.as_montgomery();

    assert!(F13::eq(&twisted.j_invariant(), &montgomery.j_invariant()));
}

#[test]
fn composed_general_weierstrass_conversion_preserves_j_invariant() {
    let twisted = TwistedEdwardsCurve::<F13>::new(F13::from_i64(3), F13::from_i64(5))
        .expect("sample twisted-Edwards curve should be non-singular");
    let general: GeneralWeierstrassCurve<F13> = twisted.as_montgomery().as_general_weierstrass();

    assert!(F13::eq(&twisted.j_invariant(), &general.j_invariant()));
}

#[test]
fn composed_short_weierstrass_conversion_preserves_j_invariant() {
    let twisted = TwistedEdwardsCurve::<F5>::new(F5::one(), F5::from_i64(2))
        .expect("sample twisted-Edwards curve should be non-singular");
    let short: ShortWeierstrassCurve<F5> = twisted
        .as_montgomery()
        .try_as_short_weierstrass()
        .expect("characteristic 5 should support the Montgomery short companion");

    assert!(F5::eq(&twisted.j_invariant(), &short.j_invariant()));
}

#[test]
fn birational_open_transport_to_montgomery_maps_a_regular_affine_point() {
    let twisted = TwistedEdwardsCurve::<F5>::new(F5::one(), F5::from_i64(2))
        .expect("sample twisted-Edwards curve should be non-singular");
    let point = twisted
        .point(F5::one(), F5::zero())
        .expect("sample point should lie on the twisted-Edwards curve");

    let image = twisted
        .try_point_to_montgomery_open(&point)
        .expect("point lies in the birational open");

    assert_eq!(image, AffinePoint::new(F5::one(), F5::one()));
}

#[test]
fn birational_open_transport_from_montgomery_maps_a_regular_affine_point() {
    let montgomery = MontgomeryCurve::<F5>::new(F5::from_i64(4), F5::one())
        .expect("sample Montgomery curve should be non-singular");
    let point = montgomery
        .point(F5::one(), F5::one())
        .expect("sample point should lie on the Montgomery curve");

    let image = montgomery
        .try_point_to_twisted_edwards_open(&point)
        .expect("point lies in the birational open");

    assert_eq!(image, AffinePoint::new(F5::one(), F5::zero()));
}

#[test]
fn birational_open_transport_roundtrips_on_the_certified_domain() {
    let twisted = TwistedEdwardsCurve::<F5>::new(F5::one(), F5::from_i64(2))
        .expect("sample twisted-Edwards curve should be non-singular");
    let montgomery = twisted.as_montgomery();
    let point = twisted
        .point(F5::one(), F5::zero())
        .expect("sample point should lie on the twisted-Edwards curve");

    let montgomery_point = twisted
        .try_point_to_montgomery_open(&point)
        .expect("point lies in the birational open");
    let roundtrip = montgomery
        .try_point_to_twisted_edwards_open(&montgomery_point)
        .expect("mapped point stays in the birational open");

    assert_eq!(roundtrip, point);
}

#[test]
fn birational_open_transport_rejects_twisted_edwards_identity() {
    let twisted = TwistedEdwardsCurve::<F5>::new(F5::one(), F5::from_i64(2))
        .expect("sample twisted-Edwards curve should be non-singular");

    assert_eq!(
        twisted.try_point_to_montgomery_open(&twisted.identity()),
        Err(TwistedEdwardsBirationalMapError::ExceptionalTwistedEdwardsPoint)
    );
}

#[test]
fn birational_open_transport_rejects_the_other_x_zero_twisted_edwards_point() {
    let twisted = TwistedEdwardsCurve::<F5>::new(F5::one(), F5::from_i64(2))
        .expect("sample twisted-Edwards curve should be non-singular");
    let point = twisted
        .point(F5::zero(), F5::from_i64(-1))
        .expect("the second x = 0 point should lie on the twisted-Edwards curve");

    assert_eq!(
        twisted.try_point_to_montgomery_open(&point),
        Err(TwistedEdwardsBirationalMapError::ExceptionalTwistedEdwardsPoint)
    );
}

#[test]
fn birational_open_transport_rejects_montgomery_infinity() {
    let montgomery = MontgomeryCurve::<F5>::new(F5::from_i64(4), F5::one())
        .expect("sample Montgomery curve should be non-singular");

    assert_eq!(
        montgomery.try_point_to_twisted_edwards_open(&AffinePoint::infinity()),
        Err(TwistedEdwardsBirationalMapError::ExceptionalMontgomeryPoint)
    );
}

#[test]
fn birational_open_transport_rejects_montgomery_points_with_y_zero() {
    let montgomery = MontgomeryCurve::<F5>::new(F5::from_i64(4), F5::one())
        .expect("sample Montgomery curve should be non-singular");
    let point = montgomery
        .point(F5::zero(), F5::zero())
        .expect("sample y = 0 point should lie on the Montgomery curve");

    assert_eq!(
        montgomery.try_point_to_twisted_edwards_open(&point),
        Err(TwistedEdwardsBirationalMapError::ExceptionalMontgomeryPoint)
    );
}

#[test]
fn birational_open_transport_rejects_montgomery_points_with_x_minus_one() {
    type F7 = Fp<7>;

    let montgomery = MontgomeryCurve::<F7>::new(F7::from_i64(3), F7::one())
        .expect("sample Montgomery curve should be non-singular");
    let point = montgomery
        .point(F7::from_i64(-1), F7::one())
        .expect("sample x = -1 point should lie on the Montgomery curve");

    assert_eq!(
        montgomery.try_point_to_twisted_edwards_open(&point),
        Err(TwistedEdwardsBirationalMapError::ExceptionalMontgomeryPoint)
    );
}
