use super::shared::{F5, F7, f5_curve};
use crate::elliptic_curves::{
    AffinePoint, MontgomeryCurve,
    models::twisted_edwards::TwistedEdwardsBirationalMapError,
    traits::{AffineCurveModel, CurveModel},
};
use crate::fields::traits::Field;

#[test]
fn birational_open_transport_to_montgomery_maps_a_regular_affine_point() {
    let twisted = f5_curve();
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
    let twisted = f5_curve();
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
    let twisted = f5_curve();

    assert_eq!(
        twisted.try_point_to_montgomery_open(&twisted.identity()),
        Err(TwistedEdwardsBirationalMapError::ExceptionalTwistedEdwardsPoint)
    );
}

#[test]
fn birational_open_transport_rejects_the_other_x_zero_twisted_edwards_point() {
    let twisted = f5_curve();
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
