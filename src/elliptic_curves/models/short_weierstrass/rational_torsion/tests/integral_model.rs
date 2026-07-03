use crate::elliptic_curves::{
    AffinePoint, ShortWeierstrassCurve,
    short_weierstrass::rational_torsion::integral_model::RationalIntegralModel, traits::CurveModel,
};
use crate::fields::Q;

use super::fixtures::{product_two_two_fixture, q, rational_scaled_fixture};

#[test]
fn integral_model_witness_stores_curve_and_scale_for_later_transport() {
    let fixture = rational_scaled_fixture();
    let witness =
        RationalIntegralModel::new(fixture.curve, q(2, 1)).expect("u = 2 should be invertible");

    assert_eq!(witness.scale(), &q(2, 1));
    assert_eq!(witness.source_curve().a(), &q(-1, 16));
    assert_eq!(witness.curve().a(), &q(-1, 1));
}

#[test]
fn integral_model_uses_unit_scale_for_integral_coefficients() {
    let fixture = product_two_two_fixture();
    let witness = RationalIntegralModel::from_curve(fixture.curve)
        .expect("integral fixture should already be an integral model");

    assert_eq!(witness.scale(), &q(1, 1));
    assert!(witness.has_integral_coefficients());
    assert_eq!(witness.curve().a(), &q(-1, 1));
    assert_eq!(witness.curve().b(), &q(0, 1));
}

#[test]
fn integral_model_computes_integer_scale_from_denominator_powers() {
    let fixture = rational_scaled_fixture();
    let witness = RationalIntegralModel::from_curve(fixture.curve)
        .expect("rational fixture should scale to an integral model");

    assert_eq!(witness.scale(), &q(2, 1));
    assert!(witness.has_integral_coefficients());
    assert_eq!(witness.curve().a(), &q(-1, 1));
    assert_eq!(witness.curve().b(), &q(0, 1));
}

#[test]
fn integral_model_transports_points_roundtrip() {
    let source_curve =
        ShortWeierstrassCurve::<Q>::new(q(0, 1), q(1, 64)).expect("valid rational source curve");
    let source_point = AffinePoint::new(q(0, 1), q(1, 8));
    let witness = RationalIntegralModel::from_curve(source_curve.clone())
        .expect("curve should scale to y² = x³ + 1");

    let integral_point = witness
        .to_integral_point(&source_point)
        .expect("source point should transport to integral model");
    assert_eq!(integral_point, AffinePoint::new(q(0, 1), q(1, 1)));
    assert!(witness.curve().contains(&integral_point));

    let roundtrip = witness
        .to_source_point(&integral_point)
        .expect("integral point should transport back to the source");
    assert_eq!(roundtrip, source_point);
    assert!(source_curve.contains(&roundtrip));
}
