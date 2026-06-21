use super::shared::{f3_curve, f5_curve};
use crate::elliptic_curves::traits::{CurveModel, EnumerableCurveModel};

#[test]
fn finite_point_enumeration_lists_exactly_the_non_identity_points() {
    let curve = f5_curve();
    let finite_points = curve.finite_points();

    assert_eq!(finite_points.len(), 7);
    assert!(
        finite_points
            .iter()
            .all(|point| curve.is_on_curve_nonzero(point))
    );
}

#[test]
fn full_point_enumeration_includes_identity_and_order() {
    let curve = f5_curve();
    let points = curve.points();

    assert_eq!(points.len(), 8);
    assert!(curve.is_identity(points.first().expect("identity should be present")));
    assert_eq!(curve.order(), 8);
}

#[test]
fn characteristic_three_curve_is_enumerable() {
    let curve = f3_curve();

    assert_eq!(curve.finite_points().len(), 3);
    assert_eq!(curve.order(), 4);
}

#[test]
fn random_point_uses_the_supplied_index_sampler() {
    let curve = f5_curve();
    let expected = curve.points()[2].clone();
    let mut sampler = |upper_bound: usize| {
        assert_eq!(upper_bound, 8);
        Some(2)
    };

    let sampled = curve
        .random_point(&mut sampler)
        .expect("sampler should choose an existing point");

    assert_eq!(sampled, expected);
}

#[test]
fn random_point_propagates_sampler_failure() {
    let curve = f5_curve();
    let mut sampler = |_upper_bound: usize| None;

    assert!(curve.random_point(&mut sampler).is_none());
}
