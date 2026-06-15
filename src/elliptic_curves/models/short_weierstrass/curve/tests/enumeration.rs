use crate::elliptic_curves::traits::{CurveModel, EnumerableCurveModel};

use super::shared::f7_curve;

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
