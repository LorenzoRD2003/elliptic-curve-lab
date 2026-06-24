use super::shared::{F5, f5_curve};
use crate::elliptic_curves::traits::{
    AffineCurveModel, CurveModel, EnumerableCurveModel, FiniteGroupCurveModel,
};
use crate::fields::traits::Field;

#[test]
fn twisted_edwards_small_field_group_axioms_hold_on_a_tiny_curve() {
    let curve = f5_curve();

    assert_eq!(curve.check_group_axioms(), Ok(()));
}

#[test]
fn enumerable_curve_model_lists_the_finite_identity_exactly_once() {
    let curve = f5_curve();
    let identity = curve.identity();
    let finite_points = curve.finite_points();
    let all_points = curve.points();

    assert!(!finite_points.iter().any(|point| curve.is_identity(point)));
    assert_eq!(all_points.first(), Some(&identity));
    assert_eq!(
        all_points
            .iter()
            .filter(|point| curve.is_identity(point))
            .count(),
        1
    );
    assert_eq!(all_points.len(), 8);
}

#[test]
fn finite_group_surfaces_match_the_montgomery_companion_on_a_small_curve() {
    let curve = f5_curve();
    let montgomery = curve.as_montgomery();

    assert_eq!(curve.order(), montgomery.order());
    assert_eq!(curve.exponent(), montgomery.exponent());
    assert_eq!(curve.order_distribution(), montgomery.order_distribution());
    assert_eq!(curve.group_structure(), montgomery.group_structure());
    assert_eq!(
        curve.describe_group_structure(),
        montgomery.describe_group_structure()
    );
}

#[test]
fn finite_group_curve_model_reports_the_expected_structure_on_a_small_curve() {
    let curve = f5_curve();
    let structure = curve.group_structure();

    assert_eq!(curve.order(), 8);
    assert_eq!(curve.exponent(), 8);
    assert_eq!(curve.order_distribution().get(&1), Some(&1));
    assert_eq!(curve.order_distribution().get(&2), Some(&1));
    assert_eq!(curve.order_distribution().get(&4), Some(&2));
    assert_eq!(curve.order_distribution().get(&8), Some(&4));
    assert_eq!(structure.order, 8);
    assert_eq!(structure.exponent, 8);
    assert!(structure.cyclic);
    assert_eq!(structure.invariant_factors, None);
}

#[test]
fn point_order_matches_the_montgomery_companion_on_the_common_open_subset() {
    let curve = f5_curve();
    let point = curve
        .point(F5::one(), F5::zero())
        .expect("sample point should lie on the twisted-Edwards curve");
    let montgomery = curve.as_montgomery();
    let montgomery_point = curve
        .try_point_to_montgomery_open(&point)
        .expect("sample point should lie in the common birational open");

    assert_eq!(
        curve.point_order(&point),
        montgomery.point_order(&montgomery_point)
    );
}
