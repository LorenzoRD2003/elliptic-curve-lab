use num_bigint::BigUint;
use std::collections::BTreeMap;

use crate::elliptic_curves::frobenius::group_order::{
    GroupOrderRoute, SmallFieldGroupOrderStrategy,
};
use crate::elliptic_curves::traits::{
    CurveModel, EnumerableCurveModel, FiniteAbelianGroupStructure, FiniteGroupCurveModel,
    FrobeniusTraceCurveModel, GroupCurveModel,
};
use crate::elliptic_curves::{AffinePoint, ShortWeierstrassCurve};

use super::shared::{F7, f7_curve, f7_point, f43_curve, small_cyclic_distribution};

#[test]
fn public_group_order_api_prefers_character_sum_in_auto_mode() {
    let curve = f43_curve();

    let report = curve
        .group_order_by_small_field(SmallFieldGroupOrderStrategy::Auto)
        .expect("automatic group order should succeed");

    assert_eq!(report.route(), GroupOrderRoute::QuadraticCharacter);
    assert_eq!(report.curve_order(), BigUint::from(curve.order() as u64));
}

#[test]
fn public_frobenius_trace_by_agrees_with_the_exhaustive_trace() {
    let curve = f43_curve();

    assert_eq!(
        curve.frobenius_trace_by_small_field(SmallFieldGroupOrderStrategy::Exhaustive),
        curve.frobenius_trace()
    );
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
    let curve =
        ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3)).expect("valid curve");
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
    let expected: BTreeMap<usize, usize> = small_cyclic_distribution();

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
