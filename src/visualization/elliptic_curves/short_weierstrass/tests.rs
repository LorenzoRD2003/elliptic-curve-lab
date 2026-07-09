use num_bigint::{BigInt, BigUint};
use num_rational::BigRational;

use super::*;
use crate::elliptic_curves::{
    AffinePoint,
    frobenius::group_order::SmallFieldGroupOrderStrategy,
    short_weierstrass::{
        group_exponent::{GroupExponentReport, GroupExponentStrategy},
        point_order::{PointOrderReport, PointOrderStrategy},
        rational_torsion::RationalTorsionStrategy,
    },
    traits::{AffineCurveModel, EnumerableCurveModel},
};
use crate::fields::Q;
use crate::visualization::Visualizable;

type F7 = crate::fields::Fp7;

fn q(numerator: i64, denominator: i64) -> BigRational {
    BigRational::new(BigInt::from(numerator), BigInt::from(denominator))
}

fn bu(value: u64) -> BigUint {
    BigUint::from(value)
}

fn f7_curve() -> crate::elliptic_curves::ShortWeierstrassCurve<F7> {
    crate::elliptic_curves::ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3))
        .expect("valid curve")
}

fn f7_point(x: i64, y: i64) -> AffinePoint<F7> {
    f7_curve()
        .point(F7::from_i64(x), F7::from_i64(y))
        .expect("point should lie on the curve")
}

#[test]
fn curve_display_and_equation_string_share_one_equation_surface() {
    let curve = f7_curve();

    assert_eq!(
        curve.to_equation_string(),
        "y^2 = x^3 + (2 (mod 7))x + (3 (mod 7))"
    );
    assert_eq!(format!("{curve}"), curve.to_equation_string());
    assert_eq!(format_curve(&curve), "y^2 = x^3 + 2x + 3");
}

#[test]
fn describe_rational_torsion_report_mentions_scaled_integral_model() {
    let curve = crate::elliptic_curves::ShortWeierstrassCurve::<Q>::new(q(-1, 16), q(0, 1))
        .expect("valid rational curve");
    let report = curve
        .rational_torsion_by(RationalTorsionStrategy::LutzNagell)
        .expect("scaled curve should have certified rational torsion");

    let description = describe_rational_torsion_report(&report);
    assert!(description.contains("Rational torsion over Q"));
    assert!(description.contains("integral transport: source curve was scaled"));
    assert!(description.contains("group: ℤ/2ℤ × ℤ/2ℤ"));
    assert!(description.contains("torsion points:"));
}

#[test]
fn point_display_uses_affine_coordinates_or_identity_symbol() {
    let point = f7_point(2, 1);
    let infinity = AffinePoint::<F7>::infinity();

    assert_eq!(point.to_coordinates_string(), "(2 (mod 7), 1 (mod 7))");
    assert_eq!(format!("{point}"), point.to_coordinates_string());
    assert_eq!(format_point(&point), point.to_coordinates_string());
    assert_eq!(format_point_compact(&point), "(2, 1)");
    assert_eq!(format_point(&infinity), "O");
    assert_eq!(format_point_compact(&infinity), "O");
}

#[test]
fn debug_output_is_more_informative_than_the_default_derives() {
    let curve = f7_curve();
    let point = f7_point(2, 1);

    assert!(format!("{curve:?}").contains("ShortWeierstrassCurve"));
    assert!(format!("{curve:?}").contains("equation"));
    assert!(format!("{point:?}").contains("AffinePoint"));
    assert!(format!("{point:?}").contains("x"));
}

#[test]
fn curve_description_mentions_invariants() {
    let description = describe_curve(&f7_curve());

    assert!(description.contains("Short-Weierstrass curve"));
    assert!(description.contains("discriminant"));
    assert!(description.contains("j-invariant"));
}

#[test]
fn point_description_mentions_identity_and_membership_status() {
    let description = describe_point(&f7_curve(), &f7_point(2, 1));

    assert!(description.contains("Curve point"));
    assert!(description.contains("identity: no"));
    assert!(description.contains("on curve: yes"));
}

#[test]
fn membership_description_shows_both_sides_of_the_equation() {
    let description = describe_membership(&f7_curve(), &f7_point(2, 1));

    assert!(description.contains("left side: y^2"));
    assert!(description.contains("right side: x^3 + ax + b"));
    assert!(description.contains("result: on curve"));
}

#[test]
fn membership_description_is_honest_about_the_point_at_infinity() {
    let description = describe_membership(&f7_curve(), &AffinePoint::<F7>::infinity());

    assert!(description.contains("point: O"));
    assert!(description.contains("convention"));
}

#[test]
fn addition_explanation_mentions_the_geometric_case_and_result() {
    let explanation =
        explain_add(&f7_curve(), &f7_point(2, 1), &f7_point(3, 1)).expect("valid addition");

    assert!(explanation.contains("Point addition"));
    assert!(explanation.contains("case: secant formula"));
    assert!(explanation.contains("result: (2 (mod 7), 6 (mod 7))"));
}

#[test]
fn point_listing_shows_group_order_and_identity() {
    let listing = list_points(&f7_curve());

    assert!(listing.contains("Curve points"));
    assert!(listing.contains("group order: 6"));
    assert!(listing.contains("0: O"));
}

#[test]
fn point_order_description_mentions_repeated_addition_method() {
    let description = describe_point_order(&f7_curve(), &f7_point(2, 1));

    assert!(description.contains("Point order"));
    assert!(description.contains("repeated addition"));
    assert!(description.contains("point order: 6"));
}

#[test]
fn point_order_description_is_honest_about_invalid_points() {
    let description = describe_point_order(
        &f7_curve(),
        &AffinePoint::<F7>::new(F7::from_i64(2), F7::from_i64(2)),
    );

    assert!(description.contains("result: point is not on the curve"));
}

#[test]
fn group_structure_description_reports_small_cyclic_example() {
    let description = describe_group_structure(&f7_curve());

    assert!(description.contains("Finite curve group structure"));
    assert!(description.contains("group order: 6"));
    assert!(description.contains("cyclic: yes"));
    assert!(description.contains("exponent: 6"));
    assert!(description.contains("invariant factors: Z/6Z"));
}

#[test]
fn compact_group_structure_summary_reports_core_invariants() {
    let summary = summarize_group_structure(&f7_curve());

    assert!(summary.contains("cyclic: yes"));
    assert!(summary.contains("exponent: 6"));
    assert!(summary.contains("invariant factors: Z/6Z"));
}

#[test]
fn order_distribution_description_lists_exact_point_orders() {
    let description = describe_order_distribution(&f7_curve());

    assert!(description.contains("Point-order distribution"));
    assert!(description.contains("order 1: 1 point(s)"));
    assert!(description.contains("order 2: 1 point(s)"));
    assert!(description.contains("order 3: 2 point(s)"));
    assert!(description.contains("order 6: 2 point(s)"));
}

#[test]
fn compact_order_distribution_summary_uses_arrow_surface() {
    let summary = summarize_order_distribution(&f7_curve());

    assert!(summary.contains("1 -> 1"));
    assert!(summary.contains("2 -> 1"));
    assert!(summary.contains("3 -> 2"));
    assert!(summary.contains("6 -> 2"));
}

#[test]
fn scalar_multiplication_description_reports_method_and_result() {
    let description =
        describe_scalar_mul(&f7_curve(), &f7_point(2, 1), 3).expect("valid scalar multiply");

    assert!(description.contains("Scalar multiplication"));
    assert!(description.contains("scalar: 3"));
    assert!(description.contains("double-and-add"));
    assert!(description.contains("result: [3]P = (6 (mod 7), 0 (mod 7))"));
}

#[test]
fn point_order_explanation_lists_successive_multiples_until_identity() {
    let description = explain_point_order(&f7_curve(), &f7_point(2, 1));

    assert!(description.contains("Point-order explanation"));
    assert!(description.contains("[1]P = (2 (mod 7), 1 (mod 7))"));
    assert!(description.contains("[6]P = O"));
    assert!(description.contains("first identity hit: [6]P = O"));
    assert!(description.contains("point order: 6"));
}

#[test]
fn point_order_from_multiple_visualization_reports_the_prime_peeling_steps() {
    let report = f7_curve()
        .point_order_from_multiple(&f7_point(6, 0), bu(6), &[(bu(2), 1), (bu(3), 1)])
        .expect("valid reduction report should build");

    assert_eq!(
        format_point_order_from_multiple_report(&report),
        "ord(P) from M = 6 is 2"
    );

    let description = describe_point_order_from_multiple_report(&report);
    assert!(description.contains("Point order from multiple"));
    assert!(description.contains("supplied multiple M: 6"));
    assert!(description.contains("exact order recovered: 2"));
    assert!(
        description
            .contains("prime 3: exponent in M = 1, removed exponent = 1, remaining multiple = 2")
    );
}

#[test]
fn unified_point_order_visualization_mentions_the_selected_strategy() {
    let report = f7_curve()
        .point_order_by(
            &f7_point(2, 1),
            PointOrderStrategy::HasseIntervalNaive {
                group_order_strategy: SmallFieldGroupOrderStrategy::Auto,
            },
        )
        .expect("Hasse-interval order recovery should succeed");

    assert_eq!(
        format_point_order_report(&report),
        "ord(P) via H(q) search = 6"
    );

    let description = describe_point_order_report(&report);
    assert!(description.contains("Point order report"));
    assert!(description.contains("strategy: naive Hasse interval"));
    assert!(description.contains("exact order: 6"));
    assert!(description.contains("group-order route: quadratic character"));
    assert!(description.contains("first H(q)-multiple annihilating P: 6"));
}

#[test]
fn exhaustive_point_order_visualization_stays_honest_about_the_route() {
    let report = f7_curve()
        .point_order_by(&f7_point(2, 1), PointOrderStrategy::Exhaustive)
        .expect("exhaustive order recovery should succeed");

    let PointOrderReport::Exhaustive(exhaustive) = report else {
        panic!("expected the exhaustive route to preserve its variant");
    };

    assert_eq!(
        describe_exhaustive_point_order_report(&exhaustive),
        exhaustive.describe()
    );
    assert!(exhaustive.describe().contains("Exhaustive point order"));
    assert!(exhaustive.describe().contains("exact order: 6"));
}

#[test]
fn group_exponent_visualization_mentions_the_selected_strategy() {
    let curve = f7_curve();
    let sampled_point = f7_point(2, 1);
    let point_index = curve
        .points()
        .iter()
        .position(|candidate| candidate == &sampled_point)
        .expect("sample point should appear in the enumerated group");
    let mut sampler = move |upper_bound: usize| (point_index < upper_bound).then_some(point_index);

    let report = curve
        .group_exponent_by(
            GroupExponentStrategy::RandomPoints {
                max_samples: 1,
                point_order_strategy: PointOrderStrategy::HasseIntervalNaive {
                    group_order_strategy: SmallFieldGroupOrderStrategy::Auto,
                },
            },
            &mut sampler,
        )
        .expect("random-point exponent accumulation should succeed");

    assert_eq!(
        format_group_exponent_report(&report),
        "group exponent lower bound after 1 sample(s) = 6"
    );

    let description = describe_group_exponent_report(&report);
    assert!(description.contains("Group exponent report"));
    assert!(description.contains("strategy: random points"));
    assert!(description.contains("exponent lower bound: 6"));
    assert!(description.contains("exact exponent: not certified by this route"));
    assert!(description.contains("point-order route: naive Hasse interval"));
}

#[test]
fn exhaustive_group_exponent_visualization_stays_honest_about_exactness() {
    let curve = f7_curve();
    let mut sampler = |_| Some(0usize);
    let report = curve
        .group_exponent_by(GroupExponentStrategy::Exhaustive, &mut sampler)
        .expect("exhaustive exponent route should succeed");

    let GroupExponentReport::Exhaustive(exact_exponent) = report else {
        panic!("expected the exhaustive group-exponent route to preserve its variant");
    };

    assert_eq!(
        describe_exhaustive_group_exponent_report(&exact_exponent),
        "Exhaustive group exponent\nexact exponent: 6\nstrategy: compute every point order in the tiny ambient group and take their lcm"
    );
    assert!(
        describe_exhaustive_group_exponent_report(&exact_exponent)
            .contains("Exhaustive group exponent")
    );
    assert!(
        describe_exhaustive_group_exponent_report(&exact_exponent).contains("exact exponent: 6")
    );
}

#[test]
fn exponent_lower_bound_group_order_verification_visualization_stays_honest_about_scope() {
    let curve = crate::elliptic_curves::ShortWeierstrassCurve::<crate::fields::Fp5>::new(
        crate::fields::Fp5::from_i64(0),
        crate::fields::Fp5::from_i64(1),
    )
    .expect("valid curve");
    let sampled_point = curve
        .point(
            crate::fields::Fp5::from_i64(2),
            crate::fields::Fp5::from_i64(2),
        )
        .expect("point should lie on the curve");
    let point_index = curve
        .points()
        .iter()
        .position(|candidate| candidate == &sampled_point)
        .expect("sample point should appear in the enumerated group");
    let mut sampler = move |upper_bound: usize| (point_index < upper_bound).then_some(point_index);

    let report = curve
        .group_exponent_by(
            GroupExponentStrategy::RandomPoints {
                max_samples: 1,
                point_order_strategy: PointOrderStrategy::Exhaustive,
            },
            &mut sampler,
        )
        .expect("random-point exponent accumulation should succeed");
    let GroupExponentReport::RandomPoints(accumulation) = report else {
        panic!("expected accumulation report");
    };
    let verification = curve
        .verify_exponent_lower_bound_by_group_order(
            &accumulation,
            SmallFieldGroupOrderStrategy::Auto,
        )
        .expect("verification should succeed");

    assert_eq!(
        format_exponent_lower_bound_group_order_verification(&verification),
        "group order verifies #E(F_q) = 6 from lower bound 6"
    );

    let description = describe_exponent_lower_bound_group_order_verification(&verification);
    assert!(description.contains("Exponent lower-bound verification by group order"));
    assert!(description.contains("exponent lower bound: 6"));
    assert!(description.contains("verified group order: 6"));
    assert!(description.contains("does not by itself certify the exponent"));
}

#[test]
fn visualizable_trait_is_hooked_up_for_curves_and_points() {
    let curve = f7_curve();
    let point = f7_point(2, 1);

    assert!(curve.describe().contains("Short-Weierstrass curve"));
    assert_eq!(point.format_compact(), format_point_compact(&point));
}

#[test]
fn curve_display_works_over_q_too() {
    let curve = crate::elliptic_curves::ShortWeierstrassCurve::<Q>::new(q(-1, 1), q(0, 1))
        .expect("valid curve");

    assert_eq!(curve.to_equation_string(), "y^2 = x^3 + (-1)x + (0)");
    assert_eq!(format!("{curve}"), curve.to_equation_string());
    assert_eq!(format_curve(&curve), "y^2 = x^3 + (-1)x");
}
