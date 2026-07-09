use num_bigint::BigUint;

use super::*;
use crate::elliptic_curves::{
    ShortWeierstrassCurve,
    frobenius::{
        AbsoluteFrobenius, FrobeniusTrace, RelativeFrobenius,
        characteristic_equation::FrobeniusCharacteristicEquationCurveModel,
        extension_counts::compare_extension_count_with_enumeration,
        group_order::{
            GroupOrderReport, MestreConfig, MestreGroupOrderReport, MestreSide, MestreStepReport,
            SmallFieldGroupOrderStrategy,
        },
        orbit::relative_frobenius_orbit,
    },
    short_weierstrass::isomorphisms::ShortWeierstrassQuadraticTwist,
    traits::{AffineCurveModel, FiniteGroupCurveModel, FrobeniusTraceCurveModel},
};
use crate::fields::traits::{EnumerableFiniteField, Field, SqrtField};
use crate::isogenies::{
    frobenius_relation::{FrobeniusComparableIsogeny, FrobeniusComparableIsogenyGraph},
    graphs::IsogenyGraphBuilder,
    scalar_multiplication::ScalarMultiplicationIsogeny,
};
use crate::proptest_support::fields::ProptestF17Sqrt3Field;
use crate::visualization::traits::Visualizable;

type F17 = crate::fields::Fp17;
type F19 = crate::fields::Fp19;
type F7 = crate::fields::Fp7;
type F41 = crate::fields::Fp41;
type F43 = crate::fields::Fp43;
type F17Squared = ProptestF17Sqrt3Field;

fn f41_curve() -> ShortWeierstrassCurve<F41> {
    ShortWeierstrassCurve::<F41>::new(F41::from_i64(2), F41::from_i64(3)).expect("valid F41 curve")
}

fn lift_f17_curve_to_f17_squared(
    curve: &ShortWeierstrassCurve<F17>,
) -> ShortWeierstrassCurve<F17Squared> {
    ShortWeierstrassCurve::<F17Squared>::new(
        F17Squared::from_base(*curve.a()),
        F17Squared::from_base(*curve.b()),
    )
    .expect("lifting a smooth F17 curve should preserve smoothness")
}

fn first_nonsquare<F>() -> F::Elem
where
    F: EnumerableFiniteField + SqrtField,
{
    F::elements()
        .into_iter()
        .find(|value| !F::is_zero(value) && !F::has_square_root(value))
        .expect("small odd prime fields should contain non-squares")
}

fn first_non_fixed_point<F>(
    curve: &ShortWeierstrassCurve<F>,
) -> crate::elliptic_curves::AffinePoint<F>
where
    F: EnumerableFiniteField + SqrtField,
    F::Elem: Clone,
{
    for x in F::elements() {
        for y in F::elements() {
            if let Ok(point) = curve.point(x.clone(), y) {
                let image = curve
                    .absolute_frobenius_power_point(&point, 1)
                    .expect("absolute Frobenius should evaluate");
                if image != point {
                    return point;
                }
            }
        }
    }

    panic!("expected a non-fixed point over the quadratic extension")
}

#[test]
fn frobenius_metadata_visualizations_keep_absolute_and_relative_distinct() {
    let absolute = AbsoluteFrobenius::for_field::<F43>(3);
    let relative = RelativeFrobenius::for_field::<F17Squared>(2);

    assert_eq!(format_absolute_frobenius(&absolute), "π_43^3");
    assert_eq!(format_relative_frobenius(&relative), "π_(17^2)^2");

    let absolute_description = describe_absolute_frobenius(&absolute);
    let relative_description = describe_relative_frobenius(&relative);

    assert!(absolute_description.contains("Absolute Frobenius"));
    assert!(absolute_description.contains("characteristic p: 43"));
    assert!(relative_description.contains("Relative Frobenius"));
    assert!(relative_description.contains("base field: F_(17^2)"));
    assert!(relative_description.contains("field order q: 289"));
}

#[test]
fn trace_polynomial_and_zeta_visualizations_share_the_same_frobenius_story() {
    let curve = ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid curve");
    let trace = curve.frobenius_trace().expect("trace should compute");
    let polynomial = trace.characteristic_polynomial();
    let zeta = trace.local_zeta_function();

    assert!(format_frobenius_trace(&trace).contains("t ="));

    let trace_description = describe_frobenius_trace(&trace);
    let polynomial_description = describe_frobenius_characteristic_polynomial(&polynomial);
    let zeta_description = describe_frobenius_local_zeta_function(&zeta);

    assert!(trace_description.contains("curve order #E(F_q)"));
    assert!(polynomial_description.contains("χ_π(T)"));
    assert!(polynomial_description.contains("discriminant t^2 - 4q"));
    assert!(zeta_description.contains("source characteristic polynomial"));
    assert!(zeta_description.contains("numerator:"));
    assert!(zeta_description.contains("denominator:"));
    assert_eq!(zeta.format_compact(), zeta.pretty());
}

#[test]
fn hasse_and_curve_type_visualizations_explain_their_exact_criteria() {
    let ordinary_curve =
        ShortWeierstrassCurve::<F43>::new(F43::zero(), F43::one()).expect("valid curve");
    let _ordinary_trace = ordinary_curve
        .frobenius_trace()
        .expect("trace should compute");
    let hasse_report = ordinary_curve
        .verify_hasse_bound()
        .expect("Hasse report should compute");

    let hasse_description = describe_hasse_bound_report(&hasse_report);

    assert!(hasse_description.contains("trace square t^2"));
    assert!(hasse_description.contains("bound square 4q"));
    assert!(hasse_description.contains("slack 4q - t^2"));
}

#[test]
fn hasse_interval_visualization_reports_discrete_search_data() {
    let interval = crate::elliptic_curves::frobenius::HasseInterval::for_q(BigUint::from(43u8))
        .expect("q = 43 should define a Hasse interval");

    assert_eq!(format_hasse_interval(&interval), "H(43) = [31 , 57]");

    let description = describe_hasse_interval(&interval);
    assert!(description.contains("field order q: 43"));
    assert!(description.contains("interval H(q): [31 , 57]"));
    assert!(description.contains("integer candidate count: 27"));
    assert!(description.contains("floor(sqrt(4q)): 13"));
}

#[test]
fn naive_hasse_multiple_search_visualization_reports_the_first_hit_and_steps() {
    let curve =
        ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3)).expect("valid curve");
    let point = curve
        .generator()
        .expect("the sample curve should be cyclic");
    let report = curve
        .find_annihilating_multiple_in_hasse_interval_naive(&point)
        .expect("naive Hasse search should succeed");

    assert_eq!(
        format_hasse_multiple_search_report(&report),
        "first H(q)-multiple annihilating P: 6"
    );

    let description = describe_hasse_multiple_search_report(&report);
    assert!(description.contains("searched interval: H(7) = [3 , 13]"));
    assert!(description.contains("tested candidates: 4"));
    assert!(description.contains("first annihilating multiple: 6"));
    assert!(description.contains("M = 6 gives [M]P = O"));
}

#[test]
fn character_sum_visualization_reports_the_counting_formula() {
    let curve = ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid curve");
    let report = curve
        .group_order_by_quadratic_character()
        .expect("character-sum count should succeed");

    assert_eq!(
        format_character_sum_point_count(&report),
        "#E(F_43) via χ-sum = 34"
    );

    let description = describe_character_sum_point_count(&report);
    assert!(description.contains("Quadratic-character point count"));
    assert!(description.contains("character sum Σ χ(f(x))"));
    assert!(description.contains("counting formula: #E(F_q) = q + 1 + Σ χ(f(x))"));
}

#[test]
fn unified_group_order_visualization_mentions_the_strategy() {
    let curve = ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid curve");
    let report = curve
        .group_order_by_small_field(SmallFieldGroupOrderStrategy::Auto)
        .expect("automatic group order should succeed");

    assert_eq!(
        format_group_order_report(&report),
        "#E(F_43) via χ-sum = 34"
    );

    let description = describe_group_order_report(&report);
    assert!(description.contains("Group order"));
    assert!(description.contains("strategy: quadratic character"));
    assert!(description.contains("curve order #E(F_q): 34"));
}

#[test]
fn mestre_visualizations_show_side_history_and_lower_bounds() {
    let base_field = crate::fields::finite_field_descriptor::FiniteFieldDescriptor::new(
        43,
        core::num::NonZeroU32::new(1).expect("1 is non-zero"),
    )
    .expect("prime field descriptor should build");
    let original = FrobeniusTrace::from_order(base_field.clone(), 52)
        .expect("original Frobenius package should build");
    let twist = FrobeniusTrace::from_order(base_field, 36).expect("twist package should build");
    let point_order_report = ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3))
        .expect("valid sample curve")
        .point_order_from_multiple(
            &ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3))
                .expect("valid sample curve")
                .point(F7::from_i64(2), F7::from_i64(1))
                .expect("sample point should lie on the curve"),
            BigUint::from(6u8),
            &[(BigUint::from(2u8), 1), (BigUint::from(3u8), 1)],
        )
        .expect("known-multiple route should recover a sample order");
    let step = MestreStepReport::new(
        MestreSide::QuadraticTwist,
        BigUint::from(45u8),
        point_order_report,
        BigUint::from(9u8),
    );
    let mestre = MestreGroupOrderReport::new(
        MestreConfig::with_iteration_cap(8),
        MestreSide::QuadraticTwist,
        original,
        twist,
        vec![step.clone()],
    );

    assert_eq!(
        format_mestre_step_report(&step),
        "quadratic twist: M = 45, ord(P) = 6, running λ lower bound = 9"
    );
    assert_eq!(
        format_mestre_group_order_report(&mestre),
        "#E(F_43) via Mestre = 52"
    );

    let step_description = describe_mestre_step_report(&step);
    let mestre_description = describe_mestre_group_order_report(&mestre);
    let unified_description =
        describe_group_order_report(&GroupOrderReport::MestreFp(Box::new(mestre.clone())));

    assert!(step_description.contains("side: quadratic twist"));
    assert!(step_description.contains("annihilating multiple in H(p): 45"));
    assert!(mestre_description.contains("resolved side: quadratic twist"));
    assert!(mestre_description.contains("shared Hasse interval: H(43)"));
    assert!(mestre_description.contains("step 1: quadratic twist: M = 45"));
    assert!(mestre_description.contains("iteration cap: 8"));
    assert!(unified_description.contains("strategy: Mestre"));
    assert!(unified_description.contains("note: the returned group order is always #E(F_p)"));
}

#[test]
fn extension_count_visualizations_show_the_derived_and_exhaustive_routes() {
    let curve = ShortWeierstrassCurve::<F17>::new(F17::from_i64(2), F17::from_i64(3))
        .expect("valid F17 curve");
    let lifted_curve = lift_f17_curve_to_f17_squared(&curve);
    let trace = curve.frobenius_trace().expect("trace should compute");
    let report =
        trace.curve_order_over_extension(core::num::NonZeroU32::new(2).expect("2 is positive"));
    let sequence = trace.curve_orders_over_extensions_through(
        core::num::NonZeroU32::new(3).expect("3 is positive"),
    );
    let comparison = compare_extension_count_with_enumeration(&lifted_curve, &trace)
        .expect("comparison should compute");

    let report_description = describe_frobenius_extension_count_report(&report);
    let sequence_description = describe_frobenius_extension_count_sequence_report(&sequence);
    let comparison_description =
        describe_frobenius_extension_enumeration_comparison_report(&comparison);

    assert!(report_description.contains("extension field: F_(17^2)"));
    assert!(report_description.contains("power sum s_n = α^n + β^n"));
    assert!(sequence_description.contains("degree 3"));
    assert!(comparison_description.contains("Frobenius-derived count"));
    assert!(comparison_description.contains("exhaustive count"));
    assert!(comparison_description.contains("agreement: yes"));
}

#[test]
fn characteristic_equation_visualizations_show_all_pointwise_terms() {
    let curve = ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid curve");
    let point = curve
        .point(F43::zero(), F43::one())
        .expect("sample point should lie on the curve");
    let polynomial = curve
        .frobenius_trace()
        .expect("trace should compute")
        .characteristic_polynomial();
    let check = curve
        .verify_frobenius_characteristic_equation_at_point(&point, &polynomial)
        .expect("pointwise check should compute");
    let exhaustive = curve
        .verify_frobenius_characteristic_equation_exhaustive()
        .expect("exhaustive report should compute");

    let check_description = describe_frobenius_characteristic_equation_check(&check);
    let exhaustive_description =
        describe_frobenius_characteristic_equation_exhaustive_report(&exhaustive);

    assert!(check_description.contains("π_q(P)"));
    assert!(check_description.contains("π_q^2(P)"));
    assert!(check_description.contains("[t]π_q(P)"));
    assert!(check_description.contains("[q]P"));
    assert!(exhaustive_description.contains("checked points"));
    assert!(exhaustive_description.contains("failed checks: 0"));
    assert!(exhaustive_description.contains("failed points: none"));
}

#[test]
fn orbit_and_torsion_visualizations_report_motion_and_periods() {
    let curve = ShortWeierstrassCurve::<F17Squared>::new(
        F17Squared::from_base(F17::zero()),
        F17Squared::from_base(F17::one()),
    )
    .expect("valid extension curve");
    let point = first_non_fixed_point(&curve);
    let orbit = curve
        .absolute_frobenius_orbit(&point, 1)
        .expect("orbit should compute");
    let torsion_report = curve
        .absolute_frobenius_on_exact_torsion(4, 1)
        .expect("torsion report should compute");

    let orbit_description = describe_frobenius_orbit(&orbit);
    let torsion_description = describe_frobenius_on_exact_torsion_report(&torsion_report);

    assert!(orbit_description.contains("period: 2"));
    assert!(orbit_description.contains("points: ["));
    assert!(torsion_description.contains("exact order n: 4"));
    assert!(torsion_description.contains("moved count:"));
    assert!(torsion_description.contains("orbit periods: ["));
    assert!(torsion_description.contains("minimal absolute-Frobenius fixing powers:"));
}

#[test]
fn quadratic_twist_and_isogeny_visualizations_report_their_invariants() {
    let twist_curve =
        ShortWeierstrassCurve::<F19>::new(F19::from_i64(2), F19::from_i64(3)).expect("valid curve");
    let twist_package = ShortWeierstrassQuadraticTwist::new(twist_curve, first_nonsquare::<F19>())
        .expect("quadratic twist should build");
    let twist_relation = twist_package
        .frobenius_relation()
        .expect("twist relation should compute");

    let isogeny =
        ScalarMultiplicationIsogeny::new(f41_curve(), 2).expect("scalar isogeny should build");
    let isogeny_relation = isogeny
        .frobenius_relation_report()
        .expect("relation should compute");

    let twist_description = describe_quadratic_twist_frobenius_relation(&twist_relation);
    let isogeny_description = describe_isogeny_frobenius_relation(&isogeny_relation);

    assert!(twist_description.contains("expected sum 2q + 2"));
    assert!(twist_description.contains("trace negation t' = -t holds: yes"));
    assert!(isogeny_description.contains("isogeny degree: 4"));
    assert!(isogeny_description.contains("same curve order: yes"));
    assert!(isogeny_description.contains("same trace: yes"));
}

#[test]
fn graph_visualization_reports_reference_and_per_node_verdicts() {
    let graph = IsogenyGraphBuilder::new(f41_curve(), 2)
        .max_depth(1)
        .build()
        .expect("small graph should build");
    let report = graph
        .frobenius_relation_report()
        .expect("graph report should compute");

    let description = describe_isogeny_graph_frobenius_report(&report);

    assert!(description.contains("reference node: 0"));
    assert!(description.contains("checked nodes:"));
    assert!(description.contains("checked edges:"));
    assert!(description.contains("per-node verdicts:"));
    assert!(description.contains("node 0: yes"));
}

#[test]
fn visualizable_trait_is_hooked_up_for_frobenius_objects() {
    let relative = RelativeFrobenius::for_field::<F17Squared>(1);
    let curve = ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid curve");
    let point = curve
        .point(F43::zero(), F43::one())
        .expect("sample point should lie on the curve");
    let orbit = relative_frobenius_orbit(&curve, &point).expect("relative orbit should compute");

    assert!(relative.format_compact().contains("π_(17^2)"));
    assert!(relative.describe().contains("field order q: 289"));
    assert!(orbit.format_compact().contains("period 1 orbit"));
    assert!(orbit.describe().contains("Frobenius orbit"));
}
