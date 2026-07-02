use crate::fields::traits::*;
use num_bigint::{BigInt, BigUint};
use proptest::prelude::*;

use crate::elliptic_curves::{
    ShortWeierstrassCurve,
    frobenius::character_sum::CharacterSumPointCount,
    frobenius::extension_counts::compare_extension_count_with_enumeration,
    frobenius::{
        AbsoluteFrobenius, RelativeFrobenius,
        group_order::{GroupOrderReport, SmallFieldGroupOrderStrategy},
    },
    short_weierstrass::isomorphisms::{ShortWeierstrassQuadraticTwist, TwistKind},
    traits::{EnumerableCurveModel, FrobeniusTraceCurveModel},
};
use crate::fields::{
    finite_field_descriptor::FiniteFieldDescriptor,
    traits::{EnumerableFiniteField, SqrtField},
};
use crate::proptest_support::{
    config::CurveStrategyConfig, elliptic_curves::arb_nonsingular_curve,
};

type F17 = crate::fields::Fp17;
type F19 = crate::fields::Fp19;
type F43 = crate::fields::Fp43;

crate::fields::extension_field::define_fp_quadratic_extension!(
    spec: F43Sqrt2Spec,
    field: F43Sqrt2,
    base: F43,
    non_residue: 2,
    name: "F43(sqrt(2))",
);

crate::fields::extension_field::define_fp_quadratic_extension!(
    spec: F17Sqrt3Spec,
    field: F17Sqrt3,
    base: F17,
    non_residue: 3,
    name: "F17(sqrt(3))",
);

fn nz(n: u32) -> core::num::NonZeroU32 {
    core::num::NonZeroU32::new(n).expect("test degrees are positive")
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

fn lift_f43_curve_to_f43_sqrt2(
    curve: &ShortWeierstrassCurve<F43>,
) -> ShortWeierstrassCurve<F43Sqrt2> {
    ShortWeierstrassCurve::<F43Sqrt2>::new(
        F43Sqrt2::from_base(*curve.a()),
        F43Sqrt2::from_base(*curve.b()),
    )
    .expect("lifting an F43 curve to F43^2 should preserve smoothness")
}

fn lift_f17_curve_to_f17_sqrt3(
    curve: &ShortWeierstrassCurve<F17>,
) -> ShortWeierstrassCurve<F17Sqrt3> {
    ShortWeierstrassCurve::<F17Sqrt3>::new(
        F17Sqrt3::from_base(*curve.a()),
        F17Sqrt3::from_base(*curve.b()),
    )
    .expect("lifting an F17 curve to F17^2 should preserve smoothness")
}

#[test]
fn character_sum_count_matches_exhaustive_order_and_trace_over_f43() {
    let curve = ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid curve");

    let report = curve
        .group_order_by_small_field(SmallFieldGroupOrderStrategy::QuadraticCharacter)
        .expect("quadratic-character route should compute");

    let GroupOrderReport::QuadraticCharacter(report) = report else {
        panic!("quadratic-character strategy should preserve its report variant");
    };

    assert_eq!(report.curve_order(), BigUint::from(curve.order() as u64));
    assert_eq!(
        report.trace(),
        BigInt::from(43_i128 + 1 - curve.order() as i128)
    );
}

#[test]
fn character_sum_report_constructor_recovers_trace_by_negating_the_sum() {
    let base_field = FiniteFieldDescriptor::new(43, nz(1)).expect("F43 descriptor should be valid");
    let report = CharacterSumPointCount::new(base_field, 4).expect("constructor should succeed");

    assert_eq!(report.trace(), BigInt::from(-4));
    assert_eq!(report.curve_order(), BigUint::from(48u8));
}

#[test]
fn extension_count_over_degree_two_matches_explicit_quadratic_extension_enumeration() {
    let base_curve =
        ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid base curve");
    let extension_curve = lift_f43_curve_to_f43_sqrt2(&base_curve);
    let trace = base_curve
        .frobenius_trace()
        .expect("Frobenius trace should compute");

    let report = trace.curve_order_over_extension(nz(2));
    let expected_power_sum = trace.trace() * trace.trace() - BigInt::from(2 * 43);

    assert_eq!(report.extension_field_order(), &BigUint::from(43u32).pow(2));
    assert_eq!(report.power_sum(), &expected_power_sum);
    assert_eq!(
        report.curve_order(),
        &BigUint::from(extension_curve.order() as u64)
    );
}

#[test]
fn extension_count_comparison_report_distinguishes_frobenius_and_enumeration_routes() {
    let base_curve =
        ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid base curve");
    let extension_curve = lift_f43_curve_to_f43_sqrt2(&base_curve);
    let trace = base_curve
        .frobenius_trace()
        .expect("Frobenius trace should compute");

    let comparison = compare_extension_count_with_enumeration(&extension_curve, &trace)
        .expect("comparison between Frobenius and exhaustive paths should compute");

    assert_eq!(comparison.relative_extension_degree(), nz(2));
    assert!(comparison.agrees());
}

#[test]
fn frobenius_metadata_record_the_expected_parameters() {
    let absolute = AbsoluteFrobenius::for_field::<F43>(3);
    let relative = RelativeFrobenius::for_field::<F17Sqrt3>(2);

    assert_eq!(absolute.characteristic(), BigUint::from(43u8));
    assert_eq!(absolute.power(), 3);
    assert_eq!(relative.base_field().characteristic, BigUint::from(17u8));
    assert_eq!(relative.base_field().extension_degree.get(), 2);
    assert_eq!(relative.power(), 2);
}

#[test]
fn absolute_frobenius_orbits_partition_quadratic_extension_points() {
    let curve = ShortWeierstrassCurve::<F43Sqrt2>::new(
        F43Sqrt2::from_base(F43::one()),
        F43Sqrt2::from_base(F43::one()),
    )
    .expect("base-defined curve should stay smooth over F43^2");

    let point = curve
        .points()
        .into_iter()
        .find(|point| {
            curve
                .absolute_frobenius_power_point(point, 1)
                .expect("absolute Frobenius should evaluate")
                != *point
        })
        .expect("expected a non-fixed point over F43^2");

    let orbit = curve
        .absolute_frobenius_orbit(&point, 1)
        .expect("absolute orbit should evaluate");
    let partition = curve
        .absolute_frobenius_orbits_on_points(1)
        .expect("absolute partition should evaluate");

    assert_eq!(orbit.period(), 2);
    assert_eq!(
        partition
            .iter()
            .map(|candidate| candidate.period())
            .sum::<usize>(),
        curve.order()
    );
}

#[test]
fn quadratic_twist_frobenius_relation_holds_for_a_nontrivial_f19_twist() {
    let curve = ShortWeierstrassCurve::<F19>::new(F19::from_i64(2), F19::from_i64(3))
        .expect("valid F19 curve");
    let factor = first_nonsquare::<F19>();
    let package = ShortWeierstrassQuadraticTwist::new(curve, factor)
        .expect("fixed non-square should define a twist package");
    let relation = package
        .frobenius_relation()
        .expect("quadratic-twist Frobenius relation should compute");

    match package.kind() {
        TwistKind::Quadratic => {
            assert!(relation.holds());
            assert!(relation.trace_negation_holds());
        }
        TwistKind::Trivial => {
            assert!(relation.same_curve_order_holds());
            assert!(relation.same_trace_holds());
        }
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(8))]

    #[test]
    fn property_extension_counts_of_degree_one_and_two_match_enumeration(
        curve in arb_nonsingular_curve::<crate::fields::Fp17>(CurveStrategyConfig::default()),
    ) {
        let trace = curve
            .frobenius_trace()
            .expect("Frobenius trace should compute on small F17 curves");
        let lifted_curve = lift_f17_curve_to_f17_sqrt3(&curve);

        let degree_one = compare_extension_count_with_enumeration(&curve, &trace)
            .expect("degree-one comparison should compute");
        let degree_two = compare_extension_count_with_enumeration(&lifted_curve, &trace)
            .expect("degree-two comparison should compute");

        prop_assert_eq!(degree_one.relative_extension_degree(), nz(1));
        prop_assert!(degree_one.agrees());
        prop_assert_eq!(degree_two.relative_extension_degree(), nz(2));
        prop_assert!(degree_two.agrees());
    }

    #[test]
    fn property_quadratic_twist_frobenius_relation_holds_over_f43(
        curve in arb_nonsingular_curve::<crate::fields::Fp43>(CurveStrategyConfig::default()),
    ) {
        let package = ShortWeierstrassQuadraticTwist::new(curve, first_nonsquare::<F43>())
            .expect("a fixed F43 non-square should define a quadratic twist package");
        let relation = package
            .frobenius_relation()
            .expect("quadratic-twist Frobenius relation should compute on small F43 curves");

        prop_assert_eq!(relation.twist_kind(), package.kind());
        match package.kind() {
            TwistKind::Quadratic => {
                prop_assert!(relation.holds());
                prop_assert_eq!(relation.twist().trace(), -relation.original().trace());
            }
            TwistKind::Trivial => {
                prop_assert!(relation.same_curve_order_holds());
                prop_assert!(relation.same_trace_holds());
            }
        }
    }
}
