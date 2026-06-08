use proptest::prelude::*;

use num_bigint::{BigInt, BigUint};
use std::collections::HashSet;

use crate::elliptic_curves::frobenius::{
    AbsoluteFrobenius, FrobeniusCharacteristicPolynomial, FrobeniusCurveType,
    FrobeniusDiscriminant, FrobeniusLocalZetaFunction, FrobeniusTrace, RelativeFrobenius,
    absolute_frobenius_on_exact_torsion, absolute_frobenius_orbit,
    absolute_frobenius_orbits_on_points, absolute_frobenius_power_point,
    compare_extension_count_with_enumeration, frobenius_twist_power,
    relative_frobenius_on_exact_torsion, relative_frobenius_orbit,
    relative_frobenius_orbits_on_points, relative_frobenius_point,
    verify_frobenius_characteristic_equation_at_point,
    verify_frobenius_characteristic_equation_exhaustive, verify_hasse_bound,
    verify_isogeny_frobenius_relation, verify_isogeny_graph_frobenius_relation,
};
use crate::elliptic_curves::{
    AffineCurveModel, AffinePoint, CurveModel, EnumerableCurveModel, FrobeniusTraceCurveModel,
    ShortWeierstrassCurve, ShortWeierstrassQuadraticTwist, TwistKind,
};
use crate::fields::{EnumerableFiniteField, Field, FiniteFieldDescriptor, Fp, SqrtField};
use crate::isogenies::graphs::IsogenyGraphBuilder;
use crate::isogenies::{Isogeny, ScalarMultiplicationIsogeny, VeluIsogeny};
use crate::proptest_support::{
    ProptestF17Sqrt3Field, curve_and_rational_point, non_singular_short_weierstrass_curve,
};

type F17 = Fp<17>;
type F43 = Fp<43>;
type F19 = Fp<19>;
type F41 = Fp<41>;
type F17Squared = ProptestF17Sqrt3Field;

crate::fields::define_fp_quadratic_extension!(
    spec: F43Sqrt2Spec,
    field: F43Sqrt2,
    base: F43,
    non_residue: 2,
    name: "F43(sqrt(2))",
);

fn alpha() -> <F43Sqrt2 as Field>::Elem {
    F43Sqrt2::element(vec![F43::zero(), F43::one()])
}

fn f41_curve() -> ShortWeierstrassCurve<F41> {
    ShortWeierstrassCurve::<F41>::new(F41::from_i64(2), F41::from_i64(3)).expect("valid F41 curve")
}

fn depth_one_f41_graph() -> crate::isogenies::graphs::IsogenyGraph<ShortWeierstrassCurve<F41>> {
    IsogenyGraphBuilder::new(f41_curve(), 2)
        .max_depth(1)
        .build()
        .expect("depth-one F41 graph should build")
}

fn f43_quadratic_twist_factor() -> <F43 as Field>::Elem {
    first_nonsquare::<F43>()
}

fn lift_f17_curve_to_f17_squared(
    curve: &ShortWeierstrassCurve<F17>,
) -> ShortWeierstrassCurve<F17Squared> {
    ShortWeierstrassCurve::<F17Squared>::new(
        F17Squared::from_base(*curve.a()),
        F17Squared::from_base(*curve.b()),
    )
    .expect("lifting a smooth F17 curve to F17^2 should preserve smoothness")
}

fn nz(n: u32) -> core::num::NonZeroU32 {
    core::num::NonZeroU32::new(n).expect("test degrees are positive")
}

fn find_non_fixed_point(curve: &ShortWeierstrassCurve<F43Sqrt2>) -> AffinePoint<F43Sqrt2> {
    for x in F43Sqrt2::elements() {
        for y in F43Sqrt2::elements() {
            if let Ok(point) = curve.point(x.clone(), y) {
                let image = absolute_frobenius_power_point(curve, &point, 1)
                    .expect("absolute Frobenius should evaluate on on-curve inputs");
                if image != point {
                    return point;
                }
            }
        }
    }

    panic!("expected a point in E(F43^2) that is not fixed by Frob_43");
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

#[test]
fn absolute_frobenius_metadata_uses_the_prime_characteristic() {
    let frobenius = AbsoluteFrobenius::for_field::<F43>(3);

    assert_eq!(frobenius.characteristic, 43);
    assert_eq!(frobenius.power, 3);
    assert!(!frobenius.is_identity());
}

#[test]
fn relative_frobenius_metadata_uses_the_full_base_field_descriptor() {
    let frobenius = RelativeFrobenius::for_field::<F43Sqrt2>(2);

    assert_eq!(frobenius.base_field.characteristic, 43);
    assert_eq!(frobenius.base_field.extension_degree.get(), 2);
    assert_eq!(frobenius.power, 2);
}

#[test]
fn frobenius_twist_can_change_extension_field_coefficients() {
    let curve = ShortWeierstrassCurve::<F43Sqrt2>::new(alpha(), F43Sqrt2::zero())
        .expect("non-zero linear coefficient should define a smooth curve");

    let first_twist = frobenius_twist_power(&curve, 1).expect("Frobenius twist should be valid");
    let second_twist =
        frobenius_twist_power(&curve, 2).expect("second Frobenius twist should be valid");

    assert_ne!(first_twist, curve);
    assert_eq!(second_twist, curve);
}

#[test]
fn absolute_frobenius_fixes_the_point_at_infinity() {
    let curve = ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid curve");
    let point = AffinePoint::<F43>::infinity();

    let image = absolute_frobenius_power_point(&curve, &point, 4)
        .expect("Frobenius should fix the distinguished identity");

    assert_eq!(image, point);
}

#[test]
fn absolute_frobenius_fixes_prime_field_rational_points() {
    let curve = ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid curve");
    let point = curve
        .point(F43::zero(), F43::one())
        .expect("(0, 1) should lie on y^2 = x^3 + x + 1 over F43");

    let image = absolute_frobenius_power_point(&curve, &point, 1)
        .expect("Frobenius should evaluate on rational points");

    assert_eq!(image, point);
}

#[test]
fn relative_frobenius_orbit_on_a_rational_point_is_a_singleton() {
    let curve = ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid curve");
    let point = curve
        .point(F43::zero(), F43::one())
        .expect("(0, 1) should lie on the curve");

    let orbit = relative_frobenius_orbit(&curve, &point).expect("relative orbit should evaluate");

    assert_eq!(orbit.start(), &point);
    assert_eq!(orbit.points(), &[point]);
    assert_eq!(orbit.period(), 1);
}

#[test]
fn relative_frobenius_orbits_on_points_are_singletons() {
    let curve = ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid curve");

    let orbits =
        relative_frobenius_orbits_on_points(&curve).expect("relative point orbits should evaluate");

    assert_eq!(orbits.len(), curve.order());
    assert!(orbits.iter().all(|orbit| orbit.period() == 1));
    assert_eq!(
        orbits.iter().map(|orbit| orbit.period()).sum::<usize>(),
        curve.order()
    );
}

#[test]
fn relative_frobenius_can_fix_a_point_that_absolute_frobenius_moves() {
    let curve = ShortWeierstrassCurve::<F43Sqrt2>::new(
        F43Sqrt2::from_base(F43::one()),
        F43Sqrt2::from_base(F43::one()),
    )
    .expect("base-defined curve should stay smooth over F43^2");
    let point = find_non_fixed_point(&curve);

    let absolute_image = absolute_frobenius_power_point(&curve, &point, 1)
        .expect("absolute Frobenius should evaluate");
    let relative_image =
        relative_frobenius_point(&curve, &point).expect("relative Frobenius should evaluate");

    assert_ne!(absolute_image, point);
    assert_eq!(relative_image, point);
}

#[test]
fn absolute_frobenius_orbit_can_have_period_two_over_f43_squared() {
    let curve = ShortWeierstrassCurve::<F43Sqrt2>::new(
        F43Sqrt2::from_base(F43::one()),
        F43Sqrt2::from_base(F43::one()),
    )
    .expect("base-defined curve should stay smooth over F43^2");
    let point = find_non_fixed_point(&curve);

    let orbit =
        absolute_frobenius_orbit(&curve, &point, 1).expect("absolute orbit should evaluate");
    let image =
        absolute_frobenius_power_point(&curve, &point, 1).expect("absolute Frobenius should apply");

    assert_eq!(orbit.start(), &point);
    assert_eq!(orbit.points(), &[point.clone(), image]);
    assert_eq!(orbit.period(), 2);
}

#[test]
fn absolute_frobenius_orbits_on_points_partition_a_quadratic_extension_curve() {
    let curve = ShortWeierstrassCurve::<F43Sqrt2>::new(
        F43Sqrt2::from_base(F43::one()),
        F43Sqrt2::from_base(F43::one()),
    )
    .expect("base-defined curve should stay smooth over F43^2");

    let orbits = absolute_frobenius_orbits_on_points(&curve, 1)
        .expect("absolute Frobenius point orbits should evaluate");

    assert_eq!(
        orbits.iter().map(|orbit| orbit.period()).sum::<usize>(),
        curve.order()
    );
    assert!(orbits.iter().any(|orbit| orbit.period() == 2));
}

#[test]
fn absolute_frobenius_orbit_rejects_curves_not_preserved_by_the_chosen_power() {
    let curve = ShortWeierstrassCurve::<F43Sqrt2>::new(alpha(), F43Sqrt2::zero())
        .expect("valid extension-defined curve");
    let point = AffinePoint::<F43Sqrt2>::infinity();

    assert_eq!(
        absolute_frobenius_orbit(&curve, &point, 1),
        Err(crate::elliptic_curves::CurveError::AbsoluteFrobeniusDoesNotPreserveCurve { power: 1 })
    );
}

#[test]
fn curves_over_enumerable_quadratic_extensions_can_list_points_and_orders() {
    let curve = ShortWeierstrassCurve::<F43Sqrt2>::new(
        F43Sqrt2::from_base(F43::one()),
        F43Sqrt2::from_base(F43::one()),
    )
    .expect("base-defined curve should stay smooth over F43^2");

    let points = curve.points();

    assert!(!points.is_empty());
    assert_eq!(curve.order(), points.len());
}

#[test]
fn frobenius_trace_matches_the_counting_formula() {
    let curve = ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid curve");

    let report = curve
        .frobenius_trace()
        .expect("Frobenius trace should compute");

    assert_eq!(report.base_field().characteristic, 43);
    assert_eq!(report.base_field().extension_degree.get(), 1);
    assert_eq!(report.field_order(), 43);
    assert_eq!(report.curve_order(), curve.order() as u64);
    assert_eq!(report.trace(), 43_i64 + 1 - curve.order() as i64);
}

#[test]
fn frobenius_trace_from_order_and_order_from_trace_roundtrip() {
    let base_field =
        FiniteFieldDescriptor::new(43, core::num::NonZeroU32::new(1).expect("1 is non-zero"))
            .expect("F43 descriptor should be valid");
    let report = FrobeniusTrace::from_order(base_field.clone(), 48)
        .expect("small Frobenius trace package should build");

    assert_eq!(report.trace(), -4);
    assert_eq!(
        FrobeniusTrace::curve_order_from_trace(base_field, report.trace()),
        Ok(report.curve_order())
    );
}

#[test]
fn frobenius_curve_order_from_trace_rejects_non_positive_orders() {
    let base_field =
        FiniteFieldDescriptor::new(43, core::num::NonZeroU32::new(1).expect("1 is non-zero"))
            .expect("F43 descriptor should be valid");
    assert_eq!(
        FrobeniusTrace::curve_order_from_trace(base_field, 100),
        Err(crate::elliptic_curves::CurveError::InvalidFrobeniusTrace { trace: 100 })
    );
}

#[test]
fn extension_count_over_degree_one_recovers_the_base_curve_order() {
    let curve = ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid curve");
    let trace = curve
        .frobenius_trace()
        .expect("Frobenius trace should compute");

    let report = trace.curve_order_over_extension(nz(1));

    assert_eq!(report.frobenius_trace(), &trace);
    assert_eq!(report.extension_degree(), nz(1));
    assert_eq!(report.extension_field_order(), &BigUint::from(43u32));
    assert_eq!(report.power_sum(), &BigInt::from(trace.trace()));
    assert_eq!(report.curve_order(), &BigUint::from(curve.order() as u64));
}

#[test]
fn extension_count_over_degree_two_matches_explicit_quadratic_extension_enumeration() {
    let base_curve =
        ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid curve");
    let extension_curve = ShortWeierstrassCurve::<F43Sqrt2>::new(
        F43Sqrt2::from_base(F43::one()),
        F43Sqrt2::from_base(F43::one()),
    )
    .expect("base-defined curve should stay smooth over F43^2");
    let trace = base_curve
        .frobenius_trace()
        .expect("Frobenius trace should compute");

    let report = trace.curve_order_over_extension(nz(2));
    let expected_power_sum =
        BigInt::from(trace.trace()) * BigInt::from(trace.trace()) - BigInt::from(2 * 43);

    assert_eq!(report.extension_degree(), nz(2));
    assert_eq!(report.extension_field_order(), &BigUint::from(43u32).pow(2));
    assert_eq!(report.power_sum(), &expected_power_sum);
    assert_eq!(
        report.curve_order(),
        &BigUint::from(extension_curve.order() as u64)
    );
}

#[test]
fn extension_count_sequence_matches_individual_reports() {
    let curve = ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid curve");
    let trace = curve
        .frobenius_trace()
        .expect("Frobenius trace should compute");

    let sequence = trace.curve_orders_over_extensions_through(nz(4));

    assert_eq!(sequence.frobenius_trace(), &trace);
    assert_eq!(sequence.reports().len(), 4);
    assert_eq!(
        sequence
            .reports()
            .iter()
            .map(|report| report.extension_degree().get())
            .collect::<Vec<_>>(),
        vec![1, 2, 3, 4]
    );

    for report in sequence.reports() {
        let individual = trace.curve_order_over_extension(report.extension_degree());
        assert_eq!(report, &individual);
    }
}

#[test]
fn extension_count_can_represent_degrees_far_beyond_u128_field_powers() {
    let curve = ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid curve");
    let trace = curve
        .frobenius_trace()
        .expect("Frobenius trace should compute");

    let report = trace.curve_order_over_extension(nz(25));

    assert!(report.extension_field_order() > &BigUint::from(u128::MAX));
    assert!(report.curve_order() > &BigUint::from(u128::MAX));
}

#[test]
fn extension_count_comparison_report_distinguishes_frobenius_and_enumeration_routes() {
    let base_curve =
        ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid base curve");
    let extension_curve = ShortWeierstrassCurve::<F43Sqrt2>::new(
        F43Sqrt2::from_base(F43::one()),
        F43Sqrt2::from_base(F43::one()),
    )
    .expect("base-defined curve should stay smooth over F43^2");
    let trace = base_curve
        .frobenius_trace()
        .expect("Frobenius trace should compute");

    let comparison = compare_extension_count_with_enumeration(&extension_curve, &trace)
        .expect("comparison between Frobenius and exhaustive paths should compute");

    assert_eq!(comparison.trace_base_field(), trace.base_field());
    assert_eq!(comparison.curve_base_field().characteristic, 43);
    assert_eq!(comparison.curve_base_field().extension_degree.get(), 2);
    assert_eq!(comparison.relative_extension_degree(), nz(2));
    assert_eq!(
        comparison.frobenius_count().curve_order(),
        comparison.exhaustive_curve_order()
    );
    assert_eq!(
        comparison.exhaustive_curve_order(),
        &BigUint::from(extension_curve.order() as u64)
    );
    assert!(comparison.agrees());
}

#[test]
fn characteristic_polynomial_is_derived_from_the_frobenius_trace() {
    let curve = ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid curve");
    let trace = curve
        .frobenius_trace()
        .expect("Frobenius trace should compute");

    let polynomial = trace.characteristic_polynomial();

    assert_eq!(polynomial.base_field(), trace.base_field());
    assert_eq!(polynomial.trace(), trace.trace());
    assert_eq!(polynomial.field_order(), trace.field_order());
}

#[test]
fn characteristic_polynomial_discriminant_and_evaluation_match_the_formula() {
    let polynomial = FrobeniusCharacteristicPolynomial::new(
        FiniteFieldDescriptor::new(43, core::num::NonZeroU32::new(1).expect("1 is non-zero"))
            .expect("F43 descriptor should be valid"),
        -4,
    );

    assert_eq!(polynomial.discriminant(), 16 - 4 * 43);
    assert_eq!(polynomial.evaluate_at_integer(0), 43);
    assert_eq!(polynomial.evaluate_at_integer(1), 48);
    assert_eq!(polynomial.evaluate_at_integer(2), 55);
}

#[test]
fn frobenius_discriminant_is_derived_from_the_trace_package() {
    let curve = ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid curve");
    let trace = curve
        .frobenius_trace()
        .expect("Frobenius trace should compute");

    let discriminant = trace.discriminant();

    assert_eq!(discriminant.frobenius_trace(), &trace);
    assert_eq!(discriminant.base_field(), trace.base_field());
    assert_eq!(discriminant.curve_order(), trace.curve_order());
    assert_eq!(discriminant.trace(), trace.trace());
}

#[test]
fn frobenius_discriminant_matches_the_quadratic_formula() {
    let base_field =
        FiniteFieldDescriptor::new(43, core::num::NonZeroU32::new(1).expect("1 is non-zero"))
            .expect("F43 descriptor should be valid");
    let trace = FrobeniusTrace::from_order(base_field, 41).expect("t = 3 should be valid over F43");

    let discriminant = trace.discriminant();

    assert_eq!(
        discriminant.quadratic_discriminant().value(),
        &num_bigint::BigInt::from(-163)
    );
    assert!(discriminant.is_negative());
    assert!(!discriminant.is_zero());
    assert!(!discriminant.is_positive());
    assert!(discriminant.is_fundamental());
}

#[test]
fn frobenius_discriminant_constructor_matches_trace_method() {
    let base_field =
        FiniteFieldDescriptor::new(43, core::num::NonZeroU32::new(1).expect("1 is non-zero"))
            .expect("F43 descriptor should be valid");
    let trace =
        FrobeniusTrace::from_order(base_field, 48).expect("t = -4 should be valid over F43");

    let from_method = trace.discriminant();
    let from_constructor = FrobeniusDiscriminant::new(trace.clone());

    assert_eq!(from_constructor, from_method);
    assert_eq!(
        from_constructor.quadratic_discriminant().value(),
        &num_bigint::BigInt::from(-156)
    );
}

#[test]
fn frobenius_discriminant_quadratic_factorization_matches_the_expected_split() {
    let base_field =
        FiniteFieldDescriptor::new(43, core::num::NonZeroU32::new(1).expect("1 is non-zero"))
            .expect("F43 descriptor should be valid");
    let trace =
        FrobeniusTrace::from_order(base_field, 48).expect("t = -4 should be valid over F43");
    let discriminant = trace.discriminant();

    let factorization = discriminant.quadratic_factorization().expect(
        "negative Frobenius discriminant should admit an imaginary quadratic factorization",
    );

    assert_eq!(
        factorization.discriminant(),
        discriminant.quadratic_discriminant()
    );
    assert_eq!(factorization.conductor(), &num_bigint::BigUint::from(2u8));
    assert_eq!(
        factorization.fundamental_discriminant(),
        &crate::elliptic_curves::QuadraticDiscriminant::new(-39)
    );
}

#[test]
fn frobenius_discriminant_can_build_the_imaginary_quadratic_order() {
    let base_field =
        FiniteFieldDescriptor::new(43, core::num::NonZeroU32::new(1).expect("1 is non-zero"))
            .expect("F43 descriptor should be valid");
    let trace =
        FrobeniusTrace::from_order(base_field, 48).expect("t = -4 should be valid over F43");
    let discriminant = trace.discriminant();

    let order = discriminant
        .frobenius_order()
        .expect("negative Frobenius discriminant should define an imaginary quadratic order");

    assert_eq!(
        order.fundamental_discriminant(),
        &crate::elliptic_curves::QuadraticDiscriminant::new(-39)
    );
    assert_eq!(order.conductor(), &num_bigint::BigUint::from(2u8));
    assert_eq!(order.discriminant(), discriminant.quadratic_discriminant());
    assert!(!order.is_maximal());
}

#[test]
fn frobenius_order_matches_the_existing_imaginary_quadratic_order_helper() {
    let base_field =
        FiniteFieldDescriptor::new(43, core::num::NonZeroU32::new(1).expect("1 is non-zero"))
            .expect("F43 descriptor should be valid");
    let trace =
        FrobeniusTrace::from_order(base_field, 48).expect("t = -4 should be valid over F43");
    let discriminant = trace.discriminant();

    let from_existing = discriminant
        .frobenius_order()
        .expect("existing helper should build the order");
    let from_frobenius_order = discriminant
        .frobenius_order()
        .expect("frobenius_order should build the same order");

    assert_eq!(from_frobenius_order, from_existing);
}

#[test]
fn frobenius_order_is_contained_in_the_maximal_order() {
    let base_field =
        FiniteFieldDescriptor::new(43, core::num::NonZeroU32::new(1).expect("1 is non-zero"))
            .expect("F43 descriptor should be valid");
    let trace =
        FrobeniusTrace::from_order(base_field, 48).expect("t = -4 should be valid over F43");
    let discriminant = trace.discriminant();
    let frobenius_order = discriminant
        .frobenius_order()
        .expect("Frobenius order should exist in the imaginary case");
    let maximal_order = discriminant
        .maximal_order()
        .expect("the same factorization should recover O_K");

    assert!(frobenius_order.is_suborder_of(&maximal_order));
    assert!(maximal_order.is_overorder_of(&frobenius_order));
}

#[test]
fn endomorphism_ring_report_from_frobenius_is_honest_about_candidates() {
    let base_field = FiniteFieldDescriptor::new(13, nz(1))
        .expect("F_13 metadata should be internally consistent");
    let trace =
        FrobeniusTrace::from_order(base_field, 18).expect("t = -4 should be valid over F13");
    let report = trace
        .discriminant()
        .endomorphism_ring_report()
        .expect("ordinary Frobenius data should produce a candidate report");

    assert_eq!(report.frobenius_discriminant().trace(), -4);
    assert!(matches!(
        report,
        crate::elliptic_curves::EndomorphismRingReport::OrdinaryQuadraticOrderCandidates { .. }
    ));
    assert_eq!(report.candidate_count(), Some(2));
    assert_eq!(report.sandwich_inclusion_holds(), Some(true));
    assert_eq!(
        report
            .candidate_orders()
            .expect("ordinary branch should expose candidate orders")
            .iter()
            .map(|order| order.conductor().clone())
            .collect::<Vec<_>>(),
        vec![BigUint::from(1u8), BigUint::from(3u8)]
    );
}

#[test]
fn supersingular_frobenius_data_uses_the_supersingular_report_branch() {
    let base_field =
        FiniteFieldDescriptor::new(5, nz(1)).expect("F_5 metadata should be internally consistent");
    let trace = FrobeniusTrace::from_order(base_field, 6).expect("t = 0 should be valid over F5");
    let report = trace
        .discriminant()
        .endomorphism_ring_report()
        .expect("supersingular Frobenius data should still produce a report");

    assert!(matches!(
        report,
        crate::elliptic_curves::EndomorphismRingReport::SupersingularQuaternionicPlaceholder { .. }
    ));
    assert_eq!(report.frobenius_discriminant().trace(), 0);
    assert_eq!(report.candidate_count(), None);
}

#[test]
fn characteristic_polynomial_pretty_and_display_use_the_same_educational_surface() {
    let positive_trace = FrobeniusCharacteristicPolynomial::new(
        FiniteFieldDescriptor::new(43, core::num::NonZeroU32::new(1).expect("1 is non-zero"))
            .expect("F43 descriptor should be valid"),
        3,
    );
    let zero_trace = FrobeniusCharacteristicPolynomial::new(
        FiniteFieldDescriptor::new(43, core::num::NonZeroU32::new(1).expect("1 is non-zero"))
            .expect("F43 descriptor should be valid"),
        0,
    );
    let negative_trace = FrobeniusCharacteristicPolynomial::new(
        FiniteFieldDescriptor::new(43, core::num::NonZeroU32::new(1).expect("1 is non-zero"))
            .expect("F43 descriptor should be valid"),
        -4,
    );

    assert_eq!(positive_trace.pretty(), "T^2 - 3T + 43");
    assert_eq!(zero_trace.pretty(), "T^2 + 43");
    assert_eq!(negative_trace.pretty(), "T^2 + 4T + 43");
    assert_eq!(format!("{negative_trace}"), negative_trace.pretty());
}

#[test]
fn local_zeta_function_is_derived_from_the_characteristic_polynomial() {
    let polynomial = FrobeniusCharacteristicPolynomial::new(
        FiniteFieldDescriptor::new(43, core::num::NonZeroU32::new(1).expect("1 is non-zero"))
            .expect("F43 descriptor should be valid"),
        -4,
    );

    let zeta = polynomial.local_zeta_function();

    assert_eq!(zeta.characteristic_polynomial(), &polynomial);
    assert_eq!(zeta.base_field(), polynomial.base_field());
    assert_eq!(zeta.field_order(), polynomial.field_order());
    assert_eq!(zeta.trace(), polynomial.trace());
}

#[test]
fn local_zeta_function_is_derived_from_the_frobenius_trace() {
    let curve = ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid curve");
    let trace = curve
        .frobenius_trace()
        .expect("Frobenius trace should compute");

    let zeta = trace.local_zeta_function();

    assert_eq!(
        zeta,
        trace.characteristic_polynomial().local_zeta_function()
    );
}

#[test]
fn local_zeta_function_pretty_and_display_use_the_same_educational_surface() {
    let prime_field_zeta = FrobeniusLocalZetaFunction::from_characteristic_polynomial(
        FrobeniusCharacteristicPolynomial::new(
            FiniteFieldDescriptor::new(17, core::num::NonZeroU32::new(1).expect("1 is non-zero"))
                .expect("F17 descriptor should be valid"),
            2,
        ),
    );
    let extension_field_zeta = FrobeniusLocalZetaFunction::from_characteristic_polynomial(
        FrobeniusCharacteristicPolynomial::new(
            FiniteFieldDescriptor::new(43, core::num::NonZeroU32::new(2).expect("2 is non-zero"))
                .expect("F43^2 descriptor should be valid"),
            -4,
        ),
    );

    assert_eq!(prime_field_zeta.numerator_pretty(), "1 - 2T + 17T²");
    assert_eq!(prime_field_zeta.denominator_pretty(), "(1 - T)(1 - 17T)");
    assert_eq!(
        prime_field_zeta.pretty(),
        "Z(E/F_17, T) = (1 - 2T + 17T²) / ((1 - T)(1 - 17T))"
    );
    assert_eq!(format!("{prime_field_zeta}"), prime_field_zeta.pretty());

    assert_eq!(
        extension_field_zeta.pretty(),
        "Z(E/F_(43^2), T) = (1 + 4T + 1849T²) / ((1 - T)(1 - 1849T))"
    );
}

#[test]
fn scalar_isogeny_frobenius_relation_holds_on_a_small_f41_curve() {
    let isogeny = ScalarMultiplicationIsogeny::new(f41_curve(), 2)
        .expect("scalar multiplication isogeny should build");

    let relation = verify_isogeny_frobenius_relation(&isogeny)
        .expect("Frobenius relation should compute on both sides of the isogeny");

    assert_eq!(relation.degree(), isogeny.degree());
    assert!(relation.same_curve_order());
    assert!(relation.same_trace());
    assert!(relation.holds());
}

#[test]
fn velu_isogeny_frobenius_relation_holds_on_the_small_f41_example() {
    let domain = f41_curve();
    let generator = domain
        .point(F41::from_i64(40), F41::from_i64(0))
        .expect("sample generator should lie on the curve");
    let isogeny =
        VeluIsogeny::from_generator(domain, generator).expect("sample Vélu isogeny should build");

    let relation = verify_isogeny_frobenius_relation(&isogeny)
        .expect("Frobenius relation should compute on the Vélu example");

    assert_eq!(relation.degree(), isogeny.degree());
    assert_eq!(
        relation.domain().curve_order(),
        relation.codomain().curve_order()
    );
    assert_eq!(relation.domain().trace(), relation.codomain().trace());
    assert!(relation.same_curve_order());
    assert!(relation.same_trace());
    assert!(relation.holds());
}

#[test]
fn isogeny_graph_frobenius_report_holds_on_the_depth_one_f41_graph() {
    let graph = depth_one_f41_graph();

    let report = verify_isogeny_graph_frobenius_relation(&graph)
        .expect("graph Frobenius relation should compute on the F41 graph");

    assert_eq!(report.checked_nodes(), graph.node_count());
    assert_eq!(report.checked_edges(), graph.edge_count());
    assert_eq!(report.reference_node().0, 0);
    assert!(report.all_same_curve_order());
    assert!(report.all_same_trace());
    assert!(report.holds());
    assert!(report.node_reports().iter().all(|node| node.holds()));
}

#[test]
fn isogeny_graph_frobenius_report_handles_a_singleton_graph() {
    let graph = IsogenyGraphBuilder::new(f41_curve(), 2)
        .max_depth(0)
        .build()
        .expect("singleton F41 graph should build");

    let report = verify_isogeny_graph_frobenius_relation(&graph)
        .expect("singleton graph Frobenius relation should compute");

    assert_eq!(report.checked_nodes(), 1);
    assert_eq!(report.checked_edges(), 0);
    assert!(report.all_same_curve_order());
    assert!(report.all_same_trace());
    assert!(report.holds());
    assert_eq!(report.node_reports().len(), 1);
    assert!(report.node_reports()[0].holds());
}

#[test]
fn quadratic_twist_frobenius_relation_holds_for_a_nontrivial_f19_twist() {
    let curve = ShortWeierstrassCurve::<F19>::new(F19::from_i64(2), F19::from_i64(3))
        .expect("valid F19 curve");
    let package = ShortWeierstrassQuadraticTwist::new(curve, first_nonsquare::<F19>())
        .expect("non-zero twist factor should define a quadratic-twist package");

    let relation = package
        .frobenius_relation()
        .expect("Frobenius relation should compute by exhaustive counting");

    assert_eq!(relation.sum_orders(), 2 * 19 + 2);
    assert_eq!(relation.expected_sum(), 2 * 19 + 2);
    assert!(relation.holds());
    assert!(relation.trace_negation_holds());
}

#[test]
fn quadratic_twist_frobenius_relation_reports_trivial_twists_without_erroring() {
    let curve = ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid F43 curve");
    let package = ShortWeierstrassQuadraticTwist::new(curve, F43::from_i64(4))
        .expect("square twist factor should still define a twist package");

    let relation = package
        .frobenius_relation()
        .expect("trivial twists should still produce a report");

    assert_eq!(
        relation.original().curve_order(),
        relation.twist().curve_order()
    );
    assert_eq!(relation.original().trace(), relation.twist().trace());
    assert!(!relation.trace_negation_holds());
    assert!(!relation.holds());
}

#[test]
fn hasse_bound_report_matches_the_frobenius_trace_data() {
    let curve = ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid curve");
    let expected_trace = 43_i64 + 1 - curve.order() as i64;
    let expected_trace_square = i128::from(expected_trace) * i128::from(expected_trace);

    let report = verify_hasse_bound(&curve).expect("Hasse bound should verify");

    assert_eq!(report.frobenius_trace().field_order(), 43);
    assert_eq!(report.frobenius_trace().curve_order(), curve.order() as u64);
    assert_eq!(report.frobenius_trace().trace(), expected_trace);
    assert_eq!(report.trace_square(), expected_trace_square);
    assert_eq!(report.bound_square(), 4 * 43);
    assert_eq!(report.slack(), 4 * 43 - expected_trace_square);
    assert!(report.holds());
}

#[test]
fn curve_type_report_classifies_an_ordinary_f43_curve() {
    let curve = ShortWeierstrassCurve::<F43>::new(F43::zero(), F43::one()).expect("valid curve");
    let trace = curve
        .frobenius_trace()
        .expect("Frobenius trace should compute");
    let report = trace.curve_type_report();

    assert_eq!(trace.trace(), 8);
    assert_eq!(trace.curve_type(), FrobeniusCurveType::Ordinary);
    assert_eq!(report.frobenius_trace(), &trace);
    assert_eq!(report.curve_type(), FrobeniusCurveType::Ordinary);
    assert_eq!(report.trace_mod_characteristic(), 8);
    assert!(!report.characteristic_divides_trace());
    assert!(report.is_ordinary());
    assert!(!report.is_supersingular());
}

#[test]
fn curve_type_report_classifies_a_supersingular_f43_curve() {
    let curve = ShortWeierstrassCurve::<F43>::new(F43::one(), F43::zero()).expect("valid curve");
    let trace = curve
        .frobenius_trace()
        .expect("Frobenius trace should compute");
    let report = trace.curve_type_report();

    assert_eq!(trace.trace(), 0);
    assert_eq!(trace.curve_type(), FrobeniusCurveType::Supersingular);
    assert_eq!(report.curve_type(), FrobeniusCurveType::Supersingular);
    assert_eq!(report.trace_mod_characteristic(), 0);
    assert!(report.characteristic_divides_trace());
    assert!(!report.is_ordinary());
    assert!(report.is_supersingular());
}

#[test]
fn relative_frobenius_on_exact_torsion_reports_fixed_f43_two_torsion() {
    let curve = ShortWeierstrassCurve::<F43>::new(F43::zero(), F43::one()).expect("valid curve");

    let report = relative_frobenius_on_exact_torsion(&curve, 2)
        .expect("relative Frobenius on exact torsion should evaluate");

    assert_eq!(report.exact_order(), 2);
    assert_eq!(report.points().len(), 3);
    assert!(report.all_fixed());
    assert_eq!(report.fixed_count(), 3);
    assert_eq!(report.moved_count(), 0);
    assert_eq!(report.prime_field_rational_count(), 3);
    assert_eq!(report.extension_only_count(), 0);
    assert_eq!(report.orbit_count(), 3);
    assert_eq!(report.orbit_periods(), vec![1, 1, 1]);
    assert_eq!(report.fixed_by_absolute_frobenius_power_count(1), None);
    assert_eq!(
        report.count_with_minimal_absolute_frobenius_fixing_power(1),
        None
    );
    assert!(
        report
            .points()
            .iter()
            .all(|point| point.fixed_by_frobenius())
    );
    assert!(
        report
            .points()
            .iter()
            .all(|point| point.point() == point.frobenius_image())
    );
}

#[test]
fn relative_frobenius_on_exact_torsion_rejects_zero_order() {
    let curve = ShortWeierstrassCurve::<F43>::new(F43::zero(), F43::one()).expect("valid curve");

    assert_eq!(
        relative_frobenius_on_exact_torsion(&curve, 0),
        Err(crate::elliptic_curves::CurveError::InvalidTorsionOrder { order: 0 })
    );
}

#[test]
fn absolute_frobenius_on_exact_four_torsion_can_move_extension_field_points() {
    let curve = ShortWeierstrassCurve::<F43Sqrt2>::new(
        F43Sqrt2::from_base(F43::zero()),
        F43Sqrt2::from_base(F43::one()),
    )
    .expect("base-defined curve should stay smooth over F43^2");

    let report = absolute_frobenius_on_exact_torsion(&curve, 4, 1)
        .expect("absolute Frobenius on exact torsion should evaluate");

    assert_eq!(report.exact_order(), 4);
    assert_eq!(report.points().len(), 12);
    assert!(!report.all_fixed());
    assert_eq!(report.fixed_count(), 0);
    assert_eq!(report.moved_count(), 12);
    assert_eq!(report.prime_field_rational_count(), 0);
    assert_eq!(report.extension_only_count(), 12);
    assert_eq!(report.orbit_count(), 6);
    assert_eq!(report.orbit_periods(), vec![2, 2, 2, 2, 2, 2]);
    assert_eq!(report.fixed_by_absolute_frobenius_power_count(1), Some(0));
    assert_eq!(report.fixed_by_absolute_frobenius_power_count(2), Some(12));
    assert_eq!(report.fixed_by_absolute_frobenius_power_count(3), Some(0));
    assert_eq!(report.fixed_by_absolute_frobenius_power_count(4), Some(12));
    assert_eq!(
        report.count_with_minimal_absolute_frobenius_fixing_power(1),
        Some(0)
    );
    assert_eq!(
        report.count_with_minimal_absolute_frobenius_fixing_power(2),
        Some(12)
    );
    assert!(
        report
            .points()
            .iter()
            .all(|point| !point.fixed_by_frobenius())
    );
    assert!(
        report
            .points()
            .iter()
            .all(|point| point.minimal_absolute_frobenius_fixing_power() == Some(2))
    );
    assert!(
        report
            .points()
            .iter()
            .all(|point| point.fixed_by_absolute_frobenius_power(1) == Some(false))
    );
    assert!(
        report
            .points()
            .iter()
            .all(|point| point.fixed_by_absolute_frobenius_power(2) == Some(true))
    );
}

#[test]
fn absolute_frobenius_squared_fixes_exact_four_torsion_over_f43_squared() {
    let curve = ShortWeierstrassCurve::<F43Sqrt2>::new(
        F43Sqrt2::from_base(F43::zero()),
        F43Sqrt2::from_base(F43::one()),
    )
    .expect("base-defined curve should stay smooth over F43^2");

    let report = absolute_frobenius_on_exact_torsion(&curve, 4, 2)
        .expect("absolute Frobenius square should evaluate");

    assert_eq!(report.exact_order(), 4);
    assert_eq!(report.points().len(), 12);
    assert!(report.all_fixed());
    assert_eq!(report.fixed_count(), 12);
    assert_eq!(report.moved_count(), 0);
    assert_eq!(report.prime_field_rational_count(), 12);
    assert_eq!(report.extension_only_count(), 0);
    assert_eq!(report.orbit_count(), 12);
    assert_eq!(report.orbit_periods(), vec![1; 12]);
    assert_eq!(report.fixed_by_absolute_frobenius_power_count(1), Some(0));
    assert_eq!(report.fixed_by_absolute_frobenius_power_count(2), Some(12));
    assert_eq!(
        report.count_with_minimal_absolute_frobenius_fixing_power(2),
        Some(12)
    );
    assert!(
        report
            .points()
            .iter()
            .all(|point| point.fixed_by_frobenius())
    );
}

#[test]
fn characteristic_equation_holds_on_a_prime_field_rational_point() {
    let curve = ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid curve");
    let point = curve
        .point(F43::zero(), F43::one())
        .expect("(0, 1) should lie on the curve");
    let characteristic_polynomial = curve
        .frobenius_trace()
        .expect("trace should compute")
        .characteristic_polynomial();

    let check = verify_frobenius_characteristic_equation_at_point(
        &curve,
        &point,
        &characteristic_polynomial,
    )
    .expect("characteristic equation check should evaluate");

    assert_eq!(check.point(), &point);
    assert_eq!(check.pi_q(), &point);
    assert_eq!(check.pi_q_squared(), &point);
    assert!(check.holds());
    assert!(curve.is_identity(check.lhs()));
}

#[test]
fn characteristic_equation_rejects_a_polynomial_from_the_wrong_base_field() {
    let curve = ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid curve");
    let point = curve
        .point(F43::zero(), F43::one())
        .expect("(0, 1) should lie on the curve");
    let wrong_polynomial = FrobeniusCharacteristicPolynomial::new(
        FiniteFieldDescriptor::new(43, core::num::NonZeroU32::new(2).expect("2 is non-zero"))
            .expect("F43^2 descriptor should be valid"),
        0,
    );

    assert_eq!(
        verify_frobenius_characteristic_equation_at_point(&curve, &point, &wrong_polynomial),
        Err(
            crate::elliptic_curves::CurveError::IncompatibleFrobeniusBaseField {
                curve_characteristic: 43,
                curve_extension_degree: 1,
                polynomial_characteristic: 43,
                polynomial_extension_degree: 2,
            }
        )
    );
}

#[test]
fn characteristic_equation_rejects_off_curve_points() {
    let curve = ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid curve");
    let off_curve = AffinePoint::<F43>::new(F43::zero(), F43::zero());
    let characteristic_polynomial = curve
        .frobenius_trace()
        .expect("trace should compute")
        .characteristic_polynomial();

    assert_eq!(
        verify_frobenius_characteristic_equation_at_point(
            &curve,
            &off_curve,
            &characteristic_polynomial,
        ),
        Err(crate::elliptic_curves::CurveError::PointNotOnCurve)
    );
}

#[test]
fn exhaustive_characteristic_equation_report_holds_on_a_small_prime_field_curve() {
    let curve = ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid curve");

    let report = verify_frobenius_characteristic_equation_exhaustive(&curve)
        .expect("exhaustive characteristic-equation check should evaluate");

    assert_eq!(report.frobenius_trace().curve_order(), curve.order() as u64);
    assert_eq!(report.checked_points(), curve.order());
    assert!(report.failed_checks().is_empty());
    assert!(report.all_hold());
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(24))]

    #[test]
    fn property_characteristic_equation_holds_for_sampled_f43_rational_points(
        (curve, point) in curve_and_rational_point::<43>(),
    ) {
        let characteristic_polynomial = curve
            .frobenius_trace()
            .expect("Frobenius trace should compute on small F43 curves")
            .characteristic_polynomial();

        let check = verify_frobenius_characteristic_equation_at_point(
            &curve,
            &point,
            &characteristic_polynomial,
        )
        .expect("characteristic equation should evaluate on enumerated rational points");

        prop_assert!(check.holds());
        prop_assert!(curve.is_identity(check.lhs()));
    }

    #[test]
    fn property_hasse_bound_holds_for_sampled_f43_curves(
        curve in non_singular_short_weierstrass_curve::<43>(),
    ) {
        let report = verify_hasse_bound(&curve)
            .expect("Hasse bound should verify on small F43 curves");

        prop_assert!(report.holds());
        prop_assert!(report.trace_square() <= report.bound_square());
        prop_assert_eq!(
            report.frobenius_trace().curve_order(),
            curve.order() as u64,
        );
    }

    #[test]
    fn property_curve_type_matches_trace_divisibility_over_f43(
        curve in non_singular_short_weierstrass_curve::<43>(),
    ) {
        let trace = curve
            .frobenius_trace()
            .expect("Frobenius trace should compute on small F43 curves");
        let report = trace.curve_type_report();
        let characteristic_divides_trace = trace.trace().rem_euclid(43) == 0;

        prop_assert_eq!(report.characteristic_divides_trace(), characteristic_divides_trace);
        prop_assert_eq!(report.trace_mod_characteristic() == 0, characteristic_divides_trace);
        prop_assert_eq!(
            report.curve_type(),
            if characteristic_divides_trace {
                FrobeniusCurveType::Supersingular
            } else {
                FrobeniusCurveType::Ordinary
            }
        );
    }

    #[test]
    fn property_relative_frobenius_on_exact_two_torsion_is_tautologically_fixed_over_f43(
        curve in non_singular_short_weierstrass_curve::<43>(),
    ) {
        let report = relative_frobenius_on_exact_torsion(&curve, 2)
            .expect("relative Frobenius on exact two-torsion should evaluate");

        prop_assert!(report.all_fixed());
        prop_assert_eq!(report.fixed_count(), report.points().len());
        prop_assert_eq!(report.moved_count(), 0);
        for point in report.points() {
            prop_assert!(point.fixed_by_frobenius());
            prop_assert_eq!(point.point(), point.frobenius_image());
        }
    }

    #[test]
    fn property_exhaustive_report_matches_pointwise_checks(
        curve in non_singular_short_weierstrass_curve::<43>(),
    ) {
        let report = verify_frobenius_characteristic_equation_exhaustive(&curve)
            .expect("exhaustive report should compute on small F43 curves");
        let characteristic_polynomial = report
            .frobenius_trace()
            .characteristic_polynomial();

        prop_assert_eq!(report.checked_points(), curve.order());
        prop_assert_eq!(report.all_hold(), report.failed_checks().is_empty());

        for point in curve.points() {
            let pointwise = verify_frobenius_characteristic_equation_at_point(
                &curve,
                &point,
                &characteristic_polynomial,
            )
            .expect("pointwise check should evaluate on enumerated rational points");

            let failed_in_report = report
                .failed_checks()
                .iter()
                .any(|check| check.point() == &point);

            prop_assert_eq!(pointwise.holds(), !failed_in_report);
        }
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(4))]

    #[test]
    fn property_quadratic_twist_frobenius_relation_holds_over_f43(
        curve in non_singular_short_weierstrass_curve::<43>(),
    ) {
        let package = ShortWeierstrassQuadraticTwist::new(curve, f43_quadratic_twist_factor())
            .expect("a fixed F43 non-square should define a quadratic twist package");
        let relation = package
            .frobenius_relation()
            .expect("quadratic-twist Frobenius relation should compute on small F43 curves");

        prop_assert_eq!(package.kind(), TwistKind::Quadratic);
        prop_assert!(relation.holds());
        prop_assert!(relation.trace_negation_holds());
        prop_assert_eq!(
            u128::from(relation.original().curve_order()) + u128::from(relation.twist().curve_order()),
            relation.expected_sum(),
        );
        prop_assert_eq!(relation.sum_orders(), relation.expected_sum());
        prop_assert_eq!(relation.twist().trace(), -relation.original().trace());
    }

    #[test]
    fn property_extension_counts_of_degree_one_and_two_match_enumeration(
        curve in non_singular_short_weierstrass_curve::<17>(),
    ) {
        let trace = curve
            .frobenius_trace()
            .expect("Frobenius trace should compute on small F17 curves");
        let lifted_curve = lift_f17_curve_to_f17_squared(&curve);

        let degree_one = compare_extension_count_with_enumeration(&curve, &trace)
            .expect("degree-one comparison should compute");
        let degree_two = compare_extension_count_with_enumeration(&lifted_curve, &trace)
            .expect("degree-two comparison should compute");

        prop_assert_eq!(degree_one.relative_extension_degree(), nz(1));
        prop_assert!(degree_one.agrees());
        prop_assert_eq!(
            degree_one.frobenius_count().curve_order(),
            &BigUint::from(curve.order() as u64),
        );

        prop_assert_eq!(degree_two.relative_extension_degree(), nz(2));
        prop_assert!(degree_two.agrees());
        prop_assert_eq!(
            degree_two.frobenius_count().curve_order(),
            &BigUint::from(lifted_curve.order() as u64),
        );
    }

    #[test]
    fn property_characteristic_polynomial_and_zeta_are_consistent_with_trace(
        curve in non_singular_short_weierstrass_curve::<43>(),
    ) {
        let trace = curve
            .frobenius_trace()
            .expect("Frobenius trace should compute on small F43 curves");
        let polynomial = trace.characteristic_polynomial();
        let zeta_from_polynomial = polynomial.local_zeta_function();
        let zeta_from_trace = trace.local_zeta_function();

        prop_assert_eq!(polynomial.base_field(), trace.base_field());
        prop_assert_eq!(polynomial.trace(), trace.trace());
        prop_assert_eq!(polynomial.field_order(), trace.field_order());
        prop_assert_eq!(polynomial.evaluate_at_integer(1), i128::from(trace.curve_order()));
        prop_assert_eq!(&zeta_from_polynomial, &zeta_from_trace);
        prop_assert_eq!(zeta_from_polynomial.base_field(), trace.base_field());
        prop_assert_eq!(zeta_from_polynomial.trace(), trace.trace());
        prop_assert_eq!(zeta_from_polynomial.field_order(), trace.field_order());
    }

    #[test]
    fn property_scalar_isogenies_preserve_curve_order_and_trace(
        curve in non_singular_short_weierstrass_curve::<43>(),
        scalar in 1u64..=4,
    ) {
        let trace = curve
            .frobenius_trace()
            .expect("Frobenius trace should compute on small F43 curves");
        let isogeny = ScalarMultiplicationIsogeny::new(curve, scalar)
            .expect("small scalar-multiplication isogenies should build");
        let relation = verify_isogeny_frobenius_relation(&isogeny)
            .expect("Frobenius relation should compute for scalar isogenies");

        prop_assert_eq!(relation.degree(), (scalar * scalar) as usize);
        prop_assert_eq!(relation.domain().curve_order(), trace.curve_order());
        prop_assert_eq!(relation.codomain().curve_order(), trace.curve_order());
        prop_assert_eq!(relation.domain().trace(), trace.trace());
        prop_assert_eq!(relation.codomain().trace(), trace.trace());
        prop_assert!(relation.same_curve_order());
        prop_assert!(relation.same_trace());
        prop_assert!(relation.holds());
    }

    #[test]
    fn property_isogeny_graph_frobenius_report_holds_on_small_f41_graphs(
        max_depth in 0usize..=1,
        deduplicate_by_base_field_isomorphism in any::<bool>(),
    ) {
        let graph = IsogenyGraphBuilder::new(f41_curve(), 2)
            .max_depth(max_depth)
            .deduplicate_by_base_field_isomorphism(deduplicate_by_base_field_isomorphism)
            .build()
            .expect("small F41 isogeny graphs should build");
        let report = verify_isogeny_graph_frobenius_relation(&graph)
            .expect("graph Frobenius relation should compute on small F41 graphs");

        prop_assert_eq!(report.checked_nodes(), graph.node_count());
        prop_assert_eq!(report.checked_edges(), graph.edge_count());
        prop_assert!(report.all_same_curve_order());
        prop_assert!(report.all_same_trace());
        prop_assert!(report.holds());
        prop_assert!(report.node_reports().iter().all(|node| node.holds()));
    }

    #[test]
    fn property_absolute_frobenius_orbits_partition_extension_points(
        curve in non_singular_short_weierstrass_curve::<17>(),
    ) {
        let lifted_curve = lift_f17_curve_to_f17_squared(&curve);
        let orbits = absolute_frobenius_orbits_on_points(&lifted_curve, 1)
            .expect("absolute Frobenius orbits should compute on lifted F17 curves");
        let mut seen_points = HashSet::new();
        let fixed_points_from_orbits = orbits
            .iter()
            .filter(|orbit| orbit.period() == 1)
            .count();

        for orbit in &orbits {
            prop_assert!(matches!(orbit.period(), 1 | 2));
            prop_assert_eq!(orbit.points().len(), orbit.period());
            prop_assert_eq!(
                absolute_frobenius_power_point(&lifted_curve, orbit.start(), orbit.period() as u32)
                    .expect("absolute Frobenius should evaluate on orbit starts"),
                orbit.start().clone(),
            );

            for point in orbit.points() {
                prop_assert!(seen_points.insert(point.clone()));
            }
        }

        prop_assert_eq!(
            orbits.iter().map(|orbit| orbit.period()).sum::<usize>(),
            lifted_curve.order(),
        );
        prop_assert_eq!(seen_points.len(), lifted_curve.order());
        prop_assert_eq!(fixed_points_from_orbits, curve.order());
    }
}
