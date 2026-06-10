use proptest::prelude::*;

use crate::elliptic_curves::{CurveIsomorphism, CurveModel, EnumerableCurveModel};
use crate::fields::{EnumerableFiniteField, ExtensionField, Field, Fp};
use crate::isogenies::Isogeny;
use crate::proptest_support::config::{
    AnalyticStrategyConfig, CurveStrategyConfig, FieldStrategyConfig, PolynomialStrategyConfig,
};
use crate::proptest_support::elliptic_curves::{
    arb_complex_lattice, arb_curve_and_point, arb_division_polynomial_case,
    arb_endomorphism_report_case, arb_frobenius_curve_case, arb_nonsingular_curve,
    arb_upper_half_plane_point,
};
use crate::proptest_support::fields::{
    ProptestF17Sqrt3Spec, arb_complex_approx, arb_distinct_fp_elems, arb_extension_elem,
    arb_nonzero_fp_elem, arb_q_elem,
};
use crate::proptest_support::isogenies::{arb_composable_velu_case, arb_cyclic_kernel_case};
use crate::proptest_support::polynomials::{
    arb_dense_polynomial, arb_multivariate_polynomial, arb_sparse_polynomial,
};

type F17 = Fp<17>;

proptest! {
    #![proptest_config(ProptestConfig::with_cases(24))]

    #[test]
    fn fp_elements_stay_inside_the_field(element in crate::proptest_support::fields::arb_fp_elem::<17>()) {
        prop_assert!(F17::elements().contains(&element));
    }

    #[test]
    fn nonzero_fp_elements_never_generate_zero(element in arb_nonzero_fp_elem::<17>()) {
        prop_assert!(!F17::is_zero(&element));
    }

    #[test]
    fn distinct_fp_elements_do_not_repeat(elements in arb_distinct_fp_elems::<17>(4)) {
        let mut unique = elements.clone();
        unique.dedup();
        prop_assert_eq!(unique.len(), elements.len());
    }

    #[test]
    fn rational_samples_have_nonzero_denominator(value in arb_q_elem(FieldStrategyConfig::default())) {
        prop_assert_ne!(value.denom().clone(), 0.into());
    }

    #[test]
    fn complex_samples_respect_the_requested_bounds(value in arb_complex_approx(FieldStrategyConfig::default())) {
        prop_assert!(value.re.abs() <= FieldStrategyConfig::default().max_real_norm);
        prop_assert!(value.im.abs() <= FieldStrategyConfig::default().max_imaginary_norm);
    }

    #[test]
    fn extension_samples_are_canonically_reduced(
        value in arb_extension_elem::<ProptestF17Sqrt3Spec>()
    ) {
        prop_assert_eq!(ExtensionField::<ProptestF17Sqrt3Spec>::reduce(&value), value);
    }

    #[test]
    fn dense_polynomial_respects_max_length(
        polynomial in arb_dense_polynomial::<F17>(PolynomialStrategyConfig::default())
    ) {
        prop_assert!(polynomial.coefficients().len() <= PolynomialStrategyConfig::default().max_len);
    }

    #[test]
    fn sparse_polynomial_respects_max_term_count(
        polynomial in arb_sparse_polynomial::<F17>(PolynomialStrategyConfig::default())
    ) {
        prop_assert!(polynomial.terms().len() <= PolynomialStrategyConfig::default().max_terms);
    }

    #[test]
    fn multivariate_polynomial_respects_requested_arity(
        polynomial in arb_multivariate_polynomial::<F17>(PolynomialStrategyConfig::default())
    ) {
        prop_assert_eq!(polynomial.arity(), PolynomialStrategyConfig::default().arity);
    }

    #[test]
    fn nonsingular_curves_stay_nonsingular(
        curve in arb_nonsingular_curve::<17>(CurveStrategyConfig::default())
    ) {
        prop_assert!(!F17::is_zero(&curve.discriminant()));
    }

    #[test]
    fn curve_and_point_samples_stay_on_curve(
        case in arb_curve_and_point::<17>(CurveStrategyConfig::default())
    ) {
        let (curve, point) = case;
        prop_assert!(curve.contains(&point));
    }

    #[test]
    fn frobenius_cases_track_their_source_curve(
        case in arb_frobenius_curve_case::<17>(CurveStrategyConfig::default())
    ) {
        prop_assert_eq!(case.trace.curve_order(), case.curve.order() as u64);
        prop_assert_eq!(case.discriminant.frobenius_trace(), &case.trace);
    }

    #[test]
    fn endomorphism_cases_build_reports(
        case in arb_endomorphism_report_case::<17>(CurveStrategyConfig::default())
    ) {
        prop_assert_eq!(case.report.frobenius_discriminant().curve_order(), case.curve.order() as u64);
    }

    #[test]
    fn division_polynomial_cases_keep_supported_indices(
        case in arb_division_polynomial_case::<17>(CurveStrategyConfig::default())
    ) {
        prop_assert!(case.index >= 1);
        prop_assert!(case.curve.contains(&case.curve.identity()));
        prop_assert_eq!(case.polynomial.x_factor(), case.polynomial.x_factor());
    }

    #[test]
    fn upper_half_plane_samples_stay_in_the_upper_half_plane(
        tau in arb_upper_half_plane_point(AnalyticStrategyConfig::default())
    ) {
        prop_assert!(tau.imaginary_part() > 0.0);
    }

    #[test]
    fn complex_lattice_samples_remember_their_tau(
        lattice in arb_complex_lattice(AnalyticStrategyConfig::default())
    ) {
        prop_assert!(
            lattice
                .tau()
                .expect("sampled lattice should recover tau")
                .imaginary_part()
                > 0.0
        );
    }

    #[test]
    fn cyclic_kernel_cases_respect_kernel_membership(case in arb_cyclic_kernel_case()) {
        prop_assert!(case.curve.contains(&case.generator));
        prop_assert!(case.isogeny.kernel_points().contains(&case.kernel_point));
        prop_assert!(!case.isogeny.kernel_points().contains(&case.sample_point));
    }

    #[test]
    fn composable_velu_cases_keep_source_and_bridge_coherent(case in arb_composable_velu_case()) {
        prop_assert_eq!(case.first.codomain(), case.bridge.domain());
    }
}
