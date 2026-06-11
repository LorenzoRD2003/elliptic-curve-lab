use proptest::prelude::*;

use crate::elliptic_curves::{CurveIsomorphism, CurveModel, EnumerableCurveModel};
use crate::fields::{AmbientField, EnumerableFiniteField, ExtensionField, Field, Fp};
use crate::isogenies::Isogeny;
use crate::proptest_support::config::{
    AnalyticStrategyConfig, CurveStrategyConfig, FieldStrategyConfig, PolynomialStrategyConfig,
};
use crate::proptest_support::elliptic_curves::{
    arb_complex_lattice, arb_curve_and_point, arb_division_polynomial_case,
    arb_endomorphism_report_case, arb_frobenius_curve_case, arb_nonsingular_curve,
    arb_short_weierstrass_function_case, arb_short_weierstrass_function_pair_case,
    arb_upper_half_plane_point,
};
use crate::proptest_support::fields::{
    ProptestF17Sqrt3Spec, arb_complex_approx, arb_distinct_fp_elems, arb_extension_elem,
    arb_nonzero_fp_elem, arb_q_elem, arb_rational_function,
};
use crate::proptest_support::isogenies::{
    arb_composable_short_weierstrass_function_field_map_case,
    arb_short_weierstrass_function_field_map_case,
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
    fn rational_function_samples_keep_nonzero_monic_denominators(
        function in arb_rational_function::<F17>(PolynomialStrategyConfig::default())
    ) {
        prop_assert!(!function.denominator().is_zero());
        prop_assert!(function.denominator().is_monic());
    }

    #[test]
    fn rational_function_samples_respect_polynomial_size_budget_after_reduction(
        function in arb_rational_function::<F17>(PolynomialStrategyConfig::default())
    ) {
        prop_assert!(function.numerator().len() <= PolynomialStrategyConfig::default().max_len);
        prop_assert!(function.denominator().len() <= PolynomialStrategyConfig::default().max_len);
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
    fn function_field_cases_keep_curve_and_function_context_coherent(
        case in arb_short_weierstrass_function_case::<17>(
            CurveStrategyConfig::default(),
            PolynomialStrategyConfig::default(),
        )
    ) {
        prop_assert!(F17::eq(case.curve.a(), case.field.curve().a()));
        prop_assert!(F17::eq(case.curve.b(), case.field.curve().b()));
        prop_assert!(F17::eq(case.curve.a(), case.function.curve().a()));
        prop_assert!(F17::eq(case.curve.b(), case.function.curve().b()));
    }

    #[test]
    fn function_field_pair_cases_support_same_curve_operations(
        case in arb_short_weierstrass_function_pair_case::<17>(
            CurveStrategyConfig::default(),
            PolynomialStrategyConfig::default(),
        )
    ) {
        prop_assert!(case.left.add(&case.right).is_ok());
        prop_assert!(case.left.mul(&case.right).is_ok());
        prop_assert!(AmbientField::add(&case.field, &case.left, &case.right).is_ok());
        prop_assert!(AmbientField::mul(&case.field, &case.left, &case.right).is_ok());
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

    #[test]
    fn function_field_map_cases_keep_pullbacks_and_ambient_fields_coherent(
        case in arb_short_weierstrass_function_field_map_case::<17>(CurveStrategyConfig::default())
    ) {
        let domain_field = case.map.domain_function_field();
        let codomain_field = case.map.codomain_function_field();
        prop_assert_eq!(case.map.domain_curve(), &case.domain_curve);
        prop_assert_eq!(case.map.codomain_curve(), &case.codomain_curve);
        prop_assert_eq!(domain_field.curve(), &case.domain_curve);
        prop_assert_eq!(codomain_field.curve(), &case.codomain_curve);

        let codomain_x = case.codomain_field.x();
        let codomain_y = case.codomain_field.y();
        let pulled_x = case.map.pullback_function(&codomain_x);
        let pulled_y = case.map.pullback_function(&codomain_y);

        prop_assert_eq!(pulled_x, Ok(case.map.x_pullback().clone()));
        prop_assert_eq!(pulled_y, Ok(case.map.y_pullback().clone()));
    }

    #[test]
    fn composable_function_field_map_cases_agree_with_explicit_generator_pullback(
        case in arb_composable_short_weierstrass_function_field_map_case::<17>(CurveStrategyConfig::default())
    ) {
        prop_assert_eq!(case.first.codomain_curve(), &case.middle_curve);
        prop_assert_eq!(case.second.domain_curve(), &case.middle_curve);
        prop_assert_eq!(case.composite.domain_curve(), &case.domain_curve);
        prop_assert_eq!(case.composite.codomain_curve(), &case.codomain_curve);

        let expected_x = case
            .first
            .pullback_function(case.second.x_pullback())
            .expect("generated composable maps should pull back x coherently");
        let expected_y = case
            .first
            .pullback_function(case.second.y_pullback())
            .expect("generated composable maps should pull back y coherently");

        prop_assert_eq!(case.composite.x_pullback(), &expected_x);
        prop_assert_eq!(case.composite.y_pullback(), &expected_y);
    }
}
