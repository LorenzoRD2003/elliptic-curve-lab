use num_complex::Complex64;
use num_rational::BigRational;
use proptest::prelude::*;

use super::{
    EisensteinSeriesQExpansion, EisensteinSeriesQExpansionApprox, EisensteinSeriesWeight,
    JInvariantQExpansion, JInvariantQExpansionApprox, ModularQExpansionApproximation,
    ModularQExpansionCoefficients, ModularQExpansionFamily, ModularQParameter,
    QExpansionTruncation, compare_j_from_eisenstein_and_q_expansion,
};
use crate::elliptic_curves::analytic::{
    AnalyticCurveError, ComplexLattice, UpperHalfPlanePoint, g4_sum, g6_sum,
};
use crate::elliptic_curves::{ApproxTolerance, LatticeSumTruncation};
use crate::fields::ComplexApprox;

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1.0e-12
}

fn q(numerator: i64, denominator: i64) -> BigRational {
    BigRational::new(numerator.into(), denominator.into())
}

fn e4_from_lattice_g4(sum: Complex64) -> Complex64 {
    sum * Complex64::new(45.0 / std::f64::consts::PI.powi(4), 0.0)
}

fn e6_from_lattice_g6(sum: Complex64) -> Complex64 {
    sum * Complex64::new(945.0 / (2.0 * std::f64::consts::PI.powi(6)), 0.0)
}

#[test]
fn q_parameter_stores_tau_and_computed_q() {
    let tau = UpperHalfPlanePoint::tau_i();
    let q = ModularQParameter::from_tau(tau.clone());

    assert_eq!(q.tau(), &tau);
    assert_eq!(q.q(), &Complex64::new((-std::f64::consts::TAU).exp(), 0.0));
}

#[test]
fn absolute_value_matches_exp_minus_two_pi_im_tau() {
    let tau = UpperHalfPlanePoint::from_re_im(0.3, 1.2).unwrap();
    let q = ModularQParameter::from_tau(tau.clone());
    let expected = (-std::f64::consts::TAU * tau.imaginary_part()).exp();

    assert!(close(q.absolute_value(), expected));
    assert!(q.absolute_value() < 1.0);
}

#[test]
fn tau_rho_produces_the_expected_real_negative_q() {
    let q = ModularQParameter::from_tau(UpperHalfPlanePoint::tau_rho());
    let expected_real = -(std::f64::consts::PI * f64::sqrt(3.0)).exp().recip();

    assert!(close(q.q().re, expected_real));
    assert!(q.q().im.abs() <= 1.0e-12);
}

#[test]
fn derived_type_is_cloneable_and_structurally_stable() {
    let q = ModularQParameter::from_tau(UpperHalfPlanePoint::tau_generic_example());
    let cloned: ModularQParameter = q.clone();

    assert_eq!(cloned, q);
}

#[test]
fn q_expansion_truncation_rejects_zero_and_accepts_positive_lengths() {
    assert_eq!(
        QExpansionTruncation::new(0),
        Err(AnalyticCurveError::InvalidSeriesPrecision)
    );

    let truncation = QExpansionTruncation::new(12).unwrap();
    assert_eq!(truncation.terms(), 12);
}

#[test]
fn q_expansion_truncation_default_is_explicit_and_j_table_has_its_own_cap() {
    assert_eq!(QExpansionTruncation::default_educational().terms(), 4);
    assert_eq!(
        JInvariantQExpansion::full_current_table_truncation().terms(),
        6
    );
    assert_eq!(JInvariantQExpansion::max_supported_terms(), 6);
}

#[test]
fn modular_q_expansion_coefficients_store_and_expose_a_small_table() {
    let coefficients = ModularQExpansionCoefficients::from_integers(-1, vec![1, 2, 3]);

    assert_eq!(coefficients.start_exponent(), -1);
    assert_eq!(coefficients.end_exponent(), Some(1));
    assert_eq!(coefficients.coefficients(), &[q(1, 1), q(2, 1), q(3, 1)]);
    assert_eq!(coefficients.coefficient_of(-1), Some(q(1, 1)));
    assert_eq!(coefficients.coefficient_of(0), Some(q(2, 1)));
    assert_eq!(coefficients.coefficient_of(1), Some(q(3, 1)));
    assert_eq!(coefficients.coefficient_of(2), None);
    assert_eq!(coefficients.len(), 3);
    assert!(!coefficients.is_empty());
}

#[test]
fn coefficient_table_evaluates_using_start_exponent() {
    let coefficients = ModularQExpansionCoefficients::from_integers(-1, vec![1, 2, 3]);
    let q = Complex64::new(0.5, 0.0);
    let truncation = QExpansionTruncation::new(3).unwrap();

    let expected =
        Complex64::new(1.0, 0.0) / q + Complex64::new(2.0, 0.0) + Complex64::new(3.0, 0.0) * q;

    assert!(
        (coefficients.evaluate_truncated_at(q, truncation).unwrap() - expected).norm() <= 1.0e-12
    );
}

#[test]
fn coefficient_table_evaluation_respects_truncation_length() {
    let coefficients = ModularQExpansionCoefficients::from_integers(0, vec![10, 20, 30]);
    let q = Complex64::new(0.25, 0.0);
    let truncation = QExpansionTruncation::new(2).unwrap();
    let expected = Complex64::new(10.0, 0.0) + Complex64::new(20.0, 0.0) * q;

    assert!(
        (coefficients.evaluate_truncated_at(q, truncation).unwrap() - expected).norm() <= 1.0e-12
    );
}

#[test]
fn coefficient_table_truncation_reports_when_request_exceeds_table_length() {
    let coefficients = ModularQExpansionCoefficients::from_integers(0, vec![10, 20, 30]);

    assert_eq!(
        coefficients.truncated(QExpansionTruncation::new(4).unwrap()),
        Err(AnalyticCurveError::InvalidSeriesPrecision)
    );
}

#[test]
fn j_invariant_coefficients_live_in_the_shared_coefficients_value_object() {
    let coefficients = ModularQExpansionCoefficients::j_invariant_nonnegative();

    assert_eq!(coefficients.start_exponent(), 0);
    assert_eq!(coefficients.end_exponent(), Some(4));
    assert_eq!(
        coefficients.coefficients(),
        &[
            q(744, 1),
            q(196_884, 1),
            q(21_493_760, 1),
            q(864_299_970, 1),
            q(20_245_856_256, 1),
        ]
    );
    assert_eq!(coefficients.coefficient_of(-1), None);
    assert_eq!(coefficients.coefficient_of(0), Some(q(744, 1)));
    assert_eq!(coefficients.coefficient_of(4), Some(q(20_245_856_256, 1)));
}

#[test]
fn exported_j_q_coefficients_match_the_nonnegative_tail_table() {
    assert_eq!(
        ModularQExpansionCoefficients::j_invariant_nonnegative().coefficients(),
        &[
            q(744, 1),
            q(196_884, 1),
            q(21_493_760, 1),
            q(864_299_970, 1),
            q(20_245_856_256, 1),
        ]
    );
}

#[test]
fn j_invariant_current_table_includes_the_principal_term() {
    let coefficients = ModularQExpansionCoefficients::j_invariant_current_table();

    assert_eq!(coefficients.start_exponent(), -1);
    assert_eq!(coefficients.end_exponent(), Some(4));
    assert_eq!(
        coefficients.coefficients(),
        &[
            q(1, 1),
            q(744, 1),
            q(196_884, 1),
            q(21_493_760, 1),
            q(864_299_970, 1),
            q(20_245_856_256, 1),
        ]
    );
    assert_eq!(coefficients.coefficient_of(-1), Some(q(1, 1)));
    assert_eq!(coefficients.coefficient_of(0), Some(q(744, 1)));
}

#[test]
fn eisenstein_e4_coefficients_follow_the_sigma_three_formula() {
    let coefficients = EisensteinSeriesQExpansion::e4()
        .coefficients(QExpansionTruncation::new(5).unwrap())
        .unwrap();

    assert_eq!(coefficients.start_exponent(), 0);
    assert_eq!(
        coefficients.coefficients(),
        &[q(1, 1), q(240, 1), q(2160, 1), q(6720, 1), q(17_520, 1)]
    );
}

#[test]
fn eisenstein_e6_coefficients_follow_the_sigma_five_formula() {
    let coefficients = EisensteinSeriesQExpansion::e6()
        .coefficients(QExpansionTruncation::new(5).unwrap())
        .unwrap();

    assert_eq!(coefficients.start_exponent(), 0);
    assert_eq!(
        coefficients.coefficients(),
        &[
            q(1, 1),
            q(-504, 1),
            q(-16_632, 1),
            q(-122_976, 1),
            q(-532_728, 1),
        ]
    );
}

#[test]
fn j_invariant_report_stores_q_parameter_value_truncation_and_term_count() {
    let approximation = JInvariantQExpansion::from_tau(
        UpperHalfPlanePoint::tau_i(),
        QExpansionTruncation::new(3).unwrap(),
    )
    .unwrap();

    assert_eq!(approximation.terms_used(), 3);
    assert_eq!(
        approximation.truncation(),
        QExpansionTruncation::new(3).unwrap()
    );
    assert_eq!(approximation.tau(), &UpperHalfPlanePoint::tau_i());
    assert_eq!(approximation.q(), approximation.q_parameter().q());
}

#[test]
fn one_term_j_expansion_matches_only_the_principal_term() {
    let tau = UpperHalfPlanePoint::tau_i();
    let q_parameter = ModularQParameter::from_tau(tau.clone());
    let expected = Complex64::new(1.0, 0.0) / *q_parameter.q();
    let approximation =
        JInvariantQExpansion::from_tau(tau, QExpansionTruncation::new(1).unwrap()).unwrap();

    assert!((*approximation.value() - expected).norm() <= 1.0e-12);
}

#[test]
fn two_term_j_expansion_adds_the_constant_744_term() {
    let tau = UpperHalfPlanePoint::tau_i();
    let q_parameter = ModularQParameter::from_tau(tau.clone());
    let expected = Complex64::new(1.0, 0.0) / *q_parameter.q() + Complex64::new(744.0, 0.0);
    let approximation =
        JInvariantQExpansion::from_tau(tau, QExpansionTruncation::new(2).unwrap()).unwrap();

    assert!((*approximation.value() - expected).norm() <= 1.0e-12);
}

#[test]
fn additional_terms_change_the_j_q_expansion_value_for_generic_tau() {
    let tau = UpperHalfPlanePoint::tau_generic_example();
    let small =
        JInvariantQExpansion::from_tau(tau.clone(), QExpansionTruncation::new(1).unwrap()).unwrap();
    let larger =
        JInvariantQExpansion::from_tau(tau, JInvariantQExpansion::full_current_table_truncation())
            .unwrap();

    assert_ne!(small.value(), larger.value());
}

#[test]
fn approximation_report_is_cloneable_and_structurally_stable() {
    let approximation = JInvariantQExpansion::from_tau(
        UpperHalfPlanePoint::tau_i(),
        QExpansionTruncation::default_educational(),
    )
    .unwrap();
    let cloned: JInvariantQExpansionApprox = approximation.clone();

    assert_eq!(cloned, approximation);
}

#[test]
fn shared_approximation_trait_exposes_tau_q_and_terms_consistently() {
    fn check_family<F: ModularQExpansionFamily>(
        family: &F,
        tau: UpperHalfPlanePoint,
        truncation: QExpansionTruncation,
    ) where
        F::Approximation: ModularQExpansionApproximation,
    {
        let approximation = family.evaluate_tau(tau.clone(), truncation).unwrap();

        assert_eq!(approximation.tau(), &tau);
        assert_eq!(approximation.q(), approximation.q_parameter().q());
        assert_eq!(approximation.terms_used(), truncation.terms());
        assert_eq!(approximation.truncation(), truncation);
        assert!(approximation.value().re.is_finite());
        assert!(approximation.value().im.is_finite());
    }

    let tau = UpperHalfPlanePoint::tau_i();
    let truncation = QExpansionTruncation::new(4).unwrap();

    check_family(&JInvariantQExpansion::new(), tau.clone(), truncation);
    check_family(&EisensteinSeriesQExpansion::e4(), tau.clone(), truncation);
    check_family(&EisensteinSeriesQExpansion::e6(), tau, truncation);
}

#[test]
fn eisenstein_runtime_family_records_the_requested_weight() {
    let family = EisensteinSeriesQExpansion::new(EisensteinSeriesWeight::new(8).unwrap());
    let approximation: EisensteinSeriesQExpansionApprox = family
        .evaluate_tau(
            UpperHalfPlanePoint::tau_i(),
            QExpansionTruncation::new(4).unwrap(),
        )
        .unwrap();

    assert_eq!(family.weight(), EisensteinSeriesWeight::new(8).unwrap());
    assert_eq!(
        approximation.weight(),
        EisensteinSeriesWeight::new(8).unwrap()
    );
    assert_eq!(approximation.k(), 8);
    assert_eq!(approximation.terms_used(), 4);
}

#[test]
fn j_invariant_comparison_report_records_both_routes_and_their_difference() {
    let tau = UpperHalfPlanePoint::tau_i();
    let lattice_truncation = LatticeSumTruncation::default_educational();
    let q_truncation = QExpansionTruncation::new(3).unwrap();
    let report = compare_j_from_eisenstein_and_q_expansion(
        tau.clone(),
        lattice_truncation,
        q_truncation,
        ApproxTolerance::loose(),
    )
    .unwrap();

    assert_eq!(report.tau(), &tau);
    assert_eq!(
        *report.difference(),
        *report.eisenstein_j() - *report.q_expansion_j()
    );
    assert_eq!(report.lattice_truncation(), lattice_truncation);
    assert_eq!(report.q_truncation(), q_truncation);
    assert_eq!(report.tolerance(), ApproxTolerance::loose());
    assert_eq!(
        report.agrees_approximately(),
        ComplexApprox::eq_with_tolerance(
            report.eisenstein_j(),
            report.q_expansion_j(),
            report.tolerance(),
        )
    );
}

#[test]
fn square_lattice_comparison_can_agree_under_richer_truncations() {
    let report = compare_j_from_eisenstein_and_q_expansion(
        UpperHalfPlanePoint::tau_i(),
        LatticeSumTruncation::larger_for_comparison(),
        JInvariantQExpansion::full_current_table_truncation(),
        ApproxTolerance::new(1.0e-1, 1.0e-8),
    )
    .unwrap();

    assert!(report.agrees_approximately());
}

#[test]
fn j_q_expansion_matches_eisenstein_j_for_tau_with_large_imaginary_part() {
    let tau = UpperHalfPlanePoint::from_re_im(0.0, 1.5).unwrap();
    let report = compare_j_from_eisenstein_and_q_expansion(
        tau,
        LatticeSumTruncation::new(64).unwrap(),
        JInvariantQExpansion::full_current_table_truncation(),
        ApproxTolerance::new(1.0, 1.0e-6),
    )
    .unwrap();

    let relative_error = report.absolute_difference() / report.eisenstein_j().norm();
    assert!(report.absolute_difference().is_finite());
    assert!(relative_error < 5.0e-4);
}

#[test]
fn e4_q_expansion_matches_lattice_eisenstein_approximately() {
    let tau = UpperHalfPlanePoint::from_re_im(0.1, 2.5).unwrap();
    let lattice = ComplexLattice::from_tau(tau.clone());
    let lattice_value = e4_from_lattice_g4(
        g4_sum(&lattice, LatticeSumTruncation::new(16).unwrap())
            .unwrap()
            .value,
    );
    let q_value = *EisensteinSeriesQExpansion::e4()
        .evaluate_tau(tau, QExpansionTruncation::new(5).unwrap())
        .unwrap()
        .value();

    assert!(ComplexApprox::eq_with_tolerance(
        &lattice_value,
        &q_value,
        ApproxTolerance::new(1.0e-3, 1.0e-3)
    ));
}

#[test]
fn e6_q_expansion_matches_lattice_eisenstein_approximately() {
    let tau = UpperHalfPlanePoint::from_re_im(0.1, 2.5).unwrap();
    let lattice = ComplexLattice::from_tau(tau.clone());
    let lattice_value = e6_from_lattice_g6(
        g6_sum(&lattice, LatticeSumTruncation::new(16).unwrap())
            .unwrap()
            .value,
    );
    let q_value = *EisensteinSeriesQExpansion::e6()
        .evaluate_tau(tau, QExpansionTruncation::new(5).unwrap())
        .unwrap()
        .value();

    assert!(ComplexApprox::eq_with_tolerance(
        &lattice_value,
        &q_value,
        ApproxTolerance::new(1.0e-3, 1.0e-3)
    ));
}

// Numerical note from the first generic-j comparison experiments:
//
// A stronger proptest that demanded approximate agreement for generic `τ`
// was too ambitious with the current short `q`-table and moderate Eisenstein
// radii. Two concrete counterexamples improved steadily once the lattice
// truncation radius `r` was increased:
//
// | τ                         | r=8 rel error | r=16 rel error | r=32 rel error |
// |---------------------------|---------------|----------------|----------------|
// | `1.5990124441568718 i`    | `2.56e-2`     | `6.79e-3`      | `1.75e-3`      |
// | `1.9655657131426314 i`    | `1.44e-1`     | `3.81e-2`      | `9.83e-3`      |
//
// So the property we keep here is the honest one for the current layer:
// generic `τ` should produce finite, self-consistent comparison reports, while
// stronger agreement claims remain reserved for selected well-behaved cases
// such as `τ = i`.
proptest! {
    #![proptest_config(ProptestConfig::with_cases(24))]

    #[test]
    fn generic_taus_produce_finite_and_self_consistent_j_comparison_reports(
        re in -0.45f64..0.45,
        im in 0.8f64..2.2,
    ) {
        let tau = UpperHalfPlanePoint::from_re_im(re, im).unwrap();
        let report = compare_j_from_eisenstein_and_q_expansion(
            tau,
            LatticeSumTruncation::new(8).unwrap(),
            JInvariantQExpansion::full_current_table_truncation(),
            ApproxTolerance::loose(),
        ).unwrap();

        prop_assert!(report.eisenstein_j().re.is_finite());
        prop_assert!(report.eisenstein_j().im.is_finite());
        prop_assert!(report.q_expansion_j().re.is_finite());
        prop_assert!(report.q_expansion_j().im.is_finite());
        prop_assert!(report.difference().re.is_finite());
        prop_assert!(report.difference().im.is_finite());
        prop_assert!(report.absolute_difference().is_finite());
        prop_assert_eq!(
            *report.difference(),
            *report.eisenstein_j() - *report.q_expansion_j()
        );
        prop_assert_eq!(
            report.agrees_approximately(),
            ComplexApprox::eq_with_tolerance(
                report.eisenstein_j(),
                report.q_expansion_j(),
                report.tolerance(),
            )
        );
    }
}
