use num_bigint::BigUint;

use crate::elliptic_curves::{
    CurveError, ShortWeierstrassCurve,
    frobenius::{
        HasseInterval, SchoofGroupOrderOutcome, SchoofTraceCrtOutcome, SchoofTraceMod2Report,
        SchoofTraceModOddPrimeOutcome, schoof::finalize_schoof_group_order_report,
    },
    short_weierstrass::division_polynomials::DivisionPolynomialError,
    traits::{EnumerableCurveModel, FrobeniusTraceCurveModel},
};
use crate::fields::{
    Fp,
    finite_field_descriptor::FiniteFieldDescriptor,
    traits::{Field, FiniteField},
};
use crate::polynomials::DensePolynomial;

type F7 = Fp<7>;
type F19 = Fp<19>;
type F43 = Fp<43>;

fn expected_trace_mod_2_from_order(curve: &ShortWeierstrassCurve<F7>) -> u8 {
    (curve.order() % 2) as u8
}

#[test]
fn schoof_mod_2_matches_exhaustive_trace_and_group_order_parity_over_f7() {
    let curve =
        ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3)).expect("valid curve");

    let report = curve.schoof_trace_mod_2();
    let exhaustive_trace = curve
        .frobenius_trace()
        .expect("small enumerable curve should supply a Frobenius trace");

    assert_eq!(report.field_order(), 7);
    assert_eq!(
        report.trace_mod_2(),
        expected_trace_mod_2_from_order(&curve)
    );
    assert_eq!(
        report.trace_mod_2(),
        exhaustive_trace.trace().rem_euclid(2) as u8
    );
    assert_eq!(
        report.group_order_is_even(),
        curve.order().is_multiple_of(2)
    );
}

#[test]
fn schoof_mod_2_detects_rational_two_torsion_when_cubic_has_root() {
    let curve = ShortWeierstrassCurve::<F19>::new(F19::from_i64(-1), F19::zero())
        .expect("valid curve with rational 2-torsion");

    let report = curve.schoof_trace_mod_2();

    assert!(report.has_rational_two_torsion());
    assert_eq!(report.trace_mod_2(), 0);
    assert!(report.group_order_is_even());
}

#[test]
fn schoof_mod_2_report_exposes_the_intermediate_polynomials() {
    let curve =
        ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3)).expect("valid curve");

    let report: SchoofTraceMod2Report<F7> = curve.schoof_trace_mod_2();
    let x = DensePolynomial::new(vec![F7::zero(), F7::one()]);

    assert_eq!(report.cubic(), &curve.to_cubic());
    assert_eq!(
        report.x_q_mod_cubic(),
        &DensePolynomial::pow_mod(&x, 7, report.cubic())
            .expect("nonzero cubic should define a quotient ring")
    );
    assert_eq!(
        report.gcd(),
        &report.cubic().gcd(&report.x_q_mod_cubic().sub(&x))
    );
}

#[test]
fn schoof_odd_prime_step_can_recover_trace_mod_three_over_f7() {
    let curve = ShortWeierstrassCurve::<F7>::new(F7::zero(), F7::from_i64(2))
        .expect("valid curve with a direct ell = 3 trace step");

    let report = curve
        .schoof_trace_mod_odd_prime(3)
        .expect("ell = 3 should be a valid odd-prime Schoof input");
    let exhaustive_trace = curve
        .frobenius_trace()
        .expect("small enumerable curve should supply a Frobenius trace");

    assert_eq!(report.field_order(), 7);
    assert_eq!(report.odd_prime(), 3);
    assert_eq!(
        report.division_polynomial(),
        curve
            .division_polynomial(3)
            .expect("psi_3")
            .stored_x_factor()
    );

    let SchoofTraceModOddPrimeOutcome::TraceFound { trace_mod_ell } = report.outcome() else {
        panic!("the chosen curve should complete the ell = 3 step without denominator refinement");
    };

    assert_eq!(
        *trace_mod_ell,
        exhaustive_trace.trace().rem_euclid(3) as usize
    );
    assert_eq!(report.trace_mod_odd_prime(), Some(*trace_mod_ell));
    assert_eq!(report.candidate_reports().len(), *trace_mod_ell + 1);
}

#[test]
fn schoof_odd_prime_step_rejects_invalid_ell_inputs() {
    let curve =
        ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3)).expect("valid curve");

    let error = curve
        .schoof_trace_mod_odd_prime(7)
        .expect_err("ell equal to the characteristic should be rejected");

    assert_eq!(
        error,
        DivisionPolynomialError::Curve(CurveError::InvalidSchoofOddPrime {
            odd_prime: 7,
            characteristic: 7,
        })
    );
}

#[test]
fn schoof_trace_crt_combines_mod_two_and_mod_three() {
    let curve = ShortWeierstrassCurve::<F7>::new(F7::zero(), F7::from_i64(2))
        .expect("valid curve with direct ell = 3 resolution");
    let exhaustive_trace = curve
        .frobenius_trace()
        .expect("small enumerable curve should supply a Frobenius trace");

    let report = curve
        .schoof_trace_crt(&[3])
        .expect("ell = 3 should be a valid Schoof prime input");

    let SchoofTraceCrtOutcome::Combined { solution } = report.outcome() else {
        panic!("the chosen curve should resolve both mod 2 and mod 3");
    };

    assert_eq!(report.resolved_congruences().len(), 2);
    assert_eq!(solution.modulus(), &BigUint::from(6u8));
    assert_eq!(
        solution.residue(),
        &BigUint::from(exhaustive_trace.trace().rem_euclid(6) as u64)
    );
    assert_eq!(
        report.mod_2_report().trace_mod_2(),
        exhaustive_trace.trace().rem_euclid(2) as u8
    );
    assert_eq!(report.odd_prime_reports().len(), 1);
}

#[test]
fn schoof_trace_crt_returns_the_partial_solution_when_an_odd_prime_blocks() {
    let curve =
        ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3)).expect("valid curve");

    let report = curve
        .schoof_trace_crt(&[3])
        .expect("ell = 3 should be a valid Schoof prime input");

    let SchoofTraceCrtOutcome::BlockedOnOddPrime {
        blocked_prime,
        partial_solution,
    } = report.outcome()
    else {
        panic!("the chosen curve should block on the odd-prime step without refinement");
    };

    assert_eq!(*blocked_prime, 3);
    assert_eq!(report.resolved_congruences().len(), 1);
    assert_eq!(report.odd_prime_reports().len(), 1);
    assert_eq!(partial_solution.modulus(), &BigUint::from(2u8));
    assert_eq!(
        partial_solution.residue(),
        &BigUint::from(report.mod_2_report().trace_mod_2())
    );
}

#[test]
fn schoof_group_order_reports_ambiguity_when_the_crt_modulus_is_too_small() {
    let curve =
        ShortWeierstrassCurve::<F7>::new(F7::zero(), F7::from_i64(2)).expect("valid small curve");
    let base_field = FiniteFieldDescriptor::new(F7::characteristic(), F7::extension_degree())
        .expect("finite field descriptor should be valid");
    let crt_report = curve
        .schoof_trace_crt(&[])
        .expect("the mod-2 Schoof step should always succeed");

    let report = finalize_schoof_group_order_report(base_field, crt_report)
        .expect("the final Hasse resolution should succeed");

    match report.outcome() {
        SchoofGroupOrderOutcome::AmbiguousTraceClass {
            first_trace,
            last_trace,
            candidate_count,
        } => {
            assert!(first_trace <= last_trace);
            assert!(*candidate_count >= 2);
        }
        other => panic!("expected an ambiguous trace class, got {other:?}"),
    }
}

#[test]
fn schoof_group_order_resolves_the_unique_trace_after_crt_and_hasse() {
    for a in -3..=3 {
        for b in -3..=3 {
            let Ok(curve) = ShortWeierstrassCurve::<F43>::new(F43::from_i64(a), F43::from_i64(b))
            else {
                continue;
            };
            let base_field =
                FiniteFieldDescriptor::new(F43::characteristic(), F43::extension_degree())
                    .expect("finite field descriptor should be valid");
            let crt_report = curve
                .schoof_trace_crt(&[3, 5])
                .expect("the requested Schoof primes should be valid inputs");
            let report = finalize_schoof_group_order_report(base_field, crt_report)
                .expect("the final Hasse resolution should succeed");

            let SchoofGroupOrderOutcome::GroupOrderFound { trace, curve_order } = report.outcome()
            else {
                continue;
            };

            let exhaustive_trace = curve
                .frobenius_trace()
                .expect("small enumerable curve should supply a Frobenius trace");
            assert_eq!(*trace, i128::from(exhaustive_trace.trace()));
            assert_eq!(*curve_order, u128::from(exhaustive_trace.curve_order()));
            assert_eq!(
                report
                    .to_frobenius_trace()
                    .expect("small resolved order should convert to FrobeniusTrace"),
                Some(exhaustive_trace)
            );
            return;
        }
    }

    panic!("expected to find one small F43 curve resolved by the current [3, 5] Schoof stage");
}

#[test]
fn schoof_trace_crt_until_hasse_uniqueness_preserves_a_blocked_prime_step() {
    let curve =
        ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3)).expect("valid curve");

    let report = curve
        .schoof_trace_crt_until_hasse_uniqueness()
        .expect("the automatic Schoof CRT driver should run");

    let SchoofTraceCrtOutcome::BlockedOnOddPrime {
        blocked_prime,
        partial_solution,
    } = report.outcome()
    else {
        panic!("the chosen curve should still block on the first odd-prime step");
    };

    assert_eq!(*blocked_prime, 3);
    assert_eq!(partial_solution.modulus(), &BigUint::from(2u8));
    assert_eq!(report.odd_prime_reports().len(), 1);
}

#[test]
fn schoof_group_order_uses_the_hasse_stopping_rule() {
    for a in -3..=3 {
        for b in -3..=3 {
            let Ok(curve) = ShortWeierstrassCurve::<F43>::new(F43::from_i64(a), F43::from_i64(b))
            else {
                continue;
            };
            let report = curve
                .schoof_group_order()
                .expect("the automatic Schoof route should run on valid curves");

            let SchoofGroupOrderOutcome::GroupOrderFound { trace, curve_order } = report.outcome()
            else {
                continue;
            };

            let exhaustive_trace = curve
                .frobenius_trace()
                .expect("small enumerable curve should supply a Frobenius trace");
            assert_eq!(*trace, i128::from(exhaustive_trace.trace()));
            assert_eq!(*curve_order, u128::from(exhaustive_trace.curve_order()));

            let combined_solution = report
                .crt_report()
                .combined_solution()
                .expect("automatic successful route should end at one combined CRT class");
            let threshold = BigUint::from(
                HasseInterval::for_q(report.field_order())
                    .expect("valid finite field order should define H(q)")
                    .trace_bound(),
            ) * BigUint::from(2u8);
            assert!(combined_solution.modulus() > &threshold);
            return;
        }
    }

    panic!("expected to find one small F43 curve resolved by the automatic Schoof route");
}
