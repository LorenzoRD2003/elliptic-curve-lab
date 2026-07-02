use num_bigint::{BigInt, BigUint};

use crate::elliptic_curves::{
    ShortWeierstrassCurve,
    group_algorithms::cyclic_roots::{
        CyclicPrimeRootBezout, CyclicPrimeRootCurveModel, CyclicPrimeRootError,
        CyclicPrimeRootInput, CyclicPrimeRootInputError, CyclicPrimeRootOutcome,
        CyclicPrimeRootReport, CyclicPrimeRootStep, CyclicPrimeRootTrace,
    },
    traits::{AffineCurveModel, CurveModel, GroupCurveModel},
};
use crate::fields::{Fp5, traits::Field};

fn bu(value: u8) -> BigUint {
    BigUint::from(value)
}

fn cyclic_f5_curve() -> ShortWeierstrassCurve<Fp5> {
    ShortWeierstrassCurve::<Fp5>::new(Fp5::zero(), Fp5::one())
        .expect("y² = x³ + 1 should be smooth over F5")
}

#[test]
fn input_records_prime_power_decomposition() {
    let input = CyclicPrimeRootInput::from_group_order_and_prime(bu(72), bu(3))
        .expect("72 = 8 * 3² should define a valid prime-root input");

    assert_eq!(input.group_order(), &bu(72));
    assert_eq!(input.root_degree(), &bu(3));
    assert_eq!(input.prime_to_root_cofactor(), &bu(8));
    assert_eq!(input.sylow_order(), &bu(9));
    assert_eq!(input.sylow_exponent(), 2);
    assert!(input.root_degree_divides_group_order());
}

#[test]
fn input_records_the_trivial_sylow_case() {
    let input = CyclicPrimeRootInput::from_group_order_and_prime(bu(35), bu(2))
        .expect("2 is prime even when it does not divide the group order");

    assert_eq!(input.group_order(), &bu(35));
    assert_eq!(input.root_degree(), &bu(2));
    assert_eq!(input.prime_to_root_cofactor(), &bu(35));
    assert_eq!(input.sylow_order(), &bu(1));
    assert_eq!(input.sylow_exponent(), 0);
    assert!(!input.root_degree_divides_group_order());
}

#[test]
fn input_rejects_zero_group_order() {
    assert_eq!(
        CyclicPrimeRootInput::from_group_order_and_prime(bu(0), bu(2)),
        Err(CyclicPrimeRootInputError::ZeroGroupOrder)
    );
}

#[test]
fn input_rejects_non_prime_root_degree() {
    assert_eq!(
        CyclicPrimeRootInput::from_group_order_and_prime(bu(72), bu(4)),
        Err(CyclicPrimeRootInputError::NonPrimeRootDegree { root_degree: bu(4) })
    );
}

#[test]
fn bezout_records_coefficients_for_root_formula() {
    let bezout = CyclicPrimeRootBezout::new(BigInt::from(-1), BigInt::from(1), bu(8), bu(27));

    assert_eq!(bezout.s(), &BigInt::from(-1));
    assert_eq!(bezout.t(), &BigInt::from(1));
    assert_eq!(bezout.cofactor(), &bu(8));
    assert_eq!(bezout.next_sylow_order(), &bu(27));
}

#[test]
fn report_preserves_target_generator_trace_and_outcome() {
    let input = CyclicPrimeRootInput::from_group_order_and_prime(bu(21), bu(2))
        .expect("2 is prime and does not divide 21");
    let step = CyclicPrimeRootStep::new(bu(1), "δ".to_string());
    let trace = CyclicPrimeRootTrace::new(
        Some("α".to_string()),
        Some("β".to_string()),
        Some(bu(1)),
        None,
        vec![step],
    );
    let outcome = CyclicPrimeRootOutcome::Root {
        root: "ρ".to_string(),
    };

    let report =
        CyclicPrimeRootReport::new(input, "γ".to_string(), "δ".to_string(), trace, outcome);

    assert_eq!(report.target(), "γ");
    assert_eq!(report.sylow_generator(), "δ");
    assert_eq!(report.input().group_order(), &bu(21));
    assert_eq!(report.trace().alpha(), Some(&"α".to_string()));
    assert_eq!(report.trace().beta(), Some(&"β".to_string()));
    assert_eq!(report.trace().discrete_log(), Some(&bu(1)));
    assert_eq!(report.trace().bezout(), None);
    assert_eq!(report.trace().steps()[0].discrete_log_candidate(), &bu(1));
    assert_eq!(report.trace().steps()[0].candidate_multiple(), "δ");
    assert_eq!(report.root(), Some(&"ρ".to_string()));
    assert!(report.outcome().found_root());
}

#[test]
fn no_root_outcome_has_no_root() {
    let outcome = CyclicPrimeRootOutcome::<String>::NoRoot;

    assert_eq!(outcome.root(), None);
    assert!(!outcome.found_root());
}

#[test]
fn cyclic_prime_root_uses_inverse_when_root_degree_does_not_divide_group_order() {
    let curve = cyclic_f5_curve();
    let generator = curve
        .point(Fp5::from_i64(2), Fp5::from_i64(2))
        .expect("(2,2) generates the cyclic order-6 curve");
    let target = curve
        .mul_scalar(&generator, bu(2))
        .expect("target should be a scalar multiple");

    let report = curve
        .cyclic_prime_root(&target, bu(5), bu(6), &curve.identity())
        .expect("5 should be invertible modulo the group order");

    let root = report.root().expect("unique root should exist");
    assert_eq!(
        curve.mul_scalar(root, bu(5)).expect("root should multiply"),
        target
    );
    assert_eq!(report.input().prime_to_root_cofactor(), &bu(6));
    assert_eq!(report.input().sylow_order(), &bu(1));
    assert!(report.trace().steps().is_empty());
}

#[test]
fn cyclic_prime_root_finds_root_inside_nontrivial_sylow_case() {
    let curve = cyclic_f5_curve();
    let generator = curve
        .point(Fp5::from_i64(2), Fp5::from_i64(2))
        .expect("(2,2) generates the cyclic order-6 curve");
    let target = curve
        .mul_scalar(&generator, bu(2))
        .expect("[2]P should be valid");
    let sylow_generator = curve
        .mul_scalar(&generator, bu(3))
        .expect("[3]P should generate the 2-Sylow subgroup");

    let report = curve
        .cyclic_prime_root(&target, bu(2), bu(6), &sylow_generator)
        .expect("target [2]P should have a square root P");

    let root = report.root().expect("root should exist");
    assert_eq!(
        curve.mul_scalar(root, bu(2)).expect("root should multiply"),
        target
    );
    assert_eq!(report.trace().discrete_log(), Some(&bu(2)));
    assert_eq!(
        report.trace().bezout().map(CyclicPrimeRootBezout::cofactor),
        Some(&bu(3))
    );
}

#[test]
fn cyclic_prime_root_reports_no_root_when_sylow_log_is_not_divisible_by_r() {
    let curve = cyclic_f5_curve();
    let generator = curve
        .point(Fp5::from_i64(2), Fp5::from_i64(2))
        .expect("(2,2) generates the cyclic order-6 curve");
    let sylow_generator = curve
        .mul_scalar(&generator, bu(3))
        .expect("[3]P should generate the 2-Sylow subgroup");

    let report = curve
        .cyclic_prime_root(&generator, bu(2), bu(6), &sylow_generator)
        .expect("valid setup should produce a no-root report");

    assert_eq!(report.root(), None);
    assert_eq!(report.outcome(), &CyclicPrimeRootOutcome::NoRoot);
    assert_eq!(report.trace().discrete_log(), Some(&bu(1)));
    assert_eq!(report.trace().bezout(), None);
}

#[test]
fn cyclic_prime_root_rejects_wrong_sylow_generator_order() {
    let curve = cyclic_f5_curve();
    let generator = curve
        .point(Fp5::from_i64(2), Fp5::from_i64(2))
        .expect("(2,2) generates the cyclic order-6 curve");

    assert_eq!(
        curve.cyclic_prime_root(&generator, bu(2), bu(6), &curve.identity()),
        Err(CyclicPrimeRootError::InvalidSylowGenerator {
            expected_order: bu(2)
        })
    );
}
