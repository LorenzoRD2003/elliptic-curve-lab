use num_bigint::{BigInt, BigUint};

use crate::elliptic_curves::group_algorithms::cyclic_roots::{
    CyclicPrimeRootBezout, CyclicPrimeRootInput, CyclicPrimeRootInputError, CyclicPrimeRootOutcome,
    CyclicPrimeRootReport, CyclicPrimeRootStep, CyclicPrimeRootTrace,
};

fn bu(value: u8) -> BigUint {
    BigUint::from(value)
}

#[test]
fn input_records_prime_power_decomposition() {
    let input = CyclicPrimeRootInput::from_group_order_and_prime(bu(72), bu(3))
        .expect("72 = 8 * 3² should define a valid prime-root input");

    assert_eq!(input.group_order(), &bu(72));
    assert_eq!(input.root_prime(), &bu(3));
    assert_eq!(input.prime_to_root_cofactor(), &bu(8));
    assert_eq!(input.sylow_order(), &bu(9));
    assert_eq!(input.sylow_exponent(), 2);
    assert!(input.root_prime_divides_group_order());
}

#[test]
fn input_records_the_trivial_sylow_case() {
    let input = CyclicPrimeRootInput::from_group_order_and_prime(bu(35), bu(2))
        .expect("2 is prime even when it does not divide the group order");

    assert_eq!(input.group_order(), &bu(35));
    assert_eq!(input.root_prime(), &bu(2));
    assert_eq!(input.prime_to_root_cofactor(), &bu(35));
    assert_eq!(input.sylow_order(), &bu(1));
    assert_eq!(input.sylow_exponent(), 0);
    assert!(!input.root_prime_divides_group_order());
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
