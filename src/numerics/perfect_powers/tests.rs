use num_bigint::{BigInt, BigUint};

use crate::numerics::{
    hensel::HenselLiftError,
    perfect_powers::{
        PerfectPowerOutcome, PerfectPowerSearchConfig, detect_perfect_power,
        exponents::prime_exponents_through,
    },
};

fn bu(value: u64) -> BigUint {
    BigUint::from(value)
}

fn detected(input: u64) -> PerfectPowerOutcome {
    detect_perfect_power(&bu(input), PerfectPowerSearchConfig::default())
        .outcome()
        .clone()
}

#[test]
fn detects_odd_square_coprime_to_six() {
    let PerfectPowerOutcome::PerfectPower { base, exponent } = detected(25) else {
        panic!("25 should be detected as a square");
    };

    assert_eq!(base, bu(5));
    assert_eq!(exponent, 2);
}

#[test]
fn detects_odd_cube_coprime_to_six() {
    assert_eq!(
        detected(125),
        PerfectPowerOutcome::PerfectPower {
            base: bu(5),
            exponent: 3
        }
    );
}

#[test]
fn detects_higher_prime_exponent() {
    assert_eq!(
        detected(17u64.pow(5)),
        PerfectPowerOutcome::PerfectPower {
            base: bu(17),
            exponent: 5
        }
    );
}

#[test]
fn reports_not_perfect_power_in_scope() {
    let report = detect_perfect_power(&bu(101), PerfectPowerSearchConfig::default());

    assert_eq!(report.input(), &bu(101));
    assert_eq!(report.outcome(), &PerfectPowerOutcome::NotPerfectPower);
    assert!(!report.candidate_reports().is_empty());
    assert!(
        report
            .candidate_reports()
            .iter()
            .any(|candidate| candidate.exponent() == 2)
    );
}

#[test]
fn reports_inputs_not_coprime_to_six_as_out_of_scope() {
    assert_eq!(detected(64), PerfectPowerOutcome::NotCoprimeToSix);
    assert_eq!(detected(27), PerfectPowerOutcome::NotCoprimeToSix);
}

#[test]
fn reports_degenerate_inputs_as_out_of_scope() {
    assert_eq!(detected(1), PerfectPowerOutcome::DegenerateInput);
}

#[test]
fn candidate_report_records_hensel_prime_choice_and_bound() {
    let report = detect_perfect_power(&bu(25), PerfectPowerSearchConfig::default());
    let square_candidate = report
        .candidate_reports()
        .iter()
        .find(|candidate| candidate.exponent() == 2)
        .expect("square exponent should be tested first");

    assert_eq!(report.config().max_seed_scan(), 10_000);
    assert_eq!(square_candidate.hensel_prime(), &bu(3));
    assert!(square_candidate.root_bound() >= &bu(5));
    assert!(
        square_candidate
            .certified_roots()
            .contains(&BigInt::from(5))
    );
}

#[test]
fn config_controls_underlying_seed_scan_limit() {
    let config = PerfectPowerSearchConfig::default().with_max_seed_scan(2);
    let report = detect_perfect_power(&bu(25), config);

    assert_eq!(
        report.outcome(),
        &PerfectPowerOutcome::HenselFailure(HenselLiftError::SeedScanLimitExceeded)
    );
}

#[test]
fn prime_exponents_stop_at_floor_log2_input() {
    assert_eq!(prime_exponents_through(1), Vec::<u32>::new());
    assert_eq!(prime_exponents_through(10), vec![2, 3, 5, 7]);
}
