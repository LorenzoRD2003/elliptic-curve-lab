use num_bigint::{BigInt, BigUint};

use super::{
    HenselLiftError, hensel_lift_simple_root, hensel_lift_simple_root_step,
    hensel_lift_square_root, hensel_lift_square_root_fast, polynomial::evaluate_polynomial,
};

fn bi(value: i64) -> BigInt {
    BigInt::from(value)
}

fn bu(value: u64) -> BigUint {
    BigUint::from(value)
}

fn polynomial(values: &[i64]) -> Vec<BigInt> {
    values.iter().copied().map(bi).collect()
}

#[test]
fn simple_step_lifts_a_square_root_modulo_the_next_prime_power() {
    let coefficients = polynomial(&[-2, 0, 1]);
    let step = hensel_lift_simple_root_step(&coefficients, &bi(3), &bu(7), 1)
        .expect("3^2 = 2 mod 7 should lift uniquely");

    assert_eq!(step.level(), 1);
    assert_eq!(step.root_before(), &bi(3));
    assert_eq!(step.lift_digit(), &bu(1));
    assert_eq!(step.root_after(), &bi(10));
    assert_eq!(
        evaluate_polynomial(&coefficients, step.root_after()) % bi(49),
        bi(0)
    );
}

#[test]
fn repeated_lift_records_each_hensel_digit() {
    let coefficients = polynomial(&[-2, 0, 1]);
    let trace = hensel_lift_simple_root(&coefficients, &bi(3), &bu(7), 4)
        .expect("simple square root should lift through p^4");

    assert_eq!(trace.prime(), &bu(7));
    assert_eq!(trace.coefficients(), coefficients.as_slice());
    assert_eq!(trace.initial_root(), &bi(3));
    assert_eq!(trace.reached_level(), 4);
    assert_eq!(trace.final_root(), &bi(2166));
    assert_eq!(trace.steps().len(), 3);
    assert_eq!(
        trace
            .steps()
            .iter()
            .map(|step| step.lift_digit().clone())
            .collect::<Vec<_>>(),
        vec![bu(1), bu(2), bu(6)]
    );
    assert_eq!(
        evaluate_polynomial(&coefficients, trace.final_root()) % bi(2401),
        bi(0)
    );
}

#[test]
fn square_root_lift_delegates_to_the_simple_polynomial_route() {
    let a = bi(2);
    let p = bi(7);
    let trace = hensel_lift_square_root(&bi(2), &bi(3), &bu(7), 4)
        .expect("3^2 = 2 mod 7 should lift through p^4");

    assert_eq!(trace.coefficients(), polynomial(&[-2, 0, 1]).as_slice());
    assert_eq!(trace.final_root(), &bi(2166));
    assert_eq!(
        (trace.final_root() * trace.final_root() - bi(2)) % bi(2401),
        bi(0)
    );

    let x = trace.initial_root();
    assert_eq!((x * x - &a) % p.pow(1), bi(0));
    for step in trace.steps() {
        let x = step.root_after();
        let e = step.level() + 1;
        assert_eq!((x * x - &a) % p.pow(e), bi(0));
    }
}

#[test]
fn square_root_lift_rejects_singular_zero_root_on_the_simple_route() {
    assert_eq!(
        hensel_lift_square_root(&bi(0), &bi(0), &bu(5), 2),
        Err(HenselLiftError::SingularDerivativeModPrime)
    );
}

#[test]
fn fast_square_root_lift_doubles_the_precision_levels() {
    let a = bi(2);
    let p = bi(7);
    let trace = hensel_lift_square_root_fast(&a, &bi(3), &bu(7), 8)
        .expect("simple square root should lift by precision doubling");

    assert_eq!(trace.prime(), &bu(7));
    assert_eq!(trace.value(), &a);
    assert_eq!(trace.initial_root(), &bi(3));
    assert_eq!(trace.target_level(), 8);
    assert_eq!(trace.reached_level(), 8);
    assert_eq!(
        trace
            .steps()
            .iter()
            .map(|step| (step.source_level(), step.target_level()))
            .collect::<Vec<_>>(),
        vec![(1, 2), (2, 4), (4, 8)]
    );

    let x = trace.initial_root();
    assert_eq!((x * x - &a) % p.pow(1), bi(0));
    let mut previous_root = trace.initial_root();
    for step in trace.steps() {
        assert_eq!(step.root_before(), previous_root);
        let x = step.root_after();
        assert_eq!((x * x - &a) % p.pow(step.target_level()), bi(0));
        previous_root = x;
    }
}

#[test]
fn fast_square_root_lift_stops_at_a_non_power_of_two_target() {
    let a = bi(2);
    let p = bi(7);
    let trace = hensel_lift_square_root_fast(&a, &bi(3), &bu(7), 5)
        .expect("final fast step should stop at the requested precision");

    assert_eq!(trace.reached_level(), 5);
    assert_eq!(
        trace
            .steps()
            .iter()
            .map(|step| (step.source_level(), step.target_level()))
            .collect::<Vec<_>>(),
        vec![(1, 2), (2, 4), (4, 5)]
    );
    assert_eq!(
        (trace.final_root() * trace.final_root() - &a) % p.pow(5),
        bi(0)
    );
}

#[test]
fn fast_square_root_lift_matches_the_linear_square_root_lift() {
    let linear =
        hensel_lift_square_root(&bi(2), &bi(3), &bu(7), 8).expect("linear lift should reach p^8");
    let fast = hensel_lift_square_root_fast(&bi(2), &bi(3), &bu(7), 8)
        .expect("fast lift should reach p^8");

    assert_eq!(fast.final_root(), linear.final_root());
}

#[test]
fn fast_square_root_lift_rejects_singular_zero_root() {
    assert_eq!(
        hensel_lift_square_root_fast(&bi(0), &bi(0), &bu(5), 2),
        Err(HenselLiftError::SingularDerivativeModPrime)
    );
}

#[test]
fn target_level_one_validates_and_returns_the_normalized_initial_root() {
    let coefficients = polynomial(&[-2, 0, 1]);
    let trace = hensel_lift_simple_root(&coefficients, &bi(-4), &bu(7), 1)
        .expect("-4 is the same root as 3 modulo 7");

    assert_eq!(trace.initial_root(), &bi(3));
    assert_eq!(trace.final_root(), &bi(3));
    assert_eq!(trace.reached_level(), 1);
    assert!(trace.steps().is_empty());
}

#[test]
fn non_root_is_rejected_before_lifting() {
    let coefficients = polynomial(&[-2, 0, 1]);

    assert_eq!(
        hensel_lift_simple_root(&coefficients, &bi(2), &bu(7), 3),
        Err(HenselLiftError::RootDoesNotSolveCurrentModulus)
    );
}

#[test]
fn singular_derivative_is_rejected_for_the_simple_route() {
    let coefficients = polynomial(&[0, 0, 1]);

    assert_eq!(
        hensel_lift_simple_root(&coefficients, &bi(0), &bu(5), 2),
        Err(HenselLiftError::SingularDerivativeModPrime)
    );
}

#[test]
fn invalid_inputs_are_reported_explicitly() {
    assert_eq!(
        hensel_lift_simple_root(&[], &bi(0), &bu(5), 2),
        Err(HenselLiftError::EmptyPolynomial)
    );
    assert_eq!(
        hensel_lift_simple_root(&polynomial(&[1]), &bi(0), &bu(5), 2),
        Err(HenselLiftError::ConstantPolynomial)
    );
    assert_eq!(
        hensel_lift_simple_root(&polynomial(&[-2, 0, 1]), &bi(3), &bu(1), 2),
        Err(HenselLiftError::NonPrimeModulus)
    );
    assert_eq!(
        hensel_lift_simple_root(&polynomial(&[-2, 0, 1]), &bi(3), &bu(9), 2),
        Err(HenselLiftError::NonPrimeModulus)
    );
    assert_eq!(
        hensel_lift_simple_root(&polynomial(&[-2, 0, 1]), &bi(3), &bu(7), 0),
        Err(HenselLiftError::ZeroTargetLevel)
    );
}
