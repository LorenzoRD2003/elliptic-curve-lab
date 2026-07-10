use crate::elliptic_curves::endomorphisms::{
    quadratic_ideals::{PrimeNormIdeal, PrimeNormIdealError},
    quadratic_orders::{ImaginaryQuadraticOrder, QuadraticDiscriminant},
};

use super::{bu, maximal_order};

#[test]
fn prime_norm_ideal_rejects_root_outside_split_pair() {
    let order = maximal_order(-23);

    let error = PrimeNormIdeal::split(order, bu(3), bu(0))
        .expect_err("0 is not a split root of -23 modulo 3");

    assert_eq!(error, PrimeNormIdealError::RootDoesNotMatchPrimeBehavior);
}

#[test]
fn prime_norm_ideal_rejects_inert_primes() {
    let order = maximal_order(-23);

    let error =
        PrimeNormIdeal::split(order, bu(5), bu(1)).expect_err("5 is inert for discriminant -23");

    assert_eq!(error, PrimeNormIdealError::NonSplitPrime);
}

#[test]
fn ramified_prime_ideal_rejects_split_primes() {
    let order = maximal_order(-23);

    let error =
        PrimeNormIdeal::ramified(order, bu(3)).expect_err("3 splits, so it is not ramified");

    assert_eq!(error, PrimeNormIdealError::PrimeDoesNotRamify);
}

#[test]
fn ramified_prime_ideal_rejects_inert_primes() {
    let order = maximal_order(-23);

    let error =
        PrimeNormIdeal::ramified(order, bu(5)).expect_err("5 is inert, so it is not ramified");

    assert_eq!(error, PrimeNormIdealError::PrimeDoesNotRamify);
}

#[test]
fn prime_norm_ideal_rejects_conductor_dividing_primes() {
    let order = ImaginaryQuadraticOrder::new(QuadraticDiscriminant::new(-23), bu(3))
        .expect("conductor 3 should define a non-maximal imaginary quadratic order");

    let error = PrimeNormIdeal::split(order, bu(3), bu(1))
        .expect_err("3 divides the conductor and is not invertible");

    assert_eq!(
        error,
        PrimeNormIdealError::NonInvertibleBecauseDividesConductor
    );
}

#[test]
fn ramified_prime_ideal_rejects_conductor_dividing_primes() {
    let order = ImaginaryQuadraticOrder::new(QuadraticDiscriminant::new(-23), bu(23))
        .expect("conductor 23 should define a non-maximal imaginary quadratic order");

    let error = PrimeNormIdeal::ramified(order, bu(23))
        .expect_err("23 divides the conductor and is not invertible");

    assert_eq!(
        error,
        PrimeNormIdealError::NonInvertibleBecauseDividesConductor
    );
}
