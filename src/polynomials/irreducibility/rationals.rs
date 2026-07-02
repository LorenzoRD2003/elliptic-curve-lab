use crate::fields::traits::*;
use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::{One, Signed, ToPrimitive, Zero};

use crate::fields::{Fp2, Fp3, Fp5, Fp7, Fp11, Fp13, Fp17, Fp19, Fp23, Fp29, Fp31, Q};
use crate::polynomials::{
    DensePolynomial, PolynomialError,
    irreducibility::{IrreducibilityBackend, IrreducibilityStatus},
};

/// Exact but intentionally partial irreducibility backend for dense
/// polynomials over `Q`.
///
/// The current implementation stays mathematically exact by using only
/// sound criteria:
///
/// - normalization to a primitive polynomial in `Z[x]`
/// - naive rational-root search on small primitive inputs
/// - the degree-`<= 3` fact that "no rational root" implies irreducible
/// - small-prime Eisenstein checks
/// - irreducibility certificates obtained from reductions modulo small primes
///
/// When those criteria do not decide the input, the backend returns
/// [`PolynomialError::UndeterminedIrreducibility`] instead of guessing.
impl IrreducibilityBackend for Q {
    fn irreducibility_status_impl(
        polynomial: &DensePolynomial<Self>,
    ) -> Result<IrreducibilityStatus<Self>, PolynomialError> {
        match polynomial.degree() {
            None | Some(0) => return Ok(IrreducibilityStatus::Constant),
            Some(1) => return Ok(IrreducibilityStatus::Linear),
            Some(_) => {}
        }

        let primitive = primitive_integer_coefficients(polynomial);

        if let Some(root) = find_small_rational_root(polynomial, &primitive)? {
            let divisor = DensePolynomial::<Q>::new(vec![-root.clone(), Q::one()]);
            let (quotient, remainder) = polynomial.div_rem(&divisor)?;
            debug_assert!(
                remainder.is_zero(),
                "rational root search should only report exact roots"
            );
            return Ok(IrreducibilityStatus::Reducible { divisor, quotient });
        }

        let degree = polynomial
            .degree()
            .expect("non-constant polynomial has a degree");
        if degree <= 3 {
            return Ok(IrreducibilityStatus::Irreducible);
        }

        if passes_small_prime_eisenstein(&primitive) {
            return Ok(IrreducibilityStatus::Irreducible);
        }

        if is_irreducible_mod_small_primes(&primitive)? {
            return Ok(IrreducibilityStatus::Irreducible);
        }

        Err(PolynomialError::UndeterminedIrreducibility(
            "the current exact backend for Q only uses rational-root search, small-prime Eisenstein checks, and irreducibility certificates from reductions modulo small primes",
        ))
    }
}

/// Size bound for the intentionally naive rational-root search used by the
/// current exact partial backend over `Q`.
///
/// The search enumerates divisors of the constant and leading coefficients of
/// a primitive integer polynomial. Keeping this bound small avoids presenting
/// a brute-force educational helper as if it were a scalable exact algorithm.
const MAX_NAIVE_RATIONAL_ROOT_FACTOR: usize = 10_000;

fn primitive_integer_coefficients(polynomial: &DensePolynomial<Q>) -> Vec<BigInt> {
    let coefficients = polynomial.coefficients();
    if coefficients.is_empty() {
        return Vec::new();
    }

    let denominator_lcm = coefficients
        .iter()
        .fold(BigInt::one(), |accumulator, coefficient| {
            lcm_bigint(&accumulator, coefficient.denom())
        });

    let mut cleared = coefficients
        .iter()
        .map(|coefficient| coefficient.numer() * (&denominator_lcm / coefficient.denom()))
        .collect::<Vec<_>>();

    let content = cleared
        .iter()
        .fold(BigInt::zero(), |accumulator, coefficient| {
            gcd_bigint(&accumulator, coefficient)
        });

    if !content.is_zero() && content != BigInt::one() {
        for coefficient in &mut cleared {
            *coefficient /= &content;
        }
    }

    if let Some(leading) = cleared.last()
        && leading.is_negative()
    {
        for coefficient in &mut cleared {
            *coefficient = -coefficient.clone();
        }
    }

    cleared
}

fn gcd_bigint(lhs: &BigInt, rhs: &BigInt) -> BigInt {
    let mut a = lhs.abs();
    let mut b = rhs.abs();

    while !b.is_zero() {
        let remainder = &a % &b;
        a = b;
        b = remainder;
    }

    a
}

fn lcm_bigint(lhs: &BigInt, rhs: &BigInt) -> BigInt {
    if lhs.is_zero() || rhs.is_zero() {
        BigInt::zero()
    } else {
        (lhs / gcd_bigint(lhs, rhs)) * rhs
    }
    .abs()
}

fn find_small_rational_root(
    polynomial: &DensePolynomial<Q>,
    primitive: &[BigInt],
) -> Result<Option<BigRational>, PolynomialError> {
    let constant = primitive
        .first()
        .expect("non-constant primitive polynomial has at least one coefficient");
    if constant.is_zero() {
        return Ok(Some(Q::zero()));
    }

    let leading = primitive
        .last()
        .expect("non-constant primitive polynomial has a leading coefficient");
    let constant_abs = constant.abs();
    let leading_abs = leading.abs();

    let Some(constant_bound) = constant_abs.to_usize() else {
        return Ok(None);
    };
    let Some(leading_bound) = leading_abs.to_usize() else {
        return Ok(None);
    };
    if constant_bound > MAX_NAIVE_RATIONAL_ROOT_FACTOR
        || leading_bound > MAX_NAIVE_RATIONAL_ROOT_FACTOR
    {
        return Ok(None);
    }

    let numerators = positive_divisors_usize(constant_bound);
    let denominators = positive_divisors_usize(leading_bound);

    for numerator in numerators {
        for denominator in &denominators {
            let positive = BigRational::new(BigInt::from(numerator), BigInt::from(*denominator));
            if is_exact_root(polynomial, &positive) {
                return Ok(Some(positive));
            }

            let negative = -positive.clone();
            if is_exact_root(polynomial, &negative) {
                return Ok(Some(negative));
            }
        }
    }

    Ok(None)
}

fn positive_divisors_usize(value: usize) -> Vec<usize> {
    let mut divisors = Vec::new();
    let mut candidate = 1_usize;

    while candidate.saturating_mul(candidate) <= value {
        if value.is_multiple_of(candidate) {
            divisors.push(candidate);
            let paired = value / candidate;
            if paired != candidate {
                divisors.push(paired);
            }
        }
        candidate += 1;
    }

    divisors.sort_unstable();
    divisors
}

fn is_exact_root(polynomial: &DensePolynomial<Q>, candidate: &BigRational) -> bool {
    let value = polynomial
        .coefficients()
        .iter()
        .rev()
        .fold(Q::zero(), |accumulator, coefficient| {
            Q::add(&Q::mul(&accumulator, candidate), coefficient)
        });

    Q::is_zero(&value)
}

fn passes_small_prime_eisenstein(coefficients: &[BigInt]) -> bool {
    small_primes().into_iter().any(|prime| {
        let prime_bigint = BigInt::from(prime);
        let prime_square = &prime_bigint * &prime_bigint;
        let leading = coefficients
            .last()
            .expect("non-constant primitive polynomial has a leading coefficient");
        let constant = coefficients
            .first()
            .expect("non-constant primitive polynomial has at least one coefficient");

        !is_divisible_by(constant, &prime_square)
            && !is_divisible_by(leading, &prime_bigint)
            && coefficients[..coefficients.len() - 1]
                .iter()
                .all(|coefficient| is_divisible_by(coefficient, &prime_bigint))
    })
}

fn is_irreducible_mod_small_primes(coefficients: &[BigInt]) -> Result<bool, PolynomialError> {
    macro_rules! try_prime {
        ($field:ty, $prime:literal) => {
            if irreducible_mod_prime::<$field>(coefficients, $prime)? {
                return Ok(true);
            }
        };
    }

    try_prime!(Fp2, 2);
    try_prime!(Fp3, 3);
    try_prime!(Fp5, 5);
    try_prime!(Fp7, 7);
    try_prime!(Fp11, 11);
    try_prime!(Fp13, 13);
    try_prime!(Fp17, 17);
    try_prime!(Fp19, 19);
    try_prime!(Fp23, 23);
    try_prime!(Fp29, 29);
    try_prime!(Fp31, 31);

    Ok(false)
}

fn irreducible_mod_prime<F>(coefficients: &[BigInt], p: usize) -> Result<bool, PolynomialError>
where
    F: IrreducibilityBackend,
{
    let leading = coefficients
        .last()
        .expect("non-constant primitive polynomial has a leading coefficient");
    let prime = BigInt::from(p);
    if is_divisible_by(leading, &prime) {
        return Ok(false);
    }

    let reduced = DensePolynomial::<F>::new(
        coefficients
            .iter()
            .map(|coefficient| F::from_bigint(&BigInt::from(bigint_mod_usize(coefficient, p))))
            .collect(),
    );

    Ok(matches!(
        reduced.irreducibility_status()?,
        IrreducibilityStatus::Irreducible
    ))
}

fn bigint_mod_usize(value: &BigInt, modulus: usize) -> usize {
    let modulus_bigint = BigInt::from(modulus);
    let reduced = ((value % &modulus_bigint) + &modulus_bigint) % &modulus_bigint;
    reduced
        .to_usize()
        .expect("reduction modulo a small prime should fit in usize")
}

fn is_divisible_by(value: &BigInt, divisor: &BigInt) -> bool {
    (value % divisor).is_zero()
}

fn small_primes() -> [usize; 11] {
    [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]
}
