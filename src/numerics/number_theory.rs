use core::fmt;
use num_bigint::{BigInt, BigUint, Sign};
use num_prime::nt_funcs::{factorize, is_prime};
use num_traits::{One, Zero};

/// Failure modes for helpers that expect a positive prime integer.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum PositivePrimeError {
    Zero,
    One,
    Composite,
}

impl fmt::Display for PositivePrimeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Zero => write!(
                f,
                "prime-dependent arithmetic requires a positive prime, not 0"
            ),
            Self::One => write!(
                f,
                "prime-dependent arithmetic requires a prime integer greater than 1"
            ),
            Self::Composite => write!(f, "prime-dependent arithmetic requires a prime integer"),
        }
    }
}

impl std::error::Error for PositivePrimeError {}

/// Returns whether an integer is squarefree.
///
/// This shared helper is purely arithmetic: it answers whether no prime square
/// divides `n`. Equivalently, in the prime factorization `n = ± Π p_i^{e_i}`,
/// the value is squarefree exactly when every exponent satisfies `e_i ≤ 1`.
///
/// By convention:
/// - `0` is not squarefree
/// - `1` and `-1` are squarefree
///
/// The current implementation factors the absolute value with `num-prime` and
/// checks whether every prime exponent is at most `1`.
///
/// Complexity: dominated by `num-prime`'s implementation.
pub fn is_squarefree(value: &BigInt) -> bool {
    match value.sign() {
        Sign::NoSign => false,
        _ => {
            let abs_value: BigUint = value.magnitude().clone();
            if abs_value == BigUint::from(1u8) {
                return true;
            }
            factorize(abs_value).values().all(|&exponent| exponent <= 1)
        }
    }
}

/// Returns the positive divisors of `n`, sorted increasingly.
///
/// The output is the finite set
/// `{ d in Z_{>0} : d | n }`, represented as a sorted vector of `BigUint`.
///
/// By convention, this helper returns the empty vector for `n = 0`, since the
/// mathematical set of positive divisors of `0` is infinite and cannot be
/// represented by a finite `Vec`.
///
/// The current implementation factors `n` with `num-prime`, enumerates all
/// divisor combinations from the prime-power data, and sorts the result.
///
/// Complexity: dominated by `num-prime`, plus `Θ(τ(n) log τ(n))` big-integer
/// comparisons for the final sort, where `τ(n)` is the number of positive
/// divisors of `n`.
pub fn positive_divisors(n: &BigUint) -> Vec<BigUint> {
    if n.is_zero() {
        return Vec::new();
    }

    let factors = factorize(n.clone());
    let mut divisors = vec![BigUint::one()];

    for (prime, exponent) in factors {
        let previous_divisors = divisors.clone();
        let mut prime_power = BigUint::one();

        for _ in 0..exponent {
            prime_power *= &prime;
            for divisor in &previous_divisors {
                divisors.push(divisor * &prime_power);
            }
        }
    }

    divisors.sort();
    divisors
}

/// Returns `⌊√n⌋` for a non-negative integer `n`.
///
/// Complexity: `Θ(log n)` exact big-integer comparisons and multiplications
/// using binary search over `[0, n + 1)`.
pub(crate) fn integer_sqrt(n: &BigUint) -> BigUint {
    if n.is_zero() || n.is_one() {
        return n.clone();
    }

    let mut low = BigUint::zero();
    let mut high = n + BigUint::one();
    while &high - &low > BigUint::one() {
        let middle = (&low + &high) >> 1usize;
        if &middle * &middle <= *n {
            low = middle;
        } else {
            high = middle;
        }
    }

    low
}

/// Returns `√n` exactly when `n` is a square in `ℤ`.
///
/// Complexity: one [`integer_sqrt`] call plus one exact square check.
pub(crate) fn exact_square_root(n: &BigUint) -> Option<BigUint> {
    let root = integer_sqrt(n);
    if &root * &root == *n {
        Some(root)
    } else {
        None
    }
}

/// Returns the `ℓ`-adic valuation `v_ℓ(n)` of a positive integer `n`.
///
/// This is the largest exponent `a >= 0` such that `ℓ^a | n`.
///
/// The input `prime` must be a genuine prime integer. For example,
/// `v_2(48) = 4` and `v_3(48) = 1`.
///
/// By convention, this helper returns `0` when `n = 0`, since the project uses
/// it only for positive conductor data where `0` should not arise and this
/// choice keeps the helper total on `BigUint`.
///
/// Complexity: prime validation is dominated by `num-prime`. After validation,
/// the implementation performs `Θ(v_ℓ(n))` exact big-integer divisions.
pub(crate) fn valuation_biguint(n: &BigUint, prime: &BigUint) -> Result<u32, PositivePrimeError> {
    validate_positive_prime(prime)?;

    if n.is_zero() {
        return Ok(0);
    }

    let mut valuation = 0u32;
    let mut residual = n.clone();

    while (&residual % prime).is_zero() {
        residual /= prime;
        valuation += 1;
    }

    Ok(valuation)
}

fn validate_positive_prime(prime: &BigUint) -> Result<(), PositivePrimeError> {
    if prime.is_zero() {
        return Err(PositivePrimeError::Zero);
    }
    if prime == &BigUint::one() {
        return Err(PositivePrimeError::One);
    }
    if !is_prime(prime, None).probably() {
        return Err(PositivePrimeError::Composite);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use crate::numerics::{
        PositivePrimeError, exact_square_root, integer_sqrt, is_squarefree, positive_divisors,
        valuation_biguint,
    };
    use num_bigint::{BigInt, BigUint};

    #[test]
    fn zero_is_not_squarefree() {
        assert!(!is_squarefree(&BigInt::from(0)));
    }

    #[test]
    fn units_are_squarefree() {
        assert!(is_squarefree(&BigInt::from(1)));
        assert!(is_squarefree(&BigInt::from(-1)));
    }

    #[test]
    fn repeated_prime_factors_are_detected() {
        assert!(!is_squarefree(&BigInt::from(12)));
        assert!(!is_squarefree(&BigInt::from(-48)));
    }

    #[test]
    fn squarefree_examples_pass() {
        assert!(is_squarefree(&BigInt::from(-163)));
        assert!(is_squarefree(&BigInt::from(30)));
    }

    #[test]
    fn positive_divisors_of_zero_are_left_empty_by_convention() {
        assert!(positive_divisors(&BigUint::from(0u8)).is_empty());
    }

    #[test]
    fn positive_divisors_of_one_are_just_one() {
        assert_eq!(
            positive_divisors(&BigUint::from(1u8)),
            vec![BigUint::from(1u8)]
        );
    }

    #[test]
    fn positive_divisors_are_sorted_and_complete() {
        assert_eq!(
            positive_divisors(&BigUint::from(12u8)),
            vec![
                BigUint::from(1u8),
                BigUint::from(2u8),
                BigUint::from(3u8),
                BigUint::from(4u8),
                BigUint::from(6u8),
                BigUint::from(12u8),
            ]
        );
    }

    #[test]
    fn integer_sqrt_returns_the_floor_root() {
        assert_eq!(integer_sqrt(&BigUint::from(0u8)), BigUint::from(0u8));
        assert_eq!(integer_sqrt(&BigUint::from(1u8)), BigUint::from(1u8));
        assert_eq!(integer_sqrt(&BigUint::from(15u8)), BigUint::from(3u8));
        assert_eq!(integer_sqrt(&BigUint::from(16u8)), BigUint::from(4u8));
        assert_eq!(integer_sqrt(&BigUint::from(17u8)), BigUint::from(4u8));
    }

    #[test]
    fn exact_square_root_distinguishes_squares_from_non_squares() {
        assert_eq!(
            exact_square_root(&BigUint::from(144u16)),
            Some(BigUint::from(12u8))
        );
        assert_eq!(exact_square_root(&BigUint::from(145u16)), None);
    }

    #[test]
    fn valuation_counts_the_prime_exponent() {
        assert_eq!(
            valuation_biguint(&BigUint::from(48u8), &BigUint::from(2u8)),
            Ok(4)
        );
        assert_eq!(
            valuation_biguint(&BigUint::from(48u8), &BigUint::from(3u8)),
            Ok(1)
        );
        assert_eq!(
            valuation_biguint(&BigUint::from(48u8), &BigUint::from(5u8)),
            Ok(0)
        );
    }

    #[test]
    fn valuation_rejects_non_prime_inputs() {
        assert_eq!(
            valuation_biguint(&BigUint::from(48u8), &BigUint::from(0u8)),
            Err(PositivePrimeError::Zero)
        );
        assert_eq!(
            valuation_biguint(&BigUint::from(48u8), &BigUint::from(1u8)),
            Err(PositivePrimeError::One)
        );
        assert_eq!(
            valuation_biguint(&BigUint::from(48u8), &BigUint::from(4u8)),
            Err(PositivePrimeError::Composite)
        );
    }

    #[test]
    fn positive_prime_error_display_messages_remain_specific() {
        assert_eq!(
            PositivePrimeError::Zero.to_string(),
            "prime-dependent arithmetic requires a positive prime, not 0"
        );
        assert_eq!(
            PositivePrimeError::One.to_string(),
            "prime-dependent arithmetic requires a prime integer greater than 1"
        );
        assert_eq!(
            PositivePrimeError::Composite.to_string(),
            "prime-dependent arithmetic requires a prime integer"
        );
    }
}
