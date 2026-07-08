use num_bigint::{BigInt, BigUint};
use num_traits::Zero;

/// Returns the Euclidean floor quotient `⌊a / b⌋` for one positive divisor `b`.
///
/// This differs from Rust's truncating integer division when `a < 0`.
///
/// Complexity: `Θ(1)` BigInt division plus one remainder inspection.
pub(crate) fn floor_div_bigint_by_positive(value: &BigInt, positive_divisor: &BigInt) -> BigInt {
    let quotient = value / positive_divisor;
    let remainder = value % positive_divisor;
    if remainder.sign() == num_bigint::Sign::Minus {
        quotient - BigInt::from(1u8)
    } else {
        quotient
    }
}

/// Returns the Euclidean ceiling quotient `⌈a / b⌉` for one positive divisor `b`.
///
/// This differs from Rust's truncating integer division when `a > 0` and
/// `a` is not divisible by `b`.
///
/// Complexity: `Θ(1)` BigInt division plus one remainder inspection.
pub(crate) fn ceil_div_bigint_by_positive(value: &BigInt, positive_divisor: &BigInt) -> BigInt {
    let quotient = value / positive_divisor;
    let remainder = value % positive_divisor;
    if remainder.sign() == num_bigint::Sign::Plus {
        quotient + BigInt::from(1u8)
    } else {
        quotient
    }
}

/// Returns the least non-negative representative of `value` modulo `modulus`.
///
/// This is the Euclidean residue in `[0, modulus)`, unlike Rust's `%` on
/// signed integers, whose remainder keeps the sign of the dividend.
///
/// Complexity: `Θ(1)` BigInt division plus normalization.
pub(crate) fn positive_mod_bigint(value: &BigInt, modulus: &BigInt) -> BigInt {
    debug_assert!(modulus > &BigInt::zero());
    ((value % modulus) + modulus) % modulus
}

/// Returns the least non-negative representative of `value` modulo `modulus`.
///
/// This is the same Euclidean residue as [`positive_mod_bigint`], specialized
/// for a positive unsigned modulus and returned as a [`BigUint`].
///
/// Complexity: `Θ(1)` BigInt division plus normalization.
pub(crate) fn positive_mod_biguint(value: &BigInt, modulus: &BigUint) -> BigUint {
    positive_mod_bigint(value, &BigInt::from(modulus.clone()))
        .to_biguint()
        .expect("positive modular representative should convert to BigUint")
}

#[cfg(test)]
mod tests {

    use super::{
        ceil_div_bigint_by_positive, floor_div_bigint_by_positive, positive_mod_bigint,
        positive_mod_biguint,
    };
    use num_bigint::{BigInt, BigUint};

    #[test]
    fn floor_division_uses_the_euclidean_convention() {
        let divisor = BigInt::from(5u8);

        assert_eq!(
            floor_div_bigint_by_positive(&BigInt::from(13i8), &divisor),
            BigInt::from(2i8)
        );
        assert_eq!(
            floor_div_bigint_by_positive(&BigInt::from(-13i8), &divisor),
            BigInt::from(-3i8)
        );
    }

    #[test]
    fn ceil_division_uses_the_euclidean_convention() {
        let divisor = BigInt::from(5u8);

        assert_eq!(
            ceil_div_bigint_by_positive(&BigInt::from(13i8), &divisor),
            BigInt::from(3i8)
        );
        assert_eq!(
            ceil_div_bigint_by_positive(&BigInt::from(-13i8), &divisor),
            BigInt::from(-2i8)
        );
    }

    #[test]
    fn positive_mod_bigint_uses_euclidean_representatives() {
        let modulus = BigInt::from(5u8);

        assert_eq!(
            positive_mod_bigint(&BigInt::from(13i8), &modulus),
            BigInt::from(3u8)
        );
        assert_eq!(
            positive_mod_bigint(&BigInt::from(-13i8), &modulus),
            BigInt::from(2u8)
        );
    }

    #[test]
    fn positive_mod_biguint_returns_unsigned_representatives() {
        assert_eq!(
            positive_mod_biguint(&BigInt::from(-13i8), &BigUint::from(5u8)),
            BigUint::from(2u8)
        );
    }
}
