use num_bigint::{BigInt, BigUint, Sign};
use num_traits::{One, Signed, Zero};

/// Returns the greatest common divisor of two nonnegative integers.
///
/// Complexity: `Θ(log min(a, b))` exact remainder steps.
#[cfg_attr(not(test), allow(dead_code))]
pub(crate) fn gcd_biguint(left: &BigUint, right: &BigUint) -> BigUint {
    let mut left = left.clone();
    let mut right = right.clone();

    while !right.is_zero() {
        let remainder = &left % &right;
        left = right;
        right = remainder;
    }

    left
}

/// Computes the extended greatest common divisor of two integers.
///
/// The returned triple `(g, x, y)` satisfies `left * x + right * y = g`,
/// where `g` is the nonnegative gcd of `left` and `right`.
#[cfg_attr(not(test), allow(dead_code))]
pub(crate) fn extended_gcd_bigint(left: BigInt, right: BigInt) -> (BigInt, BigInt, BigInt) {
    if right.is_zero() {
        return (left.abs(), BigInt::from(1u8), BigInt::zero());
    }

    let (gcd, x1, y1) = extended_gcd_bigint(right.clone(), &left % &right);
    let x = y1.clone();
    let y = x1 - (left / right) * y1;

    (gcd, x, y)
}

/// Returns the inverse of `value` modulo `modulus`, when it exists.
///
/// This returns `None` exactly when `gcd(value, modulus) != 1`.
#[cfg_attr(not(test), allow(dead_code))]
pub(crate) fn inverse_mod_biguint(value: &BigUint, modulus: &BigUint) -> Option<BigUint> {
    let (gcd, inverse, _) =
        extended_gcd_bigint(BigInt::from(value.clone()), BigInt::from(modulus.clone()));
    if gcd != BigInt::one() {
        return None;
    }

    Some(positive_bigint_modulus(inverse, modulus))
}

fn positive_bigint_modulus(value: BigInt, modulus: &BigUint) -> BigUint {
    let modulus_bigint = BigInt::from(modulus.clone());
    let reduced = ((value % &modulus_bigint) + &modulus_bigint) % &modulus_bigint;
    match reduced.sign() {
        Sign::Minus => unreachable!("the normalized residue should be nonnegative"),
        _ => reduced
            .to_biguint()
            .expect("the normalized residue should convert back to BigUint"),
    }
}

#[cfg(test)]
mod tests {
    use super::{extended_gcd_bigint, gcd_biguint, inverse_mod_biguint};
    use num_bigint::{BigInt, BigUint};

    fn bu(value: u64) -> BigUint {
        BigUint::from(value)
    }

    #[test]
    fn gcd_biguint_handles_zero_and_nontrivial_inputs() {
        assert_eq!(gcd_biguint(&bu(0), &bu(0)), bu(0));
        assert_eq!(gcd_biguint(&bu(0), &bu(42)), bu(42));
        assert_eq!(gcd_biguint(&bu(42), &bu(0)), bu(42));
        assert_eq!(gcd_biguint(&bu(84), &bu(126)), bu(42));
    }

    #[test]
    fn extended_gcd_bigint_returns_bezout_coefficients() {
        let left = BigInt::from(84);
        let right = BigInt::from(126);
        let (gcd, x, y) = extended_gcd_bigint(left.clone(), right.clone());

        assert_eq!(gcd, BigInt::from(42));
        assert_eq!(left * x + right * y, gcd);
    }

    #[test]
    fn inverse_mod_biguint_recovers_a_unit_inverse() {
        let inverse = inverse_mod_biguint(&bu(3), &bu(11)).expect("3 should be invertible mod 11");

        assert_eq!(inverse, bu(4));
    }

    #[test]
    fn inverse_mod_biguint_returns_none_for_non_units() {
        assert_eq!(inverse_mod_biguint(&bu(6), &bu(9)), None);
    }
}
