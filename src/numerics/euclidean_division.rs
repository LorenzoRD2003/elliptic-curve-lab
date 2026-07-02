use num_bigint::BigInt;

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

#[cfg(test)]
mod tests {

    use super::{ceil_div_bigint_by_positive, floor_div_bigint_by_positive};
    use num_bigint::BigInt;

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
}
