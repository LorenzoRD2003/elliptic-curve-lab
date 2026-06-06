use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::One;

/// Returns the exact Bernoulli number `B_n` using the Akiyama-Tanigawa
/// algorithm.
///
/// This implementation uses the standard Bernoulli-number convention
/// `B_0 = 1`, `B_1 = -1/2`, `B_2 = 1/6`, ...
///
/// the Akiyama-Tanigawa algorithm builds a triangular tableau of rational
/// numbers `A_j^(m)` by initializing each outer step with
/// `A_m^(m) = 1 / (m + 1)` and then updating backward through
///
/// `A_{j-1}^(m) = j * (A_{j-1}^(m) - A_j^(m))`
/// for `j = m, m-1, ..., 1`.
///
/// After the `m`-th outer step, the leftmost entry satisfies
/// `A_0^(m) = B_m^+`, where `B_1^+ = +1/2` is the “first Bernoulli number”
/// convention naturally produced by this tableau. We then translate only the
/// exceptional index `n = 1` to the more common convention
/// `B_1 = -1/2`, leaving every other `B_n` unchanged.
///
/// This is a compact exact algorithm well suited to educational and moderate
///-size symbolic computations where we want to expose the Bernoulli-number
/// definition honestly without introducing floating-point error.
///
/// Complexity:
/// - `Θ(n²)` exact rational arithmetic updates
/// - `Θ(n)` stored rational values
///
/// As usual for exact rational arithmetic, the bit-complexity also depends on
/// the size of the intermediate numerators and denominators, not only on the
/// outer index `n`.
pub fn bernoulli_number(n: usize) -> BigRational {
    let mut tableau = Vec::with_capacity(n + 1);

    for m in 0..=n {
        tableau.push(BigRational::new(BigInt::one(), BigInt::from(m + 1)));

        for j in (1..=m).rev() {
            let difference = tableau[j - 1].clone() - tableau[j].clone();
            tableau[j - 1] = BigRational::from_integer(BigInt::from(j)) * difference;
        }
    }

    if n == 1 {
        -tableau[0].clone()
    } else {
        tableau[0].clone()
    }
}

#[cfg(test)]
mod tests {
    use num_bigint::BigInt;
    use num_rational::BigRational;
    use num_traits::Zero;

    use crate::numerics::bernoulli_number;

    fn q(numerator: i64, denominator: i64) -> BigRational {
        BigRational::new(BigInt::from(numerator), BigInt::from(denominator))
    }

    #[test]
    fn first_classical_bernoulli_values_match_the_standard_convention() {
        assert_eq!(bernoulli_number(0), q(1, 1));
        assert_eq!(bernoulli_number(1), q(-1, 2));
        assert_eq!(bernoulli_number(2), q(1, 6));
        assert_eq!(bernoulli_number(4), q(-1, 30));
        assert_eq!(bernoulli_number(6), q(1, 42));
        assert_eq!(bernoulli_number(8), q(-1, 30));
        assert_eq!(bernoulli_number(10), q(5, 66));
    }

    #[test]
    fn odd_bernoulli_numbers_vanish_after_b1() {
        assert_eq!(bernoulli_number(3), BigRational::zero());
        assert_eq!(bernoulli_number(5), BigRational::zero());
        assert_eq!(bernoulli_number(7), BigRational::zero());
        assert_eq!(bernoulli_number(9), BigRational::zero());
    }

    #[test]
    fn output_is_already_normalized_as_an_exact_rational() {
        let value = bernoulli_number(12);

        assert_eq!(value, q(-691, 2730));
        assert_eq!(value.numer().clone(), BigInt::from(-691));
        assert_eq!(value.denom().clone(), BigInt::from(2730));
    }
}
