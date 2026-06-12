use num_bigint::BigUint;
use num_prime::nt_funcs::{factorize, is_prime};
use num_traits::One;

/// Failure modes for normalized prime-power factorizations.
///
/// This helper stays in `numerics` because it describes only exact integer
/// structure, not any curve-specific semantics.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) enum PrimePowerFactorizationError {
    Empty,
    TrivialInteger,
    NonPositivePrimePowerBase,
    DuplicatePrime,
    ProductMismatch,
    CompositeBase,
}

/// Cached powers `1, ℓ, ℓ², ..., ℓ^e` for one fixed prime `ℓ`.
///
/// This small exact cache is shared infrastructure:
///
/// - consumers can reuse `ℓ^k` without recomputing powers ad hoc, and
/// - local algorithms can keep the `ℓ`-power ladder explicit and inspectable.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct PrimePowerTable {
    prime: BigUint,
    powers: Vec<BigUint>,
}

impl PrimePowerTable {
    pub(crate) fn up_through(prime: &BigUint, exponent_bound: u32) -> Self {
        let mut powers = Vec::with_capacity(exponent_bound as usize + 1);
        powers.push(BigUint::one());

        for exponent in 1..=exponent_bound {
            let next = &powers[exponent as usize - 1] * prime;
            powers.push(next);
        }

        Self {
            prime: prime.clone(),
            powers,
        }
    }

    pub(crate) fn prime(&self) -> &BigUint {
        &self.prime
    }

    pub(crate) fn exponent_bound(&self) -> u32 {
        self.powers
            .len()
            .checked_sub(1)
            .expect("power table should contain at least the zeroth power") as u32
    }

    pub(crate) fn power(&self, exponent: u32) -> &BigUint {
        &self.powers[exponent as usize]
    }
}

/// Returns `base^exponent` by exponentiation by squaring.
///
/// Complexity: `Θ(log exponent)` exact integer multiplications.
pub(crate) fn pow_biguint(base: &BigUint, exponent: u32) -> BigUint {
    let mut result = BigUint::one();
    let mut power = base.clone();
    let mut exponent_bits = exponent;

    while exponent_bits > 0 {
        if exponent_bits & 1 == 1 {
            result *= &power;
        }

        exponent_bits >>= 1;
        if exponent_bits > 0 {
            power = &power * &power;
        }
    }

    result
}

/// A structurally normalized prime-power factorization `M = Π ℓᵢ^eᵢ`.
///
/// The stored factors are sorted by increasing prime, have no repeated prime,
/// have positive exponents, and multiply back to the supplied `multiple`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct NormalizedPrimePowerFactorization {
    factors: Vec<(BigUint, u32)>,
}

impl NormalizedPrimePowerFactorization {
    /// Normalizes and also certifies that each base is prime.
    pub(crate) fn checked(
        multiple: &BigUint,
        factorization: &[(BigUint, u32)],
    ) -> Result<Self, PrimePowerFactorizationError> {
        let normalized = Self::trusted(multiple, factorization)?;

        for (prime, _) in normalized.as_slice() {
            if !is_prime(prime, None).probably() {
                return Err(PrimePowerFactorizationError::CompositeBase);
            }
        }

        Ok(normalized)
    }

    /// Factors one integer and returns its normalized prime-power decomposition.
    ///
    /// This is the route to prefer when callers already have only the integer
    /// `M` and want one canonical `Π ℓᵢ^eᵢ` representation instead of passing
    /// through a second ad hoc factorization helper.
    ///
    /// The input must satisfy `M >= 2`. Values `0` and `1` are rejected
    /// because they do not admit a meaningful non-empty prime-power
    /// factorization.
    pub(crate) fn factor(value: &BigUint) -> Result<Self, PrimePowerFactorizationError> {
        if value < &BigUint::from(2u8) {
            return Err(PrimePowerFactorizationError::TrivialInteger);
        }

        let mut factors = factorize(value.clone())
            .into_iter()
            .map(|(prime, exponent)| {
                (
                    prime,
                    u32::try_from(exponent).expect(
                        "num-prime exponents should fit into the normalized prime-power surface",
                    ),
                )
            })
            .collect::<Vec<_>>();
        factors.sort_unstable_by(|left, right| left.0.cmp(&right.0));

        Ok(Self { factors })
    }

    /// Normalizes a factorization whose prime bases are already trusted.
    ///
    /// This still checks the exact arithmetic structure:
    ///
    /// - non-empty input
    /// - bases at least `2`
    /// - positive exponents
    /// - no repeated prime labels after sorting
    /// - product equal to the supplied `multiple`
    pub(crate) fn trusted(
        multiple: &BigUint,
        factorization: &[(BigUint, u32)],
    ) -> Result<Self, PrimePowerFactorizationError> {
        if factorization.is_empty() {
            return Err(PrimePowerFactorizationError::Empty);
        }

        let mut factors = factorization.to_vec();
        factors.sort_unstable_by(|left, right| left.0.cmp(&right.0));

        for window in factors.windows(2) {
            if window[0].0 == window[1].0 {
                return Err(PrimePowerFactorizationError::DuplicatePrime);
            }
        }

        let mut product = BigUint::one();
        for (prime, exponent) in &factors {
            if *exponent == 0 || prime < &BigUint::from(2u8) {
                return Err(PrimePowerFactorizationError::NonPositivePrimePowerBase);
            }
            product *= pow_biguint(prime, *exponent);
        }

        if &product != multiple {
            return Err(PrimePowerFactorizationError::ProductMismatch);
        }

        Ok(Self { factors })
    }

    pub(crate) fn as_slice(&self) -> &[(BigUint, u32)] {
        &self.factors
    }

    pub(crate) fn into_factors(self) -> Vec<(BigUint, u32)> {
        self.factors
    }
}

#[cfg(test)]
mod tests {
    use super::{
        NormalizedPrimePowerFactorization, PrimePowerFactorizationError, PrimePowerTable,
        pow_biguint,
    };
    use num_bigint::BigUint;

    fn bu(value: u64) -> BigUint {
        BigUint::from(value)
    }

    #[test]
    fn pow_biguint_uses_exponentiation_by_squaring_correctly() {
        assert_eq!(pow_biguint(&bu(7), 0), bu(1));
        assert_eq!(pow_biguint(&bu(7), 1), bu(7));
        assert_eq!(pow_biguint(&bu(7), 5), bu(16807));
    }

    #[test]
    fn prime_power_table_lists_all_powers_up_through_the_bound() {
        let table = PrimePowerTable::up_through(&bu(3), 4);

        assert_eq!(table.prime(), &bu(3));
        assert_eq!(table.exponent_bound(), 4);
        assert_eq!(table.power(0), &bu(1));
        assert_eq!(table.power(1), &bu(3));
        assert_eq!(table.power(2), &bu(9));
        assert_eq!(table.power(3), &bu(27));
        assert_eq!(table.power(4), &bu(81));
    }

    #[test]
    fn checked_prime_power_factorization_sorts_and_certifies_primes() {
        let normalized =
            NormalizedPrimePowerFactorization::checked(&bu(72), &[(bu(3), 2), (bu(2), 3)])
                .expect("72 = 2^3 * 3^2");

        assert_eq!(normalized.as_slice(), &[(bu(2), 3), (bu(3), 2)]);
    }

    #[test]
    fn trusted_prime_power_factorization_skips_primality_but_keeps_structure_checks() {
        let normalized = NormalizedPrimePowerFactorization::trusted(&bu(36), &[(bu(6), 2)])
            .expect("trusted route should keep only structural checks");

        assert_eq!(normalized.as_slice(), &[(bu(6), 2)]);
    }

    #[test]
    fn checked_prime_power_factorization_rejects_composite_bases() {
        assert_eq!(
            NormalizedPrimePowerFactorization::checked(&bu(36), &[(bu(6), 2)]),
            Err(PrimePowerFactorizationError::CompositeBase)
        );
    }

    #[test]
    fn prime_power_factorization_rejects_repeated_primes() {
        assert_eq!(
            NormalizedPrimePowerFactorization::trusted(&bu(72), &[(bu(2), 1), (bu(2), 2)]),
            Err(PrimePowerFactorizationError::DuplicatePrime)
        );
    }

    #[test]
    fn prime_power_factorization_rejects_the_wrong_product() {
        assert_eq!(
            NormalizedPrimePowerFactorization::trusted(&bu(72), &[(bu(2), 3), (bu(3), 1)]),
            Err(PrimePowerFactorizationError::ProductMismatch)
        );
    }

    #[test]
    fn factor_builds_the_canonical_prime_power_decomposition() {
        let factorization =
            NormalizedPrimePowerFactorization::factor(&bu(72)).expect("72 should factor");

        assert_eq!(factorization.as_slice(), &[(bu(2), 3), (bu(3), 2)]);
    }

    #[test]
    fn factor_rejects_zero_and_one() {
        assert_eq!(
            NormalizedPrimePowerFactorization::factor(&bu(0)),
            Err(PrimePowerFactorizationError::TrivialInteger)
        );
        assert_eq!(
            NormalizedPrimePowerFactorization::factor(&bu(1)),
            Err(PrimePowerFactorizationError::TrivialInteger)
        );
    }
}
