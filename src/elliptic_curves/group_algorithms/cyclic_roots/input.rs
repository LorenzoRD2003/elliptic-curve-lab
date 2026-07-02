use num_bigint::BigUint;

/// Input metadata for prime-degree root extraction in a finite cyclic group.
///
/// The cyclic group is written additively. Given a prime `r`, the algorithm
/// writes `|G| = a r^k`, with `gcd(a, r) = 1`, and receives an element `δ ∈ G`
/// of order `r^k`, generating the `r`-Sylow subgroup. This value object records
/// only the integer side of that setup.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct CyclicPrimeRootInput {
    group_order: BigUint,
    root_prime: BigUint,
    prime_to_root_cofactor: BigUint,
    sylow_order: BigUint,
    sylow_exponent: u32,
}

impl CyclicPrimeRootInput {
    pub(crate) fn new(
        group_order: BigUint,
        root_prime: BigUint,
        prime_to_root_cofactor: BigUint,
        sylow_order: BigUint,
        sylow_exponent: u32,
    ) -> Self {
        Self {
            group_order,
            root_prime,
            prime_to_root_cofactor,
            sylow_order,
            sylow_exponent,
        }
    }

    /// Returns the cyclic group order `|G|`.
    pub(crate) fn group_order(&self) -> &BigUint {
        &self.group_order
    }

    /// Returns the requested prime root degree `r`.
    pub(crate) fn root_prime(&self) -> &BigUint {
        &self.root_prime
    }

    /// Returns the factor `a` in `|G| = a r^k`, where `gcd(a, r) = 1`.
    pub(crate) fn prime_to_root_cofactor(&self) -> &BigUint {
        &self.prime_to_root_cofactor
    }

    /// Returns the `r`-Sylow order `r^k`.
    pub(crate) fn sylow_order(&self) -> &BigUint {
        &self.sylow_order
    }

    /// Returns the exponent `k` in the `r`-Sylow order `r^k`.
    pub(crate) fn sylow_exponent(&self) -> u32 {
        self.sylow_exponent
    }

    /// Returns whether the requested prime `r` divides `|G|`.
    pub(crate) fn root_prime_divides_group_order(&self) -> bool {
        self.sylow_exponent > 0
    }
}
