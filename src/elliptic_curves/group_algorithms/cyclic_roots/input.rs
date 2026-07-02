use num_bigint::BigUint;
use num_prime::nt_funcs::is_prime;
use num_traits::{One, Zero};

/// Validation errors for prime-degree cyclic-root input metadata.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) enum CyclicPrimeRootInputError {
    /// The group order `|G|` must be positive.
    ZeroGroupOrder,
    /// The requested root degree `r` must be prime.
    NonPrimeRootDegree { root_degree: BigUint },
    /// The exponent `k` in `|G| = a r^k` did not fit the staged report type.
    SylowExponentTooLarge { root_prime: BigUint },
}

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
    pub(crate) fn from_group_order_and_prime(
        group_order: BigUint,
        root_prime: BigUint,
    ) -> Result<Self, CyclicPrimeRootInputError> {
        if group_order.is_zero() {
            return Err(CyclicPrimeRootInputError::ZeroGroupOrder);
        }
        if !is_prime(&root_prime, None).probably() {
            return Err(CyclicPrimeRootInputError::NonPrimeRootDegree {
                root_degree: root_prime,
            });
        }

        let mut prime_to_root_cofactor = group_order.clone();
        let mut sylow_order = BigUint::one();
        let mut sylow_exponent = 0u32;

        while (&prime_to_root_cofactor % &root_prime).is_zero() {
            prime_to_root_cofactor /= &root_prime;
            sylow_order *= &root_prime;
            sylow_exponent = sylow_exponent.checked_add(1).ok_or_else(|| {
                CyclicPrimeRootInputError::SylowExponentTooLarge {
                    root_prime: root_prime.clone(),
                }
            })?;
        }

        Ok(Self {
            group_order,
            root_prime,
            prime_to_root_cofactor,
            sylow_order,
            sylow_exponent,
        })
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
