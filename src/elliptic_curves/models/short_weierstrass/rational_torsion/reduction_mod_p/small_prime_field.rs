use core::fmt;

use num_bigint::{BigInt, BigUint, Sign};
use num_prime::nt_funcs::is_prime;
use num_traits::Zero;

/// A small runtime prime used only by the reduction-mod-`p` torsion route.
///
/// This is intentionally not a field backend. It keeps just enough arithmetic
/// to enumerate an integral short-Weierstrass model modulo a good prime chosen
/// while the rational-torsion algorithm is running.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) struct ReductionPrime {
    p: u32,
}

/// A canonical residue modulo a [`ReductionPrime`].
///
/// The representative always lies in `0 ≤ value < p`; arithmetic keeps the
/// modulus as explicit context on [`ReductionPrime`] so this type does not
/// become a standalone ambient field element.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(super) struct ReductionResidue {
    value: u32,
}

/// Validation failures for a runtime reduction prime.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum ReductionPrimeError {
    /// The requested modulus was `0`.
    Zero,
    /// The requested modulus was `1`.
    One,
    /// The requested modulus was composite.
    Composite { modulus: u32 },
}

impl ReductionPrime {
    /// Validates and stores a runtime prime modulus.
    ///
    /// Complexity: dominated by `num-prime` primality testing on a `u32`
    /// integer.
    pub(super) fn new(p: u32) -> Result<Self, ReductionPrimeError> {
        match p {
            0 => Err(ReductionPrimeError::Zero),
            1 => Err(ReductionPrimeError::One),
            _ if !is_prime(&BigUint::from(p), None).probably() => {
                Err(ReductionPrimeError::Composite { modulus: p })
            }
            _ => Ok(Self { p }),
        }
    }

    pub(super) fn modulus(&self) -> u32 {
        self.p
    }

    /// Reduces an integer to its canonical residue modulo `p`.
    ///
    /// Complexity: `Θ(log |n|)` big-integer division work.
    pub(super) fn reduce_bigint(&self, n: &BigInt) -> ReductionResidue {
        let modulus = BigUint::from(self.p);
        let value = match n.sign() {
            Sign::Minus => {
                let magnitude = n.magnitude() % &modulus;
                if magnitude.is_zero() {
                    BigUint::zero()
                } else {
                    &modulus - magnitude
                }
            }
            _ => n.magnitude() % &modulus,
        };

        ReductionResidue {
            value: value
                .try_into()
                .expect("reduction modulo a u32 prime fits in u32"),
        }
    }

    /// Embeds a non-negative machine integer by reducing it modulo `p`.
    ///
    /// Complexity: `Θ(1)`.
    pub(super) fn reduce_u64(&self, n: u64) -> ReductionResidue {
        ReductionResidue {
            value: (n % u64::from(self.p)) as u32,
        }
    }

    /// Returns every residue in `ℤ/pℤ` in canonical order.
    ///
    /// Complexity: `Θ(p)`.
    pub(super) fn residues(&self) -> impl Iterator<Item = ReductionResidue> {
        (0..self.p).map(|value| ReductionResidue { value })
    }

    /// Adds two residues modulo `p`.
    ///
    /// Complexity: `Θ(1)`.
    pub(super) fn add(&self, left: ReductionResidue, right: ReductionResidue) -> ReductionResidue {
        self.reduce_u64(u64::from(left.value) + u64::from(right.value))
    }

    /// Subtracts two residues modulo `p`.
    ///
    /// Complexity: `Θ(1)`.
    #[cfg(test)]
    pub(super) fn sub(&self, left: ReductionResidue, right: ReductionResidue) -> ReductionResidue {
        let p = u64::from(self.p);
        self.reduce_u64(p + u64::from(left.value) - u64::from(right.value))
    }

    /// Multiplies two residues modulo `p`.
    ///
    /// Complexity: `Θ(1)`.
    pub(super) fn mul(&self, left: ReductionResidue, right: ReductionResidue) -> ReductionResidue {
        self.reduce_u64(u64::from(left.value) * u64::from(right.value))
    }

    /// Negates a residue modulo `p`.
    ///
    /// Complexity: `Θ(1)`.
    #[cfg(test)]
    pub(super) fn neg(&self, value: ReductionResidue) -> ReductionResidue {
        if value.is_zero() {
            value
        } else {
            ReductionResidue {
                value: self.p - value.value,
            }
        }
    }

    /// Raises a residue to a non-negative exponent by repeated squaring.
    ///
    /// Complexity: `Θ(log exponent)` modular multiplications.
    #[cfg(test)]
    pub(super) fn pow(&self, base: ReductionResidue, exponent: u32) -> ReductionResidue {
        let mut remaining = exponent;
        let mut power = base;
        let mut accumulator = self.one();

        while remaining > 0 {
            if remaining % 2 == 1 {
                accumulator = self.mul(accumulator, power);
            }
            power = self.mul(power, power);
            remaining /= 2;
        }

        accumulator
    }

    /// Returns `value⁻¹ mod p`, if `value ≠ 0`.
    ///
    /// Complexity: `Θ(log p)` modular multiplications via Fermat's theorem.
    #[cfg(test)]
    pub(super) fn inv(&self, value: ReductionResidue) -> Option<ReductionResidue> {
        if value.is_zero() {
            None
        } else {
            Some(self.pow(value, self.p - 2))
        }
    }

    pub(super) fn zero(&self) -> ReductionResidue {
        ReductionResidue { value: 0 }
    }

    #[cfg(test)]
    pub(super) fn one(&self) -> ReductionResidue {
        ReductionResidue { value: 1 }
    }
}

impl ReductionResidue {
    pub(super) fn representative(&self) -> u32 {
        self.value
    }

    pub(super) fn is_zero(&self) -> bool {
        self.value.is_zero()
    }
}

impl fmt::Display for ReductionPrimeError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Zero => write!(
                formatter,
                "reduction modulo p requires a prime modulus, not 0"
            ),
            Self::One => write!(
                formatter,
                "reduction modulo p requires a prime modulus greater than 1"
            ),
            Self::Composite { modulus } => write!(
                formatter,
                "reduction modulo p requires a prime modulus, but {modulus} is composite"
            ),
        }
    }
}

impl std::error::Error for ReductionPrimeError {}
