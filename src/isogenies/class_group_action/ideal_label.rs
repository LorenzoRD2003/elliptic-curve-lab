use core::fmt;

use num_bigint::{BigInt, BigUint};

use crate::elliptic_curves::endomorphisms::{
    binary_quadratic_forms::QuadraticClassGroup, quadratic_ideals::PrimeNormIdeal,
};
use crate::isogenies::graphs::endomorphisms::CraterReport;

/// Local behavior of a supplied prime-norm ideal in a crater label.
///
/// This mirrors the already-validated family of [`PrimeNormIdeal`]: an ideal
/// value can currently be split or ramified, but not inert or conductor
/// dividing.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum CraterIdealPrimeBehavior {
    /// The ideal is one of the two split primes above `ℓ`.
    Split,
    /// The ideal is the unique ramified prime above `ℓ`.
    Ramified,
}

/// Failure modes for labeling a crater by a supplied prime-norm ideal.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) enum CraterIdealLabelError {
    /// The crater was built at a different local prime than the ideal norm.
    PrimeNormMismatch {
        ideal_norm: BigUint,
        crater_prime: BigUint,
    },
    /// The ideal order and the class group have different discriminants.
    OrderDiscriminantMismatch {
        ideal_discriminant: BigInt,
        class_group_discriminant: BigInt,
    },
}

impl fmt::Display for CraterIdealLabelError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::PrimeNormMismatch {
                ideal_norm,
                crater_prime,
            } => write!(
                formatter,
                "prime-norm ideal has norm {ideal_norm}, but the crater prime is {crater_prime}"
            ),
            Self::OrderDiscriminantMismatch {
                ideal_discriminant,
                class_group_discriminant,
            } => write!(
                formatter,
                "ideal order has discriminant {ideal_discriminant}, but the class group has discriminant {class_group_discriminant}"
            ),
        }
    }
}

impl std::error::Error for CraterIdealLabelError {}

/// Certified local compatibility between one crater and one prime-norm ideal.
///
/// This report certifies only that a caller-supplied prime-norm ideal can label
/// a crater at the same local prime `ℓ` and in the same quadratic order as a
/// supplied class group. It does not orient crater edges, infer an ideal from
/// the graph, construct the associated binary-quadratic-form class, or compute
/// a class-group action.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct CraterIdealLabelReport {
    ideal: PrimeNormIdeal,
    crater_prime: BigUint,
    prime_behavior: CraterIdealPrimeBehavior,
}

impl CraterIdealLabelReport {
    /// Certifies local compatibility of `ideal` with `crater` and `class_group`.
    ///
    /// The input ideal is already a validated [`PrimeNormIdeal`], so this
    /// constructor does not re-run local splitting or conductor-divisibility
    /// checks. It only checks the crater prime and the class-group
    /// discriminant, then records whether the validated ideal is split or
    /// ramified.
    ///
    /// Complexity: `Θ(1)` exact integer comparisons and cloning of the stored
    /// ideal metadata.
    pub(crate) fn new(
        crater: &CraterReport,
        class_group: &QuadraticClassGroup,
        ideal: PrimeNormIdeal,
    ) -> Result<Self, CraterIdealLabelError> {
        if ideal.norm() != crater.prime() {
            return Err(CraterIdealLabelError::PrimeNormMismatch {
                ideal_norm: ideal.norm().clone(),
                crater_prime: crater.prime().clone(),
            });
        }

        let ideal_discriminant = ideal.order().discriminant().value();
        let class_group_discriminant = class_group.discriminant().value();
        if ideal_discriminant != class_group_discriminant {
            return Err(CraterIdealLabelError::OrderDiscriminantMismatch {
                ideal_discriminant: ideal_discriminant.clone(),
                class_group_discriminant: class_group_discriminant.clone(),
            });
        }

        let prime_behavior = if ideal.is_split() {
            CraterIdealPrimeBehavior::Split
        } else {
            debug_assert!(ideal.is_ramified());
            CraterIdealPrimeBehavior::Ramified
        };

        Ok(Self {
            ideal,
            crater_prime: crater.prime().clone(),
            prime_behavior,
        })
    }

    /// Returns the supplied prime-norm ideal.
    pub(crate) fn ideal(&self) -> &PrimeNormIdeal {
        &self.ideal
    }

    /// Returns the local crater prime `ℓ`.
    pub(crate) fn crater_prime(&self) -> &BigUint {
        &self.crater_prime
    }

    /// Returns whether the supplied ideal is split or ramified.
    pub(crate) fn prime_behavior(&self) -> CraterIdealPrimeBehavior {
        self.prime_behavior
    }
}
