use num_bigint::BigUint;

use crate::elliptic_curves::endomorphisms::{
    quadratic_ideals::{PrimeNormIdealError, QuadraticPrimeBehavior},
    quadratic_orders::ImaginaryQuadraticOrder,
};

/// The unique ramified prime ideal of norm `ℓ` in an imaginary quadratic order.
///
/// This crate-internal type models the invertible ramified case. If an odd
/// prime `ℓ ∤ f` ramifies in `O_f`, then there is a single prime ideal above
/// `ℓ`, and `(ℓ) = 𝔭²`. Unlike the split case, the root modulo `ℓ` is not a
/// choice between conjugate directions: it is the repeated local root of
/// `x² ≡ Δ (mod ℓ)`.
///
/// This type does *not* implement ideal multiplication, principal-ideal
/// checks, ideal classes, or any action on elliptic curves.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct RamifiedPrimeIdeal {
    order: ImaginaryQuadraticOrder,
    ell: BigUint,
    root: BigUint,
}

impl RamifiedPrimeIdeal {
    /// Builds the unique ramified prime ideal of norm `ℓ`.
    ///
    /// Validation is delegated to [`ImaginaryQuadraticOrder::prime_behavior`]:
    /// the prime must ramify in `O_f` and must not divide the conductor.
    ///
    /// Complexity: dominated by [`ImaginaryQuadraticOrder::prime_behavior`].
    pub(crate) fn new(
        order: ImaginaryQuadraticOrder,
        ell: BigUint,
    ) -> Result<Self, PrimeNormIdealError> {
        let root = Self::ramified_root(&order, &ell)?;
        Ok(Self { order, ell, root })
    }

    /// Returns the imaginary quadratic order containing the ideal.
    pub(crate) fn order(&self) -> &ImaginaryQuadraticOrder {
        &self.order
    }

    /// Returns the prime norm `ℓ`.
    pub(crate) fn norm(&self) -> &BigUint {
        &self.ell
    }

    /// Returns the repeated root of `Δ` modulo `ℓ`.
    pub(crate) fn root(&self) -> &BigUint {
        &self.root
    }

    /// Returns the conjugate ramified prime ideal.
    ///
    /// Ramified prime ideals are fixed by conjugation.
    pub(crate) fn conjugate(&self) -> Self {
        self.clone()
    }

    fn ramified_root(
        order: &ImaginaryQuadraticOrder,
        ell: &BigUint,
    ) -> Result<BigUint, PrimeNormIdealError> {
        match order.prime_behavior(ell)? {
            QuadraticPrimeBehavior::Ramified { root } => Ok(root),
            QuadraticPrimeBehavior::NonInvertibleBecauseDividesConductor => {
                Err(PrimeNormIdealError::NonInvertibleBecauseDividesConductor)
            }
            QuadraticPrimeBehavior::Split { .. } | QuadraticPrimeBehavior::Inert => {
                Err(PrimeNormIdealError::PrimeDoesNotRamify)
            }
        }
    }
}
