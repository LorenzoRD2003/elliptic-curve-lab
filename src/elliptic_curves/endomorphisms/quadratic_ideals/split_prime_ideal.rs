use num_bigint::BigUint;

use crate::elliptic_curves::endomorphisms::{
    quadratic_ideals::{PrimeNormIdealError, QuadraticPrimeBehavior},
    quadratic_orders::ImaginaryQuadraticOrder,
};

/// A selected split prime ideal of norm `ℓ` in an imaginary quadratic order.
///
/// This type models only the first explicit ideal family needed by the
/// horizontal-isogeny story. If `O_f` has discriminant `Δ` and an odd prime
/// `ℓ ∤ f` splits, then the two roots `r₁, r₂` of `x² ≡ Δ (mod ℓ)` correspond
/// to the two conjugate prime ideals above `ℓ`.
///
/// The stored `root` selects one of those two conjugate choices by recording a
/// square root `r` of the order discriminant modulo `ℓ`. The later
/// ideal-to-form bridge uses that root to choose the middle coefficient `b` of
/// a form with `b² ≡ Δ (mod 4ℓ)`.
///
/// This type does *not* implement ideal multiplication, ideal classes, or any
/// action on elliptic curves.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct SplitPrimeIdeal {
    order: ImaginaryQuadraticOrder,
    ell: BigUint,
    root: BigUint,
}

impl SplitPrimeIdeal {
    /// Builds a split prime ideal of norm `ℓ` by choosing one split root.
    ///
    /// The chosen root records one square root of the order discriminant
    /// modulo `ℓ`.
    ///
    /// Validation is delegated to [`ImaginaryQuadraticOrder::prime_behavior`]:
    /// the prime must split in `O_f`, must not divide the conductor, and the
    /// supplied `root` must be one of the two roots of `Δ` modulo `ℓ`.
    ///
    /// Complexity: dominated by [`ImaginaryQuadraticOrder::prime_behavior`].
    pub(crate) fn new(
        order: ImaginaryQuadraticOrder,
        ell: BigUint,
        root: BigUint,
    ) -> Result<Self, PrimeNormIdealError> {
        Self::validate_selected_root(&order, &ell, &root)?;
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

    /// Returns the selected split root `r`.
    pub(crate) fn root(&self) -> &BigUint {
        &self.root
    }

    /// Returns the conjugate split prime ideal above the same `ℓ`.
    ///
    /// If the stored root is `r`, the conjugate root is `ℓ - r`, using the
    /// canonical representatives returned by the prime-behavior layer.
    pub(crate) fn conjugate(&self) -> Self {
        Self {
            order: self.order.clone(),
            ell: self.ell.clone(),
            root: &self.ell - &self.root,
        }
    }

    fn validate_selected_root(
        order: &ImaginaryQuadraticOrder,
        ell: &BigUint,
        root: &BigUint,
    ) -> Result<(), PrimeNormIdealError> {
        match order.prime_behavior(ell)? {
            QuadraticPrimeBehavior::Split {
                roots: (left, right),
            } if root == &left || root == &right => Ok(()),
            QuadraticPrimeBehavior::Split { .. } => {
                Err(PrimeNormIdealError::RootDoesNotMatchPrimeBehavior)
            }
            QuadraticPrimeBehavior::NonInvertibleBecauseDividesConductor => {
                Err(PrimeNormIdealError::NonInvertibleBecauseDividesConductor)
            }
            QuadraticPrimeBehavior::Inert | QuadraticPrimeBehavior::Ramified { .. } => {
                Err(PrimeNormIdealError::NonSplitPrime)
            }
        }
    }
}
