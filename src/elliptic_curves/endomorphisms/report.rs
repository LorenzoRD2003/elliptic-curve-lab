use crate::elliptic_curves::endomorphisms::{
    EndomorphismRingCandidateSet, ImaginaryQuadraticOrder, ImaginaryQuadraticOrderError,
    QuadraticDiscriminantFactorization,
};
use crate::elliptic_curves::frobenius::{FrobeniusCurveType, FrobeniusDiscriminant};

/// Frobenius-compatible endomorphism-ring report.
///
/// The ordinary case is modeled through imaginary quadratic orders compatible
/// with the Frobenius data. The supersingular case is kept separate because it
/// does not belong to the same quadratic-order story.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum EndomorphismRingReport {
    /// Frobenius-compatible candidate data for the ordinary case.
    ///
    /// Starting from the Frobenius discriminant `Δ_π = t^2 - 4q`, this branch:
    ///
    /// - factors `Δ_π = v^2 D_K`
    /// - recovers the Frobenius-generated order `ℤ[π]`
    /// - recovers the maximal order `O_K`
    /// - enumerates every intermediate quadratic order `O_f` with `f | v`
    ///
    /// It is intentionally honest: it does not claim to compute the actual
    /// ring `End(E)`. Instead, it records the quadratic orders that remain
    /// compatible with the current ordinary Frobenius data, together with the
    /// certified sandwich `ℤ[π] ⊆ End(E) ⊆ O_K`.
    OrdinaryQuadraticOrderCandidates {
        frobenius_discriminant: FrobeniusDiscriminant,
        factorization: QuadraticDiscriminantFactorization,
        candidate_set: EndomorphismRingCandidateSet,
    },

    /// Frobenius-side placeholder for the supersingular case.
    ///
    /// This branch is also intentionally honest: the current `endomorphisms`
    /// module does not model the supersingular endomorphism ring as an
    /// imaginary quadratic order. In that case the geometric endomorphism
    /// algebra lives in a quaternion algebra, so the ordinary
    /// quadratic-order pipeline does not apply.
    SupersingularQuaternionicPlaceholder {
        frobenius_discriminant: FrobeniusDiscriminant,
    },
}

impl EndomorphismRingReport {
    /// Builds the report from one Frobenius discriminant package.
    ///
    /// Complexity:
    /// - `Θ(1)` in the supersingular branch
    /// - dominated by `num-prime` in the ordinary branch
    pub fn from_frobenius_discriminant(
        frobenius_discriminant: FrobeniusDiscriminant,
    ) -> Result<Self, ImaginaryQuadraticOrderError> {
        match frobenius_discriminant.frobenius_trace().curve_type() {
            FrobeniusCurveType::Ordinary => {
                let factorization = frobenius_discriminant
                    .quadratic_factorization()
                    .map_err(|_| ImaginaryQuadraticOrderError::NonImaginaryOrderDiscriminant)?;
                let candidate_set = factorization.endomorphism_ring_candidates()?;

                Ok(Self::OrdinaryQuadraticOrderCandidates {
                    frobenius_discriminant,
                    factorization,
                    candidate_set,
                })
            }
            FrobeniusCurveType::Supersingular => Ok(Self::SupersingularQuaternionicPlaceholder {
                frobenius_discriminant,
            }),
        }
    }

    /// Returns the underlying Frobenius discriminant package.
    pub fn frobenius_discriminant(&self) -> &FrobeniusDiscriminant {
        match self {
            Self::OrdinaryQuadraticOrderCandidates {
                frobenius_discriminant,
                ..
            }
            | Self::SupersingularQuaternionicPlaceholder {
                frobenius_discriminant,
            } => frobenius_discriminant,
        }
    }

    /// Returns whether the report is in the ordinary branch.
    pub fn is_ordinary(&self) -> bool {
        matches!(self, Self::OrdinaryQuadraticOrderCandidates { .. })
    }

    /// Returns whether the report is in the supersingular branch.
    pub fn is_supersingular(&self) -> bool {
        matches!(self, Self::SupersingularQuaternionicPlaceholder { .. })
    }

    /// Returns the canonical factorization `Δ_π = v^2 D_K` in the ordinary case.
    pub fn factorization(&self) -> Option<&QuadraticDiscriminantFactorization> {
        match self {
            Self::OrdinaryQuadraticOrderCandidates { factorization, .. } => Some(factorization),
            Self::SupersingularQuaternionicPlaceholder { .. } => None,
        }
    }

    /// Returns the candidate orders compatible with Frobenius in the ordinary case.
    pub fn candidate_set(&self) -> Option<&EndomorphismRingCandidateSet> {
        match self {
            Self::OrdinaryQuadraticOrderCandidates { candidate_set, .. } => Some(candidate_set),
            Self::SupersingularQuaternionicPlaceholder { .. } => None,
        }
    }

    /// Returns the Frobenius-generated order `ℤ[π]` in the ordinary case.
    pub fn frobenius_order(&self) -> Option<&ImaginaryQuadraticOrder> {
        self.candidate_set()
            .map(EndomorphismRingCandidateSet::frobenius_order)
    }

    /// Returns the maximal order `O_K` in the ordinary case.
    pub fn maximal_order(&self) -> Option<&ImaginaryQuadraticOrder> {
        self.candidate_set()
            .map(EndomorphismRingCandidateSet::maximal_order)
    }

    /// Returns every candidate order compatible with Frobenius in the ordinary case.
    ///
    /// This is a candidate list, not a certification of which one equals
    /// `End(E)`.
    pub fn candidate_orders(&self) -> Option<&[ImaginaryQuadraticOrder]> {
        self.candidate_set()
            .map(EndomorphismRingCandidateSet::candidate_orders)
    }

    /// Returns the number of Frobenius-compatible candidate orders in the ordinary case.
    pub fn candidate_count(&self) -> Option<usize> {
        self.candidate_set().map(EndomorphismRingCandidateSet::len)
    }

    /// Returns whether the certified sandwich `ℤ[π] ⊆ O_K` holds in the
    /// ordinary case.
    ///
    /// This is the order-theoretic part that the current Frobenius-side
    /// arithmetic really certifies before any later algorithm singles out one
    /// candidate for `End(E)`.
    pub fn sandwich_inclusion_holds(&self) -> Option<bool> {
        match (self.frobenius_order(), self.maximal_order()) {
            (Some(frobenius_order), Some(maximal_order)) => {
                Some(frobenius_order.is_suborder_of(maximal_order))
            }
            _ => None,
        }
    }
}
