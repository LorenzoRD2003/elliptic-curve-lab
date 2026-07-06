use num_bigint::BigUint;

use crate::elliptic_curves::endomorphisms::candidate_sets::EndomorphismRingCandidateSet;
use crate::numerics::PositivePrimeError;

/// Tentative local endomorphism-side relation attached to one isogeny edge.
///
/// This enum is intentionally not a definitive classification of the edge.
/// It records only what remains possible after comparing the `ℓ`-local
/// candidate levels attached to the source and target endomorphism-order
/// candidate sets.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum IsogenyEdgeEndomorphismTentativeRelation {
    PossiblyHorizontal,
    PossiblyAscending,
    PossiblyDescending,
    Ambiguous,
    Unsupported,
}

impl IsogenyEdgeEndomorphismTentativeRelation {
    /// Returns this relation when it has a single tentative direction.
    ///
    /// `Ambiguous` and `Unsupported` are intentionally excluded because the
    /// conservative refinement pass must not eliminate candidates from evidence
    /// that does not point to one local edge shape.
    pub(crate) fn as_unambiguous(&self) -> Option<Self> {
        match self {
            Self::PossiblyHorizontal | Self::PossiblyAscending | Self::PossiblyDescending => {
                Some(self.clone())
            }
            Self::Ambiguous | Self::Unsupported => None,
        }
    }

    /// Returns whether source and target local levels are compatible with this
    /// tentative relation.
    ///
    /// The levels are arithmetic conductor levels `v_ℓ(f)` on the source and
    /// target candidate orders. Ambiguous or unsupported relations return
    /// `false`; callers that want conservative behavior should first use
    /// [`Self::as_unambiguous`].
    pub(crate) fn allows_levels(&self, source_level: u32, target_level: u32) -> bool {
        match self {
            Self::PossiblyHorizontal => source_level == target_level,
            Self::PossiblyAscending => source_level.checked_sub(1) == Some(target_level),
            Self::PossiblyDescending => target_level.checked_sub(1) == Some(source_level),
            Self::Ambiguous | Self::Unsupported => false,
        }
    }
}

/// Tentative arithmetic report comparing the source and target candidate
/// endomorphism orders of one edge at a fixed prime `ℓ`.
///
/// This report does **not** certify the exact endomorphism ring of either
/// endpoint, nor does it definitively classify the edge. It only compares the
/// locally possible candidate levels derived from Frobenius-compatible
/// imaginary quadratic orders.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct IsogenyEdgeEndomorphismReport {
    prime: BigUint,
    source_possible_levels: Vec<u32>,
    target_possible_levels: Vec<u32>,
    relation: IsogenyEdgeEndomorphismTentativeRelation,
}

impl IsogenyEdgeEndomorphismReport {
    /// Builds the tentative local edge report from the source and target
    /// candidate endomorphism-order sets.
    ///
    /// The current rule looks only at level differences `0`, `+1`, and `-1`
    /// between the source and target possible local levels:
    ///
    /// - `0` contributes a horizontal possibility
    /// - `+1` contributes a descending possibility
    /// - `-1` contributes an ascending possibility
    ///
    /// If more than one of those possibilities remains, the result is
    /// [`IsogenyEdgeEndomorphismRelation::Ambiguous`]. If none remain, the
    /// current arithmetic layer marks the relation as
    /// [`IsogenyEdgeEndomorphismRelation::Unsupported`].
    ///
    /// Complexity: dominated by `num-prime`, plus `Θ(st)` integer comparisons,
    /// where `s` and `t` are the numbers of distinct possible source and target
    /// levels.
    pub(crate) fn from_candidate_sets(
        prime: &BigUint,
        source: &EndomorphismRingCandidateSet,
        target: &EndomorphismRingCandidateSet,
    ) -> Result<Self, PositivePrimeError> {
        let source_possible_levels = source.distinct_levels(prime)?;
        let target_possible_levels = target.distinct_levels(prime)?;
        let relation = if source.fundamental_discriminant() != target.fundamental_discriminant() {
            IsogenyEdgeEndomorphismTentativeRelation::Unsupported
        } else {
            infer_relation(&source_possible_levels, &target_possible_levels)
        };

        Ok(Self {
            prime: prime.clone(),
            source_possible_levels,
            target_possible_levels,
            relation,
        })
    }

    /// Returns the chosen prime `ℓ`.
    pub fn prime(&self) -> &BigUint {
        &self.prime
    }

    /// Returns the distinct source levels compatible with the source candidate set.
    pub fn source_possible_levels(&self) -> &[u32] {
        &self.source_possible_levels
    }

    /// Returns the distinct target levels compatible with the target candidate set.
    pub fn target_possible_levels(&self) -> &[u32] {
        &self.target_possible_levels
    }

    /// Returns the tentative local endomorphism-side relation for the edge.
    pub fn relation(&self) -> &IsogenyEdgeEndomorphismTentativeRelation {
        &self.relation
    }
}

impl EndomorphismRingCandidateSet {
    /// Builds the tentative local edge report from `self` to `target`.
    ///
    /// Complexity: dominated by `num-prime`, plus `Θ(st)` integer comparisons,
    /// where `s` and `t` are the numbers of distinct possible source and target
    /// levels.
    pub(crate) fn tentative_edge_endomorphism_report(
        &self,
        prime: &BigUint,
        target: &EndomorphismRingCandidateSet,
    ) -> Result<IsogenyEdgeEndomorphismReport, PositivePrimeError> {
        IsogenyEdgeEndomorphismReport::from_candidate_sets(prime, self, target)
    }

    fn distinct_levels(&self, prime: &BigUint) -> Result<Vec<u32>, PositivePrimeError> {
        let mut levels = self
            .volcanic_level_candidates_at(prime)?
            .into_iter()
            .map(|candidate| candidate.level())
            .collect::<Vec<_>>();
        levels.sort_unstable();
        levels.dedup();
        Ok(levels)
    }
}

fn infer_relation(
    source_possible_levels: &[u32],
    target_possible_levels: &[u32],
) -> IsogenyEdgeEndomorphismTentativeRelation {
    let mut horizontal = false;
    let mut ascending = false;
    let mut descending = false;

    for &source_level in source_possible_levels {
        for &target_level in target_possible_levels {
            if source_level == target_level {
                horizontal = true;
            } else if source_level == target_level + 1 {
                ascending = true;
            } else if target_level == source_level + 1 {
                descending = true;
            }
        }
    }

    let possibility_count =
        usize::from(horizontal) + usize::from(ascending) + usize::from(descending);

    match possibility_count {
        0 => IsogenyEdgeEndomorphismTentativeRelation::Unsupported,
        1 if horizontal => IsogenyEdgeEndomorphismTentativeRelation::PossiblyHorizontal,
        1 if ascending => IsogenyEdgeEndomorphismTentativeRelation::PossiblyAscending,
        1 if descending => IsogenyEdgeEndomorphismTentativeRelation::PossiblyDescending,
        _ => IsogenyEdgeEndomorphismTentativeRelation::Ambiguous,
    }
}

#[cfg(test)]
mod tests;
