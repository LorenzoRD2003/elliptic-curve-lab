use crate::elliptic_curves::endomorphisms::candidate_sets::EndomorphismRingCandidateSet;
use crate::numerics::PositivePrimeError;
use num_bigint::BigUint;

/// Tentative local endomorphism-side relation attached to one isogeny edge.
///
/// This enum is intentionally not a definitive classification of the edge.
/// It records only what remains possible after comparing the `ℓ`-local
/// candidate levels attached to the source and target endomorphism-order
/// candidate sets.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) enum IsogenyEdgeEndomorphismRelation {
    PossiblyHorizontal,
    PossiblyAscending,
    PossiblyDescending,
    Ambiguous,
    Unsupported,
}

/// Tentative arithmetic report comparing the source and target candidate
/// endomorphism orders of one edge at a fixed prime `ℓ`.
///
/// This report does **not** certify the exact endomorphism ring of either
/// endpoint, nor does it definitively classify the edge. It only compares the
/// locally possible candidate levels derived from Frobenius-compatible
/// imaginary quadratic orders.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct IsogenyEdgeEndomorphismReport {
    prime: BigUint,
    source_possible_levels: Vec<u32>,
    target_possible_levels: Vec<u32>,
    relation: IsogenyEdgeEndomorphismRelation,
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
            IsogenyEdgeEndomorphismRelation::Unsupported
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
    pub(crate) fn prime(&self) -> &BigUint {
        &self.prime
    }

    /// Returns the distinct source levels compatible with the source candidate set.
    pub(crate) fn source_possible_levels(&self) -> &[u32] {
        &self.source_possible_levels
    }

    /// Returns the distinct target levels compatible with the target candidate set.
    pub(crate) fn target_possible_levels(&self) -> &[u32] {
        &self.target_possible_levels
    }

    /// Returns the tentative local endomorphism-side relation for the edge.
    pub(crate) fn relation(&self) -> &IsogenyEdgeEndomorphismRelation {
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
) -> IsogenyEdgeEndomorphismRelation {
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
        0 => IsogenyEdgeEndomorphismRelation::Unsupported,
        1 if horizontal => IsogenyEdgeEndomorphismRelation::PossiblyHorizontal,
        1 if ascending => IsogenyEdgeEndomorphismRelation::PossiblyAscending,
        1 if descending => IsogenyEdgeEndomorphismRelation::PossiblyDescending,
        _ => IsogenyEdgeEndomorphismRelation::Ambiguous,
    }
}

#[cfg(test)]
mod tests {

    use num_bigint::BigUint;

    use super::infer_relation;
    use crate::elliptic_curves::endomorphisms::{
        candidate_sets::EndomorphismRingCandidateSet, quadratic_orders::QuadraticDiscriminant,
    };
    use crate::isogenies::graphs::{
        IsogenyEdgeEndomorphismRelation, IsogenyEdgeEndomorphismReport,
    };

    fn candidate_set(discriminant: i64) -> EndomorphismRingCandidateSet {
        QuadraticDiscriminant::new(discriminant)
            .factorization()
            .expect("test discriminant should factor canonically")
            .endomorphism_ring_candidates()
            .expect("candidate orders should construct")
    }

    #[test]
    fn identical_singleton_levels_are_possibly_horizontal() {
        let source = candidate_set(-4);
        let target = candidate_set(-4);

        let report = IsogenyEdgeEndomorphismReport::from_candidate_sets(
            &BigUint::from(2u8),
            &source,
            &target,
        )
        .expect("report should build");

        assert_eq!(report.source_possible_levels(), &[0]);
        assert_eq!(report.target_possible_levels(), &[0]);
        assert_eq!(
            report.relation(),
            &IsogenyEdgeEndomorphismRelation::PossiblyHorizontal
        );
    }

    #[test]
    fn different_quadratic_fields_are_unsupported() {
        let source = candidate_set(-16);
        let target = candidate_set(-3);

        let report = source
            .tentative_edge_endomorphism_report(&BigUint::from(2u8), &target)
            .expect("report should build");

        assert_eq!(report.source_possible_levels(), &[0, 1]);
        assert_eq!(report.target_possible_levels(), &[0]);
        assert_eq!(
            report.relation(),
            &IsogenyEdgeEndomorphismRelation::Unsupported
        );
    }

    #[test]
    fn mixed_horizontal_and_vertical_possibilities_are_ambiguous() {
        let source = candidate_set(-16);
        let target = candidate_set(-16);

        let report = IsogenyEdgeEndomorphismReport::from_candidate_sets(
            &BigUint::from(2u8),
            &source,
            &target,
        )
        .expect("report should build");

        assert_eq!(report.source_possible_levels(), &[0, 1]);
        assert_eq!(report.target_possible_levels(), &[0, 1]);
        assert_eq!(
            report.relation(),
            &IsogenyEdgeEndomorphismRelation::Ambiguous
        );
    }

    #[test]
    fn level_classifier_recovers_all_tentative_variants() {
        assert_eq!(
            infer_relation(&[0], &[0]),
            IsogenyEdgeEndomorphismRelation::PossiblyHorizontal
        );
        assert_eq!(
            infer_relation(&[1], &[0]),
            IsogenyEdgeEndomorphismRelation::PossiblyAscending
        );
        assert_eq!(
            infer_relation(&[0], &[1]),
            IsogenyEdgeEndomorphismRelation::PossiblyDescending
        );
        assert_eq!(
            infer_relation(&[0, 1], &[0]),
            IsogenyEdgeEndomorphismRelation::Ambiguous
        );
        assert_eq!(
            infer_relation(&[0], &[2]),
            IsogenyEdgeEndomorphismRelation::Unsupported
        );
    }
}
