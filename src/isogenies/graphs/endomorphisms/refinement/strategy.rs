/// Policy for deciding which local graph evidence may eliminate candidates.
///
/// A strategy controls how aggressively a refinement run shrinks
/// `C₀ = {O_f : ℤ[π] ⊆ O_f ⊆ O_K}`. The default is intentionally conservative:
/// ambiguous graph evidence is recorded but not used to discard candidates.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum CandidateRefinementStrategy {
    /// Uses local node evidence and only unequivocal incident-edge evidence.
    ///
    /// This strategy never eliminates an order using an `Ambiguous` or
    /// `Unsupported` tentative edge relation. It is the default because the
    /// current graph-side relations are educational evidence, not a
    /// certification of `End(E)`.
    #[default]
    Conservative,
    /// Uses only the node's own local `ℓ`-level evidence.
    ///
    /// This ignores incident edge relations, so it is useful as the smallest
    /// explainable refinement pass and as a baseline for tests.
    NodeLocalLevelsOnly,
    /// Uses node-local evidence and all incident edge relations with a single
    /// tentative direction.
    ///
    /// `PossiblyHorizontal`, `PossiblyAscending`, and `PossiblyDescending`
    /// may eliminate candidates. `Ambiguous` and `Unsupported` still do not.
    IncidentUnambiguousEdges,
}
