/// Policy for deciding which local graph evidence may eliminate candidates.
///
/// A strategy controls how aggressively a refinement run shrinks
/// `C₀ = {O_f : ℤ[π] ⊆ O_f ⊆ O_K}`. The default is intentionally conservative:
/// ambiguous graph evidence is recorded but not used to discard candidates.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum CandidateRefinementStrategy {
    /// Uses local node evidence and only unequivocal graph evidence.
    ///
    /// This strategy never eliminates an order using an `Ambiguous` or
    /// `Unsupported` tentative edge relation. Today this includes conservative
    /// observed node-level evidence and incident edges with a single tentative
    /// direction. It is the default because the current graph-side relations
    /// are educational evidence, not a certification of `End(E)`.
    #[default]
    Conservative,
    /// Uses only the node's own local `ℓ`-level evidence.
    ///
    /// This ignores incident edge relations, so it is useful as the smallest
    /// explainable refinement pass and as a baseline for tests.
    NodeLocalLevelsOnly,
}
