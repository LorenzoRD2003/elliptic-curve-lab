/// Kind of evidence supporting a candidate-refinement result.
///
/// This is not a probability score and does not certify the exact
/// endomorphism ring `End(E)`. It records how strong the evidence source is
/// within the current educational graph pipeline.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RefinementConfidence {
    /// Only the original Frobenius-compatible arithmetic candidate set `C₀`
    /// has been recorded.
    ArithmeticOnly,
    /// Candidates have survived conservative local evidence at the chosen
    /// prime `ℓ`, without eliminations from ambiguous graph relations.
    ConservativeLocalEvidence,
    /// Candidates have survived tentative graph evidence beyond the purely
    /// arithmetic `C₀` construction.
    TentativeGraphEvidence,
    /// Candidates have survived tentative graph evidence after monotone
    /// fixed-point propagation across the graph.
    ///
    /// This is stronger than one-hop evidence but still does not certify the
    /// exact endomorphism ring `End(E)`.
    PropagatedTentativeGraphEvidence,
}
