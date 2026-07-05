use std::collections::BTreeSet;

use num_bigint::BigUint;

use crate::isogenies::graphs::{
    IsogenyGraphEdgeId, IsogenyGraphNodeId, endomorphisms::IsogenyEdgeEndomorphismTentativeRelation,
};

/// Provenance of a local endomorphism constraint used in a refinement run.
///
/// The provenance is intentionally structural rather than textual, so examples
/// and visualization code can explain where a constraint came from without
/// parsing a free-form string.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ConstraintSource {
    /// The constraint was read from one stored node report.
    NodeReport { node_id: IsogenyGraphNodeId },
    /// The constraint was read from one stored edge report.
    EdgeReport { edge_id: IsogenyGraphEdgeId },
}

/// One local constraint used to filter endomorphism-ring candidates.
///
/// A constraint is local to the chosen prime `ℓ` and to the evidence already
/// present in the graph-side endomorphism report. It narrows candidates by
/// compatibility, but it is not a certificate that a surviving order equals
/// `End(E)`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum LocalEndomorphismConstraint {
    /// The candidate order must have local conductor level
    /// `v_ℓ(f) ∈ allowed_levels`.
    NodeLevel {
        /// The chosen local prime `ℓ`.
        ell: BigUint,
        /// Allowed local levels for `O_f = ℤ + fO_K`.
        allowed_levels: BTreeSet<u32>,
        /// Where this constraint came from in the graph report.
        provenance: ConstraintSource,
    },
    /// The source and target candidate levels must be compatible with one
    /// tentative edge relation.
    EdgeRelation {
        /// The chosen local prime `ℓ`.
        ell: BigUint,
        /// Source node of the observed edge.
        source_node: IsogenyGraphNodeId,
        /// Target node of the observed edge.
        target_node: IsogenyGraphNodeId,
        /// Tentative relation observed on the edge.
        relation: IsogenyEdgeEndomorphismTentativeRelation,
        /// Where this constraint came from in the graph report.
        provenance: ConstraintSource,
    },
}
