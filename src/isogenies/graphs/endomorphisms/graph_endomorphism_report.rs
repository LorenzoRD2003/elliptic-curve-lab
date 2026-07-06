//! Aggregate graph-side endomorphism reports.
//!
//! This module keeps the public report vocabulary small while delegating the
//! construction, refinement, observed-volcano evidence, and tests to focused
//! sibling files.

mod build;
mod edge;
mod node;
mod observed_graph_evidence;
mod observed_volcano_evidence;
mod refinement;

#[cfg(test)]
mod tests;

use num_bigint::BigUint;

pub use edge::IsogenyGraphEndomorphismEdgeReport;
pub use node::IsogenyGraphEndomorphismNodeReport;

use crate::isogenies::graphs::IsogenyGraphNodeId;

/// Endomorphism-side report for an entire educational isogeny graph at one chosen prime `ℓ`.
///
/// This aggregate report is still conservative. It packages:
///
/// - automatic Frobenius-compatible candidate-order data for each node
/// - the corresponding `ℓ`-local candidate levels at each node
/// - graph-observed volcano evidence from certified floor distances when
///   available, with a weak-BFS fallback for older educational cases
/// - tentative edge relations derived from those node-wise candidate sets
///
/// It does **not** certify exact endomorphism rings or definitive edge types.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct IsogenyGraphEndomorphismReport {
    prime: BigUint,
    nodes: Vec<IsogenyGraphEndomorphismNodeReport>,
    edges: Vec<IsogenyGraphEndomorphismEdgeReport>,
}

impl IsogenyGraphEndomorphismReport {
    pub(crate) fn new(
        prime: BigUint,
        nodes: Vec<IsogenyGraphEndomorphismNodeReport>,
        edges: Vec<IsogenyGraphEndomorphismEdgeReport>,
    ) -> Self {
        Self {
            prime,
            nodes,
            edges,
        }
    }

    /// Returns the chosen prime `ℓ`.
    pub fn prime(&self) -> &BigUint {
        &self.prime
    }

    /// Returns the node reports in dense node-id order.
    pub fn nodes(&self) -> &[IsogenyGraphEndomorphismNodeReport] {
        &self.nodes
    }

    /// Returns the edge reports in stored edge order.
    pub fn edges(&self) -> &[IsogenyGraphEndomorphismEdgeReport] {
        &self.edges
    }

    /// Returns the node report for the requested id when present.
    pub fn node_report(
        &self,
        node_id: IsogenyGraphNodeId,
    ) -> Option<&IsogenyGraphEndomorphismNodeReport> {
        self.nodes
            .get(node_id.0)
            .filter(|report| report.node_id() == node_id)
    }
}
