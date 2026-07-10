use std::collections::{BTreeSet, HashMap, HashSet};

use crate::isogenies::{
    class_group_action::crater_walk::CraterWalkTermination,
    graphs::{
        IsogenyGraphNodeId,
        endomorphisms::{CraterReport, HorizontalEdgeReport},
    },
};

pub(super) struct CraterWalkRun {
    visited: Vec<IsogenyGraphNodeId>,
    start_in_crater: bool,
    termination: CraterWalkTermination,
}

impl CraterWalkRun {
    pub(super) fn into_parts(self) -> (Vec<IsogenyGraphNodeId>, bool, CraterWalkTermination) {
        (self.visited, self.start_in_crater, self.termination)
    }
}

pub(super) struct CraterWalk {
    start: IsogenyGraphNodeId,
    crater_nodes: HashSet<IsogenyGraphNodeId>,
    outgoing: HashMap<IsogenyGraphNodeId, Vec<HorizontalEdgeReport>>,
}

impl CraterWalk {
    pub(super) fn from_crater_report(crater: &CraterReport, start: IsogenyGraphNodeId) -> Self {
        let crater_nodes = crater.nodes().iter().copied().collect::<HashSet<_>>();
        let mut outgoing = crater.certified_internal_outgoing_edge_map();

        for edges in outgoing.values_mut() {
            edges.sort_by_key(|edge| (edge.target(), edge.edge_id()));
        }

        Self {
            start,
            crater_nodes,
            outgoing,
        }
    }

    pub(super) fn run(&mut self) -> CraterWalkRun {
        if !self.crater_nodes.contains(&self.start) {
            return CraterWalkRun {
                visited: vec![self.start],
                start_in_crater: false,
                termination: CraterWalkTermination::StartOutsideCrater,
            };
        }

        let mut visited = vec![self.start];
        let mut seen = BTreeSet::from([self.start]);
        let mut previous = None;
        let mut current = self.start;

        loop {
            let Some(edges) = self.outgoing.get(&current) else {
                return CraterWalkRun {
                    visited,
                    start_in_crater: true,
                    termination: CraterWalkTermination::NoCertifiedOutgoingEdge,
                };
            };
            let Some(next) =
                Self::choose_forward_crater_edge(edges, previous).map(|edge| edge.target())
            else {
                return CraterWalkRun {
                    visited,
                    start_in_crater: true,
                    termination: CraterWalkTermination::NoCertifiedOutgoingEdge,
                };
            };

            visited.push(next);

            if next == self.start {
                return CraterWalkRun {
                    visited,
                    start_in_crater: true,
                    termination: CraterWalkTermination::ClosedCycle,
                };
            }

            if !seen.insert(next) {
                return CraterWalkRun {
                    visited,
                    start_in_crater: true,
                    termination: CraterWalkTermination::RepeatedNonStartNode,
                };
            }

            previous = Some(current);
            current = next;
        }
    }

    /// Chooses the next certified crater edge in the deterministic local direction.
    ///
    /// The input edges are expected to have already been sorted by target node
    /// and edge id. The helper first chooses the earliest edge whose target is
    /// not the previous node, so a walk through a crater keeps moving forward
    /// when the local graph offers a non-backtracking choice. If every
    /// available edge returns to the previous node, it chooses the first edge,
    /// which lets degenerate two-vertex craters close as `v₀ → v₁ → v₀`.
    fn choose_forward_crater_edge(
        edges: &[HorizontalEdgeReport],
        previous: Option<IsogenyGraphNodeId>,
    ) -> Option<&HorizontalEdgeReport> {
        edges
            .iter()
            .find(|edge| Some(edge.target()) != previous)
            .or_else(|| edges.first())
    }
}
