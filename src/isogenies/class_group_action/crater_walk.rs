use std::collections::{BTreeSet, HashMap, HashSet};

use crate::elliptic_curves::endomorphisms::quadratic_ideals::PrimeNormIdeal;
use crate::isogenies::graphs::{
    IsogenyGraphNodeId,
    endomorphisms::{CraterReport, CraterShape, HorizontalEdgeReport},
};

/// Deterministic walk on a certified crater, labeled by a prime-norm ideal.
///
/// The report records the graph-theoretic horizontal cycle seen in an
/// `ℓ`-volcano after the caller supplies an ideal of norm `ℓ`. The walk starts
/// at `start`, follows certified horizontal crater edges in a deterministic
/// local direction, and records the visited path. The local direction is
/// chosen from graph data: outgoing crater edges are ordered by target node and
/// edge id, and the walk avoids immediately backtracking when another outgoing
/// crater edge is available.
///
/// When the path returns to `start`, [`Self::cycle_length`] records the number
/// of horizontal steps. In that case [`Self::visited`] includes the closing
/// copy of `start`, so a 2-cycle is stored as `v₀, v₁, v₀`.
/// If the walk cannot start or cannot close a cycle, the report keeps the
/// maximal path found and [`Self::cycle_length`] returns `None`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CraterWalkReport {
    ideal: PrimeNormIdeal,
    crater_shape: CraterShape,
    visited: Vec<IsogenyGraphNodeId>,
}

impl CraterWalkReport {
    pub(crate) fn from_crater_report(
        crater: &CraterReport,
        ideal: PrimeNormIdeal,
        start: IsogenyGraphNodeId,
    ) -> Self {
        let mut walk = CraterWalk::from_crater_report(crater, start);
        let visited = walk.run();

        Self {
            ideal,
            crater_shape: crater.shape(),
            visited,
        }
    }

    /// Returns the ideal labeling this crater walk.
    pub fn ideal(&self) -> &PrimeNormIdeal {
        &self.ideal
    }

    /// Returns the requested starting node.
    pub fn start(&self) -> IsogenyGraphNodeId {
        self.visited[0]
    }

    /// Returns the certified crater shape that supplied the walk context.
    pub fn crater_shape(&self) -> CraterShape {
        self.crater_shape
    }

    /// Returns the path visited by the deterministic crater walk.
    ///
    /// If the walk closes a cycle, the final entry repeats [`Self::start`].
    pub fn visited(&self) -> &[IsogenyGraphNodeId] {
        &self.visited
    }

    /// Returns the certified cycle length when the walk returns to `start`.
    pub fn cycle_length(&self) -> Option<usize> {
        let closes_at_start =
            self.visited.len() > 1 && self.visited.last().copied() == Some(self.start());
        closes_at_start.then_some(self.visited.len() - 1)
    }

    /// Returns whether the recorded path closes back at its starting node.
    pub fn is_closed_cycle(&self) -> bool {
        self.cycle_length().is_some()
    }
}

struct CraterWalk {
    start: IsogenyGraphNodeId,
    crater_nodes: HashSet<IsogenyGraphNodeId>,
    outgoing: HashMap<IsogenyGraphNodeId, Vec<HorizontalEdgeReport>>,
}

impl CraterWalk {
    fn from_crater_report(crater: &CraterReport, start: IsogenyGraphNodeId) -> Self {
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

    fn run(&mut self) -> Vec<IsogenyGraphNodeId> {
        if !self.crater_nodes.contains(&self.start) {
            return vec![self.start];
        }

        let mut visited = vec![self.start];
        let mut seen = BTreeSet::from([self.start]);
        let mut previous = None;
        let mut current = self.start;

        loop {
            let Some(edges) = self.outgoing.get(&current) else {
                return visited;
            };
            let Some(next) =
                Self::choose_forward_crater_edge(edges, previous).map(|edge| edge.target())
            else {
                return visited;
            };

            visited.push(next);

            if next == self.start {
                return visited;
            }

            if !seen.insert(next) {
                return visited;
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
