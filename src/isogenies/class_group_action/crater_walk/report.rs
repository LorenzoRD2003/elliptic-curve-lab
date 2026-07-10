use crate::elliptic_curves::endomorphisms::quadratic_ideals::PrimeNormIdeal;
use crate::isogenies::{
    class_group_action::crater_walk::{CraterWalkTermination, engine::CraterWalk},
    graphs::{
        IsogenyGraphNodeId,
        endomorphisms::{CraterReport, CraterShape},
    },
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
///
/// Complexity: building the local outgoing-edge map is linear in the certified
/// crater evidence, sorting each local edge list costs `O(E_c log E_c)` in the
/// worst case for `E_c` certified internal crater edges, and the deterministic
/// walk itself is `O(L)` for the observed path length `L`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CraterWalkReport {
    ideal: PrimeNormIdeal,
    crater_shape: CraterShape,
    visited: Vec<IsogenyGraphNodeId>,
    start_in_crater: bool,
    termination: CraterWalkTermination,
}

impl CraterWalkReport {
    pub(crate) fn from_crater_report(
        crater: &CraterReport,
        ideal: PrimeNormIdeal,
        start: IsogenyGraphNodeId,
    ) -> Self {
        let mut walk = CraterWalk::from_crater_report(crater, start);
        let run = walk.run();
        let (visited, start_in_crater, termination) = run.into_parts();

        Self {
            ideal,
            crater_shape: crater.shape(),
            visited,
            start_in_crater,
            termination,
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

    /// Returns whether the requested start node belonged to the certified crater.
    pub fn start_in_crater(&self) -> bool {
        self.start_in_crater
    }

    /// Returns the graph-theoretic reason why the deterministic walk stopped.
    pub fn termination(&self) -> CraterWalkTermination {
        self.termination
    }

    /// Returns the certified cycle length when the walk returns to `start`.
    pub fn cycle_length(&self) -> Option<usize> {
        (self.termination == CraterWalkTermination::ClosedCycle).then_some(self.visited.len() - 1)
    }

    /// Returns whether the recorded path closes back at its starting node.
    pub fn is_closed_cycle(&self) -> bool {
        self.termination == CraterWalkTermination::ClosedCycle
    }
}
