use std::collections::{BTreeMap, BTreeSet, HashSet};

use crate::isogenies::{
    class_group_action::CraterOrientationWitnessError,
    graphs::{IsogenyGraphNodeId, endomorphisms::CraterReport},
};

/// User-supplied orientation of one certified crater cycle.
///
/// The witness records only graph orientation data: for each certified crater
/// node it names the next node in the chosen positive direction. Construction
/// validates that every declared step follows a certified internal horizontal
/// crater edge and that the successors form one closed cycle. It does not
/// certify that this direction is the arithmetic direction of `𝔭` rather than
/// `\bar{𝔭}`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CraterOrientationWitness {
    positive_successor: BTreeMap<IsogenyGraphNodeId, IsogenyGraphNodeId>,
}

impl CraterOrientationWitness {
    /// Validates a user-supplied crater orientation against certified evidence.
    ///
    /// Complexity: `O(V_c + E_c)` set/map construction plus `O(V_c)` cycle
    /// validation, where `V_c` and `E_c` are certified crater nodes and
    /// certified internal horizontal crater edges.
    pub fn new(
        crater: &CraterReport,
        positive_successor: BTreeMap<IsogenyGraphNodeId, IsogenyGraphNodeId>,
    ) -> Result<Self, CraterOrientationWitnessError> {
        Self::validate(crater, &positive_successor)?;
        Ok(Self { positive_successor })
    }

    /// Returns the declared positive successor of `node`.
    pub fn successor(&self, node: IsogenyGraphNodeId) -> Option<IsogenyGraphNodeId> {
        self.positive_successor.get(&node).copied()
    }

    /// Returns whether `node` belongs to this oriented crater cycle.
    pub fn contains_node(&self, node: IsogenyGraphNodeId) -> bool {
        self.positive_successor.contains_key(&node)
    }

    /// Returns the declared predecessor of `node` in the positive orientation.
    pub fn predecessor(&self, node: IsogenyGraphNodeId) -> Option<IsogenyGraphNodeId> {
        self.positive_successor
            .iter()
            .find_map(|(source, target)| (*target == node).then_some(*source))
    }

    /// Returns the closed oriented cycle starting at `start`, if available.
    ///
    /// The returned path repeats `start` as its last node.
    pub fn oriented_cycle_from(
        &self,
        start: IsogenyGraphNodeId,
    ) -> Option<Vec<IsogenyGraphNodeId>> {
        let mut visited = vec![start];
        let mut seen = BTreeSet::from([start]);
        let mut current = start;

        loop {
            let next = self.successor(current)?;
            visited.push(next);

            if next == start {
                return Some(visited);
            }

            if !seen.insert(next) {
                return None;
            }

            current = next;
        }
    }

    /// Builds the opposite orientation and validates it against the same crater.
    pub fn inverse(&self, crater: &CraterReport) -> Result<Self, CraterOrientationWitnessError> {
        let inverse_successor = self
            .positive_successor
            .iter()
            .map(|(source, target)| (*target, *source))
            .collect::<BTreeMap<_, _>>();

        Self::new(crater, inverse_successor)
    }

    fn validate(
        crater: &CraterReport,
        positive_successor: &BTreeMap<IsogenyGraphNodeId, IsogenyGraphNodeId>,
    ) -> Result<(), CraterOrientationWitnessError> {
        if crater.nodes().is_empty() {
            return Err(CraterOrientationWitnessError::EmptyCrater);
        }

        let crater_nodes = crater.nodes().iter().copied().collect::<HashSet<_>>();
        for source in positive_successor.keys().copied() {
            if !crater_nodes.contains(&source) {
                return Err(CraterOrientationWitnessError::SourceOutsideCrater { source });
            }
        }

        let outgoing = crater.certified_internal_outgoing_edge_map();
        for source in crater.nodes().iter().copied() {
            let Some(target) = positive_successor.get(&source).copied() else {
                return Err(CraterOrientationWitnessError::MissingSuccessor { source });
            };
            if !crater_nodes.contains(&target) {
                return Err(CraterOrientationWitnessError::TargetOutsideCrater { source, target });
            }
            let has_certified_edge = outgoing
                .get(&source)
                .into_iter()
                .flatten()
                .any(|edge| edge.target() == target);
            if !has_certified_edge {
                return Err(
                    CraterOrientationWitnessError::MissingCertifiedHorizontalEdge {
                        source,
                        target,
                    },
                );
            }
        }

        Self::validate_single_cycle(crater.nodes()[0], crater.nodes().len(), positive_successor)
    }

    fn validate_single_cycle(
        start: IsogenyGraphNodeId,
        expected_len: usize,
        positive_successor: &BTreeMap<IsogenyGraphNodeId, IsogenyGraphNodeId>,
    ) -> Result<(), CraterOrientationWitnessError> {
        let mut seen = BTreeSet::from([start]);
        let mut current = start;

        for step in 0..expected_len {
            let next = positive_successor
                .get(&current)
                .copied()
                .ok_or(CraterOrientationWitnessError::MissingSuccessor { source: current })?;

            if next == start {
                if step + 1 == expected_len {
                    return Ok(());
                }
                return Err(CraterOrientationWitnessError::DoesNotCloseCycle { start });
            }

            if !seen.insert(next) {
                return Err(CraterOrientationWitnessError::RepeatsBeforeClosing { node: next });
            }

            current = next;
        }

        Err(CraterOrientationWitnessError::DoesNotCloseCycle { start })
    }
}
