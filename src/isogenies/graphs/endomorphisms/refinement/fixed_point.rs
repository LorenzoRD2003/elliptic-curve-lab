use std::collections::BTreeSet;

use crate::elliptic_curves::endomorphisms::quadratic_orders::ImaginaryQuadraticOrder;
use crate::isogenies::graphs::{
    IsogenyGraphNodeId,
    endomorphisms::{
        IsogenyGraphEndomorphismReport,
        refinement::{
            CandidateElimination, CandidateEliminationReason, CandidateRefinementEdgeDirection,
            CandidateRefinementError, CandidateRefinementStrategy, ConstraintSource,
            EndomorphismCandidateRefinement, IsogenyGraphCandidateRefinementReport,
            LocalEndomorphismConstraint, RefinementConfidence,
        },
    },
};

/// Fixed-point refinement outcome plus internal convergence diagnostics.
pub(crate) struct FixedPointCandidateRefinement {
    report: IsogenyGraphCandidateRefinementReport,
    rounds_with_eliminations: usize,
}

impl FixedPointCandidateRefinement {
    #[cfg(test)]
    pub(crate) fn report(&self) -> &IsogenyGraphCandidateRefinementReport {
        &self.report
    }

    pub(crate) fn into_report(self) -> IsogenyGraphCandidateRefinementReport {
        let Self {
            report,
            rounds_with_eliminations,
        } = self;
        let _ = rounds_with_eliminations;
        report
    }

    #[cfg(test)]
    pub(crate) fn rounds_with_eliminations(&self) -> usize {
        self.rounds_with_eliminations
    }
}

impl IsogenyGraphEndomorphismReport {
    /// Refines candidate endomorphism orders globally until the survivor sets
    /// reach a fixed point.
    ///
    /// This is a monotone propagation pass. Each round computes eliminations
    /// from the survivor sets of the previous round and applies all removals
    /// simultaneously, so the iteration order of nodes and edges does not
    /// affect the result. A node may end with no survivors; that represents
    /// incompatibility of the observed tentative evidence under the selected
    /// strategy, not a construction error.
    ///
    /// The result still does not certify `End(E)`.
    pub fn refine_candidates_to_fixed_point(
        &self,
        strategy: CandidateRefinementStrategy,
    ) -> Result<IsogenyGraphCandidateRefinementReport, CandidateRefinementError> {
        Ok(self
            .fixed_point_candidate_refinement(strategy)?
            .into_report())
    }

    pub(crate) fn fixed_point_candidate_refinement(
        &self,
        strategy: CandidateRefinementStrategy,
    ) -> Result<FixedPointCandidateRefinement, CandidateRefinementError> {
        if strategy == CandidateRefinementStrategy::NodeLocalLevelsOnly {
            return Ok(FixedPointCandidateRefinement {
                report: self.refine_candidates(strategy)?,
                rounds_with_eliminations: 0,
            });
        }

        let node_levels = self.node_candidate_levels()?;
        let mut survivors = node_levels
            .iter()
            .map(|levels| (0..levels.len()).collect::<BTreeSet<_>>())
            .collect::<Vec<_>>();
        let mut elimination_reasons = node_levels
            .iter()
            .map(|levels| vec![None; levels.len()])
            .collect::<Vec<Vec<Option<CandidateEliminationReason>>>>();
        let mut rounds_with_eliminations = 0;

        loop {
            let previous = survivors.clone();
            let mut next = previous.clone();
            let mut changed = false;

            for node in self.nodes() {
                let node_index = node.node_id().0;
                for &candidate_index in &previous[node_index] {
                    let candidate_level = node_levels[node_index][candidate_index];
                    if let Some(reason) = self.fixed_point_elimination_reason(
                        node.node_id(),
                        candidate_level,
                        &previous,
                        &node_levels,
                    ) {
                        next[node_index].remove(&candidate_index);
                        if elimination_reasons[node_index][candidate_index].is_none() {
                            elimination_reasons[node_index][candidate_index] = Some(reason);
                        }
                        changed = true;
                    }
                }
            }

            if !changed {
                break;
            }

            rounds_with_eliminations += 1;
            survivors = next;
        }

        Ok(FixedPointCandidateRefinement {
            report: self.fixed_point_report_from_survivors(
                strategy,
                survivors,
                elimination_reasons,
            )?,
            rounds_with_eliminations,
        })
    }

    fn node_candidate_levels(&self) -> Result<Vec<Vec<u32>>, CandidateRefinementError> {
        self.nodes()
            .iter()
            .map(|node| {
                node.candidate_set()
                    .candidate_orders()
                    .iter()
                    .map(|candidate| {
                        candidate
                            .volcanic_level_at(self.prime())
                            .map_err(|_| CandidateRefinementError::InvalidLocalPrime)
                            .map(|level| level.level())
                    })
                    .collect()
            })
            .collect()
    }

    fn fixed_point_elimination_reason(
        &self,
        node_id: IsogenyGraphNodeId,
        candidate_level: u32,
        previous_survivors: &[BTreeSet<usize>],
        node_levels: &[Vec<u32>],
    ) -> Option<CandidateEliminationReason> {
        let allowed_levels = self.node_report(node_id)?.refinement_allowed_levels();

        if !allowed_levels.contains(&candidate_level) {
            return Some(CandidateEliminationReason::IncompatibleLocalLevel {
                ell: self.prime().clone(),
                candidate_level,
                allowed_levels,
            });
        }

        for edge in self.edges() {
            let (direction, adjacent_node) = if edge.source() == node_id {
                (CandidateRefinementEdgeDirection::Outgoing, edge.target())
            } else if edge.target() == node_id {
                (CandidateRefinementEdgeDirection::Incoming, edge.source())
            } else {
                continue;
            };

            let Some(relation) = edge.refinement_relation() else {
                continue;
            };

            let adjacent_levels = previous_survivors[adjacent_node.0]
                .iter()
                .map(|&index| node_levels[adjacent_node.0][index])
                .collect::<BTreeSet<_>>();

            let compatible = adjacent_levels.iter().any(|&adjacent_level| {
                let (source_level, target_level) = match direction {
                    CandidateRefinementEdgeDirection::Outgoing => (candidate_level, adjacent_level),
                    CandidateRefinementEdgeDirection::Incoming => (adjacent_level, candidate_level),
                };
                relation.allows_levels(source_level, target_level)
            });

            if !compatible {
                return Some(
                    CandidateEliminationReason::IncompatibleIncidentEdgeRelation {
                        ell: self.prime().clone(),
                        direction,
                        adjacent_node,
                        candidate_level,
                        compatible_adjacent_levels: adjacent_levels,
                        expected_relation: relation,
                    },
                );
            }
        }

        None
    }

    fn fixed_point_report_from_survivors(
        &self,
        strategy: CandidateRefinementStrategy,
        survivors: Vec<BTreeSet<usize>>,
        elimination_reasons: Vec<Vec<Option<CandidateEliminationReason>>>,
    ) -> Result<IsogenyGraphCandidateRefinementReport, CandidateRefinementError> {
        let mut node_refinements = Vec::with_capacity(self.nodes().len());

        for node in self.nodes() {
            let node_index = node.node_id().0;
            let candidate_orders = node.candidate_set().candidate_orders();
            let surviving_candidates = survivors[node_index]
                .iter()
                .map(|&index| candidate_orders[index].clone())
                .collect::<Vec<ImaginaryQuadraticOrder>>();
            let eliminated_candidates = candidate_orders
                .iter()
                .enumerate()
                .filter(|(index, _)| !survivors[node_index].contains(index))
                .map(|(index, candidate)| {
                    let reason = elimination_reasons[node_index][index]
                        .clone()
                        .expect("every eliminated candidate should record its first reason");
                    CandidateElimination::new(candidate.clone(), reason)
                })
                .collect();

            node_refinements.push(EndomorphismCandidateRefinement::new(
                node.node_id(),
                node.candidate_set().clone(),
                surviving_candidates,
                eliminated_candidates,
                self.fixed_point_constraints_for_node(node.node_id())?,
                self.fixed_point_confidence_for_node(node.node_id()),
            ));
        }

        Ok(IsogenyGraphCandidateRefinementReport::new(
            self.prime().clone(),
            strategy,
            node_refinements,
        ))
    }

    fn fixed_point_constraints_for_node(
        &self,
        node_id: IsogenyGraphNodeId,
    ) -> Result<Vec<LocalEndomorphismConstraint>, CandidateRefinementError> {
        let node = self
            .node_report(node_id)
            .ok_or(CandidateRefinementError::NodeNotFound { node_id })?;
        let mut constraints = vec![LocalEndomorphismConstraint::NodeLevel {
            ell: self.prime().clone(),
            allowed_levels: node.refinement_allowed_levels(),
            provenance: ConstraintSource::NodeReport { node_id },
        }];

        for edge in self.edges() {
            if edge.source() != node_id && edge.target() != node_id {
                continue;
            }
            let Some(relation) = edge.refinement_relation() else {
                continue;
            };
            constraints.push(LocalEndomorphismConstraint::EdgeRelation {
                ell: self.prime().clone(),
                source_node: edge.source(),
                target_node: edge.target(),
                relation,
                provenance: ConstraintSource::EdgeReport {
                    edge_id: edge.edge_id(),
                },
            });
        }

        Ok(constraints)
    }

    fn fixed_point_confidence_for_node(&self, node_id: IsogenyGraphNodeId) -> RefinementConfidence {
        if self.edges().iter().any(|edge| {
            (edge.source() == node_id || edge.target() == node_id)
                && edge.refinement_relation().is_some()
        }) {
            RefinementConfidence::PropagatedTentativeGraphEvidence
        } else {
            RefinementConfidence::ConservativeLocalEvidence
        }
    }
}
