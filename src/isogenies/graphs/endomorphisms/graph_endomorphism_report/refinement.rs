use crate::isogenies::graphs::{
    IsogenyGraphNodeId,
    endomorphisms::{
        graph_endomorphism_report::IsogenyGraphEndomorphismReport,
        refinement::{
            CandidateRefinementError, CandidateRefinementStrategy, EndomorphismCandidateRefinement,
            IncidentEdgeRefinementConstraint, IsogenyGraphCandidateRefinementReport,
        },
    },
};

impl IsogenyGraphEndomorphismReport {
    /// Refines the candidate endomorphism orders for one node using evidence
    /// already recorded in this graph report.
    ///
    /// [`CandidateRefinementStrategy::NodeLocalLevelsOnly`] uses only the
    /// node-local conductor-level constraint `v_ℓ(f) ∈ L`, where `L` may be
    /// narrowed by endpoint-role evidence observed in the graph.
    /// [`CandidateRefinementStrategy::Conservative`] and
    /// [`CandidateRefinementStrategy::IncidentUnambiguousEdges`] also use
    /// one-hop incident edge constraints, but only when either the observed
    /// graph relation or the arithmetic candidate-set relation is unequivocal.
    /// `Ambiguous` and `Unsupported` edge evidence is ignored by these
    /// conservative strategies.
    ///
    /// The result is not a certification of `End(E)`: even a unique survivor is
    /// only the unique candidate compatible with the evidence used by the
    /// selected strategy.
    pub fn refine_candidates_for_node(
        &self,
        node_id: IsogenyGraphNodeId,
        strategy: CandidateRefinementStrategy,
    ) -> Result<EndomorphismCandidateRefinement, CandidateRefinementError> {
        let node_report = self
            .node_report(node_id)
            .ok_or(CandidateRefinementError::NodeNotFound { node_id })?;

        match strategy {
            CandidateRefinementStrategy::NodeLocalLevelsOnly => {
                EndomorphismCandidateRefinement::from_node_local_levels(node_report, self.prime())
            }
            CandidateRefinementStrategy::Conservative
            | CandidateRefinementStrategy::IncidentUnambiguousEdges => {
                EndomorphismCandidateRefinement::from_incident_unambiguous_edges(
                    node_report,
                    self.prime(),
                    self.incident_edge_constraints_for_node(node_id)?,
                )
            }
        }
    }

    /// Refines candidate endomorphism orders independently for every stored node.
    ///
    /// This aggregate pass is intentionally not a fixed-point propagation
    /// algorithm. Each node refinement uses the evidence already present in
    /// this report, and incident-edge checks compare against the neighbor's
    /// original allowed levels from the report rather than against recursively
    /// refined survivor sets.
    pub fn refine_candidates(
        &self,
        strategy: CandidateRefinementStrategy,
    ) -> Result<IsogenyGraphCandidateRefinementReport, CandidateRefinementError> {
        let node_refinements = self
            .nodes()
            .iter()
            .map(|node| self.refine_candidates_for_node(node.node_id(), strategy))
            .collect::<Result<Vec<_>, _>>()?;

        Ok(IsogenyGraphCandidateRefinementReport::new(
            self.prime().clone(),
            strategy,
            node_refinements,
        ))
    }

    fn incident_edge_constraints_for_node(
        &self,
        node_id: IsogenyGraphNodeId,
    ) -> Result<Vec<IncidentEdgeRefinementConstraint>, CandidateRefinementError> {
        let mut constraints = Vec::new();

        for edge in self.edges() {
            let adjacent_node = if edge.source() == node_id {
                edge.target()
            } else if edge.target() == node_id {
                edge.source()
            } else {
                continue;
            };

            let adjacent_levels = self
                .node_report(adjacent_node)
                .ok_or(CandidateRefinementError::NodeNotFound {
                    node_id: adjacent_node,
                })?
                .refinement_allowed_levels();

            if let Some(constraint) =
                IncidentEdgeRefinementConstraint::from_edge_report(node_id, edge, adjacent_levels)
            {
                constraints.push(constraint);
            }
        }

        Ok(constraints)
    }
}
