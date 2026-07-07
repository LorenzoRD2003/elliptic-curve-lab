use std::collections::HashSet;
use std::fmt;
use std::hash::Hash;

use num_bigint::BigUint;

use crate::isogenies::graphs::{
    GraphCurveModel, IsogenyGraph, IsogenyGraphNodeId, VolcanoRole,
    endomorphisms::{
        VolcanoSearchError,
        volcano_structure::{
            CraterReport, CraterShape, HorizontalEdgeReport, HorizontalEdgeStatus,
        },
    },
};

impl<C: GraphCurveModel> IsogenyGraph<C>
where
    C::Point: Clone + Eq + Hash,
    C::IsomorphismWitness: Clone + fmt::Debug,
{
    /// Builds a crater report from certified altitude and weak surface evidence.
    ///
    /// The report first builds [`Self::volcano_structure_report`]. Its certified
    /// crater is the surface `V₀` from that report. Stored edges between
    /// certified crater nodes are horizontal by altitude. When altitude is not
    /// available, weak surface-like evidence may mark an edge as suspected, and
    /// partial endpoints are reported as not certifiable.
    ///
    /// Complexity: one volcano-structure pass, plus linear scans over stored
    /// nodes and edges.
    pub fn volcano_crater_report(
        &self,
        prime: &BigUint,
    ) -> Result<CraterReport, VolcanoSearchError> {
        let structure = self.volcano_structure_report(prime)?;
        let crater_nodes = structure
            .crater()
            .map(|level| level.nodes().to_vec())
            .unwrap_or_default();
        let crater_node_set = crater_nodes.iter().copied().collect::<HashSet<_>>();
        let weak_surface_nodes = self.weak_surface_nodes();
        let partial_nodes = structure
            .partial_uncertified_node_ids()
            .into_iter()
            .collect::<HashSet<_>>();

        let horizontal_edges = self
            .edges()
            .iter()
            .filter_map(|edge| {
                let status = HorizontalEdgeStatus::from_evidence(
                    edge.source(),
                    edge.target(),
                    &crater_node_set,
                    &weak_surface_nodes,
                    &partial_nodes,
                )?;
                Some(HorizontalEdgeReport::new(
                    edge.id(),
                    edge.source(),
                    edge.target(),
                    status,
                ))
            })
            .collect::<Vec<_>>();

        let shape = CraterShape::from_crater_evidence(&crater_nodes, &horizontal_edges);

        Ok(CraterReport::new(
            prime.clone(),
            structure,
            crater_nodes,
            horizontal_edges,
            shape,
        ))
    }

    fn weak_surface_nodes(&self) -> HashSet<IsogenyGraphNodeId> {
        if self.node(IsogenyGraphNodeId(0)).is_none() {
            return HashSet::new();
        }

        self.infer_volcano_like_layers(IsogenyGraphNodeId(0))
            .roles()
            .iter()
            .filter_map(|(node_id, role)| (*role == VolcanoRole::Surface).then_some(*node_id))
            .collect()
    }
}
