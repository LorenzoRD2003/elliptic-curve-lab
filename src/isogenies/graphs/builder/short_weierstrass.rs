use std::collections::VecDeque;
use std::hash::Hash;

use super::graph::{IsogenyGraph, IsogenyGraphBuilder};
use crate::elliptic_curves::{ShortWeierstrassCurve, ShortWeierstrassIsomorphism};
use crate::fields::{EnumerableFiniteField, SqrtField};
use crate::isogenies::graphs::{
    EdgeTargetWitness, IsogenyGraphEdge, IsogenyGraphEdgeId, IsogenyGraphError, IsogenyGraphNode,
    IsogenyGraphNodeId, cyclic_kernels_of_order,
};
use crate::isogenies::{Isogeny, VeluIsogeny};

struct ResolvedTarget<I> {
    id: IsogenyGraphNodeId,
    witness: EdgeTargetWitness<I>,
    inserted_new_node: bool,
}

impl<F> IsogenyGraph<ShortWeierstrassCurve<F>>
where
    F: EnumerableFiniteField + SqrtField + Clone,
    F::Elem: Clone + Eq + Hash + PartialEq,
{
    /// Searches the current node set for a stored short-Weierstrass
    /// representative that is base-field isomorphic to `curve`.
    ///
    /// - first, filter by equal `j`-invariant as a cheap necessary condition
    /// - then, run the existing exhaustive base-field witness search on the
    ///   surviving short-Weierstrass representatives
    pub(super) fn find_isomorphic_node(
        nodes: &[IsogenyGraphNode<ShortWeierstrassCurve<F>>],
        curve: &ShortWeierstrassCurve<F>,
    ) -> Option<(IsogenyGraphNodeId, ShortWeierstrassIsomorphism<F>)> {
        let target_j = curve.j_invariant();
        for node in nodes {
            if !F::eq(&node.j_invariant(), &target_j) {
                continue;
            }
            if let Some(witness) = curve.find_isomorphism_to(node.representative()) {
                return Some((node.id(), witness));
            }
        }
        None
    }

    /// Resolves a raw Vélu codomain onto one stored target-node representative.
    pub(super) fn resolve_target_node(
        &mut self,
        raw_codomain: ShortWeierstrassCurve<F>,
    ) -> (
        IsogenyGraphNodeId,
        EdgeTargetWitness<ShortWeierstrassIsomorphism<F>>,
    ) {
        let resolved = self.resolve_target_node_with_options(raw_codomain, true);
        (resolved.id, resolved.witness)
    }

    fn resolve_target_node_with_options(
        &mut self,
        raw_codomain: ShortWeierstrassCurve<F>,
        deduplicate_by_base_field_isomorphism: bool,
    ) -> ResolvedTarget<ShortWeierstrassIsomorphism<F>> {
        if let Some(target_id) = self
            .nodes
            .iter()
            .find(|node| node.representative() == &raw_codomain)
            .map(IsogenyGraphNode::id)
        {
            return ResolvedTarget {
                id: target_id,
                witness: EdgeTargetWitness::Identity,
                inserted_new_node: false,
            };
        }

        if deduplicate_by_base_field_isomorphism
            && let Some((target_id, witness)) =
                Self::find_isomorphic_node(&self.nodes, &raw_codomain)
        {
            return ResolvedTarget {
                id: target_id,
                witness: EdgeTargetWitness::Explicit(witness),
                inserted_new_node: false,
            };
        }

        let target_id = IsogenyGraphNodeId(self.nodes.len());
        self.nodes
            .push(IsogenyGraphNode::new(target_id, raw_codomain));
        ResolvedTarget {
            id: target_id,
            witness: EdgeTargetWitness::Identity,
            inserted_new_node: true,
        }
    }

    /// Appends all outgoing degree-`ell` Vélu edges from the requested source node.
    ///
    /// - `ell` must be prime
    /// - kernels must be rational cyclic subgroups of exact order `ell`
    /// - target nodes are deduplicated up to base-field isomorphism
    ///
    /// The returned ids are exactly the newly stored edge ids, in append order.
    pub fn outgoing_velu_edges_from_node(
        &mut self,
        source: IsogenyGraphNodeId,
        ell: usize,
    ) -> Result<Vec<IsogenyGraphEdgeId>, IsogenyGraphError> {
        if ell < 2 {
            return Err(IsogenyGraphError::InvalidDegree);
        }

        if !is_prime(ell) {
            return Err(IsogenyGraphError::DegreeMustBePrimeForThisBuilder { degree: ell });
        }

        let source_curve = self
            .node(source)
            .ok_or(IsogenyGraphError::MissingSourceNode(source))?
            .representative()
            .clone();

        let kernels = cyclic_kernels_of_order(&source_curve, ell)?;
        if kernels.is_empty() {
            return Err(IsogenyGraphError::NonRationalKernelForRequestedDegree { degree: ell });
        }

        let mut new_edge_ids = Vec::with_capacity(kernels.len());

        for kernel in kernels {
            let generator =
                kernel.points().get(1).cloned().ok_or(
                    IsogenyGraphError::NonRationalKernelForRequestedDegree { degree: ell },
                )?;
            let velu = VeluIsogeny::from_generator(source_curve.clone(), generator)?;
            let (target, target_witness) = self.resolve_target_node(velu.codomain().clone());

            let edge_id = IsogenyGraphEdgeId(self.edges.len());
            self.edges.push(IsogenyGraphEdge::new(
                edge_id,
                source,
                target,
                velu.kernel().clone(),
                target_witness,
            ));
            new_edge_ids.push(edge_id);
        }

        Ok(new_edge_ids)
    }
}

impl<F> IsogenyGraphBuilder<ShortWeierstrassCurve<F>>
where
    F: EnumerableFiniteField + SqrtField + Clone,
    F::Elem: Clone + Eq + Hash + PartialEq,
{
    /// Builds a small breadth-first `l`-isogeny graph from the configured start curve.
    /// `l` must be prime.
    pub fn build(self) -> Result<IsogenyGraph<ShortWeierstrassCurve<F>>, IsogenyGraphError> {
        if self.ell < 2 {
            return Err(IsogenyGraphError::InvalidDegree);
        }

        if !is_prime(self.ell) {
            return Err(IsogenyGraphError::DegreeMustBePrimeForThisBuilder { degree: self.ell });
        }

        let mut graph = IsogenyGraph {
            nodes: vec![IsogenyGraphNode::new(
                IsogenyGraphNodeId(0),
                self.start_curve,
            )],
            edges: vec![],
        };
        let mut queue = VecDeque::from([(IsogenyGraphNodeId(0), 0usize)]);

        while let Some((source, depth)) = queue.pop_front() {
            if depth >= self.max_depth {
                continue;
            }

            let source_curve = graph
                .node(source)
                .expect("queued node ids should stay valid")
                .representative()
                .clone();
            let kernels = cyclic_kernels_of_order(&source_curve, self.ell)?;

            for kernel in kernels {
                let generator = kernel.points().get(1).cloned().ok_or(
                    IsogenyGraphError::NonRationalKernelForRequestedDegree { degree: self.ell },
                )?;
                let velu = VeluIsogeny::from_generator(source_curve.clone(), generator)?;
                let resolved = graph.resolve_target_node_with_options(
                    velu.codomain().clone(),
                    self.deduplicate_by_base_field_isomorphism,
                );

                let edge_id = IsogenyGraphEdgeId(graph.edges.len());
                graph.edges.push(IsogenyGraphEdge::new(
                    edge_id,
                    source,
                    resolved.id,
                    velu.kernel().clone(),
                    resolved.witness,
                ));

                if resolved.inserted_new_node {
                    queue.push_back((resolved.id, depth + 1));
                }
            }
        }

        Ok(graph)
    }
}

fn is_prime(n: usize) -> bool {
    if n < 2 {
        return false;
    }
    if n == 2 {
        return true;
    }
    if n.is_multiple_of(2) {
        return false;
    }

    let mut divisor = 3;
    while divisor * divisor <= n {
        if n.is_multiple_of(divisor) {
            return false;
        }
        divisor += 2;
    }
    true
}
