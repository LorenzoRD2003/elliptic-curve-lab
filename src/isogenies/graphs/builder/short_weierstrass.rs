use num_prime::nt_funcs::is_prime;
use std::collections::VecDeque;
use std::hash::Hash;

use crate::elliptic_curves::{
    ShortWeierstrassCurve, short_weierstrass::isomorphisms::ShortWeierstrassIsomorphism,
};
use crate::fields::traits::{EnumerableFiniteField, Field, SqrtField};
use crate::isogenies::graphs::edge::EdgeTargetWitness;
use crate::isogenies::{
    graphs::{
        GraphTorsionCurveModel, IsogenyGraphEdge, IsogenyGraphEdgeId, IsogenyGraphError,
        IsogenyGraphNode, IsogenyGraphNodeId,
    },
    traits::Isogeny,
};

use crate::isogenies::graphs::builder::{IsogenyGraph, IsogenyGraphBuilder};

struct ResolvedTarget<I> {
    id: IsogenyGraphNodeId,
    witness: EdgeTargetWitness<I>,
    inserted_new_node: bool,
}

fn degree_is_prime(degree: usize) -> bool {
    is_prime(&(degree as u64), None).probably()
}

impl<F> IsogenyGraph<ShortWeierstrassCurve<F>>
where
    F: Field + EnumerableFiniteField + SqrtField + Clone,
    F::Elem: Clone + Eq + Hash + PartialEq,
{
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

    #[cfg(test)]
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

    #[cfg(test)]
    pub(crate) fn outgoing_velu_edges_from_node(
        &mut self,
        source: IsogenyGraphNodeId,
        ell: usize,
    ) -> Result<Vec<IsogenyGraphEdgeId>, IsogenyGraphError> {
        if ell < 2 {
            return Err(IsogenyGraphError::InvalidDegree);
        }

        if !degree_is_prime(ell) {
            return Err(IsogenyGraphError::DegreeMustBePrimeForThisBuilder { degree: ell });
        }

        let source_curve = self
            .node(source)
            .ok_or(IsogenyGraphError::MissingSourceNode(source))?
            .representative()
            .clone();

        let kernels = source_curve.cyclic_kernels_of_order(ell)?;
        if kernels.is_empty() {
            return Err(IsogenyGraphError::NonRationalKernelForRequestedDegree { degree: ell });
        }

        let mut new_edge_ids = Vec::with_capacity(kernels.len());

        for kernel in kernels {
            let generator =
                kernel.points().get(1).cloned().ok_or(
                    IsogenyGraphError::NonRationalKernelForRequestedDegree { degree: ell },
                )?;
            let velu = source_curve.velu_isogeny_from_generator(generator)?;
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
    F: Field + EnumerableFiniteField + SqrtField + Clone,
    F::Elem: Clone + Eq + Hash + PartialEq,
{
    pub fn build(self) -> Result<IsogenyGraph<ShortWeierstrassCurve<F>>, IsogenyGraphError> {
        if self.ell < 2 {
            return Err(IsogenyGraphError::InvalidDegree);
        }

        if !degree_is_prime(self.ell) {
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
            let kernels = source_curve.cyclic_kernels_of_order(self.ell)?;

            for kernel in kernels {
                let generator = kernel.points().get(1).cloned().ok_or(
                    IsogenyGraphError::NonRationalKernelForRequestedDegree { degree: self.ell },
                )?;
                let velu = source_curve.velu_isogeny_from_generator(generator)?;
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
