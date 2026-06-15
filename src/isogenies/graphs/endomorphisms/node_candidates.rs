use std::hash::Hash;

use crate::elliptic_curves::{
    endomorphisms::{
        candidate_sets::EndomorphismRingCandidateSet,
        quadratic_orders::{
            ImaginaryQuadraticOrderError, QuadraticDiscriminantFactorization,
            QuadraticDiscriminantFactorizationError,
        },
    },
    frobenius::FrobeniusTraceCurveModel,
};
use crate::fields::{traits::EnumerableFiniteField, traits::FiniteField, traits::SqrtField};
use crate::isogenies::graphs::{
    GraphCurveModel, IsogenyGraph, IsogenyGraphError, IsogenyGraphNode, IsogenyGraphNodeId,
};

impl<C: GraphCurveModel + FrobeniusTraceCurveModel> IsogenyGraphNode<C>
where
    C::BaseField: EnumerableFiniteField<Elem = C::Elem> + SqrtField<Elem = C::Elem> + FiniteField,
    C::Point: Clone + PartialEq,
{
    /// Derives the Frobenius-compatible candidate endomorphism orders for the stored node.
    ///
    /// This uses the representative curve already stored on the node:
    ///
    /// `representative -> FrobeniusTrace -> Δ_π -> candidate quadratic orders`.
    ///
    /// The result is still only a Frobenius-compatible candidate set. It does
    /// not certify the exact geometric endomorphism ring of the node's curve.
    ///
    /// Complexity: this is the cost of [`FrobeniusTraceCurveModel::frobenius_trace`],
    /// and then dominated by `num-prime`.
    pub fn endomorphism_ring_candidates(
        &self,
    ) -> Result<EndomorphismRingCandidateSet, IsogenyGraphError> {
        let discriminant = self.representative().frobenius_trace()?.discriminant();
        let factorization =
            QuadraticDiscriminantFactorization::from_frobenius_discriminant(&discriminant)
                .map_err(IsogenyGraphError::from)?;

        factorization
            .endomorphism_ring_candidates()
            .map_err(IsogenyGraphError::from)
    }
}

impl<C: GraphCurveModel + FrobeniusTraceCurveModel> IsogenyGraph<C>
where
    C::BaseField: EnumerableFiniteField<Elem = C::Elem> + SqrtField<Elem = C::Elem> + FiniteField,
    C::Point: Clone + Eq + Hash + PartialEq,
    C::IsomorphismWitness: Clone + core::fmt::Debug,
{
    /// Derives the Frobenius-compatible candidate endomorphism orders for one stored node.
    ///
    /// Complexity: this is the cost of [`FrobeniusTraceCurveModel::frobenius_trace`],
    /// and then dominated by `num-prime`.
    pub fn node_endomorphism_candidates(
        &self,
        node_id: IsogenyGraphNodeId,
    ) -> Result<EndomorphismRingCandidateSet, IsogenyGraphError> {
        let node = self
            .node(node_id)
            .ok_or(IsogenyGraphError::MissingSourceNode(node_id))?;
        node.endomorphism_ring_candidates()
    }

    /// Derives the Frobenius-compatible candidate endomorphism orders for every stored node.
    ///
    /// The result is returned in dense node-id order, matching the graph's node
    /// storage convention.
    ///
    /// Complexity: one exhaustive Frobenius-trace computation per node, plus
    /// arithmetic dominated by `num-prime` for each node.
    pub fn graph_endomorphism_candidates(
        &self,
    ) -> Result<Vec<(IsogenyGraphNodeId, EndomorphismRingCandidateSet)>, IsogenyGraphError> {
        self.nodes()
            .iter()
            .map(|node| {
                node.endomorphism_ring_candidates()
                    .map(|candidate_set| (node.id(), candidate_set))
            })
            .collect()
    }
}

impl From<ImaginaryQuadraticOrderError> for IsogenyGraphError {
    fn from(error: ImaginaryQuadraticOrderError) -> Self {
        Self::Endomorphism(error)
    }
}

impl From<QuadraticDiscriminantFactorizationError> for IsogenyGraphError {
    fn from(error: QuadraticDiscriminantFactorizationError) -> Self {
        Self::EndomorphismFactorization(error)
    }
}

#[cfg(test)]
mod tests {
    use crate::elliptic_curves::ShortWeierstrassCurve;
    use crate::fields::{Fp, traits::Field};
    use crate::isogenies::graphs::{IsogenyGraphBuilder, IsogenyGraphError, IsogenyGraphNodeId};

    type F41 = Fp<41>;
    type Curve41 = ShortWeierstrassCurve<F41>;

    fn f41_curve() -> Curve41 {
        Curve41::new(F41::from_i64(2), F41::from_i64(3)).expect("valid curve")
    }

    #[test]
    fn node_can_derive_its_candidate_endomorphism_orders() {
        let graph = IsogenyGraphBuilder::new(f41_curve(), 2)
            .max_depth(0)
            .build()
            .expect("depth-zero graph should build");
        let node = graph
            .node(IsogenyGraphNodeId(0))
            .expect("root node should exist");

        let candidate_set = node
            .endomorphism_ring_candidates()
            .expect("candidate endomorphism orders should derive");

        assert!(!candidate_set.is_empty());
        assert_eq!(
            candidate_set.maximal_order().conductor(),
            &num_bigint::BigUint::from(1u8)
        );
    }

    #[test]
    fn graph_can_derive_candidates_for_one_requested_node() {
        let graph = IsogenyGraphBuilder::new(f41_curve(), 2)
            .max_depth(1)
            .build()
            .expect("depth-one graph should build");

        let candidate_set = graph
            .node_endomorphism_candidates(IsogenyGraphNodeId(0))
            .expect("root candidate set should derive");

        assert!(!candidate_set.is_empty());
    }

    #[test]
    fn graph_can_derive_candidates_for_every_stored_node() {
        let graph = IsogenyGraphBuilder::new(f41_curve(), 2)
            .max_depth(1)
            .build()
            .expect("depth-one graph should build");

        let candidate_sets = graph
            .graph_endomorphism_candidates()
            .expect("all node candidate sets should derive");

        assert_eq!(candidate_sets.len(), graph.node_count());
        assert_eq!(candidate_sets[0].0, IsogenyGraphNodeId(0));
        assert!(
            candidate_sets
                .iter()
                .all(|(_, candidate_set)| !candidate_set.is_empty())
        );
    }

    #[test]
    fn missing_node_reports_a_graph_error() {
        let graph = IsogenyGraphBuilder::new(f41_curve(), 2)
            .max_depth(0)
            .build()
            .expect("depth-zero graph should build");

        assert_eq!(
            graph.node_endomorphism_candidates(IsogenyGraphNodeId(99)),
            Err(IsogenyGraphError::MissingSourceNode(IsogenyGraphNodeId(99)))
        );
    }
}
