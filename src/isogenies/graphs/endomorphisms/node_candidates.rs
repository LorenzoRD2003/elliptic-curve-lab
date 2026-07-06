use std::hash::Hash;

use crate::elliptic_curves::{
    endomorphisms::{
        candidate_sets::EndomorphismRingCandidateSet,
        quadratic_orders::{
            ImaginaryQuadraticOrderError, QuadraticDiscriminantFactorization,
            QuadraticDiscriminantFactorizationError,
        },
    },
    traits::FrobeniusTraceCurveModel,
};
use crate::fields::traits::{EnumerableFiniteField, FiniteField, SqrtField};
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
mod tests;
