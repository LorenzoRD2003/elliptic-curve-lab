use std::hash::Hash;

use super::GraphCurveModel;
use super::node::IsogenyGraphNodeId;
use crate::isogenies::IsogenyKernel;

/// Stable identifier for one stored graph edge.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct IsogenyGraphEdgeId(pub usize);

/// How an edge's raw codomain relates to the chosen target-node representative.
///
/// When the target representative is exactly the raw codomain curve produced by
/// the edge construction, [`Self::Identity`] is enough. Otherwise, the graph
/// stores an explicit same-family isomorphism witness carrying that codomain to
/// the chosen representative.
#[derive(Clone, Debug)]
pub enum EdgeTargetWitness<I> {
    Identity,
    Explicit(I),
}

/// One explicit directed edge in the educational isogeny graph. The edge stores:
///
/// - source and target node identifiers
/// - one explicit rational kernel subgroup
/// - an optional witness transporting the raw codomain onto the target
///   representative
///
/// It deliberately does not store an entire `VeluIsogeny<C>` object. The goal
/// is to keep the graph light and explicit while preserving the mathematically
/// meaningful data needed to reconstruct later educational summaries.
#[derive(Clone, Debug)]
pub struct IsogenyGraphEdge<C: GraphCurveModel>
where
    C::Point: Clone + Eq + Hash,
{
    id: IsogenyGraphEdgeId,
    source: IsogenyGraphNodeId,
    target: IsogenyGraphNodeId,
    kernel: IsogenyKernel<C>,
    target_witness: EdgeTargetWitness<C::IsomorphismWitness>,
}

impl<C: GraphCurveModel> IsogenyGraphEdge<C>
where
    C::Point: Clone + Eq + Hash,
{
    /// Builds one explicit graph edge from kernel data and target transport.
    pub fn new(
        id: IsogenyGraphEdgeId,
        source: IsogenyGraphNodeId,
        target: IsogenyGraphNodeId,
        kernel: IsogenyKernel<C>,
        target_witness: EdgeTargetWitness<C::IsomorphismWitness>,
    ) -> Self {
        Self {
            id,
            source,
            target,
            kernel,
            target_witness,
        }
    }

    /// Returns the edge identifier.
    pub fn id(&self) -> IsogenyGraphEdgeId {
        self.id
    }

    /// Returns the source node identifier.
    pub fn source(&self) -> IsogenyGraphNodeId {
        self.source
    }

    /// Returns the target node identifier.
    pub fn target(&self) -> IsogenyGraphNodeId {
        self.target
    }

    /// Returns the validated rational kernel subgroup stored on the edge.
    pub fn kernel(&self) -> &IsogenyKernel<C> {
        &self.kernel
    }

    /// Returns the separable degree inferred from the stored kernel size.
    pub fn degree(&self) -> usize {
        self.kernel.degree()
    }

    /// Returns the explicit transport from the raw codomain to the target
    /// representative.
    pub fn target_witness(&self) -> &EdgeTargetWitness<C::IsomorphismWitness> {
        &self.target_witness
    }
}

#[cfg(test)]
mod tests {
    use std::collections::HashSet;

    use super::{EdgeTargetWitness, IsogenyGraphEdge, IsogenyGraphEdgeId};
    use crate::elliptic_curves::{
        AffineCurveModel, CurveIsomorphism, CurveModel, ShortWeierstrassCurve,
        ShortWeierstrassIsomorphism,
    };
    use crate::fields::{Field, Fp};
    use crate::isogenies::IsogenyKernel;
    use crate::isogenies::graphs::IsogenyGraphNodeId;

    type F41 = Fp<41>;

    /// y^2 = x^3 + 2x + 3
    fn f41_curve() -> ShortWeierstrassCurve<F41> {
        ShortWeierstrassCurve::<F41>::new(F41::from_i64(2), F41::from_i64(3)).expect("valid curve")
    }

    fn sample_kernel(
        curve: &ShortWeierstrassCurve<F41>,
    ) -> IsogenyKernel<ShortWeierstrassCurve<F41>> {
        let two_torsion = curve
            .point(F41::from_i64(40), F41::from_i64(0))
            .expect("sample point should lie on the curve");

        IsogenyKernel::new(curve, HashSet::from([curve.identity(), two_torsion]))
            .expect("two-torsion subgroup should be valid")
    }

    #[test]
    fn edge_exposes_ids_kernel_degree_and_identity_witness() {
        let curve = f41_curve();
        let kernel = sample_kernel(&curve);
        let edge = IsogenyGraphEdge::new(
            IsogenyGraphEdgeId(5),
            IsogenyGraphNodeId(1),
            IsogenyGraphNodeId(2),
            kernel,
            EdgeTargetWitness::Identity,
        );

        assert_eq!(edge.id(), IsogenyGraphEdgeId(5));
        assert_eq!(edge.source(), IsogenyGraphNodeId(1));
        assert_eq!(edge.target(), IsogenyGraphNodeId(2));
        assert_eq!(edge.degree(), 2);
        assert_eq!(edge.kernel().points().len(), 2);
        assert!(matches!(edge.target_witness(), EdgeTargetWitness::Identity));
    }

    #[test]
    fn edge_can_store_an_explicit_target_witness_on_f41() {
        let curve = f41_curve();
        let kernel = sample_kernel(&curve);
        let expected_codomain = curve
            .scaled_by(F41::from_i64(3))
            .expect("non-zero scale should define a codomain");
        let witness = ShortWeierstrassIsomorphism::new(curve.clone(), F41::from_i64(3))
            .expect("non-zero scale should define an isomorphism");
        let edge = IsogenyGraphEdge::new(
            IsogenyGraphEdgeId(8),
            IsogenyGraphNodeId(10),
            IsogenyGraphNodeId(11),
            kernel,
            EdgeTargetWitness::Explicit(witness),
        );

        match edge.target_witness() {
            EdgeTargetWitness::Identity => panic!("expected an explicit witness"),
            EdgeTargetWitness::Explicit(witness) => {
                assert_eq!(witness.domain(), &curve);
                assert_eq!(witness.codomain(), &expected_codomain);
            }
        }
    }
}
