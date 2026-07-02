use crate::isogenies::graphs::GraphCurveModel;
use core::hash::Hash;

/// Stable identifier for one stored graph node.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct IsogenyGraphNodeId(pub usize);

/// One stored representative curve in the educational isogeny graph.
///
/// Nodes deliberately store a single representative of the underlying
/// base-field isomorphism class. Later summaries can recover the
/// `j`-invariant on demand from that representative, without caching a second
/// source of truth in the node itself.
#[derive(Clone, Debug)]
pub struct IsogenyGraphNode<C: GraphCurveModel>
where
    C::Point: Clone,
{
    id: IsogenyGraphNodeId,
    representative: C,
}

impl<C: GraphCurveModel> IsogenyGraphNode<C>
where
    C::Point: Clone,
{
    /// Builds a graph node from one chosen representative curve.
    pub(crate) fn new(id: IsogenyGraphNodeId, representative: C) -> Self {
        Self { id, representative }
    }

    /// Returns the node identifier.
    pub fn id(&self) -> IsogenyGraphNodeId {
        self.id
    }

    /// Returns the stored representative curve.
    pub(crate) fn representative(&self) -> &C {
        &self.representative
    }

    /// Computes the `j`-invariant from the stored representative on demand.
    pub fn j_invariant(&self) -> C::Elem {
        self.representative.j_invariant()
    }
}

#[cfg(test)]
mod tests {

    use crate::elliptic_curves::ShortWeierstrassCurve;
    use crate::fields::traits::Field;
    use crate::isogenies::graphs::{IsogenyGraphNode, IsogenyGraphNodeId};

    type F41 = crate::fields::Fp41;

    fn f41_curve() -> ShortWeierstrassCurve<F41> {
        ShortWeierstrassCurve::<F41>::new(F41::from_i64(2), F41::from_i64(3)).expect("valid curve")
    }

    #[test]
    fn node_exposes_id_representative_and_j_invariant() {
        let curve = f41_curve();
        let node = IsogenyGraphNode::new(IsogenyGraphNodeId(3), curve.clone());

        assert_eq!(node.id(), IsogenyGraphNodeId(3));
        assert_eq!(node.representative(), &curve);
        assert!(F41::eq(&node.j_invariant(), &curve.j_invariant()));
    }
}
