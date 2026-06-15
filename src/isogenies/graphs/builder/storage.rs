use std::hash::Hash;

use crate::isogenies::graphs::{
    GraphCurveModel, IsogenyGraphEdge, IsogenyGraphNode, IsogenyGraphNodeId,
};

/// Small stored `ℓ`-isogeny graph built from explicit node representatives and
/// explicit edge kernels.
///
/// - nodes are stored in a dense vector
/// - edges are stored in a flat vector
/// - incoming and outgoing adjacency is recovered by filtering the edge list
#[derive(Clone, Debug)]
pub struct IsogenyGraph<C: GraphCurveModel>
where
    C::Point: Clone + Eq + Hash,
    C::IsomorphismWitness: Clone + core::fmt::Debug,
{
    pub(crate) nodes: Vec<IsogenyGraphNode<C>>,
    pub(crate) edges: Vec<IsogenyGraphEdge<C>>,
}

impl<C: GraphCurveModel> IsogenyGraph<C>
where
    C::Point: Clone + Eq + Hash,
    C::IsomorphismWitness: Clone + core::fmt::Debug,
{
    /// Returns every stored node in dense id order.
    pub fn nodes(&self) -> &[IsogenyGraphNode<C>] {
        &self.nodes
    }

    /// Returns how many nodes are stored in the graph.
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }

    /// Returns every stored edge.
    pub fn edges(&self) -> &[IsogenyGraphEdge<C>] {
        &self.edges
    }

    /// Returns how many edges are stored in the graph.
    pub fn edge_count(&self) -> usize {
        self.edges.len()
    }

    /// Returns the node with the given id when that dense id is present.
    ///
    /// The current representation expects `node.id().0` to match the node's
    /// vector index. The extra equality check keeps the lookup honest even if a
    /// malformed internal graph is ever constructed during future refactors.
    pub fn node(&self, id: IsogenyGraphNodeId) -> Option<&IsogenyGraphNode<C>> {
        self.nodes.get(id.0).filter(|node| node.id() == id)
    }

    /// Iterates over every stored outgoing edge from `source`.
    pub fn outgoing_edges(
        &self,
        source: IsogenyGraphNodeId,
    ) -> impl Iterator<Item = &IsogenyGraphEdge<C>> + '_ {
        self.edges
            .iter()
            .filter(move |edge| edge.source() == source)
    }

    /// Iterates over every stored incoming edge to `target`.
    pub fn incoming_edges(
        &self,
        target: IsogenyGraphNodeId,
    ) -> impl Iterator<Item = &IsogenyGraphEdge<C>> + '_ {
        self.edges
            .iter()
            .filter(move |edge| edge.target() == target)
    }

    /// Returns the number of outgoing edges from `node`.
    pub fn out_degree(&self, node: IsogenyGraphNodeId) -> usize {
        self.outgoing_edges(node).count()
    }

    /// Returns the number of incoming edges to `node`.
    pub fn in_degree(&self, node: IsogenyGraphNodeId) -> usize {
        self.incoming_edges(node).count()
    }

    #[cfg(test)]
    pub(crate) fn push_edge(&mut self, edge: IsogenyGraphEdge<C>) {
        self.edges.push(edge);
    }
}
