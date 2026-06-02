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
    pub(super) nodes: Vec<IsogenyGraphNode<C>>,
    pub(super) edges: Vec<IsogenyGraphEdge<C>>,
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
}

/// Small BFS-style builder for educational `ℓ`-isogeny graphs.
#[derive(Clone, Debug)]
pub struct IsogenyGraphBuilder<C: GraphCurveModel> {
    pub(super) start_curve: C,
    pub(super) ell: usize,
    pub(super) max_depth: usize,
    pub(super) deduplicate_by_base_field_isomorphism: bool,
}

impl<C: GraphCurveModel> IsogenyGraphBuilder<C> {
    /// Starts a graph build from one chosen representative and one prime degree `ℓ`.
    pub fn new(start_curve: C, degree: usize) -> Self {
        Self {
            start_curve,
            ell: degree,
            max_depth: 1,
            deduplicate_by_base_field_isomorphism: true,
        }
    }

    /// Sets the maximum BFS depth measured in edge traversals from the start node.
    pub fn max_depth(mut self, max_depth: usize) -> Self {
        self.max_depth = max_depth;
        self
    }

    /// Enables or disables base-field-isomorphism deduplication for newly
    /// discovered codomain curves.
    ///
    /// Exact representative equality is still deduplicated even when this flag
    /// is `false`.
    pub fn deduplicate_by_base_field_isomorphism(mut self, yes: bool) -> Self {
        self.deduplicate_by_base_field_isomorphism = yes;
        self
    }
}
