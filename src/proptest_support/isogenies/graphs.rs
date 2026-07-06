use std::collections::HashSet;

use num_bigint::BigUint;
use proptest::prelude::*;

use crate::elliptic_curves::{
    ShortWeierstrassCurve,
    traits::{AffineCurveModel, CurveModel},
};
use crate::fields::traits::*;
use crate::isogenies::{
    graphs::{
        IsogenyGraph, IsogenyGraphEdge, IsogenyGraphEdgeId, IsogenyGraphNode, IsogenyGraphNodeId,
        edge::EdgeTargetWitness,
    },
    kernel::IsogenyKernel,
};

type F41 = crate::fields::Fp41;
type Curve41 = ShortWeierstrassCurve<F41>;

/// Small structural ordinary `2`-volcano fixture for graph-side floor search.
///
/// The generated graph has one surface vertex, no horizontal surface edges, and
/// a complete regular downward tree of the requested depth. Every non-floor
/// vertex has outgoing degree `ℓ + 1 = 3`, and every floor vertex has outgoing
/// degree `1`, namely the edge back to its parent.
#[derive(Clone, Debug)]
pub struct VolcanicFloorSearchCase {
    graph: IsogenyGraph<Curve41>,
    prime: BigUint,
    start: IsogenyGraphNodeId,
    floor_nodes: Vec<IsogenyGraphNodeId>,
    depth: usize,
}

impl VolcanicFloorSearchCase {
    /// Returns the structural volcano graph.
    pub fn graph(&self) -> &IsogenyGraph<Curve41> {
        &self.graph
    }

    /// Returns the volcano prime `ℓ`.
    pub fn prime(&self) -> &BigUint {
        &self.prime
    }

    /// Returns the start node for floor search.
    pub fn start(&self) -> IsogenyGraphNodeId {
        self.start
    }

    /// Returns the floor nodes.
    pub fn floor_nodes(&self) -> &[IsogenyGraphNodeId] {
        &self.floor_nodes
    }

    /// Returns the structural volcano depth.
    pub fn depth(&self) -> usize {
        self.depth
    }
}

/// Generates small complete structural `2`-volcanoes for floor-search tests.
pub fn arb_volcanic_floor_search_case() -> BoxedStrategy<VolcanicFloorSearchCase> {
    (1usize..=3)
        .prop_map(build_complete_binary_two_volcano)
        .boxed()
}

fn build_complete_binary_two_volcano(depth: usize) -> VolcanicFloorSearchCase {
    let curve = sample_curve();
    let kernel = sample_two_torsion_kernel(&curve);
    let mut nodes = vec![IsogenyGraphNode::new(IsogenyGraphNodeId(0), curve.clone())];
    let mut levels = vec![vec![IsogenyGraphNodeId(0)]];

    for level in 1..=depth {
        let previous = levels[level - 1].clone();
        let mut current_level = Vec::new();
        let child_count_per_parent = if level == 1 { 3 } else { 2 };

        for _parent in previous {
            for _ in 0..child_count_per_parent {
                let node_id = IsogenyGraphNodeId(nodes.len());
                nodes.push(IsogenyGraphNode::new(node_id, curve.clone()));
                current_level.push(node_id);
            }
        }
        levels.push(current_level);
    }

    let mut edges = Vec::new();
    for level in 0..depth {
        let mut next_index = 0usize;
        let child_count_per_parent = if level == 0 { 3 } else { 2 };

        for parent in &levels[level] {
            for _ in 0..child_count_per_parent {
                let child = levels[level + 1][next_index];
                next_index += 1;
                push_structural_edge(&mut edges, *parent, child, &kernel);
                push_structural_edge(&mut edges, child, *parent, &kernel);
            }
        }
    }

    VolcanicFloorSearchCase {
        graph: IsogenyGraph {
            fully_expanded_nodes: vec![true; nodes.len()],
            nodes,
            edges,
        },
        prime: BigUint::from(2u8),
        start: IsogenyGraphNodeId(0),
        floor_nodes: levels[depth].clone(),
        depth,
    }
}

fn push_structural_edge(
    edges: &mut Vec<IsogenyGraphEdge<Curve41>>,
    source: IsogenyGraphNodeId,
    target: IsogenyGraphNodeId,
    kernel: &IsogenyKernel<Curve41>,
) {
    let edge_id = IsogenyGraphEdgeId(edges.len());
    edges.push(IsogenyGraphEdge::new(
        edge_id,
        source,
        target,
        kernel.clone(),
        EdgeTargetWitness::Identity,
    ));
}

fn sample_curve() -> Curve41 {
    Curve41::new(F41::from_i64(2), F41::from_i64(3)).expect("shared F41 curve should stay valid")
}

fn sample_two_torsion_kernel(curve: &Curve41) -> IsogenyKernel<Curve41> {
    let point = curve
        .point(F41::from_i64(40), F41::zero())
        .expect("sample point should lie on the curve");

    IsogenyKernel::new(curve, HashSet::from([curve.identity(), point]))
        .expect("two-torsion subgroup should be valid")
}

pub(crate) fn touch_volcanic_floor_search_case_fields() {
    let _ = |case: VolcanicFloorSearchCase| {
        let _ = case.graph;
        let _ = case.prime;
        let _ = case.start;
        let _ = case.floor_nodes;
        let _ = case.depth;
    };
}
