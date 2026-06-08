//! The current graph layer is intentionally small. It focuses on:
//!
//! - one representative curve per node
//! - one explicit rational kernel per edge
//! - an optional explicit witness that transports a raw codomain curve onto
//!   the chosen representative of the target node
//!
//! This keeps the graph honest about the difference between:
//!
//! - “the Vélu codomain produced by this kernel”
//! - “the representative curve we chose to store for this isomorphism class”
//!
//! while still leaving room for future curve families to plug in their own
//! isomorphism witnesses through [`GraphCurveModel`].

use crate::elliptic_curves::invariants::HasJInvariant;
use crate::elliptic_curves::isomorphisms::{CurveIsomorphism, ShortWeierstrassIsomorphism};
use crate::elliptic_curves::short_weierstrass::ShortWeierstrassCurve;
use crate::elliptic_curves::traits::GroupCurveModel;
use crate::fields::Field;

pub mod builder;
pub mod connected_components;
pub mod cycles;
pub mod edge;
pub mod endomorphisms;
pub mod error;
pub mod node;
pub mod torsion;
pub mod verification;
pub mod volcano;

pub use builder::{IsogenyGraph, IsogenyGraphBuilder};
pub use connected_components::weakly_connected_components;
pub use cycles::{find_small_directed_cycles, has_directed_cycle};
pub use edge::{EdgeTargetWitness, IsogenyGraphEdge, IsogenyGraphEdgeId};
pub use endomorphisms::{
    EndomorphismVolcanoReport, IsogenyEdgeEndomorphismRelation, IsogenyEdgeEndomorphismReport,
    IsogenyGraphEndomorphismEdgeReport, IsogenyGraphEndomorphismNodeReport,
    IsogenyGraphEndomorphismReport, VolcanoHeuristicComparison,
};
pub use error::IsogenyGraphError;
pub use node::{IsogenyGraphNode, IsogenyGraphNodeId};
pub use torsion::cyclic_kernels_of_order;
pub use verification::{IsogenyGraphVerificationReport, ReverseEdgeStatus};
pub use volcano::{VolcanoLayering, VolcanoLikeLayering, VolcanoRole, infer_volcano_like_layers};

/// This trait packages the small collection of capabilities the current graph
/// representation needs:
///
/// - group-curve structure so kernels remain honest subgroup objects
/// - a `j`-invariant for later graph summaries
/// - one explicit same-family isomorphism witness type so edges can record how
///   a raw codomain is transported onto the chosen node representative
pub trait GraphCurveModel: GroupCurveModel + HasJInvariant + Clone
where
    Self::Point: Clone,
{
    /// Explicit same-family isomorphism witness used to transport codomain
    /// curves onto stored node representatives.
    type IsomorphismWitness: CurveIsomorphism<Domain = Self, Codomain = Self>;
}

impl<F: Field> GraphCurveModel for ShortWeierstrassCurve<F> {
    type IsomorphismWitness = ShortWeierstrassIsomorphism<F>;
}
