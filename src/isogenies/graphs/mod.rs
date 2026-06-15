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
//! internal graph-model witnesses through one hidden capability trait.

pub(crate) mod builder;
pub(crate) mod connected_components;
pub(crate) mod cycles;
pub(crate) mod edge;
pub(crate) mod endomorphisms;
pub(crate) mod error;
pub(crate) mod node;
pub(crate) mod torsion;
pub(crate) mod verification;
pub(crate) mod volcano;
mod weak_graph;

pub use builder::{IsogenyGraph, IsogenyGraphBuilder};
pub use edge::{IsogenyGraphEdge, IsogenyGraphEdgeId};
pub use error::IsogenyGraphError;
pub use node::{IsogenyGraphNode, IsogenyGraphNodeId};
pub(crate) use torsion::GraphTorsionCurveModel;
pub use verification::{IsogenyGraphVerificationReport, ReverseEdgeStatus};
pub use volcano::{VolcanoLikeLayering, VolcanoRole};

#[allow(unused_imports)]
pub(crate) use endomorphisms::{
    EndomorphismVolcanoReport, IsogenyEdgeEndomorphismRelation, IsogenyEdgeEndomorphismReport,
    IsogenyGraphEndomorphismEdgeReport, IsogenyGraphEndomorphismNodeReport,
    IsogenyGraphEndomorphismReport, VolcanoHeuristicComparison,
};

use crate::elliptic_curves::{
    ShortWeierstrassCurve,
    short_weierstrass::isomorphisms::ShortWeierstrassIsomorphism,
    traits::{CurveIsomorphism, GroupCurveModel, HasJInvariant},
};
use crate::fields::traits::Field;

/// Hidden infrastructure trait packaging the capabilities the current graph
/// representation needs.
///
/// Public graph values still depend on this bound today, so the trait remains
/// public for type-checking purposes, but it is intentionally not part of the
/// documented user-facing graph story.
///
/// It packages:
///
/// - group-curve structure so kernels remain honest subgroup objects
/// - a `j`-invariant for later graph summaries
/// - one explicit same-family isomorphism witness type so edges can record how
///   a raw codomain is transported onto the chosen node representative
#[doc(hidden)]
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
