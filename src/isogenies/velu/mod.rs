use std::hash::Hash;

use crate::elliptic_curves::GroupCurveModel;
use crate::isogenies::{IsogenyError, IsogenyKernel};

mod core;
mod short_weierstrass;

#[cfg(test)]
mod tests;

pub use core::VeluIsogeny;
pub use short_weierstrass::{
    DualVeluIsogeny, verify_left_dual_relation, verify_right_dual_relation,
};

/// Private capability trait that fills in the model-specific parts of Vélu's construction.
///
/// The generic core in [`VeluIsogeny`] handles the parts that do not depend on
/// a concrete curve model:
///
/// - validating or building the finite kernel subgroup
/// - storing the domain, codomain, kernel, and degree
/// - applying the universal evaluation logic "off-curve is an error, kernel
///   points map to the codomain identity"
///
/// Each concrete curve model that wants to support Vélu must provide the
/// model-specific pieces through this trait:
///
/// - how to compute the codomain curve from the domain and kernel
/// - how to evaluate a non-kernel point once the universal cases have been
///   excluded by the generic core
pub(super) trait SupportsVeluConstruction: GroupCurveModel + Clone
where
    Self::Point: Clone + Eq + Hash,
{
    fn velu_codomain_from_kernel(
        domain: &Self,
        kernel: &IsogenyKernel<Self>,
    ) -> Result<Self, IsogenyError>
    where
        Self: Sized;

    fn velu_evaluate_non_kernel_point(
        isogeny: &VeluIsogeny<Self>,
        point: &Self::Point,
    ) -> Result<Self::Point, IsogenyError>
    where
        Self: Sized;
}
