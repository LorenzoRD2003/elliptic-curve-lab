use crate::elliptic_curves::traits::CurveModel;

use crate::isogenies::kernel::description::{
    MixedKernelDescription, NonReducedKernelDescription, ReducedKernelDescription,
};

/// Public description of an isogeny kernel.
///
/// This enum is intentionally broader than [`crate::isogenies::IsogenyKernel`].
/// The latter models a concrete reduced finite subgroup given explicitly by
/// points. That is the right representation for separable small-field
/// constructions such as Vélu isogenies, but it is not mathematically honest
/// for inseparable maps.
///
/// For example, absolute Frobenius in characteristic `p` has geometric kernel
/// of degree `p`, yet that kernel may contribute no nontrivial rational points
/// at all. A point-only API would therefore misleadingly suggest a trivial
/// kernel. `KernelDescription` makes that distinction explicit.
#[derive(Clone, Debug)]
pub enum KernelDescription<C: CurveModel> {
    /// A fully reduced kernel description.
    Reduced(ReducedKernelDescription<C>),
    /// A purely infinitesimal nonreduced kernel with no reduced point data.
    NonReduced(NonReducedKernelDescription),
    /// A kernel with both reduced and infinitesimal contributions.
    Mixed(MixedKernelDescription<C>),
    /// The current implementation does not have an honest kernel description.
    Unknown,
}

impl<C: CurveModel> KernelDescription<C> {
    /// Returns the total kernel degree when the current description knows it.
    pub fn degree(&self) -> Option<usize> {
        match self {
            Self::Reduced(description) => Some(description.degree()),
            Self::NonReduced(description) => Some(description.degree()),
            Self::Mixed(description) => description
                .reduced_degree()
                .checked_mul(description.infinitesimal_degree()),
            Self::Unknown => None,
        }
    }

    /// Returns the reduced degree when the current description knows it.
    pub fn reduced_degree(&self) -> Option<usize> {
        match self {
            Self::Reduced(description) => Some(description.degree()),
            Self::NonReduced(_) => Some(0),
            Self::Mixed(description) => Some(description.reduced_degree()),
            Self::Unknown => None,
        }
    }

    /// Returns the infinitesimal degree when the current description knows it.
    pub fn infinitesimal_degree(&self) -> Option<usize> {
        match self {
            Self::Reduced(_) => Some(1),
            Self::NonReduced(description) => Some(description.degree()),
            Self::Mixed(description) => Some(description.infinitesimal_degree()),
            Self::Unknown => None,
        }
    }

    /// Returns the explicit rational points currently visible in the reduced
    /// part of the kernel description, if any.
    pub fn rational_points(&self) -> Option<&[C::Point]> {
        match self {
            Self::Reduced(description) => Some(description.points()),
            Self::Mixed(description) => Some(description.reduced_points()),
            Self::NonReduced(_) | Self::Unknown => None,
        }
    }

    /// Returns whether the current description is fully reduced.
    pub fn is_fully_reduced(&self) -> bool {
        matches!(self, Self::Reduced(_))
    }
}
