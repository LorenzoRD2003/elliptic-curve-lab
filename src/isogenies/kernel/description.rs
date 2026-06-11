use crate::elliptic_curves::traits::CurveModel;
use crate::isogenies::IsogenyKernel;

/// Public description of an isogeny kernel.
///
/// This enum is intentionally broader than [`IsogenyKernel`]. The latter
/// models a concrete reduced finite subgroup given explicitly by points. That
/// is the right representation for separable small-field constructions such as
/// Vélu isogenies, but it is not mathematically honest for inseparable maps.
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

/// Reduced kernel data that can be exposed through explicit points.
#[derive(Clone, Debug)]
pub enum ReducedKernelDescription<C: CurveModel> {
    /// A validated explicit finite subgroup.
    RationalPointSubgroup(IsogenyKernel<C>),
    /// A finite reduced subgroup visible through explicit points, but not
    /// currently reified as an [`IsogenyKernel`].
    FiniteSubgroupSchemeVisibleAsPoints {
        /// Explicit points currently visible on the chosen base field.
        points: Vec<C::Point>,
        /// Degree contributed by the reduced subgroup.
        degree: usize,
    },
}

/// Purely infinitesimal kernel data.
#[derive(Clone, Debug)]
pub struct NonReducedKernelDescription {
    /// Total degree of the infinitesimal kernel contribution.
    degree: usize,
    /// Human-readable mathematical label for the current construction.
    label: String,
}

/// Mixed kernel data with both reduced and infinitesimal parts.
#[derive(Clone, Debug)]
pub struct MixedKernelDescription<C: CurveModel> {
    /// Explicit reduced points currently visible on the chosen base field.
    reduced_points: Vec<C::Point>,
    /// Degree contributed by the reduced part.
    reduced_degree: usize,
    /// Degree contributed by the infinitesimal part.
    infinitesimal_degree: usize,
    /// Optional human-readable label for the infinitesimal side.
    label: Option<String>,
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

impl<C: CurveModel> ReducedKernelDescription<C> {
    /// Returns the reduced degree carried by this description.
    pub fn degree(&self) -> usize {
        match self {
            Self::RationalPointSubgroup(kernel) => kernel.degree(),
            Self::FiniteSubgroupSchemeVisibleAsPoints { degree, .. } => *degree,
        }
    }

    /// Returns the explicit points carried by this reduced description.
    pub fn points(&self) -> &[C::Point] {
        match self {
            Self::RationalPointSubgroup(kernel) => kernel.points(),
            Self::FiniteSubgroupSchemeVisibleAsPoints { points, .. } => points.as_slice(),
        }
    }
}

impl NonReducedKernelDescription {
    /// Builds a purely infinitesimal kernel description.
    pub fn new(degree: usize, label: impl Into<String>) -> Self {
        Self {
            degree,
            label: label.into(),
        }
    }

    /// Returns the infinitesimal degree.
    pub fn degree(&self) -> usize {
        self.degree
    }

    /// Returns the human-readable label.
    pub fn label(&self) -> &str {
        &self.label
    }
}

impl<C: CurveModel> MixedKernelDescription<C> {
    /// Builds a mixed kernel description with reduced and infinitesimal parts.
    pub fn new(
        reduced_points: Vec<C::Point>,
        reduced_degree: usize,
        infinitesimal_degree: usize,
        label: Option<String>,
    ) -> Self {
        Self {
            reduced_points,
            reduced_degree,
            infinitesimal_degree,
            label,
        }
    }

    /// Returns the visible reduced points.
    pub fn reduced_points(&self) -> &[C::Point] {
        self.reduced_points.as_slice()
    }

    /// Returns the reduced degree.
    pub fn reduced_degree(&self) -> usize {
        self.reduced_degree
    }

    /// Returns the infinitesimal degree.
    pub fn infinitesimal_degree(&self) -> usize {
        self.infinitesimal_degree
    }

    /// Returns the optional label for the infinitesimal side.
    pub fn label(&self) -> Option<&str> {
        self.label.as_deref()
    }
}
