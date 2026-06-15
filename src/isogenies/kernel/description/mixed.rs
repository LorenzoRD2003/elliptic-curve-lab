use crate::elliptic_curves::traits::CurveModel;

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
