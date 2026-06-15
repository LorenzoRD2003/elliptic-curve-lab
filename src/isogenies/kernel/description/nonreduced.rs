/// Purely infinitesimal kernel data.
#[derive(Clone, Debug)]
pub struct NonReducedKernelDescription {
    /// Total degree of the infinitesimal kernel contribution.
    degree: usize,
    /// Human-readable mathematical label for the current construction.
    label: String,
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
