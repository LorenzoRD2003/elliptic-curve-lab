use crate::visualization::Visualizable;

/// Shared educational visualization surface for polynomial representations.
pub trait VisualizablePolynomial: Visualizable {
    /// Returns a compact human-readable polynomial string.
    fn format_polynomial(&self) -> String {
        self.format_compact()
    }

    /// Returns a richer educational description of the polynomial.
    fn describe_polynomial(&self) -> String {
        self.describe()
    }
}
