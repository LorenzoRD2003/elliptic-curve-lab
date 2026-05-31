use crate::visualization::Visualizable;

/// Shared educational visualization surface for field elements and closely
/// related field-domain values.
pub trait VisualizableField: Visualizable + Sized {
    /// Returns a compact field-element-style representation.
    fn format_elem(&self) -> String {
        self.format_compact()
    }

    /// Returns a compact visualization of the multiplicative inverse when it
    /// exists.
    fn inverse(&self) -> Option<String> {
        None
    }

    /// Explains addition when that concept is meaningful for the type.
    fn explain_add(_lhs: &Self, _rhs: &Self) -> Option<String> {
        None
    }

    /// Explains multiplication when that concept is meaningful for the type.
    fn explain_mul(_lhs: &Self, _rhs: &Self) -> Option<String> {
        None
    }

    /// Explains division when that concept is meaningful for the type.
    fn explain_div(_lhs: &Self, _rhs: &Self) -> Option<String> {
        None
    }
}
