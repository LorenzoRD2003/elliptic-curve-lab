/// Shared educational visualization surface for algebraic values.
///
/// The goal of this trait is not to replace every specialized helper in the
/// `visualization` module. Instead, it captures the smallest common interface
/// that is useful across multiple mathematical domains in this crate:
///
/// - a compact textual form for an element-like value
/// - a richer description
/// - an optional compact view of the inverse
/// - optional explanations of addition and multiplication
/// - an optional explanation of division
///
/// Types that do not have meaningful arithmetic explanations yet may rely on
/// the default `None` implementations for the optional methods.
pub trait Visualizable: Sized {
    /// Returns a compact human-readable representation of the value.
    fn format_elem(&self) -> String;

    /// Returns a compatibility alias for [`Self::format_elem`].
    fn format_compact(&self) -> String {
        self.format_elem()
    }

    /// Returns a richer educational description of the value.
    fn describe(&self) -> String;

    /// Returns a compact human-readable representation of the inverse when that
    /// concept is available and the value is invertible.
    fn inverse(&self) -> Option<String> {
        None
    }

    /// Explains an addition involving two values of the same kind when that
    /// concept is supported by the visualization layer.
    fn explain_add(_lhs: &Self, _rhs: &Self) -> Option<String> {
        None
    }

    /// Explains a multiplication involving two values of the same kind when
    /// that concept is supported by the visualization layer.
    fn explain_mul(_lhs: &Self, _rhs: &Self) -> Option<String> {
        None
    }

    /// Explains a division involving two values of the same kind when that
    /// concept is supported by the visualization layer.
    fn explain_div(_lhs: &Self, _rhs: &Self) -> Option<String> {
        None
    }
}
