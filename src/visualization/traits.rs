/// Small shared visualization surface for educational text output.
///
/// Domain-specific visualization traits can extend this base trait with richer
/// operations, but every visualizable object should at least provide:
///
/// - a compact representation
/// - a richer description
pub trait Visualizable {
    /// Returns a compact human-readable representation.
    fn format_compact(&self) -> String;

    /// Returns a richer educational description.
    fn describe(&self) -> String;
}
