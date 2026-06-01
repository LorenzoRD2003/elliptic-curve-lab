use super::IsogenyError;
use crate::elliptic_curves::CurveModel;

/// Minimal shared interface for explicit elliptic-curve isogeny objects.
///
/// This trait is intentionally austere. An isogeny is exposed only through:
///
/// - its domain curve
/// - its codomain curve
/// - its degree
/// - point evaluation
/// - the explicit kernel points used in the small finite educational setting
pub trait Isogeny<Domain, Codomain>
where
    Domain: CurveModel,
    Codomain: CurveModel,
{
    /// Returns the domain curve.
    fn domain(&self) -> &Domain;

    /// Returns the codomain curve.
    fn codomain(&self) -> &Codomain;

    /// Returns the degree of the isogeny.
    fn degree(&self) -> usize;

    /// Evaluates the isogeny at a point of the domain and returns a point of
    /// the codomain.
    ///
    /// In concrete implementations this will send every kernel point to
    /// the identity of the codomain and identify points that differ by a
    /// kernel element.
    fn evaluate(&self, point: &Domain::Point) -> Result<Codomain::Point, IsogenyError>;

    /// Returns the explicit kernel points used by the isogeny.
    ///
    /// This deliberately exposes the small finite representation
    /// instead of hiding the kernel behind a more opaque quotient object.
    fn kernel_points(&self) -> &[Domain::Point];
}
