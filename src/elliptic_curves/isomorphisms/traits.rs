use crate::elliptic_curves::isomorphisms::CurveIsomorphismError;
use crate::elliptic_curves::traits::CurveModel;

/// Explicit curve isomorphisms between two concrete curve models.
///
/// This trait is intentionally small. An educational curve isomorphism is
/// currently exposed through:
///
/// - its domain curve
/// - its codomain curve
/// - point evaluation
///
pub trait CurveIsomorphism {
    type Domain: CurveModel;
    type Codomain: CurveModel;

    /// Returns the domain curve.
    fn domain(&self) -> &Self::Domain;

    /// Returns the codomain curve.
    fn codomain(&self) -> &Self::Codomain;

    /// Evaluates the isomorphism at one point of the domain curve.
    fn evaluate(
        &self,
        point: &<Self::Domain as CurveModel>::Point,
    ) -> Result<<Self::Codomain as CurveModel>::Point, CurveIsomorphismError>;
}
