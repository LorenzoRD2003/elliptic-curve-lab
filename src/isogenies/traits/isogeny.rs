use crate::elliptic_curves::traits::CurveModel;
use crate::isogenies::{
    error::IsogenyError,
    kernel::{KernelDescription, ReducedKernelDescription},
};

/// Minimal shared interface for explicit elliptic-curve isogeny objects.
///
/// This trait is intentionally austere. An isogeny is exposed only through:
///
/// - its domain curve
/// - its codomain curve
/// - its degree
/// - point evaluation
/// - an honest public description of the kernel
pub trait Isogeny<Domain: CurveModel, Codomain: CurveModel> {
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

    /// Returns the current public kernel description for this isogeny.
    ///
    /// This is broader than a slice of rational points: inseparable maps may
    /// have nonreduced geometric kernels that are not honestly visible through
    /// explicit point data alone.
    fn kernel_description(&self) -> KernelDescription<Domain>;

    /// Returns the explicit rational kernel points currently visible in the
    /// kernel description.
    ///
    /// This is a convenience helper for small reduced examples. It should not
    /// be read as “the full geometric kernel” in inseparable settings.
    fn kernel_points(&self) -> Vec<Domain::Point> {
        match self.kernel_description() {
            KernelDescription::Reduced(ReducedKernelDescription::RationalPointSubgroup(kernel)) => {
                kernel.points().to_vec()
            }
            KernelDescription::Reduced(
                ReducedKernelDescription::FiniteSubgroupSchemeVisibleAsPoints { points, .. },
            ) => points,
            KernelDescription::Mixed(description) => description.reduced_points().to_vec(),
            KernelDescription::NonReduced(_) | KernelDescription::Unknown => Vec::new(),
        }
    }
}
