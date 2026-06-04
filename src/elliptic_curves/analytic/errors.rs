/// Typed error surface for the educational complex-analytic elliptic-curve
/// milestone.
#[derive(Clone, Debug, PartialEq)]
pub enum AnalyticCurveError {
    TauNotInUpperHalfPlane,
    DegenerateLattice,
    NonPositiveLatticeOrientation,
    InvalidTruncationRadius,
    InvalidSeriesPrecision,
    NearlySingularAnalyticCurve,
    PointTooCloseToLatticePoint,
    PointNotInFundamentalParallelogram,
    InvalidModularMatrix,
    NonPositiveImaginaryPartAfterModularAction,
    UnsupportedNormalization,
    NumericalComparisonFailed,
}
