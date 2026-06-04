use crate::fields::ComplexApprox;

use super::{super::AnalyticCurveError, FundamentalParallelogramCoordinate};

impl FundamentalParallelogramCoordinate {
    /// Builds a validated coordinate pair `(u, v)`.
    ///
    /// This constructor accepts exactly the half-open unit-square region
    /// `0 ≤ u < 1`, `0 ≤ v < 1`. Non-finite inputs return
    /// [`AnalyticCurveError::NumericalComparisonFailed`]. Finite but
    /// out-of-range inputs return
    /// [`AnalyticCurveError::PointNotInFundamentalParallelogram`].
    pub fn new(u: f64, v: f64) -> Result<Self, AnalyticCurveError> {
        if !u.is_finite() || !v.is_finite() {
            return Err(AnalyticCurveError::NumericalComparisonFailed);
        }

        if Self::raw_is_in_half_open_unit_square(u, v) {
            Ok(Self { u, v })
        } else {
            Err(AnalyticCurveError::PointNotInFundamentalParallelogram)
        }
    }

    /// Returns the coordinate along `ω₁`.
    pub fn u(&self) -> f64 {
        self.u
    }

    /// Returns the coordinate along `ω₂`.
    pub fn v(&self) -> f64 {
        self.v
    }

    /// Returns whether this pair lies in the half-open unit square.
    ///
    /// Successful constructor calls preserve this invariant, so valid
    /// `FundamentalParallelogramCoordinate` values always satisfy it.
    pub fn is_in_half_open_unit_square(&self) -> bool {
        Self::raw_is_in_half_open_unit_square(self.u, self.v)
    }

    /// Reduces raw coordinates modulo the unit square to a canonical representative.
    ///
    /// The reduction uses `rem_euclid(1)` on each coordinate, then snaps
    /// values numerically very close to either `0` or `1` back to `0` using
    /// [`ComplexApprox::default_tolerance`].
    pub fn reduce_mod_unit_square(u: f64, v: f64) -> Result<Self, AnalyticCurveError> {
        Self::new(reduce_unit_coordinate(u), reduce_unit_coordinate(v))
    }

    fn raw_is_in_half_open_unit_square(u: f64, v: f64) -> bool {
        (0.0..1.0).contains(&u) && (0.0..1.0).contains(&v)
    }
}

fn reduce_unit_coordinate(value: f64) -> f64 {
    let reduced = value.rem_euclid(1.0);
    let tolerance = ComplexApprox::default_tolerance();

    if !reduced.is_finite() {
        reduced
    } else if tolerance.real_close(reduced, 0.0) || tolerance.real_close(reduced, 1.0) {
        0.0
    } else {
        reduced
    }
}
