use num_complex::Complex64;

use crate::elliptic_curves::analytic::{AnalyticCurveError, ComplexLattice, UpperHalfPlanePoint};

/// One recovered ordered period basis for an analytic elliptic curve.
///
/// This type intentionally wraps one validated [`ComplexLattice`] instead of
/// storing `ω₁`, `ω₂`, and `τ` as parallel fields. That keeps the non-
/// degeneracy and positive-orientation invariants in one place.
#[derive(Clone, Debug, PartialEq)]
pub struct RecoveredPeriodBasis {
    lattice: ComplexLattice,
}

impl RecoveredPeriodBasis {
    /// Builds a recovered period basis from explicit periods.
    pub fn new(omega1: Complex64, omega2: Complex64) -> Result<Self, AnalyticCurveError> {
        Ok(Self {
            lattice: ComplexLattice::new(omega1, omega2)?,
        })
    }

    /// Wraps an already validated complex lattice as a recovered period basis.
    pub fn from_lattice(lattice: ComplexLattice) -> Self {
        Self { lattice }
    }

    /// Returns the first recovered period `ω₁`.
    pub fn omega1(&self) -> &Complex64 {
        self.lattice.omega1()
    }

    /// Returns the second recovered period `ω₂`.
    pub fn omega2(&self) -> &Complex64 {
        self.lattice.omega2()
    }

    /// Returns the recovered period ratio `τ = ω₂ / ω₁`.
    ///
    /// Because this type only stores validated lattices, the ratio is expected
    /// to lie in the upper half-plane.
    pub fn tau(&self) -> UpperHalfPlanePoint {
        self.lattice
            .tau()
            .expect("validated recovered period basis must have tau in the upper half-plane")
    }

    /// Returns the underlying validated lattice.
    pub fn lattice(&self) -> &ComplexLattice {
        &self.lattice
    }

    /// Consumes the wrapper and returns the underlying validated lattice.
    pub fn into_lattice(self) -> ComplexLattice {
        self.lattice
    }

    /// Returns the oriented area of the recovered period parallelogram.
    pub fn oriented_area(&self) -> f64 {
        self.lattice.oriented_area()
    }

    /// Returns the covolume of the recovered period lattice.
    pub fn covolume(&self) -> f64 {
        self.lattice.covolume()
    }
}
