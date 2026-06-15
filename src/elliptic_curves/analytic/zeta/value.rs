use num_complex::Complex64;

use crate::elliptic_curves::analytic::{
    EllipticFunctionTruncation, elliptic_functions::traits::HasPoleDistance,
};

/// One truncated approximation to the Weierstrass `ζ`-function.
#[derive(Clone, Debug, PartialEq)]
pub struct WeierstrassZetaApprox {
    z: Complex64,
    value: Complex64,
    truncation: EllipticFunctionTruncation,
    terms_used: usize,
    pole_distance: f64,
}

impl WeierstrassZetaApprox {
    pub(crate) fn new(
        z: Complex64,
        value: Complex64,
        truncation: EllipticFunctionTruncation,
        terms_used: usize,
        pole_distance: f64,
    ) -> Self {
        Self {
            z,
            value,
            truncation,
            terms_used,
            pole_distance,
        }
    }

    /// Returns the original evaluation point supplied by the caller.
    pub fn z(&self) -> &Complex64 {
        &self.z
    }

    /// Returns the approximate complex value produced by the truncation.
    pub fn value(&self) -> &Complex64 {
        &self.value
    }

    /// Returns the truncation policy used for this approximation.
    pub fn truncation(&self) -> EllipticFunctionTruncation {
        self.truncation
    }

    /// Returns the number of nonzero lattice terms that were summed.
    pub fn terms_used(&self) -> usize {
        self.terms_used
    }

    /// Returns the smallest inspected distance to a lattice pole.
    pub fn pole_distance(&self) -> f64 {
        self.pole_distance
    }
}

impl HasPoleDistance for WeierstrassZetaApprox {
    fn pole_distance(&self) -> f64 {
        self.pole_distance
    }
}
