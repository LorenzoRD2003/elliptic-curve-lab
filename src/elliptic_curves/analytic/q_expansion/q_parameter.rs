use num_complex::Complex64;

use crate::elliptic_curves::analytic::UpperHalfPlanePoint;

/// The classical modular parameter `q = e^{2π i τ}`
/// attached to a point `τ` in the upper half-plane.
///
/// This is the small complex parameter used in Fourier and `q`-expansions of
/// modular forms and elliptic functions. Since `Im(τ) > 0`, its absolute value
/// satisfies `|q| = e^{-2π Im(τ)} < 1`,
/// so `q` always lies strictly inside the open unit disc.
#[derive(Clone, Debug, PartialEq)]
pub struct ModularQParameter {
    tau: UpperHalfPlanePoint,
    q: Complex64,
}

impl ModularQParameter {
    /// Builds the modular parameter `q = e^{2π i τ}` from a validated
    /// upper-half-plane point `τ`.
    ///
    /// This stores both `τ` and the derived `q` so later `q`-expansion
    /// routines can report their analytic input and the actual small
    /// parameter they used.
    pub fn from_tau(tau: UpperHalfPlanePoint) -> ModularQParameter {
        let q = (Complex64::new(0.0, std::f64::consts::TAU) * *tau.tau()).exp();

        ModularQParameter { tau, q }
    }

    /// Returns the upper-half-plane parameter `τ` used to build this `q`.
    pub fn tau(&self) -> &UpperHalfPlanePoint {
        &self.tau
    }

    /// Returns the complex modular parameter `q = e^{2π i τ}`.
    pub fn q(&self) -> &Complex64 {
        &self.q
    }

    /// Returns the complex absolute value `|q|`.
    ///
    /// For `τ ∈ ℍ`, this equals `e^{-2π Im(τ)}` and is therefore always
    /// strictly smaller than `1`.
    pub fn absolute_value(&self) -> f64 {
        self.q.norm()
    }
}
