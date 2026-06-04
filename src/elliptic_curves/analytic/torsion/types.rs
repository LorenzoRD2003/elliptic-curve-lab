use num_complex::Complex64;

use super::super::{
    AnalyticCurveMembershipReport, AnalyticCurvePoint, ApproxTolerance,
    FundamentalParallelogramCoordinate,
};
/// A reduced index class for a torus `n`-torsion point.
///
/// Relative to a lattice basis `Λ = ℤω₁ + ℤω₂`, the torus `n`-torsion subgroup is
/// `E[n] ≅ (1/n)Λ / Λ`.
///
/// Concretely, every `n`-torsion class has a representative of the form
/// `z = (a/n)ω₁ + (b/n)ω₂ mod Λ` with `0 ≤ a, b < n`. This type stores that
/// reduced pair `(a, b)` together with its common denominator `n`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct TorusTorsionIndex {
    pub(super) a: usize,
    pub(super) b: usize,
    pub(super) n: usize,
}

impl TorusTorsionIndex {
    /// Builds a validated reduced torus-torsion index.
    ///
    /// The current torus model accepts exactly those triples with `n > 0`,
    /// `0 ≤ a < n`, and `0 ≤ b < n`.
    pub fn new(a: usize, b: usize, n: usize) -> Result<Self, super::super::AnalyticCurveError> {
        if n == 0 || a >= n || b >= n {
            return Err(super::super::AnalyticCurveError::InvalidTorusTorsionIndex);
        }

        Ok(Self { a, b, n })
    }

    /// Returns the reduced coefficient of `ω₁`.
    pub fn a(&self) -> usize {
        self.a
    }

    /// Returns the reduced coefficient of `ω₂`.
    pub fn b(&self) -> usize {
        self.b
    }

    /// Returns the common torsion denominator `n`.
    pub fn n(&self) -> usize {
        self.n
    }

    /// Returns whether this reduced class is the torus identity `(0, 0; n)`.
    pub fn is_identity_class(&self) -> bool {
        self.a == 0 && self.b == 0
    }

    /// Returns whether this reduced class has exact torus order `n`.
    ///
    /// In the current torus model, the point indexed by `(a, b; n)` is
    /// primitive exactly when `gcd(a, b, n) = 1`.
    pub fn is_primitive(&self) -> bool {
        gcd3(self.a, self.b, self.n) == 1
    }
}

/// One explicit torus `n`-torsion point.
///
/// The stored data is intentionally redundant:
/// - [`TorusTorsionIndex`] records the arithmetic class `(a, b; n)`
/// - [`FundamentalParallelogramCoordinate`] records the canonical reduced
///   torus coordinate `(a/n, b/n)`
/// - `z` records the concrete complex representative `(a/n)ω₁ + (b/n)ω₂`
#[derive(Clone, Debug, PartialEq)]
pub struct TorusTorsionPoint {
    pub(super) index: TorusTorsionIndex,
    pub(super) coordinate: FundamentalParallelogramCoordinate,
    pub(super) z: Complex64,
}

impl TorusTorsionPoint {
    /// Returns the reduced arithmetic index `(a, b; n)`.
    pub fn index(&self) -> &TorusTorsionIndex {
        &self.index
    }

    /// Returns the canonical reduced torus coordinate `(a/n, b/n)`.
    pub fn coordinate(&self) -> &FundamentalParallelogramCoordinate {
        &self.coordinate
    }

    /// Returns the explicit complex representative `(a/n)ω₁ + (b/n)ω₂`.
    pub fn z(&self) -> &Complex64 {
        &self.z
    }

    /// Returns whether this torus point is the identity class.
    pub fn is_identity_class(&self) -> bool {
        self.index.is_identity_class()
    }
}

/// One approximate image of a torus `n`-torsion point on the analytic
/// Weierstrass curve attached to the same lattice.
///
/// This keeps both sides of the correspondence visible at once:
/// - `torus_point` stores the reduced class in `(1/n)Λ / Λ`
/// - `curve_point` stores its image under `z ↦ (℘(z), ℘′(z))`
/// - `membership_report` records the approximate residual in
///   `y² = 4x³ - g₂x - g₃`
///
/// The distinguished torus identity class `(a, b) = (0, 0)` maps to
/// [`AnalyticCurvePoint::Infinity`], since `℘` and `℘′` have poles at lattice
/// points and the complex torus compactifies to the cubic by sending those
/// poles to infinity.
#[derive(Clone, Debug, PartialEq)]
pub struct AnalyticTorsionPointApprox {
    pub(super) torus_point: TorusTorsionPoint,
    pub(super) curve_point: AnalyticCurvePoint,
    pub(super) membership_report: AnalyticCurveMembershipReport,
}

impl AnalyticTorsionPointApprox {
    /// Returns the torus-side `n`-torsion point.
    pub fn torus_point(&self) -> &TorusTorsionPoint {
        &self.torus_point
    }

    /// Returns the curve-side image under `z ↦ (℘(z), ℘′(z))`.
    pub fn curve_point(&self) -> &AnalyticCurvePoint {
        &self.curve_point
    }

    /// Returns the approximate curve-membership report for the image point.
    pub fn membership_report(&self) -> &AnalyticCurveMembershipReport {
        &self.membership_report
    }

    /// Returns whether the mapped point was accepted as lying on the cubic.
    pub fn lies_on_curve(&self) -> bool {
        self.membership_report.is_on_curve
    }
}

/// Approximate outcome of comparing analytic torus torsion against a complex
/// division-polynomial `x`-criterion.
#[derive(Clone, Debug, PartialEq)]
pub enum AnalyticDivisionPolynomialComparisonStatus {
    /// The torus identity class maps to the point at infinity, so no finite
    /// `x = ℘(z)` value is available to evaluate.
    PoleAtIdentity,
    /// The chosen `x`-criterion value is approximately zero.
    VanishesApproximately,
    /// The chosen `x`-criterion value is not approximately zero.
    DoesNotVanishApproximately,
}

/// Which branch of the even division-polynomial vanishing condition is active.
///
/// For even indices, the full curve-side condition is
/// `ψ_n(P) = y(P) ε_n(x(P)) = 0`, so vanishing can happen through the `y`
/// branch, the `ε_n(x)` branch, or both.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum EvenDivisionPolynomialVanishingBranch {
    /// `y(P) ≈ 0` while `ε_n(x(P))` is not approximately zero.
    YApproxZero,
    /// `ε_n(x(P)) ≈ 0` while `y(P)` is not approximately zero.
    XCriterionApproxZero,
    /// Both `y(P) ≈ 0` and `ε_n(x(P)) ≈ 0`.
    BothBranches,
    /// Neither branch is approximately zero.
    NeitherBranch,
}

/// Structured comparison between one analytic torsion point and the odd
/// division polynomial `ψ_n(x)` at `x = ℘(z)`.
#[derive(Clone, Debug, PartialEq)]
pub struct AnalyticOddDivisionPolynomialReport {
    pub(super) torsion_point: AnalyticTorsionPointApprox,
    pub(super) x_value: Complex64,
    pub(super) psi_n_x: Complex64,
    pub(super) absolute_value: f64,
    pub(super) status: AnalyticDivisionPolynomialComparisonStatus,
    pub(super) tolerance: ApproxTolerance,
}

impl AnalyticOddDivisionPolynomialReport {
    /// Returns the mapped analytic torsion point being compared.
    pub fn torsion_point(&self) -> &AnalyticTorsionPointApprox {
        &self.torsion_point
    }

    /// Returns the finite `x = ℘(z)` value.
    pub fn x_value(&self) -> &Complex64 {
        &self.x_value
    }

    /// Returns the evaluated odd division polynomial `ψ_n(x)`.
    pub fn psi_n_x(&self) -> &Complex64 {
        &self.psi_n_x
    }

    /// Returns `|ψ_n(x)|`.
    pub fn absolute_value(&self) -> f64 {
        self.absolute_value
    }

    /// Returns the comparison status.
    pub fn status(&self) -> &AnalyticDivisionPolynomialComparisonStatus {
        &self.status
    }

    /// Returns the tolerance used for the approximate zero check.
    pub fn tolerance(&self) -> ApproxTolerance {
        self.tolerance
    }
}

/// Structured comparison between one analytic torsion point and the stripped
/// even division-polynomial factor `ε_n(x)` at `x = ℘(z)`.
///
/// For even indices, the full curve-side vanishing condition is
/// `ψ_n(P) = y(P) ε_n(x(P)) = 0`, so this report stores both the
/// `x`-criterion value `ε_n(x(P))` and the active vanishing branch.
#[derive(Clone, Debug, PartialEq)]
pub struct AnalyticEvenDivisionPolynomialReport {
    pub(super) torsion_point: AnalyticTorsionPointApprox,
    pub(super) x_value: Complex64,
    pub(super) epsilon_n_x: Complex64,
    pub(super) absolute_value: f64,
    pub(super) branch: EvenDivisionPolynomialVanishingBranch,
    pub(super) status: AnalyticDivisionPolynomialComparisonStatus,
    pub(super) tolerance: ApproxTolerance,
}

impl AnalyticEvenDivisionPolynomialReport {
    /// Returns the mapped analytic torsion point being compared.
    pub fn torsion_point(&self) -> &AnalyticTorsionPointApprox {
        &self.torsion_point
    }

    /// Returns the finite `x = ℘(z)` value.
    pub fn x_value(&self) -> &Complex64 {
        &self.x_value
    }

    /// Returns the evaluated stripped even factor `ε_n(x)`.
    pub fn epsilon_n_x(&self) -> &Complex64 {
        &self.epsilon_n_x
    }

    /// Returns `|ε_n(x)|`.
    pub fn absolute_value(&self) -> f64 {
        self.absolute_value
    }

    /// Returns which branch of `y(P) ε_n(x(P)) ≈ 0` is active.
    pub fn branch(&self) -> &EvenDivisionPolynomialVanishingBranch {
        &self.branch
    }

    /// Returns the comparison status for the full even-index vanishing test.
    pub fn status(&self) -> &AnalyticDivisionPolynomialComparisonStatus {
        &self.status
    }

    /// Returns the tolerance used for the approximate zero checks.
    pub fn tolerance(&self) -> ApproxTolerance {
        self.tolerance
    }
}

/// Typed analytic comparison case for one torus torsion class.
///
/// This enum avoids mixing three mathematically distinct situations behind one
/// struct with many optional fields:
///
/// - [`Self::Pole`] for the identity class, which maps to infinity
/// - [`Self::Odd`] for odd `n`, where the tested criterion is `ψ_n(x)`
/// - [`Self::Even`] for even `n`, where the tested criterion is `ε_n(x)` and
///   the report also records the active `y(P) ε_n(x(P)) ≈ 0` branch
#[derive(Clone, Debug, PartialEq)]
pub enum AnalyticDivisionPolynomialComparisonCase {
    /// The torus identity class maps to the point at infinity, so no finite
    /// `x = ℘(z)` value is available.
    Pole {
        torsion_point: AnalyticTorsionPointApprox,
        tolerance: ApproxTolerance,
    },
    /// Odd-index comparison through `ψ_n(x)`.
    Odd(AnalyticOddDivisionPolynomialReport),
    /// Even-index comparison through `ε_n(x)` plus branch information.
    Even(AnalyticEvenDivisionPolynomialReport),
}

impl AnalyticDivisionPolynomialComparisonCase {
    /// Returns the mapped analytic torsion point being compared.
    pub fn torsion_point(&self) -> &AnalyticTorsionPointApprox {
        match self {
            Self::Pole { torsion_point, .. } => torsion_point,
            Self::Odd(report) => report.torsion_point(),
            Self::Even(report) => report.torsion_point(),
        }
    }

    /// Returns the tolerance used for this comparison.
    pub fn tolerance(&self) -> ApproxTolerance {
        match self {
            Self::Pole { tolerance, .. } => *tolerance,
            Self::Odd(report) => report.tolerance(),
            Self::Even(report) => report.tolerance(),
        }
    }

    /// Returns the comparison status.
    pub fn status(&self) -> AnalyticDivisionPolynomialComparisonStatus {
        match self {
            Self::Pole { .. } => AnalyticDivisionPolynomialComparisonStatus::PoleAtIdentity,
            Self::Odd(report) => report.status().clone(),
            Self::Even(report) => report.status().clone(),
        }
    }
}

fn gcd3(a: usize, b: usize, c: usize) -> usize {
    gcd(gcd(a, b), c)
}

fn gcd(a: usize, b: usize) -> usize {
    let mut x = a;
    let mut y = b;

    while y != 0 {
        let remainder = x % y;
        x = y;
        y = remainder;
    }

    x
}
