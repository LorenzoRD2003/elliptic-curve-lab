use num_complex::Complex64;

use crate::numerics::ComplexApproxComparison;

/// A rank-two complex lattice `Λ = ℤω₁ + ℤω₂`.
///
/// This type stores an explicit ordered basis in `ℂ`. The current
/// implementation requires that basis to be non-degenerate and
/// positively oriented at construction time, so the associated ratio
/// `τ = ω₂ / ω₁` lives in the upper half-plane. Degenerate bases and bases
/// with non-positive orientation are rejected as distinct construction
/// failures.
#[derive(Clone, Debug, PartialEq)]
pub struct ComplexLattice {
    pub(super) omega1: Complex64,
    pub(super) omega2: Complex64,
}

/// One enumerated lattice point together with its integer coordinates.
///
/// The `value` field stores the concrete complex number `mω₁ + nω₂`, while
/// `m` and `n` record which integer combination produced it.
#[derive(Clone, Debug, PartialEq)]
pub struct LatticeIndexPoint {
    /// Integer coefficient of `ω₁`.
    pub m: i64,
    /// Integer coefficient of `ω₂`.
    pub n: i64,
    /// The corresponding lattice value `mω₁ + nω₂`.
    pub value: Complex64,
}

/// Coordinates in the standard half-open unit square `[0, 1) × [0, 1)`.
///
/// Relative to a lattice basis `ω₁, ω₂`, a pair `(u, v)` represents the
/// complex point `uω₁ + vω₂`. When `0 ≤ u < 1` and `0 ≤ v < 1`, that point
/// lies in the chosen half-open fundamental parallelogram for the quotient
/// `ℂ / Λ`.
///
/// This type validates its inputs at construction time. A successful value is
/// guaranteed to satisfy `0 ≤ u < 1` and `0 ≤ v < 1`, so it already lives in
/// the chosen canonical representative region.
#[derive(Clone, Debug, PartialEq)]
pub struct FundamentalParallelogramCoordinate {
    /// Real coordinate along `ω₁`.
    pub(super) u: f64,
    /// Real coordinate along `ω₂`.
    pub(super) v: f64,
}

/// A canonical point of the complex torus `ℂ / Λ`.
///
/// The quotient identifies complex numbers that differ by a lattice vector:
/// `z ∼ z + λ` for `λ ∈ Λ`. This type stores one canonical representative of
/// that equivalence class as reduced coordinates in the chosen half-open
/// fundamental parallelogram.
///
/// The lattice itself is not stored inside the point. That keeps the value
/// small and educational, but it also means the point only has meaning
/// relative to the `ComplexLattice` that created it.
#[derive(Clone, Debug)]
pub struct ComplexTorusPoint {
    pub(super) coordinate: FundamentalParallelogramCoordinate,
}

/// Approximate comparison between two complex representatives modulo one
/// ambient lattice.
///
/// Given two complex numbers `z_left` and `z_right`, the comparison searches
/// over lattice shifts
///
/// `z_right + mω₁ + nω₂`
///
/// inside the finite box `-r ≤ m, n ≤ r` and keeps the shift producing the
/// smallest residual norm against `z_left`.
///
/// The stored approximate verdict therefore answers:
///
/// “does `z_left` agree approximately with `z_right` modulo this lattice,
/// within the searched shift budget?”
#[derive(Clone, Debug, PartialEq)]
pub struct ComplexModuloLatticeComparison {
    pub(super) original_right: Complex64,
    pub(super) best_shift: LatticeIndexPoint,
    pub(super) comparison: ComplexApproxComparison,
    pub(super) search_radius: usize,
}
