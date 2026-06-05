use num_complex::Complex64;

use super::WeierstrassCubicRoots;
use crate::fields::ComplexApprox;
use crate::numerics::ApproxTolerance;

/// Coarse geometric configuration of a recovered Weierstrass cubic.
///
/// This classification is intentionally separate from near-collision
/// diagnostics. A triple may be geometrically “three approximately real” or
/// “one approximately real plus an approximately conjugate pair” while still
/// being numerically close to a repeated-root regime.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CubicRootConfiguration {
    /// All three roots have approximately vanishing imaginary part.
    ThreeApproximatelyReal,
    /// Exactly one root is approximately real and the other two are
    /// approximately complex conjugates.
    OneApproximatelyRealTwoApproximatelyConjugate,
    /// No simpler real/conjugate pattern was detected under the tolerance.
    GenericComplex,
}

/// Pairwise-separation status for a recovered cubic-root triple.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CubicRootSeparation {
    /// The roots remain numerically separated under the requested tolerance.
    WellSeparated,
    /// Some pair becomes approximately repeated under the requested tolerance.
    NearlyRepeated,
}

/// Structured classification report for one cubic-root triple.
///
/// The report stores both a geometric classification and a separate
/// near-collision status so callers can distinguish “shape” from “distance to
/// singularity”.
#[derive(Clone, Debug, PartialEq)]
pub struct CubicRootConfigurationReport {
    roots: WeierstrassCubicRoots,
    configuration: CubicRootConfiguration,
    separation: CubicRootSeparation,
    tolerance: ApproxTolerance,
    min_pairwise_distance: f64,
    conjugate_pair_residual: Option<f64>,
}

impl CubicRootConfigurationReport {
    /// Builds a classification report for `roots` under `tolerance`.
    ///
    /// Complexity: `Θ(1)`.
    pub fn new(roots: WeierstrassCubicRoots, tolerance: ApproxTolerance) -> Self {
        let min_pairwise_distance = roots.min_pairwise_distance();
        let separation = if roots.has_repeated_root_approx(tolerance) {
            CubicRootSeparation::NearlyRepeated
        } else {
            CubicRootSeparation::WellSeparated
        };

        let all_real = all_roots_are_approximately_real(&roots, tolerance);
        let conjugate_candidate = best_real_split_conjugate_candidate(&roots, tolerance);
        let conjugate_pair_residual = if all_real {
            None
        } else {
            conjugate_candidate.map(|candidate| candidate.residual)
        };
        let configuration = if all_real {
            CubicRootConfiguration::ThreeApproximatelyReal
        } else if conjugate_candidate.is_some_and(|candidate| candidate.is_close) {
            CubicRootConfiguration::OneApproximatelyRealTwoApproximatelyConjugate
        } else {
            CubicRootConfiguration::GenericComplex
        };

        Self {
            roots,
            configuration,
            separation,
            tolerance,
            min_pairwise_distance,
            conjugate_pair_residual,
        }
    }

    /// Returns the classified roots.
    pub fn roots(&self) -> &WeierstrassCubicRoots {
        &self.roots
    }

    /// Returns the coarse geometric configuration.
    pub fn configuration(&self) -> CubicRootConfiguration {
        self.configuration
    }

    /// Returns the pairwise-separation status.
    pub fn separation(&self) -> CubicRootSeparation {
        self.separation
    }

    /// Returns the tolerance used for the classification.
    pub fn tolerance(&self) -> ApproxTolerance {
        self.tolerance
    }

    /// Returns the smallest pairwise root distance.
    pub fn min_pairwise_distance(&self) -> f64 {
        self.min_pairwise_distance
    }

    /// Returns the best conjugate-pair residual among splits with an
    /// approximately real leftover root, when such a split exists.
    pub fn conjugate_pair_residual(&self) -> Option<f64> {
        self.conjugate_pair_residual
    }
}

/// Returns the coarse geometric configuration of `roots`.
///
/// Complexity: `Θ(1)`.
pub fn classify_cubic_root_configuration(
    roots: &WeierstrassCubicRoots,
    tolerance: ApproxTolerance,
) -> CubicRootConfiguration {
    cubic_root_configuration_report(roots, tolerance).configuration()
}

/// Builds a structured cubic-root configuration report.
///
/// Complexity: `Θ(1)`.
pub fn cubic_root_configuration_report(
    roots: &WeierstrassCubicRoots,
    tolerance: ApproxTolerance,
) -> CubicRootConfigurationReport {
    CubicRootConfigurationReport::new(roots.clone(), tolerance)
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct ConjugatePairCandidate {
    residual: f64,
    is_close: bool,
}

fn all_roots_are_approximately_real(
    roots: &WeierstrassCubicRoots,
    tolerance: ApproxTolerance,
) -> bool {
    roots
        .roots()
        .into_iter()
        .all(|root| is_approximately_real(root, tolerance))
}

fn best_real_split_conjugate_candidate(
    roots: &WeierstrassCubicRoots,
    tolerance: ApproxTolerance,
) -> Option<ConjugatePairCandidate> {
    let roots = roots.roots();
    let splits = [(0usize, 1usize, 2usize), (1, 0, 2), (2, 0, 1)];
    let mut best: Option<ConjugatePairCandidate> = None;

    for (real_index, first_pair, second_pair) in splits {
        if !is_approximately_real(roots[real_index], tolerance) {
            continue;
        }

        let pair_residual = conjugate_pair_residual(roots[first_pair], roots[second_pair]);
        let candidate = ConjugatePairCandidate {
            residual: pair_residual,
            is_close: ComplexApprox::eq_with_tolerance(
                roots[first_pair],
                &roots[second_pair].conj(),
                tolerance,
            ),
        };

        match best {
            None => best = Some(candidate),
            Some(previous) if candidate.residual < previous.residual => best = Some(candidate),
            _ => {}
        }
    }

    best
}

fn is_approximately_real(root: &Complex64, tolerance: ApproxTolerance) -> bool {
    ComplexApprox::eq_with_tolerance(root, &root.conj(), tolerance)
}

fn conjugate_pair_residual(first: &Complex64, second: &Complex64) -> f64 {
    (*first - second.conj()).norm()
}
