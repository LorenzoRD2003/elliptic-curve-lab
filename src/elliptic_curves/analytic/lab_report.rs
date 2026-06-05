use core::fmt;

use num_complex::Complex64;

use super::{
    AnalyticCurveError, AnalyticInvariants, AnalyticWeierstrassCurve, ComplexLattice,
    EllipticFunctionTruncation, LatticeSumTruncation, ModularQParameter, TorusToCurveMapResult,
    UpperHalfPlanePoint, analytic_invariants, map_torus_point_to_curve,
};
use crate::elliptic_curves::ShortWeierstrassCurve;
use crate::fields::ComplexApprox;

/// Shared access to the ambient upper-half-plane parameter and its standard
/// lattice.
pub trait HasAnalyticLatticeContext {
    /// Returns the upper-half-plane parameter `τ`.
    fn tau(&self) -> &UpperHalfPlanePoint;

    /// Returns the associated lattice `Λ_τ = ℤ + ℤτ`.
    fn lattice(&self) -> &ComplexLattice;
}

/// Distinguishes the classical special upper-half-plane parameters used
/// throughout the analytic examples.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SpecialTauKind {
    /// The square-lattice point `τ = i`.
    TauI,
    /// The equianharmonic point `τ = ρ = -1/2 + (√3/2)i`.
    TauRho,
    /// Any other currently non-special sample point.
    Generic,
}

/// Distinguishes the classical special `j`-invariant regimes visible in the
/// first analytic experiments.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SpecialJKind {
    /// `j` is numerically close to `1728`, as for the square lattice.
    Near1728,
    /// `j` is numerically close to `0`, as for the equianharmonic lattice.
    NearZero,
    /// No special cancellation is currently detected.
    Generic,
}

/// Small value object for the short-Weierstrass companion
/// `y² = x³ + ax + b` attached to an analytic cubic.
///
/// This keeps the short-form coefficients and their `j`-invariant visible
/// without exposing the full generic algebraic curve type inside a
/// higher-level educational report.
#[derive(Clone, Debug, PartialEq)]
pub struct AnalyticShortWeierstrassModel {
    a: Complex64,
    b: Complex64,
    j_invariant: Complex64,
}

impl AnalyticShortWeierstrassModel {
    /// Builds the short-Weierstrass companion attached to an analytic cubic.
    pub fn from_analytic_curve(curve: &AnalyticWeierstrassCurve) -> Self {
        let short_curve: ShortWeierstrassCurve<ComplexApprox> = curve.as_short_weierstrass();

        Self {
            a: *short_curve.a(),
            b: *short_curve.b(),
            j_invariant: short_curve.j_invariant(),
        }
    }

    /// Returns the short-form coefficient `a`.
    pub fn a(&self) -> &Complex64 {
        &self.a
    }

    /// Returns the short-form coefficient `b`.
    pub fn b(&self) -> &Complex64 {
        &self.b
    }

    /// Returns the short-model `j`-invariant.
    pub fn j_invariant(&self) -> &Complex64 {
        &self.j_invariant
    }
}

impl fmt::Display for AnalyticShortWeierstrassModel {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "y^2 = x^3 + ({:?})x + ({:?})", self.a, self.b)
    }
}

/// Aggregated educational report for one complex-analytic lattice experiment.
///
/// This bundles together the main objects computed from a single parameter
/// `τ`:
///
/// - the standard lattice `Λ_τ = ℤ + ℤτ`
/// - the modular cusp parameter `q = e^{2π i τ}`
/// - the approximate analytic invariants `g₂`, `g₃`, `Δ`, and `j`
/// - the analytic cubic `y² = 4x³ - g₂x - g₃`
/// - its short-Weierstrass companion `y² = x³ + ax + b`
/// - lightweight classifications of the special `τ` and `j` cases
///
/// The point of this report is not to add new mathematics, but to keep the
/// first analytic milestone coherent as one inspectable bundle.
#[derive(Clone, Debug, PartialEq)]
pub struct ComplexAnalyticCurveLabReport {
    tau: UpperHalfPlanePoint,
    lattice: ComplexLattice,
    q_parameter: ModularQParameter,
    invariants: AnalyticInvariants,
    analytic_curve: AnalyticWeierstrassCurve,
    short_model: AnalyticShortWeierstrassModel,
    special_tau_kind: SpecialTauKind,
    special_j_kind: SpecialJKind,
}

impl ComplexAnalyticCurveLabReport {
    /// Builds the aggregated analytic report attached to one `τ`.
    ///
    /// The only numerical knob here is the lattice truncation used to compute
    /// the approximate invariants; the `q`-parameter itself is exact at the
    /// current `Complex64` level.
    pub fn from_tau(
        tau: UpperHalfPlanePoint,
        truncation: LatticeSumTruncation,
    ) -> Result<Self, AnalyticCurveError> {
        let lattice = ComplexLattice::from_tau(tau.clone());
        let q_parameter = ModularQParameter::from_tau(tau.clone());
        let invariants = analytic_invariants(&lattice, truncation)?;
        let analytic_curve = AnalyticWeierstrassCurve::new(invariants.g2, invariants.g3)?;
        let short_model = AnalyticShortWeierstrassModel::from_analytic_curve(&analytic_curve);
        let special_tau_kind = classify_special_tau(&tau);
        let special_j_kind = classify_special_j(&invariants.j_invariant);

        Ok(Self {
            tau,
            lattice,
            q_parameter,
            invariants,
            analytic_curve,
            short_model,
            special_tau_kind,
            special_j_kind,
        })
    }

    /// Returns the original upper-half-plane parameter `τ`.
    pub fn tau(&self) -> &UpperHalfPlanePoint {
        &self.tau
    }

    /// Returns the standard lattice `Λ_τ = ℤ + ℤτ`.
    pub fn lattice(&self) -> &ComplexLattice {
        &self.lattice
    }

    /// Returns the derived modular parameter `q = e^{2π i τ}`.
    pub fn q_parameter(&self) -> &ModularQParameter {
        &self.q_parameter
    }

    /// Returns the approximate analytic invariants.
    pub fn invariants(&self) -> &AnalyticInvariants {
        &self.invariants
    }

    /// Returns the analytic cubic `y² = 4x³ - g₂x - g₃`.
    pub fn analytic_curve(&self) -> &AnalyticWeierstrassCurve {
        &self.analytic_curve
    }

    /// Returns the short-Weierstrass companion model.
    pub fn short_model(&self) -> &AnalyticShortWeierstrassModel {
        &self.short_model
    }

    /// Returns the lightweight classification of the supplied `τ`.
    pub fn special_tau_kind(&self) -> SpecialTauKind {
        self.special_tau_kind
    }

    /// Returns the lightweight classification of the computed `j`.
    pub fn special_j_kind(&self) -> SpecialJKind {
        self.special_j_kind
    }
}

impl HasAnalyticLatticeContext for ComplexAnalyticCurveLabReport {
    fn tau(&self) -> &UpperHalfPlanePoint {
        &self.tau
    }

    fn lattice(&self) -> &ComplexLattice {
        &self.lattice
    }
}

/// Aggregated experiment for the analytic uniformization map
/// `ℂ / Λ → E(ℂ)`, `z ↦ (℘(z), ℘′(z))`.
///
/// This report keeps one lattice/curve pair together with several sampled
/// torus representatives and the corresponding mapped curve points.
#[derive(Clone, Debug, PartialEq)]
pub struct UniformizationExperimentReport {
    tau: UpperHalfPlanePoint,
    lattice: ComplexLattice,
    curve: AnalyticWeierstrassCurve,
    sampled_points: Vec<TorusToCurveMapResult>,
    all_points_lie_on_curve: bool,
}

impl UniformizationExperimentReport {
    /// Builds one uniformization experiment from explicit complex
    /// representatives in `ℂ`.
    ///
    /// The sampled points may include genuine lattice points; those are mapped
    /// to [`super::AnalyticCurvePoint::Infinity`] through the same pole convention
    /// already used by [`map_torus_point_to_curve`].
    pub fn from_sample_points(
        tau: UpperHalfPlanePoint,
        sample_points: Vec<Complex64>,
        invariant_truncation: LatticeSumTruncation,
        function_truncation: EllipticFunctionTruncation,
        tolerance: crate::numerics::ApproxTolerance,
    ) -> Result<Self, AnalyticCurveError> {
        let lattice = ComplexLattice::from_tau(tau.clone());
        let curve = AnalyticWeierstrassCurve::from_lattice(&lattice, invariant_truncation)?;
        let sampled_points = sample_points
            .into_iter()
            .map(|z| {
                map_torus_point_to_curve(
                    &lattice,
                    z,
                    invariant_truncation,
                    function_truncation,
                    tolerance,
                )
            })
            .collect::<Result<Vec<_>, _>>()?;
        let all_points_lie_on_curve = sampled_points.iter().all(|point| point.lies_on_curve());

        Ok(Self {
            tau,
            lattice,
            curve,
            sampled_points,
            all_points_lie_on_curve,
        })
    }

    /// Returns the original upper-half-plane parameter.
    pub fn tau(&self) -> &UpperHalfPlanePoint {
        &self.tau
    }

    /// Returns the underlying lattice.
    pub fn lattice(&self) -> &ComplexLattice {
        &self.lattice
    }

    /// Returns the common analytic cubic used in the experiment.
    pub fn curve(&self) -> &AnalyticWeierstrassCurve {
        &self.curve
    }

    /// Returns the sampled torus-to-curve map results.
    pub fn sampled_points(&self) -> &[TorusToCurveMapResult] {
        &self.sampled_points
    }

    /// Returns whether every sampled point was accepted as lying on the curve.
    pub fn all_points_lie_on_curve(&self) -> bool {
        self.all_points_lie_on_curve
    }
}

impl HasAnalyticLatticeContext for UniformizationExperimentReport {
    fn tau(&self) -> &UpperHalfPlanePoint {
        &self.tau
    }

    fn lattice(&self) -> &ComplexLattice {
        &self.lattice
    }
}

fn classify_special_tau(tau: &UpperHalfPlanePoint) -> SpecialTauKind {
    if ComplexApprox::eq_with_tolerance(
        tau.tau(),
        UpperHalfPlanePoint::tau_i().tau(),
        ComplexApprox::default_tolerance(),
    ) {
        return SpecialTauKind::TauI;
    }

    if ComplexApprox::eq_with_tolerance(
        tau.tau(),
        UpperHalfPlanePoint::tau_rho().tau(),
        ComplexApprox::default_tolerance(),
    ) {
        return SpecialTauKind::TauRho;
    }

    SpecialTauKind::Generic
}

fn classify_special_j(j_invariant: &Complex64) -> SpecialJKind {
    if ComplexApprox::eq_with_tolerance(
        j_invariant,
        &Complex64::new(1728.0, 0.0),
        crate::numerics::ApproxTolerance::new(1.0e-3, 1.0e-6),
    ) {
        return SpecialJKind::Near1728;
    }

    if ComplexApprox::eq_with_tolerance(
        j_invariant,
        &Complex64::new(0.0, 0.0),
        crate::numerics::ApproxTolerance::new(1.0e-4, 1.0e-6),
    ) {
        return SpecialJKind::NearZero;
    }

    SpecialJKind::Generic
}

#[cfg(test)]
mod tests {
    use num_complex::Complex64;

    use super::{
        ComplexAnalyticCurveLabReport, SpecialJKind, SpecialTauKind,
        UniformizationExperimentReport, classify_special_j,
    };
    use crate::elliptic_curves::analytic::{
        EllipticFunctionTruncation, LatticeSumTruncation, TorusToCurveValues, UpperHalfPlanePoint,
    };
    use crate::fields::ComplexApprox;

    #[test]
    fn lab_report_keeps_tau_lattice_and_q_parameter_consistent() {
        let tau = UpperHalfPlanePoint::tau_generic_example();
        let report = ComplexAnalyticCurveLabReport::from_tau(
            tau.clone(),
            LatticeSumTruncation::new(12).unwrap(),
        )
        .unwrap();

        assert_eq!(report.tau(), &tau);
        assert_eq!(
            report.lattice(),
            &crate::elliptic_curves::analytic::ComplexLattice::from_tau(tau.clone())
        );
        assert_eq!(report.q_parameter().tau(), &tau);
        assert_eq!(report.lattice().tau().unwrap(), tau);
    }

    #[test]
    fn lab_report_keeps_j_consistent_across_all_model_surfaces() {
        let report = ComplexAnalyticCurveLabReport::from_tau(
            UpperHalfPlanePoint::tau_generic_example(),
            LatticeSumTruncation::new(12).unwrap(),
        )
        .unwrap();

        assert!(ComplexApprox::eq_with_tolerance(
            &report.invariants().j_invariant,
            &report.analytic_curve().j_invariant().unwrap(),
            crate::numerics::ApproxTolerance::new(1.0e-6, 1.0e-6),
        ));
        assert!(ComplexApprox::eq_with_tolerance(
            &report.invariants().j_invariant,
            report.short_model().j_invariant(),
            crate::numerics::ApproxTolerance::new(1.0e-6, 1.0e-6),
        ));
    }

    #[test]
    fn special_tau_kind_distinguishes_i_rho_and_generic_cases() {
        let tau_i_report = ComplexAnalyticCurveLabReport::from_tau(
            UpperHalfPlanePoint::tau_i(),
            LatticeSumTruncation::new(12).unwrap(),
        )
        .unwrap();
        let tau_rho_report = ComplexAnalyticCurveLabReport::from_tau(
            UpperHalfPlanePoint::tau_rho(),
            LatticeSumTruncation::new(12).unwrap(),
        )
        .unwrap();
        let generic_report = ComplexAnalyticCurveLabReport::from_tau(
            UpperHalfPlanePoint::tau_generic_example(),
            LatticeSumTruncation::new(12).unwrap(),
        )
        .unwrap();

        assert_eq!(tau_i_report.special_tau_kind(), SpecialTauKind::TauI);
        assert_eq!(tau_rho_report.special_tau_kind(), SpecialTauKind::TauRho);
        assert_eq!(generic_report.special_tau_kind(), SpecialTauKind::Generic);
    }

    #[test]
    fn special_j_kind_reflects_the_classical_special_values() {
        let tau_i_report = ComplexAnalyticCurveLabReport::from_tau(
            UpperHalfPlanePoint::tau_i(),
            LatticeSumTruncation::new(12).unwrap(),
        )
        .unwrap();
        let tau_rho_report = ComplexAnalyticCurveLabReport::from_tau(
            UpperHalfPlanePoint::tau_rho(),
            LatticeSumTruncation::new(12).unwrap(),
        )
        .unwrap();
        let generic_report = ComplexAnalyticCurveLabReport::from_tau(
            UpperHalfPlanePoint::tau_generic_example(),
            LatticeSumTruncation::new(12).unwrap(),
        )
        .unwrap();

        assert_eq!(tau_i_report.special_j_kind(), SpecialJKind::Near1728);
        assert_eq!(tau_rho_report.special_j_kind(), SpecialJKind::NearZero);
        assert_eq!(generic_report.special_j_kind(), SpecialJKind::Generic);
    }

    #[test]
    fn j_classifier_is_generic_away_from_the_special_values() {
        assert_eq!(
            classify_special_j(&Complex64::new(120.0, -1692.0)),
            SpecialJKind::Generic
        );
    }

    #[test]
    fn uniformization_report_derives_global_curve_membership_from_samples() {
        let report = UniformizationExperimentReport::from_sample_points(
            UpperHalfPlanePoint::tau_i(),
            vec![
                Complex64::new(0.0, 0.0),
                Complex64::new(0.3, 0.2),
                Complex64::new(0.5, 0.0),
            ],
            LatticeSumTruncation::new(16).unwrap(),
            EllipticFunctionTruncation::new(14).unwrap(),
            crate::numerics::ApproxTolerance::new(1.0e-2, 1.0e-2),
        )
        .unwrap();

        assert_eq!(
            report.all_points_lie_on_curve(),
            report
                .sampled_points()
                .iter()
                .all(|point| point.lies_on_curve())
        );
    }

    #[test]
    fn uniformization_report_keeps_the_same_curve_across_all_samples() {
        let report = UniformizationExperimentReport::from_sample_points(
            UpperHalfPlanePoint::tau_i(),
            vec![Complex64::new(0.3, 0.2), Complex64::new(0.5, 0.0)],
            LatticeSumTruncation::new(16).unwrap(),
            EllipticFunctionTruncation::new(14).unwrap(),
            crate::numerics::ApproxTolerance::new(1.0e-2, 1.0e-2),
        )
        .unwrap();

        assert!(
            report
                .sampled_points()
                .iter()
                .all(|point| point.curve() == report.curve())
        );
    }

    #[test]
    fn uniformization_report_can_include_both_finite_points_and_a_pole() {
        let report = UniformizationExperimentReport::from_sample_points(
            UpperHalfPlanePoint::tau_i(),
            vec![Complex64::new(0.0, 0.0), Complex64::new(0.3, 0.2)],
            LatticeSumTruncation::new(16).unwrap(),
            EllipticFunctionTruncation::new(14).unwrap(),
            crate::numerics::ApproxTolerance::new(1.0e-2, 1.0e-2),
        )
        .unwrap();

        assert!(matches!(
            report.sampled_points()[0].values(),
            TorusToCurveValues::Pole
        ));
        assert!(matches!(
            report.sampled_points()[1].values(),
            TorusToCurveValues::FiniteValues { .. }
        ));
    }
}
