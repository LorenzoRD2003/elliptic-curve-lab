use num_complex::Complex64;

use super::{
    AnalyticCurveError, ComplexLattice, LatticeSumTruncation, UpperHalfPlanePoint, g4_sum, g6_sum,
};
use crate::fields::ComplexApprox;

/// Approximate classical analytic invariants attached to a complex lattice `Λ`.
/// Here every quantity is only an approximation, because `G₄` and `G₆` are
/// themselves computed by finite square-box truncation.
///
/// - `g₂(Λ) = 60 G₄(Λ)`
/// - `g₃(Λ) = 140 G₆(Λ)`
/// - `Δ(Λ) = g₂(Λ)^3 - 27 g₃(Λ)^2`
/// - `j(Λ) = 1728 g₂(Λ)^3 / Δ(Λ)`
///
/// The values `g₂`, `g₃`, and `Δ` depend on the scaling of the lattice basis,
/// while `j` is the classical homothety-invariant quantity.
#[derive(Clone, Debug, PartialEq)]
pub struct AnalyticInvariants {
    /// Approximation to `g₂(Λ) = 60 G₄(Λ)`.
    pub g2: Complex64,
    /// Approximation to `g₃(Λ) = 140 G₆(Λ)`.
    pub g3: Complex64,
    /// Approximation to `Δ(Λ) = g₂(Λ)^3 - 27 g₃(Λ)^2`.
    pub discriminant: Complex64,
    /// Approximation to `j(Λ) = 1728 g₂(Λ)^3 / Δ(Λ)`.
    pub j_invariant: Complex64,
    /// Truncation policy used to compute `g₂` and `g₃`.
    pub truncation: LatticeSumTruncation,
}

/// Approximates the classical analytic invariant `g₂(Λ)`.
///
/// This implementation uses the relation `g₂(Λ) = 60 G₄(Λ)`,
/// where `G₄(Λ)` is approximated by a finite punctured square-box lattice sum.
pub fn analytic_g2(
    lattice: &ComplexLattice,
    truncation: LatticeSumTruncation,
) -> Result<Complex64, AnalyticCurveError> {
    let g4 = g4_sum(lattice, truncation)?;
    Ok(Complex64::new(60.0, 0.0) * g4.value)
}

/// Approximates the classical analytic invariant `g₃(Λ)`.
///
/// This implementation uses the relation `g₃(Λ) = 140 G₆(Λ)`,
/// where `G₆(Λ)` is approximated by a finite punctured square-box lattice sum.
pub fn analytic_g3(
    lattice: &ComplexLattice,
    truncation: LatticeSumTruncation,
) -> Result<Complex64, AnalyticCurveError> {
    let g6 = g6_sum(lattice, truncation)?;
    Ok(Complex64::new(140.0, 0.0) * g6.value)
}

/// Computes the classical discriminant expression `Δ = g₂^3 - 27 g₃^2`.
pub fn analytic_discriminant(g2: &Complex64, g3: &Complex64) -> Complex64 {
    g2.powu(3) - Complex64::new(27.0, 0.0) * g3.powu(2)
}

/// Computes the classical analytic `j`-invariant from `g₂` and `g₃`.
///
/// The formula is `j = 1728 g₂^3 / Δ`, where `Δ = g₂^3 - 27 g₃^2`.
///
/// If the discriminant is numerically too close to zero under the default
/// `ComplexApprox` tolerance policy, the function returns
/// [`AnalyticCurveError::NearlySingularAnalyticCurve`] instead of dividing by a
/// value that is too unstable to interpret honestly.
pub fn analytic_j_invariant(
    g2: &Complex64,
    g3: &Complex64,
) -> Result<Complex64, AnalyticCurveError> {
    let discriminant = analytic_discriminant(g2, g3);

    if ComplexApprox::is_zero_with_tolerance(&discriminant, ComplexApprox::default_tolerance()) {
        return Err(AnalyticCurveError::NearlySingularAnalyticCurve);
    }

    Ok(Complex64::new(1728.0, 0.0) * g2.powu(3) / discriminant)
}

/// Computes the approximate analytic invariants attached to a complex lattice.
///
/// This bundles together the truncated approximations to `g₂`, `g₃`, `Δ`, and
/// `j` produced from one common truncation policy.
pub fn analytic_invariants(
    lattice: &ComplexLattice,
    truncation: LatticeSumTruncation,
) -> Result<AnalyticInvariants, AnalyticCurveError> {
    let g2 = analytic_g2(lattice, truncation)?;
    let g3 = analytic_g3(lattice, truncation)?;
    let discriminant = analytic_discriminant(&g2, &g3);
    let j_invariant = analytic_j_invariant(&g2, &g3)?;

    Ok(AnalyticInvariants {
        g2,
        g3,
        discriminant,
        j_invariant,
        truncation,
    })
}

/// Computes the approximate analytic invariants attached to the standard
/// lattice `Λ_τ = ℤ + ℤτ`.
///
/// This is a convenience wrapper around [`ComplexLattice::from_tau`] followed
/// by [`analytic_invariants`].
pub fn analytic_invariants_from_tau(
    tau: &UpperHalfPlanePoint,
    truncation: LatticeSumTruncation,
) -> Result<AnalyticInvariants, AnalyticCurveError> {
    let lattice = ComplexLattice::from_tau(tau.clone());
    analytic_invariants(&lattice, truncation)
}

#[cfg(test)]
mod tests {
    use num_complex::Complex64;

    use super::{
        AnalyticInvariants, analytic_discriminant, analytic_g2, analytic_g3, analytic_invariants,
        analytic_invariants_from_tau, analytic_j_invariant,
    };
    use crate::elliptic_curves::analytic::{
        AnalyticCurveError, ComplexLattice, LatticeSumTruncation, UpperHalfPlanePoint, g4_sum,
        g6_sum,
    };
    use crate::fields::{ComplexApprox, Field};

    fn standard_square_lattice() -> ComplexLattice {
        ComplexLattice::from_tau(UpperHalfPlanePoint::tau_i())
    }

    fn standard_hexagonal_lattice() -> ComplexLattice {
        ComplexLattice::from_tau(UpperHalfPlanePoint::tau_rho())
    }

    #[test]
    fn analytic_g2_matches_sixty_times_g4() {
        let lattice = standard_square_lattice();
        let truncation = LatticeSumTruncation::default_educational();

        let g2 = analytic_g2(&lattice, truncation).unwrap();
        let g4 = g4_sum(&lattice, truncation).unwrap();

        assert!(ComplexApprox::eq(
            &g2,
            &(Complex64::new(60.0, 0.0) * g4.value)
        ));
    }

    #[test]
    fn analytic_g3_matches_one_hundred_forty_times_g6() {
        let lattice = standard_hexagonal_lattice();
        let truncation = LatticeSumTruncation::default_educational();

        let g3 = analytic_g3(&lattice, truncation).unwrap();
        let g6 = g6_sum(&lattice, truncation).unwrap();

        assert!(ComplexApprox::eq(
            &g3,
            &(Complex64::new(140.0, 0.0) * g6.value)
        ));
    }

    #[test]
    fn analytic_discriminant_uses_the_classical_formula() {
        let g2 = Complex64::new(3.0, -1.0);
        let g3 = Complex64::new(-2.0, 4.0);

        let discriminant = analytic_discriminant(&g2, &g3);
        let expected = g2.powu(3) - Complex64::new(27.0, 0.0) * g3.powu(2);

        assert_eq!(discriminant, expected);
    }

    #[test]
    fn analytic_j_invariant_rejects_nearly_singular_input() {
        let g2 = Complex64::new(0.0, 0.0);
        let g3 = Complex64::new(0.0, 0.0);

        assert_eq!(
            analytic_j_invariant(&g2, &g3),
            Err(AnalyticCurveError::NearlySingularAnalyticCurve)
        );
    }

    #[test]
    fn analytic_j_invariant_matches_the_classical_formula_when_stable() {
        let g2 = Complex64::new(12.0, 1.0);
        let g3 = Complex64::new(4.0, -2.0);

        let j = analytic_j_invariant(&g2, &g3).unwrap();
        let discriminant = analytic_discriminant(&g2, &g3);
        let expected = Complex64::new(1728.0, 0.0) * g2.powu(3) / discriminant;

        assert!(ComplexApprox::eq(&j, &expected));
    }

    #[test]
    fn analytic_invariants_bundle_all_derived_quantities() {
        let lattice = standard_square_lattice();
        let truncation = LatticeSumTruncation::default_educational();

        let invariants = analytic_invariants(&lattice, truncation).unwrap();
        let g2 = analytic_g2(&lattice, truncation).unwrap();
        let g3 = analytic_g3(&lattice, truncation).unwrap();
        let discriminant = analytic_discriminant(&g2, &g3);
        let j_invariant = analytic_j_invariant(&g2, &g3).unwrap();

        assert_eq!(
            invariants,
            AnalyticInvariants {
                g2,
                g3,
                discriminant,
                j_invariant,
                truncation,
            }
        );
    }

    #[test]
    fn analytic_invariants_from_tau_matches_standard_lattice_construction() {
        let tau = UpperHalfPlanePoint::tau_rho();
        let truncation = LatticeSumTruncation::larger_for_comparison();

        let from_tau = analytic_invariants_from_tau(&tau, truncation).unwrap();
        let lattice = ComplexLattice::from_tau(tau);
        let from_lattice = analytic_invariants(&lattice, truncation).unwrap();

        assert_eq!(from_tau, from_lattice);
    }
}
