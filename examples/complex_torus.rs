use num_complex::Complex64;

use elliptic_algorithms_lab::elliptic_curves::analytic::{
    AnalyticWeierstrassCurve, ComplexLattice, LatticeSumTruncation, UpperHalfPlanePoint,
};
use elliptic_algorithms_lab::visualization::Visualizable;

#[derive(Clone, Copy)]
enum SpecialExpectation {
    Square,
    Hexagonal,
    Generic,
}

fn print_tau_block(
    label: &str,
    tau: &UpperHalfPlanePoint,
    truncation: LatticeSumTruncation,
    expectation: SpecialExpectation,
) -> Result<(), String> {
    let lattice = ComplexLattice::from_tau(tau.clone());
    let g4 = lattice
        .g4_sum(truncation)
        .map_err(|error| format!("{error:?}"))?;
    let g6 = lattice
        .g6_sum(truncation)
        .map_err(|error| format!("{error:?}"))?;
    let invariants = lattice
        .analytic_invariants(truncation)
        .map_err(|error| format!("{error:?}"))?;
    let analytic_curve = AnalyticWeierstrassCurve::from_lattice(&lattice, truncation)
        .map_err(|error| format!("{error:?}"))?;
    let short_curve = analytic_curve.as_short_weierstrass();

    println!("{label}");
    println!("{}", "-".repeat(label.len()));
    println!("τ = {}", tau.tau().format_compact());
    println!("{}", lattice.describe());
    println!();
    println!("{}", g4.describe());
    println!();
    println!("{}", g6.describe());
    println!();
    println!("{}", invariants.describe());
    println!("analytic cubic model: {}", analytic_curve.format_compact());
    println!("short-Weierstrass model: {}", short_curve.format_compact());
    match expectation {
        SpecialExpectation::Square => {
            println!(
                "expected checks for the square lattice: |g₃| ≈ {:.6e}, |j - 1728| ≈ {:.6e}",
                invariants.g3().norm(),
                (*invariants.j_invariant() - Complex64::new(1728.0, 0.0)).norm()
            );
        }
        SpecialExpectation::Hexagonal => {
            println!(
                "expected checks for the equianharmonic lattice: |g₂| ≈ {:.6e}, |j| ≈ {:.6e}",
                invariants.g2().norm(),
                invariants.j_invariant().norm()
            );
        }
        SpecialExpectation::Generic => {
            println!("qualitative check: no special cancellation is expected in g₂, g₃, or j.");
        }
    }
    println!();

    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let truncation = LatticeSumTruncation::new(12).unwrap();
    let tau_i = UpperHalfPlanePoint::tau_i();
    let tau_rho = UpperHalfPlanePoint::tau_rho();
    let tau_generic = UpperHalfPlanePoint::new(Complex64::new(0.3, 1.2))?;

    println!("Complex tori and analytic invariants");
    println!("===============================================");
    println!();
    println!(
        "Note: we use a slightly larger truncation (r = {}) so the special cases τ = i and τ = ρ appear numerically closer to their expected invariants.",
        truncation.radius()
    );
    println!();

    print_tau_block("τ = i", &tau_i, truncation, SpecialExpectation::Square)?;
    print_tau_block("τ = ρ", &tau_rho, truncation, SpecialExpectation::Hexagonal)?;
    print_tau_block(
        "τ = 0.3 + 1.2i",
        &tau_generic,
        truncation,
        SpecialExpectation::Generic,
    )?;

    Ok(())
}
