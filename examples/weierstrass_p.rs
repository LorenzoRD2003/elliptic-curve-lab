use elliptic_algorithms_lab::elliptic_curves::analytic::{
    TorusToCurveValues, verify_weierstrass_differential_equation,
};
use elliptic_algorithms_lab::visualization::fields::format_complex;
use elliptic_algorithms_lab::visualization::{
    describe_complex_lattice, describe_weierstrass_differential_equation,
    describe_weierstrass_p_approx, describe_weierstrass_p_derivative_approx,
};
use elliptic_algorithms_lab::{
    ApproxTolerance, ComplexLattice, EllipticFunctionTruncation, LatticeSumTruncation,
    UpperHalfPlanePoint, weierstrass_p, weierstrass_p_derivative,
};
use num_complex::Complex64;

fn show_point(
    lattice: &ComplexLattice,
    z: Complex64,
    invariant_truncation: LatticeSumTruncation,
    function_truncation: EllipticFunctionTruncation,
    tolerance: ApproxTolerance,
) -> Result<(), Box<dyn std::error::Error>> {
    let p = weierstrass_p(lattice, z, function_truncation)?;
    let p_prime = weierstrass_p_derivative(lattice, z, function_truncation)?;
    let report = verify_weierstrass_differential_equation(
        lattice,
        z,
        invariant_truncation,
        function_truncation,
        tolerance,
    )?;

    println!("{}", describe_weierstrass_p_approx(&p));
    println!("{}", describe_weierstrass_p_derivative_approx(&p_prime));
    println!("{}", describe_weierstrass_differential_equation(&report));
    println!();

    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let tau = UpperHalfPlanePoint::tau_i();
    let lattice = ComplexLattice::from_tau(tau.clone());
    let invariant_truncation = LatticeSumTruncation::new(4)?;
    let function_truncation = EllipticFunctionTruncation::new(6)?;
    let tolerance = ApproxTolerance::new(1.0e-4, 1.0e-2);

    let points = [
        Complex64::new(0.20, 0.15),
        Complex64::new(0.35, 0.40),
        Complex64::new(0.45, 0.10),
        Complex64::new(0.15, 0.70),
    ];

    println!("Weierstrass ℘ and the differential equation");
    println!("=========================================================");
    println!();
    println!("τ = {}", format_complex(tau.tau()));
    println!("{}", describe_complex_lattice(&lattice));
    println!(
        "invariant truncation: r = {}",
        invariant_truncation.radius()
    );
    println!(
        "elliptic-function truncation: r = {}",
        function_truncation.radius()
    );
    println!(
        "verification tolerance: abs = {:.3e}, rel = {:.3e}",
        tolerance.absolute, tolerance.relative
    );
    println!();
    println!("Chosen z-points in the fundamental parallelogram:");
    println!();

    for z in points {
        show_point(
            &lattice,
            z,
            invariant_truncation,
            function_truncation,
            tolerance,
        )?;
    }

    let pole_report = verify_weierstrass_differential_equation(
        &lattice,
        Complex64::new(0.0, 0.0),
        invariant_truncation,
        function_truncation,
        tolerance,
    )?;

    println!("Pole case at z = 0:");
    println!("  values = {:?}", pole_report.values());
    println!("  verdict = {:?}", pole_report.status());

    if let TorusToCurveValues::Pole = pole_report.values() {
        println!("  interpretation: z ∈ Λ, so the map lands at the point at infinity.");
    }

    Ok(())
}
