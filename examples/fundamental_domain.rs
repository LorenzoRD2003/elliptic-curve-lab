use num_complex::Complex64;

use elliptic_algorithms_lab::elliptic_curves::analytic::{
    LatticeSumTruncation, UpperHalfPlanePoint,
};
use elliptic_algorithms_lab::numerics::ApproxTolerance;
use elliptic_algorithms_lab::visualization::Visualizable;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let tau = UpperHalfPlanePoint::new(Complex64::new(3.7, 0.2))?;
    let max_steps = 16;
    // NOTE: with radius = 1000 it is invariant under tolerance
    let truncation = LatticeSumTruncation::new(100)?;
    let tolerance = ApproxTolerance::new(1.0e-4, 1.0e-4);

    let reduction = tau.reduce_to_standard_fundamental_domain(max_steps)?;
    let modular_check = reduction.accumulated_matrix().verify_j_invariance_at(
        tau.clone(),
        truncation,
        tolerance,
    )?;

    println!("Reduction to the standard fundamental domain");
    println!("========================================================");
    println!();
    println!("original τ = {}", tau.tau().format_compact());
    println!("max steps = {max_steps}");
    println!("j-invariant truncation radius = {}", truncation.radius());
    println!(
        "comparison tolerance = abs {:.3e}, rel {:.3e}",
        tolerance.absolute, tolerance.relative
    );
    println!(
        "Note: the reduction itself is exact at the modular-action level, but the truncated analytic j-value is still quite sensitive for this messy τ."
    );
    println!();
    println!("{}", reduction.describe());
    println!();
    println!("Reduction steps");
    println!("---------------");
    if reduction.steps().is_empty() {
        println!("No modular step was needed.");
    } else {
        for (index, step) in reduction.steps().iter().enumerate() {
            println!("Step {}", index + 1);
            println!("{}", step.describe());
            println!();
        }
    }
    println!("Accumulated matrix");
    println!("------------------");
    println!("{}", reduction.accumulated_matrix().describe());
    println!();
    println!("j before and after reduction");
    println!("----------------------------");
    println!("{}", modular_check.describe());

    Ok(())
}
