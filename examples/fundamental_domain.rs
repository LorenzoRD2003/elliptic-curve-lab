use elliptic_algorithms_lab::elliptic_curves::analytic::verify_j_modular_invariance;
use elliptic_algorithms_lab::visualization::fields::format_complex;
use elliptic_algorithms_lab::visualization::{
    describe_fundamental_domain_reduction_report, describe_fundamental_domain_reduction_step,
    describe_modular_invariance_report, describe_modular_matrix,
};
use elliptic_algorithms_lab::{
    ApproxTolerance, LatticeSumTruncation, UpperHalfPlanePoint,
    reduce_tau_to_standard_fundamental_domain,
};
use num_complex::Complex64;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let tau = UpperHalfPlanePoint::new(Complex64::new(3.7, 0.2))?;
    let max_steps = 16;
    // NOTE: with radius = 1000 it is invariant under tolerance
    let truncation = LatticeSumTruncation::new(100)?;
    let tolerance = ApproxTolerance::new(1.0e-4, 1.0e-4);

    let reduction = reduce_tau_to_standard_fundamental_domain(tau.clone(), max_steps)?;
    let modular_check = verify_j_modular_invariance(
        tau.clone(),
        reduction.accumulated_matrix(),
        truncation,
        tolerance,
    )?;

    println!("Reduction to the standard fundamental domain");
    println!("========================================================");
    println!();
    println!("original τ = {}", format_complex(tau.tau()));
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
    println!(
        "{}",
        describe_fundamental_domain_reduction_report(&reduction)
    );
    println!();
    println!("Reduction steps");
    println!("---------------");
    if reduction.steps().is_empty() {
        println!("No modular step was needed.");
    } else {
        for (index, step) in reduction.steps().iter().enumerate() {
            println!("Step {}", index + 1);
            println!("{}", describe_fundamental_domain_reduction_step(step));
            println!();
        }
    }
    println!("Accumulated matrix");
    println!("------------------");
    println!(
        "{}",
        describe_modular_matrix(&reduction.accumulated_matrix())
    );
    println!();
    println!("j before and after reduction");
    println!("----------------------------");
    println!("{}", describe_modular_invariance_report(&modular_check));

    Ok(())
}
