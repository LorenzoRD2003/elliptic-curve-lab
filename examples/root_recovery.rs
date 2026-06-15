use num_complex::Complex64;

use elliptic_algorithms_lab::elliptic_curves::analytic::AnalyticWeierstrassCurve;
use elliptic_algorithms_lab::elliptic_curves::analytic::periods::{
    PeriodRecoveryConfig, WeierstrassCubicRoots,
};
use elliptic_algorithms_lab::fields::complex_approx::ComplexApprox;
use elliptic_algorithms_lab::numerics::ApproxTolerance;
use elliptic_algorithms_lab::visualization::{Visualizable, format_analytic_cubic_model};

fn c(re: f64, im: f64) -> Complex64 {
    Complex64::new(re, im)
}

fn yes_no(value: bool) -> &'static str {
    if value { "yes" } else { "no" }
}

fn indent(text: &str, spaces: usize) -> String {
    let prefix = " ".repeat(spaces);
    text.lines()
        .map(|line| format!("{prefix}{line}"))
        .collect::<Vec<_>>()
        .join("\n")
}

fn reconstruction_matches_curve(
    roots: &WeierstrassCubicRoots,
    curve: &AnalyticWeierstrassCurve,
    tolerance: ApproxTolerance,
) -> bool {
    ComplexApprox::eq_with_tolerance(&roots.g2(), curve.g2(), tolerance)
        && ComplexApprox::eq_with_tolerance(&roots.g3(), curve.g3(), tolerance)
}

fn print_root_recovery_case(
    title: &str,
    source_roots: WeierstrassCubicRoots,
    config: PeriodRecoveryConfig,
) -> Result<(), Box<dyn std::error::Error>> {
    let source_classification = source_roots.configuration_report(config.tolerance());
    let curve = AnalyticWeierstrassCurve::new(source_roots.g2(), source_roots.g3())?;
    let from_curve = curve.recover_weierstrass_cubic_roots(config)?;
    let from_invariants = WeierstrassCubicRoots::from_invariants(curve.g2(), curve.g3(), config)?;
    let report = curve.recover_weierstrass_cubic_roots_with_report(config)?;
    let recovered_classification = report
        .roots()
        .configuration_report(report.metadata().tolerance());

    println!("{title}");
    println!("{}", "=".repeat(title.len()));
    println!();
    println!("curve:");
    println!("  {}", format_analytic_cubic_model(&curve));
    println!();
    println!("source roots:");
    println!("{}", indent(&source_roots.describe(), 2));
    println!();
    println!("source classification:");
    println!("{}", indent(&source_classification.describe(), 2));
    println!();
    println!("recovered roots from curve:");
    println!("{}", indent(&from_curve.describe(), 2));
    println!();
    println!("recovered roots from invariants:");
    println!("{}", indent(&from_invariants.describe(), 2));
    println!();
    println!("recovery report:");
    println!("{}", indent(&report.describe(), 2));
    println!();
    println!("recovered classification:");
    println!("{}", indent(&recovered_classification.describe(), 2));
    println!();
    println!("metadata:");
    println!("{}", indent(&report.metadata().describe(), 2));
    println!();
    println!(
        "curve-level and invariant-level reconstruction agree with the curve? {} / {}",
        yes_no(reconstruction_matches_curve(
            &from_curve,
            &curve,
            config.tolerance()
        )),
        yes_no(reconstruction_matches_curve(
            &from_invariants,
            &curve,
            config.tolerance()
        ))
    );
    println!(
        "report says reconstructed coefficients agree with the curve? {}",
        yes_no(report.reconstruction_agrees())
    );
    println!();

    Ok(())
}

fn print_noisy_invariants_case(
    title: &str,
    source_roots: WeierstrassCubicRoots,
    g2_noise: Complex64,
    g3_noise: Complex64,
    config: PeriodRecoveryConfig,
) -> Result<(), Box<dyn std::error::Error>> {
    let exact_curve = AnalyticWeierstrassCurve::new(source_roots.g2(), source_roots.g3())?;
    let noisy_curve =
        AnalyticWeierstrassCurve::new(source_roots.g2() + g2_noise, source_roots.g3() + g3_noise)?;
    let report = noisy_curve.recover_weierstrass_cubic_roots_with_report(config)?;
    let strict_classification = report
        .roots()
        .configuration_report(report.metadata().tolerance());
    let loose_classification = report
        .roots()
        .configuration_report(ApproxTolerance::loose());

    println!("{title}");
    println!("{}", "=".repeat(title.len()));
    println!();
    println!("exact source curve:");
    println!("  {}", format_analytic_cubic_model(&exact_curve));
    println!("noisy inverse input:");
    println!("  {}", format_analytic_cubic_model(&noisy_curve));
    println!("  Δg₂ = {}", g2_noise);
    println!("  Δg₃ = {}", g3_noise);
    println!();
    println!("recovery report:");
    println!("{}", indent(&report.describe(), 2));
    println!();
    println!("classification after noisy recovery under strict tolerance:");
    println!("{}", indent(&strict_classification.describe(), 2));
    println!();
    println!("classification after noisy recovery under loose tolerance:");
    println!("{}", indent(&loose_classification.describe(), 2));
    println!();

    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let config = PeriodRecoveryConfig::strict();

    println!("Cubic-root recovery");
    println!("================================");
    println!();
    println!("active config:");
    println!("{}", indent(&config.describe(), 2));
    println!();
    println!(
        "All root triples below are shown in stored order; that order is implementation-stable but not mathematically canonical."
    );
    println!();

    print_root_recovery_case(
        "Case 1: three approximately real roots",
        WeierstrassCubicRoots::new(c(1.0, 0.0), c(2.0, 0.0), c(-3.0, 0.0), config.tolerance())?,
        config,
    )?;

    print_root_recovery_case(
        "Case 2: one real root plus a conjugate pair",
        WeierstrassCubicRoots::new(c(1.0, 1.0), c(-2.0, 0.0), c(1.0, -1.0), config.tolerance())?,
        config,
    )?;

    print_root_recovery_case(
        "Case 3: generic complex roots",
        WeierstrassCubicRoots::new(c(1.0, 1.0), c(-0.2, 0.1), c(-0.8, -1.1), config.tolerance())?,
        config,
    )?;

    print_root_recovery_case(
        "Case 4: larger complex roots where Newton helps",
        WeierstrassCubicRoots::new(
            c(8.813789020059971, -6.296193572032816),
            c(-5.70258988712044, -4.026550473696494),
            c(-3.1111991329395314, 10.32274404572931),
            config.tolerance(),
        )?,
        config,
    )?;

    print_noisy_invariants_case(
        "Case 5: inverse recovery from slightly noisy invariants",
        WeierstrassCubicRoots::new(c(1.0, 0.0), c(2.0, 0.0), c(-3.0, 0.0), config.tolerance())?,
        c(1.0e-8, -3.0e-9),
        c(-2.0e-8, 4.0e-9),
        config,
    )?;
    Ok(())
}
