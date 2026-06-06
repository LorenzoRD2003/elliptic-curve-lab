use elliptic_algorithms_lab::visualization::fields::format_complex;
use elliptic_algorithms_lab::visualization::{Visualizable, format_analytic_cubic_model};
use elliptic_algorithms_lab::{
    AnalyticWeierstrassCurve, ComplexApprox, PeriodRecoveryConfig,
    recover_canonical_tau_from_curve, recover_period_basis, recover_tau_from_curve,
    reduce_tau_to_standard_fundamental_domain,
};
use num_complex::Complex64;

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

fn print_case(
    title: &str,
    curve: &AnalyticWeierstrassCurve,
    config: PeriodRecoveryConfig,
) -> Result<(), Box<dyn std::error::Error>> {
    let period_basis_report = recover_period_basis(curve, config)?;
    let tau_report = recover_tau_from_curve(curve, config)?;
    let canonical_tau_report = recover_canonical_tau_from_curve(curve, config)?;
    let matrix_image = canonical_tau_report
        .accumulated_matrix()
        .apply(&canonical_tau_report.original_tau())?;
    let one_step_config = PeriodRecoveryConfig::new(
        config.tolerance(),
        config.newton_max_iterations(),
        config.agm_max_iterations(),
        config.abel_jacobi_integration_steps(),
        config.branch_lattice_search_radius(),
        1,
    )?;
    let one_step_canonical = recover_canonical_tau_from_curve(curve, one_step_config)?;
    let zero_step_reduction =
        reduce_tau_to_standard_fundamental_domain(canonical_tau_report.original_tau(), 0)?;

    println!("{title}");
    println!("{}", "=".repeat(title.len()));
    println!();
    println!("curve:");
    println!("  {}", format_analytic_cubic_model(curve));
    println!();
    println!("full period-basis recovery report:");
    println!("{}", indent(&period_basis_report.describe(), 2));
    println!();
    println!("Legendre-to-period-basis transport summary:");
    println!(
        "  {}",
        period_basis_report
            .legendre_reduction()
            .parameter()
            .format_compact()
            .replace("lambda", "λ")
    );
    println!(
        "  invariant differential scale ≈ {}",
        format_complex(
            &period_basis_report
                .basis_report()
                .invariant_differential_scale()
        )
    );
    println!(
        "  {}",
        period_basis_report
            .periods()
            .describe()
            .replace('\n', "\n  ")
    );
    println!();
    println!("natural τ summary:");
    println!("  τ ≈ {}", format_complex(tau_report.tau().tau()));
    println!("  this is the τ induced directly by the recovered basis (ω₁, ω₂)");
    println!();
    println!("canonical τ summary:");
    println!(
        "  canonical τ ≈ {}",
        format_complex(canonical_tau_report.canonical_tau().tau())
    );
    println!(
        "  γ = [[{}, {}], [{}, {}]]",
        canonical_tau_report.accumulated_matrix().a(),
        canonical_tau_report.accumulated_matrix().b(),
        canonical_tau_report.accumulated_matrix().c(),
        canonical_tau_report.accumulated_matrix().d()
    );
    println!(
        "  canonical τ = γ(original τ)? {}",
        yes_no(ComplexApprox::eq_with_tolerance(
            matrix_image.tau(),
            canonical_tau_report.canonical_tau().tau(),
            config.tolerance(),
        ))
    );
    println!(
        "  reduction status = {:?}",
        canonical_tau_report.fundamental_domain_reduction().status()
    );
    println!();
    println!("consistency checks:");
    println!(
        "  τ from full report agrees with τ from τ-wrapper? {}",
        yes_no(ComplexApprox::eq_with_tolerance(
            period_basis_report.tau().tau(),
            tau_report.tau().tau(),
            config.tolerance(),
        ))
    );
    println!(
        "  natural τ matches the original τ stored in the canonical report? {}",
        yes_no(ComplexApprox::eq_with_tolerance(
            tau_report.tau().tau(),
            canonical_tau_report.original_tau().tau(),
            config.tolerance(),
        ))
    );
    println!(
        "  canonical τ lies in the standard fundamental domain? {}",
        yes_no(
            canonical_tau_report
                .fundamental_domain_reduction()
                .is_reduced()
        )
    );
    println!();
    println!(
        "  caveat: the recovered basis (ω₁, ω₂) and the resulting natural τ are valid before any canonical SL₂(ℤ)-normalization."
    );
    println!();
    println!("tiny-budget modular reduction demo:");
    println!(
        "  with fundamental-domain reduction budget = 1, canonical recovery still succeeds here? {}",
        yes_no(
            one_step_canonical
                .fundamental_domain_reduction()
                .is_reduced()
        )
    );
    println!(
        "  steps actually used in that one-step run = {}",
        one_step_canonical
            .fundamental_domain_reduction()
            .steps()
            .len()
    );
    println!(
        "  if we allow 0 direct modular steps after recovering natural τ, the reduction status becomes {:?}",
        zero_step_reduction.status()
    );
    println!();

    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let config = PeriodRecoveryConfig::strict();

    println!("End-to-end period recovery from a curve");
    println!("====================================================");
    println!();
    println!("active config:");
    println!("{}", indent(&config.describe(), 2));
    println!();
    println!("The recovery pipeline below starts from an analytic curve y² = 4x³ - g₂x - g₃.");
    println!(
        "It then recovers cubic roots, chooses a Legendre representative, computes K(λ) and K(1-λ), transports periods back to the curve, and exposes the resulting τ."
    );
    println!();

    let real_split_curve = AnalyticWeierstrassCurve::new(c(28.0, 0.0), c(-24.0, 0.0))?;
    print_case("Case 1: real split cubic", &real_split_curve, config)?;

    let harder_complex_curve = AnalyticWeierstrassCurve::new(
        c(-85.06447698350209, -258.5906723468124),
        c(923.8388816938035, -3127.3197765358514),
    )?;
    print_case(
        "Case 2: harder generic complex curve",
        &harder_complex_curve,
        config,
    )?;

    Ok(())
}
