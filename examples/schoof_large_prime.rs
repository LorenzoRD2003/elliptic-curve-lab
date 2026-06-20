use elliptic_algorithms_lab::elliptic_curves::{
    ShortWeierstrassCurve,
    frobenius::{
        group_order::FiniteFieldGroupOrderStrategy, schoof::SchoofTraceModOddPrimeOutcome,
    },
};
use elliptic_algorithms_lab::fields::{Fp, traits::Field};
use elliptic_algorithms_lab::visualization::{Visualizable, format_curve};

type F = Fp<1_000_000_007>;

fn heading(title: &str) {
    println!("{title}");
    println!("{}", "-".repeat(title.len()));
}

enum ExampleMode {
    DefaultCurve,
    ExplicitCurve(ShortWeierstrassCurve<F>),
}

fn parse_mode() -> Result<ExampleMode, Box<dyn std::error::Error>> {
    let mut args = std::env::args().skip(1);
    let Some(first) = args.next() else {
        return Ok(ExampleMode::DefaultCurve);
    };

    let Some(second) = args.next() else {
        return Err("expected both a and b when passing explicit curve coefficients".into());
    };
    if args.next().is_some() {
        return Err("expected at most two positional arguments: a b".into());
    }

    let a: i64 = first.parse()?;
    let b: i64 = second.parse()?;
    Ok(ExampleMode::ExplicitCurve(ShortWeierstrassCurve::<F>::new(
        F::from_i64(a),
        F::from_i64(b),
    )?))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let curve = match parse_mode()? {
        ExampleMode::DefaultCurve => {
            ShortWeierstrassCurve::<F>::new(F::from_i64(-16), F::from_i64(-16))?
        }
        ExampleMode::ExplicitCurve(curve) => curve,
    };
    let report = curve.schoof_group_order()?;

    println!("Automatic Schoof over Fp<10^9 + 7>");
    println!("======================================================");
    println!();
    println!("Curve: {}", format_curve(&curve));
    println!();

    heading("Detailed report");
    println!("{}", report.describe());
    println!();

    heading("Odd-prime step log");
    for odd_prime_report in report.crt_report().odd_prime_reports() {
        match odd_prime_report.outcome() {
            SchoofTraceModOddPrimeOutcome::TraceFound { trace_mod_ell } => {
                println!(
                    "  ℓ = {} resolved: t ≡ {} (mod {})",
                    odd_prime_report.odd_prime(),
                    trace_mod_ell,
                    odd_prime_report.odd_prime()
                );
            }
            SchoofTraceModOddPrimeOutcome::NonUnitDenominator {
                candidate_trace_mod_ell,
                witness_gcd,
            } => {
                println!(
                    "  ℓ = {} skipped: non-unit denominator at candidate {} with gcd witness {:?}",
                    odd_prime_report.odd_prime(),
                    candidate_trace_mod_ell,
                    witness_gcd
                );
            }
            SchoofTraceModOddPrimeOutcome::ExhaustedCandidates => {
                println!(
                    "  ℓ = {} gave no usable trace residue before refinement",
                    odd_prime_report.odd_prime()
                );
            }
        }
    }
    println!();

    heading("Summary");
    match (
        report.trace(),
        report.curve_order(),
        report.crt_report().combined_solution(),
    ) {
        (Some(trace), Some(curve_order), Some(solution)) => {
            println!("  resolved trace t = {trace}");
            println!("  resolved group order #E(F_q) = {curve_order}");
            println!(
                "  final CRT class: t ≡ {} (mod {})",
                solution.residue(),
                solution.modulus()
            );
        }
        _ => {
            println!("  automatic Schoof did not resolve the final trace");
        }
    }

    let skipped = report
        .crt_report()
        .odd_prime_reports()
        .iter()
        .filter(|odd_prime_report| {
            matches!(
                odd_prime_report.outcome(),
                SchoofTraceModOddPrimeOutcome::NonUnitDenominator { .. }
            )
        })
        .count();
    println!("  skipped odd primes due to non-unit denominators: {skipped}");

    let integrated = curve.group_order_by(FiniteFieldGroupOrderStrategy::Schoof)?;
    println!(
        "  integrated group-order route: {}",
        integrated.format_compact()
    );

    Ok(())
}
