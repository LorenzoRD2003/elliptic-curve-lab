use elliptic_algorithms_lab::fields::traits::*;
use num_bigint::BigUint;
use std::time::Instant;

use elliptic_algorithms_lab::elliptic_curves::{
    ShortWeierstrassCurve,
    frobenius::{HasseInterval, schoof::SchoofTraceModOddPrimeOutcome},
};
use elliptic_algorithms_lab::fields::{Fp1000000007, traits::FiniteField};
use elliptic_algorithms_lab::visualization::Visualizable;

type F = Fp1000000007;

fn heading(title: &str) {
    println!("{title}");
    println!("{}", "-".repeat(title.len()));
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let curve = ShortWeierstrassCurve::<F>::new(F::one(), F::from_i64(3))?;
    let started_at = Instant::now();
    let report = curve.schoof_group_order()?;
    let elapsed = started_at.elapsed();

    println!("Automatic Schoof over Fp<10^9 + 7>");
    println!("======================================================");
    println!();
    println!("Curve: {}", curve.format_compact());
    println!("Elapsed: {:.3?}", elapsed);
    println!();

    heading("Prime-by-prime progress");
    let threshold = BigUint::from(2u8) * HasseInterval::for_q(F::order())?.trace_bound();
    println!(
        "  start: mod 2 gives t ≡ {} (mod 2), Hasse threshold = {}",
        report.crt_report().mod_2_report().trace_mod_2(),
        threshold
    );
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
                    "  ℓ = {} skipped: non-unit denominator at candidate {} with gcd witness {}",
                    odd_prime_report.odd_prime(),
                    candidate_trace_mod_ell,
                    witness_gcd.format_compact()
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
    if let Some(solution) = report.crt_report().combined_solution() {
        println!(
            "  stop: CRT modulus {} is now above the Hasse threshold {}",
            solution.modulus(),
            threshold
        );
    }
    println!();

    heading("Detailed report");
    println!("{}", report.describe());
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

    Ok(())
}
