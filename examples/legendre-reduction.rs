use elliptic_algorithms_lab::{
    AnalyticWeierstrassCurve, ApproxTolerance, ComplexApprox, LegendreReduction, Visualizable,
    WeierstrassCubicRoots, classify_legendre_parameter_conditioning,
    cubic_root_configuration_report, format_analytic_cubic_model, legendre_reduction_report,
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

fn with_lambda_symbol(text: &str) -> String {
    text.replace("lambda", "λ")
}

fn approx_eq(left: &Complex64, right: &Complex64, tolerance: ApproxTolerance) -> bool {
    ComplexApprox::eq_with_tolerance(left, right, tolerance)
}

fn evaluate_curve_rhs(curve: &AnalyticWeierstrassCurve, x: Complex64) -> Complex64 {
    Complex64::new(4.0, 0.0) * x.powu(3) - *curve.g2() * x - *curve.g3()
}

fn print_legendre_sanity_checks(
    reduction: &LegendreReduction,
    curve: &AnalyticWeierstrassCurve,
    tolerance: ApproxTolerance,
) {
    let [e1, e2, e3] = reduction.selected_root_triple();
    let legendre_x_e1 = reduction.legendre_x_from_original_x(*e1);
    let legendre_x_e2 = reduction.legendre_x_from_original_x(*e2);
    let legendre_x_e3 = reduction.legendre_x_from_original_x(*e3);
    let sample_legendre_x = c(0.375, -0.25);
    let sample_original_x = reduction.original_x_from_legendre_x(sample_legendre_x);
    let direct_curve_rhs = evaluate_curve_rhs(curve, sample_original_x);
    let reduced_curve_rhs = reduction.evaluate_original_cubic_from_legendre_x(sample_legendre_x);
    let y_scale_squared = reduction.legendre_y_scale().powu(2);
    let rhs_scale = reduction.legendre_rhs_scale_factor();
    let differential_identity =
        reduction.invariant_differential_scale() * reduction.legendre_y_scale();

    println!("sanity checks:");
    println!(
        "  selected roots map to X = 1, 0, λ? {} / {} / {}",
        yes_no(approx_eq(
            &legendre_x_e1,
            &Complex64::new(1.0, 0.0),
            tolerance,
        )),
        yes_no(approx_eq(
            &legendre_x_e2,
            &Complex64::new(0.0, 0.0),
            tolerance,
        )),
        yes_no(approx_eq(
            &legendre_x_e3,
            reduction.parameter().lambda(),
            tolerance,
        ))
    );
    println!(
        "  principal y-scale squares to rhs scale factor? {}",
        yes_no(approx_eq(&y_scale_squared, &rhs_scale, tolerance))
    );
    println!(
        "  differential scale times y-scale recovers x-scale? {}",
        yes_no(approx_eq(
            &differential_identity,
            &reduction.x_scale(),
            tolerance
        ))
    );
    println!(
        "  original cubic agrees with scaled Legendre cubic at sample X? {}",
        yes_no(approx_eq(&direct_curve_rhs, &reduced_curve_rhs, tolerance))
    );
}

fn print_case(
    title: &str,
    roots: WeierstrassCubicRoots,
    tolerance: ApproxTolerance,
) -> Result<(), Box<dyn std::error::Error>> {
    let curve = AnalyticWeierstrassCurve::new(roots.g2(), roots.g3())?;
    let root_configuration = cubic_root_configuration_report(&roots, tolerance);
    let reduction = LegendreReduction::from_roots(&roots, tolerance)?;
    let report = legendre_reduction_report(&roots, tolerance)?;
    let direct_conditioning =
        classify_legendre_parameter_conditioning(reduction.parameter(), tolerance);
    let loose_conditioning =
        classify_legendre_parameter_conditioning(reduction.parameter(), ApproxTolerance::loose());

    println!("{title}");
    println!("{}", "=".repeat(title.len()));
    println!();
    println!("curve:");
    println!("  {}", format_analytic_cubic_model(&curve));
    println!();
    println!("source roots:");
    println!("{}", indent(&roots.describe(), 2));
    println!();
    println!("cubic-root configuration:");
    println!("{}", indent(&root_configuration.describe(), 2));
    println!();
    println!("chosen Legendre parameter:");
    println!(
        "{}",
        indent(&with_lambda_symbol(&reduction.parameter().describe()), 2)
    );
    println!();
    println!("full S3 orbit:");
    println!(
        "{}",
        indent(&with_lambda_symbol(&reduction.orbit().describe()), 2)
    );
    println!();
    println!("explicit reduction:");
    println!("{}", indent(&with_lambda_symbol(&reduction.describe()), 2));
    println!();
    println!("reduction report:");
    println!("{}", indent(&with_lambda_symbol(&report.describe()), 2));
    println!();
    println!(
        "conditioning from direct classifier / report: {:?} / {:?}",
        direct_conditioning,
        report.conditioning()
    );
    println!(
        "conditioning under loose tolerance: {:?}",
        loose_conditioning
    );
    println!(
        "near singular locus according to report? {}",
        yes_no(report.is_near_singular())
    );
    println!(
        "singularity distance score = {:.6e}",
        report.singularity_distance()
    );
    println!();
    print_legendre_sanity_checks(&reduction, &curve, tolerance);
    println!();

    Ok(())
}

fn print_permutation_invariance_case(
    title: &str,
    original: WeierstrassCubicRoots,
    permuted: WeierstrassCubicRoots,
    tolerance: ApproxTolerance,
) -> Result<(), Box<dyn std::error::Error>> {
    let original_reduction = LegendreReduction::from_roots(&original, tolerance)?;
    let permuted_reduction = LegendreReduction::from_roots(&permuted, tolerance)?;
    let original_report = legendre_reduction_report(&original, tolerance)?;
    let permuted_report = legendre_reduction_report(&permuted, tolerance)?;

    println!("{title}");
    println!("{}", "=".repeat(title.len()));
    println!();
    println!("original stored order:");
    println!("{}", indent(&original.describe(), 2));
    println!();
    println!("permuted stored order:");
    println!("{}", indent(&permuted.describe(), 2));
    println!();
    println!("what changes under permutation:");
    println!(
        "  input-order-relative orbit label: {:?} -> {:?}",
        original_report.selected_orbit_element_relative_to_input_order(),
        permuted_report.selected_orbit_element_relative_to_input_order()
    );
    println!(
        "  selected permutation relative to caller order: {:?} -> {:?}",
        original_reduction.selected_permutation(),
        permuted_reduction.selected_permutation()
    );
    println!();
    println!("what does not change:");
    println!(
        "  chosen lambda agrees? {}",
        yes_no(approx_eq(
            original_reduction.parameter().lambda(),
            permuted_reduction.parameter().lambda(),
            tolerance,
        ))
    );
    println!(
        "  conditioning agrees? {}",
        yes_no(original_report.conditioning() == permuted_report.conditioning())
    );
    println!(
        "  singularity distance agrees? {}",
        yes_no(tolerance.real_close(
            original_report.singularity_distance(),
            permuted_report.singularity_distance()
        ))
    );
    println!(
        "  full Legendre orbit agrees? {}",
        yes_no(
            original_reduction
                .orbit()
                .values()
                .into_iter()
                .zip(permuted_reduction.orbit().values())
                .all(|(lhs, rhs)| approx_eq(&lhs, &rhs, tolerance))
        )
    );
    println!();

    Ok(())
}

fn print_controlled_rejection_case(
    title: &str,
    tolerance: ApproxTolerance,
) -> Result<(), Box<dyn std::error::Error>> {
    println!("{title}");
    println!("{}", "=".repeat(title.len()));
    println!();
    println!("attempt 1: repeated-root input at the root-constructor boundary");
    match WeierstrassCubicRoots::new(
        c(1.0, 0.0),
        c(1.0 + 5.0e-13, 0.0),
        c(-2.0 - 5.0e-13, 0.0),
        tolerance,
    ) {
        Ok(_) => println!("  unexpected success"),
        Err(error) => {
            println!("  rejected with: {error}");
            println!(
                "  why this is mathematically sane: Legendre reduction would divide by e1 - e2, and here that denominator is within the chosen repeated-root tolerance."
            );
        }
    }
    println!();
    println!("attempt 2: exactly singular analytic invariants");
    match AnalyticWeierstrassCurve::new(c(12.0, 0.0), c(-8.0, 0.0)) {
        Ok(_) => println!("  unexpected success"),
        Err(error) => {
            println!("  rejected with: {error}");
            println!(
                "  why this is mathematically sane: these invariants come from 4(x-1)^2(x+2), whose discriminant is zero, so there is no nonsingular elliptic curve or honest Legendre parameter behind them."
            );
        }
    }
    println!();

    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let tolerance = ApproxTolerance::strict();

    println!("Milestone 9: Legendre reduction");
    println!("===============================");
    println!();
    println!(
        "This example starts from stored cubic roots for 4(x-e1)(x-e2)(x-e3), chooses one deterministic Legendre representative, and then checks the affine normalization numerically."
    );
    println!();
    println!(
        "All stored root orders below are implementation-level input orderings, not canonical mathematical orderings."
    );
    println!();

    print_case(
        "Case 1: three approximately real roots",
        WeierstrassCubicRoots::new(c(1.0, 0.0), c(2.0, 0.0), c(-3.0, 0.0), tolerance)?,
        tolerance,
    )?;

    print_case(
        "Case 2: one approximately real root plus a conjugate pair",
        WeierstrassCubicRoots::new(c(1.0, 1.0), c(-2.0, 0.0), c(1.0, -1.0), tolerance)?,
        tolerance,
    )?;

    print_case(
        "Case 3: generic complex roots",
        WeierstrassCubicRoots::new(c(1.0, 1.0), c(-0.2, 0.1), c(-0.8, -1.1), tolerance)?,
        tolerance,
    )?;

    print_case(
        "Case 4: nearly singular but still distinct roots",
        WeierstrassCubicRoots::new(
            c(1.0, 0.0),
            c(1.0 + 1.0e-7, 0.0),
            c(-2.0 - 1.0e-7, 0.0),
            tolerance,
        )?,
        tolerance,
    )?;

    print_permutation_invariance_case(
        "Case 5: what changes and what does not under root permutation",
        WeierstrassCubicRoots::new(c(1.0, 0.0), c(2.0, 0.0), c(-3.0, 0.0), tolerance)?,
        WeierstrassCubicRoots::new(c(2.0, 0.0), c(-3.0, 0.0), c(1.0, 0.0), tolerance)?,
        tolerance,
    )?;

    print_controlled_rejection_case(
        "Case 6: controlled rejection near the singular boundary",
        tolerance,
    )?;

    Ok(())
}
