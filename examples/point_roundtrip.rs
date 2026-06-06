use elliptic_algorithms_lab::elliptic_curves::analytic::{
    PointRoundTripValidationConfig, validate_point_inverse_uniformization_roundtrip_with_periods,
};
use elliptic_algorithms_lab::{
    AbelJacobiConfig, AnalyticCurvePoint, AnalyticWeierstrassCurve, ApproxTolerance, ComplexApprox,
    ComplexLattice, EllipticFunctionTruncation, LatticeSumTruncation, PeriodRecoveryConfig,
    UpperHalfPlanePoint, Visualizable, format_analytic_cubic_model, format_complex,
    format_point_compact, map_torus_point_to_curve, recover_period_basis,
    reduce_tau_to_standard_fundamental_domain,
};
use num_complex::Complex64;

fn c(re: f64, im: f64) -> Complex64 {
    Complex64::new(re, im)
}

fn yes_no(value: bool) -> &'static str {
    if value { "yes" } else { "no" }
}

fn print_prefixed_block(prefix: &str, text: &str) {
    for line in text.lines() {
        println!("{prefix}{line}");
    }
}

#[derive(Clone)]
struct PointRoundTripCase {
    title: &'static str,
    tau: UpperHalfPlanePoint,
    source_z: Complex64,
    validation_config: PointRoundTripValidationConfig,
    expectation: &'static str,
}

struct PointRoundTripOutcome {
    title: &'static str,
    agrees_approximately: bool,
    same_torus_class_as_source_approximately: bool,
    x_residual_norm: f64,
    y_residual_norm: f64,
}

fn build_source_point(
    lattice: &ComplexLattice,
    source_z: Complex64,
) -> Result<AnalyticCurvePoint, Box<dyn std::error::Error>> {
    let point = map_torus_point_to_curve(
        lattice,
        source_z,
        LatticeSumTruncation::new(18)?,
        EllipticFunctionTruncation::new(16)?,
        ApproxTolerance::strict(),
    )?
    .point()
    .clone();

    Ok(point)
}

fn print_case(
    case: &PointRoundTripCase,
    period_config: PeriodRecoveryConfig,
    abel_config: AbelJacobiConfig,
) -> Result<PointRoundTripOutcome, Box<dyn std::error::Error>> {
    let torus_tolerance = ApproxTolerance::new(1.0e-2, 1.0e-2);
    let lattice = ComplexLattice::from_tau(case.tau.clone());
    let curve = AnalyticWeierstrassCurve::from_tau(&case.tau, LatticeSumTruncation::new(18)?)?;
    let source_point = build_source_point(&lattice, case.source_z)?;
    let period_basis_report = recover_period_basis(&curve, period_config)?;
    let source_canonical_tau = reduce_tau_to_standard_fundamental_domain(
        case.tau.clone(),
        period_config.fundamental_domain_reduction_max_steps(),
    )?;
    let recovered_canonical_tau = reduce_tau_to_standard_fundamental_domain(
        period_basis_report.tau(),
        period_config.fundamental_domain_reduction_max_steps(),
    )?;
    let report = validate_point_inverse_uniformization_roundtrip_with_periods(
        &curve,
        &source_point,
        period_basis_report.periods(),
        abel_config,
        case.validation_config,
    )?;
    let torus_comparison = period_basis_report
        .periods()
        .lattice()
        .compare_complex_points_mod_lattice_approx(
            *report.reduced_representative(),
            case.source_z,
            1,
            torus_tolerance,
        )?;

    println!("{}", case.title);
    println!("{}", "=".repeat(case.title.len()));
    println!("curve:");
    println!("  {}", format_analytic_cubic_model(&curve));
    println!("ambient lattice comparison:");
    println!("  source τ ≈ {}", format_complex(case.tau.tau()));
    println!(
        "  recovered ω₁ ≈ {}, recovered ω₂ ≈ {}",
        format_complex(period_basis_report.periods().omega1()),
        format_complex(period_basis_report.periods().omega2())
    );
    println!(
        "  τ_recovered ≈ {}",
        format_complex(period_basis_report.tau().tau())
    );
    println!(
        "  canonical τ(source) ≈ {}, canonical τ(recovered) ≈ {}",
        format_complex(source_canonical_tau.reduced_tau().tau()),
        format_complex(recovered_canonical_tau.reduced_tau().tau())
    );
    println!(
        "  same modular class after canonicalization? {}",
        yes_no(ComplexApprox::eq_with_tolerance(
            source_canonical_tau.reduced_tau().tau(),
            recovered_canonical_tau.reduced_tau().tau(),
            ApproxTolerance::new(1.0e-6, 1.0e-6),
        ))
    );
    println!(
        "  covolume ratio recovered/source ≈ {:.6e}",
        period_basis_report.periods().covolume() / lattice.covolume()
    );
    println!("source torus seed and point:");
    println!("  source z ≈ {}", format_complex(&case.source_z));
    println!("  source point P = {}", format_point_compact(&source_point));
    println!("inverse-uniformization result:");
    println!(
        "  recovered z_P ≈ {}",
        format_complex(report.reduced_representative())
    );
    println!(
        "  |z_P - z_source| ≈ {:.6e}",
        (*report.reduced_representative() - case.source_z).norm()
    );
    println!(
        "  same torus class as source z approximately? {}",
        yes_no(torus_comparison.agrees_approximately())
    );
    println!(
        "  best lattice shift = ({}, {}), shifted residual ≈ {:.6e}",
        torus_comparison.best_shift().m,
        torus_comparison.best_shift().n,
        torus_comparison.shifted_difference_norm()
    );
    println!(
        "  contour = {:?}, tail length ≈ {:.6e}",
        report
            .point_recovery_report()
            .contour()
            .legendre_contour_strategy(),
        report.point_recovery_report().contour().tail_length()
    );
    println!(
        "  numerics: {} Simpson steps, {} branch adjustments",
        report
            .point_recovery_report()
            .metadata()
            .integration_steps_used(),
        report
            .point_recovery_report()
            .metadata()
            .branch_adjustments_used()
    );
    println!(
        "  recovered point via (wp, wp') = {}",
        format_point_compact(report.recovered_curve_point())
    );
    println!("validation:");
    println!("  x residual norm = {:.6e}", report.x_residual_norm());
    println!("  y residual norm = {:.6e}", report.y_residual_norm());
    println!(
        "  final verdict = {}",
        if report.agrees_approximately() {
            "agrees"
        } else {
            "does not agree"
        }
    );
    println!("  expectation = {}", case.expectation);
    println!(
        "  validation config = {}",
        case.validation_config.format_compact()
    );
    println!(
        "  torus-class comparison tolerance = abs {:.1e}, rel {:.1e}",
        torus_tolerance.absolute, torus_tolerance.relative
    );
    println!();

    Ok(PointRoundTripOutcome {
        title: case.title,
        agrees_approximately: report.agrees_approximately(),
        same_torus_class_as_source_approximately: torus_comparison.agrees_approximately(),
        x_residual_norm: report.x_residual_norm(),
        y_residual_norm: report.y_residual_norm(),
    })
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let period_config = PeriodRecoveryConfig::strict();
    let abel_config = AbelJacobiConfig::strict();

    println!("Point roundtrip validation");
    println!("=======================================");
    println!();
    println!("active configs:");
    print_prefixed_block("  ", &period_config.describe());
    println!(
        "  Abel-Jacobi config: contour = {:?}, steps = {}, segment samples = {}, ray samples = {}, max branch adjustments = {}",
        abel_config.legendre_contour_strategy,
        abel_config.integration_steps,
        abel_config.segment_samples,
        abel_config.ray_samples,
        abel_config.max_branch_adjustments
    );
    println!();
    println!("Each case starts from a known torus sample z, maps it forward to a curve point P,");
    println!("then runs the inverse experiment P -> z_P mod Λ -> (wp(z_P), wp'(z_P)).");
    println!();

    let cases = vec![
        PointRoundTripCase {
            title: "Case 1: square lattice, generic interior point",
            tau: UpperHalfPlanePoint::tau_i(),
            source_z: c(0.20, 0.15),
            validation_config: PointRoundTripValidationConfig::new(
                LatticeSumTruncation::new(24)?,
                EllipticFunctionTruncation::new(22)?,
                ApproxTolerance::new(5.0e-3, 5.0e-3),
            ),
            expectation: "this is the cleanest sanity check: both the torus class and the final roundtrip should look healthy",
        },
        PointRoundTripCase {
            title: "Case 2: square lattice, point near a half-period",
            tau: UpperHalfPlanePoint::tau_i(),
            source_z: c(0.49, 0.02),
            validation_config: PointRoundTripValidationConfig::new(
                LatticeSumTruncation::new(24)?,
                EllipticFunctionTruncation::new(22)?,
                ApproxTolerance::new(5.0e-3, 5.0e-3),
            ),
            expectation: "this keeps the same tolerance as the easy case, but the geometry near a semiperiod should make the forward check noticeably harder",
        },
        PointRoundTripCase {
            title: "Case 3: generic non-CM lattice, generic interior point",
            tau: UpperHalfPlanePoint::tau_generic_example(),
            source_z: c(0.18, 0.11),
            validation_config: PointRoundTripValidationConfig::new(
                LatticeSumTruncation::new(24)?,
                EllipticFunctionTruncation::new(22)?,
                ApproxTolerance::new(5.0e-3, 5.0e-3),
            ),
            expectation: "this shows that the same pedagogical tolerance can still certify a generic non-CM case",
        },
        PointRoundTripCase {
            title: "Case 4: same generic lattice, near a semiperiod, with deliberately weak forward validation",
            tau: UpperHalfPlanePoint::tau_generic_example(),
            source_z: c(0.49, 0.02),
            validation_config: PointRoundTripValidationConfig::new(
                LatticeSumTruncation::new(1)?,
                EllipticFunctionTruncation::new(1)?,
                ApproxTolerance::new(5.0e-3, 5.0e-3),
            ),
            expectation: "the torus recovery may still look reasonable, but now the same tolerance should fail because the forward validation budget is intentionally too small for this harder geometry",
        },
    ];

    let mut outcomes = Vec::new();

    for case in &cases {
        outcomes.push(print_case(case, period_config, abel_config)?);
    }

    println!("Summary");
    println!("=======");
    for outcome in &outcomes {
        println!(
            "- {}: verdict = {}, same torus class = {}, x residual = {:.3e}, y residual = {:.3e}",
            outcome.title,
            if outcome.agrees_approximately {
                "agrees"
            } else {
                "does not agree"
            },
            yes_no(outcome.same_torus_class_as_source_approximately),
            outcome.x_residual_norm,
            outcome.y_residual_norm
        );
    }

    Ok(())
}
