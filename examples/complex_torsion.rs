use elliptic_algorithms_lab::{
    AnalyticCurvePoint, AnalyticDivisionPolynomialComparisonCase, AnalyticWeierstrassCurve,
    ApproxTolerance, ComplexLattice, EllipticFunctionApproximation, EllipticFunctionTruncation,
    LatticeSumTruncation, UpperHalfPlanePoint, analytic_invariants,
    compare_primitive_analytic_torsion_with_division_polynomial,
    describe_analytic_division_polynomial_comparison, describe_analytic_invariants,
    describe_analytic_torsion_point_approx, describe_complex_lattice, format_analytic_cubic_model,
    format_complex, format_short_weierstrass_over_complex, map_primitive_torus_torsion_to_curve,
    map_torus_point_to_curve, weierstrass_p, weierstrass_p_derivative,
};
use num_complex::Complex64;

fn eval_monic_cubic(x: Complex64, a1: Complex64, a0: Complex64) -> Complex64 {
    x * x * x + a1 * x + a0
}

fn durand_kerner_monic_cubic(a1: Complex64, a0: Complex64) -> [Complex64; 3] {
    let mut roots = [
        Complex64::new(1.0, 0.0),
        Complex64::new(-0.5, 0.866_025_403_784_438_6),
        Complex64::new(-0.5, -0.866_025_403_784_438_6),
    ];

    for _ in 0..100 {
        let old = roots;
        let mut max_delta = 0.0_f64;

        for i in 0..3 {
            let mut denominator = Complex64::new(1.0, 0.0);
            for j in 0..3 {
                if i != j {
                    denominator *= old[i] - old[j];
                }
            }

            let delta = eval_monic_cubic(old[i], a1, a0) / denominator;
            roots[i] = old[i] - delta;
            max_delta = max_delta.max(delta.norm());
        }

        if max_delta < 1.0e-12 {
            break;
        }
    }

    roots
}

fn primitive_point_z(
    lattice: &ComplexLattice,
    n: usize,
    a: usize,
    b: usize,
) -> Result<Complex64, Box<dyn std::error::Error>> {
    let points = elliptic_algorithms_lab::primitive_torus_n_torsion_points(lattice, n)?;
    let point = points
        .into_iter()
        .find(|point| point.index().a() == a && point.index().b() == b)
        .ok_or_else(|| format!("missing primitive torus point ({a}, {b}; {n})"))?;
    Ok(*point.z())
}

fn primitive_mapped_point(
    lattice: &ComplexLattice,
    n: usize,
    a: usize,
    b: usize,
    invariant_truncation: LatticeSumTruncation,
    function_truncation: EllipticFunctionTruncation,
    tolerance: ApproxTolerance,
) -> Result<elliptic_algorithms_lab::AnalyticTorsionPointApprox, Box<dyn std::error::Error>> {
    let points = map_primitive_torus_torsion_to_curve(
        lattice,
        n,
        invariant_truncation,
        function_truncation,
        tolerance,
    )?;
    points
        .into_iter()
        .find(|point| point.torus_point().index().a() == a && point.torus_point().index().b() == b)
        .ok_or_else(|| format!("missing mapped primitive point ({a}, {b}; {n})").into())
}

fn primitive_division_report(
    lattice: &ComplexLattice,
    n: usize,
    a: usize,
    b: usize,
    invariant_truncation: LatticeSumTruncation,
    function_truncation: EllipticFunctionTruncation,
    tolerance: ApproxTolerance,
) -> Result<AnalyticDivisionPolynomialComparisonCase, Box<dyn std::error::Error>> {
    let reports = compare_primitive_analytic_torsion_with_division_polynomial(
        lattice,
        n,
        invariant_truncation,
        function_truncation,
        tolerance,
    )?;
    reports
        .into_iter()
        .find(|report| {
            let index = report.torsion_point().torus_point().index();
            index.a() == a && index.b() == b
        })
        .ok_or_else(|| format!("missing division report for ({a}, {b}; {n})").into())
}

fn point_distance(left: &AnalyticCurvePoint, right: &AnalyticCurvePoint) -> Option<f64> {
    match (left, right) {
        (AnalyticCurvePoint::Infinity, AnalyticCurvePoint::Infinity) => Some(0.0),
        (
            AnalyticCurvePoint::Finite { x: x1, y: y1 },
            AnalyticCurvePoint::Finite { x: x2, y: y2 },
        ) => {
            let dx = (*x1 - *x2).norm();
            let dy = (*y1 - *y2).norm();
            Some((dx * dx + dy * dy).sqrt())
        }
        _ => None,
    }
}

fn print_n2_experiment(
    lattice: &ComplexLattice,
    invariant_truncation: LatticeSumTruncation,
    tolerance: ApproxTolerance,
) -> Result<(), Box<dyn std::error::Error>> {
    let invariants = analytic_invariants(lattice, invariant_truncation)?;
    let a1 = -invariants.g2 / Complex64::new(4.0, 0.0);
    let a0 = -invariants.g3 / Complex64::new(4.0, 0.0);
    let cubic_roots = durand_kerner_monic_cubic(a1, a0);

    println!("n = 2");
    println!("-----");
    println!(
        "For non-trivial 2-torsion, the useful signal is ℘′(z) ≈ 0 and x = ℘(z) near a root of 4x^3 - g₂x - g₃."
    );
    println!("Approximate cubic roots:");
    for (i, root) in cubic_roots.iter().enumerate() {
        println!("  root[{i}] ≈ {root}");
    }
    println!();

    for &(a, b) in &[(0_usize, 1_usize), (1, 0), (1, 1)] {
        let z = primitive_point_z(lattice, 2, a, b)?;
        println!("point ({a}, {b}; 2), z = {}", format_complex(&z));

        for radius in [6_usize, 10, 14] {
            let function_truncation = EllipticFunctionTruncation::new(radius)?;
            let mapped = primitive_mapped_point(
                lattice,
                2,
                a,
                b,
                invariant_truncation,
                function_truncation,
                tolerance,
            )?;
            let p_prime = weierstrass_p_derivative(lattice, z, function_truncation)?;

            match mapped.curve_point() {
                AnalyticCurvePoint::Finite { x, .. } => {
                    let nearest_root_distance = cubic_roots
                        .iter()
                        .map(|root| (*x - *root).norm())
                        .fold(f64::INFINITY, f64::min);

                    println!(
                        "  r_fun = {:>2}: |℘′(z)| = {:.6e}, curve residual = {:.6e}, distance to nearest cubic root = {:.6e}",
                        radius,
                        p_prime.value().norm(),
                        mapped.membership_report().absolute_error(),
                        nearest_root_distance
                    );
                }
                AnalyticCurvePoint::Infinity => {
                    println!("  r_fun = {:>2}: unexpected pole", radius);
                }
            }
        }

        let report = primitive_division_report(
            lattice,
            2,
            a,
            b,
            invariant_truncation,
            EllipticFunctionTruncation::new(14)?,
            tolerance,
        )?;
        println!("  secondary division-polynomial diagnostic at r_fun = 14:");
        println!(
            "{}",
            indent(
                &describe_analytic_division_polynomial_comparison(&report),
                4
            )
        );
        println!();
    }

    Ok(())
}

fn print_n3_experiment(
    lattice: &ComplexLattice,
    invariant_truncation: LatticeSumTruncation,
    tolerance: ApproxTolerance,
) -> Result<(), Box<dyn std::error::Error>> {
    println!("n = 3");
    println!("-----");
    println!(
        "Here the main experiment is convergence: ℘(z), ℘′(z), and the cubic residual stabilize as r_fun grows."
    );
    println!();

    for &(a, b) in &[(0_usize, 1_usize), (1, 1), (1, 2)] {
        let z = primitive_point_z(lattice, 3, a, b)?;
        let mapped = primitive_mapped_point(
            lattice,
            3,
            a,
            b,
            invariant_truncation,
            EllipticFunctionTruncation::new(14)?,
            tolerance,
        )?;
        println!("{}", describe_analytic_torsion_point_approx(&mapped));

        let mut p_values = Vec::new();
        let mut p_prime_values = Vec::new();
        let mut residuals = Vec::new();

        for radius in [6_usize, 10, 14] {
            let truncation = EllipticFunctionTruncation::new(radius)?;
            let p = weierstrass_p(lattice, z, truncation)?;
            let p_prime = weierstrass_p_derivative(lattice, z, truncation)?;
            let mapped_point =
                map_torus_point_to_curve(lattice, z, invariant_truncation, truncation, tolerance)?;

            p_values.push(*p.value());
            p_prime_values.push(*p_prime.value());
            residuals.push(mapped_point.membership_report().absolute_error());
        }

        println!(
            "  Δ℘: 6→10 = {:.6e}, 10→14 = {:.6e}",
            (p_values[1] - p_values[0]).norm(),
            (p_values[2] - p_values[1]).norm()
        );
        println!(
            "  Δ℘′: 6→10 = {:.6e}, 10→14 = {:.6e}",
            (p_prime_values[1] - p_prime_values[0]).norm(),
            (p_prime_values[2] - p_prime_values[1]).norm()
        );
        println!(
            "  residuals: r=6 -> {:.6e}, r=10 -> {:.6e}, r=14 -> {:.6e}",
            residuals[0], residuals[1], residuals[2]
        );

        let report = primitive_division_report(
            lattice,
            3,
            a,
            b,
            invariant_truncation,
            EllipticFunctionTruncation::new(14)?,
            tolerance,
        )?;
        println!("  secondary division-polynomial diagnostic at r_fun = 14:");
        println!(
            "{}",
            indent(
                &describe_analytic_division_polynomial_comparison(&report),
                4
            )
        );
        println!();
    }

    Ok(())
}

fn print_n6_experiment(
    lattice: &ComplexLattice,
    tolerance: ApproxTolerance,
) -> Result<(), Box<dyn std::error::Error>> {
    println!("n = 6");
    println!("-----");
    println!(
        "For order 6, the useful test is structural: doubling should land near 3-torsion and tripling should land near 2-torsion."
    );
    println!();

    for (r_inv, r_fun) in [(16_usize, 14_usize), (24, 28)] {
        let invariant_truncation = LatticeSumTruncation::new(r_inv)?;
        let function_truncation = EllipticFunctionTruncation::new(r_fun)?;
        let six = map_primitive_torus_torsion_to_curve(
            lattice,
            6,
            invariant_truncation,
            function_truncation,
            tolerance,
        )?;
        let three = map_primitive_torus_torsion_to_curve(
            lattice,
            3,
            invariant_truncation,
            function_truncation,
            tolerance,
        )?;
        let two = map_primitive_torus_torsion_to_curve(
            lattice,
            2,
            invariant_truncation,
            function_truncation,
            tolerance,
        )?;

        let mut min_double_to_three = f64::INFINITY;
        let mut max_double_to_three = 0.0_f64;
        let mut min_triple_to_two = f64::INFINITY;
        let mut max_triple_to_two = 0.0_f64;

        for point in &six {
            let z = *point.torus_point().z();
            let double = map_torus_point_to_curve(
                lattice,
                z * 2.0,
                invariant_truncation,
                function_truncation,
                tolerance,
            )?;
            let triple = map_torus_point_to_curve(
                lattice,
                z * 3.0,
                invariant_truncation,
                function_truncation,
                tolerance,
            )?;

            let nearest_three = three
                .iter()
                .filter_map(|target| point_distance(double.point(), target.curve_point()))
                .fold(f64::INFINITY, f64::min);
            let nearest_two = two
                .iter()
                .filter_map(|target| point_distance(triple.point(), target.curve_point()))
                .fold(f64::INFINITY, f64::min);

            min_double_to_three = min_double_to_three.min(nearest_three);
            max_double_to_three = max_double_to_three.max(nearest_three);
            min_triple_to_two = min_triple_to_two.min(nearest_two);
            max_triple_to_two = max_triple_to_two.max(nearest_two);
        }

        println!(
            "r_inv = {r_inv}, r_fun = {r_fun}: [2]P to 3-torsion distance in [{min_double_to_three:.6e}, {max_double_to_three:.6e}], [3]P to 2-torsion distance in [{min_triple_to_two:.6e}, {max_triple_to_two:.6e}]"
        );
    }

    println!();

    let representative = primitive_mapped_point(
        lattice,
        6,
        1,
        1,
        LatticeSumTruncation::new(24)?,
        EllipticFunctionTruncation::new(28)?,
        tolerance,
    )?;
    println!("Representative primitive 6-torsion image at the larger truncation:");
    println!(
        "{}",
        describe_analytic_torsion_point_approx(&representative)
    );

    let report = primitive_division_report(
        lattice,
        6,
        1,
        1,
        LatticeSumTruncation::new(24)?,
        EllipticFunctionTruncation::new(28)?,
        tolerance,
    )?;
    println!("Secondary division-polynomial diagnostic:");
    println!(
        "{}",
        indent(
            &describe_analytic_division_polynomial_comparison(&report),
            2
        )
    );
    println!();

    Ok(())
}

fn print_n7_experiment(
    lattice: &ComplexLattice,
    invariant_truncation: LatticeSumTruncation,
    tolerance: ApproxTolerance,
) -> Result<(), Box<dyn std::error::Error>> {
    println!("n = 7");
    println!("-----");
    println!("For order 7, the best signal is convergence of ℘ and ℘′ as r_fun grows.");
    println!();

    for &(a, b) in &[(0_usize, 1_usize), (0, 2), (0, 3), (0, 4)] {
        let z = primitive_point_z(lattice, 7, a, b)?;
        println!("point ({a}, {b}; 7), z = {}", format_complex(&z));

        let mut p_values = Vec::new();
        let mut p_prime_values = Vec::new();

        for radius in [10_usize, 14, 20, 28] {
            let truncation = EllipticFunctionTruncation::new(radius)?;
            p_values.push(*weierstrass_p(lattice, z, truncation)?.value());
            p_prime_values.push(*weierstrass_p_derivative(lattice, z, truncation)?.value());
        }

        println!(
            "  Δ℘: 10→14 = {:.6e}, 14→20 = {:.6e}, 20→28 = {:.6e}",
            (p_values[1] - p_values[0]).norm(),
            (p_values[2] - p_values[1]).norm(),
            (p_values[3] - p_values[2]).norm()
        );
        println!(
            "  Δ℘′: 10→14 = {:.6e}, 14→20 = {:.6e}, 20→28 = {:.6e}",
            (p_prime_values[1] - p_prime_values[0]).norm(),
            (p_prime_values[2] - p_prime_values[1]).norm(),
            (p_prime_values[3] - p_prime_values[2]).norm()
        );
    }

    println!();
    let report = primitive_division_report(
        lattice,
        7,
        0,
        1,
        invariant_truncation,
        EllipticFunctionTruncation::new(28)?,
        tolerance,
    )?;
    println!("Secondary division-polynomial diagnostic for (0, 1; 7) at r_fun = 28:");
    println!(
        "{}",
        indent(
            &describe_analytic_division_polynomial_comparison(&report),
            2
        )
    );
    println!();

    Ok(())
}

fn indent(text: &str, spaces: usize) -> String {
    let padding = " ".repeat(spaces);
    text.lines()
        .map(|line| format!("{padding}{line}"))
        .collect::<Vec<_>>()
        .join("\n")
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let tau = UpperHalfPlanePoint::tau_i();
    let lattice = ComplexLattice::from_tau(tau.clone());
    let invariant_truncation = LatticeSumTruncation::new(16)?;
    let tolerance = ApproxTolerance::new(1.0e-4, 1.0e-2);
    let invariants = analytic_invariants(&lattice, invariant_truncation)?;
    let analytic_curve = AnalyticWeierstrassCurve::from_lattice(&lattice, invariant_truncation)?;
    let short_curve = analytic_curve.as_short_weierstrass();

    println!("Milestone 8 + 7 bridge: complex torus torsion experiments");
    println!("=========================================================");
    println!();
    println!("τ = {}", format_complex(tau.tau()));
    println!("{}", describe_complex_lattice(&lattice));
    println!();
    println!("{}", describe_analytic_invariants(&invariants));
    println!(
        "analytic cubic model: {}",
        format_analytic_cubic_model(&analytic_curve)
    );
    println!(
        "short-Weierstrass model: {}",
        format_short_weierstrass_over_complex(&short_curve)
    );
    println!(
        "base tolerance: abs = {:.3e}, rel = {:.3e}",
        tolerance.absolute, tolerance.relative
    );
    println!();

    print_n2_experiment(&lattice, invariant_truncation, tolerance)?;
    print_n3_experiment(&lattice, invariant_truncation, tolerance)?;
    print_n6_experiment(&lattice, tolerance)?;
    print_n7_experiment(&lattice, invariant_truncation, tolerance)?;

    Ok(())
}
