use elliptic_algorithms_lab::elliptic_curves::division_polynomials::{
    DivisionPolynomialForm, compare_division_polynomial_torsion_with_enumeration,
    division_polynomial, exact_n_torsion_points_from_division_polynomial,
    rational_x_candidates_for_division_polynomial, torsion_candidates_from_division_polynomial,
    torsion_points_from_division_polynomial,
};
use elliptic_algorithms_lab::{
    AffinePoint, Field, Fp, GroupCurveModel, ShortWeierstrassCurve, format_curve,
    format_dense_polynomial, format_point_compact,
};

type F = Fp<11>;

fn format_points(points: &[AffinePoint<F>]) -> String {
    if points.is_empty() {
        "[]".to_string()
    } else {
        points
            .iter()
            .map(format_point_compact::<F>)
            .collect::<Vec<_>>()
            .join(", ")
    }
}

fn format_xs(xs: &[<F as Field>::Elem]) -> String {
    if xs.is_empty() {
        "[]".to_string()
    } else {
        xs.iter()
            .map(ToString::to_string)
            .collect::<Vec<_>>()
            .join(", ")
    }
}

fn format_division_polynomial(form: &DivisionPolynomialForm<F>) -> String {
    match form {
        DivisionPolynomialForm::InX(polynomial) => format_dense_polynomial(polynomial),
        DivisionPolynomialForm::YTimes(polynomial) => {
            format!("y * ({})", format_dense_polynomial(polynomial))
        }
    }
}

fn show_division_polynomial_walkthrough(
    curve: &ShortWeierstrassCurve<F>,
    n: usize,
) -> Result<(), Box<dyn std::error::Error>> {
    let polynomial = division_polynomial(curve, n)?;
    let rational_roots = rational_x_candidates_for_division_polynomial(curve, n)?;
    let lifted_points = torsion_candidates_from_division_polynomial(curve, n)?;
    let n_torsion_points = torsion_points_from_division_polynomial(curve, n)?;
    let exact_order_points = exact_n_torsion_points_from_division_polynomial(curve, n)?;
    let comparison = compare_division_polynomial_torsion_with_enumeration(curve, n)?;

    println!("n = {n}");
    println!("-----");
    println!("ψ_{n} = {}", format_division_polynomial(&polynomial));
    println!("raíces racionales de ψ_{n}: {}", format_xs(&rational_roots));
    println!(
        "puntos racionales levantados desde esas raíces: {}",
        format_points(&lifted_points)
    );
    println!("verificación de [{n}]P = O sobre los puntos levantados:");
    if lifted_points.is_empty() {
        println!("  no hay puntos racionales candidatos para verificar");
    } else {
        for point in &lifted_points {
            let is_n_torsion = curve.is_torsion_point(point, n as u64);
            println!("  {} -> {}", format_point_compact(point), is_n_torsion);
        }
    }
    println!(
        "puntos que satisfacen [{n}]P = O: {}",
        format_points(&n_torsion_points)
    );
    println!(
        "puntos de orden exacto {n}: {}",
        format_points(&exact_order_points)
    );
    println!("comparación contra enumeración:");
    println!(
        "  candidatos por polinomio: {}",
        comparison.polynomial_candidate_count
    );
    println!(
        "  puntos {n}-torsión por polinomio: {}",
        comparison.polynomial_n_torsion_count
    );
    println!(
        "  puntos {n}-torsión por enumeración: {}",
        comparison.enumerated_n_torsion_count
    );
    println!(
        "  puntos de orden exacto {n} por polinomio: {}",
        comparison.exact_order_polynomial_count
    );
    println!(
        "  puntos de orden exacto {n} por enumeración: {}",
        comparison.exact_order_enumerated_count
    );
    println!(
        "  faltantes desde el polinomio: {}",
        format_points(&comparison.missing_from_polynomial)
    );
    println!(
        "  extras desde el polinomio: {}",
        format_points(&comparison.extra_from_polynomial)
    );
    println!();

    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let curve = ShortWeierstrassCurve::<F>::new(F::from_i64(1), F::from_i64(7))?;

    println!("Seventh milestone: division polynomials and rational torsion");
    println!("============================================================");
    println!();
    println!("curva pequeña E/Fp:");
    println!("  {}", format_curve(&curve));
    println!();

    show_division_polynomial_walkthrough(&curve, 3)?;
    show_division_polynomial_walkthrough(&curve, 5)?;

    Ok(())
}
