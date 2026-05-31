use elliptic_algorithms_lab::{
    CurveError, EnumerableCurveModel, Field, Fp, ShortWeierstrassCurve, format_curve, format_point,
};

type F = Fp<101>;

fn main() -> Result<(), CurveError> {
    let curve = ShortWeierstrassCurve::<F>::new(F::from_i64(2), F::from_i64(3))?;
    let points = curve.points();
    let order = curve.order();

    println!("First milestone: curve order over a small prime field");
    println!("====================================================");
    println!();
    println!("curve: {}", format_curve(&curve));
    println!("#E(F_101) = {order}");
    println!("sample of enumerated points:");

    for point in points.iter().take(10) {
        println!("  {}", format_point(point));
    }

    if points.len() > 10 {
        println!("  ...");
    }

    Ok(())
}
