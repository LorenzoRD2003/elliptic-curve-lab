use elliptic_algorithms_lab::{
    AffineCurveModel, CurveError, EnumerableCurveModel, Field, FiniteGroupCurveModel, Fp,
    ShortWeierstrassCurve, format_point_compact, summarize_group_structure,
    summarize_order_distribution,
};

// Change prime to test other finite fields: Fp<101>
type F = Fp<101>;

fn main() -> Result<(), CurveError> {
    // E(F): y^2 = x^3 + 2x + 3
    let curve = ShortWeierstrassCurve::<F>::new(F::from_i64(2), F::from_i64(3))?;
    let point = curve.point(F::from_i64(1), F::from_i64(39))?;
    let point_order = curve
        .point_order(&point)
        .expect("chosen sample point should lie on the curve");

    println!("Second milestone: finite elliptic curve group structure");
    println!("======================================================");
    println!();
    println!("curve: y^2 = x^3 + 2x + 3 over F_101");
    println!("#E(F_101) = {}", curve.order());
    println!();
    println!("identity:");
    println!("  O");
    println!();
    println!("sample point:");
    println!("  P = {}", format_point_compact(&point));
    println!("  order(P) = {point_order}");
    println!();
    println!("order distribution:");
    for line in summarize_order_distribution(&curve).lines() {
        println!("  {line}");
    }
    println!();
    println!("group summary:");
    for line in summarize_group_structure(&curve).lines() {
        println!("  {line}");
    }

    Ok(())
}
