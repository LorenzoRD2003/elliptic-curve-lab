use std::collections::BTreeMap;

use elliptic_algorithms_lab::elliptic_curves::{
    CurveError, ShortWeierstrassCurve,
    traits::{AffineCurveModel, EnumerableCurveModel, FiniteGroupCurveModel},
};
use elliptic_algorithms_lab::visualization::Visualizable;

// Change prime to test other finite fields: elliptic_algorithms_lab::fields::Fp101
type F = elliptic_algorithms_lab::fields::Fp101;

fn main() -> Result<(), CurveError> {
    // E(F): y^2 = x^3 + 2x + 3
    let curve = ShortWeierstrassCurve::<F>::new(F::from_i64(2), F::from_i64(3))?;
    let point = curve.point(F::from_i64(1), F::from_i64(39))?;
    let point_order = curve
        .point_order(&point)
        .expect("chosen sample point should lie on the curve");

    println!("Finite elliptic curve group structure");
    println!("======================================================");
    println!();
    println!("curve: y^2 = x^3 + 2x + 3 over F_101");
    println!("#E(F_101) = {}", curve.order());
    println!();
    println!("identity:");
    println!("  O");
    println!();
    println!("sample point:");
    println!("  P = {}", point.format_compact());
    println!("  order(P) = {point_order}");
    println!();
    println!("order distribution:");
    for (order, count) in order_distribution(&curve) {
        println!("  {order} -> {count}");
    }
    println!();
    println!("group summary:");
    for line in curve.group_structure().describe().lines() {
        println!("  {line}");
    }

    Ok(())
}

fn order_distribution(curve: &ShortWeierstrassCurve<F>) -> BTreeMap<usize, usize> {
    let mut distribution = BTreeMap::new();
    for point in curve.points() {
        let order = curve
            .point_order(&point)
            .expect("enumerated points have finite order");
        *distribution.entry(order).or_insert(0) += 1;
    }
    distribution
}
