use elliptic_algorithms_lab::elliptic_curves::{
    GeneralWeierstrassCurve,
    traits::{CurveModelConversion, EnumerableCurveModel, GroupCurveModel},
};
use elliptic_algorithms_lab::fields::traits::*;
use elliptic_algorithms_lab::visualization::Visualizable;

type F = elliptic_algorithms_lab::fields::Fp5;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let curve =
        GeneralWeierstrassCurve::<F>::new(F::one(), F::one(), F::one(), F::one(), F::zero())?;
    let finite_points = curve.finite_points();
    let left = finite_points
        .first()
        .cloned()
        .expect("the chosen example curve should have a finite point");
    let right = finite_points
        .get(1)
        .cloned()
        .unwrap_or_else(|| left.clone());
    let general_sum = curve.add(&left, &right)?;
    let conversion = curve.conversion_to_short_weierstrass()?;
    let short_left = conversion.map_source_point(&left)?;
    let short_right = conversion.map_source_point(&right)?;
    let short_sum = conversion.target().add(&short_left, &short_right)?;
    let transported_back = conversion.map_target_point(&short_sum)?;

    println!("General Weierstrass educational walkthrough");
    println!("======================================================");
    println!();
    println!("{}", curve.describe());
    println!();
    println!("Short-Weierstrass companion");
    println!("---------------------------");
    println!("{}", conversion.target().describe());
    println!();
    println!("Sample calculation");
    println!("------------------");
    println!("P  = {}", left.format_compact());
    println!("Q  = {}", right.format_compact());
    println!("P' = {}", short_left.format_compact());
    println!("Q' = {}", short_right.format_compact());
    println!();
    println!(
        "P + Q on the general model       = {}",
        general_sum.format_compact()
    );
    println!(
        "P' + Q' on the short companion   = {}",
        short_sum.format_compact()
    );
    println!(
        "transport back(P' + Q')          = {}",
        transported_back.format_compact()
    );
    println!(
        "agreement                        = {}",
        if general_sum == transported_back {
            "yes"
        } else {
            "no"
        }
    );

    Ok(())
}
