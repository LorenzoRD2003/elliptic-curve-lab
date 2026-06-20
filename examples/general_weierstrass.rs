use elliptic_algorithms_lab::elliptic_curves::{
    GeneralWeierstrassCurve,
    traits::{CurveModelConversion, EnumerableCurveModel, GroupCurveModel},
};
use elliptic_algorithms_lab::fields::{Fp, traits::Field};
use elliptic_algorithms_lab::visualization::{
    describe_general_weierstrass_curve, describe_general_weierstrass_short_reduction,
    format_point_compact,
};

type F = Fp<5>;

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
    println!("{}", describe_general_weierstrass_curve(&curve));
    println!();
    println!("{}", describe_general_weierstrass_short_reduction(&curve));
    println!();
    println!("Sample calculation");
    println!("------------------");
    println!("P  = {}", format_point_compact(&left));
    println!("Q  = {}", format_point_compact(&right));
    println!("P' = {}", format_point_compact(&short_left));
    println!("Q' = {}", format_point_compact(&short_right));
    println!();
    println!(
        "P + Q on the general model       = {}",
        format_point_compact(&general_sum)
    );
    println!(
        "P' + Q' on the short companion   = {}",
        format_point_compact(&short_sum)
    );
    println!(
        "transport back(P' + Q')          = {}",
        format_point_compact(&transported_back)
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
