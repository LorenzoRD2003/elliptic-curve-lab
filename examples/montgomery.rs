use elliptic_algorithms_lab::elliptic_curves::{
    AffinePoint, MontgomeryCurve,
    montgomery::MontgomeryXzPoint,
    traits::{CurveModelConversion, GroupCurveModel, LiftXCoordinate},
};
use elliptic_algorithms_lab::visualization::Visualizable;

type F89 = elliptic_algorithms_lab::fields::Fp89;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let curve = MontgomeryCurve::<F89>::new(F89::from_i64(3), F89::from_i64(2))?;
    let base_x = F89::from_i64(2);
    let point = curve
        .point_from_x(base_x.clone())?
        .expect("the chosen x-coordinate should lift to a Montgomery point");
    let scalar = 37u64;
    let base_x = match &point {
        AffinePoint::Finite { x, .. } => x.clone(),
        AffinePoint::Infinity => unreachable!("lifting a finite x should not return infinity"),
    };
    let ladder_report = curve.try_ladder_x_report(base_x.clone(), scalar)?;
    let affine_multiple = curve.mul_scalar(&point, scalar)?;
    let affine_multiple_x = MontgomeryXzPoint::from_affine_point(&affine_multiple);
    let general = curve.as_general_weierstrass();

    println!("Montgomery ladder educational walkthrough");
    println!("======================================================");
    println!();
    println!("{}", curve.describe());
    println!();
    println!("{}", ladder_report.describe());
    println!();
    println!("Short-Weierstrass companion");
    println!("---------------------------");
    match curve.conversion_to_short_weierstrass() {
        Ok(conversion) => println!("{}", conversion.target().describe()),
        Err(error) => println!("unavailable: {error}"),
    }
    println!();
    println!("General-Weierstrass embedding");
    println!("-----------------------------");
    println!("{}", general.describe());
    println!();
    println!("General companion view");
    println!("---------------------");
    println!("{}", general.describe());
    println!();
    println!("x-only sample calculation");
    println!("-------------------------");
    println!("P              = {}", point.format_compact());
    println!("x(P)           = {}", base_x);
    println!("n              = {}", scalar);
    println!();
    println!(
        "ladder x([n]P) = {}",
        ladder_report.multiple_x().format_compact()
    );
    println!(
        "ladder x([n+1]P) = {}",
        ladder_report.next_multiple_x().format_compact()
    );
    println!("affine [n]P    = {}", affine_multiple.format_compact());
    println!(
        "x([n]P) from affine validation = {}",
        affine_multiple_x.format_compact()
    );
    println!(
        "agreement     = {}",
        if ladder_report.multiple_x() == &affine_multiple_x {
            "yes"
        } else {
            "no"
        }
    );
    println!();
    println!(
        "note: the ladder result is x-only; recovering the sign of y([n]P) would require extra data."
    );

    Ok(())
}
