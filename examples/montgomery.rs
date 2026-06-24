use elliptic_algorithms_lab::elliptic_curves::{
    AffinePoint, MontgomeryCurve,
    montgomery::MontgomeryXzPoint,
    traits::{GroupCurveModel, LiftXCoordinate},
};
use elliptic_algorithms_lab::fields::{Fp, traits::Field};
use elliptic_algorithms_lab::visualization::elliptic_curves::{
    describe_montgomery_curve, describe_montgomery_general_embedding,
    describe_montgomery_ladder_report, describe_montgomery_short_reduction,
    format_montgomery_xz_point,
};
use elliptic_algorithms_lab::visualization::{
    describe_general_weierstrass_curve, format_point_compact,
};

type F89 = Fp<89>;

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
    println!("{}", describe_montgomery_curve(&curve));
    println!();
    println!(
        "{}",
        describe_montgomery_ladder_report(&curve, &ladder_report)
    );
    println!();
    println!("{}", describe_montgomery_short_reduction(&curve));
    println!();
    println!("{}", describe_montgomery_general_embedding(&curve));
    println!();
    println!("General companion view");
    println!("---------------------");
    println!("{}", describe_general_weierstrass_curve(&general));
    println!();
    println!("x-only sample calculation");
    println!("-------------------------");
    println!("P              = {}", format_point_compact(&point));
    println!("x(P)           = {}", base_x);
    println!("n              = {}", scalar);
    println!();
    println!(
        "ladder x([n]P) = {}",
        format_montgomery_xz_point(ladder_report.multiple_x())
    );
    println!(
        "ladder x([n+1]P) = {}",
        format_montgomery_xz_point(ladder_report.next_multiple_x())
    );
    println!(
        "affine [n]P    = {}",
        format_point_compact(&affine_multiple)
    );
    println!(
        "x([n]P) from affine validation = {}",
        format_montgomery_xz_point(&affine_multiple_x)
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
