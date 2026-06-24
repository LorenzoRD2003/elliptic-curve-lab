use elliptic_algorithms_lab::elliptic_curves::{
    TwistedEdwardsCurve,
    traits::{AffineCurveModel, CurveModel, GroupCurveModel},
};
use elliptic_algorithms_lab::fields::{Fp, traits::Field};
use elliptic_algorithms_lab::visualization::elliptic_curves::{
    describe_montgomery_curve, describe_twisted_edwards_birational_transport,
    describe_twisted_edwards_curve, describe_twisted_edwards_montgomery_companion,
};
use elliptic_algorithms_lab::visualization::format_point_compact;

type F5 = Fp<5>;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let curve = TwistedEdwardsCurve::<F5>::new(F5::one(), F5::from_i64(2))?;
    let montgomery = curve.as_montgomery();
    let left = curve.point(F5::one(), F5::zero())?;
    let right = curve.point(F5::from_i64(2), F5::from_i64(2))?;
    let sum = curve.add(&left, &right)?;
    let identity = curve.identity();
    let second_x_zero = curve.point(F5::zero(), F5::from_i64(-1))?;

    let left_m = curve.try_point_to_montgomery_open(&left)?;
    let right_m = curve.try_point_to_montgomery_open(&right)?;
    let sum_m = montgomery.add(&left_m, &right_m)?;
    let sum_roundtrip = montgomery.try_point_to_twisted_edwards_open(&sum_m)?;
    let identity_m = curve.point_to_montgomery(&identity)?;
    let second_x_zero_m = curve.point_to_montgomery(&second_x_zero)?;

    println!("Twisted Edwards educational walkthrough");
    println!("======================================================");
    println!();
    println!("{}", describe_twisted_edwards_curve(&curve));
    println!();
    println!("{}", describe_twisted_edwards_montgomery_companion(&curve));
    println!();
    println!("{}", describe_twisted_edwards_birational_transport(&curve));
    println!();
    println!("Montgomery companion view");
    println!("------------------------");
    println!("{}", describe_montgomery_curve(&montgomery));
    println!();
    println!("Transported sample calculation");
    println!("-----------------------------");
    println!("P                = {}", format_point_compact(&left));
    println!("Q                = {}", format_point_compact(&right));
    println!("P + Q (native)   = {}", format_point_compact(&sum));
    println!();
    println!("phi(P)           = {}", format_point_compact(&left_m));
    println!("phi(Q)           = {}", format_point_compact(&right_m));
    println!("phi(P) + phi(Q)  = {}", format_point_compact(&sum_m));
    println!(
        "phi^(-1)(sum)    = {}",
        format_point_compact(&sum_roundtrip)
    );
    println!(
        "agreement        = {}",
        if sum == sum_roundtrip { "yes" } else { "no" }
    );
    println!();
    println!("Extended exceptional-point transport");
    println!("-----------------------------------");
    println!("(0, 1)           = {}", format_point_compact(&identity));
    println!("phi(0, 1)        = {}", format_point_compact(&identity_m));
    println!(
        "(0, -1)          = {}",
        format_point_compact(&second_x_zero)
    );
    println!(
        "phi(0, -1)       = {}",
        format_point_compact(&second_x_zero_m)
    );
    println!();
    println!(
        "note: Edwards -> Montgomery is now total on affine Edwards points, but the reverse affine transport is still only partial."
    );
    Ok(())
}
