use elliptic_algorithms_lab::elliptic_curves::{
    TwistedEdwardsCurve,
    traits::{AffineCurveModel, CurveModel, GroupCurveModel},
};
use elliptic_algorithms_lab::fields::traits::*;
use elliptic_algorithms_lab::visualization::Visualizable;

type F5 = elliptic_algorithms_lab::fields::Fp5;

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
    println!("{}", curve.describe());
    println!();
    println!("Montgomery companion");
    println!("--------------------");
    println!("{}", montgomery.describe());
    println!();
    println!("Montgomery companion view");
    println!("------------------------");
    println!("{}", montgomery.describe());
    println!();
    println!("Transported sample calculation");
    println!("-----------------------------");
    println!("P                = {}", left.format_compact());
    println!("Q                = {}", right.format_compact());
    println!("P + Q (native)   = {}", sum.format_compact());
    println!();
    println!("phi(P)           = {}", left_m.format_compact());
    println!("phi(Q)           = {}", right_m.format_compact());
    println!("phi(P) + phi(Q)  = {}", sum_m.format_compact());
    println!("phi^(-1)(sum)    = {}", sum_roundtrip.format_compact());
    println!(
        "agreement        = {}",
        if sum == sum_roundtrip { "yes" } else { "no" }
    );
    println!();
    println!("Extended exceptional-point transport");
    println!("-----------------------------------");
    println!("(0, 1)           = {}", identity.format_compact());
    println!("phi(0, 1)        = {}", identity_m.format_compact());
    println!("(0, -1)          = {}", second_x_zero.format_compact());
    println!("phi(0, -1)       = {}", second_x_zero_m.format_compact());
    println!();
    println!(
        "note: Edwards -> Montgomery is now total on affine Edwards points, but the reverse affine transport is still only partial."
    );
    Ok(())
}
