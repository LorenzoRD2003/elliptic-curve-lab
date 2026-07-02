use elliptic_algorithms_lab::elliptic_curves::{
    AffinePoint, ProjectivePoint, ShortWeierstrassCurve,
    short_weierstrass::projective::{
        ShortWeierstrassProjectiveOperationCost, ShortWeierstrassProjectiveOperationKind,
    },
    traits::{AffineCurveModel, HasProjectiveModel, ProjectiveGroupCurveModel},
};
use elliptic_algorithms_lab::fields::traits::*;
use elliptic_algorithms_lab::visualization::{
    describe_projective_affine_roundtrip, describe_projective_normalization,
    describe_short_weierstrass_projective_cost, format_curve, format_point_compact,
    format_projective_point,
};

type F = elliptic_algorithms_lab::fields::Fp7;

fn scaled_projective(point: &AffinePoint<F>, scale: i64) -> ProjectivePoint<F> {
    match point {
        AffinePoint::Infinity => ProjectivePoint::Infinity,
        AffinePoint::Finite { x, y } => {
            let lambda = F::from_i64(scale);
            ProjectivePoint::new(F::mul(x, &lambda), F::mul(y, &lambda), lambda)
        }
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let curve = ShortWeierstrassCurve::<F>::new(F::from_i64(2), F::from_i64(3))?;
    let left = curve.point(F::from_i64(2), F::from_i64(1))?;
    let right = curve.point(F::from_i64(3), F::from_i64(1))?;
    let left_projective = scaled_projective(&left, 3);
    let sum = curve.add_projective(&left_projective, &curve.to_projective(&right)?)?;
    let mixed_sum = curve.mixed_add_projective(&left_projective, &right)?;

    println!("Projective short-Weierstrass walkthrough");
    println!("=======================================");
    println!("curve: {}", format_curve(&curve));
    println!("P (affine): {}", format_point_compact(&left));
    println!("Q (affine): {}", format_point_compact(&right));
    println!(
        "P (projective representative): {}",
        format_projective_point(&left_projective)
    );
    println!();
    println!("{}", describe_projective_normalization(&left_projective)?);
    println!();
    println!(
        "{}",
        describe_projective_affine_roundtrip(&left_projective)?
    );
    println!();
    println!(
        "P + Q (projective baseline): {}",
        format_projective_point(&sum)
    );
    println!(
        "mixed P + Q (projective + affine): {}",
        format_projective_point(&mixed_sum)
    );
    println!();
    println!(
        "{}",
        describe_short_weierstrass_projective_cost(
            &ShortWeierstrassProjectiveOperationCost::for_kind(
                ShortWeierstrassProjectiveOperationKind::Normalize,
            ),
        )
    );
    println!();
    println!(
        "{}",
        describe_short_weierstrass_projective_cost(
            &ShortWeierstrassProjectiveOperationCost::for_kind(
                ShortWeierstrassProjectiveOperationKind::Add,
            ),
        )
    );
    println!();
    println!(
        "{}",
        describe_short_weierstrass_projective_cost(
            &ShortWeierstrassProjectiveOperationCost::for_kind(
                ShortWeierstrassProjectiveOperationKind::MixedAdd,
            ),
        )
    );

    Ok(())
}
