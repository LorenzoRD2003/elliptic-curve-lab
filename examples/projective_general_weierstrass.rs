use elliptic_algorithms_lab::elliptic_curves::{
    AffinePoint, GeneralWeierstrassCurve, ProjectivePoint,
    general_weierstrass::projective::{
        GeneralWeierstrassProjectiveOperationCost, GeneralWeierstrassProjectiveOperationKind,
    },
    traits::{
        AffineCurveModel, CurveModelConversion, GroupCurveModel, HasProjectiveModel,
        ProjectiveGroupCurveModel,
    },
};
use elliptic_algorithms_lab::fields::{Fp, traits::Field};
use elliptic_algorithms_lab::visualization::{
    describe_general_weierstrass_projective_cost, describe_general_weierstrass_short_reduction,
    describe_projective_affine_roundtrip, describe_projective_normalization,
    format_general_weierstrass_curve, format_point_compact, format_projective_point,
};

type F = Fp<5>;

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
    let curve =
        GeneralWeierstrassCurve::<F>::new(F::one(), F::one(), F::one(), F::one(), F::zero())?;
    let left = curve.point(F::zero(), F::zero())?;
    let right = curve.point(F::from_i64(2), F::one())?;
    let left_projective = scaled_projective(&left, 2);
    let sum = curve.add_projective(&left_projective, &curve.to_projective(&right)?)?;
    let conversion = curve.conversion_to_short_weierstrass()?;
    let short_sum = conversion.target().add(
        &conversion.map_source_point(&left)?,
        &conversion.map_source_point(&right)?,
    )?;

    println!("Projective general-Weierstrass walkthrough");
    println!("=========================================");
    println!("curve: {}", format_general_weierstrass_curve(&curve));
    println!("P (affine): {}", format_point_compact(&left));
    println!("Q (affine): {}", format_point_compact(&right));
    println!(
        "P (projective representative): {}",
        format_projective_point(&left_projective)
    );
    println!();
    println!("{}", describe_general_weierstrass_short_reduction(&curve));
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
        "P + Q via short companion: {}",
        format_point_compact(&conversion.map_target_point(&short_sum)?)
    );
    println!();
    println!(
        "{}",
        describe_general_weierstrass_projective_cost(
            &GeneralWeierstrassProjectiveOperationCost::for_kind(
                GeneralWeierstrassProjectiveOperationKind::Normalize,
            ),
        )
    );
    println!();
    println!(
        "{}",
        describe_general_weierstrass_projective_cost(
            &GeneralWeierstrassProjectiveOperationCost::for_kind(
                GeneralWeierstrassProjectiveOperationKind::Add,
            ),
        )
    );

    Ok(())
}
