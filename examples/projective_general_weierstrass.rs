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
use elliptic_algorithms_lab::fields::traits::*;
use elliptic_algorithms_lab::visualization::Visualizable;

type F = elliptic_algorithms_lab::fields::Fp5;

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
    println!("curve: {}", curve.format_compact());
    println!("P (affine): {}", left.format_compact());
    println!("Q (affine): {}", right.format_compact());
    println!(
        "P (projective representative): {}",
        left_projective.format_compact()
    );
    println!();
    println!("Short-Weierstrass companion");
    println!("---------------------------");
    println!("{}", conversion.target().describe());
    println!();
    println!("{}", describe_projective_normalization(&left_projective)?);
    println!();
    println!(
        "{}",
        describe_projective_affine_roundtrip(&left_projective)?
    );
    println!();
    println!("P + Q (projective baseline): {}", sum.format_compact());
    println!(
        "P + Q via short companion: {}",
        conversion.map_target_point(&short_sum)?.format_compact()
    );
    println!();
    println!(
        "{}",
        GeneralWeierstrassProjectiveOperationCost::for_kind(
            GeneralWeierstrassProjectiveOperationKind::Normalize,
        )
        .describe()
    );
    println!();
    println!(
        "{}",
        GeneralWeierstrassProjectiveOperationCost::for_kind(
            GeneralWeierstrassProjectiveOperationKind::Add,
        )
        .describe()
    );

    Ok(())
}

fn describe_projective_normalization(
    point: &ProjectivePoint<F>,
) -> Result<String, Box<dyn std::error::Error>> {
    Ok([
        "Projective normalization".to_string(),
        format!("input: {}", point.format_compact()),
        format!(
            "already normalized: {}",
            if point.is_normalized() { "yes" } else { "no" }
        ),
        format!(
            "normalized representative: {}",
            point.normalize()?.format_compact()
        ),
        "rule: finite points are rescaled to the Z = 1 chart".to_string(),
    ]
    .join("\n"))
}

fn describe_projective_affine_roundtrip(
    point: &ProjectivePoint<F>,
) -> Result<String, Box<dyn std::error::Error>> {
    let affine = point.to_affine()?;
    let lifted_back = ProjectivePoint::from_affine(&affine);

    Ok([
        "Affine/projective roundtrip".to_string(),
        format!("projective input: {}", point.format_compact()),
        format!("affine chart point: {}", affine.format_compact()),
        format!("lifted projective point: {}", lifted_back.format_compact()),
    ]
    .join("\n"))
}
