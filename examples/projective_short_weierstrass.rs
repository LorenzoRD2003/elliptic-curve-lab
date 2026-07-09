use elliptic_algorithms_lab::elliptic_curves::{
    AffinePoint, ProjectivePoint, ShortWeierstrassCurve,
    short_weierstrass::projective::{
        ShortWeierstrassProjectiveOperationCost, ShortWeierstrassProjectiveOperationKind,
    },
    traits::{AffineCurveModel, HasProjectiveModel, ProjectiveGroupCurveModel},
};
use elliptic_algorithms_lab::fields::traits::*;
use elliptic_algorithms_lab::visualization::Visualizable;

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
    println!("curve: {}", curve.format_compact());
    println!("P (affine): {}", left.format_compact());
    println!("Q (affine): {}", right.format_compact());
    println!(
        "P (projective representative): {}",
        left_projective.format_compact()
    );
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
        "mixed P + Q (projective + affine): {}",
        mixed_sum.format_compact()
    );
    println!();
    println!(
        "{}",
        ShortWeierstrassProjectiveOperationCost::for_kind(
            ShortWeierstrassProjectiveOperationKind::Normalize,
        )
        .describe()
    );
    println!();
    println!(
        "{}",
        ShortWeierstrassProjectiveOperationCost::for_kind(
            ShortWeierstrassProjectiveOperationKind::Add
        )
        .describe()
    );
    println!();
    println!(
        "{}",
        ShortWeierstrassProjectiveOperationCost::for_kind(
            ShortWeierstrassProjectiveOperationKind::MixedAdd,
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
