use crate::visualization::*;
use core::fmt;

use crate::elliptic_curves::{
    CurveError, ProjectivePoint,
    general_weierstrass::projective::GeneralWeierstrassProjectiveOperationCost,
    short_weierstrass::projective::ShortWeierstrassProjectiveOperationCost,
};
use crate::fields::traits::Field;
use crate::visualization::{
    Visualizable, VisualizableField,
    shared::{format_field_elem as format_elem, yes_no},
};

/// Formats a projective point compactly.
fn format_projective_point<F: Field>(point: &ProjectivePoint<F>) -> String
where
    F::Elem: VisualizableField + fmt::Display,
{
    match point {
        ProjectivePoint::Infinity => "O".to_string(),
        ProjectivePoint::Finite { x, y, z } => format!(
            "({} : {} : {})",
            format_elem::<F>(x),
            format_elem::<F>(y),
            format_elem::<F>(z)
        ),
    }
}

/// Describes one projective point as stored in the current homogeneous chart.
fn describe_projective_point<F: Field>(point: &ProjectivePoint<F>) -> String
where
    F::Elem: VisualizableField + fmt::Display,
{
    match point {
        ProjectivePoint::Infinity => [
            "Projective point".to_string(),
            "point: O".to_string(),
            "role: distinguished point at infinity".to_string(),
            "normalized: yes".to_string(),
        ]
        .join("\n"),
        ProjectivePoint::Finite { x, y, z } => [
            "Projective point".to_string(),
            format!("point: {}", format_projective_point(point)),
            format!("X-coordinate: {}", format_elem::<F>(x)),
            format!("Y-coordinate: {}", format_elem::<F>(y)),
            format!("Z-coordinate: {}", format_elem::<F>(z)),
            format!("normalized: {}", yes_no(point.is_normalized())),
        ]
        .join("\n"),
    }
}

/// Explains the normalization step for one projective point.
fn describe_projective_normalization<F: Field>(
    point: &ProjectivePoint<F>,
) -> Result<String, CurveError>
where
    F::Elem: VisualizableField + fmt::Display,
{
    let normalized = point.normalize()?;
    Ok([
        "Projective normalization".to_string(),
        format!("input: {}", format_projective_point(point)),
        format!("already normalized: {}", yes_no(point.is_normalized())),
        format!(
            "normalized representative: {}",
            format_projective_point(&normalized)
        ),
        "rule: finite points are rescaled to the Z = 1 chart".to_string(),
    ]
    .join("\n"))
}

/// Explains the affine/projective roundtrip for one stored representative.
fn describe_projective_affine_roundtrip<F: Field>(
    point: &ProjectivePoint<F>,
) -> Result<String, CurveError>
where
    F::Elem: VisualizableField + fmt::Display,
{
    let affine = point.to_affine()?;
    let lifted_back = ProjectivePoint::from_affine(&affine);
    Ok([
        "Projective-affine roundtrip".to_string(),
        format!("projective input: {}", format_projective_point(point)),
        format!("affine image: {}", affine.format_compact()),
        format!(
            "lifted-back normalized representative: {}",
            format_projective_point(&lifted_back)
        ),
        "note: the roundtrip chooses the normalized projective chart Z = 1".to_string(),
    ]
    .join("\n"))
}

/// Describes one educational cost model for the current short-Weierstrass
/// projective baseline.
fn describe_short_weierstrass_projective_cost(
    cost: &ShortWeierstrassProjectiveOperationCost,
) -> String {
    [
        "Short-Weierstrass projective operation cost".to_string(),
        format!("operation: {:?}", cost.kind()),
        format!(
            "representation cost: additions={}, multiplications={}, squarings={}, inversions={}",
            cost.representation_cost().additions(),
            cost.representation_cost().multiplications(),
            cost.representation_cost().squarings(),
            cost.representation_cost().inversions()
        ),
        format!("delegated affine additions: {}", cost.affine_additions()),
        format!("delegated affine doublings: {}", cost.affine_doublings()),
        format!("note: {}", cost.note()),
    ]
    .join("\n")
}

impl Visualizable for ShortWeierstrassProjectiveOperationCost {
    fn format_compact(&self) -> String {
        format!(
            "short-Weierstrass projective {:?}: {}M, {}S, {}I",
            self.kind(),
            self.representation_cost().multiplications(),
            self.representation_cost().squarings(),
            self.representation_cost().inversions()
        )
    }

    fn describe(&self) -> String {
        describe_short_weierstrass_projective_cost(self)
    }
}

/// Describes one educational cost model for the current general-Weierstrass
/// projective baseline.
fn describe_general_weierstrass_projective_cost(
    cost: &GeneralWeierstrassProjectiveOperationCost,
) -> String {
    [
        "General-Weierstrass projective operation cost".to_string(),
        format!("operation: {:?}", cost.kind()),
        format!(
            "representation cost: additions={}, multiplications={}, squarings={}, inversions={}",
            cost.representation_cost().additions(),
            cost.representation_cost().multiplications(),
            cost.representation_cost().squarings(),
            cost.representation_cost().inversions()
        ),
        format!("delegated affine additions: {}", cost.affine_additions()),
        format!("delegated affine doublings: {}", cost.affine_doublings()),
        format!("note: {}", cost.note()),
    ]
    .join("\n")
}

impl Visualizable for GeneralWeierstrassProjectiveOperationCost {
    fn format_compact(&self) -> String {
        format!(
            "general-Weierstrass projective {:?}: {}M, {}S, {}I",
            self.kind(),
            self.representation_cost().multiplications(),
            self.representation_cost().squarings(),
            self.representation_cost().inversions()
        )
    }

    fn describe(&self) -> String {
        describe_general_weierstrass_projective_cost(self)
    }
}

impl<F: Field> Visualizable for ProjectivePoint<F>
where
    F::Elem: VisualizableField + fmt::Display,
{
    fn format_compact(&self) -> String {
        format_projective_point(self)
    }

    fn describe(&self) -> String {
        describe_projective_point(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::elliptic_curves::{
        ProjectivePoint,
        general_weierstrass::projective::{
            GeneralWeierstrassProjectiveOperationCost, GeneralWeierstrassProjectiveOperationKind,
        },
        short_weierstrass::projective::{
            ShortWeierstrassProjectiveOperationCost, ShortWeierstrassProjectiveOperationKind,
        },
    };
    use crate::fields::traits::Field;
    use crate::visualization::traits::Visualizable;

    type F7 = crate::fields::Fp7;

    #[test]
    fn compact_formatter_uses_colon_separated_coordinates() {
        let point = ProjectivePoint::<F7>::new(F7::from_i64(2), F7::from_i64(5), F7::one());

        assert_eq!(format_projective_point(&point), "(2 : 5 : 1)");
        assert_eq!(ProjectivePoint::<F7>::infinity().format_compact(), "O");
    }

    #[test]
    fn projective_description_mentions_normalization_status() {
        let point = ProjectivePoint::<F7>::new(F7::from_i64(6), F7::from_i64(1), F7::from_i64(3));
        let description = describe_projective_point(&point);

        assert!(description.contains("Projective point"));
        assert!(description.contains("normalized: no"));
    }

    #[test]
    fn normalization_and_roundtrip_descriptions_surface_the_expected_story() {
        let point = ProjectivePoint::<F7>::new(F7::from_i64(6), F7::from_i64(1), F7::from_i64(3));

        let normalization = describe_projective_normalization(&point)
            .expect("nonzero z should admit a normalization explanation");
        let roundtrip = describe_projective_affine_roundtrip(&point)
            .expect("nonzero z should admit an affine roundtrip explanation");

        assert!(normalization.contains("normalized representative"));
        assert!(roundtrip.contains("affine image: (2, 5)"));
    }

    #[test]
    fn cost_description_lists_representation_work_and_auxiliary_counts() {
        let cost = ShortWeierstrassProjectiveOperationCost::for_kind(
            ShortWeierstrassProjectiveOperationKind::Add,
        );
        let normalization = ShortWeierstrassProjectiveOperationCost::for_kind(
            ShortWeierstrassProjectiveOperationKind::Normalize,
        );
        let description = describe_short_weierstrass_projective_cost(&cost);

        assert!(description.contains("Short-Weierstrass projective operation cost"));
        assert!(description.contains("delegated affine additions: 0"));
        assert_eq!(normalization.representation_cost().inversions(), 1);
    }

    #[test]
    fn general_cost_description_uses_the_general_title() {
        let cost = GeneralWeierstrassProjectiveOperationCost::for_kind(
            GeneralWeierstrassProjectiveOperationKind::Neg,
        );
        let description = describe_general_weierstrass_projective_cost(&cost);

        assert!(description.contains("General-Weierstrass projective operation cost"));
        assert!(description.contains("operation: Neg"));
    }
}
