use crate::elliptic_curves::traits::{EnumerableCurveModel, GroupCurveModel};
use crate::elliptic_curves::{AffinePoint, ShortWeierstrassCurve};
use crate::fields::traits::*;
use crate::fields::traits::{EnumerableFiniteField, SqrtField};

pub(super) fn same_x_set<F: Field>(left: &[F::Elem], right: &[F::Elem]) -> bool {
    left.len() == right.len()
        && left.iter().all(|x| right.iter().any(|y| F::eq(x, y)))
        && right.iter().all(|x| left.iter().any(|y| F::eq(x, y)))
}

pub(super) fn same_point_set<F: Field>(left: &[AffinePoint<F>], right: &[AffinePoint<F>]) -> bool {
    left.len() == right.len()
        && left
            .iter()
            .all(|point| right.iter().any(|other| point == other))
        && right
            .iter()
            .all(|point| left.iter().any(|other| point == other))
}

pub(super) fn unique_x_coordinates_of_rational_n_torsion<F: EnumerableFiniteField + SqrtField>(
    curve: &ShortWeierstrassCurve<F>,
    n: usize,
) -> Vec<F::Elem> {
    let mut xs = Vec::new();
    for point in curve.points() {
        let AffinePoint::Finite { ref x, .. } = point else {
            continue;
        };
        if !curve.is_torsion_point(&point, n) {
            continue;
        }
        if xs.iter().any(|existing| F::eq(existing, x)) {
            continue;
        }
        xs.push(x.clone());
    }
    xs
}
