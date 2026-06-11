use crate::elliptic_curves::traits::{CurveModel, FiniteGroupCurveModel};
use crate::fields::{EnumerableFiniteField, SqrtField};

use crate::isogenies::{Isogeny, IsogenyError, IsogenyMapError};

/// Exhaustively compares two explicit maps on a small finite domain curve.
///
/// This helper is intentionally educational and brute-force. It first checks
/// that both maps have the same concrete domain and codomain curves, then
/// evaluates both maps on every point of `C(F_q)` and returns whether the
/// images agree pointwise.
///
/// It is meant for the same small finite settings as [`FiniteGroupCurveModel`]
/// and provides a simple confidence-building notion of map equality:
///
/// `left = right` exactly when `left(P) = right(P)` for every enumerated
/// domain point `P`.
pub fn maps_equal_exhaustively<I, J, C, D>(left: &I, right: &J) -> Result<bool, IsogenyError>
where
    C: FiniteGroupCurveModel + PartialEq,
    C::BaseField: EnumerableFiniteField<Elem = C::Elem> + SqrtField<Elem = C::Elem>,
    C::Point: Clone + PartialEq,
    D: CurveModel + PartialEq,
    D::Point: PartialEq,
    I: Isogeny<C, D>,
    J: Isogeny<C, D>,
{
    if left.domain() != right.domain() || left.codomain() != right.codomain() {
        return Err(IsogenyError::Map(
            IsogenyMapError::MapComparisonDomainCodomainMismatch,
        ));
    }

    for point in left.domain().points() {
        if left.evaluate(&point)? != right.evaluate(&point)? {
            return Ok(false);
        }
    }

    Ok(true)
}

#[cfg(test)]
mod tests {
    use crate::elliptic_curves::{AffineCurveModel, ShortWeierstrassCurve};
    use crate::fields::{Field, Fp};
    use crate::isogenies::{
        IsogenyError, IsogenyMapError, ScalarMultiplicationIsogeny, VeluIsogeny,
        maps_equal_exhaustively,
    };

    type F41 = Fp<41>;
    type Curve = ShortWeierstrassCurve<F41>;

    fn curve_a() -> Curve {
        Curve::new(F41::from_i64(2), F41::from_i64(3)).expect("valid curve")
    }

    fn curve_b() -> Curve {
        Curve::new(F41::from_i64(18), F41::from_i64(38)).expect("valid curve")
    }

    #[test]
    fn identical_scalar_multiplication_maps_are_equal_exhaustively() {
        let left = ScalarMultiplicationIsogeny::new(curve_a(), 2).expect("map should build");
        let right = ScalarMultiplicationIsogeny::new(curve_a(), 2).expect("map should build");

        assert_eq!(
            maps_equal_exhaustively::<_, _, Curve, Curve>(&left, &right),
            Ok(true)
        );
    }

    #[test]
    fn different_scalar_multiplication_maps_are_not_equal_exhaustively() {
        let left = ScalarMultiplicationIsogeny::new(curve_a(), 2).expect("map should build");
        let right = ScalarMultiplicationIsogeny::new(curve_a(), 3).expect("map should build");

        assert_eq!(
            maps_equal_exhaustively::<_, _, Curve, Curve>(&left, &right),
            Ok(false)
        );
    }

    #[test]
    fn map_comparison_rejects_different_concrete_domain_or_codomain_curves() {
        let left = ScalarMultiplicationIsogeny::new(curve_a(), 2).expect("map should build");
        let right = ScalarMultiplicationIsogeny::new(curve_b(), 2).expect("map should build");

        assert_eq!(
            maps_equal_exhaustively::<_, _, Curve, Curve>(&left, &right),
            Err(IsogenyError::Map(
                IsogenyMapError::MapComparisonDomainCodomainMismatch
            ))
        );
    }

    #[test]
    fn identical_velu_isogenies_are_equal_exhaustively() {
        let domain = curve_a();
        let generator = domain
            .point(F41::from_i64(40), F41::from_i64(0))
            .expect("sample generator should lie on the curve");
        let left = VeluIsogeny::from_generator(domain.clone(), generator.clone())
            .expect("Vélu map should build");
        let right = VeluIsogeny::from_generator(domain, generator).expect("Vélu map should build");

        assert_eq!(
            maps_equal_exhaustively::<_, _, Curve, Curve>(&left, &right),
            Ok(true)
        );
    }
}
