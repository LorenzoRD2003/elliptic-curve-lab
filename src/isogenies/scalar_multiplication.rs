use crate::elliptic_curves::{FiniteGroupCurveModel, GroupCurveModel};
use crate::fields::{EnumerableFiniteField, SqrtField};

use super::{Isogeny, IsogenyError};

/// Scalar-multiplication isogeny `[n] : E -> E` on a small explicit curve group.
///
/// For a non-zero integer `n`, elliptic-curve multiplication by `n`
///
/// `[n](P) = P + P + ... + P`
///
/// is a standard example of an isogeny from a curve to itself.
///
/// In this educational implementation:
///
/// - the domain and codomain are the same curve value
/// - the degree is reported as `n^2`
/// - `kernel_points()` exposes the points of `E(F_q)` killed by `[n]`
pub struct ScalarMultiplicationIsogeny<C>
where
    C: GroupCurveModel,
    C::Point: Clone,
{
    curve: C,
    scalar: u64,
    kernel_points: Vec<C::Point>,
}

impl<C> ScalarMultiplicationIsogeny<C>
where
    C: FiniteGroupCurveModel,
    C::BaseField: EnumerableFiniteField<Elem = C::Elem> + SqrtField<Elem = C::Elem>,
    C::Point: Clone + PartialEq,
{
    /// Builds the scalar-multiplication isogeny `[n]`.
    ///
    /// The current constructor is intentionally restricted to small finite
    /// curve groups so it can materialize the explicit rational kernel
    ///
    /// `{ P in E(F_q) : [n]P = O }`.
    ///
    /// Scalar `0` is rejected because this crate reserves the isogeny surface
    /// for the usual non-constant multiplication-by-`n` maps.
    pub fn new(curve: C, scalar: u64) -> Result<Self, IsogenyError> {
        if scalar == 0 {
            return Err(IsogenyError::ZeroScalarIsNotIsogeny);
        }

        let identity = curve.identity();
        let kernel_points = curve
            .points()
            .into_iter()
            .map(|point| -> Result<_, IsogenyError> {
                let image = curve.mul_scalar(&point, scalar)?;
                Ok((point, image == identity))
            })
            .collect::<Result<Vec<_>, IsogenyError>>()?
            .into_iter()
            .filter_map(|(point, kills_point)| kills_point.then_some(point))
            .collect();

        Ok(Self {
            curve,
            scalar,
            kernel_points,
        })
    }

    /// Returns the underlying scalar `n` in `[n]`.
    pub fn scalar(&self) -> u64 {
        self.scalar
    }
}

impl<C> Isogeny<C, C> for ScalarMultiplicationIsogeny<C>
where
    C: GroupCurveModel,
    C::Point: Clone,
{
    fn domain(&self) -> &C {
        &self.curve
    }

    fn codomain(&self) -> &C {
        &self.curve
    }

    fn degree(&self) -> usize {
        usize::try_from(u128::from(self.scalar) * u128::from(self.scalar))
            .expect("educational scalar-multiplication degrees should fit in usize")
    }

    fn evaluate(&self, p: &C::Point) -> Result<C::Point, IsogenyError> {
        self.curve.mul_scalar(p, self.scalar).map_err(Into::into)
    }

    fn kernel_points(&self) -> &[C::Point] {
        &self.kernel_points
    }
}

#[cfg(test)]
mod tests {
    use crate::elliptic_curves::{
        AffineCurveModel, CurveModel, EnumerableCurveModel, FiniteGroupCurveModel, GroupCurveModel,
        ShortWeierstrassCurve,
    };
    use crate::fields::{Field, Fp};
    use crate::isogenies::{Isogeny, IsogenyError, ScalarMultiplicationIsogeny, VerifiableIsogeny};

    type F41 = Fp<41>;
    type Curve = ShortWeierstrassCurve<F41>;

    fn curve() -> Curve {
        Curve::new(F41::from_i64(2), F41::from_i64(3)).expect("valid curve")
    }

    #[test]
    fn constructor_rejects_zero_scalar() {
        assert!(matches!(
            ScalarMultiplicationIsogeny::new(curve(), 0),
            Err(IsogenyError::ZeroScalarIsNotIsogeny)
        ));
    }

    #[test]
    fn degree_of_multiplication_by_two_is_four() {
        let isogeny =
            ScalarMultiplicationIsogeny::new(curve(), 2).expect("scalar isogeny should build");

        assert_eq!(isogeny.degree(), 4);
        assert_eq!(isogeny.scalar(), 2);
    }

    #[test]
    fn evaluation_matches_group_scalar_multiplication() {
        let curve = curve();
        let point = curve
            .point(F41::from_i64(3), F41::from_i64(6))
            .expect("sample point should lie on the curve");
        let isogeny = ScalarMultiplicationIsogeny::new(curve.clone(), 3)
            .expect("scalar isogeny should build");

        assert_eq!(
            isogeny.evaluate(&point),
            curve.mul_scalar(&point, 3).map_err(Into::into)
        );
    }

    #[test]
    fn scalar_one_is_identity_map() {
        let curve = curve();
        let isogeny = ScalarMultiplicationIsogeny::new(curve.clone(), 1)
            .expect("scalar-one isogeny should build");

        for point in curve.points() {
            assert_eq!(
                isogeny
                    .evaluate(&point)
                    .expect("scalar-one isogeny should evaluate"),
                point
            );
        }
    }

    #[test]
    fn kernel_points_match_the_rational_two_torsion_plus_identity() {
        let curve = curve();
        let isogeny = ScalarMultiplicationIsogeny::new(curve.clone(), 2)
            .expect("scalar isogeny should build");

        let mut expected = vec![curve.identity()];
        expected.extend(curve.points_of_order(2));

        assert_eq!(isogeny.kernel_points(), expected.as_slice());
    }

    #[test]
    fn exhaustive_verifier_passes_for_multiplication_by_two() {
        let isogeny =
            ScalarMultiplicationIsogeny::new(curve(), 2).expect("scalar isogeny should build");

        assert_eq!(isogeny.verify_maps_domain_to_codomain(), Ok(()));
        assert_eq!(isogeny.verify_maps_kernel_to_identity(), Ok(()));
        assert_eq!(isogeny.verify_homomorphism(), Ok(()));
        assert_eq!(isogeny.verify_kernel_exactness(), Ok(()));
    }
}
