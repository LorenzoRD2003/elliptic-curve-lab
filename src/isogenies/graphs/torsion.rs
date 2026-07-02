use core::hash::Hash;

use crate::elliptic_curves::traits::FiniteGroupCurveModel;
use crate::fields::traits::{EnumerableFiniteField, SqrtField};
use crate::isogenies::{
    graphs::{GraphCurveModel, IsogenyGraphError},
    kernel::IsogenyKernel,
};

/// Graph-side torsion helpers for small finite curve models.
///
/// This trait stays narrower than the base curve-model traits. It packages the
/// graph-specific torsion operation we actually need when building educational
/// `\ell`-isogeny graphs:
///
/// - enumerate rational points of exact order `\ell`
/// - collapse generators `P` and `-P` to the same cyclic subgroup
/// - return explicit finite kernels suitable for later Vélu construction
pub(crate) trait GraphTorsionCurveModel: FiniteGroupCurveModel + GraphCurveModel
where
    Self::BaseField: EnumerableFiniteField<Elem = Self::Elem> + SqrtField<Elem = Self::Elem>,
    Self::Point: Clone + Eq + Hash,
{
    /// Returns the distinct cyclic kernels of exact order `ell`.
    ///
    /// This is the graph-layer wrapper over the generic exact-order helpers on
    /// [`crate::elliptic_curves::traits::FiniteGroupCurveModel`].
    fn cyclic_kernels_of_order(
        &self,
        ell: usize,
    ) -> Result<Vec<IsogenyKernel<Self>>, IsogenyGraphError>
    where
        Self: Sized,
    {
        let points = self
            .points_of_exact_order(ell)
            .map_err(|error| match error {
                crate::elliptic_curves::CurveError::InvalidTorsionOrder { .. } => {
                    IsogenyGraphError::InvalidDegree
                }
                other => IsogenyGraphError::from(other),
            })?;

        let mut kernels = Vec::new();

        for point in points {
            let kernel = IsogenyKernel::cyclic(self, &point)?;

            if kernel.order() != ell {
                continue;
            }

            if !kernels.iter().any(|existing| existing == &kernel) {
                kernels.push(kernel);
            }
        }

        Ok(kernels)
    }
}

impl<C> GraphTorsionCurveModel for C
where
    C: FiniteGroupCurveModel + GraphCurveModel,
    C::BaseField: EnumerableFiniteField<Elem = C::Elem> + SqrtField<Elem = C::Elem>,
    C::Point: Clone + Eq + Hash,
{
}

#[cfg(test)]
mod tests {

    use crate::elliptic_curves::{
        ShortWeierstrassCurve,
        traits::{AffineCurveModel, CurveModel, FiniteGroupCurveModel},
    };
    use crate::fields::traits::Field;
    use crate::isogenies::{
        graphs::{GraphTorsionCurveModel, IsogenyGraphError},
        kernel::IsogenyKernel,
        traits::{Isogeny, VerifiableIsogeny},
        velu::VeluIsogeny,
    };

    type F5 = crate::fields::Fp5;
    type F7 = crate::fields::Fp7;
    type F41 = crate::fields::Fp41;
    type Curve5 = ShortWeierstrassCurve<F5>;
    type Curve7 = ShortWeierstrassCurve<F7>;
    type Curve = ShortWeierstrassCurve<F41>;

    fn f5_noncyclic_curve() -> Curve5 {
        Curve5::new(F5::from_i64(-1), F5::zero()).expect("valid curve")
    }

    fn f7_curve() -> Curve7 {
        Curve7::new(F7::from_i64(2), F7::from_i64(3)).expect("valid curve")
    }

    fn f41_curve() -> Curve {
        Curve::new(F41::from_i64(2), F41::from_i64(3)).expect("valid curve")
    }

    #[test]
    fn cyclic_kernels_reject_zero_degree() {
        let curve = f41_curve();

        assert_eq!(
            curve.cyclic_kernels_of_order(0),
            Err(IsogenyGraphError::InvalidDegree)
        );
    }

    #[test]
    fn exact_order_points_generate_kernels_of_size_l() {
        let curve = f41_curve();
        let kernels = curve
            .cyclic_kernels_of_order(2)
            .expect("order two kernels should exist");

        assert_eq!(kernels.len(), 1);
        assert!(kernels.iter().all(|kernel| kernel.order() == 2));
    }

    #[test]
    fn p_and_minus_p_generate_same_kernel_once() {
        let curve = f7_curve();
        let kernels = curve
            .cyclic_kernels_of_order(3)
            .expect("order three kernels should exist");

        assert_eq!(kernels.len(), 1);
    }

    #[test]
    fn cyclic_kernels_of_order_l_have_size_l() {
        let curve = f7_curve();
        let kernels = curve
            .cyclic_kernels_of_order(3)
            .expect("order three kernels should exist");

        assert!(kernels.iter().all(|kernel| kernel.order() == 3));
    }

    #[test]
    fn number_of_kernels_matches_expected_small_example() {
        let curve = f7_curve();
        let kernels = curve
            .cyclic_kernels_of_order(3)
            .expect("order three kernels should exist");

        assert_eq!(kernels.len(), 1);
    }

    #[test]
    fn full_rational_two_torsion_has_l_plus_one_kernels() {
        let curve = f5_noncyclic_curve();
        let kernels = curve
            .cyclic_kernels_of_order(2)
            .expect("order two kernels should exist");

        assert_eq!(curve.points_of_order(2).len(), 3);
        assert_eq!(kernels.len(), 3);
    }

    #[test]
    fn cyclic_kernel_matches_known_order_two_generator() {
        let curve = f41_curve();
        let expected_point = curve
            .point(F41::from_i64(40), F41::from_i64(0))
            .expect("sample point should lie on the curve");
        let kernels = curve
            .cyclic_kernels_of_order(2)
            .expect("order two kernels should exist");

        assert_eq!(kernels[0].points(), &[curve.identity(), expected_point]);
    }

    #[test]
    fn graph_cyclic_kernels_of_order_three_agree_with_division_polynomial_torsion() {
        let curve = f7_curve();
        let exact_points = curve
            .exact_n_torsion_points_from_division_polynomial(3)
            .expect("division-polynomial exact three-torsion should exist");
        let graph_kernels = curve
            .cyclic_kernels_of_order(3)
            .expect("graph kernels of order three should exist");
        let mut polynomial_kernels = Vec::new();

        for point in exact_points {
            let kernel = IsogenyKernel::cyclic(&curve, &point)
                .expect("exact three-torsion point should generate a cyclic kernel");
            if !polynomial_kernels
                .iter()
                .any(|existing| existing == &kernel)
            {
                polynomial_kernels.push(kernel);
            }
        }

        assert_eq!(polynomial_kernels.len(), graph_kernels.len());
        assert!(
            polynomial_kernels
                .iter()
                .all(|kernel| graph_kernels.iter().any(|other| other == kernel))
        );
    }

    #[test]
    fn division_polynomial_torsion_can_feed_velu_kernel_construction() {
        let curve = f7_curve();
        let generator = curve
            .exact_n_torsion_points_from_division_polynomial(3)
            .expect("division-polynomial exact three-torsion should exist")
            .into_iter()
            .next()
            .expect("the sample curve should have a rational three-torsion generator");

        let isogeny =
            VeluIsogeny::from_generator(curve.clone(), generator).expect("Vélu should build");

        assert_eq!(isogeny.degree(), 3);
        assert_eq!(isogeny.kernel().order(), 3);
        assert_eq!(isogeny.verify_maps_kernel_to_identity(), Ok(()));
        assert_eq!(isogeny.verify_kernel_exactness(), Ok(()));
    }
}
