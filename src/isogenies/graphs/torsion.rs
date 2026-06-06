use std::hash::Hash;

use crate::elliptic_curves::torsion::points_of_exact_order;
use crate::elliptic_curves::traits::FiniteGroupCurveModel;
use crate::fields::{EnumerableFiniteField, SqrtField};
use crate::isogenies::IsogenyKernel;
use crate::isogenies::graphs::GraphCurveModel;
use crate::isogenies::graphs::IsogenyGraphError;

/// Returns the distinct cyclic kernels of exact order `ell`.
///
/// This is the graph-layer wrapper over the generic elliptic-curve torsion
/// helpers in `elliptic_curves::torsion`. It:
///
/// - asks the curve layer for the rational points of exact order `ell`
/// - turns each such point into an explicit cyclic subgroup
/// - deduplicates generators such as `P` and `-P` when they define the same
///   kernel
pub fn cyclic_kernels_of_order<C>(
    curve: &C,
    ell: usize,
) -> Result<Vec<IsogenyKernel<C>>, IsogenyGraphError>
where
    C: FiniteGroupCurveModel + GraphCurveModel,
    C::BaseField: EnumerableFiniteField<Elem = C::Elem> + SqrtField<Elem = C::Elem>,
    C::Point: Clone + Eq + Hash,
{
    let points = points_of_exact_order(curve, ell).map_err(|error| match error {
        crate::elliptic_curves::CurveError::InvalidTorsionOrder { .. } => {
            IsogenyGraphError::InvalidDegree
        }
        other => IsogenyGraphError::from(other),
    })?;

    let mut kernels = Vec::new();

    for point in points {
        let kernel = IsogenyKernel::cyclic(curve, &point)?;

        if kernel.order() != ell {
            continue;
        }

        if !kernels.iter().any(|existing| existing == &kernel) {
            kernels.push(kernel);
        }
    }

    Ok(kernels)
}

#[cfg(test)]
mod tests {
    use crate::elliptic_curves::{
        AffineCurveModel, CurveModel, FiniteGroupCurveModel, ShortWeierstrassCurve,
        division_polynomials::exact_n_torsion_points_from_division_polynomial,
    };
    use crate::fields::{Field, Fp};
    use crate::isogenies::graphs::IsogenyGraphError;
    use crate::isogenies::graphs::cyclic_kernels_of_order;
    use crate::isogenies::{Isogeny, VeluIsogeny, VerifiableIsogeny};

    type F5 = Fp<5>;
    type F7 = Fp<7>;
    type F41 = Fp<41>;
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
            cyclic_kernels_of_order(&curve, 0),
            Err(IsogenyGraphError::InvalidDegree)
        );
    }

    #[test]
    fn exact_order_points_generate_kernels_of_size_l() {
        let curve = f41_curve();
        let kernels = cyclic_kernels_of_order(&curve, 2).expect("order two kernels should exist");

        assert_eq!(kernels.len(), 1);
        assert!(kernels.iter().all(|kernel| kernel.order() == 2));
    }

    #[test]
    fn p_and_minus_p_generate_same_kernel_once() {
        let curve = f7_curve();
        let kernels = cyclic_kernels_of_order(&curve, 3).expect("order three kernels should exist");

        assert_eq!(kernels.len(), 1);
    }

    #[test]
    fn cyclic_kernels_of_order_l_have_size_l() {
        let curve = f7_curve();
        let kernels = cyclic_kernels_of_order(&curve, 3).expect("order three kernels should exist");

        assert!(kernels.iter().all(|kernel| kernel.order() == 3));
    }

    #[test]
    fn number_of_kernels_matches_expected_small_example() {
        let curve = f7_curve();
        let kernels = cyclic_kernels_of_order(&curve, 3).expect("order three kernels should exist");

        assert_eq!(kernels.len(), 1);
    }

    #[test]
    fn full_rational_two_torsion_has_l_plus_one_kernels() {
        let curve = f5_noncyclic_curve();
        let kernels = cyclic_kernels_of_order(&curve, 2).expect("order two kernels should exist");

        assert_eq!(curve.points_of_order(2).len(), 3);
        assert_eq!(kernels.len(), 3);
    }

    #[test]
    fn cyclic_kernel_matches_known_order_two_generator() {
        let curve = f41_curve();
        let expected_point = curve
            .point(F41::from_i64(40), F41::from_i64(0))
            .expect("sample point should lie on the curve");
        let kernels = cyclic_kernels_of_order(&curve, 2).expect("order two kernels should exist");

        assert_eq!(kernels[0].points(), &[curve.identity(), expected_point]);
    }

    #[test]
    fn graph_cyclic_kernels_of_order_three_agree_with_division_polynomial_torsion() {
        let curve = f7_curve();
        let exact_points = exact_n_torsion_points_from_division_polynomial(&curve, 3)
            .expect("division-polynomial exact three-torsion should exist");
        let graph_kernels =
            cyclic_kernels_of_order(&curve, 3).expect("graph kernels of order three should exist");
        let mut polynomial_kernels = Vec::new();

        for point in exact_points {
            let kernel = crate::isogenies::IsogenyKernel::cyclic(&curve, &point)
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
        let generator = exact_n_torsion_points_from_division_polynomial(&curve, 3)
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
