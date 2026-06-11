use std::hash::Hash;

use crate::elliptic_curves::isomorphisms::{CurveIsomorphism, ShortWeierstrassIsomorphism};
use crate::elliptic_curves::short_weierstrass::ShortWeierstrassCurve;
use crate::elliptic_curves::traits::{CurveModel, EnumerableCurveModel, FiniteGroupCurveModel};
use crate::fields::{EnumerableFiniteField, Field, SqrtField};
use crate::isogenies::{
    DualIsogenyError, Isogeny, IsogenyError, KernelDescription, ScalarMultiplicationIsogeny,
    VeluIsogeny,
};

/// Exhaustively searched dual of a short-Weierstrass Vélu isogeny.
///
/// This is represented as `alpha ∘ psi`, where:
///
/// - `psi` is a Vélu isogeny built from a candidate kernel on `E'(F_q)`
/// - `alpha` is a base-field short-Weierstrass isomorphism from the Vélu
///   codomain back to the original domain curve `E`
pub struct DualVeluIsogeny<F: Field> {
    velu_part: VeluIsogeny<ShortWeierstrassCurve<F>>,
    codomain_to_original: ShortWeierstrassIsomorphism<F>,
}

impl<F> DualVeluIsogeny<F>
where
    F: Field,
{
    pub fn velu_part(&self) -> &VeluIsogeny<ShortWeierstrassCurve<F>> {
        &self.velu_part
    }

    pub fn codomain_to_original(&self) -> &ShortWeierstrassIsomorphism<F> {
        &self.codomain_to_original
    }
}

impl<F> Isogeny<ShortWeierstrassCurve<F>, ShortWeierstrassCurve<F>> for DualVeluIsogeny<F>
where
    F: Field + Clone,
    F::Elem: Clone + Eq + Hash,
{
    fn domain(&self) -> &ShortWeierstrassCurve<F> {
        self.velu_part.domain()
    }

    fn codomain(&self) -> &ShortWeierstrassCurve<F> {
        self.codomain_to_original.codomain()
    }

    fn degree(&self) -> usize {
        self.velu_part.degree()
    }

    fn evaluate(
        &self,
        point: &<ShortWeierstrassCurve<F> as CurveModel>::Point,
    ) -> Result<<ShortWeierstrassCurve<F> as CurveModel>::Point, IsogenyError> {
        let mid = self.velu_part.evaluate(point)?;
        self.codomain_to_original.evaluate(&mid).map_err(Into::into)
    }

    fn kernel_description(&self) -> KernelDescription<ShortWeierstrassCurve<F>> {
        self.velu_part.kernel_description()
    }
}

impl<F> VeluIsogeny<ShortWeierstrassCurve<F>>
where
    F: Field + EnumerableFiniteField + SqrtField + Clone,
    F::Elem: Clone + Eq + Hash,
{
    /// Searches exhaustively for a dual isogeny on a small finite curve.
    ///
    /// For a separable Vélu isogeny
    /// `phi : E -> E'`
    /// of degree `n`, the current educational strategy is:
    ///
    /// 1. enumerate points of exact order `n` on `E'(F_q)`
    /// 2. build the corresponding Vélu isogenies `psi : E' -> C`
    /// 3. search for a base-field isomorphism `alpha : C -> E`
    /// 4. form the candidate dual `alpha ∘ psi`
    /// 5. verify the duality relations exhaustively on rational points
    pub fn find_dual_exhaustively(&self) -> Result<DualVeluIsogeny<F>, IsogenyError> {
        let degree = self.degree();
        let mut saw_isomorphic_candidate = false;

        for generator in self.codomain().points_of_order(degree) {
            let psi = VeluIsogeny::from_generator(self.codomain().clone(), generator)?;
            if psi.degree() != degree {
                return Err(IsogenyError::Dual(DualIsogenyError::DegreeMismatch));
            }

            let alphas = exhaustive_isomorphisms_to_target(psi.codomain(), self.domain());
            if alphas.is_empty() {
                continue;
            }
            saw_isomorphic_candidate = true;

            for alpha in alphas {
                let candidate_dual = DualVeluIsogeny {
                    velu_part: psi.clone(),
                    codomain_to_original: alpha,
                };

                if dual_relations_hold_exhaustively(self, &candidate_dual)? {
                    if candidate_dual.degree() != degree {
                        return Err(IsogenyError::Dual(DualIsogenyError::DegreeMismatch));
                    }

                    return Ok(candidate_dual);
                }
            }
        }

        if saw_isomorphic_candidate {
            Err(IsogenyError::Dual(DualIsogenyError::DualRelationViolation))
        } else {
            Err(IsogenyError::Dual(DualIsogenyError::DualNotFound))
        }
    }
}

/// Verifies exhaustively that `\hat{\varphi} \circ \varphi = [n]_E` on
/// rational points of the original domain curve.
///
/// This helper is intentionally small-field and educational: it enumerates all
/// points of `E(F_q)` and checks that applying the candidate dual after `phi`
/// agrees pointwise with multiplication by `n = deg(phi)`.
pub fn verify_left_dual_relation<F>(
    phi: &VeluIsogeny<ShortWeierstrassCurve<F>>,
    dual: &DualVeluIsogeny<F>,
) -> Result<(), IsogenyError>
where
    F: Field + EnumerableFiniteField + SqrtField + Clone,
    F::Elem: Clone + Eq + Hash,
{
    let degree_u64 =
        u64::try_from(phi.degree()).expect("tiny educational isogeny degrees should fit in u64");
    let scalar = ScalarMultiplicationIsogeny::new(phi.domain().clone(), degree_u64)?;

    for point in phi.domain().points() {
        let left = dual.evaluate(&phi.evaluate(&point)?)?;
        let right = scalar.evaluate(&point)?;
        if left != right {
            return Err(IsogenyError::Dual(DualIsogenyError::DualRelationViolation));
        }
    }

    Ok(())
}

/// Verifies exhaustively that `\varphi \circ \hat{\varphi} = [n]_{E'}`
/// on rational points of the codomain curve.
///
/// This helper is intentionally small-field and educational: it enumerates all
/// points of `E'(F_q)` and checks that applying `phi` after the candidate dual
/// agrees pointwise with multiplication by `n = deg(phi)`.
pub fn verify_right_dual_relation<F>(
    phi: &VeluIsogeny<ShortWeierstrassCurve<F>>,
    dual: &DualVeluIsogeny<F>,
) -> Result<(), IsogenyError>
where
    F: Field + EnumerableFiniteField + SqrtField + Clone,
    F::Elem: Clone + Eq + Hash,
{
    let degree_u64 =
        u64::try_from(phi.degree()).expect("tiny educational isogeny degrees should fit in u64");
    let scalar = ScalarMultiplicationIsogeny::new(phi.codomain().clone(), degree_u64)?;

    for point in phi.codomain().points() {
        let left = phi.evaluate(&dual.evaluate(&point)?)?;
        let right = scalar.evaluate(&point)?;
        if left != right {
            return Err(IsogenyError::Dual(DualIsogenyError::DualRelationViolation));
        }
    }

    Ok(())
}

fn dual_relations_hold_exhaustively<F>(
    phi: &VeluIsogeny<ShortWeierstrassCurve<F>>,
    candidate_dual: &DualVeluIsogeny<F>,
) -> Result<bool, IsogenyError>
where
    F: Field + EnumerableFiniteField + SqrtField + Clone,
    F::Elem: Clone + Eq + Hash,
{
    match verify_left_dual_relation(phi, candidate_dual) {
        Ok(()) => {}
        Err(IsogenyError::Dual(DualIsogenyError::DualRelationViolation)) => return Ok(false),
        Err(other) => return Err(other),
    }

    match verify_right_dual_relation(phi, candidate_dual) {
        Ok(()) => {}
        Err(IsogenyError::Dual(DualIsogenyError::DualRelationViolation)) => return Ok(false),
        Err(other) => return Err(other),
    }

    Ok(true)
}

fn exhaustive_isomorphisms_to_target<F>(
    source: &ShortWeierstrassCurve<F>,
    target: &ShortWeierstrassCurve<F>,
) -> Vec<ShortWeierstrassIsomorphism<F>>
where
    F: Field + EnumerableFiniteField + Clone,
    F::Elem: Clone + Eq + Hash,
{
    F::elements()
        .into_iter()
        .filter_map(|u| ShortWeierstrassIsomorphism::new(source.clone(), u).ok())
        .filter(|isomorphism| isomorphism.codomain() == target)
        .collect()
}

#[cfg(test)]
mod tests {
    use crate::elliptic_curves::{AffineCurveModel, EnumerableCurveModel, ShortWeierstrassCurve};
    use crate::fields::{Field, Fp};
    use crate::isogenies::velu::short_weierstrass::{
        DualVeluIsogeny, verify_left_dual_relation, verify_right_dual_relation,
    };
    use crate::isogenies::{Isogeny, ScalarMultiplicationIsogeny, VeluIsogeny};

    type F41 = Fp<41>;
    type F29 = Fp<29>;
    type Curve = ShortWeierstrassCurve<F41>;
    type CurveF29 = ShortWeierstrassCurve<F29>;
    type Dual = DualVeluIsogeny<F41>;
    type DualF29 = DualVeluIsogeny<F29>;

    fn curve() -> Curve {
        Curve::new(F41::from_i64(2), F41::from_i64(3)).expect("valid curve")
    }

    fn sample_phi() -> VeluIsogeny<Curve> {
        let curve = curve();
        let generator = curve
            .point(F41::from_i64(40), F41::from_i64(0))
            .expect("sample generator should lie on the curve");
        VeluIsogeny::from_generator(curve, generator).expect("sample Vélu isogeny should build")
    }

    fn curve_f29() -> CurveF29 {
        CurveF29::new(F29::from_i64(2), F29::from_i64(2)).expect("valid F29 curve")
    }

    fn sample_degree_three_phi() -> VeluIsogeny<CurveF29> {
        let curve = curve_f29();
        let generator = curve
            .point(F29::from_i64(10), F29::from_i64(23))
            .expect("sample degree-three generator should lie on the curve");
        VeluIsogeny::from_generator(curve, generator)
            .expect("sample degree-three Vélu isogeny should build")
    }

    fn left_relation_holds(phi: &VeluIsogeny<Curve>, dual: &Dual) -> bool {
        let scalar =
            ScalarMultiplicationIsogeny::new(phi.domain().clone(), phi.degree() as u64).unwrap();

        phi.domain().points().into_iter().all(|point| {
            dual.evaluate(&phi.evaluate(&point).unwrap()).unwrap()
                == scalar.evaluate(&point).unwrap()
        })
    }

    fn right_relation_holds(phi: &VeluIsogeny<Curve>, dual: &Dual) -> bool {
        let scalar =
            ScalarMultiplicationIsogeny::new(phi.codomain().clone(), phi.degree() as u64).unwrap();

        phi.codomain().points().into_iter().all(|point| {
            phi.evaluate(&dual.evaluate(&point).unwrap()).unwrap()
                == scalar.evaluate(&point).unwrap()
        })
    }

    #[test]
    fn dual_search_finds_a_degree_matching_candidate_on_the_f41_example() {
        let phi = sample_phi();
        let dual = phi
            .find_dual_exhaustively()
            .expect("dual should be found on the sample curve");

        assert_eq!(dual.degree(), phi.degree());
    }

    #[test]
    fn dual_search_candidate_satisfies_both_duality_relations_on_rational_points() {
        let phi = sample_phi();
        let dual = phi
            .find_dual_exhaustively()
            .expect("dual should be found on the sample curve");

        assert!(left_relation_holds(&phi, &dual));
        assert!(right_relation_holds(&phi, &dual));
    }

    #[test]
    fn public_left_dual_relation_verifier_accepts_the_enumerated_dual() {
        let phi = sample_phi();
        let dual = phi
            .find_dual_exhaustively()
            .expect("dual should be found on the sample curve");

        verify_left_dual_relation(&phi, &dual)
            .expect("the exhaustively found dual should satisfy the left relation");
    }

    #[test]
    fn public_right_dual_relation_verifier_accepts_the_enumerated_dual() {
        let phi = sample_phi();
        let dual = phi
            .find_dual_exhaustively()
            .expect("dual should be found on the sample curve");

        verify_right_dual_relation(&phi, &dual)
            .expect("the exhaustively found dual should satisfy the right relation");
    }

    #[test]
    fn find_dual_by_enumeration_finds_dual_for_degree_3_example() {
        let phi = sample_degree_three_phi();

        let dual = phi
            .find_dual_exhaustively()
            .expect("dual should be found for the degree-three sample");

        assert_eq!(dual.degree(), 3);
    }

    #[test]
    fn dual_kernel_has_expected_order() {
        let phi = sample_degree_three_phi();
        let dual: DualF29 = phi
            .find_dual_exhaustively()
            .expect("dual should be found for the degree-three sample");

        assert_eq!(dual.kernel_points().len(), phi.degree());
    }
}
