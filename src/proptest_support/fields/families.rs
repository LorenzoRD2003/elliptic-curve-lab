use crate::fields::traits::*;
use proptest::prelude::*;

use crate::fields::extension_field::{ExtensionField, ExtensionFieldElement, ExtensionFieldSpec};
use crate::fields::polynomial_field::PolynomialModulus;
use crate::proptest_support::fields::{arb_fp_elem, arb_semantic_extension_elem};

type F17 = crate::fields::Fp17;

/// Canonical quadratic extension fixture used across the repository's tests.
pub struct ProptestF17Sqrt3Spec;

impl ExtensionFieldSpec for ProptestF17Sqrt3Spec {
    type Base = F17;

    const NAME: &'static str = "proptest F17(sqrt(3))";

    fn defining_modulus() -> PolynomialModulus<Self::Base> {
        PolynomialModulus::<Self::Base>::new(vec![
            <Self::Base as Field>::from_i64(-3),
            <Self::Base as Field>::zero(),
            <Self::Base as Field>::one(),
        ])
        .expect("x^2 - 3 should be a valid structural modulus")
    }

    fn check_field_conditions() -> Result<(), crate::fields::FieldError> {
        Self::defining_modulus().check_field_modulus_requirements()
    }
}

/// Canonical degree-three tower step above `F17(sqrt(3))`.
pub struct ProptestF17TowerSpec;

impl ExtensionFieldSpec for ProptestF17TowerSpec {
    type Base = ProptestF17Sqrt3Field;

    const NAME: &'static str = "proptest F17(sqrt(3))(u)";

    fn defining_modulus() -> PolynomialModulus<Self::Base> {
        PolynomialModulus::<ProptestF17Sqrt3Field>::new(vec![
            ProptestF17Sqrt3Field::one(),
            ProptestF17Sqrt3Field::one(),
            ProptestF17Sqrt3Field::zero(),
            ProptestF17Sqrt3Field::one(),
        ])
        .expect("tower modulus should be structurally valid")
    }

    fn check_field_conditions() -> Result<(), crate::fields::FieldError> {
        Ok(())
    }
}

pub type ProptestF17Sqrt3Field = ExtensionField<ProptestF17Sqrt3Spec>;
pub type ProptestF17TowerField = ExtensionField<ProptestF17TowerSpec>;
pub type ProptestF17Sqrt3Elem = ExtensionFieldElement<ProptestF17Sqrt3Spec>;
pub type ProptestF17TowerElem = ExtensionFieldElement<ProptestF17TowerSpec>;

/// Coupled fixture that exposes compatible base, quadratic, and tower values.
#[derive(Clone, Debug)]
pub struct TowerElementCase {
    pub base_left: crate::fields::Fp17Elem,
    pub base_right: crate::fields::Fp17Elem,
    pub quadratic_left: ProptestF17Sqrt3Elem,
    pub quadratic_right: ProptestF17Sqrt3Elem,
    pub tower_left: ProptestF17TowerElem,
    pub tower_right: ProptestF17TowerElem,
}

/// Returns a shrink-friendly coupled tower fixture.
pub fn arb_tower_element_case() -> BoxedStrategy<TowerElementCase> {
    (
        arb_fp_elem::<crate::fields::Fp17>(),
        arb_fp_elem::<crate::fields::Fp17>(),
        arb_semantic_extension_elem::<ProptestF17Sqrt3Spec>(),
        arb_semantic_extension_elem::<ProptestF17Sqrt3Spec>(),
        arb_semantic_extension_elem::<ProptestF17TowerSpec>(),
        arb_semantic_extension_elem::<ProptestF17TowerSpec>(),
    )
        .prop_map(
            |(base_left, base_right, quadratic_left, quadratic_right, tower_left, tower_right)| {
                TowerElementCase {
                    base_left,
                    base_right,
                    quadratic_left,
                    quadratic_right,
                    tower_left,
                    tower_right,
                }
            },
        )
        .boxed()
}

pub(crate) fn touch_family_case_fields() {
    let _ = |case: TowerElementCase| {
        let _ = case.base_left;
        let _ = case.base_right;
        let _ = case.quadratic_left;
        let _ = case.quadratic_right;
        let _ = case.tower_left;
        let _ = case.tower_right;
    };
}
