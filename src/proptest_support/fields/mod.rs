//! Field-oriented strategies and reusable canonical field fixtures.

use crate::fields::traits::*;
pub mod complex;
pub mod extension;
pub mod families;
pub mod prime;
pub mod rational;
pub mod rational_function;

pub use complex::arb_complex_approx;
pub use extension::{arb_extension_elem, arb_semantic_extension_elem};
pub use families::{
    ProptestF17Sqrt3Elem, ProptestF17Sqrt3Field, ProptestF17Sqrt3Spec, ProptestF17TowerElem,
    ProptestF17TowerField, ProptestF17TowerSpec, TowerElementCase, arb_tower_element_case,
};
pub use prime::{arb_distinct_fp_elems, arb_fp_elem, arb_nonzero_fp_elem};
pub use rational::arb_q_elem;
pub use rational_function::arb_rational_function;

pub(crate) fn touch_field_inventory() {
    let _ = arb_complex_approx(crate::proptest_support::config::FieldStrategyConfig::default());
    let _ = arb_extension_elem::<ProptestF17Sqrt3Spec>();
    let _ = arb_semantic_extension_elem::<ProptestF17TowerSpec>();
    let _ = arb_fp_elem::<crate::fields::Fp17>();
    let _ = arb_nonzero_fp_elem::<crate::fields::Fp17>();
    let _ = arb_distinct_fp_elems::<crate::fields::Fp17>(2);
    let _ = arb_q_elem(crate::proptest_support::config::FieldStrategyConfig::default());
    let _ = arb_rational_function::<crate::fields::Fp17>(
        crate::proptest_support::config::PolynomialStrategyConfig::default(),
    );
    let _ = core::mem::size_of::<ProptestF17Sqrt3Field>();
    let _ = core::mem::size_of::<ProptestF17TowerField>();
    let _ = core::mem::size_of::<ProptestF17Sqrt3Elem>();
    let _ = core::mem::size_of::<ProptestF17TowerElem>();
    let _ = core::mem::size_of::<TowerElementCase>();
    let _ = arb_tower_element_case();
    families::touch_family_case_fields();
}
