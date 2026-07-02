//! Isogeny-oriented fixtures and reusable tiny finite-field cases.

pub mod composition;
pub mod function_field_maps;
pub mod isomorphisms;
pub mod kernels;
pub mod velu;

pub use composition::{ComposableVeluCase, arb_composable_velu_case};
pub use function_field_maps::{
    ComposableFunctionFieldMapCase, FunctionFieldMapCase,
    arb_composable_short_weierstrass_function_field_map_case,
    arb_short_weierstrass_function_field_map_case,
};
pub use isomorphisms::{ShortWeierstrassIsomorphismCase, arb_short_weierstrass_isomorphism_case};
pub use velu::{CyclicKernelCase, arb_cyclic_kernel_case};

pub(crate) fn touch_isogeny_inventory() {
    let config = crate::proptest_support::config::CurveStrategyConfig::default();
    let _ = arb_cyclic_kernel_case();
    let _ = arb_composable_velu_case();
    let _ = arb_short_weierstrass_function_field_map_case::<crate::fields::Fp17>(config);
    let _ = arb_composable_short_weierstrass_function_field_map_case::<crate::fields::Fp17>(config);
    let _ = arb_short_weierstrass_isomorphism_case::<crate::fields::Fp17>(config);
    let _ = core::mem::size_of::<CyclicKernelCase>();
    let _ = core::mem::size_of::<ComposableVeluCase>();
    let _ = core::mem::size_of::<FunctionFieldMapCase<crate::fields::Fp17>>();
    let _ = core::mem::size_of::<ComposableFunctionFieldMapCase<crate::fields::Fp17>>();
    let _ = core::mem::size_of::<ShortWeierstrassIsomorphismCase<crate::fields::Fp17>>();
    velu::touch_cyclic_kernel_case_fields();
    function_field_maps::touch_function_field_map_case_fields();
    isomorphisms::touch_isomorphism_case_fields();
}
