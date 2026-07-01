//! Exact numerics-oriented property-test strategies.

mod square_roots;

pub(crate) use square_roots::{
    ModularSquareRootCase, arb_modular_square_root_case, brute_force_square_roots_mod_u64,
};

pub(crate) fn touch_numerics_inventory() {
    let _ = arb_modular_square_root_case(
        crate::proptest_support::config::NumericsStrategyConfig::default(),
    );
    let _ = brute_force_square_roots_mod_u64(1, 8);
    let _ = core::mem::size_of::<ModularSquareRootCase>();
}
