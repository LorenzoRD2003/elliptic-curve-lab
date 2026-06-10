//! Internal property-testing infrastructure for the repository.
//!
//! This module is intentionally organized by mathematical domain so generators
//! and reusable fixtures live near the concepts they exercise instead of
//! accumulating in one flat helper file.

pub mod combinators;
pub mod config;
pub mod elliptic_curves;
pub mod fields;
pub mod isogenies;
pub mod polynomials;

#[cfg(test)]
mod tests;

fn touch_inventory() {
    combinators::touch_combinator_inventory();
    config::touch_config_inventory();
    fields::touch_field_inventory();
    polynomials::touch_polynomial_inventory();
    elliptic_curves::touch_elliptic_curve_inventory();
    isogenies::touch_isogeny_inventory();
}

const _: fn() = touch_inventory;
