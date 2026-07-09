use num_bigint::BigUint;
use num_complex::Complex64;

use crate::elliptic_curves::endomorphisms::quadratic_orders::ImaginaryQuadraticOrder;
use crate::fields::traits::Field;
use crate::polynomials::irreducibility::ReducibilityReason;
use crate::visualization::{Visualizable, VisualizableField};

pub(crate) fn format_field_elem<F: Field>(value: &F::Elem) -> String
where
    F::Elem: VisualizableField,
{
    value.format_elem()
}

pub(crate) fn parenthesize_if_needed(text: &str) -> String {
    if text.contains(' ') || text.starts_with('-') || text.contains('/') {
        format!("({text})")
    } else {
        text.to_string()
    }
}

pub(crate) fn yes_no(value: bool) -> &'static str {
    if value { "yes" } else { "no" }
}

pub(crate) fn comma_list<I>(items: I) -> String
where
    I: IntoIterator<Item = String>,
{
    items.into_iter().collect::<Vec<_>>().join(", ")
}

pub(crate) fn compact_visualizable_list<'a, I, T>(items: I) -> String
where
    I: IntoIterator<Item = &'a T>,
    T: Visualizable + 'a,
{
    comma_list(items.into_iter().map(Visualizable::format_compact))
}

pub(crate) fn is_small_real(value: f64) -> bool {
    value.abs() <= 1.0e-12
}

pub(crate) fn is_small_complex(value: &Complex64) -> bool {
    value.norm() <= 1.0e-12
}

pub(crate) fn format_reducibility_reason(reason: ReducibilityReason) -> &'static str {
    match reason {
        ReducibilityReason::AlgebraicallyClosed => {
            "the base field is algebraically closed, so every degree >= 2 polynomial factors non-trivially"
        }
    }
}

pub(crate) fn format_order_conductor_label(conductor: &BigUint) -> String {
    format!("O_{conductor}")
}

pub(crate) fn format_imaginary_quadratic_order_label(order: &ImaginaryQuadraticOrder) -> String {
    format!(
        "O_{} (Δ = {})",
        order.conductor(),
        order.discriminant().value()
    )
}
