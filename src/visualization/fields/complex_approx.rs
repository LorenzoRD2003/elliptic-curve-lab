use num_complex::Complex64;

use crate::fields::{complex_approx::ComplexApprox, traits::Field};
use crate::visualization::{
    Visualizable, VisualizableField,
    shared::{is_small_complex, is_small_real},
};

fn is_negligible_component(component: f64, other_component: f64) -> bool {
    let scale = other_component.abs().max(1.0);
    component.abs() <= 1.0e-9 * scale
}

fn format_decimal_compact(value: f64) -> String {
    let mut text = format!("{value:.6}");

    if text.contains('.') {
        while text.ends_with('0') {
            text.pop();
        }
        if text.ends_with('.') {
            text.pop();
        }
    }

    if text == "-0" { "0".to_string() } else { text }
}

/// Formats a complex number as `a + bi` or `a - bi`.
pub(crate) fn format_complex(z: &Complex64) -> String {
    let imag_sign = if z.im < 0.0 { '-' } else { '+' };
    format!("{:.6} {} {:.6}i", z.re, imag_sign, z.im.abs())
}

/// Formats a complex number compactly, suppressing numerically negligible
/// real or imaginary parts while still showing exact zero as `0`.
pub(crate) fn format_complex_compact(z: &Complex64) -> String {
    if is_small_complex(z) {
        return "0".to_string();
    }

    if is_small_real(z.im) || is_negligible_component(z.im, z.re) {
        return format_decimal_compact(z.re);
    }

    if is_small_real(z.re) || is_negligible_component(z.re, z.im) {
        return match format_decimal_compact(z.im) {
            value if value == "1" => "i".to_string(),
            value if value == "-1" => "-i".to_string(),
            value => format!("{value}i"),
        };
    }

    let real = format_decimal_compact(z.re);
    let imag = match format_decimal_compact(z.im.abs()) {
        value if value == "1" => "i".to_string(),
        value => format!("{value}i"),
    };
    let imag_sign = if z.im < 0.0 { '-' } else { '+' };

    format!("{real} {imag_sign} {imag}")
}

/// Returns a short textual description of a complex number.
fn describe_complex(z: &Complex64) -> String {
    format!(
        "z = {}\n|z| = {:.6}\narg(z) = {:.6} rad\napprox zero: {}",
        format_complex(z),
        z.norm(),
        z.arg(),
        ComplexApprox::is_zero(z)
    )
}

impl Visualizable for Complex64 {
    fn format_compact(&self) -> String {
        format_complex_compact(self)
    }

    fn describe(&self) -> String {
        describe_complex(self)
    }
}

impl VisualizableField for Complex64 {
    fn format_elem(&self) -> String {
        format_complex_compact(self)
    }

    fn inverse(&self) -> Option<String> {
        ComplexApprox::inv(self).map(|value| format_complex(&value))
    }

    fn explain_add(lhs: &Self, rhs: &Self) -> Option<String> {
        let result = ComplexApprox::add(lhs, rhs);
        Some(format!(
            "Addition in C (approx)\n\
             lhs: {}\n\
             rhs: {}\n\
             result: {}",
            format_complex(lhs),
            format_complex(rhs),
            format_complex(&result)
        ))
    }

    fn explain_mul(lhs: &Self, rhs: &Self) -> Option<String> {
        let result = ComplexApprox::mul(lhs, rhs);
        Some(format!(
            "Multiplication in C (approx)\n\
             lhs: {}\n\
             rhs: {}\n\
             result: {}",
            format_complex(lhs),
            format_complex(rhs),
            format_complex(&result)
        ))
    }

    fn explain_div(lhs: &Self, rhs: &Self) -> Option<String> {
        let reciprocal = ComplexApprox::inv(rhs)?;
        let result = ComplexApprox::mul(lhs, &reciprocal);

        Some(format!(
            "Division in C (approx)\n\
             lhs: {}\n\
             rhs: {}\n\
             reciprocal of rhs: {}\n\
             reduction to multiplication: {} * {} = {}",
            format_complex(lhs),
            format_complex(rhs),
            format_complex(&reciprocal),
            format_complex(lhs),
            format_complex(&reciprocal),
            format_complex(&result)
        ))
    }
}

#[cfg(test)]
mod tests {
    use num_complex::Complex64;

    use super::*;
    use crate::visualization::{Visualizable, VisualizableField};

    #[test]
    fn complex_formatting_is_human_readable() {
        let z = Complex64::new(2.5, -3.0);
        assert_eq!(format_complex(&z), "2.500000 - 3.000000i");
    }

    #[test]
    fn complex_description_contains_key_quantities() {
        let z = Complex64::new(3.0, 4.0);
        let description = describe_complex(&z);
        assert!(description.contains("z = 3.000000 + 4.000000i"));
        assert!(description.contains("|z| = 5.000000"));
        assert!(description.contains("approx zero: false"));
    }

    #[test]
    fn complex_visualizable_trait_reuses_core_helpers() {
        let lhs = Complex64::new(1.0, 2.0);
        let rhs = Complex64::new(3.0, -1.0);

        assert_eq!(lhs.format_compact(), "1 + 2i");
        assert!(lhs.describe().contains("|z|"));
        assert_eq!(
            lhs.inverse().expect("complex inverse should exist"),
            "0.200000 - 0.400000i"
        );

        let add =
            Complex64::explain_add(&lhs, &rhs).expect("complex addition should be explainable");
        assert!(add.contains("Addition in C (approx)"));
        assert!(add.contains("4.000000 + 1.000000i"));

        let mul = Complex64::explain_mul(&lhs, &rhs)
            .expect("complex multiplication should be explainable");
        assert!(mul.contains("Multiplication in C (approx)"));

        let div =
            Complex64::explain_div(&lhs, &rhs).expect("complex division should be explainable");
        assert!(div.contains("Division in C (approx)"));
    }

    #[test]
    fn compact_complex_formatter_suppresses_negligible_real_or_imaginary_noise() {
        assert_eq!(
            format_complex_compact(&Complex64::new(-1.0e-15, -24.15094)),
            "-24.15094i"
        );
        assert_eq!(
            format_complex_compact(&Complex64::new(188.795905, 1.0e-15)),
            "188.795905"
        );
        assert_eq!(format_complex_compact(&Complex64::new(0.0, 0.0)), "0");
    }

    #[test]
    fn compact_complex_formatter_drops_pointless_integer_decimals_and_unit_imaginary_coefficients()
    {
        assert_eq!(format_complex_compact(&Complex64::new(1.0, 1.0)), "1 + i");
        assert_eq!(format_complex_compact(&Complex64::new(1.0, -1.0)), "1 - i");
        assert_eq!(format_complex_compact(&Complex64::new(0.0, 1.0)), "i");
        assert_eq!(format_complex_compact(&Complex64::new(0.0, -1.0)), "-i");
        assert_eq!(format_complex_compact(&Complex64::new(2.0, 0.0)), "2");
        assert_eq!(format_complex_compact(&Complex64::new(2.5, 0.0)), "2.5");
    }
}
