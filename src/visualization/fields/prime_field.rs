use crate::fields::{
    errors::FieldError,
    prime_field::{Fp, FpElem},
    traits::Field,
};
use crate::visualization::{Visualizable, VisualizableField};

/// Returns a short textual description of the prime field `GF(P)`.
pub fn format_prime_field<const P: u64>() -> Result<String, FieldError> {
    Fp::<P>::validate_modulus()?;
    Ok(format!(
        "GF({P})\ncharacteristic: {P}\nextension degree: 1\ncardinality: {P}"
    ))
}

/// Formats a prime-field element using its canonical representative.
pub fn format_fp_elem<const P: u64>(elem: &FpElem<P>) -> String {
    format!("{} (mod {P})", elem.value())
}

impl<const P: u64> Visualizable for FpElem<P> {
    fn format_compact(&self) -> String {
        format_fp_elem(self)
    }

    fn describe(&self) -> String {
        format!(
            "element: {}\nrepresentative: {}\nfield: GF({P})",
            format_fp_elem(self),
            self.value()
        )
    }
}

impl<const P: u64> VisualizableField for FpElem<P> {
    fn format_elem(&self) -> String {
        format_fp_elem(self)
    }

    fn inverse(&self) -> Option<String> {
        Fp::<P>::inv(self).map(|value| format_fp_elem(&value))
    }

    fn explain_add(lhs: &Self, rhs: &Self) -> Option<String> {
        let raw_sum = u128::from(lhs.value()) + u128::from(rhs.value());
        let reduced = Fp::<P>::add(lhs, rhs);

        Some(format!(
            "Addition in GF({P})\n\
             canonical lhs: {}\n\
             canonical rhs: {}\n\
             raw sum: {} + {} = {raw_sum}\n\
             reduction: {raw_sum} mod {P} = {}\n\
             result: {}",
            format_fp_elem(lhs),
            format_fp_elem(rhs),
            lhs.value(),
            rhs.value(),
            reduced.value(),
            format_fp_elem(&reduced)
        ))
    }

    fn explain_mul(lhs: &Self, rhs: &Self) -> Option<String> {
        let raw_product = u128::from(lhs.value()) * u128::from(rhs.value());
        let reduced = Fp::<P>::mul(lhs, rhs);

        Some(format!(
            "Multiplication in GF({P})\n\
             canonical lhs: {}\n\
             canonical rhs: {}\n\
             raw product: {} * {} = {raw_product}\n\
             reduction: {raw_product} mod {P} = {}\n\
             result: {}",
            format_fp_elem(lhs),
            format_fp_elem(rhs),
            lhs.value(),
            rhs.value(),
            reduced.value(),
            format_fp_elem(&reduced)
        ))
    }

    fn explain_div(lhs: &Self, rhs: &Self) -> Option<String> {
        let reciprocal = Fp::<P>::inv(rhs)?;
        let result = Fp::<P>::mul(lhs, &reciprocal);

        Some(format!(
            "Division in GF({P})\n\
             lhs: {}\n\
             rhs: {}\n\
             inverse of rhs: {}\n\
             reduction to multiplication: {} * {} mod {P} = {}",
            format_fp_elem(lhs),
            format_fp_elem(rhs),
            format_fp_elem(&reciprocal),
            lhs.value(),
            reciprocal.value(),
            result.value()
        ))
    }
}

/// Explains a modular addition step by step in `GF(P)`.
pub fn explain_add<const P: u64>(lhs: u64, rhs: u64) -> Result<String, FieldError> {
    Fp::<P>::validate_modulus()?;
    let left = Fp::<P>::new_elem(lhs)?;
    let right = Fp::<P>::new_elem(rhs)?;
    Ok(FpElem::<P>::explain_add(&left, &right).expect("prime-field elements explain addition"))
}

/// Explains a modular multiplication step by step in `GF(P)`.
pub fn explain_mul<const P: u64>(lhs: u64, rhs: u64) -> Result<String, FieldError> {
    Fp::<P>::validate_modulus()?;
    let left = Fp::<P>::new_elem(lhs)?;
    let right = Fp::<P>::new_elem(rhs)?;
    Ok(FpElem::<P>::explain_mul(&left, &right)
        .expect("prime-field elements explain multiplication"))
}

/// Explains how the multiplicative inverse of an element behaves in `GF(P)`.
pub fn explain_inverse<const P: u64>(value: u64) -> Result<String, FieldError> {
    Fp::<P>::validate_modulus()?;
    let element = Fp::<P>::new_elem(value)?;

    if Fp::<P>::is_zero(&element) {
        return Err(FieldError::DivisionByZero);
    }

    let inverse = Fp::<P>::inverse(&element)?;
    let verification = Fp::<P>::mul(&element, &inverse);

    Ok(format!(
        "Inverse in GF({P})\n\
         element: {}\n\
         inverse: {}\n\
         verification: {} * {} mod {P} = {}",
        format_fp_elem(&element),
        format_fp_elem(&inverse),
        element.value(),
        inverse.value(),
        verification.value()
    ))
}

/// Builds the full addition table for `GF(P)` as aligned plain text.
pub fn addition_table<const P: u64>() -> Result<String, FieldError> {
    Fp::<P>::validate_modulus()?;
    render_binary_operation_table::<P>("Addition table", Fp::<P>::add)
}

/// Builds the full multiplication table for `GF(P)` as aligned plain text.
pub fn multiplication_table<const P: u64>() -> Result<String, FieldError> {
    Fp::<P>::validate_modulus()?;
    render_binary_operation_table::<P>("Multiplication table", Fp::<P>::mul)
}

/// Builds a table of multiplicative inverses for the non-zero elements of `GF(P)`.
pub fn inverses_table<const P: u64>() -> Result<String, FieldError> {
    Fp::<P>::validate_modulus()?;

    let mut lines = vec![
        format!("Inverse table for GF({P})"),
        "a | a^-1 | check".to_string(),
        "----------------".to_string(),
    ];

    for value in 1..P {
        let element = Fp::<P>::new_elem(value)?;
        let inverse = Fp::<P>::inverse(&element)?;
        let check = Fp::<P>::mul(&element, &inverse);
        lines.push(format!(
            "{:>2} | {:>4} | {} * {} = {}",
            element.value(),
            inverse.value(),
            element.value(),
            inverse.value(),
            check.value()
        ));
    }

    Ok(lines.join("\n"))
}

/// Renders a binary operation table for `GF(P)`.
fn render_binary_operation_table<const P: u64>(
    title: &str,
    operation: fn(&FpElem<P>, &FpElem<P>) -> FpElem<P>,
) -> Result<String, FieldError> {
    let width = cell_width(P);
    let mut lines = vec![format!("{title} for GF({P})")];

    let header = table_header::<P>(width);
    lines.push(header);
    lines.push("-".repeat(lines.last().map_or(0, String::len)));

    for row in 0..P {
        let row_elem = Fp::<P>::new_elem(row)?;
        let mut line = format!("{:>width$} |", row_elem.value(), width = width);

        for col in 0..P {
            let col_elem = Fp::<P>::new_elem(col)?;
            let value = operation(&row_elem, &col_elem);
            line.push_str(&format!(" {:>width$}", value.value(), width = width));
        }

        lines.push(line);
    }

    Ok(lines.join("\n"))
}

/// Builds the common header row for operation tables.
fn table_header<const P: u64>(width: usize) -> String {
    let mut header = format!("{:>width$} |", "", width = width);
    for value in 0..P {
        header.push_str(&format!(" {:>width$}", value, width = width));
    }
    header
}

/// Computes a stable column width for plain-text tables.
fn cell_width(modulus: u64) -> usize {
    modulus.saturating_sub(1).to_string().len().max(1)
}

#[cfg(test)]
mod tests {
    use super::{
        addition_table, explain_add, explain_inverse, explain_mul, format_fp_elem,
        format_prime_field, inverses_table, multiplication_table,
    };
    use crate::fields::{FieldError, Fp, FpElem};
    use crate::visualization::{Visualizable, VisualizableField};

    type F17 = Fp<17>;
    type E17 = FpElem<17>;

    #[test]
    fn prime_field_summary_is_readable() {
        let summary = format_prime_field::<17>().expect("field should be valid");
        assert!(summary.contains("GF(17)"));
        assert!(summary.contains("characteristic: 17"));
        assert!(summary.contains("cardinality: 17"));
    }

    #[test]
    fn prime_field_element_format_is_compact() {
        let element = E17::new(20).expect("element should be created");
        assert_eq!(format_fp_elem(&element), "3 (mod 17)");
    }

    #[test]
    fn addition_explanation_shows_reduction() {
        let explanation = explain_add::<17>(13, 9).expect("explanation should succeed");
        assert!(explanation.contains("raw sum: 13 + 9 = 22"));
        assert!(explanation.contains("reduction: 22 mod 17 = 5"));
        assert!(explanation.contains("result: 5 (mod 17)"));
    }

    #[test]
    fn multiplication_explanation_shows_reduction() {
        let explanation = explain_mul::<17>(5, 7).expect("explanation should succeed");
        assert!(explanation.contains("raw product: 5 * 7 = 35"));
        assert!(explanation.contains("reduction: 35 mod 17 = 1"));
        assert!(explanation.contains("result: 1 (mod 17)"));
    }

    #[test]
    fn inverse_explanation_shows_verification() {
        let explanation = explain_inverse::<17>(3).expect("inverse should exist");
        assert!(explanation.contains("inverse: 6 (mod 17)"));
        assert!(explanation.contains("verification: 3 * 6 mod 17 = 1"));
    }

    #[test]
    fn inverse_explanation_rejects_zero() {
        let error = explain_inverse::<17>(0).expect_err("zero should not be invertible");
        assert!(matches!(error, FieldError::DivisionByZero));
    }

    #[test]
    fn addition_table_contains_expected_entries() {
        let table = addition_table::<17>().expect("table should render");
        assert!(table.contains("Addition table for GF(17)"));
        assert!(table.contains("13 |"));
        assert!(table.contains(" 5"));
    }

    #[test]
    fn multiplication_table_contains_expected_entries() {
        let table = multiplication_table::<17>().expect("table should render");
        assert!(table.contains("Multiplication table for GF(17)"));
        assert!(table.contains(" 5 |"));
        assert!(table.contains(" 1"));
    }

    #[test]
    fn inverse_table_contains_expected_inverse_pair() {
        let table = inverses_table::<17>().expect("table should render");
        assert!(table.contains("Inverse table for GF(17)"));
        assert!(table.contains(" 3 |    6 | 3 * 6 = 1"));
    }

    #[test]
    fn invalid_prime_field_is_rejected_by_visualizers() {
        let error = format_prime_field::<15>().expect_err("GF(15) should be rejected");
        assert!(matches!(error, FieldError::InvalidModulus { modulus: 15 }));
    }

    #[test]
    fn formatting_matches_runtime_values() {
        let element = F17::new_elem(34).expect("element should be created");
        assert_eq!(format_fp_elem(&element), "0 (mod 17)");
    }

    #[test]
    fn prime_field_visualizable_trait_reuses_core_helpers() {
        let lhs = E17::new(13).expect("lhs should exist");
        let rhs = E17::new(9).expect("rhs should exist");

        assert_eq!(lhs.format_compact(), "13 (mod 17)");
        assert_eq!(lhs.format_elem(), "13 (mod 17)");
        assert!(lhs.describe().contains("field: GF(17)"));
        assert_eq!(rhs.inverse().expect("inverse should exist"), "2 (mod 17)");

        let add = E17::explain_add(&lhs, &rhs).expect("prime-field addition should be explainable");
        assert!(add.contains("Addition in GF(17)"));
        assert!(add.contains("result: 5 (mod 17)"));

        let mul =
            E17::explain_mul(&lhs, &rhs).expect("prime-field multiplication should be explainable");
        assert!(mul.contains("Multiplication in GF(17)"));

        let div = E17::explain_div(&lhs, &rhs).expect("prime-field division should be explainable");
        assert!(div.contains("Division in GF(17)"));
    }
}
