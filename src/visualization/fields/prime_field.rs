use crate::visualization::*;
use num_bigint::{BigInt, BigUint};

use crate::fields::{
    FieldError, Fp, Fp2, Fp2Elem, FpElem,
    traits::{EnumerableFiniteField, FiniteField},
};
use crate::visualization::VisualizableField;
use crate::visualization::traits::Visualizable;
use crypto_bigint::modular::ConstPrimeMontyParams;

/// Returns a short textual description of a finite prime field.
pub fn format_prime_field<F>() -> Result<String, FieldError>
where
    F: FiniteField,
{
    F::check_structure()?;
    let characteristic = F::characteristic()
        .to_positive_biguint()
        .expect("finite fields have positive characteristic");
    Ok(format!(
        "GF({characteristic})\ncharacteristic: {characteristic}\nextension degree: {}\ncardinality: {}",
        F::extension_degree(),
        F::cardinality_biguint()
    ))
}

/// Formats a field element using its compact canonical representation.
pub fn format_fp_elem<E>(elem: &E) -> String
where
    E: VisualizableField,
{
    elem.format_elem()
}

impl<M, const LIMBS: usize> Visualizable for FpElem<M, LIMBS>
where
    M: ConstPrimeMontyParams<LIMBS>,
{
    fn format_compact(&self) -> String {
        self.to_biguint().to_string()
    }

    fn describe(&self) -> String {
        format!(
            "element: {}\nrepresentative: {}\nfield: GF({})",
            self.format_compact(),
            self.to_biguint(),
            Fp::<M, LIMBS>::modulus_biguint()
        )
    }
}

impl<M, const LIMBS: usize> VisualizableField for FpElem<M, LIMBS>
where
    M: ConstPrimeMontyParams<LIMBS>,
{
    fn format_elem(&self) -> String {
        self.format_compact()
    }

    fn inverse(&self) -> Option<String> {
        Fp::<M, LIMBS>::inv(self).map(|value| value.format_compact())
    }

    fn explain_add(lhs: &Self, rhs: &Self) -> Option<String> {
        let modulus = Fp::<M, LIMBS>::modulus_biguint();
        let raw_sum = lhs.to_biguint() + rhs.to_biguint();
        let reduced = Fp::<M, LIMBS>::add(lhs, rhs);

        Some(format!(
            "Addition in GF({modulus})\n\
             canonical lhs: {}\n\
             canonical rhs: {}\n\
             raw sum: {} + {} = {raw_sum}\n\
             reduction: {raw_sum} mod {modulus} = {}\n\
             result: {}",
            lhs.format_compact(),
            rhs.format_compact(),
            lhs.to_biguint(),
            rhs.to_biguint(),
            reduced.to_biguint(),
            reduced.format_compact()
        ))
    }

    fn explain_mul(lhs: &Self, rhs: &Self) -> Option<String> {
        let modulus = Fp::<M, LIMBS>::modulus_biguint();
        let raw_product = lhs.to_biguint() * rhs.to_biguint();
        let reduced = Fp::<M, LIMBS>::mul(lhs, rhs);

        Some(format!(
            "Multiplication in GF({modulus})\n\
             canonical lhs: {}\n\
             canonical rhs: {}\n\
             raw product: {} * {} = {raw_product}\n\
             reduction: {raw_product} mod {modulus} = {}\n\
             result: {}",
            lhs.format_compact(),
            rhs.format_compact(),
            lhs.to_biguint(),
            rhs.to_biguint(),
            reduced.to_biguint(),
            reduced.format_compact()
        ))
    }

    fn explain_div(lhs: &Self, rhs: &Self) -> Option<String> {
        let modulus = Fp::<M, LIMBS>::modulus_biguint();
        let reciprocal = Fp::<M, LIMBS>::inv(rhs)?;
        let result = Fp::<M, LIMBS>::mul(lhs, &reciprocal);

        Some(format!(
            "Division in GF({modulus})\n\
             lhs: {}\n\
             rhs: {}\n\
             inverse of rhs: {}\n\
             reduction to multiplication: {} * {} mod {modulus} = {}",
            lhs.format_compact(),
            rhs.format_compact(),
            reciprocal.format_compact(),
            lhs.to_biguint(),
            reciprocal.to_biguint(),
            result.to_biguint()
        ))
    }
}

impl Visualizable for Fp2Elem {
    fn format_compact(&self) -> String {
        fp2_value(self).to_string()
    }

    fn describe(&self) -> String {
        format!(
            "element: {}\nrepresentative: {}\nfield: GF(2)",
            self.format_compact(),
            fp2_value(self)
        )
    }
}

impl VisualizableField for Fp2Elem {
    fn format_elem(&self) -> String {
        self.format_compact()
    }

    fn inverse(&self) -> Option<String> {
        Fp2::inv(self).map(|value| value.format_compact())
    }

    fn explain_add(lhs: &Self, rhs: &Self) -> Option<String> {
        let reduced = Fp2::add(lhs, rhs);
        Some(format!(
            "Addition in GF(2)\n\
             canonical lhs: {}\n\
             canonical rhs: {}\n\
             reduction: {} + {} mod 2 = {}\n\
             result: {}",
            lhs.format_compact(),
            rhs.format_compact(),
            fp2_value(lhs),
            fp2_value(rhs),
            fp2_value(&reduced),
            reduced.format_compact()
        ))
    }

    fn explain_mul(lhs: &Self, rhs: &Self) -> Option<String> {
        let reduced = Fp2::mul(lhs, rhs);
        Some(format!(
            "Multiplication in GF(2)\n\
             canonical lhs: {}\n\
             canonical rhs: {}\n\
             reduction: {} * {} mod 2 = {}\n\
             result: {}",
            lhs.format_compact(),
            rhs.format_compact(),
            fp2_value(lhs),
            fp2_value(rhs),
            fp2_value(&reduced),
            reduced.format_compact()
        ))
    }
}

fn fp2_value(element: &Fp2Elem) -> u8 {
    element.value()
}

/// Explains a modular addition step by step in a finite field.
pub fn explain_add<F>(lhs: &BigUint, rhs: &BigUint) -> Result<String, FieldError>
where
    F: FiniteField,
    F::Elem: VisualizableField,
{
    let left = elem_from_biguint::<F>(lhs)?;
    let right = elem_from_biguint::<F>(rhs)?;
    Ok(F::Elem::explain_add(&left, &right).expect("field elements explain addition"))
}

/// Explains a modular multiplication step by step in a finite field.
pub fn explain_mul<F>(lhs: &BigUint, rhs: &BigUint) -> Result<String, FieldError>
where
    F: FiniteField,
    F::Elem: VisualizableField,
{
    let left = elem_from_biguint::<F>(lhs)?;
    let right = elem_from_biguint::<F>(rhs)?;
    Ok(F::Elem::explain_mul(&left, &right).expect("field elements explain multiplication"))
}

/// Explains how the multiplicative inverse of an element behaves in a finite field.
pub fn explain_inverse<F>(value: &BigUint) -> Result<String, FieldError>
where
    F: FiniteField,
    F::Elem: VisualizableField,
{
    let element = elem_from_biguint::<F>(value)?;

    if F::is_zero(&element) {
        return Err(FieldError::DivisionByZero);
    }

    let inverse = F::inverse(&element)?;
    let verification = F::mul(&element, &inverse);

    Ok(format!(
        "Inverse in GF({})\n\
         element: {}\n\
         inverse: {}\n\
         verification: {} * {} = {}",
        F::characteristic()
            .to_positive_biguint()
            .expect("finite fields have positive characteristic"),
        element.format_elem(),
        inverse.format_elem(),
        element.format_elem(),
        inverse.format_elem(),
        verification.format_elem()
    ))
}

fn elem_from_biguint<F>(value: &BigUint) -> Result<F::Elem, FieldError>
where
    F: FiniteField,
{
    F::check_structure()?;
    Ok(F::from_bigint(&BigInt::from(value.clone())))
}

/// Builds the full addition table for a small enumerable finite field.
pub fn addition_table<F>() -> Result<String, FieldError>
where
    F: FiniteField + EnumerableFiniteField,
    F::Elem: VisualizableField,
{
    render_binary_operation_table::<F>("Addition table", F::add)
}

/// Builds the full multiplication table for a small enumerable finite field.
pub fn multiplication_table<F>() -> Result<String, FieldError>
where
    F: FiniteField + EnumerableFiniteField,
    F::Elem: VisualizableField,
{
    render_binary_operation_table::<F>("Multiplication table", F::mul)
}

/// Builds a table of multiplicative inverses for the non-zero elements of a field.
pub fn inverses_table<F>() -> Result<String, FieldError>
where
    F: FiniteField + EnumerableFiniteField,
    F::Elem: VisualizableField,
{
    F::check_structure()?;
    let mut lines = vec![
        format!(
            "Inverse table for GF({})",
            F::characteristic().to_positive_biguint().unwrap()
        ),
        "a | a^-1 | check".to_string(),
        "----------------".to_string(),
    ];

    for element in F::elements()
        .into_iter()
        .filter(|element| !F::is_zero(element))
    {
        let inverse = F::inverse(&element)?;
        let check = F::mul(&element, &inverse);
        lines.push(format!(
            "{} | {} | {} * {} = {}",
            element.format_elem(),
            inverse.format_elem(),
            element.format_elem(),
            inverse.format_elem(),
            check.format_elem()
        ));
    }

    Ok(lines.join("\n"))
}

fn render_binary_operation_table<F>(
    title: &str,
    operation: fn(&F::Elem, &F::Elem) -> F::Elem,
) -> Result<String, FieldError>
where
    F: FiniteField + EnumerableFiniteField,
    F::Elem: VisualizableField,
{
    F::check_structure()?;
    let elements = F::elements();
    let width = elements
        .iter()
        .map(|element| element.format_elem().len())
        .max()
        .unwrap_or(1);
    let mut lines = vec![format!(
        "{title} for GF({})",
        F::characteristic().to_positive_biguint().unwrap()
    )];

    let mut header = format!("{:>width$} |", "", width = width);
    for element in &elements {
        header.push_str(&format!(
            " {:>width$}",
            element.format_elem(),
            width = width
        ));
    }
    lines.push(header);
    lines.push("-".repeat(lines.last().map_or(0, String::len)));

    for row in &elements {
        let mut line = format!("{:>width$} |", row.format_elem(), width = width);
        for col in &elements {
            let value = operation(row, col);
            line.push_str(&format!(" {:>width$}", value.format_elem(), width = width));
        }
        lines.push(line);
    }

    Ok(lines.join("\n"))
}
