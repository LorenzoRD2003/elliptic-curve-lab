use num_complex::Complex64;

use crate::fields::{
    FieldError,
    complex_approx::ComplexApprox,
    polynomial_field::{PolynomialFieldElement, PolynomialModulus},
    traits::Field,
};
use crate::polynomials::{
    DensePolynomial, PolynomialError,
    irreducibility::{IrreducibilityBackend, IrreducibilityStatus},
};
use crate::visualization::{
    Visualizable, VisualizableField,
    fields::complex_approx::format_complex,
    shared::{comma_list, format_reducibility_reason, yes_no},
};

/// Formats a polynomial over `GF(P)` from coefficients stored in ascending degree order.
fn format_prime_polynomial<F: Field>(coefficients: &[F::Elem]) -> String
where
    F::Elem: VisualizableField,
{
    let mut terms = Vec::new();

    for (power, coefficient) in coefficients.iter().enumerate().rev() {
        if F::is_zero(coefficient) {
            continue;
        }

        let value = coefficient.format_elem();
        let term = match power {
            0 => value.to_string(),
            1 if value == "1" => "x".to_string(),
            1 => format!("{value}*x"),
            _ if value == "1" => format!("x^{power}"),
            _ => format!("{value}*x^{power}"),
        };
        terms.push(term);
    }

    if terms.is_empty() {
        "0".to_string()
    } else {
        terms.join(" + ")
    }
}

/// Formats a polynomial over the approximate complex backend.
///
/// Coefficients are interpreted in ascending degree order, exactly like the
/// dense storage used elsewhere in the crate.
fn format_complex_polynomial(coefficients: &[Complex64]) -> String {
    let mut terms = Vec::new();

    for (power, coefficient) in coefficients.iter().enumerate().rev() {
        if ComplexApprox::is_zero(coefficient) {
            continue;
        }

        let coeff = format_complex(coefficient);
        let term = match power {
            0 => coeff,
            1 => format!("({coeff})*x"),
            _ => format!("({coeff})*x^{power}"),
        };
        terms.push(term);
    }

    if terms.is_empty() {
        "0".to_string()
    } else {
        terms.join(" + ")
    }
}

fn format_coefficients<F: Field>(coefficients: &[F::Elem]) -> String
where
    F::Elem: VisualizableField,
{
    format!(
        "[{}]",
        comma_list(coefficients.iter().map(VisualizableField::format_elem))
    )
}

/// Explains how the coefficient vector maps to a polynomial over `GF(P)`.
fn explain_prime_polynomial_storage<F: Field>(coefficients: &[F::Elem]) -> String
where
    F::Elem: VisualizableField,
{
    let mut lines = vec![
        format!("Polynomial over GF({})", F::characteristic()),
        "coefficients are stored in ascending degree order".to_string(),
        format!("polynomial: {}", format_prime_polynomial::<F>(coefficients)),
        "storage mapping:".to_string(),
    ];

    if coefficients.is_empty() {
        lines.push("- empty vector represents the zero polynomial".to_string());
        return lines.join("\n");
    }

    for (power, coefficient) in coefficients.iter().enumerate() {
        lines.push(format!(
            "- index {power}: coefficient {} multiplies x^{power}",
            coefficient.format_elem()
        ));
    }

    lines.join("\n")
}

/// Formats a modulus polynomial used in a quotient construction over `GF(P)`.
fn format_prime_polynomial_modulus<F>(modulus: &PolynomialModulus<F>) -> String
where
    F: Field,
    F::Elem: VisualizableField,
{
    format!(
        "m(x) = {}",
        format_prime_polynomial::<F>(modulus.coefficients())
    )
}

/// Returns a short textual description of a modulus polynomial over `GF(P)`.
fn describe_prime_polynomial_modulus<F>(modulus: &PolynomialModulus<F>) -> String
where
    F: Field,
    F::Elem: VisualizableField,
{
    format!(
        "Polynomial modulus over GF({})\n\
         degree: {}\n\
         raw coefficients (ascending): {}\n\
         expression: {}",
        F::characteristic(),
        modulus.degree(),
        format_coefficients::<F>(modulus.coefficients()),
        format_prime_polynomial_modulus::<F>(modulus)
    )
}

/// Describes whether a modulus polynomial over `GF(P)` is suitable as the
/// defining polynomial of a quotient field.
///
/// This helper bridges the `fields` and `polynomials` educational layers:
///
/// - it starts from the field-oriented notion of a modulus polynomial
/// - it reuses the polynomial irreducibility checker currently implemented for
///   dense polynomials over prime fields
/// - it explains the practical consequence for the quotient construction
///
/// The underlying irreducibility check currently uses the library's
/// exhaustive educational baseline algorithm.
fn describe_prime_polynomial_modulus_as_field_modulus<F: Field + IrreducibilityBackend>(
    modulus: &PolynomialModulus<F>,
) -> Result<String, PolynomialError>
where
    F::Elem: VisualizableField,
{
    let dense_modulus = DensePolynomial::<F>::new(modulus.coefficients().to_vec());
    let status = dense_modulus.irreducibility_status()?;

    let suitability = match &status {
        IrreducibilityStatus::Irreducible => {
            "suitable for a quotient field: yes, this modulus is irreducible over the base field"
        }
        IrreducibilityStatus::Linear => {
            "suitable for a quotient field: technically yes, but the quotient is just the base field again"
        }
        IrreducibilityStatus::Reducible { .. }
        | IrreducibilityStatus::ReducibleWithoutWitness { .. } => {
            "suitable for a quotient field: no, a reducible modulus gives only a quotient algebra in general"
        }
        IrreducibilityStatus::Constant => {
            "suitable for a quotient field: no, a constant polynomial cannot define a meaningful field extension"
        }
    };

    Ok(format!(
        "Field-modulus check over GF({})\n\
         expression: {}\n\
         base field algebraically closed: {}\n\
         irreducibility status: {}\n\
         {}",
        F::characteristic(),
        format_prime_polynomial_modulus::<F>(modulus),
        yes_no(F::IS_ALGEBRAICALLY_CLOSED),
        format_irreducibility_status::<F>(&status),
        suitability
    ))
}

/// Explains the irreducibility result for a modulus polynomial over `GF(P)`.
///
/// This helper is aimed at the field-construction use case: it explains not
/// just whether the polynomial is reducible, but what that means for the
/// quotient `GF(P)[x] / (m(x))`.
fn explain_prime_polynomial_modulus_irreducibility<F: IrreducibilityBackend>(
    modulus: &PolynomialModulus<F>,
) -> Result<String, PolynomialError>
where
    F::Elem: VisualizableField,
{
    let dense_modulus = DensePolynomial::<F>::new(modulus.coefficients().to_vec());
    let status = dense_modulus.irreducibility_status()?;

    let mut lines = vec![
        format!(
            "Irreducibility check for a field modulus over GF({})",
            F::characteristic()
        ),
        format!("modulus: {}", format_prime_polynomial_modulus::<F>(modulus)),
        format!(
            "base field algebraically closed: {}",
            yes_no(F::IS_ALGEBRAICALLY_CLOSED)
        ),
    ];

    if !F::IS_ALGEBRAICALLY_CLOSED {
        lines.push(
            "note: the base field is not algebraically closed, so higher-degree irreducible polynomials may exist".to_string(),
        );
    }

    match status {
        IrreducibilityStatus::Constant => {
            lines.push("status: constant".to_string());
            lines.push(
                "consequence: a constant polynomial is not a valid field-extension modulus"
                    .to_string(),
            );
        }
        IrreducibilityStatus::Linear => {
            lines.push("status: linear".to_string());
            lines.push(
                "consequence: a linear modulus is irreducible, but the quotient does not create a genuine new extension"
                    .to_string(),
            );
        }
        IrreducibilityStatus::Irreducible => {
            lines.push("status: irreducible".to_string());
            lines.push(
                "consequence: this modulus is suitable for a quotient-field construction over the base field".to_string(),
            );
        }
        IrreducibilityStatus::Reducible { divisor, quotient } => {
            lines.push("status: reducible".to_string());
            lines.push(format!(
                "witness divisor: {}",
                format_prime_polynomial::<F>(divisor.coefficients())
            ));
            lines.push(format!(
                "witness quotient: {}",
                format_prime_polynomial::<F>(quotient.coefficients())
            ));
            lines.push(format!(
                "factorization: {} = ({}) * ({})",
                format_prime_polynomial::<F>(modulus.coefficients()),
                format_prime_polynomial::<F>(divisor.coefficients()),
                format_prime_polynomial::<F>(quotient.coefficients())
            ));
            lines.push(
                "consequence: a reducible modulus does not define a field extension in general"
                    .to_string(),
            );
        }
        IrreducibilityStatus::ReducibleWithoutWitness { reason } => {
            lines.push("status: reducible".to_string());
            lines.push(format!(
                "current explanation: {}",
                format_reducibility_reason(reason)
            ));
            lines.push(
                "consequence: a reducible modulus does not define a field extension in general"
                    .to_string(),
            );
        }
    }

    Ok(lines.join("\n"))
}

/// Describes whether a modulus polynomial over the approximate complex
/// backend can define a non-trivial quotient field.
///
/// Because `ComplexApprox` models an algebraically closed field, every
/// polynomial of degree at least `2` is reducible. The current irreducibility
/// API reports that conclusion without fabricating a numerical factorization
/// witness.
fn describe_complex_polynomial_modulus_as_field_modulus(
    modulus: &PolynomialModulus<ComplexApprox>,
) -> Result<String, PolynomialError> {
    let dense_modulus = DensePolynomial::<ComplexApprox>::new(modulus.coefficients().to_vec());
    let status = dense_modulus.irreducibility_status()?;

    let suitability = match &status {
        IrreducibilityStatus::Irreducible => {
            "suitable for a quotient field: yes, this modulus is irreducible over the base field"
        }
        IrreducibilityStatus::Linear => {
            "suitable for a quotient field: technically yes, but the quotient is just the base field again"
        }
        IrreducibilityStatus::Reducible { .. }
        | IrreducibilityStatus::ReducibleWithoutWitness { .. } => {
            "suitable for a quotient field: no, a reducible modulus gives only a quotient algebra in general"
        }
        IrreducibilityStatus::Constant => {
            "suitable for a quotient field: no, a constant polynomial cannot define a meaningful field extension"
        }
    };

    Ok(format!(
        "Field-modulus check over C (approx)\n\
         expression: m(x) = {}\n\
         base field algebraically closed: yes\n\
         irreducibility status: {}\n\
         {}",
        format_complex_polynomial(modulus.coefficients()),
        format_complex_irreducibility_status(&status),
        suitability
    ))
}

/// Explains the irreducibility result for a modulus polynomial over the
/// approximate complex backend.
fn explain_complex_polynomial_modulus_irreducibility(
    modulus: &PolynomialModulus<ComplexApprox>,
) -> Result<String, PolynomialError> {
    let dense_modulus = DensePolynomial::<ComplexApprox>::new(modulus.coefficients().to_vec());
    let status = dense_modulus.irreducibility_status()?;

    let mut lines = vec![
        "Irreducibility check for a field modulus over C (approx)".to_string(),
        format!("modulus: m(x) = {}", format_complex_polynomial(modulus.coefficients())),
        "base field algebraically closed: yes".to_string(),
        "note: in an algebraically closed field, every polynomial of degree at least two is reducible".to_string(),
    ];

    match status {
        IrreducibilityStatus::Constant => {
            lines.push("status: constant".to_string());
            lines.push(
                "consequence: a constant polynomial is not a valid field-extension modulus"
                    .to_string(),
            );
        }
        IrreducibilityStatus::Linear => {
            lines.push("status: linear".to_string());
            lines.push(
                "consequence: a linear modulus is irreducible, but the quotient does not create a genuine new extension"
                    .to_string(),
            );
        }
        IrreducibilityStatus::Irreducible => {
            lines.push("status: irreducible".to_string());
            lines.push(
                "consequence: this modulus is suitable for a quotient-field construction over the base field".to_string(),
            );
        }
        IrreducibilityStatus::Reducible { divisor, quotient } => {
            lines.push("status: reducible".to_string());
            lines.push(format!(
                "witness divisor: {}",
                format_complex_polynomial(divisor.coefficients())
            ));
            lines.push(format!(
                "witness quotient: {}",
                format_complex_polynomial(quotient.coefficients())
            ));
            lines.push(format!(
                "factorization: {} = ({}) * ({})",
                format_complex_polynomial(modulus.coefficients()),
                format_complex_polynomial(divisor.coefficients()),
                format_complex_polynomial(quotient.coefficients())
            ));
            lines.push(
                "consequence: a reducible modulus does not define a field extension in general"
                    .to_string(),
            );
        }
        IrreducibilityStatus::ReducibleWithoutWitness { reason } => {
            lines.push("status: reducible".to_string());
            lines.push(format!(
                "current explanation: {}",
                format_reducibility_reason(reason)
            ));
            lines.push(
                "consequence: a reducible modulus does not define a field extension in general"
                    .to_string(),
            );
        }
    }

    Ok(lines.join("\n"))
}

fn format_irreducibility_status<F: Field>(status: &IrreducibilityStatus<F>) -> String
where
    F::Elem: VisualizableField,
{
    match status {
        IrreducibilityStatus::Constant => "constant".to_string(),
        IrreducibilityStatus::Linear => "linear".to_string(),
        IrreducibilityStatus::Irreducible => "irreducible".to_string(),
        IrreducibilityStatus::Reducible { divisor, quotient } => format!(
            "reducible; witness: {} = ({}) * ({})",
            format_prime_polynomial_modulus(
                &PolynomialModulus::<F>::new(quotient.mul(divisor).coefficients().to_vec())
                    .expect("product of non-trivial factors is a valid non-constant modulus")
            ),
            format_prime_polynomial::<F>(divisor.coefficients()),
            format_prime_polynomial::<F>(quotient.coefficients())
        ),
        IrreducibilityStatus::ReducibleWithoutWitness { reason } => {
            format!("reducible; reason: {}", format_reducibility_reason(*reason))
        }
    }
}

fn format_complex_irreducibility_status(status: &IrreducibilityStatus<ComplexApprox>) -> String {
    match status {
        IrreducibilityStatus::Constant => "constant".to_string(),
        IrreducibilityStatus::Linear => "linear".to_string(),
        IrreducibilityStatus::Irreducible => "irreducible".to_string(),
        IrreducibilityStatus::Reducible { divisor, quotient } => {
            let product = quotient.mul(divisor);
            format!(
                "reducible; witness: {} = ({}) * ({})",
                format_complex_polynomial(product.coefficients()),
                format_complex_polynomial(divisor.coefficients()),
                format_complex_polynomial(quotient.coefficients())
            )
        }
        IrreducibilityStatus::ReducibleWithoutWitness { reason } => {
            format!("reducible; reason: {}", format_reducibility_reason(*reason))
        }
    }
}

/// Formats a quotient representative over `GF(P)` together with its modulus.
fn format_prime_polynomial_field_element<F: Field>(element: &PolynomialFieldElement<F>) -> String
where
    F::Elem: VisualizableField,
{
    format!(
        "[{}] mod ({})",
        format_prime_polynomial::<F>(element.coefficients()),
        format_prime_polynomial::<F>(element.modulus().coefficients())
    )
}

/// Returns a short educational description of a quotient element over `GF(P)`.
fn describe_prime_polynomial_field_element<F: Field>(element: &PolynomialFieldElement<F>) -> String
where
    F::Elem: VisualizableField,
{
    let reduced = element
        .reduce()
        .expect("prime-field quotient reduction should succeed for non-zero modulus");

    format!(
        "Quotient element over GF({})\n\
         representative coefficients (ascending): {}\n\
         representative polynomial: {}\n\
         reduced representative: {}\n\
         already reduced: {}\n\
         reduced degree: {}\n\
         modulus polynomial: {}\n\
         note: arithmetic is interpreted modulo the defining polynomial",
        F::characteristic(),
        format_coefficients::<F>(element.coefficients()),
        format_prime_polynomial::<F>(element.coefficients()),
        format_prime_polynomial::<F>(reduced.coefficients()),
        yes_no(
            element
                .is_reduced()
                .expect("prime-field quotient reduction should succeed"),
        ),
        reduced
            .degree()
            .map_or_else(|| "none (zero)".to_string(), |degree| degree.to_string()),
        format_prime_polynomial_modulus::<F>(element.modulus())
    )
}

/// Explains quotient reduction for an element over `GF(P)`.
fn explain_prime_polynomial_field_reduction<F: Field>(
    element: &PolynomialFieldElement<F>,
) -> Result<String, FieldError>
where
    F::Elem: VisualizableField,
{
    let reduced = element.reduce()?;

    Ok(format!(
        "Reduction in GF({})[x] / (m(x))\n\
         raw representative: {}\n\
         modulus: {}\n\
         reduced representative: {}\n\
         note: the current backend computes the Euclidean remainder modulo the defining polynomial",
        F::characteristic(),
        format_prime_polynomial::<F>(element.coefficients()),
        format_prime_polynomial_modulus::<F>(element.modulus()),
        format_prime_polynomial::<F>(reduced.coefficients())
    ))
}

/// Explains quotient addition over `GF(P)`.
fn explain_prime_polynomial_field_add<F: Field>(
    left: &PolynomialFieldElement<F>,
    right: &PolynomialFieldElement<F>,
) -> Result<String, FieldError>
where
    F::Elem: VisualizableField,
{
    let result = left.add(right)?;
    Ok(format!(
        "Addition in GF({})[x] / (m(x))\n\
         lhs: {}\n\
         rhs: {}\n\
         raw sum in GF({})[x]: ({}) + ({})\n\
         reduced result: {}",
        F::characteristic(),
        format_prime_polynomial_field_element::<F>(left),
        format_prime_polynomial_field_element::<F>(right),
        F::characteristic(),
        format_prime_polynomial::<F>(left.coefficients()),
        format_prime_polynomial::<F>(right.coefficients()),
        format_prime_polynomial_field_element::<F>(&result)
    ))
}

/// Explains quotient multiplication over `GF(P)`.
fn explain_prime_polynomial_field_mul<F: Field>(
    left: &PolynomialFieldElement<F>,
    right: &PolynomialFieldElement<F>,
) -> Result<String, FieldError>
where
    F::Elem: VisualizableField,
{
    let result = left.mul(right)?;
    Ok(format!(
        "Multiplication in GF({})[x] / (m(x))\n\
         lhs: {}\n\
         rhs: {}\n\
         raw product in GF({})[x]: ({}) * ({})\n\
         reduced result: {}\n\
         note: multiplication happens in the polynomial ring first, then the product is reduced modulo m(x)",
        F::characteristic(),
        format_prime_polynomial_field_element::<F>(left),
        format_prime_polynomial_field_element::<F>(right),
        F::characteristic(),
        format_prime_polynomial::<F>(left.coefficients()),
        format_prime_polynomial::<F>(right.coefficients()),
        format_prime_polynomial_field_element::<F>(&result)
    ))
}

/// Explains quotient inversion over `GF(P)`.
fn explain_prime_polynomial_field_inverse<F: Field>(
    element: &PolynomialFieldElement<F>,
) -> Result<String, FieldError>
where
    F::Elem: VisualizableField,
{
    let inverse = element.inverse()?;
    let check = element.mul(&inverse)?;
    Ok(format!(
        "Inverse in GF({})[x] / (m(x))\n\
         element: {}\n\
         inverse: {}\n\
         verification: {} * {} = {}\n\
         note: invertibility depends on the quotient; reducible moduli can admit non-zero non-units",
        F::characteristic(),
        format_prime_polynomial_field_element::<F>(element),
        format_prime_polynomial_field_element::<F>(&inverse),
        format_prime_polynomial::<F>(element.coefficients()),
        format_prime_polynomial::<F>(inverse.coefficients()),
        format_prime_polynomial_field_element::<F>(&check)
    ))
}

impl<F: Field> Visualizable for PolynomialModulus<F>
where
    F::Elem: VisualizableField,
{
    fn format_compact(&self) -> String {
        format_prime_polynomial_modulus::<F>(self)
    }

    fn describe(&self) -> String {
        describe_prime_polynomial_modulus::<F>(self)
    }
}

impl<F: Field> Visualizable for PolynomialFieldElement<F>
where
    F::Elem: VisualizableField,
{
    fn format_compact(&self) -> String {
        format_prime_polynomial_field_element::<F>(self)
    }

    fn describe(&self) -> String {
        describe_prime_polynomial_field_element::<F>(self)
    }
}

impl<F: Field> VisualizableField for PolynomialFieldElement<F>
where
    F::Elem: VisualizableField,
{
    fn format_elem(&self) -> String {
        format_prime_polynomial_field_element::<F>(self)
    }

    fn inverse(&self) -> Option<String> {
        self.inverse()
            .ok()
            .map(|value| format_prime_polynomial_field_element::<F>(&value))
    }

    fn explain_add(lhs: &Self, rhs: &Self) -> Option<String> {
        explain_prime_polynomial_field_add::<F>(lhs, rhs).ok()
    }

    fn explain_mul(lhs: &Self, rhs: &Self) -> Option<String> {
        explain_prime_polynomial_field_mul::<F>(lhs, rhs).ok()
    }

    fn explain_div(lhs: &Self, rhs: &Self) -> Option<String> {
        let reciprocal = rhs.inverse().ok()?;
        let result = lhs.div(rhs).ok()?;

        Some(format!(
            "Division in GF({})[x] / (m(x))\n\
             lhs: {}\n\
             rhs: {}\n\
             reciprocal of rhs: {}\n\
             reduced result: {}",
            F::characteristic(),
            format_prime_polynomial_field_element::<F>(lhs),
            format_prime_polynomial_field_element::<F>(rhs),
            format_prime_polynomial_field_element::<F>(&reciprocal),
            format_prime_polynomial_field_element::<F>(&result)
        ))
    }
}

#[cfg(test)]
mod tests;
