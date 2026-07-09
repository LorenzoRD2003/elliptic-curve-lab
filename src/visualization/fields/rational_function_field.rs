use crate::fields::{rational_function_field::RationalFunction, traits::Field};
use crate::visualization::{
    Visualizable, VisualizableField, polynomials::dense::format_dense_polynomial, shared::yes_no,
};

/// Formats a rational function as a compact quotient of dense polynomials.
pub(crate) fn format_rational_function<F: Field>(function: &RationalFunction<F>) -> String
where
    F::Elem: VisualizableField,
{
    if function.denominator().degree() == Some(0) && function.denominator().is_monic() {
        return format_dense_polynomial(function.numerator());
    }

    format!(
        "({}) / ({})",
        format_dense_polynomial(function.numerator()),
        format_dense_polynomial(function.denominator())
    )
}

/// Returns a richer educational description of a rational function value.
fn describe_rational_function<F: Field>(function: &RationalFunction<F>) -> String
where
    F::Elem: VisualizableField,
{
    format!(
        "Rational function over a field\n\
         numerator: {}\n\
         denominator: {}\n\
         compact form: {}\n\
         zero: {}\n\
         one: {}\n\
         note: the stored representation is canonical, gcd-reduced, and uses a monic denominator",
        format_dense_polynomial(function.numerator()),
        format_dense_polynomial(function.denominator()),
        format_rational_function(function),
        yes_no(function.is_zero()),
        yes_no(function.is_one())
    )
}

/// Returns a short educational description of the rational function field
/// family `F(x)`.
fn describe_rational_function_field<F: Field>() -> String
where
    F::Elem: VisualizableField,
{
    format!(
        "Rational function field\n\
         presentation: F(x) over a base field of characteristic {}\n\
         distinguished indeterminate: {}\n\
         algebraically closed: {}\n\
         note: values are stored as gcd-reduced quotients of dense polynomials with monic denominator",
        crate::fields::rational_function_field::RationalFunctionField::<F>::characteristic(),
        format_rational_function(&crate::fields::rational_function_field::RationalFunctionField::<F>::indeterminate()),
        yes_no(crate::fields::rational_function_field::RationalFunctionField::<F>::IS_ALGEBRAICALLY_CLOSED)
    )
}

/// Explains addition of two rational functions.
fn explain_rational_function_add<F: Field>(
    lhs: &RationalFunction<F>,
    rhs: &RationalFunction<F>,
) -> String
where
    F::Elem: VisualizableField,
{
    let result = lhs.add(rhs);
    format!(
        "Addition in a rational function field\n\
         lhs: {}\n\
         rhs: {}\n\
         cross-multiplied numerator: ({}) * ({}) + ({}) * ({})\n\
         common denominator before reduction: ({}) * ({})\n\
         reduced result: {}",
        format_rational_function(lhs),
        format_rational_function(rhs),
        format_dense_polynomial(lhs.numerator()),
        format_dense_polynomial(rhs.denominator()),
        format_dense_polynomial(rhs.numerator()),
        format_dense_polynomial(lhs.denominator()),
        format_dense_polynomial(lhs.denominator()),
        format_dense_polynomial(rhs.denominator()),
        format_rational_function(&result)
    )
}

/// Explains multiplication of two rational functions.
fn explain_rational_function_mul<F: Field>(
    lhs: &RationalFunction<F>,
    rhs: &RationalFunction<F>,
) -> String
where
    F::Elem: VisualizableField,
{
    let result = lhs.mul(rhs);
    format!(
        "Multiplication in a rational function field\n\
         lhs: {}\n\
         rhs: {}\n\
         raw numerator product: ({}) * ({})\n\
         raw denominator product: ({}) * ({})\n\
         reduced result: {}",
        format_rational_function(lhs),
        format_rational_function(rhs),
        format_dense_polynomial(lhs.numerator()),
        format_dense_polynomial(rhs.numerator()),
        format_dense_polynomial(lhs.denominator()),
        format_dense_polynomial(rhs.denominator()),
        format_rational_function(&result)
    )
}

/// Explains multiplicative inversion of a rational function.
fn explain_rational_function_inverse<F: Field>(function: &RationalFunction<F>) -> Option<String>
where
    F::Elem: VisualizableField,
{
    let inverse = function.inverse().ok()?;
    Some(format!(
        "Inverse in a rational function field\n\
         element: {}\n\
         inverse: {}\n\
         verification: {} * {} = {}",
        format_rational_function(function),
        format_rational_function(&inverse),
        format_rational_function(function),
        format_rational_function(&inverse),
        format_rational_function(&function.mul(&inverse))
    ))
}

/// Explains division of two rational functions.
fn explain_rational_function_div<F: Field>(
    lhs: &RationalFunction<F>,
    rhs: &RationalFunction<F>,
) -> Option<String>
where
    F::Elem: VisualizableField,
{
    let reciprocal = rhs.inverse().ok()?;
    let result = lhs.div(rhs).ok()?;
    Some(format!(
        "Division in a rational function field\n\
         lhs: {}\n\
         rhs: {}\n\
         reciprocal of rhs: {}\n\
         reduced result: {}",
        format_rational_function(lhs),
        format_rational_function(rhs),
        format_rational_function(&reciprocal),
        format_rational_function(&result)
    ))
}

/// Explains formal differentiation of a rational function.
fn explain_rational_function_derivative<F: Field>(function: &RationalFunction<F>) -> String
where
    F::Elem: VisualizableField,
{
    let result = function.derivative();
    format!(
        "Derivative in a rational function field\n\
         function: {}\n\
         quotient rule numerator: ({})' * ({}) - ({}) * ({})'\n\
         quotient rule denominator: ({})^2\n\
         reduced derivative: {}\n\
         note: the backend stores the derivative in canonical reduced form",
        format_rational_function(function),
        format_dense_polynomial(function.numerator()),
        format_dense_polynomial(function.denominator()),
        format_dense_polynomial(function.numerator()),
        format_dense_polynomial(function.denominator()),
        format_dense_polynomial(function.denominator()),
        format_rational_function(&result)
    )
}

impl<F> Visualizable for RationalFunction<F>
where
    F: Field,
    F::Elem: VisualizableField,
{
    fn format_compact(&self) -> String {
        format_rational_function(self)
    }

    fn describe(&self) -> String {
        describe_rational_function(self)
    }
}

impl<F: Field> VisualizableField for RationalFunction<F>
where
    F::Elem: VisualizableField,
{
    fn format_elem(&self) -> String {
        format_rational_function(self)
    }

    fn inverse(&self) -> Option<String> {
        self.inverse()
            .ok()
            .map(|value| format_rational_function(&value))
    }

    fn explain_add(lhs: &Self, rhs: &Self) -> Option<String> {
        Some(explain_rational_function_add(lhs, rhs))
    }

    fn explain_mul(lhs: &Self, rhs: &Self) -> Option<String> {
        Some(explain_rational_function_mul(lhs, rhs))
    }

    fn explain_div(lhs: &Self, rhs: &Self) -> Option<String> {
        explain_rational_function_div(lhs, rhs)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fields::{
        Q,
        rational_function_field::{RationalFunction, RationalFunctionField},
    };
    use crate::polynomials::DensePolynomial;
    use crate::visualization::{Visualizable, VisualizableField};

    type F17 = crate::fields::Fp17;

    fn f17_dense(values: &[u64]) -> DensePolynomial<F17> {
        DensePolynomial::<F17>::new(values.iter().copied().map(F17::from_i64).collect())
    }

    fn q_dense(values: &[(i64, i64)]) -> DensePolynomial<Q> {
        DensePolynomial::<Q>::new(
            values
                .iter()
                .map(|&(numerator, denominator)| {
                    Q::div(&Q::from_i64(numerator), &Q::from_i64(denominator))
                        .expect("denominator should be non-zero")
                })
                .collect(),
        )
    }

    #[test]
    fn rational_function_formatter_handles_polynomials_and_true_quotients() {
        let polynomial = RationalFunction::<F17>::from_polynomial(f17_dense(&[2, 5, 1]));
        let quotient = RationalFunction::<F17>::new(f17_dense(&[1]), f17_dense(&[1, 1]))
            .expect("quotient should exist");

        assert_eq!(format_rational_function(&polynomial), "x^2 + 5*x + 2");
        assert_eq!(format_rational_function(&quotient), "(1) / (x + 1)");
    }

    #[test]
    fn rational_function_description_mentions_canonical_parts() {
        let function = RationalFunction::<F17>::new(f17_dense(&[1, 1]), f17_dense(&[1, 2]))
            .expect("function should exist");
        let description = describe_rational_function(&function);

        assert!(description.contains("Rational function over a field"));
        assert!(description.contains("numerator: 9*x + 9"));
        assert!(description.contains("denominator: x + 9"));
        assert!(description.contains("monic denominator"));
    }

    #[test]
    fn rational_function_field_description_mentions_x_and_characteristic() {
        let description = describe_rational_function_field::<F17>();

        assert!(description.contains("Rational function field"));
        assert!(description.contains("characteristic 17"));
        assert!(description.contains("distinguished indeterminate: x"));
        assert!(description.contains("algebraically closed: no"));
    }

    #[test]
    fn rational_function_operation_explanations_mention_core_story() {
        let left =
            RationalFunction::<F17>::new(f17_dense(&[1]), f17_dense(&[1, 1])).expect("left exists");
        let right = RationalFunction::<F17>::new(f17_dense(&[1]), f17_dense(&[1, 16]))
            .expect("right exists");

        let add = explain_rational_function_add(&left, &right);
        let mul = explain_rational_function_mul(&left, &right);
        let inv = explain_rational_function_inverse(&left).expect("left is invertible");
        let div = explain_rational_function_div(&left, &right).expect("division should work");

        assert!(add.contains("cross-multiplied numerator"));
        assert!(add.contains("reduced result"));
        assert!(mul.contains("raw numerator product"));
        assert!(inv.contains("verification:"));
        assert!(div.contains("reciprocal of rhs"));
    }

    #[test]
    fn rational_function_derivative_explanation_mentions_quotient_rule() {
        let function =
            RationalFunction::<Q>::new(q_dense(&[(0, 1), (1, 1)]), q_dense(&[(1, 1), (1, 1)]))
                .expect("exists");
        let explanation = explain_rational_function_derivative(&function);

        assert!(explanation.contains("Derivative in a rational function field"));
        assert!(explanation.contains("quotient rule numerator"));
        assert!(explanation.contains("reduced derivative:"));
    }

    #[test]
    fn rational_function_visualizable_traits_reuse_field_helpers() {
        let function = RationalFunction::<F17>::new(f17_dense(&[1]), f17_dense(&[1, 1]))
            .expect("function should exist");

        assert_eq!(function.format_compact(), "(1) / (x + 1)");
        assert!(
            function
                .describe()
                .contains("Rational function over a field")
        );
        assert_eq!(function.format_elem(), "(1) / (x + 1)");
        assert!(
            RationalFunction::<F17>::explain_add(&function, &function)
                .expect("explanation should exist")
                .contains("Addition in a rational function field")
        );
    }

    #[test]
    fn rational_function_visualizable_division_returns_none_for_zero_divisor() {
        let lhs = RationalFunction::<F17>::from_polynomial(f17_dense(&[1]));
        let rhs = RationalFunction::<F17>::constant(F17::zero());

        assert!(RationalFunction::<F17>::explain_div(&lhs, &rhs).is_none());
    }

    #[test]
    fn rational_function_field_indeterminate_is_visualized_compactly() {
        let x = RationalFunctionField::<F17>::indeterminate();
        assert_eq!(format_rational_function(&x), "x");
    }
}
