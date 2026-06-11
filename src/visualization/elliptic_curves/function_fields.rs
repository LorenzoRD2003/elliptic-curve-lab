use core::fmt;

use crate::elliptic_curves::{ShortWeierstrassFunction, ShortWeierstrassFunctionField};
use crate::fields::Field;
use crate::visualization::elliptic_curves::short_weierstrass::format_curve;
use crate::visualization::fields::format_rational_function;
use crate::visualization::fields::traits::VisualizableField;
use crate::visualization::traits::Visualizable;

/// Formats one short-Weierstrass function-field element compactly.
///
/// The current presentation is specific to the basis `1, y` over `F(x)`, so
/// the compact formatter writes the element as
///
/// `A(x) + y*B(x)`.
pub fn format_short_weierstrass_function<F>(function: &ShortWeierstrassFunction<F>) -> String
where
    F: Field,
    F::Elem: VisualizableField,
{
    let a_text = format_rational_function(function.a_part());
    let b_text = format_rational_function(function.b_part());

    if function.b_part().is_zero() {
        return a_text;
    }

    if function.a_part().is_zero() && function.b_part().is_one() {
        return "y".to_string();
    }

    if function.a_part().is_zero() {
        return format!("y*({b_text})");
    }

    if function.b_part().is_one() {
        return format!("{a_text} + y");
    }

    format!("{a_text} + y*({b_text})")
}

/// Returns a richer educational description of one function-field element.
pub fn describe_short_weierstrass_function<F>(function: &ShortWeierstrassFunction<F>) -> String
where
    F: Field,
    F::Elem: VisualizableField + fmt::Display,
{
    format!(
        "Short-Weierstrass function-field element\n\
         curve: {}\n\
         basis over F(x): 1, y\n\
         a-part A(x): {}\n\
         b-part B(x): {}\n\
         compact form: {}\n\
         zero: {}\n\
         one: {}\n\
         note: this represents A(x) + yB(x) and uses the specific short-Weierstrass relation y^2 = x^3 + ax + b",
        format_curve(function.curve()),
        format_rational_function(function.a_part()),
        format_rational_function(function.b_part()),
        format_short_weierstrass_function(function),
        if function.is_zero() { "yes" } else { "no" },
        if function.is_one() { "yes" } else { "no" }
    )
}

/// Returns a short educational description of the ambient field `F(E)` of one
/// concrete short-Weierstrass curve.
pub fn describe_short_weierstrass_function_field<F>(
    field: &ShortWeierstrassFunctionField<F>,
) -> String
where
    F: Field,
    F::Elem: VisualizableField + fmt::Display,
{
    format!(
        "Short-Weierstrass function field\n\
         curve: {}\n\
         vector-space presentation: F(E) = F(x) ⊕ yF(x)\n\
         x-generator: {}\n\
         y-generator: {}\n\
         note: this ambient field is runtime-dependent because the multiplication rule uses the concrete curve equation y^2 = x^3 + ax + b",
        format_curve(field.curve()),
        format_short_weierstrass_function(&field.x()),
        format_short_weierstrass_function(&field.y())
    )
}

/// Explains the conjugation involution `y ↦ -y`.
pub fn explain_short_weierstrass_function_conjugate<F>(
    function: &ShortWeierstrassFunction<F>,
) -> String
where
    F: Field,
    F::Elem: VisualizableField + fmt::Display,
{
    let conjugate = function.conjugate();

    format!(
        "Conjugation in F(E)\n\
         curve: {}\n\
         element: {}\n\
         rule: A(x) + yB(x) ↦ A(x) - yB(x)\n\
         conjugate: {}",
        format_curve(function.curve()),
        format_short_weierstrass_function(function),
        format_short_weierstrass_function(&conjugate)
    )
}

/// Explains the norm `A^2 - fB^2`.
pub fn explain_short_weierstrass_function_norm<F>(function: &ShortWeierstrassFunction<F>) -> String
where
    F: Field,
    F::Elem: VisualizableField + fmt::Display,
{
    let norm = function.norm();

    format!(
        "Norm in F(E)\n\
         curve: {}\n\
         element: {}\n\
         formula: N(A, B) = A(x)^2 - f(x)B(x)^2 with f(x) = x^3 + ax + b\n\
         computed norm: {}",
        format_curve(function.curve()),
        format_short_weierstrass_function(function),
        format_rational_function(&norm)
    )
}

/// Explains addition in the basis `1, y`.
pub fn explain_short_weierstrass_function_add<F>(
    left: &ShortWeierstrassFunction<F>,
    right: &ShortWeierstrassFunction<F>,
) -> Result<String, crate::elliptic_curves::CurveError>
where
    F: Field,
    F::Elem: VisualizableField + fmt::Display,
{
    let result = left.add(right)?;

    Ok(format!(
        "Addition in F(E)\n\
         curve: {}\n\
         lhs: {}\n\
         rhs: {}\n\
         component rule: (A, B) + (C, D) = (A + C, B + D)\n\
         reduced result: {}",
        format_curve(left.curve()),
        format_short_weierstrass_function(left),
        format_short_weierstrass_function(right),
        format_short_weierstrass_function(&result)
    ))
}

/// Explains multiplication using the short-Weierstrass reduction rule.
pub fn explain_short_weierstrass_function_mul<F>(
    left: &ShortWeierstrassFunction<F>,
    right: &ShortWeierstrassFunction<F>,
) -> Result<String, crate::elliptic_curves::CurveError>
where
    F: Field,
    F::Elem: VisualizableField + fmt::Display,
{
    let result = left.mul(right)?;

    Ok(format!(
        "Multiplication in F(E)\n\
         curve: {}\n\
         lhs: {}\n\
         rhs: {}\n\
         relation used: y^2 = x^3 + ax + b\n\
         basis rule: (A, B)(C, D) = (AC + fBD, AD + BC)\n\
         reduced result: {}",
        format_curve(left.curve()),
        format_short_weierstrass_function(left),
        format_short_weierstrass_function(right),
        format_short_weierstrass_function(&result)
    ))
}

/// Explains the inverse via conjugate over norm.
pub fn explain_short_weierstrass_function_inverse<F>(
    function: &ShortWeierstrassFunction<F>,
) -> Result<String, crate::elliptic_curves::CurveError>
where
    F: Field,
    F::Elem: VisualizableField + fmt::Display,
{
    let conjugate = function.conjugate();
    let norm = function.norm();
    let inverse = function.inverse()?;

    Ok(format!(
        "Inverse in F(E)\n\
         curve: {}\n\
         element: {}\n\
         conjugate: {}\n\
         norm: {}\n\
         formula: (A, B)^(-1) = (A / (A^2 - fB^2), -B / (A^2 - fB^2))\n\
         inverse: {}",
        format_curve(function.curve()),
        format_short_weierstrass_function(function),
        format_short_weierstrass_function(&conjugate),
        format_rational_function(&norm),
        format_short_weierstrass_function(&inverse)
    ))
}

/// Explains the derivative in the basis `1, y` over `F(x)`.
pub fn explain_short_weierstrass_function_derivative<F>(
    function: &ShortWeierstrassFunction<F>,
) -> String
where
    F: Field,
    F::Elem: VisualizableField + fmt::Display,
{
    let derivative = function.derivative();

    format!(
        "Derivative in F(E)\n\
         curve: {}\n\
         element: {}\n\
         implicit relation: y^2 = x^3 + ax + b = f(x)\n\
         derived rule: (A, B)' = (A', B' + f'(x)B(x)/(2f(x)))\n\
         basis interpretation: d/dx(A(x) + yB(x)) = A'(x) + y*(B'(x) + f'(x)B(x)/(2f(x)))\n\
         reduced derivative: {}",
        format_curve(function.curve()),
        format_short_weierstrass_function(function),
        format_short_weierstrass_function(&derivative)
    )
}

impl<F> Visualizable for ShortWeierstrassFunction<F>
where
    F: Field,
    F::Elem: VisualizableField + fmt::Display,
{
    fn format_compact(&self) -> String {
        format_short_weierstrass_function(self)
    }

    fn describe(&self) -> String {
        describe_short_weierstrass_function(self)
    }
}

#[cfg(test)]
mod tests {
    use crate::elliptic_curves::{
        ShortWeierstrassCurve, ShortWeierstrassFunction, ShortWeierstrassFunctionField,
    };
    use crate::fields::{Field, Fp, RationalFunction};
    use crate::polynomials::DensePolynomial;
    use crate::visualization::elliptic_curves::{
        describe_short_weierstrass_function, describe_short_weierstrass_function_field,
        explain_short_weierstrass_function_add, explain_short_weierstrass_function_conjugate,
        explain_short_weierstrass_function_derivative, explain_short_weierstrass_function_inverse,
        explain_short_weierstrass_function_mul, explain_short_weierstrass_function_norm,
        format_short_weierstrass_function,
    };
    use crate::visualization::traits::Visualizable;

    type F17 = Fp<17>;

    fn f17_dense(values: &[u64]) -> DensePolynomial<F17> {
        DensePolynomial::<F17>::new(values.iter().copied().map(F17::elem_from_u64).collect())
    }

    fn curve() -> ShortWeierstrassCurve<F17> {
        ShortWeierstrassCurve::<F17>::new(F17::elem_from_u64(2), F17::elem_from_u64(3))
            .expect("curve should be nonsingular")
    }

    #[test]
    fn function_formatter_handles_a_part_only_y_and_full_pair() {
        let ambient = ShortWeierstrassFunctionField::<F17>::new(curve());
        let polynomial_only = ShortWeierstrassFunction::<F17>::from_rational_function(
            curve(),
            RationalFunction::<F17>::from_polynomial(f17_dense(&[1, 1])),
        );
        let y = ambient.y();
        let full = ShortWeierstrassFunction::<F17>::new(
            curve(),
            RationalFunction::<F17>::from_polynomial(f17_dense(&[1])),
            RationalFunction::<F17>::from_polynomial(f17_dense(&[1, 1])),
        );

        assert_eq!(format_short_weierstrass_function(&polynomial_only), "x + 1");
        assert_eq!(format_short_weierstrass_function(&y), "y");
        assert_eq!(format_short_weierstrass_function(&full), "1 + y*(x + 1)");
    }

    #[test]
    fn function_description_mentions_basis_and_curve_specific_relation() {
        let function = ShortWeierstrassFunction::<F17>::new(
            curve(),
            RationalFunction::<F17>::from_polynomial(f17_dense(&[1])),
            RationalFunction::<F17>::from_polynomial(f17_dense(&[0, 1])),
        );
        let description = describe_short_weierstrass_function(&function);

        assert!(description.contains("Short-Weierstrass function-field element"));
        assert!(description.contains("basis over F(x): 1, y"));
        assert!(description.contains("specific short-Weierstrass relation"));
    }

    #[test]
    fn field_description_mentions_fx_plus_yfx_and_generators() {
        let field = ShortWeierstrassFunctionField::<F17>::new(curve());
        let description = describe_short_weierstrass_function_field(&field);

        assert!(description.contains("F(E) = F(x) ⊕ yF(x)"));
        assert!(description.contains("x-generator: x"));
        assert!(description.contains("y-generator: y"));
        assert!(description.contains("runtime-dependent"));
    }

    #[test]
    fn conjugation_norm_add_mul_and_inverse_explanations_mention_core_formulas() {
        let left = ShortWeierstrassFunction::<F17>::new(
            curve(),
            RationalFunction::<F17>::from_polynomial(f17_dense(&[1])),
            RationalFunction::<F17>::from_polynomial(f17_dense(&[1])),
        );
        let right = ShortWeierstrassFunction::<F17>::new(
            curve(),
            RationalFunction::<F17>::from_polynomial(f17_dense(&[1, 1])),
            RationalFunction::<F17>::from_polynomial(f17_dense(&[0, 1])),
        );

        assert!(explain_short_weierstrass_function_conjugate(&left).contains("A(x) + yB(x)"));
        assert!(explain_short_weierstrass_function_norm(&left).contains("A(x)^2 - f(x)B(x)^2"));
        assert!(
            explain_short_weierstrass_function_add(&left, &right)
                .expect("same-curve add")
                .contains("(A, B) + (C, D) = (A + C, B + D)")
        );
        assert!(
            explain_short_weierstrass_function_mul(&left, &right)
                .expect("same-curve mul")
                .contains("(A, B)(C, D) = (AC + fBD, AD + BC)")
        );
        assert!(
            explain_short_weierstrass_function_inverse(&left)
                .expect("left should be invertible")
                .contains("(A, B)^(-1)")
        );
    }

    #[test]
    fn derivative_explanation_mentions_implicit_rule_and_basis_formula() {
        let function = ShortWeierstrassFunction::<F17>::new(
            curve(),
            RationalFunction::<F17>::from_polynomial(f17_dense(&[1])),
            RationalFunction::<F17>::from_polynomial(f17_dense(&[0, 1])),
        );
        let explanation = explain_short_weierstrass_function_derivative(&function);

        assert!(explanation.contains("Derivative in F(E)"));
        assert!(explanation.contains("implicit relation: y^2 = x^3 + ax + b = f(x)"));
        assert!(explanation.contains("(A, B)' = (A', B' + f'(x)B(x)/(2f(x)))"));
        assert!(explanation.contains("reduced derivative:"));
    }

    #[test]
    fn visualizable_trait_reuses_function_field_helpers() {
        let function = ShortWeierstrassFunction::<F17>::new(
            curve(),
            RationalFunction::<F17>::from_polynomial(f17_dense(&[1])),
            RationalFunction::<F17>::from_polynomial(f17_dense(&[1])),
        );

        assert_eq!(function.format_compact(), "1 + y");
        assert!(
            function
                .describe()
                .contains("Short-Weierstrass function-field element")
        );
    }
}
