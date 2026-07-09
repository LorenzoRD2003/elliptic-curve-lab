use core::fmt;

use crate::elliptic_curves::{
    CurveError,
    short_weierstrass::function_fields::{ShortWeierstrassFunction, ShortWeierstrassFunctionField},
};
use crate::fields::traits::{Field, FiniteField, PthRootExtraction};
use crate::visualization::{
    Visualizable, VisualizableField, elliptic_curves::short_weierstrass::format_curve,
    fields::rational_function_field::format_rational_function, shared::yes_no,
};

/// Formats one short-Weierstrass function-field element compactly.
///
/// The current presentation is specific to the basis `1, y` over `F(x)`, so
/// the compact formatter writes the element as
///
/// `A(x) + y*B(x)`.
pub(crate) fn format_short_weierstrass_function<F: Field>(
    function: &ShortWeierstrassFunction<F>,
) -> String
where
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
fn describe_short_weierstrass_function<F: Field>(function: &ShortWeierstrassFunction<F>) -> String
where
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
        yes_no(function.is_zero()),
        yes_no(function.is_one())
    )
}

/// Returns a short educational description of the ambient field `F(E)` of one
/// concrete short-Weierstrass curve.
pub(crate) fn describe_short_weierstrass_function_field<F: Field>(
    field: &ShortWeierstrassFunctionField<F>,
) -> String
where
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
fn explain_short_weierstrass_function_conjugate<F: Field>(
    function: &ShortWeierstrassFunction<F>,
) -> String
where
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
fn explain_short_weierstrass_function_norm<F: Field>(
    function: &ShortWeierstrassFunction<F>,
) -> String
where
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
fn explain_short_weierstrass_function_add<F: Field>(
    left: &ShortWeierstrassFunction<F>,
    right: &ShortWeierstrassFunction<F>,
) -> Result<String, crate::elliptic_curves::CurveError>
where
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
fn explain_short_weierstrass_function_mul<F: Field>(
    left: &ShortWeierstrassFunction<F>,
    right: &ShortWeierstrassFunction<F>,
) -> Result<String, CurveError>
where
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
fn explain_short_weierstrass_function_inverse<F: Field>(
    function: &ShortWeierstrassFunction<F>,
) -> Result<String, crate::elliptic_curves::CurveError>
where
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
fn explain_short_weierstrass_function_derivative<F: Field>(
    function: &ShortWeierstrassFunction<F>,
) -> String
where
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

/// Explains `p`-th-root extraction in the short-Weierstrass function field
/// `F(E) = F(x) ⊕ yF(x)`.
fn explain_short_weierstrass_function_pth_root<F: FiniteField>(
    function: &ShortWeierstrassFunction<F>,
) -> String
where
    F::Elem: PthRootExtraction + VisualizableField + fmt::Display,
{
    let p = F::characteristic().to_biguint();

    match function.pth_root() {
        Some(root) => format!(
            "Characteristic-p root extraction in F(E)\n\
             curve: {}\n\
             characteristic: {}\n\
             function: {}\n\
             root: {}\n\
             criterion: if u = A(x) + yB(x), then u^p = A(x)^p + y*f(x)^((p-1)/2)B(x)^p with f(x) = x^3 + ax + b\n\
             interpretation: the root exists because A(x) and B(x)/f(x)^((p-1)/2) both admit p-th roots in F(x)",
            format_curve(function.curve()),
            p,
            format_short_weierstrass_function(function),
            format_short_weierstrass_function(&root),
        ),
        None => format!(
            "Characteristic-p root extraction in F(E)\n\
             curve: {}\n\
             characteristic: {}\n\
             function: {}\n\
             p-th root: does not exist in F(E)\n\
             criterion: if u = A(x) + yB(x), then a root would require A(x) to be a p-th power in F(x) and B(x)/f(x)^((p-1)/2) to be a p-th power in F(x)",
            format_curve(function.curve()),
            p,
            format_short_weierstrass_function(function),
        ),
    }
}

impl<F: Field> Visualizable for ShortWeierstrassFunction<F>
where
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
    use crate::fields::traits::*;

    use super::{
        describe_short_weierstrass_function, describe_short_weierstrass_function_field,
        explain_short_weierstrass_function_add, explain_short_weierstrass_function_conjugate,
        explain_short_weierstrass_function_derivative, explain_short_weierstrass_function_inverse,
        explain_short_weierstrass_function_mul, explain_short_weierstrass_function_norm,
        explain_short_weierstrass_function_pth_root, format_short_weierstrass_function,
    };
    use crate::elliptic_curves::{
        ShortWeierstrassCurve, short_weierstrass::function_fields::ShortWeierstrassFunction,
    };
    use crate::fields::rational_function_field::RationalFunction;
    use crate::polynomials::DensePolynomial;
    use crate::visualization::traits::Visualizable;

    type F17 = crate::fields::Fp17;

    fn f17_dense(values: &[u64]) -> DensePolynomial<F17> {
        DensePolynomial::<F17>::new(values.iter().copied().map(F17::from_i64).collect())
    }

    fn curve() -> ShortWeierstrassCurve<F17> {
        ShortWeierstrassCurve::<F17>::new(F17::from_i64(2), F17::from_i64(3))
            .expect("curve should be nonsingular")
    }

    #[test]
    fn function_formatter_handles_a_part_only_y_and_full_pair() {
        let ambient = crate::elliptic_curves::short_weierstrass::function_fields::ShortWeierstrassFunctionField::<F17>::new(curve());
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
        let field = crate::elliptic_curves::short_weierstrass::function_fields::ShortWeierstrassFunctionField::<F17>::new(curve());
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
    fn pth_root_explanation_mentions_the_characteristic_p_formula_and_failure_story() {
        let field = crate::elliptic_curves::short_weierstrass::function_fields::ShortWeierstrassFunctionField::<F17>::new(curve());
        let rhs = RationalFunction::<F17>::from_polynomial(f17_dense(&[3, 2, 0, 1]));
        let mut rhs_to_the_eighth = RationalFunction::<F17>::constant(F17::one());
        for _ in 0..8 {
            rhs_to_the_eighth = rhs_to_the_eighth.mul(&rhs);
        }
        let y_to_the_p = ShortWeierstrassFunction::<F17>::new(
            field.curve().clone(),
            RationalFunction::<F17>::constant(F17::zero()),
            rhs_to_the_eighth,
        );
        let success = explain_short_weierstrass_function_pth_root(&y_to_the_p);
        let failure = explain_short_weierstrass_function_pth_root(&field.x());

        assert!(success.contains("u^p = A(x)^p + y*f(x)^((p-1)/2)B(x)^p"));
        assert!(success.contains("root:"));
        assert!(failure.contains("does not exist"));
        assert!(failure.contains("B(x)/f(x)^((p-1)/2)"));
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
