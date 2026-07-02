use crate::fields::traits::*;
use crate::fields::{
    FieldError,
    extension_field::{ExtensionField, ExtensionFieldElement, ExtensionFieldSpec},
};
use crate::visualization::VisualizableField;
use crate::visualization::traits::Visualizable;

/// Formats the quotient presentation of a statically specified extension
/// field.
pub fn format_extension_field<S>() -> String
where
    S: ExtensionFieldSpec,
    BaseElem<S>: VisualizableField,
{
    format!(
        "{} = F[x] / ({})",
        ExtensionField::<S>::name(),
        format_extension_polynomial::<S>(ExtensionField::<S>::modulus().coefficients())
    )
}

/// Returns a richer educational description of an extension field family.
pub fn describe_extension_field<S>() -> String
where
    S: ExtensionFieldSpec,
    BaseElem<S>: VisualizableField,
{
    format!(
        "Extension field\n\
         name: {}\n\
         presentation: {}\n\
         degree over immediate base field: {}\n\
         algebraically closed: {}\n\
         note: representatives are reduced modulo the defining polynomial stored in the type-level specification",
        ExtensionField::<S>::name(),
        format_extension_field::<S>(),
        ExtensionField::<S>::extension_degree().get(),
        if <ExtensionField<S> as Field>::IS_ALGEBRAICALLY_CLOSED {
            "yes"
        } else {
            "no"
        }
    )
}

/// Formats a canonical quotient representative using the ambient static
/// extension specification.
pub fn format_extension_field_element<S>(element: &ExtensionFieldElement<S>) -> String
where
    S: ExtensionFieldSpec,
    BaseElem<S>: VisualizableField,
{
    format!(
        "[{}] mod ({})",
        format_extension_polynomial::<S>(element.coefficients()),
        format_extension_polynomial::<S>(ExtensionField::<S>::modulus().coefficients())
    )
}

/// Returns a richer educational description of an extension-field element.
pub fn describe_extension_field_element<S>(element: &ExtensionFieldElement<S>) -> String
where
    S: ExtensionFieldSpec,
    BaseElem<S>: VisualizableField,
{
    let reduced = ExtensionField::<S>::reduce(element);

    format!(
        "Extension-field element\n\
         ambient field: {}\n\
         raw representative: {}\n\
         reduced representative: {}\n\
         reduced degree: {}\n\
         zero representative: {}",
        ExtensionField::<S>::name(),
        format_extension_polynomial::<S>(element.coefficients()),
        format_extension_field_element::<S>(&reduced),
        reduced
            .degree()
            .map_or_else(|| "none (zero)".to_string(), |degree| degree.to_string()),
        if reduced.is_zero() { "yes" } else { "no" }
    )
}

/// Explains quotient reduction modulo the defining polynomial.
pub fn explain_extension_field_reduction<S>(element: &ExtensionFieldElement<S>) -> String
where
    S: ExtensionFieldSpec,
    BaseElem<S>: VisualizableField,
{
    let reduced = ExtensionField::<S>::reduce(element);

    format!(
        "Reduction in an extension-field quotient\n\
         ambient field: {}\n\
         raw representative: {}\n\
         modulus relation: {} = 0 in the quotient\n\
         reduced representative: {}\n\
         note: the backend computes the Euclidean remainder modulo the defining polynomial",
        format_extension_field::<S>(),
        format_extension_polynomial::<S>(element.coefficients()),
        format_extension_polynomial::<S>(ExtensionField::<S>::modulus().coefficients()),
        format_extension_field_element::<S>(&reduced)
    )
}

/// Explains addition inside the quotient.
pub fn explain_extension_field_add<S>(
    left: &ExtensionFieldElement<S>,
    right: &ExtensionFieldElement<S>,
) -> String
where
    S: ExtensionFieldSpec,
    BaseElem<S>: VisualizableField,
{
    let left_reduced = ExtensionField::<S>::reduce(left);
    let right_reduced = ExtensionField::<S>::reduce(right);
    let result = ExtensionField::<S>::add(left, right);

    format!(
        "Addition in an extension-field quotient\n\
         lhs: {}\n\
         rhs: {}\n\
         canonical lhs: {}\n\
         canonical rhs: {}\n\
         coefficient-wise polynomial sum: {}\n\
         reduced result: {}",
        format_extension_field_element::<S>(left),
        format_extension_field_element::<S>(right),
        format_extension_field_element::<S>(&left_reduced),
        format_extension_field_element::<S>(&right_reduced),
        format_extension_polynomial_addition::<S>(&left_reduced, &right_reduced),
        format_extension_field_element::<S>(&result)
    )
}

/// Explains multiplication inside the quotient.
pub fn explain_extension_field_mul<S>(
    left: &ExtensionFieldElement<S>,
    right: &ExtensionFieldElement<S>,
) -> String
where
    S: ExtensionFieldSpec,
    BaseElem<S>: VisualizableField,
{
    let left_reduced = ExtensionField::<S>::reduce(left);
    let right_reduced = ExtensionField::<S>::reduce(right);
    let result = ExtensionField::<S>::mul(left, right);

    format!(
        "Multiplication in an extension-field quotient\n\
         lhs: {}\n\
         rhs: {}\n\
         canonical lhs: {}\n\
         canonical rhs: {}\n\
         raw polynomial product: {}\n\
         reduced result: {}\n\
         note: multiplication happens in the polynomial ring over the base field and is then reduced modulo the defining polynomial",
        format_extension_field_element::<S>(left),
        format_extension_field_element::<S>(right),
        format_extension_field_element::<S>(&left_reduced),
        format_extension_field_element::<S>(&right_reduced),
        format_extension_polynomial_product::<S>(&left_reduced, &right_reduced),
        format_extension_field_element::<S>(&result)
    )
}

/// Explains multiplicative inversion inside the quotient.
pub fn explain_extension_field_inverse<S>(
    element: &ExtensionFieldElement<S>,
) -> Result<String, FieldError>
where
    S: ExtensionFieldSpec,
    BaseElem<S>: VisualizableField,
{
    let reduced = ExtensionField::<S>::reduce(element);
    let inverse = ExtensionField::<S>::inverse(element)?;
    let check = ExtensionField::<S>::mul(&reduced, &inverse);

    Ok(format!(
        "Inverse in an extension-field quotient\n\
         element: {}\n\
         canonical representative: {}\n\
         inverse: {}\n\
         verification: {} * {} = {}\n\
         note: the current backend computes inverses through the extended Euclidean algorithm in the polynomial ring over the base field",
        format_extension_field_element::<S>(element),
        format_extension_field_element::<S>(&reduced),
        format_extension_field_element::<S>(&inverse),
        format_extension_polynomial::<S>(reduced.coefficients()),
        format_extension_polynomial::<S>(inverse.coefficients()),
        format_extension_field_element::<S>(&check)
    ))
}

impl<S> Visualizable for ExtensionFieldElement<S>
where
    S: ExtensionFieldSpec,
    BaseElem<S>: VisualizableField,
{
    fn format_compact(&self) -> String {
        format_extension_field_element::<S>(self)
    }

    fn describe(&self) -> String {
        describe_extension_field_element::<S>(self)
    }
}

impl<S> core::fmt::Display for ExtensionFieldElement<S>
where
    S: ExtensionFieldSpec,
    BaseElem<S>: VisualizableField,
{
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            formatter,
            "{}",
            format_extension_field_element_compact::<S>(self)
        )
    }
}

impl<S> VisualizableField for ExtensionFieldElement<S>
where
    S: ExtensionFieldSpec,
    BaseElem<S>: VisualizableField,
{
    fn format_elem(&self) -> String {
        format_extension_field_element_compact::<S>(self)
    }

    fn inverse(&self) -> Option<String> {
        ExtensionField::<S>::inv(self).map(|value| format_extension_field_element::<S>(&value))
    }

    fn explain_add(lhs: &Self, rhs: &Self) -> Option<String> {
        Some(explain_extension_field_add::<S>(lhs, rhs))
    }

    fn explain_mul(lhs: &Self, rhs: &Self) -> Option<String> {
        Some(explain_extension_field_mul::<S>(lhs, rhs))
    }

    fn explain_div(lhs: &Self, rhs: &Self) -> Option<String> {
        let reciprocal = ExtensionField::<S>::inv(rhs)?;
        let result = ExtensionField::<S>::mul(lhs, &reciprocal);

        Some(format!(
            "Division in an extension-field quotient\n\
             lhs: {}\n\
             rhs: {}\n\
             reciprocal of rhs: {}\n\
             reduction to multiplication: {} * {} = {}",
            format_extension_field_element::<S>(lhs),
            format_extension_field_element::<S>(rhs),
            format_extension_field_element::<S>(&reciprocal),
            format_extension_field_element::<S>(lhs),
            format_extension_field_element::<S>(&reciprocal),
            format_extension_field_element::<S>(&result)
        ))
    }
}

fn format_extension_field_element_compact<S>(element: &ExtensionFieldElement<S>) -> String
where
    S: ExtensionFieldSpec,
    BaseElem<S>: VisualizableField,
{
    format_extension_polynomial_with_generator::<S>(element.coefficients(), "α")
}

fn format_extension_polynomial<S>(coefficients: &[BaseElem<S>]) -> String
where
    S: ExtensionFieldSpec,
    BaseElem<S>: VisualizableField,
{
    format_extension_polynomial_with_generator::<S>(coefficients, "x")
}

fn format_extension_polynomial_with_generator<S>(
    coefficients: &[BaseElem<S>],
    generator: &str,
) -> String
where
    S: ExtensionFieldSpec,
    BaseElem<S>: VisualizableField,
{
    let mut terms = Vec::new();

    for (power, coefficient) in coefficients.iter().enumerate().rev() {
        if <S::Base as Field>::is_zero(coefficient) {
            continue;
        }

        let coefficient_text = coefficient.format_elem();
        let term = match power {
            0 => parenthesize_if_needed(&coefficient_text),
            1 if <S::Base as Field>::eq(coefficient, &<S::Base as Field>::one()) => {
                generator.to_string()
            }
            1 => format!("{}*{generator}", parenthesize_if_needed(&coefficient_text)),
            _ if <S::Base as Field>::eq(coefficient, &<S::Base as Field>::one()) => {
                format!("{generator}^{power}")
            }
            _ => format!(
                "{}*{generator}^{power}",
                parenthesize_if_needed(&coefficient_text)
            ),
        };
        terms.push(term);
    }

    if terms.is_empty() {
        "0".to_string()
    } else {
        terms.join(" + ")
    }
}

fn format_extension_polynomial_addition<S>(
    left: &ExtensionFieldElement<S>,
    right: &ExtensionFieldElement<S>,
) -> String
where
    S: ExtensionFieldSpec,
    BaseElem<S>: VisualizableField,
{
    format!(
        "({}) + ({})",
        format_extension_polynomial::<S>(left.coefficients()),
        format_extension_polynomial::<S>(right.coefficients())
    )
}

fn format_extension_polynomial_product<S>(
    left: &ExtensionFieldElement<S>,
    right: &ExtensionFieldElement<S>,
) -> String
where
    S: ExtensionFieldSpec,
    BaseElem<S>: VisualizableField,
{
    format!(
        "({}) * ({})",
        format_extension_polynomial::<S>(left.coefficients()),
        format_extension_polynomial::<S>(right.coefficients())
    )
}

fn parenthesize_if_needed(text: &str) -> String {
    if text.contains(' ') || text.starts_with('-') || text.contains('/') {
        format!("({text})")
    } else {
        text.to_string()
    }
}

type BaseElem<S> = <<S as ExtensionFieldSpec>::Base as Field>::Elem;

#[cfg(test)]
mod tests {
    use crate::fields::traits::*;

    use crate::fields::extension_field::{
        ExtensionField, ExtensionFieldElement, ExtensionFieldSpec,
    };
    use crate::fields::{Q, polynomial_field::PolynomialModulus};
    use crate::visualization::fields::{
        describe_extension_field, describe_extension_field_element, explain_extension_field_add,
        explain_extension_field_inverse, explain_extension_field_mul,
        explain_extension_field_reduction, format_extension_field, format_extension_field_element,
    };
    use crate::visualization::{Visualizable, VisualizableField};

    crate::fields::extension_field::define_q_quadratic_extension!(
        spec: QSqrt2Spec,
        field: QSqrt2,
        radicand: 2,
        name: "Q(sqrt(2))",
    );

    struct QSqrt2ISpec;

    impl ExtensionFieldSpec for QSqrt2ISpec {
        type Base = QSqrt2;

        const NAME: &'static str = "Q(sqrt(2), i)";

        fn defining_modulus() -> PolynomialModulus<Self::Base> {
            PolynomialModulus::<QSqrt2>::new(vec![QSqrt2::one(), QSqrt2::zero(), QSqrt2::one()])
                .expect("x^2 + 1 should be structurally valid")
        }

        fn check_field_conditions() -> Result<(), crate::fields::FieldError> {
            Ok(())
        }
    }

    type QSqrt2I = ExtensionField<QSqrt2ISpec>;

    #[test]
    fn extension_field_summary_mentions_static_quotient_structure() {
        let formatted = format_extension_field::<QSqrt2Spec>();
        let description = describe_extension_field::<QSqrt2Spec>();

        assert!(formatted.contains("Q(sqrt(2))"));
        assert!(formatted.contains("x^2"));
        assert!(description.contains("degree over immediate base field: 2"));
        assert!(description.contains("reduced modulo the defining polynomial"));
    }

    #[test]
    fn extension_field_element_format_and_description_are_useful() {
        let element = QSqrt2::element(vec![Q::from_i64(3), Q::one()]);

        let formatted = format_extension_field_element::<QSqrt2Spec>(&element);
        let description = describe_extension_field_element::<QSqrt2Spec>(&element);

        assert!(formatted.contains("[x + 3]"));
        assert!(formatted.contains("mod (x^2 + (-2))"));
        assert!(description.contains("raw representative: x + 3"));
        assert!(description.contains("reduced degree: 1"));
    }

    #[test]
    fn extension_field_visualizable_field_surface_uses_compact_alpha_notation() {
        let element = QSqrt2::element(vec![Q::one(), Q::one()]);

        assert_eq!(element.format_elem(), "α + 1");
    }

    #[test]
    fn extension_field_reduction_explanation_handles_q_sqrt2_relation() {
        let element =
            ExtensionFieldElement::<QSqrt2Spec>::new(vec![Q::zero(), Q::zero(), Q::one()]);

        let explanation = explain_extension_field_reduction::<QSqrt2Spec>(&element);
        assert!(explanation.contains("raw representative: x^2"));
        assert!(explanation.contains("reduced representative: [2] mod (x^2 + (-2))"));
    }

    #[test]
    fn extension_field_addition_and_multiplication_explanations_are_readable() {
        let left = QSqrt2::element(vec![Q::one(), Q::one()]);
        let right = QSqrt2::element(vec![Q::from_i64(3), Q::from_i64(-1)]);

        let add = explain_extension_field_add::<QSqrt2Spec>(&left, &right);
        let mul = explain_extension_field_mul::<QSqrt2Spec>(&left, &right);

        assert!(add.contains("Addition in an extension-field quotient"));
        assert!(add.contains("reduced result: [4] mod (x^2 + (-2))"));

        assert!(mul.contains("Multiplication in an extension-field quotient"));
        assert!(mul.contains("raw polynomial product:"));
        assert!(mul.contains("reduced result: [2*x + 1] mod (x^2 + (-2))"));
    }

    #[test]
    fn extension_field_inverse_explanation_is_exact_for_q_sqrt2() {
        let element = QSqrt2::element(vec![Q::one(), Q::one()]);
        let explanation =
            explain_extension_field_inverse::<QSqrt2Spec>(&element).expect("inverse should exist");

        assert!(explanation.contains("Inverse in an extension-field quotient"));
        assert!(explanation.contains("inverse: [x + (-1)] mod (x^2 + (-2))"));
        assert!(explanation.contains("verification:"));
    }

    #[test]
    fn extension_field_elements_are_visualizable_and_work_inside_towers() {
        let i = QSqrt2I::element(vec![QSqrt2::zero(), QSqrt2::one()]);

        assert_eq!(i.format_compact(), "[x] mod (x^2 + 1)");
        assert!(i.to_string().contains("α"));
        assert!(i.describe().contains("ambient field: Q(sqrt(2), i)"));
        assert!(i.inverse().is_some());

        let add = ExtensionFieldElement::<QSqrt2ISpec>::explain_add(&i, &QSqrt2I::one())
            .expect("tower addition should be explainable");
        assert!(add.contains("Addition in an extension-field quotient"));
    }
}
