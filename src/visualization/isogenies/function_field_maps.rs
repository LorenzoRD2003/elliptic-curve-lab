use core::fmt;

use crate::elliptic_curves::{ShortWeierstrassFunction, ShortWeierstrassFunctionField};
use crate::fields::{Field, RationalFunction};
use crate::isogenies::ShortWeierstrassFunctionFieldMap;
use crate::polynomials::DensePolynomial;
use crate::visualization::elliptic_curves::{
    describe_short_weierstrass_function_field, format_curve, format_short_weierstrass_function,
};
use crate::visualization::fields::format_rational_function;
use crate::visualization::fields::traits::VisualizableField;
use crate::visualization::traits::Visualizable;

/// Formats one short-Weierstrass function-field pullback map compactly.
pub fn format_short_weierstrass_function_field_map<F>(
    map: &ShortWeierstrassFunctionFieldMap<F>,
) -> String
where
    F: Field,
    F::Elem: VisualizableField + PartialEq,
{
    format!(
        "x' -> {}, y' -> {}",
        format_short_weierstrass_function(map.x_pullback()),
        format_short_weierstrass_function(map.y_pullback())
    )
}

/// Returns a richer educational description of one function-field pullback map.
pub fn describe_short_weierstrass_function_field_map<F>(
    map: &ShortWeierstrassFunctionFieldMap<F>,
) -> String
where
    F: Field,
    F::Elem: VisualizableField + fmt::Display + PartialEq,
{
    format!(
        "Short-Weierstrass function-field pullback\n\
         domain curve E: {}\n\
         codomain curve E': {}\n\
         induced direction: F(E') -> F(E)\n\
         stored generators:\n\
         - phi*(x') = {}\n\
         - phi*(y') = {}\n\
         compact form: {}\n\
         note: this map is determined by the images of x' and y' in the basis 1, y over F(x), and the constructor already checked that these images satisfy the codomain equation after substitution",
        format_curve(map.domain_curve()),
        format_curve(map.codomain_curve()),
        format_short_weierstrass_function(map.x_pullback()),
        format_short_weierstrass_function(map.y_pullback()),
        format_short_weierstrass_function_field_map(map)
    )
}

/// Explains the pullback of a polynomial in the codomain `x'`-coordinate.
pub fn explain_short_weierstrass_function_field_map_pullback_polynomial<F>(
    map: &ShortWeierstrassFunctionFieldMap<F>,
    polynomial: &DensePolynomial<F>,
) -> Result<String, crate::isogenies::IsogenyError>
where
    F: Field,
    F::Elem: VisualizableField + fmt::Display + PartialEq,
{
    let result = map.pullback_polynomial(polynomial)?;

    Ok(format!(
        "Polynomial pullback through phi*\n\
         codomain polynomial p(x'): {}\n\
         substitution rule: x' -> {}\n\
         computed pullback p(phi*(x')): {}",
        crate::visualization::polynomials::format_dense_polynomial(polynomial),
        format_short_weierstrass_function(map.x_pullback()),
        format_short_weierstrass_function(&result)
    ))
}

/// Explains the pullback of a rational function in the codomain `x'`-coordinate.
pub fn explain_short_weierstrass_function_field_map_pullback_rational_function<F>(
    map: &ShortWeierstrassFunctionFieldMap<F>,
    function: &RationalFunction<F>,
) -> Result<String, crate::isogenies::IsogenyError>
where
    F: Field,
    F::Elem: VisualizableField + fmt::Display + PartialEq,
{
    let result = map.pullback_rational_function(function)?;

    Ok(format!(
        "Rational-function pullback through phi*\n\
         codomain rational function r(x'): {}\n\
         substitution rule: x' -> {}\n\
         computed pullback r(phi*(x')): {}",
        format_rational_function(function),
        format_short_weierstrass_function(map.x_pullback()),
        format_short_weierstrass_function(&result)
    ))
}

/// Explains the pullback of a full function-field element `A(x') + y'B(x')`.
pub fn explain_short_weierstrass_function_field_map_pullback_function<F>(
    map: &ShortWeierstrassFunctionFieldMap<F>,
    function: &ShortWeierstrassFunction<F>,
) -> Result<String, crate::isogenies::IsogenyError>
where
    F: Field,
    F::Elem: VisualizableField + fmt::Display + PartialEq,
{
    let result = map.pullback_function(function)?;

    Ok(format!(
        "Function-field pullback through phi*\n\
         codomain function: {}\n\
         basis rule: phi*(A(x') + y'B(x')) = A(phi*(x')) + phi*(y')*B(phi*(x'))\n\
         phi*(x'): {}\n\
         phi*(y'): {}\n\
         computed pullback: {}",
        format_short_weierstrass_function(function),
        format_short_weierstrass_function(map.x_pullback()),
        format_short_weierstrass_function(map.y_pullback()),
        format_short_weierstrass_function(&result)
    ))
}

/// Explains the contravariant composition of two pullback maps.
pub fn explain_short_weierstrass_function_field_map_composition<F>(
    first: &ShortWeierstrassFunctionFieldMap<F>,
    second: &ShortWeierstrassFunctionFieldMap<F>,
) -> Result<String, crate::isogenies::IsogenyError>
where
    F: Field,
    F::Elem: VisualizableField + fmt::Display + PartialEq,
{
    let composite = first.compose(second)?;

    Ok(format!(
        "Composition of function-field pullbacks\n\
         first map: {}\n\
         second map: {}\n\
         contravariant rule: (psi o phi)* = phi* o psi*\n\
         middle curve agreement: codomain(first) = domain(second)\n\
         composite generators:\n\
         - x'' -> {}\n\
         - y'' -> {}\n\
         compact composite: {}",
        format_short_weierstrass_function_field_map(first),
        format_short_weierstrass_function_field_map(second),
        format_short_weierstrass_function(composite.x_pullback()),
        format_short_weierstrass_function(composite.y_pullback()),
        format_short_weierstrass_function_field_map(&composite)
    ))
}

/// Returns a compact description of the ambient fields attached to one pullback map.
pub fn describe_short_weierstrass_function_field_map_ambient_fields<F>(
    map: &ShortWeierstrassFunctionFieldMap<F>,
) -> String
where
    F: Field,
    F::Elem: VisualizableField + fmt::Display + PartialEq,
{
    let domain_field = ShortWeierstrassFunctionField::<F>::new(map.domain_curve().clone());
    let codomain_field = ShortWeierstrassFunctionField::<F>::new(map.codomain_curve().clone());

    format!(
        "Ambient fields around phi*\n\
         codomain source:\n{}\n\
         \n\
         domain target:\n{}",
        describe_short_weierstrass_function_field(&codomain_field),
        describe_short_weierstrass_function_field(&domain_field)
    )
}

impl<F> Visualizable for ShortWeierstrassFunctionFieldMap<F>
where
    F: Field,
    F::Elem: VisualizableField + fmt::Display + PartialEq,
{
    fn format_compact(&self) -> String {
        format_short_weierstrass_function_field_map(self)
    }

    fn describe(&self) -> String {
        describe_short_weierstrass_function_field_map(self)
    }
}

#[cfg(test)]
mod tests {
    use crate::elliptic_curves::{ShortWeierstrassCurve, ShortWeierstrassFunctionField};
    use crate::fields::{Field, Fp, RationalFunction};
    use crate::isogenies::ShortWeierstrassFunctionFieldMap;
    use crate::polynomials::DensePolynomial;
    use crate::visualization::isogenies::{
        describe_short_weierstrass_function_field_map,
        describe_short_weierstrass_function_field_map_ambient_fields,
        explain_short_weierstrass_function_field_map_composition,
        explain_short_weierstrass_function_field_map_pullback_function,
        explain_short_weierstrass_function_field_map_pullback_polynomial,
        explain_short_weierstrass_function_field_map_pullback_rational_function,
        format_short_weierstrass_function_field_map,
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

    fn identity_map() -> ShortWeierstrassFunctionFieldMap<F17> {
        let field = ShortWeierstrassFunctionField::<F17>::new(curve());
        ShortWeierstrassFunctionFieldMap::new(curve(), curve(), field.x(), field.y())
            .expect("identity map should validate")
    }

    fn negation_map() -> ShortWeierstrassFunctionFieldMap<F17> {
        let field = ShortWeierstrassFunctionField::<F17>::new(curve());
        ShortWeierstrassFunctionFieldMap::new(curve(), curve(), field.x(), field.y().neg())
            .expect("negation map should validate")
    }

    #[test]
    fn formatter_and_description_mention_generators_and_direction() {
        let map = negation_map();
        let compact = format_short_weierstrass_function_field_map(&map);
        let description = describe_short_weierstrass_function_field_map(&map);

        assert!(compact.contains("x' -> x"));
        assert!(compact.contains("y' ->"));
        assert!(description.contains("Short-Weierstrass function-field pullback"));
        assert!(description.contains("F(E') -> F(E)"));
        assert!(description.contains("phi*(x')"));
        assert!(description.contains("phi*(y')"));
    }

    #[test]
    fn polynomial_and_rational_pullback_explanations_show_substitution_rule() {
        let map = identity_map();
        let polynomial = f17_dense(&[3, 2, 1]);
        let rational = RationalFunction::<F17>::new(f17_dense(&[1, 0, 1]), f17_dense(&[1, 1]))
            .expect("denominator should be non-zero");

        let polynomial_text =
            explain_short_weierstrass_function_field_map_pullback_polynomial(&map, &polynomial)
                .expect("pullback should work");
        let rational_text =
            explain_short_weierstrass_function_field_map_pullback_rational_function(
                &map, &rational,
            )
            .expect("pullback should work");

        assert!(polynomial_text.contains("Polynomial pullback through phi*"));
        assert!(polynomial_text.contains("substitution rule: x' -> x"));
        assert!(rational_text.contains("Rational-function pullback through phi*"));
        assert!(rational_text.contains("r(phi*(x'))"));
    }

    #[test]
    fn full_function_pullback_and_composition_explanations_show_basis_and_contravariance() {
        let field = ShortWeierstrassFunctionField::<F17>::new(curve());
        let function = field
            .x()
            .add(&field.y())
            .expect("same-curve addition should work");
        let pullback_text = explain_short_weierstrass_function_field_map_pullback_function(
            &negation_map(),
            &function,
        )
        .expect("pullback should work");
        let composition_text = explain_short_weierstrass_function_field_map_composition(
            &negation_map(),
            &negation_map(),
        )
        .expect("composition should work");

        assert!(pullback_text.contains("A(x') + y'B(x')"));
        assert!(pullback_text.contains("phi*(y')"));
        assert!(composition_text.contains("contravariant rule: (psi o phi)* = phi* o psi*"));
        assert!(composition_text.contains("compact composite: x' -> x, y' -> y"));
    }

    #[test]
    fn ambient_field_description_reuses_function_field_summaries() {
        let text = describe_short_weierstrass_function_field_map_ambient_fields(&identity_map());

        assert!(text.contains("Ambient fields around phi*"));
        assert!(text.contains("Short-Weierstrass function field"));
        assert!(text.contains("codomain source"));
        assert!(text.contains("domain target"));
    }

    #[test]
    fn visualizable_trait_reuses_function_field_map_helpers() {
        let map = negation_map();

        assert!(map.format_compact().contains("x' -> x"));
        assert!(
            map.describe()
                .contains("Short-Weierstrass function-field pullback")
        );
    }
}
