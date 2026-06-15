use crate::fields::{
    rational_function_field::RationalFunction, traits::Field, traits::FiniteField,
    traits::PthRootExtraction,
};
use crate::polynomials::{DensePolynomial, SparsePolynomial};
use crate::visualization::fields::format_rational_function;
use crate::visualization::fields::traits::VisualizableField;
use crate::visualization::polynomials::{format_dense_polynomial, format_sparse_polynomial};

/// Explains `p`-th-root extraction for one finite-field element.
pub fn explain_finite_field_pth_root<F>(element: &F::Elem) -> String
where
    F: FiniteField,
    F::Elem: PthRootExtraction + VisualizableField,
{
    let root = element
        .pth_root()
        .expect("finite-field elements always admit a p-th root");

    format!(
        "Characteristic-p root extraction in a finite field\n\
         characteristic: {}\n\
         element: {}\n\
         root: {}\n\
         note: in F_(p^n), Frobenius x -> x^p is an automorphism, so every element has a unique p-th root",
        F::characteristic(),
        element.format_elem(),
        root.format_elem(),
    )
}

/// Explains `p`-th-root extraction for a dense polynomial over a finite field.
pub fn explain_dense_polynomial_pth_root<F>(polynomial: &DensePolynomial<F>) -> String
where
    F: FiniteField,
    F::Elem: PthRootExtraction + VisualizableField,
{
    let offending_degrees = dense_offending_degrees::<F>(polynomial);

    match polynomial.pth_root() {
        Some(root) => format!(
            "Characteristic-p root extraction for a dense polynomial\n\
             characteristic: {}\n\
             polynomial: {}\n\
             root: {}\n\
             criterion: every non-zero term degree is divisible by p, and each surviving coefficient admits a p-th root",
            F::characteristic(),
            format_dense_polynomial(polynomial),
            format_dense_polynomial(&root),
        ),
        None => format!(
            "Characteristic-p root extraction for a dense polynomial\n\
             characteristic: {}\n\
             polynomial: {}\n\
             p-th root: does not exist in F[x]\n\
             failing degrees: {}\n\
             criterion: a non-zero term a_i*x^i can belong to a p-th power only when i is divisible by p",
            F::characteristic(),
            format_dense_polynomial(polynomial),
            format_degree_list(&offending_degrees),
        ),
    }
}

/// Explains `p`-th-root extraction for a sparse polynomial over a finite field.
pub fn explain_sparse_polynomial_pth_root<F>(polynomial: &SparsePolynomial<F>) -> String
where
    F: FiniteField,
    F::Elem: PthRootExtraction + VisualizableField,
{
    let offending_degrees = sparse_offending_degrees::<F>(polynomial);

    match polynomial.pth_root() {
        Some(root) => format!(
            "Characteristic-p root extraction for a sparse polynomial\n\
             characteristic: {}\n\
             polynomial: {}\n\
             root: {}\n\
             criterion: every stored non-zero term degree is divisible by p, and each stored coefficient admits a p-th root",
            F::characteristic(),
            format_sparse_polynomial(polynomial),
            format_sparse_polynomial(&root),
        ),
        None => format!(
            "Characteristic-p root extraction for a sparse polynomial\n\
             characteristic: {}\n\
             polynomial: {}\n\
             p-th root: does not exist in F[x]\n\
             failing degrees: {}\n\
             criterion: a stored term a_i*x^i can belong to a p-th power only when i is divisible by p",
            F::characteristic(),
            format_sparse_polynomial(polynomial),
            format_degree_list(&offending_degrees),
        ),
    }
}

/// Explains `p`-th-root extraction for a rational function over a finite field.
pub fn explain_rational_function_pth_root<F>(function: &RationalFunction<F>) -> String
where
    F: FiniteField,
    F::Elem: PthRootExtraction + VisualizableField,
{
    match function.pth_root() {
        Some(root) => format!(
            "Characteristic-p root extraction in F(x)\n\
             characteristic: {}\n\
             function: {}\n\
             root: {}\n\
             criterion: in the canonical reduced presentation P(x)/Q(x), both P and Q must admit p-th roots in F[x]",
            F::characteristic(),
            format_rational_function(function),
            format_rational_function(&root),
        ),
        None => format!(
            "Characteristic-p root extraction in F(x)\n\
             characteristic: {}\n\
             function: {}\n\
             p-th root: does not exist in F(x)\n\
             criterion: in the canonical reduced presentation P(x)/Q(x), both P and Q must admit p-th roots in F[x]",
            F::characteristic(),
            format_rational_function(function),
        ),
    }
}

fn dense_offending_degrees<F: Field>(polynomial: &DensePolynomial<F>) -> Vec<usize> {
    let characteristic = F::characteristic();
    polynomial
        .coefficients()
        .iter()
        .enumerate()
        .filter_map(|(degree, coefficient)| {
            if !F::is_zero(coefficient) && (u64::try_from(degree).ok()? % characteristic != 0) {
                Some(degree)
            } else {
                None
            }
        })
        .collect()
}

fn sparse_offending_degrees<F: Field>(polynomial: &SparsePolynomial<F>) -> Vec<usize> {
    let characteristic = F::characteristic();
    polynomial
        .terms()
        .iter()
        .filter_map(|term| {
            if u64::try_from(term.degree).ok()? % characteristic != 0 {
                Some(term.degree)
            } else {
                None
            }
        })
        .collect()
}

fn format_degree_list(degrees: &[usize]) -> String {
    if degrees.is_empty() {
        "none".to_string()
    } else {
        degrees
            .iter()
            .map(usize::to_string)
            .collect::<Vec<_>>()
            .join(", ")
    }
}

#[cfg(test)]
mod tests {
    use crate::fields::{Fp, rational_function_field::RationalFunction, traits::Field};
    use crate::polynomials::sparse::SparsePolynomialTerm;
    use crate::polynomials::{DensePolynomial, SparsePolynomial};
    use crate::visualization::fields::{
        explain_dense_polynomial_pth_root, explain_finite_field_pth_root,
        explain_rational_function_pth_root, explain_sparse_polynomial_pth_root,
    };

    type F17 = Fp<17>;

    fn f17_dense(values: &[u64]) -> DensePolynomial<F17> {
        DensePolynomial::<F17>::new(values.iter().copied().map(F17::elem_from_u64).collect())
    }

    fn f17_sparse_term(coefficient: u64, degree: usize) -> SparsePolynomialTerm<F17> {
        SparsePolynomialTerm {
            coefficient: F17::elem_from_u64(coefficient),
            degree,
        }
    }

    #[test]
    fn finite_field_pth_root_explanation_mentions_frobenius_story() {
        let text = explain_finite_field_pth_root::<F17>(&F17::elem_from_u64(5));

        assert!(text.contains("Characteristic-p root extraction"));
        assert!(text.contains("Frobenius"));
        assert!(text.contains("unique p-th root"));
    }

    #[test]
    fn dense_polynomial_pth_root_explanation_mentions_success_and_failure_criteria() {
        let success = explain_dense_polynomial_pth_root(&DensePolynomial::<F17>::new({
            let mut coefficients = vec![F17::zero(); 18];
            coefficients[17] = F17::one();
            coefficients
        }));
        let failure = explain_dense_polynomial_pth_root(&f17_dense(&[0, 1]));

        assert!(success.contains("root: x"));
        assert!(success.contains("divisible by p"));
        assert!(failure.contains("does not exist"));
        assert!(failure.contains("failing degrees: 1"));
    }

    #[test]
    fn sparse_polynomial_pth_root_explanation_mentions_success_and_failure_criteria() {
        let success = explain_sparse_polynomial_pth_root(&SparsePolynomial::<F17>::new(vec![
            f17_sparse_term(1, 17),
        ]));
        let failure = explain_sparse_polynomial_pth_root(&SparsePolynomial::<F17>::new(vec![
            f17_sparse_term(1, 1),
        ]));

        assert!(success.contains("root: x"));
        assert!(failure.contains("does not exist"));
        assert!(failure.contains("failing degrees: 1"));
    }

    #[test]
    fn rational_function_pth_root_explanation_mentions_reduced_presentation() {
        let success = explain_rational_function_pth_root(
            &RationalFunction::<F17>::new(
                DensePolynomial::<F17>::new({
                    let mut coefficients = vec![F17::zero(); 18];
                    coefficients[17] = F17::one();
                    coefficients
                }),
                DensePolynomial::<F17>::new({
                    let mut coefficients = vec![F17::zero(); 18];
                    coefficients[0] = F17::one();
                    coefficients[17] = F17::one();
                    coefficients
                }),
            )
            .expect("example function should exist"),
        );
        let failure = explain_rational_function_pth_root(
            &RationalFunction::<F17>::new(f17_dense(&[0, 1]), f17_dense(&[1]))
                .expect("x should exist"),
        );

        assert!(success.contains("root: (x) / (x + 1)"));
        assert!(success.contains("canonical reduced presentation"));
        assert!(failure.contains("does not exist"));
    }
}
