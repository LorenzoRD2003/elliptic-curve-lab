mod element;
mod field;
mod macros;
mod spec;
mod traits;

#[cfg(test)]
mod tests;

use crate::fields::Field;
use crate::polynomials::DensePolynomial;

pub use element::ExtensionFieldElement;
pub use field::ExtensionField;
pub use spec::ExtensionFieldSpec;

type BaseElem<S> = <<S as ExtensionFieldSpec>::Base as Field>::Elem;
type DenseTriple<S> = (
    DensePolynomial<<S as ExtensionFieldSpec>::Base>,
    DensePolynomial<<S as ExtensionFieldSpec>::Base>,
    DensePolynomial<<S as ExtensionFieldSpec>::Base>,
);
