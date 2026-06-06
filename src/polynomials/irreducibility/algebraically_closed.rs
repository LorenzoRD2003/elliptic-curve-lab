use crate::fields::{ComplexApprox, Field};
use crate::polynomials::{DensePolynomial, PolynomialError};

use crate::polynomials::irreducibility::{
    IrreducibilityBackend, IrreducibilityStatus, ReducibilityReason,
};

pub(super) fn classify_in_algebraically_closed_field<F: Field>(
    polynomial: &DensePolynomial<F>,
) -> IrreducibilityStatus<F> {
    debug_assert!(
        F::IS_ALGEBRAICALLY_CLOSED,
        "this helper should be used only for algebraically closed backends"
    );

    match polynomial.degree() {
        None | Some(0) => IrreducibilityStatus::Constant,
        Some(1) => IrreducibilityStatus::Linear,
        Some(_) => IrreducibilityStatus::ReducibleWithoutWitness {
            reason: ReducibilityReason::AlgebraicallyClosedField,
        },
    }
}

impl IrreducibilityBackend for ComplexApprox {
    fn irreducibility_status_impl(
        polynomial: &DensePolynomial<Self>,
    ) -> Result<IrreducibilityStatus<Self>, PolynomialError> {
        Ok(classify_in_algebraically_closed_field(polynomial))
    }
}
