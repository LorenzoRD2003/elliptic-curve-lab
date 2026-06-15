use crate::elliptic_curves::frobenius::{
    FrobeniusTrace,
    torsion::matrix::{FrobeniusTorsionMatrixError, ModNMatrix2, NTorsionBasis},
};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FrobeniusTorsionMatrixReport<P> {
    frobenius_trace: FrobeniusTrace,
    basis: NTorsionBasis<P>,
    matrix: ModNMatrix2,
    trace_matches_mod_n: bool,
    determinant_matches_mod_n: bool,
}

impl<P> FrobeniusTorsionMatrixReport<P> {
    pub(crate) fn from_matrix_and_trace(
        frobenius_trace: FrobeniusTrace,
        basis: NTorsionBasis<P>,
        matrix: ModNMatrix2,
    ) -> Result<Self, FrobeniusTorsionMatrixError> {
        let trace_matches_mod_n = frobenius_trace.trace_matches_torsion_matrix_mod_n(&matrix)?;
        let determinant_matches_mod_n =
            frobenius_trace.determinant_matches_torsion_matrix_mod_n(&matrix)?;

        Ok(Self {
            frobenius_trace,
            basis,
            matrix,
            trace_matches_mod_n,
            determinant_matches_mod_n,
        })
    }

    pub fn frobenius_trace(&self) -> &FrobeniusTrace {
        &self.frobenius_trace
    }

    pub fn basis(&self) -> &NTorsionBasis<P> {
        &self.basis
    }

    pub fn matrix(&self) -> &ModNMatrix2 {
        &self.matrix
    }

    pub fn trace_matches_mod_n(&self) -> bool {
        self.trace_matches_mod_n
    }

    pub fn determinant_matches_mod_n(&self) -> bool {
        self.determinant_matches_mod_n
    }
}
