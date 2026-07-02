use crate::elliptic_curves::frobenius::torsion::matrix::FrobeniusTorsionMatrixError;
use num_bigint::BigInt;
use num_traits::ToPrimitive;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ModNMatrix2 {
    modulus: usize,
    entries: [[usize; 2]; 2],
}

impl ModNMatrix2 {
    pub fn from_columns(
        modulus: usize,
        first_column: [usize; 2],
        second_column: [usize; 2],
    ) -> Result<Self, FrobeniusTorsionMatrixError> {
        Self::new(
            modulus,
            [
                [first_column[0], second_column[0]],
                [first_column[1], second_column[1]],
            ],
        )
    }

    pub(crate) fn new(
        modulus: usize,
        entries: [[usize; 2]; 2],
    ) -> Result<Self, FrobeniusTorsionMatrixError> {
        for row in entries {
            for entry in row {
                if entry >= modulus {
                    return Err(FrobeniusTorsionMatrixError::CoefficientOutOfRange {
                        value: entry,
                        modulus,
                    });
                }
            }
        }
        Ok(Self { modulus, entries })
    }

    pub fn modulus(&self) -> usize {
        self.modulus
    }

    pub fn entries(&self) -> [[usize; 2]; 2] {
        self.entries
    }

    pub fn trace_mod_n(&self) -> usize {
        (self.entries[0][0] + self.entries[1][1]) % self.modulus
    }

    pub fn determinant_mod_n(&self) -> Result<usize, FrobeniusTorsionMatrixError> {
        let modulus = BigInt::from(self.modulus);
        let a11 = BigInt::from(self.entries[0][0]);
        let a12 = BigInt::from(self.entries[0][1]);
        let a21 = BigInt::from(self.entries[1][0]);
        let a22 = BigInt::from(self.entries[1][1]);
        let determinant = a11 * a22 - a12 * a21;
        let reduced = ((determinant % &modulus) + &modulus) % &modulus;
        reduced
            .to_usize()
            .ok_or(FrobeniusTorsionMatrixError::DeterminantOverflow {
                modulus: self.modulus,
            })
    }
}
