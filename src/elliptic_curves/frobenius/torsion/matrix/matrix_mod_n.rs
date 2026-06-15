use crate::elliptic_curves::frobenius::torsion::matrix::FrobeniusTorsionMatrixError;

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
        let modulus = i128::try_from(self.modulus).map_err(|_| {
            FrobeniusTorsionMatrixError::DeterminantOverflow {
                modulus: self.modulus,
            }
        })?;
        let a11 = i128::try_from(self.entries[0][0]).unwrap();
        let a12 = i128::try_from(self.entries[0][1]).unwrap();
        let a21 = i128::try_from(self.entries[1][0]).unwrap();
        let a22 = i128::try_from(self.entries[1][1]).unwrap();
        let determinant = a11 * a22 - a12 * a21;
        let reduced = ((determinant % modulus) + modulus) % modulus;
        usize::try_from(reduced).map_err(|_| FrobeniusTorsionMatrixError::DeterminantOverflow {
            modulus: self.modulus,
        })
    }
}
