mod matrix;
mod points_report;

#[cfg(test)]
mod tests;

pub use matrix::{
    FrobeniusTorsionMatrixError, FrobeniusTorsionMatrixReport, ModNMatrix2, NTorsionBasis,
};

pub(crate) use matrix::TorsionCoordinateTable;
pub(crate) use points_report::{FrobeniusOnExactTorsionPoint, FrobeniusOnExactTorsionReport};
