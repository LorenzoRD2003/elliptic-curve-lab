mod curve_model;
mod report;

#[cfg(test)]
mod tests;

pub use curve_model::FrobeniusCharacteristicEquationCurveModel;
pub use report::{
    FrobeniusCharacteristicEquationCheck, FrobeniusCharacteristicEquationExhaustiveReport,
    FrobeniusCharacteristicEquationTerms,
};
