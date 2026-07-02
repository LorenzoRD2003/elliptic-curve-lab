#![cfg_attr(not(test), allow(dead_code))]

mod bezout;
mod input;
mod outcome;
mod report;
mod step;
mod trace;

#[cfg(test)]
mod tests;

pub(crate) use bezout::CyclicPrimeRootBezout;
pub(crate) use input::CyclicPrimeRootInput;
#[cfg(test)]
pub(crate) use input::CyclicPrimeRootInputError;
pub(crate) use outcome::CyclicPrimeRootOutcome;
#[cfg(test)]
pub(crate) use report::CyclicPrimeRootReport;
pub(crate) use step::CyclicPrimeRootStep;
pub(crate) use trace::CyclicPrimeRootTrace;
