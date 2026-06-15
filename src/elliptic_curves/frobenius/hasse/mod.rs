#[cfg(test)]
mod benchmarks;
mod bound;
mod interval;
mod search_bsgs;
mod search_naive;
mod search_report;
#[cfg(test)]
mod tests;

pub use bound::{HasseBoundReport, verify_hasse_bound};
pub use interval::HasseInterval;
pub use search_report::{HasseMultipleSearchReport, HasseMultipleSearchStep};

pub(crate) use search_bsgs::{
    HasseBsgsConfig, HasseBsgsParity, find_annihilating_multiple_in_interval_bsgs,
    find_annihilating_multiple_in_interval_bsgs_with_config,
};
pub(crate) use search_naive::find_annihilating_multiple_in_interval_naive_report;
pub(crate) use search_report::hasse_multiple_search_report;
