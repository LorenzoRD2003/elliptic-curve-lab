//! Integer-root recovery from simple p-adic Hensel lifts.

mod config;
mod lift;
mod report;
mod search;

pub(crate) use config::HenselIntegerRootSearchConfig;
pub(crate) use lift::hensel_lift_integer_root;
pub(crate) use report::{HenselIntegerRootSearchReport, HenselIntegerRootTrace};
pub(crate) use search::find_integer_roots_by_hensel;
