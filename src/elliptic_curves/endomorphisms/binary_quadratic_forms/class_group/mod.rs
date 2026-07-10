mod cayley_table;
mod concordant;
mod dirichlet;
mod enumeration;
mod group;
mod inverse;
mod membership;

#[cfg(test)]
pub(crate) mod equivalence;

pub use cayley_table::QuadraticClassGroupCayleyTable;
pub use group::QuadraticClassGroup;
