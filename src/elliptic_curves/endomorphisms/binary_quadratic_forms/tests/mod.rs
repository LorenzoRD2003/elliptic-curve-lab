use num_bigint::BigInt;

mod cayley_table;
mod class_group_enumeration;
mod concordant_composition;
mod form_basics;
mod gp_tables;
mod group_laws;
mod reduction;

fn z(value: i64) -> BigInt {
    BigInt::from(value)
}
