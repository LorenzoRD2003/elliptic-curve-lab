use std::collections::HashMap;

use num_bigint::BigInt;

use crate::elliptic_curves::endomorphisms::binary_quadratic_forms::{
    BinaryQuadraticForm, BinaryQuadraticFormError, QuadraticClassGroup,
};
use crate::elliptic_curves::endomorphisms::quadratic_orders::QuadraticDiscriminant;

/// Cayley table for a small enumerated quadratic-form class group.
///
/// The table stores one canonical payload:
///
/// - the group discriminant `D`;
/// - the reduced representatives in enumeration order;
/// - a product matrix whose entry `(i,j)` is the index of
///   `representatives[i] * representatives[j]`.
///
/// Everything else, such as row labels, inverses, or textual layout, should be
/// derived from this data so summaries cannot drift out of sync.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct QuadraticClassGroupCayleyTable {
    discriminant: QuadraticDiscriminant,
    representatives: Vec<BinaryQuadraticForm>,
    products: Vec<Vec<usize>>,
}

impl QuadraticClassGroup {
    /// Builds the Cayley table for the enumerated reduced representatives.
    ///
    /// This is an educational/exploratory helper for small class groups. It
    /// first calls [`Self::enumerate_reduced_forms`], then composes every pair
    /// with [`Self::compose`] and records the index of the reduced product.
    ///
    /// Complexity: `Θ(h(D))` index construction plus `Θ(h(D)²)` compositions
    /// and expected constant-time product-index lookups, where `h(D)` is the
    /// number of reduced representatives enumerated for the discriminant `D`.
    pub fn cayley_table(&self) -> Result<QuadraticClassGroupCayleyTable, BinaryQuadraticFormError> {
        let representatives = self.enumerate_reduced_forms();
        let representative_indices: HashMap<_, _> = representatives
            .iter()
            .enumerate()
            .map(|(index, representative)| (form_key(representative), index))
            .collect();
        let mut products = Vec::with_capacity(representatives.len());

        for left in &representatives {
            let mut row = Vec::with_capacity(representatives.len());
            for right in &representatives {
                let product = self.compose(left, right)?;
                let product_index = representative_indices
                    .get(&form_key(&product))
                    .copied()
                    .expect("class-group composition should return an enumerated representative");
                row.push(product_index);
            }
            products.push(row);
        }

        Ok(QuadraticClassGroupCayleyTable {
            discriminant: self.discriminant().clone(),
            representatives,
            products,
        })
    }
}

fn form_key(form: &BinaryQuadraticForm) -> (BigInt, BigInt, BigInt) {
    (form.a().clone(), form.b().clone(), form.c().clone())
}

impl QuadraticClassGroupCayleyTable {
    /// Returns the fixed negative quadratic-order discriminant `D`.
    pub fn discriminant(&self) -> &QuadraticDiscriminant {
        &self.discriminant
    }

    /// Returns the reduced representatives indexing the table rows and columns.
    pub fn representatives(&self) -> &[BinaryQuadraticForm] {
        &self.representatives
    }

    /// Returns the product-index matrix.
    pub fn products(&self) -> &[Vec<usize>] {
        &self.products
    }

    /// Returns the number of classes represented in the table.
    pub fn class_number(&self) -> usize {
        self.representatives.len()
    }

    /// Returns the representative at `index`, if present.
    pub fn representative(&self, index: usize) -> Option<&BinaryQuadraticForm> {
        self.representatives.get(index)
    }

    /// Returns the product index for row `left` and column `right`, if present.
    pub fn product_index(&self, left: usize, right: usize) -> Option<usize> {
        self.products
            .get(left)
            .and_then(|row| row.get(right))
            .copied()
    }
}
