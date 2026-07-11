use std::collections::BTreeMap;

use num_bigint::{BigInt, BigUint};

use crate::elliptic_curves::endomorphisms::{
    binary_quadratic_forms::BinaryQuadraticForm, quadratic_orders::QuadraticDiscriminant,
};
use crate::isogenies::class_group_action::{
    ClassGroupActionPlanFactor, ClassGroupIsogenyActionError, OrientedLabeledCraterWalkReport,
};

/// Indexed collection of oriented local crater witnesses for one action plan.
///
/// The current public execution API accepts a slice of local witnesses. This
/// crate-private index validates that the slice is unambiguous for the target
/// discriminant, then provides logarithmic lookup by the local ideal/form label
/// used by [`ClassGroupActionPlanFactor`].
///
/// Complexity: construction costs `O(w log w)` comparisons for `w` witnesses.
/// Each factor lookup costs `O(log w)`.
#[derive(Clone, Debug)]
pub(crate) struct LocalCraterWitnessSet<'a> {
    witnesses: BTreeMap<LocalCraterWitnessKey, IndexedLocalCraterWitness<'a>>,
}

impl<'a> LocalCraterWitnessSet<'a> {
    pub(crate) fn new(
        discriminant: &QuadraticDiscriminant,
        local_witnesses: &'a [OrientedLabeledCraterWalkReport],
    ) -> Result<Self, ClassGroupIsogenyActionError> {
        let mut witnesses: BTreeMap<LocalCraterWitnessKey, IndexedLocalCraterWitness<'a>> =
            BTreeMap::new();

        for (witness_index, witness) in local_witnesses.iter().enumerate() {
            let key = LocalCraterWitnessKey::from_witness(discriminant, witness_index, witness)?;
            let generator_form = witness.labeled_walk().form_label().reduced_form().clone();
            match witnesses.get(&key) {
                Some(existing) if existing.witness.orientation() == witness.orientation() => {
                    return Err(ClassGroupIsogenyActionError::DuplicateLocalWitness {
                        first_witness_index: existing.index,
                        duplicate_witness_index: witness_index,
                        ideal_norm: key.ideal_norm.clone(),
                        generator_form,
                    });
                }
                Some(existing) => {
                    return Err(
                        ClassGroupIsogenyActionError::ConflictingLocalWitnessOrientation {
                            first_witness_index: existing.index,
                            conflicting_witness_index: witness_index,
                            ideal_norm: key.ideal_norm.clone(),
                            generator_form,
                        },
                    );
                }
                None => {
                    witnesses.insert(
                        key,
                        IndexedLocalCraterWitness {
                            index: witness_index,
                            witness,
                        },
                    );
                }
            }
        }

        Ok(Self { witnesses })
    }

    pub(crate) fn get(
        &self,
        factor_index: usize,
        factor: &ClassGroupActionPlanFactor,
    ) -> Result<&'a OrientedLabeledCraterWalkReport, ClassGroupIsogenyActionError> {
        let key = LocalCraterWitnessKey::from_factor(factor);
        self.witnesses
            .get(&key)
            .map(|indexed| indexed.witness)
            .ok_or_else(|| ClassGroupIsogenyActionError::MissingLocalWitness {
                factor_index,
                ideal_norm: factor.ideal().norm().clone(),
                generator_form: factor.generator_form().clone(),
            })
    }
}

#[derive(Clone, Debug)]
struct IndexedLocalCraterWitness<'a> {
    index: usize,
    witness: &'a OrientedLabeledCraterWalkReport,
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct LocalCraterWitnessKey {
    discriminant: BigInt,
    ideal_norm: BigUint,
    ideal_root_mod_ell: BigUint,
    generator_form_key: (BigInt, BigInt, BigInt),
}

impl LocalCraterWitnessKey {
    fn from_witness(
        plan_discriminant: &QuadraticDiscriminant,
        witness_index: usize,
        witness: &OrientedLabeledCraterWalkReport,
    ) -> Result<Self, ClassGroupIsogenyActionError> {
        let ideal = witness.labeled_walk().local_label().ideal();
        let witness_discriminant = ideal.order().discriminant();
        if witness_discriminant != plan_discriminant {
            return Err(
                ClassGroupIsogenyActionError::LocalWitnessDiscriminantMismatch {
                    witness_index,
                    witness_discriminant: witness_discriminant.value().clone(),
                    plan_discriminant: plan_discriminant.value().clone(),
                },
            );
        }

        Ok(Self {
            discriminant: witness_discriminant.value().clone(),
            ideal_norm: ideal.norm().clone(),
            ideal_root_mod_ell: ideal.root_mod_ell().clone(),
            generator_form_key: form_key(witness.labeled_walk().form_label().reduced_form()),
        })
    }

    fn from_factor(factor: &ClassGroupActionPlanFactor) -> Self {
        Self {
            discriminant: factor.ideal().order().discriminant().value().clone(),
            ideal_norm: factor.ideal().norm().clone(),
            ideal_root_mod_ell: factor.ideal().root_mod_ell().clone(),
            generator_form_key: form_key(factor.generator_form()),
        }
    }
}

fn form_key(form: &BinaryQuadraticForm) -> (BigInt, BigInt, BigInt) {
    let (a, b, c) = form.coefficients();
    (a.clone(), b.clone(), c.clone())
}
