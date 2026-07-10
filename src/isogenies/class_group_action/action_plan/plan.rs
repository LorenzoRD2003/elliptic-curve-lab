use std::collections::VecDeque;

use crate::elliptic_curves::endomorphisms::{
    binary_quadratic_forms::{BinaryQuadraticForm, BinaryQuadraticFormError, QuadraticClassGroup},
    quadratic_ideals::{IdealFormCorrespondence, PrimeNormIdeal},
    quadratic_orders::QuadraticDiscriminant,
};
use crate::isogenies::class_group_action::{ClassGroupActionPlanError, ClassGroupActionPlanFactor};

/// Algebraic plan for a later class-group action computation.
///
/// A plan factors a target reduced form class as a product of local ideal
/// classes:
///
/// [𝔞] = [𝔭₁]^{e₁} ··· [𝔭ᵣ]^{eᵣ}.
///
/// This type records the algebraic factorization needed by the geometric
/// layer. Each nonzero factor requires its own compatible crater-orientation
/// witness before it can move a graph node.
///
/// Note: This type does **not** execute isogenies and does **not** certify an
/// arithmetic crater orientation.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ClassGroupActionPlan {
    discriminant: QuadraticDiscriminant,
    target_form: BinaryQuadraticForm,
    factors: Vec<ClassGroupActionPlanFactor>,
    generated_subgroup_order: usize,
    ambient_class_number: usize,
}

impl ClassGroupActionPlan {
    /// Builds an algebraic action plan from local prime-norm ideals.
    ///
    /// The target must be a reduced representative in `class_group`. Each local
    /// ideal is translated to a reduced form class, the generated subgroup is
    /// closed algebraically, and a BFS records one exponent vector whose product
    /// is the target class.
    ///
    /// Complexity: `O(m · h(D) · C)`, where `m` is the number of local ideals,
    /// `h(D)` is the class number, and `C` is the cost of one form composition.
    pub fn from_local_ideals(
        class_group: &QuadraticClassGroup,
        local_ideals: &[PrimeNormIdeal],
        target_form: &BinaryQuadraticForm,
    ) -> Result<Self, ClassGroupActionPlanError> {
        if local_ideals.is_empty() {
            return Err(ClassGroupActionPlanError::EmptyLocalGeneratorSet);
        }
        class_group
            .validate_reduced_member(target_form)
            .map_err(|source| match source {
                BinaryQuadraticFormError::ClassGroupDiscriminantMismatch => {
                    ClassGroupActionPlanError::TargetDiscriminantMismatch
                }
                _ => ClassGroupActionPlanError::TargetNotReducedMember,
            })?;

        let generator_forms = Self::local_generator_forms(class_group, local_ideals)?;
        let generated_subgroup = class_group
            .generated_subgroup_by_set(&generator_forms)
            .map_err(ClassGroupActionPlanError::GeneratedSubgroup)?;
        if !generated_subgroup.contains_reduced_form(target_form) {
            return Err(ClassGroupActionPlanError::TargetOutsideGeneratedSubgroup);
        }

        let exponents = Self::find_target_exponents(
            class_group,
            &generator_forms,
            generated_subgroup.elements()[0].clone(),
            target_form,
            generated_subgroup.order(),
        )?;
        let factors = ClassGroupActionPlanFactor::from_nonzero_exponents(
            local_ideals,
            generator_forms,
            exponents,
        );

        Ok(Self {
            discriminant: class_group.discriminant().clone(),
            target_form: target_form.clone(),
            factors,
            generated_subgroup_order: generated_subgroup.order(),
            ambient_class_number: generated_subgroup.class_number(),
        })
    }

    /// Returns the discriminant of the ambient class group.
    pub fn discriminant(&self) -> &QuadraticDiscriminant {
        &self.discriminant
    }

    /// Returns the target reduced form class represented by this plan.
    pub fn target_form(&self) -> &BinaryQuadraticForm {
        &self.target_form
    }

    /// Returns the number of nonzero local factors in the plan.
    pub fn nonzero_factor_count(&self) -> usize {
        self.factors.len()
    }

    /// Returns the nonzero local factors in the plan.
    pub(crate) fn factors(&self) -> &[ClassGroupActionPlanFactor] {
        &self.factors
    }

    /// Returns whether the plan represents the principal class.
    pub fn is_trivial(&self) -> bool {
        self.factors.is_empty()
    }

    /// Returns the order of the generated subgroup used by the search.
    pub fn generated_subgroup_order(&self) -> usize {
        self.generated_subgroup_order
    }

    /// Returns the class number of the ambient class group.
    pub fn ambient_class_number(&self) -> usize {
        self.ambient_class_number
    }

    fn local_generator_forms(
        class_group: &QuadraticClassGroup,
        local_ideals: &[PrimeNormIdeal],
    ) -> Result<Vec<BinaryQuadraticForm>, ClassGroupActionPlanError> {
        let mut generator_forms = Vec::with_capacity(local_ideals.len());
        for (index, ideal) in local_ideals.iter().enumerate() {
            if ideal.order().discriminant() != class_group.discriminant() {
                return Err(
                    ClassGroupActionPlanError::LocalGeneratorDiscriminantMismatch { index },
                );
            }
            let correspondence =
                IdealFormCorrespondence::from_prime_norm_ideal(ideal).map_err(|source| {
                    ClassGroupActionPlanError::LocalGeneratorForm { index, source }
                })?;
            generator_forms.push(correspondence.reduced_form().clone());
        }

        Ok(generator_forms)
    }

    fn find_target_exponents(
        class_group: &QuadraticClassGroup,
        generator_forms: &[BinaryQuadraticForm],
        identity: BinaryQuadraticForm,
        target_form: &BinaryQuadraticForm,
        subgroup_order: usize,
    ) -> Result<Vec<usize>, ClassGroupActionPlanError> {
        if &identity == target_form {
            return Ok(vec![0; generator_forms.len()]);
        }

        let mut visited = vec![identity.clone()];
        let mut frontier = VecDeque::from([(identity, vec![0usize; generator_forms.len()])]);

        while let Some((form, exponents)) = frontier.pop_front() {
            for (index, generator) in generator_forms.iter().enumerate() {
                let product = class_group
                    .compose(&form, generator)
                    .map_err(ClassGroupActionPlanError::GeneratedSubgroup)?;
                if visited.contains(&product) {
                    continue;
                }

                let mut product_exponents = exponents.clone();
                product_exponents[index] += 1;
                if &product == target_form {
                    return Ok(product_exponents);
                }

                if visited.len() >= subgroup_order {
                    return Err(ClassGroupActionPlanError::InternalPlanningInvariantViolation);
                }
                visited.push(product.clone());
                frontier.push_back((product, product_exponents));
            }
        }

        Err(ClassGroupActionPlanError::InternalPlanningInvariantViolation)
    }
}
