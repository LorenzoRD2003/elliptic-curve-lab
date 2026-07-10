use num_traits::Zero;

use crate::isogenies::class_group_action::{ClassGroupActionPlan, ClassGroupActionPlanFactor};
use crate::visualization::Visualizable;

impl Visualizable for ClassGroupActionPlan {
    fn format_compact(&self) -> String {
        format!(
            "algebraic class-group action plan for {} using {} local factor(s)",
            self.target_form().format_compact(),
            self.nonzero_factor_count()
        )
    }

    fn describe(&self) -> String {
        describe_class_group_action_plan(self)
    }
}

fn describe_class_group_action_plan(plan: &ClassGroupActionPlan) -> String {
    let mut lines = vec![
        "Algebraic class-group action plan".to_string(),
        "---------------------------------".to_string(),
        format!("discriminant D: {}", plan.discriminant().value()),
        format!("target form class: {}", plan.target_form().format_compact()),
        format!("nonzero local factors: {}", plan.nonzero_factor_count()),
        format!(
            "generated subgroup order: {}",
            plan.generated_subgroup_order()
        ),
        format!("ambient class number h(D): {}", plan.ambient_class_number()),
        "factors:".to_string(),
    ];

    if plan.factors().is_empty() {
        lines.push("  none (principal class)".to_string());
    } else {
        for (index, factor) in plan.factors().iter().enumerate() {
            lines.push(format!("  {}. {}", index + 1, format_factor(factor)));
        }
    }

    lines.push(
        "interpretation: this is an algebraic factorization plan; executing it geometrically still requires one certified crater orientation per local ideal."
            .to_string(),
    );
    lines.join("\n")
}

fn format_factor(factor: &ClassGroupActionPlanFactor) -> String {
    format!(
        "{} with form {} and exponent {}",
        format_prime_ideal(factor),
        factor.generator_form().format_compact(),
        factor.exponent()
    )
}

fn format_prime_ideal(factor: &ClassGroupActionPlanFactor) -> String {
    if factor.ideal().root_mod_ell().is_zero() {
        format!("𝔭 = ({}, ω)", factor.ideal().norm())
    } else {
        format!(
            "𝔭 = ({}, ω - {})",
            factor.ideal().norm(),
            factor.ideal().root_mod_ell()
        )
    }
}

#[cfg(test)]
mod tests {
    use num_bigint::{BigInt, BigUint};

    use super::*;
    use crate::elliptic_curves::endomorphisms::{
        binary_quadratic_forms::{BinaryQuadraticForm, QuadraticClassGroup},
        quadratic_ideals::PrimeNormIdeal,
        quadratic_orders::{ImaginaryQuadraticOrder, QuadraticDiscriminant},
    };

    #[test]
    fn class_group_action_plan_description_stays_algebraic() {
        let class_group = QuadraticClassGroup::new(QuadraticDiscriminant::new(-84))
            .expect("D = -84 should define a class group");
        let order =
            ImaginaryQuadraticOrder::new(QuadraticDiscriminant::new(-84), BigUint::from(1u8))
                .expect("D = -84 should define an imaginary quadratic order");
        let first = PrimeNormIdeal::split(order.clone(), BigUint::from(11u8), BigUint::from(2u8))
            .expect("11 should split in D = -84");
        let second = PrimeNormIdeal::ramified(order, BigUint::from(3u8))
            .expect("3 should ramify in D = -84");
        let target = BinaryQuadraticForm::new(BigInt::from(5), BigInt::from(4), BigInt::from(5));

        let plan = ClassGroupActionPlan::from_local_ideals(&class_group, &[first, second], &target)
            .expect("Klein product should have an algebraic plan");

        let text = plan.describe();

        assert!(text.contains("Algebraic class-group action plan"));
        assert!(text.contains("target form class: (5,4,5)"));
        assert!(text.contains("𝔭 = (11, ω - 2) with form (2,2,11) and exponent 1"));
        assert!(text.contains("𝔭 = (3, ω) with form (3,0,7) and exponent 1"));
        assert!(text.contains("executing it geometrically still requires"));
        assert!(!text.contains("computed isogeny action"));
    }
}
