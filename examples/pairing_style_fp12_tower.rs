use elliptic_algorithms_lab::fields::traits::*;
use elliptic_algorithms_lab::fields::{
    FieldError,
    extension_field::{ExtensionField, ExtensionFieldSpec},
    polynomial_field::PolynomialModulus,
};
use elliptic_algorithms_lab::visualization::{Visualizable, VisualizableField};

type Fp19 = elliptic_algorithms_lab::fields::Fp19;

struct Fp2Spec;

impl ExtensionFieldSpec for Fp2Spec {
    type Base = Fp19;

    const NAME: &'static str = "Fp2 over F19";

    fn defining_modulus() -> PolynomialModulus<Self::Base> {
        PolynomialModulus::<Fp19>::new(vec![Fp19::one(), Fp19::zero(), Fp19::one()])
            .expect("u^2 + 1 should be structurally valid")
    }

    fn check_field_conditions() -> Result<(), FieldError> {
        Self::defining_modulus().check_field_modulus_requirements()
    }
}

type Fp2 = ExtensionField<Fp2Spec>;

struct Fp6Spec;

impl ExtensionFieldSpec for Fp6Spec {
    type Base = Fp2;

    const NAME: &'static str = "Fp6 over Fp2";

    fn defining_modulus() -> PolynomialModulus<Self::Base> {
        // Pairing-style shape: Fp6 = Fp2[v] / (v^3 - xi), with xi = 1 + u.
        let xi = Fp2::element(vec![Fp19::one(), Fp19::one()]);

        PolynomialModulus::<Fp2>::new(vec![Fp2::neg(&xi), Fp2::zero(), Fp2::zero(), Fp2::one()])
            .expect("v^3 - xi should be structurally valid")
    }

    fn check_field_conditions() -> Result<(), FieldError> {
        // TODO: replace this manual acceptance with a generic irreducibility
        // backend over algebraic extension bases once the crate supports it.
        Ok(())
    }
}

type Fp6 = ExtensionField<Fp6Spec>;

struct Fp12Spec;

impl ExtensionFieldSpec for Fp12Spec {
    type Base = Fp6;

    const NAME: &'static str = "Fp12 over Fp6";

    fn defining_modulus() -> PolynomialModulus<Self::Base> {
        // Pairing-style shape: Fp12 = Fp6[w] / (w^2 - v), where v is the
        // degree-one generator of Fp6 over Fp2.
        let v = Fp6::element(vec![Fp2::zero(), Fp2::one()]);

        PolynomialModulus::<Fp6>::new(vec![Fp6::neg(&v), Fp6::zero(), Fp6::one()])
            .expect("w^2 - v should be structurally valid")
    }

    fn check_field_conditions() -> Result<(), FieldError> {
        // TODO: replace this manual acceptance with a generic irreducibility
        // backend over algebraic extension bases once the crate supports it.
        Ok(())
    }
}

type Fp12 = ExtensionField<Fp12Spec>;

fn main() -> Result<(), FieldError> {
    Fp2::check_structure()?;
    Fp6::check_structure()?;
    Fp12::check_structure()?;

    println!("Pairing-style educational tower");
    println!("===============================");
    println!();

    println!("Base field:");
    println!("  GF(19)");
    println!();

    println!("First extension:");
    println!("{}", describe_extension::<Fp2Spec>());
    println!();

    println!("Second extension:");
    println!("{}", describe_extension::<Fp6Spec>());
    println!();

    println!("Third extension:");
    println!("{}", describe_extension::<Fp12Spec>());
    println!();

    let u = Fp2::element(vec![Fp19::zero(), Fp19::one()]);
    let xi = Fp2::element(vec![Fp19::one(), Fp19::one()]);
    let v = Fp6::element(vec![Fp2::zero(), Fp2::one()]);
    let w = Fp12::element(vec![Fp6::zero(), Fp6::one()]);

    println!("Tower generators:");
    println!("  u in Fp2   = {}", u.format_compact());
    println!("  v in Fp6   = {}", v.format_compact());
    println!("  w in Fp12  = {}", w.format_compact());
    println!();

    let u_squared = Fp2::mul(&u, &u);
    let v_cubed = Fp6::mul(&Fp6::mul(&v, &v), &v);
    let w_squared = Fp12::mul(&w, &w);

    println!("Tower shorthand:");
    println!("  xi   = 1 + u");
    println!("  u^2  = -1");
    println!("  v^3  = xi");
    println!("  w^2  = v");
    println!();

    println!("Defining relations after quotient reduction:");
    println!("  xi   = {}", xi.format_compact());
    println!("  u^2  = {}", u_squared.format_compact());
    println!("  v^3  = {}", v_cubed.format_compact());
    println!("  w^2  = {}", w_squared.format_compact());
    println!();

    let element = Fp12::element(vec![
        Fp6::element(vec![Fp2::from_i64(3), Fp2::from_i64(1)]),
        Fp6::element(vec![Fp2::from_i64(2)]),
    ]);
    let conjugate_like = Fp12::element(vec![
        Fp6::element(vec![Fp2::from_i64(3), Fp2::from_i64(1)]),
        Fp6::neg(&Fp6::element(vec![Fp2::from_i64(2)])),
    ]);
    let product = Fp12::mul(&element, &conjugate_like);

    println!("Sample Fp12 element:");
    println!("{}", element.describe());
    println!();

    println!("A useful multiplication trace:");
    println!(
        "{}",
        VisualizableField::explain_mul(&element, &conjugate_like)
            .expect("Fp12 multiplication is visualizable")
    );
    println!();

    println!("Result of e * e':");
    println!("  {}", product.format_compact());

    Ok(())
}

fn describe_extension<S>() -> String
where
    S: ExtensionFieldSpec,
    S::Base: Field,
    <S::Base as Field>::Elem: VisualizableField + std::fmt::Display,
{
    [
        format!("Extension field: {}", S::NAME),
        "defining modulus:".to_string(),
        S::defining_modulus().describe(),
    ]
    .join("\n")
}
