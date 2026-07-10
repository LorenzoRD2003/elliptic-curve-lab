use elliptic_algorithms_lab::elliptic_curves::endomorphisms::{
    BinaryQuadraticForm, QuadraticClassGroup, quadratic_orders::QuadraticDiscriminant,
};
use elliptic_algorithms_lab::visualization::Visualizable;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("Binary quadratic form class-group composition");
    println!("=============================================");
    println!();
    println!("Each case fixes a negative discriminant D, enumerates the reduced");
    println!("primitive positive-definite forms of discriminant D, composes a few");
    println!("classes by Dirichlet/Gauss composition, and prints the Cayley table.");
    println!();

    for discriminant in [-20, -23, -31, -84] {
        describe_case(discriminant)?;
    }

    Ok(())
}

fn describe_case(discriminant: i64) -> Result<(), Box<dyn std::error::Error>> {
    let class_group = QuadraticClassGroup::new(QuadraticDiscriminant::new(discriminant))?;
    let representatives = class_group.enumerate_reduced_forms();
    let table = class_group.cayley_table()?;

    println!("D = {discriminant}");
    println!("{}", "-".repeat(discriminant.to_string().len() + 4));
    println!("reduced representatives:");
    for (index, representative) in representatives.iter().enumerate() {
        println!("  f{index} = {}", representative.format_compact());
    }
    println!();

    println!("sample products:");
    println!(
        "{}",
        SampleProduct::new(&class_group, &representatives, 0, 0)?.describe()
    );
    if representatives.len() > 1 {
        println!(
            "{}",
            SampleProduct::new(&class_group, &representatives, 0, 1)?.describe()
        );
        println!(
            "{}",
            SampleProduct::new(&class_group, &representatives, 1, 1)?.describe()
        );
    }
    if representatives.len() > 2 {
        println!(
            "{}",
            SampleProduct::new(&class_group, &representatives, 1, 2)?.describe()
        );
    }
    println!();

    println!("{}", table.describe());
    println!();

    Ok(())
}

struct SampleProduct {
    left: usize,
    right: usize,
    product_index: usize,
    product: BinaryQuadraticForm,
}

impl Visualizable for SampleProduct {
    fn format_compact(&self) -> String {
        format!("f{} · f{} = f{}", self.left, self.right, self.product_index)
    }

    fn describe(&self) -> String {
        format!(
            "  {} = {}",
            self.format_compact(),
            self.product.format_compact()
        )
    }
}

impl SampleProduct {
    fn new(
        class_group: &QuadraticClassGroup,
        representatives: &[BinaryQuadraticForm],
        left: usize,
        right: usize,
    ) -> Result<Self, Box<dyn std::error::Error>> {
        let product = class_group.compose(&representatives[left], &representatives[right])?;
        let product_index = representatives
            .iter()
            .position(|representative| representative == &product)
            .expect("sample product should be one of the enumerated representatives");

        Ok(Self {
            left,
            right,
            product_index,
            product,
        })
    }
}
