use num_bigint::{BigInt, BigUint};

/// Inputs accepted by non-negative scalar-multiplication APIs.
///
/// This keeps educational call sites ergonomic for literals such as `2` while
/// normalizing the executable scalar to arbitrary-precision [`BigUint`].
pub trait ScalarInput {
    /// Converts the input into one non-negative scalar.
    fn into_biguint_scalar(self) -> BigUint;
}

impl ScalarInput for BigUint {
    fn into_biguint_scalar(self) -> BigUint {
        self
    }
}

impl ScalarInput for &BigUint {
    fn into_biguint_scalar(self) -> BigUint {
        self.clone()
    }
}

macro_rules! impl_unsigned_scalar_input {
    ($($type:ty),* $(,)?) => {
        $(
            impl ScalarInput for $type {
                fn into_biguint_scalar(self) -> BigUint {
                    BigUint::from(self)
                }
            }
        )*
    };
}

macro_rules! impl_signed_scalar_input {
    ($($type:ty),* $(,)?) => {
        $(
            impl ScalarInput for $type {
                fn into_biguint_scalar(self) -> BigUint {
                    assert!(
                        !self.is_negative(),
                        "non-negative scalar multiplication received a negative scalar"
                    );
                    BigInt::from(self)
                        .to_biguint()
                        .expect("non-negative signed scalar should convert to BigUint")
                }
            }
        )*
    };
}

impl_unsigned_scalar_input!(u8, u16, u32, u64, u128, usize);
impl_signed_scalar_input!(i8, i16, i32, i64, i128, isize);
