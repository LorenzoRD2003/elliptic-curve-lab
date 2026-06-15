use core::ops::{Add, Mul};
use num_traits::{One, Zero};

/// Exact linear recurrence of order `2`
///
/// `s_n = a s_{n-1} + b s_{n-2}`
///
/// together with the initial values `s_0` and `s_1`.
///
/// This value object is intended for small shared exact recurrences that arise
/// in several mathematical domains. Its main use is that the same recurrence
/// can be evaluated in two complementary ways:
///
/// - one isolated term `s_n` via binary exponentiation of the companion matrix
/// - the whole prefix `s_0, s_1, ..., s_N` by one forward linear pass
///
/// Mathematically, if `M = [[a, b], [1, 0]]`, then for every `n >= 1`,
/// `[[s_n], [s_{n-1}]] = M^(n-1) [[s_1], [s_0]]`.
///
/// The `nth_term(...)` method uses this identity with repeated squaring.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct OrderTwoLinearRecurrence<T> {
    /// The initial value `s_0`.
    initial_zero: T,
    /// The initial value `s_1`.
    initial_one: T,
    /// The coefficient `a`.
    coefficient_prev1: T,
    /// The coefficient `b`.
    coefficient_prev2: T,
}

impl<T> OrderTwoLinearRecurrence<T>
where
    T: Clone + Add<Output = T> + Mul<Output = T> + Zero + One,
{
    /// Builds the recurrence
    ///
    /// `s_0 = initial_zero`,
    /// `s_1 = initial_one`,
    /// `s_n = coefficient_prev1 * s_{n-1} + coefficient_prev2 * s_{n-2}`.
    ///
    /// Complexity: `Θ(1)`.
    pub fn new(
        initial_zero: T,
        initial_one: T,
        coefficient_prev1: T,
        coefficient_prev2: T,
    ) -> Self {
        Self {
            initial_zero,
            initial_one,
            coefficient_prev1,
            coefficient_prev2,
        }
    }

    /// Returns the initial value `s_0`.
    pub fn initial_zero(&self) -> &T {
        &self.initial_zero
    }

    /// Returns the initial value `s_1`.
    pub fn initial_one(&self) -> &T {
        &self.initial_one
    }

    /// Returns the coefficient multiplying `s_{n-1}`.
    pub fn coefficient_prev1(&self) -> &T {
        &self.coefficient_prev1
    }

    /// Returns the coefficient multiplying `s_{n-2}`.
    pub fn coefficient_prev2(&self) -> &T {
        &self.coefficient_prev2
    }

    /// Returns the term `s_n`.
    ///
    /// This uses binary exponentiation of the companion matrix
    /// `[[a, b], [1, 0]]` attached to the recurrence
    /// `s_n = a s_{n-1} + b s_{n-2}`.
    ///
    /// Complexity: `Θ(log n)` matrix multiplications of fixed `2 x 2` size,
    /// under the usual unit-cost arithmetic model for `T`.
    pub fn nth_term(&self, index: u64) -> T {
        match index {
            0 => self.initial_zero.clone(),
            1 => self.initial_one.clone(),
            _ => {
                let companion = self.companion_matrix();
                let power = matrix_pow(&companion, index - 1);
                let (value, _) = power.apply_to_pair(&self.initial_one, &self.initial_zero);
                value
            }
        }
    }

    /// Returns the whole prefix `s_0, s_1, ..., s_max_index`.
    ///
    /// This is the preferred surface when a caller needs many consecutive
    /// values instead of one isolated term.
    ///
    /// Complexity: `Θ(max_index)` scalar additions and multiplications, and
    /// `Θ(max_index)`output storage.
    pub fn terms_through(&self, max_index: usize) -> Vec<T> {
        if max_index == 0 {
            return vec![self.initial_zero.clone()];
        }

        let mut terms = Vec::with_capacity(max_index + 1);
        terms.push(self.initial_zero.clone());
        terms.push(self.initial_one.clone());

        let mut previous = self.initial_zero.clone();
        let mut current = self.initial_one.clone();

        for _ in 2..=max_index {
            let next = self.next_term_from_pair(&previous, &current);
            terms.push(next.clone());
            previous = current;
            current = next;
        }

        terms
    }

    fn companion_matrix(&self) -> Matrix2<T> {
        Matrix2::new(
            self.coefficient_prev1.clone(),
            self.coefficient_prev2.clone(),
            T::one(),
            T::zero(),
        )
    }

    fn next_term_from_pair(&self, previous: &T, current: &T) -> T {
        self.coefficient_prev1.clone() * current.clone()
            + self.coefficient_prev2.clone() * previous.clone()
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct Matrix2<T> {
    a11: T,
    a12: T,
    a21: T,
    a22: T,
}

impl<T> Matrix2<T>
where
    T: Clone + Add<Output = T> + Mul<Output = T> + Zero + One,
{
    fn new(a11: T, a12: T, a21: T, a22: T) -> Self {
        Self { a11, a12, a21, a22 }
    }

    fn identity() -> Self {
        Self::new(T::one(), T::zero(), T::zero(), T::one())
    }

    fn multiply(&self, other: &Self) -> Self {
        Self::new(
            self.a11.clone() * other.a11.clone() + self.a12.clone() * other.a21.clone(),
            self.a11.clone() * other.a12.clone() + self.a12.clone() * other.a22.clone(),
            self.a21.clone() * other.a11.clone() + self.a22.clone() * other.a21.clone(),
            self.a21.clone() * other.a12.clone() + self.a22.clone() * other.a22.clone(),
        )
    }

    fn apply_to_pair(&self, top: &T, bottom: &T) -> (T, T) {
        (
            self.a11.clone() * top.clone() + self.a12.clone() * bottom.clone(),
            self.a21.clone() * top.clone() + self.a22.clone() * bottom.clone(),
        )
    }
}

fn matrix_pow<T>(matrix: &Matrix2<T>, exponent: u64) -> Matrix2<T>
where
    T: Clone + Add<Output = T> + Mul<Output = T> + Zero + One,
{
    let mut result = Matrix2::identity();
    let mut base = matrix.clone();
    let mut exponent = exponent;

    while exponent > 0 {
        if exponent % 2 == 1 {
            result = result.multiply(&base);
        }
        exponent /= 2;
        if exponent > 0 {
            base = base.multiply(&base);
        }
    }

    result
}

#[cfg(test)]
mod tests {
    use crate::numerics::OrderTwoLinearRecurrence;

    #[test]
    fn fibonacci_values_match_the_classical_sequence() {
        let recurrence = OrderTwoLinearRecurrence::new(0i128, 1i128, 1i128, 1i128);

        assert_eq!(recurrence.nth_term(0), 0);
        assert_eq!(recurrence.nth_term(1), 1);
        assert_eq!(recurrence.nth_term(2), 1);
        assert_eq!(recurrence.nth_term(10), 55);
    }

    #[test]
    fn lucas_values_match_the_classical_sequence() {
        let recurrence = OrderTwoLinearRecurrence::new(2i128, 1i128, 1i128, 1i128);

        assert_eq!(recurrence.nth_term(0), 2);
        assert_eq!(recurrence.nth_term(1), 1);
        assert_eq!(recurrence.nth_term(5), 11);
        assert_eq!(recurrence.nth_term(8), 47);
    }

    #[test]
    fn terms_through_matches_repeated_nth_queries() {
        let recurrence = OrderTwoLinearRecurrence::new(2i128, 3i128, 4i128, -5i128);
        let terms = recurrence.terms_through(12);

        assert_eq!(terms.len(), 13);
        for (index, term) in terms.iter().enumerate() {
            assert_eq!(*term, recurrence.nth_term(index as u64));
        }
    }

    #[test]
    fn frobenius_style_recurrence_is_supported() {
        let trace = 4i128;
        let field_order = 43i128;
        let recurrence = OrderTwoLinearRecurrence::new(2i128, trace, trace, -field_order);

        assert_eq!(recurrence.terms_through(4), vec![2, 4, -70, -452, 1202]);
    }
}
