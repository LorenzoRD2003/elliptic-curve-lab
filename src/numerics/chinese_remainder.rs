use core::fmt;

use crate::numerics::{gcd_biguint, inverse_mod_biguint};
use num_bigint::{BigInt, BigUint};
use num_traits::{One, Zero};

/// One congruence `x ≡ a (mod m)` with positive modulus `m ≥ 2`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Congruence {
    residue: BigUint,
    modulus: BigUint,
}

impl Congruence {
    /// Builds one normalized congruence `x ≡ residue (mod modulus)`.
    ///
    /// The stored residue is reduced modulo `modulus`, so
    /// `Congruence::new(17, 5)` stores the canonical representative `2 mod 5`.
    pub fn new(residue: BigUint, modulus: BigUint) -> Result<Self, ChineseRemainderError> {
        if modulus.is_zero() {
            return Err(ChineseRemainderError::ZeroModulus);
        } else if modulus.is_one() {
            return Err(ChineseRemainderError::ModulusOne);
        }

        Ok(Self {
            residue: &residue % &modulus,
            modulus,
        })
    }

    /// Returns the stored normalized residue `a`.
    pub fn residue(&self) -> &BigUint {
        &self.residue
    }

    /// Returns the positive modulus `m`.
    pub fn modulus(&self) -> &BigUint {
        &self.modulus
    }
}

/// One Chinese-remainder class `x ≡ a (mod M)`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ChineseRemainderSolution {
    residue: BigUint,
    modulus: BigUint,
}

impl ChineseRemainderSolution {
    pub(crate) fn new(residue: BigUint, modulus: BigUint) -> Self {
        Self { residue, modulus }
    }

    /// Returns the canonical residue `a`.
    pub fn residue(&self) -> &BigUint {
        &self.residue
    }

    /// Returns the combined modulus `M`.
    pub fn modulus(&self) -> &BigUint {
        &self.modulus
    }

    /// Returns whether `value` lies in this congruence class.
    pub fn contains(&self, value: &BigUint) -> bool {
        value % &self.modulus == self.residue
    }
}

/// Failure modes for the current coprime Chinese-remainder solver.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ChineseRemainderError {
    EmptySystem,
    ZeroModulus,
    ModulusOne,
    NonCoprimeModuli {
        left: BigUint,
        right: BigUint,
        gcd: BigUint,
    },
}

impl fmt::Display for ChineseRemainderError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptySystem => write!(
                formatter,
                "Chinese remainder reconstruction requires at least one congruence"
            ),
            Self::ZeroModulus => write!(
                formatter,
                "Chinese remainder reconstruction requires a positive modulus, not 0"
            ),
            Self::ModulusOne => write!(
                formatter,
                "Chinese remainder reconstruction requires a modulus at least 2"
            ),
            Self::NonCoprimeModuli { left, right, gcd } => write!(
                formatter,
                "Chinese remainder reconstruction currently requires pairwise coprime moduli, but gcd({left}, {right}) = {gcd}"
            ),
        }
    }
}

impl std::error::Error for ChineseRemainderError {}

/// Combines one existing congruence class with one more coprime congruence.
///
/// - `left` represents `x ≡ a (mod M)`
/// - `right` represents `x ≡ b (mod N)`
///
/// with `gcd(M, m) = 1`. This returns the unique class modulo `MN`
/// satisfying both constraints.
///
/// Complexity: `Θ(1)` Chinese-remainder combinations plus one exact extended-gcd
/// computation on the moduli.
pub fn combine_coprime_congruences(
    left: &ChineseRemainderSolution,
    right: &Congruence,
) -> Result<ChineseRemainderSolution, ChineseRemainderError> {
    let gcd = gcd_biguint(left.modulus(), right.modulus());
    if gcd != BigUint::one() {
        return Err(ChineseRemainderError::NonCoprimeModuli {
            left: left.modulus().clone(),
            right: right.modulus().clone(),
            gcd,
        });
    }

    let modulus_inverse = inverse_mod_biguint(left.modulus(), right.modulus())
        .expect("coprime moduli should have a modular inverse");
    let left_residue = BigInt::from(left.residue().clone());
    let right_residue = BigInt::from(right.residue().clone());
    let right_modulus = right.modulus().clone();
    let delta = right_residue - left_residue;
    let delta_mod_right = positive_bigint_modulus(delta, &right_modulus);
    let lift = (&delta_mod_right * &modulus_inverse) % &right_modulus;
    let combined_modulus = left.modulus() * right.modulus();
    let combined_residue = (left.residue() + (left.modulus() * &lift)) % &combined_modulus;

    Ok(ChineseRemainderSolution::new(
        combined_residue,
        combined_modulus,
    ))
}

/// Solves a finite system of pairwise coprime congruences.
///
/// The input slice represents
///
/// `x ≡ a_i (mod m_i)` for `i = 1, ..., n`
///
/// and the current implementation requires the moduli `m_i` to be pairwise
/// coprime. The returned value is the unique solution class modulo `M = Π_i m_i`.
///
/// Complexity: `Θ(n^2)` gcd checks across the input moduli, followed by
/// `Θ(n)` incremental Chinese-remainder combinations.
pub fn solve_coprime_congruences(
    congruences: &[Congruence],
) -> Result<ChineseRemainderSolution, ChineseRemainderError> {
    let (first, rest) = congruences
        .split_first()
        .ok_or(ChineseRemainderError::EmptySystem)?;

    ensure_pairwise_coprime(congruences)?;

    let mut solution =
        ChineseRemainderSolution::new(first.residue().clone(), first.modulus().clone());

    for congruence in rest {
        solution = combine_coprime_congruences(&solution, congruence)?;
    }
    Ok(solution)
}

fn ensure_pairwise_coprime(congruences: &[Congruence]) -> Result<(), ChineseRemainderError> {
    for (index, left) in congruences.iter().enumerate() {
        for right in &congruences[index + 1..] {
            let gcd = gcd_biguint(left.modulus(), right.modulus());
            if gcd != BigUint::one() {
                return Err(ChineseRemainderError::NonCoprimeModuli {
                    left: left.modulus().clone(),
                    right: right.modulus().clone(),
                    gcd,
                });
            }
        }
    }
    Ok(())
}

fn positive_bigint_modulus(value: BigInt, modulus: &BigUint) -> BigUint {
    let modulus_bigint = BigInt::from(modulus.clone());
    (((value % &modulus_bigint) + &modulus_bigint) % &modulus_bigint)
        .to_biguint()
        .expect("the normalized residue should convert back to BigUint")
}

#[cfg(test)]
mod tests {

    use super::{
        ChineseRemainderError, ChineseRemainderSolution, Congruence, combine_coprime_congruences,
        solve_coprime_congruences,
    };
    use num_bigint::BigUint;

    fn bu(value: u64) -> BigUint {
        BigUint::from(value)
    }

    #[test]
    fn congruence_constructor_normalizes_the_residue() {
        let congruence = Congruence::new(bu(17), bu(5)).expect("5 should be a valid modulus");

        assert_eq!(congruence.residue(), &bu(2));
        assert_eq!(congruence.modulus(), &bu(5));
    }

    #[test]
    fn congruence_constructor_rejects_zero_and_one_moduli() {
        assert_eq!(
            Congruence::new(bu(0), bu(0)),
            Err(ChineseRemainderError::ZeroModulus)
        );
        assert_eq!(
            Congruence::new(bu(0), bu(1)),
            Err(ChineseRemainderError::ModulusOne)
        );
    }

    #[test]
    fn combine_coprime_congruences_returns_the_expected_class() {
        let left = ChineseRemainderSolution::new(bu(2), bu(3));
        let right = Congruence::new(bu(3), bu(5)).expect("5 should be a valid modulus");

        let combined = combine_coprime_congruences(&left, &right).expect("3 and 5 are coprime");

        assert_eq!(combined.residue(), &bu(8));
        assert_eq!(combined.modulus(), &bu(15));
        assert!(combined.contains(&bu(8)));
        assert!(combined.contains(&bu(23)));
        assert!(!combined.contains(&bu(9)));
    }

    #[test]
    fn combine_coprime_congruences_rejects_non_coprime_moduli() {
        let left = ChineseRemainderSolution::new(bu(1), bu(6));
        let right = Congruence::new(bu(3), bu(9)).expect("9 should be a valid modulus");

        assert_eq!(
            combine_coprime_congruences(&left, &right),
            Err(ChineseRemainderError::NonCoprimeModuli {
                left: bu(6),
                right: bu(9),
                gcd: bu(3),
            })
        );
    }

    #[test]
    fn solve_coprime_congruences_combines_an_entire_system() {
        let congruences = vec![
            Congruence::new(bu(2), bu(3)).expect("3 should be a valid modulus"),
            Congruence::new(bu(3), bu(5)).expect("5 should be a valid modulus"),
            Congruence::new(bu(2), bu(7)).expect("7 should be a valid modulus"),
        ];

        let solution =
            solve_coprime_congruences(&congruences).expect("3, 5, 7 are pairwise coprime");

        assert_eq!(solution.residue(), &bu(23));
        assert_eq!(solution.modulus(), &bu(105));
        assert!(solution.contains(&bu(23)));
        assert!(solution.contains(&bu(128)));
    }

    #[test]
    fn solve_coprime_congruences_rejects_the_empty_system() {
        assert_eq!(
            solve_coprime_congruences(&[]),
            Err(ChineseRemainderError::EmptySystem)
        );
    }

    #[test]
    fn solve_coprime_congruences_rejects_non_coprime_inputs() {
        let congruences = vec![
            Congruence::new(bu(1), bu(6)).expect("6 should be a valid modulus"),
            Congruence::new(bu(2), bu(9)).expect("9 should be a valid modulus"),
        ];

        assert_eq!(
            solve_coprime_congruences(&congruences),
            Err(ChineseRemainderError::NonCoprimeModuli {
                left: bu(6),
                right: bu(9),
                gcd: bu(3),
            })
        );
    }
}
