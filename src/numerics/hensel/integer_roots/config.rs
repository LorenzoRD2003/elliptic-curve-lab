use num_bigint::BigUint;

const DEFAULT_MAX_SEED_SCAN: u64 = 10_000;

/// Configuration for scanning simple roots modulo one prime and recovering
/// integer roots inside a bounded interval.
///
/// The search is intentionally a finite, explicit scan over residues
/// `0, 1, …, p − 1`. This keeps the first reusable surface educational and
/// predictable: a caller chooses a prime `p`, an integer-root bound `B₀`, and a
/// maximum prime size worth scanning. The bound means that any certified root
/// `r ∈ ℤ` must satisfy `|r| ≤ B₀`.
///
/// If the scan completes, an empty result means:
///
/// `No integer root with |r| ≤ B₀ was certified from a simple root modulo p.`
///
/// This is not yet the same as a complete theorem that no bounded integer root
/// exists, because a genuine integer root can reduce to a singular modular root
/// for an unlucky prime. Later stages can try more primes or add singular
/// lifting. For the simple-root case, the conclusion is exact: every simple
/// seed modulo this prime was tested, lifted to precision `pᵉ > 2B₀`, centered
/// to the only possible bounded integer representative, and checked by exact
/// integer evaluation.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct HenselIntegerRootSearchConfig {
    prime: BigUint,
    root_bound: BigUint,
    max_seed_scan: u64,
}

impl HenselIntegerRootSearchConfig {
    /// Creates a search config with a conservative default residue-scan limit.
    ///
    /// Complexity: `Θ(1)` big-integer moves.
    pub(crate) fn new(prime: BigUint, root_bound: BigUint) -> Self {
        Self {
            prime,
            root_bound,
            max_seed_scan: DEFAULT_MAX_SEED_SCAN,
        }
    }

    /// Sets the maximum prime value for which the search will enumerate every
    /// residue modulo `p`.
    ///
    /// Complexity: `Θ(1)`.
    pub(crate) fn with_max_seed_scan(mut self, max_seed_scan: u64) -> Self {
        self.max_seed_scan = max_seed_scan;
        self
    }

    /// Returns the prime used for the modular seed scan.
    pub(crate) fn prime(&self) -> &BigUint {
        &self.prime
    }

    /// Returns the certified absolute root bound `B₀`.
    pub(crate) fn root_bound(&self) -> &BigUint {
        &self.root_bound
    }

    /// Returns the largest prime value the config will scan exhaustively.
    pub(crate) fn max_seed_scan(&self) -> u64 {
        self.max_seed_scan
    }
}
