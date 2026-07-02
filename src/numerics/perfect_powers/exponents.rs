use num_bigint::BigUint;

/// Lists the only exponents that need to be tested.
///
/// If `N = a^e` with composite `e = rs`, then `N = (a^r)^s`; therefore it is
/// enough to test prime exponents. The caller supplies the bound
/// `max_exponent = ⌊log₂ N⌋`, since `a ≥ 2` implies `e ≤ ⌊log₂ N⌋`.
///
/// A small Eratosthenes sieve is a better fit than independent primality tests.
///
/// Complexity: `Θ(n log log n)` time and `Θ(n)` memory.
pub(super) fn prime_exponents_through(max_exponent: u32) -> Vec<u32> {
    let Ok(limit) = usize::try_from(max_exponent) else {
        return Vec::new();
    };
    if limit < 2 {
        return Vec::new();
    }

    let mut is_prime = vec![true; limit + 1];
    is_prime[0] = false;
    is_prime[1] = false;

    let mut candidate = 2usize;
    while candidate <= limit / candidate {
        if is_prime[candidate] {
            let mut multiple = candidate * candidate;
            while multiple <= limit {
                is_prime[multiple] = false;
                multiple += candidate;
            }
        }
        candidate += 1;
    }

    is_prime
        .iter()
        .enumerate()
        .filter_map(|(candidate, is_prime)| {
            if *is_prime {
                u32::try_from(candidate).ok()
            } else {
                None
            }
        })
        .collect()
}

/// Chooses the Hensel prime used for the candidate exponent `q`.
///
/// The staged perfect-power route assumes `gcd(N, 6) = 1`. For `q = 2`, it uses
/// `p = 3`; for every odd prime exponent, it uses `p = 2`. In both cases
/// `p ∤ q`, and any genuine base `a` is a unit modulo `p`. Thus the derivative
/// of `x^q − N` at a genuine root, namely `q·a^(q−1)`, is non-zero modulo `p`.
///
/// Complexity: `Θ(1)`.
pub(super) fn hensel_prime_for_exponent(exponent: u32) -> BigUint {
    if exponent == 2 {
        BigUint::from(3u8)
    } else {
        BigUint::from(2u8)
    }
}
