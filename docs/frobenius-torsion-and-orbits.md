# Frobenius, Torsion, and Orbits

Source files:

- [src/elliptic_curves/frobenius/torsion.rs](../src/elliptic_curves/frobenius/torsion.rs)
- [src/elliptic_curves/frobenius/orbit.rs](../src/elliptic_curves/frobenius/orbit.rs)

This note explains the mathematics behind the current Frobenius-on-torsion and
Frobenius-orbit helpers.

The intended reader is someone who is comfortable with the idea of an elliptic
curve equation such as $y^2 = x^3 + ax + b$, but is new to finite fields, torsion, and Frobenius.

## Why These Modules Exist

Over a finite field, an elliptic curve has only finitely many rational points,
so two natural questions immediately appear:

1. which points have order dividing a given integer $n$?
2. how does Frobenius move those points around?

The file `torsion.rs` studies the first question together with the pointwise
action of Frobenius on exact-$n$ torsion.

The file `orbit.rs` studies the second question at the orbit level: instead of
looking at one point at a time, it groups points into cycles under a chosen
Frobenius map.

Those two viewpoints are closely related:

- the torsion report tells us what happens to each torsion point
- the orbit report tells us how those pointwise motions fit together into
  cycles

## 1. The Basic Finite-Field Setting

Let $p$ be a prime and let $q = p^r$.

We often work with an elliptic curve $E$ over a finite field $\mathbb{F}_q$.
That means the coefficients of the defining equation lie in $\mathbb{F}_q$,
and we are interested in the rational point set $E(\mathbb{F}_q).$

This set is finite and forms an abelian group, with identity element $O$, the
point at infinity. So when we say that a point $P$ has order $n$, we mean:$
[n]P = O$, and $n$ is the smallest positive integer with that property.

## 2. What “Torsion” Means Here

For any integer $n \ge 1$, the $n$-torsion subgroup is

$$E[n] = \{P \in E(\overline{\mathbb{F}}_p) : [n]P = O\},$$

where $\overline{\mathbb{F}}_p$ is an algebraic closure.

This is the mathematically largest place where all finite-field torsion lives.
But the current codebase is intentionally more concrete and more modest:

- it works with a represented finite field such as $\mathbb{F}_p$ or
  $\mathbb{F}_{p^2}$
- it enumerates the rational points that are visible in that represented field

So the current torsion helpers do not try to compute all of
$E[n](\overline{\mathbb{F}}_p)$.
They compute the exact-$n$ torsion points that are already visible in the
current rational point set. In other words, the practical object is:

$$\{P \in E(\mathbb{F}_q) : P \text{ has exact order } n\}.$$

That distinction is important for beginners:

- $E[n]$ means torsion over the algebraic closure
- the current code usually means exact-$n$ torsion inside one represented
  finite field

## 3. Relative and Absolute Frobenius

There are two Frobenius maps that matter here.

### Relative Frobenius

If the current base field is $\mathbb{F}_q$, the relative Frobenius is

$$\pi_q : (x, y) \mapsto (x^q, y^q).$$

For a point already in $E(\mathbb{F}_q)$, we have $x^q = x$ and $y^q = y$.
So on rational points,

$$\pi_q(P) = P \qquad \forall P \in E(\mathbb{F}_q).$$

This is why the current relative-Frobenius torsion report is mathematically
tautological: every listed rational point is fixed. That is not useless. It provides:

- a clean API surface
- notation consistent with later extension work
- a direct contrast with the absolute Frobenius story

### Absolute Frobenius

If the characteristic is $p$, the absolute Frobenius is

$$\pi_p : (x, y) \mapsto (x^p, y^p).$$

When we represent points in a larger field such as $\mathbb{F}_{p^r}$, this
map need not fix every point. Instead, it can permute them.

This is the first genuinely interesting finite-field phenomenon in the current
Frobenius layer.

For example, a point may lie in $E(\mathbb{F}_{p^2})$ but not in
$E(\mathbb{F}_p)$. Then:

- it is visible in the represented field
- but $\pi_p(P)$ may be a different point
- only after applying $\pi_p$ twice do we return to $P$

So $\pi_p$ detects whether a point is already defined over the prime field or
only appears over an extension.

## 4. When Does Absolute Frobenius Stay on the Same Curve?

This subtlety matters a lot.

If a curve is defined over $\mathbb{F}_p$, then applying $\pi_p$ to its point
coordinates keeps us on the same curve.

But if the coefficients of the curve live only in $\mathbb{F}_{p^r}$, then
applying $\pi_p$ to a point may land on the Frobenius twist of the curve
rather than on the original curve itself.

So there are really two different statements:

1. Frobenius acts on coordinates.
2. Frobenius defines an endomorphism of the current curve model.

Those are not the same.

The current `orbit.rs` code is honest about this: the absolute-orbit helpers
first check that the chosen power $\pi_p^k$ preserves the current curve model.
If it does not, the helper returns an error instead of pretending the map is an
endomorphism of the current curve.

## 5. Exact-Order Torsion and Frobenius

Suppose now that we have enumerated the points of exact order $n$ in
$E(\mathbb{F}_q)$.

The torsion report records, for each such point $P$:

- the point $P$
- its Frobenius image
- whether the chosen Frobenius fixes it
- in the absolute case, the smallest positive power that fixes it

This is the meaning of the struct
`FrobeniusOnExactTorsionPoint<P>`.

### Fixed versus moved

If $\pi(P) = P$, the point is fixed by the chosen Frobenius action.

If $\pi(P) \ne P$, the point is moved.

For the relative Frobenius on $E(\mathbb{F}_q)$, every point is fixed.

For the absolute Frobenius on a curve viewed over $\mathbb{F}_{p^r}$, some
points may be fixed and others moved. This is the first signal that some
torsion points already descend to $\mathbb{F}_p$ while others only become
visible over the extension.

### Prime-field rational versus extension-only

If a base-defined curve is viewed over $\mathbb{F}_{p^r}$, then:

- points fixed by $\pi_p$ are already $\mathbb{F}_p$-rational
- points moved by $\pi_p$ are visible in the chosen extension field but do not
  yet descend to $\mathbb{F}_p$.

## 6. The Minimal Fixing Power

For the absolute Frobenius, one very useful invariant is the smallest positive
integer $d$ such that

$$
\pi_p^d(P) = P.
$$

This is stored in the report as the `minimal_absolute_frobenius_fixing_power`.
Why is this useful? Because it tells us the first extension degree over which the point becomes
visible as a rational point.

Very roughly:

- $d = 1$ means the point is already defined over $\mathbb{F}_p$
- $d = 2$ means it first stabilizes over an $\mathbb{F}_{p^2}$ viewpoint
- more generally, $d$ measures the Frobenius period of the point

Inside the represented field, this is the most concrete way to say:
“How far away from the prime field does this point really live?”

The helper `fixed_by_absolute_frobenius_power(k)` is then based on the simple
divisibility fact:

$$\pi_p^k(P) = P \quad \Longleftrightarrow \quad d \mid k,$$

where $d$ is the minimal fixing power.

## 7. Frobenius Orbits

Now we move from individual points to cycles.

Given a map $\pi$ and a starting point $P$, the Frobenius orbit of $P$ is

$$\mathcal{O}_\pi(P) = \{P, \pi(P), \pi^2(P), \pi^3(P), \dots\}.$$

Since the relevant sets are finite in the current codebase, this process
eventually repeats, so the orbit becomes a cycle. The orbit _period_ is the
smallest positive integer $m$ such that $\pi^m(P) = P$.

This is exactly what `FrobeniusOrbit<P>` stores:

- a chosen start point
- the distinct points in cyclic order
- the period, which is just the orbit length

## 8. Relative Orbits versus Absolute Orbits

### Relative orbits

For the relative Frobenius on $E(\mathbb{F}_q)$, every point is fixed, so every
orbit is a singleton:

$$\mathcal{O}_{\pi_q}(P) = \{P\}.$$

The current relative-orbit API makes explicit that the relative Frobenius
is trivial on the currently represented rational point set.

### Absolute orbits

For the absolute Frobenius on a curve over $\mathbb{F}_{p^r}$, orbits can have
size greater than $1$.

The smallest nontrivial example is usually an orbit of size $2$ over
$\mathbb{F}_{p^2}$:

$$P \mapsto \pi_p(P) \mapsto P$$

This means:

- $P$ is not fixed by $\pi_p$
- but it is fixed by $\pi_p^2$
- so it is naturally an extension-field point rather than a prime-field point

This is exactly the kind of phenomenon that `absolute_frobenius_orbit(...)`
and `absolute_frobenius_orbits_on_points(...)` are meant to expose.

## 9. Why Orbits Partition the Point Set

Whenever a map sends a finite set to itself, every point belongs to exactly one
orbit. So the full point set decomposes into disjoint Frobenius orbits.

That means orbit data is not just decorative. It gives a structural picture of
the whole rational point set:

- fixed points produce singleton orbits
- points defined only over extensions produce larger cycles

In the torsion setting, this becomes even more meaningful:

- the exact-$n$ torsion set is stable under Frobenius
- so we can partition exact-$n$ torsion into Frobenius orbits

That is why `FrobeniusOnExactTorsionReport::orbits()` is such a natural bridge
between `torsion.rs` and `orbit.rs`.

## 10. Torsion Reports and Orbit Reports Are Telling the Same Story

The pointwise report says:

- here is each exact-$n$ torsion point
- here is where Frobenius sends it

The orbit report says:

- if you keep applying Frobenius, these points cycle together

So the torsion report and orbit report are two resolutions of the same
information:

- one at the point level
- one at the cycle level

This is why the current implementation can compute torsion orbits directly from
the stored point-to-image action, without recomputing Frobenius from the curve
again.

Mathematically, that is exactly what should happen: once you know the action of
Frobenius on each point of a finite set, the orbit decomposition is determined.

## 11. The Key Beginner Intuition

A useful way to remember the whole story is:

- torsion asks “which points die under multiplication by $n$?”
- Frobenius asks “over which finite field are these points really defined?”

The first question is group-theoretic.
The second question is field-theoretic.

Putting them together gives a refined picture:

- not just “this is an exact-$n$ torsion point”
- but also “this torsion point is already visible over $\mathbb{F}_p$” or
  “this torsion point only appears over $\mathbb{F}_{p^2}$”, and so on.

## 12. How the Current API Fits the Theory

### In `torsion.rs`

- `relative_frobenius_on_exact_torsion(curve, n)`
  - computes the exact-$n$ rational torsion points in the current represented
    field
  - applies $\pi_q$
  - currently yields a tautological fixed-point report

- `absolute_frobenius_on_exact_torsion(curve, n, k)`
  - computes the exact-$n$ rational torsion points in the current represented
    field
  - applies $\pi_p^k$
  - records which points are fixed or moved
  - records the minimal positive fixing power

- `FrobeniusOnExactTorsionReport::orbits()`
  - recovers the cycle decomposition of the stored Frobenius action on the
    exact-$n$ torsion set

### In `orbit.rs`

- `relative_frobenius_orbit(curve, point)`
  - returns the singleton orbit of a rational point under $\pi_q$

- `relative_frobenius_orbits_on_points(curve)`
  - partitions the full rational point set into singleton orbits

- `absolute_frobenius_orbit(curve, point, k)`
  - returns the orbit of one point under $\pi_p^k$

- `absolute_frobenius_orbits_on_points(curve, k)`
  - partitions the full rational point set into absolute-Frobenius orbits
