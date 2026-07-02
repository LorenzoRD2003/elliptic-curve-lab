# Plan Temporal: Campos Primos Grandes Con `crypto-bigint`

Este documento es un plan de migración temporal. La idea es borrarlo cuando el
milestone quede convertido en tareas chicas o cuando decidamos otra dirección.

## Objetivo

Reemplazar gradualmente el límite actual de `Fp<const P: u64>` por un backend de
campo primo de módulo constante y tamaño arbitrario basado en
[`crypto-bigint`](https://docs.rs/crypto-bigint/latest/crypto_bigint/).

El objetivo no es reintroducir `BigPrimeField` ni una familia runtime/ambient.
La dirección correcta es seguir usando tipos de campo estáticos, pero permitir
módulos grandes representados por `Uint<LIMBS>` y aritmética Montgomery
constante.

## API De `crypto-bigint` Que Importa

La superficie relevante del crate es:

- `crypto_bigint::modular::ConstMontyForm<M, LIMBS>`
- `crypto_bigint::modular::ConstMontyParams<LIMBS>`
- `crypto_bigint::modular::ConstPrimeMontyParams<LIMBS>`
- macros `const_monty_params!` y `const_prime_monty_params!`

`ConstMontyForm<M, LIMBS>` representa un residuo en forma Montgomery con módulo
constante. Provee operaciones como suma, resta, multiplicación, negación,
inversión, potencia y recuperación del representante canónico. La raíz cuadrada
requiere el bound más fuerte `M: ConstPrimeMontyParams<LIMBS>`.

## Diseño Recomendado

Mantener separadas la familia de campo y el valor:

```rust
pub struct Fp<M, const LIMBS: usize>(core::marker::PhantomData<M>);

pub struct FpElem<M, const LIMBS: usize>
where
    M: ConstMontyParams<LIMBS>,
{
    value: ConstMontyForm<M, LIMBS>,
}
```

Luego:

```rust
impl<M, const LIMBS: usize> Field for Fp<M, LIMBS>
where
    M: ConstMontyParams<LIMBS>,
{
    type Elem = FpElem<M, LIMBS>;
}
```

Esto conserva la semántica actual del repo:

- `F` es la familia de campo.
- `F::Elem` es el valor almacenado en polinomios, puntos y coeficientes.
- Las curvas siguen siendo `ShortWeierstrassCurve<F>`.

La alternativa:

```rust
pub struct Fp<M, const LIMBS: usize> {
    value: ConstMontyForm<M, LIMBS>,
}
```

también puede compilar si `Field for Fp<M, LIMBS>` define `type Elem = Self`,
pero mezcla familia y elemento en un solo tipo. Esa mezcla va contra el estilo
del repo y probablemente complique extensión, tests, `CurveModel::BaseField`,
descriptores y documentación. Preferiría evitarla salvo que encontremos una
ventaja muy concreta.

## Cambios Necesarios, Uno Por Uno

### 1. `Field::characteristic() -> u64` Ya No Sirve

Archivo principal:

- `src/fields/traits/field.rs`

Hoy:

```rust
fn characteristic() -> u64;
```

Problema:

- Un módulo de 256 bits no entra en `u64`.
- Varias capas comparan `F::characteristic()` con `2`, `3`, o lo convierten a
  `u128`.

Cambio propuesto:

```rust
pub enum FieldCharacteristic {
    Zero,
    Positive(num_bigint::BigUint),
}
```

o una newtype equivalente con helpers:

- `is_zero()`
- `is_small(value: u64)`
- `as_u64() -> Option<u64>`
- `to_biguint() -> Option<BigUint>`

Nueva API posible:

```rust
fn characteristic() -> FieldCharacteristic;

fn characteristic_u64() -> Option<u64> {
    Self::characteristic().as_u64()
}

fn has_characteristic(value: u64) -> bool {
    Self::characteristic().is_small(value)
}
```

Impacto directo:

- `Q` y `ComplexApprox` devuelven `Zero`.
- `Fp<17>` o el backend chico devuelven `Positive(17)`.
- `Fp<P256, 4>` devuelve `Positive(p)` con `p` grande.

### 2. `FiniteFieldDescriptor` Debe Dejar De Guardar `u64`

Archivo:

- `src/fields/finite_field_descriptor.rs`

Hoy:

```rust
pub characteristic: u64
pub fn cardinality(&self) -> Result<u128, FieldError>
```

Problema:

- La característica y el orden de campo pueden ser mayores que `u128`.
- Frobenius, Schoof y Hasse usan `FiniteFieldDescriptor` como fuente de verdad.

Cambio propuesto:

```rust
pub struct FiniteFieldDescriptor {
    pub characteristic: BigUint,
    pub extension_degree: NonZeroU32,
}
```

Agregar helpers:

- `characteristic_u64() -> Option<u64>`
- `cardinality_biguint() -> BigUint`
- `cardinality_u128() -> Result<u128, FieldError>`
- `is_prime_field()`

Mantener el formato educativo `F_17` / `F_(43^2)` y permitir `F_p` grande.

### 3. `FiniteField::cardinality()` Y `FiniteField::order()` Deben Separar Grande Y Chico

Archivo:

- `src/fields/traits/finite_field.rs`

Hoy:

```rust
fn cardinality() -> Option<u128>
fn order() -> Result<u128, FieldError>
fn descriptor() -> Result<FiniteFieldDescriptor, FieldError>
```

Problema:

- `order()` se usa en Hasse, Schoof, Mestre y ejemplos.
- Para primos grandes puede no entrar en `u128`.

Cambio propuesto:

```rust
fn cardinality_biguint() -> BigUint;

fn cardinality_u128() -> Result<u128, FieldError> {
    Self::cardinality_biguint()
        .try_into()
        .map_err(|_| FieldError::CardinalityOverflow)
}
```

Decisión de compatibilidad:

- Mantener `order()` como alias de `cardinality_u128()` sólo durante migración.
- Migrar call sites explícitamente a `cardinality_biguint()` o
  `cardinality_u128()` según lo que matemáticamente necesiten.

### 4. `FieldError` Necesita Errores De Conversión Grandes

Archivo:

- `src/fields/error.rs`

Hoy ya existe:

- `InvalidModulus { modulus: u64 }`
- `CardinalityOverflow`

Problema:

- Un módulo grande inválido no cabe en `u64`.
- Varios sitios van a necesitar decir “esta ruta sólo acepta característica
  chica”.

Cambios propuestos:

- `InvalidModulus { modulus: String }` o `InvalidModulusBig { modulus: String }`
- `CharacteristicOverflow`
- posiblemente `UnsupportedLargeCharacteristic(&'static str)`

Evitar volver a un `InvalidBigModulus` atado a un backend específico. El error
debe hablar de la restricción matemática/API, no de `BigPrimeField`.

### 5. `Fp<const P: u64>` No Debe Desaparecer De Golpe

Archivo:

- `src/fields/prime_field.rs`

Hoy `Fp<P>` sirve muy bien para:

- ejemplos chicos,
- tests exhaustivos,
- `EnumerableFiniteField`,
- `proptest`,
- documentación pedagógica.

Plan recomendado:

1. Mantener `Fp<const P: u64>` como backend chico.
2. Agregar un nuevo backend, por ejemplo:
   - `MontgomeryFp<M, LIMBS>`
   - o reutilizar el nombre `Fp<M, LIMBS>` sólo después de una migración más
     grande.
3. Cuando el backend grande esté maduro, decidir si:
   - `FpSmall<P>` / `Fp<M,LIMBS>` conviven,
   - o `Fp<P>` se vuelve un alias/macro de conveniencia para primos chicos.

No conviene hacer un rename masivo de entrada.

### 6. `FpElem::value() -> u64` Tiene Que Volverse Backend-Specific

Archivo:

- `src/fields/prime_field.rs`

Problema:

- Muchos tests usan `.value()` esperando `u64`.
- `ConstMontyForm::retrieve()` devuelve un `Uint<LIMBS>`, no un `u64`.

Plan:

- Para `FpElem<const P: u64>`, mantener `value() -> u64`.
- Para `MontgomeryFpElem<M,LIMBS>`, exponer:
  - `to_uint() -> Uint<LIMBS>`
  - `to_biguint() -> BigUint` si necesitamos interoperabilidad con el resto del
    repo
  - `to_u64() -> Option<u64>` para tests pequeños o compatibilidad.

No imponer un único `value()` genérico en el trait `Field`.

### 7. `from_i64` Y `elem_from_u64` Siguen Sirviendo, Pero No Alcanzan

Archivo:

- `src/fields/traits/field.rs`

Hoy:

```rust
fn from_i64(n: i64) -> Self::Elem;
fn elem_from_u64(value: u64) -> Self::Elem;
```

Estos métodos siguen siendo útiles para coeficientes chicos (`0`, `1`, `2`,
`3`, `1728`, etc.). No hace falta borrarlos.

Pero para construir elementos grandes falta algo como:

```rust
fn elem_from_biguint(value: &BigUint) -> Result<Self::Elem, FieldError>;
```

o, si queremos evitar `BigUint` en el trait base:

```rust
trait BigIntegerField: Field {
    type BigRepr;
    fn elem_from_repr(value: &Self::BigRepr) -> Result<Self::Elem, FieldError>;
}
```

Recomendación inicial:

- No tocar el trait base todavía para esto.
- Agregar constructores inherentes en el backend Montgomery.
- Generalizar sólo si aparece un caller genérico real.

### 8. `EnumerableFiniteField` Debe Quedar Para Campos Chicos

Archivo:

- `src/fields/traits/enumerative_finite_field.rs`

Esta separación ya está bien. No hay que implementar
`EnumerableFiniteField` para un primo de 256 bits.

Call sites que dependen de enumeración deben seguir siendo explícitamente
“small field”:

- `src/elliptic_curves/models/traits/finite.rs`
- `src/elliptic_curves/models/short_weierstrass/enumerable.rs`
- `src/elliptic_curves/models/short_weierstrass/group_order/quadratic_character.rs`
- `src/elliptic_curves/models/short_weierstrass/division_polynomials/torsion_search/*`
- `src/isogenies/*` para verificación exhaustiva, grafos y kernels chicos
- `src/proptest_support/*`

No migrar estas rutas a grandes primos. En su lugar, sus nombres/docs deberían
decir “small/enumerable field”.

### 9. Polinomios En Característica Positiva Usan `u64`

Archivos:

- `src/polynomials/dense/trait_impls.rs`
- `src/polynomials/sparse/trait_impls.rs`

Problema:

- El cálculo de p-th roots de polinomios usa `F::characteristic()` como `u64`.
- Divide grados por la característica y convierte grados a `u64`.

Plan:

- Si la característica no entra en `usize`/`u64`, esa ruta debe devolver
  `None` o un error explícito.
- Agregar helpers tipo:
  - `F::characteristic_usize() -> Option<usize>`
  - `F::characteristic_u64() -> Option<u64>`
- Mantener estos algoritmos como rutas de característica chica hasta que haya
  una razón para soportar grados enormes.

### 10. `PthRootExtraction` Usa Característica Como Exponente

Archivo:

- `src/fields/traits/pth_root_extraction.rs`

Problemas:

- `BigUint::from(F::characteristic())`
- `F::pow(x, F::characteristic())`
- `impl<const P: u64> PthRootExtraction for FpElem<P>`

Plan:

- Separar `pow_u64` de `pow_biguint`.
- Agregar una ruta genérica `finite_field_pow_biguint`.
- Para Montgomery, usar `ConstMontyForm::pow` si acepta el formato de exponente
  necesario, o convertir el exponente a una representación que `crypto-bigint`
  acepte.
- Mantener `PthRootExtraction` inicialmente sólo para backends chicos y
  extensiones enumerables si migrarlo amenaza con expandir demasiado el alcance.

### 11. `SqrtField` Y `QuadraticCharacterFiniteField`

Archivos:

- `src/fields/traits/sqrt_field.rs`
- `src/fields/traits/quadratic_character.rs`
- `src/fields/prime_field.rs`

`crypto-bigint` ayuda mucho acá:

- `ConstMontyForm` tiene inversión y potencia.
- `sqrt` requiere `ConstPrimeMontyParams<LIMBS>`.

Plan:

- Implementar `Field` primero con `M: ConstMontyParams<LIMBS>`.
- Implementar `SqrtField` después con `M: ConstPrimeMontyParams<LIMBS>`.
- Implementar `QuadraticCharacterFiniteField` usando el default si ya depende
  de `FiniteField::cardinality_biguint`/pow grande, o con una implementación
  específica del backend Montgomery.

El caso `char = 2` debe seguir siendo explícito.

### 12. `CurveError` Y Errores De Característica También Usan `u64`

Archivo:

- `src/elliptic_curves/error.rs`

Variantes actuales guardan `characteristic: u64` y `field_order: u128`.

Ejemplos:

- `UnsupportedCharacteristic { characteristic: u64 }`
- `InvalidFrobeniusBaseField { characteristic: u64, ... }`
- `UnsupportedFrobeniusFieldOrder { field_order: u128 }`
- `InvalidHasseIntervalFieldOrder { field_order: u128 }`
- `MestrePrimeTooSmall { characteristic: u64 }`

Plan:

- Crear tipos compartidos para metadata:
  - `FieldCharacteristic`
  - `FieldOrder`
- O guardar strings/BigUint sólo en los errores que pueden cruzar el límite
  chico/grande.
- Mantener errores `u64` en rutas explícitamente small-field si eso simplifica,
  pero nombrarlas como tal.

### 13. Modelos De Curvas Comparan Característica Con `2` Y `3`

Archivos representativos:

- `src/elliptic_curves/models/short_weierstrass/curve/type_definition.rs`
- `src/elliptic_curves/models/general_weierstrass/reduction.rs`
- `src/elliptic_curves/models/general_weierstrass/y_fiber/*`
- `src/elliptic_curves/models/montgomery/type_definition.rs`
- `src/elliptic_curves/models/montgomery/reduction.rs`
- `src/elliptic_curves/models/twisted_edwards/type_definition.rs`

Plan:

- Reemplazar `F::characteristic() == 2` por `F::has_characteristic(2)`.
- Reemplazar `matches!(F::characteristic(), 2 | 3)` por helpers semánticos.
- Los errores deben reportar una característica posiblemente grande.

Este cambio es mecánico y debería hacerse temprano.

### 14. Frobenius, Hasse, Schoof Y Mestre Usan `u128`

Archivos representativos:

- `src/elliptic_curves/frobenius/hasse/interval.rs`
- `src/elliptic_curves/frobenius/group_order.rs`
- `src/elliptic_curves/frobenius/character_sum.rs`
- `src/elliptic_curves/frobenius/metadata.rs`
- `src/elliptic_curves/frobenius/extension_counts.rs`
- `src/elliptic_curves/frobenius/schoof/report/*`
- `src/elliptic_curves/models/short_weierstrass/schoof/api.rs`
- `src/elliptic_curves/models/short_weierstrass/group_order/mestre/*`

Problemas:

- `HasseInterval::for_q` acepta `u128`.
- Muchos reports guardan `field_order: u128` y `curve_order: u128`.
- Schoof convierte entre `u128`, `i128`, `BigUint`, `BigInt`.
- Mestre y Hasse-search tienen partes small-field y partes que podrían ser
  grandes, mezcladas.

Plan:

1. Mantener `HasseInterval` chico inicialmente.
2. Crear `BigHasseInterval` o migrar internamente a `BigUint`/`BigInt`.
3. No bloquear el primer backend Montgomery en Schoof grande.
4. Primer objetivo de curva grande: construcción, grupo, lift-x, scalar
   multiplication, no conteo avanzado.
5. Migrar Schoof/Hasse después como milestone propio.

### 15. `GroupCurveModel::mul_scalar` Usa `u64`

Archivos:

- `src/elliptic_curves/models/traits/group.rs` o equivalente de traits
- implementaciones de short-Weierstrass, Montgomery, proyectivo, etc.

Problema:

- Para grandes primos, scalar multiplication con escalares grandes es
  necesario para workflows reales.

Plan:

- Mantener `mul_scalar(..., u64)` para ejemplos y compatibilidad.
- Agregar una capacidad nueva:

```rust
trait BigScalarGroupCurveModel: GroupCurveModel {
    fn mul_scalar_biguint(&self, point: &Self::Point, scalar: &BigUint)
        -> Result<Self::Point, CurveError>;
}
```

Ya hay rutas de point-order que usan `BigUint`; aprovechar esa dirección.

### 16. Isogenias Y Verificación Exhaustiva Deben Quedar En Campo Chico

Archivos representativos:

- `src/isogenies/comparison.rs`
- `src/isogenies/traits/verifiable.rs`
- `src/isogenies/graphs/*`
- `src/isogenies/scalar_multiplication/kernel_structure.rs`
- `src/isogenies/scalar_multiplication/function_fields/verschiebung.rs`

Problemas:

- Mucho código requiere `EnumerableFiniteField`.
- Varios lugares usan `F::characteristic()` como scalar `u64`.

Plan:

- No migrar estas rutas en el primer milestone.
- Renombrar/documentar como small-field donde haga falta.
- Para Verschiebung/Frobenius scalar `[p]`, permitir sólo si `p` entra en
  `u64`, o migrar esa llamada a scalar grande cuando exista.

### 17. `proptest_support` Y Generadores

Archivos:

- `src/proptest_support/fields/prime.rs`
- `src/proptest_support/fields/families.rs`
- `src/proptest_support/polynomials/*`
- `src/proptest_support/elliptic_curves/*`

Plan:

- Mantener generadores exhaustivos para `Fp<17>`, `Fp<41>`, etc.
- Agregar generadores específicos para Montgomery sólo con valores chicos
  embebidos en un módulo grande si se necesitan.
- No intentar generar elementos uniformes de campos de 256 bits al principio.

### 18. Ejemplos

Ejemplos actuales usan aliases como:

- `type F = Fp<101>`
- `type F = Fp<1_000_000_007>`

Plan:

- Mantener ejemplos chicos intactos.
- Agregar un ejemplo nuevo para un módulo grande con macro de `crypto-bigint`.
- No convertir todos los ejemplos a Montgomery.

Ejemplo objetivo:

```rust
crypto_bigint::const_prime_monty_params!(
    Secp256k1Prime,
    crypto_bigint::U256,
    "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F"
);

type F = MontyFp<Secp256k1Prime, 4>;
```

La sintaxis exacta de macro debe confirmarse en una prueba mínima al iniciar la
fase de implementación.

## Orden De Implementación Recomendado

### Fase 0: Spike Mínimo Con `crypto-bigint`

- Agregar dependencia.
- Crear un test aislado que use `ConstMontyForm` y una macro de módulo.
- Verificar suma, multiplicación, recuperación de representante, inversión y
  raíz cuadrada.
- No tocar traits del repo todavía.

### Fase 1: Metadata Grande

- Introducir `FieldCharacteristic`.
- Cambiar `Field::characteristic()`.
- Agregar helpers `has_characteristic`, `characteristic_u64`.
- Migrar comparaciones con `2` y `3`.
- Mantener tests existentes verdes.

### Fase 2: Orden De Campo Grande

- Migrar `FiniteFieldDescriptor` a `BigUint`.
- Agregar `cardinality_biguint` y `cardinality_u128`.
- Marcar rutas small-field cuando pidan `u128`.
- Mantener `EnumerableFiniteField` separado.

### Fase 3: Backend Montgomery

- Agregar `MontyFp<M, LIMBS>` y `MontyFpElem<M, LIMBS>`.
- Implementar `Field`.
- Implementar `FiniteField` sin `EnumerableFiniteField`.
- Implementar `SqrtField` con `ConstPrimeMontyParams`.
- Agregar tests contra módulos chicos y uno grande.

### Fase 4: Curvas Básicas Sobre Montgomery

- Construcción de `ShortWeierstrassCurve<F>`.
- Membership.
- `LiftXCoordinate`.
- Grupo afín/proyectivo existente.
- Scalar multiplication `u64` primero.
- Agregar scalar grande si el trait ya existe.

### Fase 5: Algoritmos No Exhaustivos

- Revisar qué rutas no necesitan enumeración:
  - Schoof
  - Hasse BSGS
  - point order from known multiple
- Migrarlas una por una a `BigUint`/`BigInt`.
- No mezclar esto con el primer backend Montgomery.

## Riesgos

- Cambiar `Field::characteristic()` toca muchas capas, pero es el cambio
  correcto. Parchearlo con `u64` opcional tarde o temprano vuelve a romper.
- `FiniteFieldDescriptor` con `BigUint` puede generar mucho ruido de clones.
  Conviene usar referencias en reports grandes.
- `crypto-bigint` tiene APIs constant-time que devuelven wrappers sutiles
  (`CtOption`, choice types, etc.). Hay que mapearlas con cuidado a nuestros
  `Option`/`Result`.
- Las rutas de enumeración no deben accidentalmente aceptar campos grandes.
- Los reports de Frobenius/Schoof tienen muchos `u128`; migrarlos merece otro
  milestone.

## Criterio De Éxito Del Primer Milestone

El primer milestone estaría bien cerrado si:

- existe un backend estático Montgomery para un primo grande,
- implementa `Field`, `FiniteField` y `SqrtField`,
- puede construir una short-Weierstrass curve existente sin duplicar lógica,
- puede validar puntos y ejecutar group law/scalar multiplication existente,
- no implementa `EnumerableFiniteField`,
- no reintroduce `BigPrimeField`, `AmbientShortWeierstrassCurve` ni feature de
  ambient fields.

## Decisión De Diseño Principal

La clave es esta:

> El backend grande tiene que entrar por la abstracción estática existente,
> aunque haya que agrandar la metadata de `Field`/`FiniteField`.

No queremos otro carril paralelo de curvas. Queremos que las curvas existentes
usen un campo más capaz.
