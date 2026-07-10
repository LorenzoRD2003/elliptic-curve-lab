# Plan de implementación: `QuadraticClassGroup::compose()`

Fecha: 2026-07-10

## Resumen

Queremos implementar composición de clases de formas cuadráticas binarias
primitivas positivas definidas de discriminante `D < 0`:

```text
(a,b,c) * (a′,b′,c′) = reduce(A,B,C)
```

La API final deseada es:

```rust
impl QuadraticClassGroup {
    pub fn compose(
        &self,
        left: &BinaryQuadraticForm,
        right: &BinaryQuadraticForm,
    ) -> Result<BinaryQuadraticForm, BinaryQuadraticFormError>;
}
```

El resultado debe ser siempre el representante reducido positivo definido de la
clase compuesta.

La implementación no debe usar NUCOMP todavía. La prioridad es una composición
clásica correcta, legible y testeada sobre ejemplos educativos.

## Estado actual

Ya existe:

1. `BinaryQuadraticForm` como terna integral `(a,b,c)`.
2. `BinaryQuadraticForm::reduce_positive_definite()`.
3. `QuadraticClassGroup::enumerate_reduced_forms()`.
4. Validación interna de representantes reducidos del grupo.
5. Errores explícitos de membresía: discriminante incorrecto, no primitiva, no
   positiva definida, no reducida.
6. Un helper concordante testeado bajo `#[cfg(test)]`, para el caso
   `gcd(a,a′,(b+b′)/2) = 1`.
7. CRT compatible para módulos no coprimos, también todavía como scaffolding de
   test.

No existe todavía una transformación general y documentada que convierta dos
representantes arbitrarios en representantes equivalentes concordantes. Por eso
no conviene exponer `compose()` público hasta cerrar esa capa.

## Principios de diseño

1. **No usar una búsqueda acotada con constante mágica.**
   Una constante tipo `CONCORDANT_SEARCH_RADIUS = 16` no es una API honesta para
   composición general.

2. **Mantener `compose()` como operación del grupo.**
   La composición depende del discriminante fijo y de la pertenencia al grupo,
   por lo que debe vivir en `QuadraticClassGroup`, no como método libre de
   `BinaryQuadraticForm`.

3. **Aceptar solo representantes reducidos en la API pública.**
   Las formas no reducidas pueden aparecer internamente por equivalencia, pero
   la frontera pública debe ser austera:

   ```text
   input reducido + input reducido -> output reducido
   ```

4. **Mantener los helpers privados o `pub(crate)` según consumidor real.**
   La API pública estable inicial debe ser solo `compose()` y los errores
   necesarios para explicar fallos alcanzables.

5. **Separar composición concordante de reducción a concordancia.**
   Son responsabilidades distintas:

   - el helper concordante compone cuando la precondición ya vale;
   - otra capa transforma representantes equivalentes hasta satisfacer esa
     precondición.

## Etapa 1. Promover infraestructura exacta necesaria

Objetivo: dejar disponibles en build normal los helpers aritméticos que
realmente usa la composición.

Trabajo:

1. Promover `combine_compatible_congruences(...)` fuera de `#[cfg(test)]` cuando
   tenga un consumidor normal.
2. Mantenerlo `pub(crate)` en `numerics::chinese_remainder`.
3. Confirmar que usa `positive_mod_biguint(...)` y no duplica normalización.
4. Agregar tests de CRT compatible si no quedaron ya suficientes:
   módulos no coprimos compatibles, incompatibles y anidados.

Aceptación:

1. `RUSTC_WRAPPER= cargo test -q --all-features chinese_remainder` pasa.
2. `RUSTC_WRAPPER= cargo check -q --all-features --lib` no deja warnings por
   código muerto.

## Etapa 2. Hacer permanente el helper concordante

Objetivo: mover la composición concordante desde scaffolding de test a código
normal, sin hacerla pública.

Trabajo:

1. Compilar `class_group/concordant.rs` en build normal.
2. Mantener `compose_concordant_forms(...)` privado al submódulo de class group.
3. Aceptar internamente representantes primitivos positivos definidos del mismo
   discriminante, no necesariamente reducidos.
4. Documentar la fórmula concordante:

   ```text
   gcd(a,a′,(b+b′)/2) = 1
   A = aa′
   B ≡ b  (mod 2a)
   B ≡ b′ (mod 2a′)
   B² ≡ D (mod 4A)
   C = (B² − D)/(4A)
   ```

5. Reducir el resultado con `reduce_positive_definite()`.

Aceptación:

1. Tests con `D = -23` cubren identidad y un producto concordante no trivial.
2. Un par no concordante devuelve `NotConcordantForms`.
3. El resultado siempre satisface:
   - mismo discriminante;
   - primitiva;
   - reducida positiva definida.

## Etapa 3. Implementar equivalencia propia de formas

Objetivo: representar cambios unimodulares propios de variables sin mezclar esa
lógica con la composición.

Trabajo:

1. Crear un helper interno para aplicar una matriz

   ```text
   M = [[p,q],[r,s]], det(M) = 1
   ```

   a una forma:

   ```text
   f(x,y) -> f(px+qy, rx+sy)
   ```

2. Validar o construir la matriz de forma que `det(M) = 1`.
3. Testear que la transformación preserva:
   - discriminante;
   - primitividad;
   - positividad definida.
4. Testear transformaciones pequeñas contra ejemplos calculados a mano.

Aceptación:

1. La transformación vive en un archivo separado, por ejemplo
   `class_group/equivalence.rs`.
2. No usa búsqueda acotada como parte de la API.
3. No cambia la forma reducida salvo cuando se llama explícitamente desde una
   capa de composición.

## Etapa 4. Reducción efectiva a representantes concordantes

Objetivo: reemplazar una de las dos formas por una representante propiamente
equivalente que sea concordante con la otra, usando un algoritmo matemático
determinado, no una búsqueda con radio fijo.

Trabajo:

1. Elegir y documentar la versión clásica de reducción a concordancia que se va
   a implementar.
2. Implementarla como helper interno:

   ```rust
   fn concordant_representatives(
       &self,
       left: &BinaryQuadraticForm,
       right: &BinaryQuadraticForm,
   ) -> Result<(BinaryQuadraticForm, BinaryQuadraticForm), BinaryQuadraticFormError>
   ```

3. La salida debe satisfacer:

   ```text
   gcd(a,a′,(b+b′)/2) = 1
   ```

4. Agregar errores solo si son alcanzables desde inputs públicos válidos.
5. Evitar constantes mágicas. Si algún algoritmo requiere una cota derivada,
   documentarla y calcularla desde los coeficientes.

Aceptación:

1. Para todos los representantes reducidos enumerados de discriminantes chicos,
   la capa produce un par concordante o falla con un error matemáticamente
   justificado.
2. No hay `CONCORDANT_SEARCH_RADIUS` ni equivalente.
3. Los tests incluyen al menos un par inicialmente no concordante.

## Etapa 5. Exponer `QuadraticClassGroup::compose()`

Objetivo: agregar la API pública mínima.

Trabajo:

1. Implementar:

   ```rust
   pub fn compose(
       &self,
       left: &BinaryQuadraticForm,
       right: &BinaryQuadraticForm,
   ) -> Result<BinaryQuadraticForm, BinaryQuadraticFormError>
   ```

2. Validar `left` y `right` con `validate_reduced_member(...)`.
3. Convertirlos internamente a representantes concordantes.
4. Llamar al helper concordante.
5. Devolver el resultado reducido.

Rustdocs:

1. Decir qué hace:

   ```text
   compone dos representantes reducidos primitivos positivos definidos del
   discriminante del grupo y devuelve el representante reducido de la clase
   compuesta
   ```

2. No explicar extensamente lo que no hace.
3. Mencionar que usa composición clásica de Dirichlet/Gauss.
4. Usar una complejidad en `Θ(...)` solo si puede expresarse honestamente sin
   ruido; si no, omitirla hasta que el algoritmo esté estabilizado.

Aceptación:

1. `compose()` es la única API pública nueva para composición.
2. No se exponen helpers concordantes, matrices ni estrategias internas.
3. `RUSTC_WRAPPER= cargo check -q --all-features --lib` no deja warnings.

## Etapa 6. Tests de ley de grupo en discriminantes chicos

Objetivo: asegurar que la operación compone clases, no solo ejemplos sueltos.

Trabajo:

1. Usar `enumerate_reduced_forms()` para construir todos los representantes de
   discriminantes pequeños.
2. Testear cierre: el producto de dos formas enumeradas vuelve a estar en la
   lista.
3. Testear identidad:

   ```text
   principal * f = f
   f * principal = f
   ```

4. Testear inversos con `conjugate()`.
5. Testear asociatividad exhaustiva en discriminantes pequeños donde la lista
   sea corta.
6. Incluir ejemplos:
   - `D = -3`;
   - `D = -4`;
   - `D = -20`;
   - `D = -23`.

Aceptación:

1. Las leyes de grupo pasan exhaustivamente para esos discriminantes.
2. Los tests no dependen del orden accidental de helpers internos.
3. Si se compara contra una tabla esperada, la tabla se documenta en el test.

## Etapa 7. Limpieza de errores y visibilidad

Objetivo: dejar una API pública austera.

Trabajo:

1. Revisar `BinaryQuadraticFormError`.
2. Mantener públicos solo los errores alcanzables desde APIs públicas.
3. Bajar a privado o eliminar errores internos que hayan quedado obsoletos.
4. Revisar que los helpers de:
   - CRT compatible;
   - equivalencia propia;
   - concordancia;
   - composición concordante;

   no sean públicos.

Aceptación:

1. `rg "pub "` en el submódulo no muestra helpers accidentales.
2. `compose()` y los value objects principales son la única superficie estable
   nueva.

## Etapa 8. Ejemplo educativo opcional

Objetivo: mostrar la composición en un caso pequeño sin convertirlo en una
tabla gigante.

Trabajo:

1. Agregar un ejemplo que enumere las formas reducidas de `D = -23`.
2. Mostrar:
   - la forma principal;
   - una forma generadora;
   - uno o dos productos;
   - el resultado reducido.
3. Usar la capa de visualización si ya existe una `describe()` adecuada para
   formas o listas de formas.

Aceptación:

1. El ejemplo compila.
2. La salida no promete NUCOMP ni acción sobre curvas.
3. La narrativa queda enfocada en composición de clases de formas.

## Orden recomendado

1. Etapa 1: CRT compatible normal.
2. Etapa 2: helper concordante normal.
3. Etapa 3: equivalencia propia.
4. Etapa 4: reducción determinística a concordancia.
5. Etapa 5: `compose()` público.
6. Etapa 6: leyes de grupo.
7. Etapa 7: visibilidad/errores.
8. Etapa 8: ejemplo.

El punto crítico es la Etapa 4. Hasta que esa etapa tenga un algoritmo claro y
testeado, no conviene exponer `compose()` público.
