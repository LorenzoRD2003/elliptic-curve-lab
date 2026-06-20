# Plan Por Etapas: `GeneralWeierstrassCurve<F>`

## Objetivo

Agregar un nuevo modelo algebraico

`GeneralWeierstrassCurve<F>`

para curvas de la forma

`y^2 + a1*x*y + a3*y = x^3 + a2*x^2 + a4*x + a6`

sin desarmar prematuramente la infraestructura existente alrededor de

`ShortWeierstrassCurve<F>`.

La estrategia recomendada es:

- introducir `GeneralWeierstrassCurve<F>` como un modelo hermano
- reutilizar los traits genéricos ya existentes
- delegar las primeras capacidades avanzadas a una reducción explícita hacia
  `ShortWeierstrassCurve<F>` cuando esa reducción sea honesta
- posponer la generalización profunda de subsistemas short-specific hasta que el
  nuevo modelo ya sea útil y estable

## Resumen Del Estado Actual

El repo ya tiene una buena base abstracta para convivir con más de un modelo:

- `CurveModel`
- `AffineCurveModel`
- `GroupCurveModel`
- `LiftXCoordinate`
- `EnumerableCurveModel`
- `FiniteGroupCurveModel`
- `RelativeFrobeniusCurveModel`
- `FrobeniusTraceCurveModel`

Esos traits ya encapsulan bastante comportamiento reusable para curvas afines y
grupos pequeños.

Sin embargo, gran parte del valor actual vive en inherent methods y módulos
directamente acoplados a `ShortWeierstrassCurve<F>`:

- invariantes y discriminante
- ley de grupo
- isomorfismos y twists
- división polinómica
- Schoof
- function fields
- Vélu e isogenias short-specific
- visualización detallada
- `proptest_support`
- ejemplos

Conclusión práctica:

- sí conviene reutilizar la capa de traits
- no conviene intentar convertir toda la infraestructura actual en una capa
  completamente genérica en la primera iteración

## Principios De Diseño

1. `GeneralWeierstrassCurve<F>` debe entrar como tipo nuevo, no como refactor
   destructivo de `ShortWeierstrassCurve<F>`.
2. La primera versión debe priorizar corrección y compatibilidad conceptual,
   no cobertura total de features.
3. La reducción a short debe ser explícita y rastreable, no una conversión
   implícita escondida.
4. Los algoritmos que solo son honestos en el modelo corto deben seguir
   declarándose como tales, aunque temporalmente se usen desde el modelo
   general a través de un puente.
5. La documentación debe distinguir con claridad:
   - capacidades nativas del modelo general
   - capacidades obtenidas vía companion short

## Qué Conviene Preservar Desde El Día Uno

Idealmente, una vez que el modelo general implemente los traits correctos,
deberían seguir funcionando temprano:

- validación de puntos afines
- punto identidad
- pertenencia al modelo
- suma, resta, duplicación y multiplicación escalar
- enumeración exhaustiva sobre campos finitos pequeños
- orden de punto por recorrido exhaustivo
- orden de grupo por enumeración o por rutas ya abstractas
- búsquedas de Hasse basadas en `GroupCurveModel`
- kernels explícitos e `IsogenyKernel::cyclic`

En cambio, conviene dejar para etapas posteriores:

- división polinómica general
- Schoof nativo para el modelo general
- function fields del modelo general
- Vélu nativo sobre ecuación general
- graph layer con testigos propios
- visualización exhaustiva al nivel short actual

## Etapa 0: Preparación E Inventario

### Objetivo

Preparar el terreno sin tocar todavía la semántica del modelo corto.

### Trabajo

- crear el módulo hermano:
  - `src/elliptic_curves/models/general_weierstrass/`
- exportarlo desde:
  - `src/elliptic_curves/models/mod.rs`
  - `src/elliptic_curves/mod.rs`
- revisar nombres de error y docs para evitar que todo siga sonando
  short-specific cuando el error en realidad sea más general

### Resultado esperado

El crate puede compilar con ambos nombres de modelo expuestos aunque el nuevo
tipo todavía tenga una superficie pequeña.

## Etapa 1: Núcleo Del Modelo General

### Objetivo

Tener un descriptor algebraico válido y autocontenido para la ecuación general.

### Trabajo

- definir:
  - `pub struct GeneralWeierstrassCurve<F: Field>`
- almacenar:
  - `a1, a2, a3, a4, a6`
- implementar:
  - constructor validado
  - accessors para coeficientes
  - formateo de ecuación
  - discriminante general
  - invariantes clásicos del modelo general (`b2`, `b4`, `b6`, `b8`, `c4`,
    `c6`, `j`)
- documentar cuidadosamente el criterio de no singularidad usado

### Notas

- esta etapa no debe asumir `char(F) != 2,3)` como restricción global
- la restricción a característica distinta de `2` y `3` debe aparecer solo
  donde sea necesaria para pasar a short

### Resultado esperado

El modelo general ya representa honestamente una curva elíptica afín en forma
de Weierstrass general.

## Etapa 2: Reducción Explícita A Short

### Objetivo

Construir el puente más importante del proyecto:

`GeneralWeierstrassCurve<F> -> ShortWeierstrassCurve<F>`

cuando la característica lo permita.

### Trabajo

- crear un objeto explícito, por ejemplo:
  - `GeneralWeierstrassReduction<F>`
  - o `GeneralToShortWeierstrass<F>`
- ese objeto debe guardar:
  - la curva general original
  - la curva short companion
  - los datos del cambio de variables
  - las funciones de transporte de puntos en ambos sentidos
- implementar algo como:
  - `try_as_short_weierstrass()`
  - `reduction_to_short()`

### Recomendación de diseño

No devolver solo la curva short. Conviene devolver una estructura que preserve
el cambio de coordenadas como dato de primer orden. Eso evita perder
información y permite delegar ley de grupo y validaciones con ida y vuelta
explícita.

### Resultado esperado

El proyecto gana un “companion short” verificable para una gran clase de
curvas generales, reutilizando gran parte del trabajo actual sin duplicarlo.

## Etapa 3: Traits Básicos Del Modelo General

### Objetivo

Hacer que `GeneralWeierstrassCurve<F>` participe de la infraestructura genérica
del crate.

### Trabajo

- implementar `CurveModel`
- implementar `AffineCurveModel`
- implementar `HasJInvariant`
- evaluar si conviene implementar `LiftXCoordinate`

### Nota importante sobre `LiftXCoordinate`

Para la ecuación general, fijado `x`, la variable `y` sigue satisfaciendo una
ecuación cuadrática. Por lo tanto:

- conceptualmente sí puede haber una ruta de “lift de `x`”
- pero ya no es la simple historia `y^2 = rhs(x)`

Mi recomendación es:

- no forzar `LiftXCoordinate` en la primera etapa salvo que la API se extienda
  con honestidad
- si se quiere soportar enumeración exhaustiva pronto, introducir un trait
  nuevo más general para “resolver `y` a partir de `x`” o implementar un helper
  interno específico del modelo general

### Resultado esperado

El nuevo modelo deja de ser un simple contenedor y empieza a integrarse con la
arquitectura de curvas del crate.

## Etapa 4: Ley De Grupo Reutilizando El Companion Short

### Objetivo

Conseguir una ley de grupo correcta pronto, evitando introducir fórmulas
afines generales complejas demasiado temprano.

### Trabajo

- implementar `GroupCurveModel` para `GeneralWeierstrassCurve<F>`
- usar la reducción explícita a short en cada operación:
  - validar puntos en el modelo general
  - transportar al companion short
  - sumar o duplicar ahí
  - transportar de vuelta

### Atención especial

La negación ya no es `(x, -y)`. Para la ecuación general debe codificarse la
involución correcta del modelo.

Esto es importante porque hoy `AffinePoint<F>::neg()` está documentado como
conveniente para modelos simétricos bajo `y -> -y`, y no debe reutilizarse como
semántica de negación del modelo general.

### Resultado esperado

El modelo general obtiene:

- suma
- resta
- duplicación
- multiplicación escalar
- chequeos de torsión

sin necesidad de una segunda implementación algebraica pesada desde el día uno.

## Etapa 5: Compatibilidad Con La Capa Genérica De Grupos Finitos

### Objetivo

Hacer que el modelo general reciba gratis la mayor cantidad posible de helpers
ya existentes.

### Trabajo

Una vez que el modelo general tenga:

- `CurveModel`
- `AffineCurveModel`
- `GroupCurveModel`
- y una ruta honesta para obtener puntos desde cada `x` en campos pequeños

entonces se puede integrar con:

- `EnumerableCurveModel`
- `FiniteGroupCurveModel`
- `FrobeniusTraceCurveModel`
- búsquedas Hasse internas
- `IsogenyKernel`

### Recomendación

Esta etapa puede apoyarse temporalmente en helpers internos delegados al
companion short, pero las interfaces públicas deben seguir hablando en términos
de `GeneralWeierstrassCurve<F>`.

### Resultado esperado

Empiezan a “seguir funcionando” varias superficies educativas ya valiosas sin
duplicar algoritmos.

## Etapa 6: Wrappers Ergonómicos De Alto Valor

### Objetivo

Exponer una superficie útil para usuarios del crate sin prometer todavía que
todo es nativo del modelo general.

### Trabajo

Agregar en `GeneralWeierstrassCurve<F>` solo los wrappers con mejor retorno:

- `group_order_by(...)`
- `group_order_by_small_field(...)`
- `point_order_by(...)`
- `group_exponent_by(...)`
- `frobenius_trace_by(...)`

### Regla

Cada wrapper debe documentar si:

- el cálculo es nativo del modelo general
- o si se hace vía reducción a `ShortWeierstrassCurve<F>`

### Resultado esperado

Los callers obtienen APIs familiares y consistentes con el resto del repo.

## Etapa 7: Test Suite De Compatibilidad

### Objetivo

Blindar el diseño antes de expandirlo a subsistemas más complejos.

### Tests mínimos recomendados

- construcción válida e inválida
- singularidad detectada correctamente
- invariantes del modelo general consistentes con la reducción a short
- roundtrip de puntos general -> short -> general
- `contains` preservado por el transporte
- suma compatible con transporte
- duplicación compatible con transporte
- multiplicación escalar compatible con transporte
- orden de punto compatible con transporte en campos pequeños
- orden de grupo compatible con transporte en campos pequeños

### Property tests recomendados

- dos puntos aleatorios sobre una curva reducible a short
- comparar `P + Q` en ambos modelos
- comparar `nP` en ambos modelos
- comparar `j`

### Resultado esperado

La reducción deja de ser una idea informal y pasa a ser una pieza certificada
del diseño.

## Etapa 8: Visualización, Ejemplos Y `proptest_support`

### Objetivo

Volver visible el nuevo modelo en la experiencia educativa del repo.

### Trabajo

- agregar un formatter compacto y un descriptor pedagógico para el modelo
  general
- sumar ejemplos en `examples/`
- crear generadores de curvas generales no singulares en `proptest_support`
- mostrar al menos un ejemplo donde:
  - se vea la ecuación general
  - se vea la reducción a short
  - se compare un cálculo en ambos modelos

### Resultado esperado

El nuevo modelo deja de ser una pieza “solo interna” y entra en la narrativa
del proyecto.

## Etapa 9: Generalización Profunda Opcional

### Objetivo

Decidir, con evidencia, qué partes del stack deben volverse realmente
model-agnostic y cuáles deben seguir siendo short-specific con puente.

### Candidatos a evaluar más adelante

- ley de grupo nativa en coordenadas generales
- function fields del modelo general
- isomorfismos generales
- Vélu o isogenias nativas del modelo general
- graph layer con witness propio
- división polinómica general
- Schoof nativo general

### Recomendación

No abrir esta etapa hasta que las etapas 1 a 8 estén estabilizadas. Antes de
eso, el riesgo de abstraer de más es alto.

## Orden Recomendado De Implementación

1. módulo nuevo y export surface
2. tipo y constructor del modelo general
3. invariantes y discriminante
4. objeto explícito de reducción a short
5. `CurveModel`, `AffineCurveModel`, `HasJInvariant`
6. `GroupCurveModel` vía companion short
7. compatibilidad con grupos finitos pequeños
8. wrappers ergonómicos
9. tests de roundtrip y compatibilidad
10. visualización, ejemplos y generators
11. generalización profunda, solo si sigue teniendo sentido

## Riesgos Principales

### Riesgo 1: abstraer demasiado pronto

Intentar sacar una “supertrait universal de curvas Weierstrass” antes de tener
el modelo general funcionando probablemente va a dispersar mucho el esfuerzo.

Mitigación:

- priorizar tipo nuevo + puente explícito

### Riesgo 2: semántica equivocada de negación

Reutilizar la intuición `P -> (x, -y)` sería incorrecto para el modelo general.

Mitigación:

- implementar la negación del modelo general como operación propia

### Riesgo 3: esconder demasiada delegación al modelo short

Si todo se delega de forma implícita, después cuesta distinguir qué es
matemáticamente nativo del modelo general.

Mitigación:

- documentar explícitamente los wrappers “via short companion”

### Riesgo 4: arrastrar restricciones innecesarias de característica

El modelo general no debería nacer limitado artificialmente por restricciones
propias del modelo corto.

Mitigación:

- separar validación del modelo general de validación de la reducción a short

## Criterio De Éxito Del Primer Milestone

Consideraría exitoso el primer milestone si el repo llega a tener:

- `GeneralWeierstrassCurve<F>` expuesto públicamente
- construcción y validación honestas
- invariantes clásicos básicos
- reducción explícita a `ShortWeierstrassCurve<F>` cuando corresponda
- ley de grupo correcta vía esa reducción
- tests de compatibilidad entre ambos modelos
- al menos un ejemplo educativo de ida y vuelta

Eso ya convertiría al nuevo modelo en una pieza real del proyecto sin exigir
la migración masiva de todo el árbol short-specific.

## Recomendación Final

La mejor estrategia para este repo es:

- agregar primero un modelo general útil
- usar `ShortWeierstrassCurve<F>` como backend algebraico reutilizable cuando
  la reducción exista
- generalizar subsistemas uno por uno solo después de que el nuevo modelo
  demuestre valor real

En otras palabras:

primero coexistencia y compatibilidad;
luego ergonomía;
recién después, generalización profunda.
