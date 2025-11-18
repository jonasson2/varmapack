---
name: varmapack-style
description: Coding standards and project architecture for the varmapack C library (varmasim). Use when working with varmapack source code, refactoring, writing new functions, or any C development tasks within the varmapack project. Enforces strict ISO C11 style with 2-space indentation, Stroustrup braces, and varmapack_ prefixed public API.
---

# Varmapack Style

Development guidelines for the varmapack C library (varmasim).

## C Coding Style (Strict)

- **Standard**: ISO C11 only. No compiler extensions or attributes.
- **Indentation**: 2 spaces (no tabs).
- **Comments**: `//` only.
- **Pointer style**: `int *x` (space before pointer).
- **Brace style**: Stroustrup; `else` on its own line.
- **Max line width**: 90 characters.
- **Loop increment**: Prefer `i++` to `++i` in for-loops.
- **Multiplication/Division**: No spaces around `*` or `/` operators.
- **Boolean**: Use `bool`/`true`/`false`; prefer `0` to `NULL`.
- **Integer types**: Use `int` for all integers. Never use `size_t`.
- **Public API**: `snake_case` with `varmapack_` prefix.
- **Internal helpers**: `static` or `static inline`, use `CamelCase`, no prefix.
- **Local variables**: Prefer short or shortish names.
- **Bool control variables**: Use ALL_UPPERCASE.
- **No min/max macros**: Use `fmin`/`fmax` for doubles; define `static inline` `imin`/`imax` for ints.
- **For-loops**: Use compact expressions with no spaces: `for (i=0; i<n; i++)`.
- **Macros**: Only for include guards. No other macros.
- **Line length**: Maximum line length 90 characters.

## Project Architecture

- **Public API**: Declared in `varmapack.h` only (with C++ guards).
- **Internal code**: Lives in `src/` with internal headers (`static inline` only).
- **No visibility attributes**: Library is pure C.
- **Directory layout**:
  - `src/` - implementation + internal headers
  - `tests/` - unit tests
  - `examples/` - sample programs
- **Build system**: Meson.
- **Exports**: Only `varmapack_*` symbols are exported.

## Project Status

- R interface is on hold.
- Goal: Clean, portable C library with consistent naming and internal structure.

## Code Examples

Function definition (public API):
```c
// Calculate mean of array
double varmapack_mean(const double *data, int n) {
  if (n == 0) {
    return 0.0;
  }
  double sum = 0.0;
  for (int i=0; i<n; i++) {
    sum += data[i];
  }
  return sum/n;
}
```

Function definition (internal helper with CamelCase):
```c
// Internal helper for bounds checking
static inline bool IsValidIndex(int idx, int max) {
  return idx < max;
}
```

Bool control variable (ALL_UPPERCASE):
```c
bool CONVERGED = false;
bool VERBOSE = true;
```

Conditional with Stroustrup bracing:
```c
if (condition) {
  // action
}
else {
  // alternative
}
```

Compact for-loop with short variable names:
```c
for (int i=0; i<n; i++) {
  for (int j=0; j<m; j++) {
    matrix[i*m+j] = i*j;
  }
}
```

## Output Requirements

- Follow style exactly: 2-space indentation, Stroustrup braces, `//` comments.
- All public functions must use `varmapack_` prefix with snake_case.
- All internal helpers must be `static` or `static inline` with CamelCase.
- Use `int` for all integers; never `size_t`.
- Bool control variables in ALL_UPPERCASE.
- Short variable names for locals.
- Compact for-loops: `for (i=0; i<n; i++)`.
- No spaces around `*` or `/`.
- Keep lines under 90 characters.
- Tone: precise, technical, concise.
