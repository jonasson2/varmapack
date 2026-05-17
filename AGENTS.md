Always obey the rules below.

---------------------------
CORE C STYLE
---------------------------
- Platform agnostic; should compile on ARM and Intel with gcc, clang,
  MSVC (cl.exe) and icx
- Indent 2 spaces
- No blank lines inside functions (allowed between code blocks outside functions)
- Stroustrup braces
- Open brace ({) on same line as function definition
- else on a new line
- // comments.
- Pointer form: int *x. Prefer i++ in loops.
- Max line length 90; no spaces around * or /, spaces around comparison operators
- Use bool/true/false; prefer 0 over NULL and '\0'
- Don't use any casts unless necessary
- Prefer brackets around sizeof argument
- Never use const unless absolutely necessary
- Use "int" for index and size variables when bounds are known to be small
- Always use int constants for whole number doubles and let C convert them to double

---------------------------
ARCHITECTURE
---------------------------
- Public header: varmapack.h
- Public C API uses the varmapack_ prefix
- randompack is an external dependency, not an internal sublibrary
- Sources under src/, tests live in tests/, examples in examples/
- R package files live in r-package/
- Do not manually edit generated R source copies once sync scripts exist
- In tests, use utilities declared/defined in Tests.h, xCheck.h, and ExtraUtil.h
- Meson/Ninja drive builds
- Build and verify newly written code in the `release` Meson directory by default
- Julia and Fortran interfaces are out of scope for now

---------------------------
EDITING
---------------------------
- Before editing a file, check whether it has changed and reread it from disk.
  Do not rely on stale context, since user refactoring may have happened between
  Codex turns.

---------------------------
OUTPUT
---------------------------
- Keep responses concise, technical, and style-compliant.
- Do not restate these instructions.
- Provide concrete answers or code following the rules above.
- Stroustrup braces also on function signature lines
- Use 0 instead of NULL
- No unnecessary casts
- Prefer 0, 1, 2, 3 to 0.0, 1.0, 2.0, 3.0 in double initializations.
