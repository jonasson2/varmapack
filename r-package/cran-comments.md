## Test environments

- macOS Sequoia 15.7.4 (aarch64), R 4.5.2
- Windows Server 2022 (x86_64-w64-mingw32), R 4.6.1 and R-devel,
  checked with Win-Builder

## R CMD check results

0 errors | 0 warnings | 1 note

The note reports that this is a new submission.

The Win-Builder R-devel check reports one warning during installation. It is
limited to two gfortran obsolescence messages from the bundled SLICOT source:
an old-style character-length declaration in `MB01RD` and a computed `GOTO`
in `SB04PX`. These pre-existing constructs remain in the current SLICOT 5.9.1
sources. All subsequent checks, including examples, tests, and vignette
rebuilding, pass.

## Reverse dependencies

There are currently no reverse dependencies on CRAN.

## Submission

This is a new submission.
