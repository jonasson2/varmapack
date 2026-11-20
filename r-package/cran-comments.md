## Test environments

- macOS Sequoia 15.7.4 (aarch64), R 4.5.2
- Ubuntu 26.04.1 (x86_64), R 4.6.2 (submission tarball)
- Windows Server 2022 x64, R-devel r90926 (UCRT), Win-Builder (submission tarball)

## R CMD check results

Local macOS `R CMD check --as-cran`: 0 errors | 0 warnings | 1 note.

On macOS, R could not verify the current time over the network. The installation,
compiled code, tests, vignettes, and PDF and HTML manuals passed their checks.

Linux `R CMD check --as-cran`: 0 errors | 0 warnings | 1 note.
The flagged options (`-Wdate-time`, `-Werror=format-security`, `-Wformat`)
come from the container's `R CMD config CFLAGS`, not the package.

Windows R-devel: 0 errors | 0 warnings | 1 note.
The spelling check flags Ornstein, Uhlenbeck, and Slepian in DESCRIPTION.
These are the surnames in the model names and are correctly spelled.

## Reverse dependencies

There are currently no reverse dependencies on CRAN.

## Submission

This is an update to varmapack (version 0.2.0).
