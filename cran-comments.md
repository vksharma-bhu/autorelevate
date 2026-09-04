## R CMD check results

0 errors | 0 warnings | 1 note

The note is:

    checking for future file timestamps ... NOTE
    unable to verify current time

This is an environment-level check (R CMD check trying to reach an NTP
time server) unrelated to package content, and is expected to be absent
on CRAN's own build machines.

## Test environments

* local: Windows 11 x64, R 4.5.2 (devtools::check())
* (add win-builder / R-hub results here if you run them before submitting)

## Downstream dependencies

This is a new package; there are no downstream dependencies.

## Notes to CRAN maintainers

This is a first submission. The package implements the "autorelevated"
(self-relevation) family of distributions described in Krakowski (1973),
Dileepkumar and Sankaran (2022), and Dileep Kumar, Shabeer, and Sankaran
(2025); see `inst/CITATION` and the roxygen `@references` fields for
full details. The bundled `bladder_cancer` dataset reproduces the
widely-used Lee and Wang (2003) benchmark dataset already distributed
(in the same or equivalent form) in numerous other CRAN packages in
this field.

