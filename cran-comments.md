## R CMD check results

0 errors | 0 warnings | 0 notes (local); 0 errors | 0 warnings | 1 note
(win-builder, both R-release and R-devel -- see below)

## Test environments

* local: Windows 11 x64, R 4.6.1 (2026-06-24 ucrt) (devtools::check())
* win-builder (R-release, R 4.6.1): 0 errors, 0 warnings, 1 note
* win-builder (R-devel, R Under development (2026-09-03 r90483 ucrt)):
  0 errors, 0 warnings, 1 note
* GitHub Actions (windows-latest, macos-latest, ubuntu-latest x3 R
  versions): all passing

The win-builder note is the standard "New submission" notice plus a
list of words flagged by the spell-checker that are all proper nouns
(author names, e.g. Krakowski, Sankaran, Dileepkumar) or standard
statistical/distribution terminology and acronyms (e.g. Gompertz,
Lomax, Lindley, CAIC, HQIC, TTT) -- not actual misspellings.

## Downstream dependencies

This is a new package; there are no downstream dependencies.

## Notes to CRAN maintainers

This is a first submission. The package implements the "autorelevated"
family of distributions described in Krakowski (1973),
Dileepkumar and Sankaran (2022), and Dileep Kumar, Shabeer, and Sankaran
(2025); see `inst/CITATION` and the roxygen `@references` fields for
full details. The bundled `bladder_cancer` dataset reproduces the
widely-used Lee and Wang (2003) benchmark dataset already distributed
(in the same or equivalent form) in numerous other CRAN packages in
this field.

