# Versioned protein-level mapping profiles

Each JSON file declares one vendor export header set as `profile@version`.
The configuration names the exact profile (`input.profile`, for example
`maxquant_proteingroups@2.x-lfq-headers-1`) and the profile file
(`input.mapping`). The adapter never infers a vendor or version from
similar-looking names. Any header that is neither required, a declared
annotation, an ignored pattern nor a declared intensity column is
rejected with `E_MAPPING_HEADER`. Peptide, precursor or site columns are
rejected with `E_UNSUPPORTED_SCOPE`.

These profiles were authored from public documentation of the header
sets and are qualified only against tiny synthetic exports (V016). They do
not claim compatibility with any specific vendor release. A new header
set needs a new versioned profile.

`legacy-crosswalk-method-a.json` and `legacy-crosswalk-method-b.json` are
the semantic Method A/Method B crosswalk templates for the historical
workbook importer (V017). Contrast identity is the biological
numerator/denominator, never the D1/T1 position.
