# Releasing autoboot

`VERSION` is the canonical software version. Stable releases follow
[Semantic Versioning](https://semver.org/): backward-compatible fixes increment
PATCH, backward-compatible features increment MINOR, and incompatible changes
to the documented public API increment MAJOR.

The documented group constructors, representation labels, operator registration,
crossing-equation and export interfaces form the public API. Private
implementation details do not.

## Release procedure

1. Update `VERSION`, the software `version` in `CITATION.cff`, and
   `.zenodo.json` together. Set `date-released` in `CITATION.cff`.
2. Run `make test` and `python3 test/check-version.py --tag vX.Y.Z`.
   The full suite requires an activated Wolfram Engine; record any checks
   that could not be completed in the release PR.
3. Merge the release preparation PR into `master` after reviewing CI results.
4. Create the `vX.Y.Z` tag on the merged commit and publish a GitHub release.
   Published tags must not be moved. Tag CI checks metadata consistency.
5. With Zenodo's GitHub integration enabled, check the resulting software
   record, author ORCID, version DOI, and concept DOI. Use the concept DOI
   for a badge representing all versions; use a version DOI to identify a
   particular archived release.

`CITATION.cff` supplies GitHub's software citation information and references the
accompanying paper. `.zenodo.json` supplies the software metadata for Zenodo
archiving and takes precedence over `CITATION.cff` there. Keep both consistent.
