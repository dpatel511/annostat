# Releasing Annostat

## Release preparation

Use a focused branch and PR. Update `annostat/__init__.py`, `pyproject.toml`, the
release workflow's manual tag default, README version and dated changelog together.
Remove development-only wording for shipped changes; retain historical validation
records as records of their original audit. Document gene schema and calculation
policy changes separately from the package version.

Require green main CI after merge. Build/check distributions and exercise the
installed wheel outside the checkout, including
`examples/check_installed_gene_comparison.py`. This checks package assets and gene
reports rather than relying on imports/help alone. Public-data regression steps
are in [public-data-check.md](public-data-check.md).

## PyPI and GitHub

Publishing a GitHub release at the matching `vX.Y.Z` tag triggers
`.github/workflows/release.yml`. The workflow checks tag/version agreement, runs
the release tests, builds/checks wheel and source distributions, publishes through
the existing PyPI Trusted Publisher and attaches the distributions to GitHub.
Manual dispatch requires an existing tag and can retry publication. Do not create
a second release or upload a different artifact under an existing PyPI version.

Verify the Actions run completed, both distribution types exist on PyPI, and a
fresh installation outside the checkout reports the expected version and passes
`pip check`. A GitHub tag alone does not establish successful PyPI publication.

## Bioconda automation

Use Bioconda's existing autobump service rather than a duplicate Annostat workflow
with a cross-repository personal access token. Its scheduler checks upstream
releases and can open recipe-update PRs. Annostat's PyPI source URL is compatible
with that process. The recipe should explicitly retain:

```yaml
extra:
  autobump:
    enable: true
```

Keep `recipe-maintainers` alongside that block. Current bioconda-utils code uses
`extra.autobump`; older documentation also mentions `extra.watch` for subrecipes.
Annostat is a top-level recipe, so explicit enablement documents intent rather
than creating a new scheduler.

The normal sequence is PyPI publication, autobump detection/recipe PR, Bioconda
lint/build tests, merge, then channel publication. Scans may be delayed and are
not an immediate-release guarantee. Review Python/runtime dependency changes;
automated version/checksum updates do not establish their correctness. Do not
report Conda availability until the built version is actually in the channel.

If an existing maintainer PR already covers Annostat, update that PR to the latest
published version instead of creating duplicates. For a manual fallback, update
`recipes/annostat/meta.yaml` in the separate Bioconda fork: set the version, reset
the build number to 0, download the actual released sdist, verify its SHA-256
against PyPI, and check requirements against the released package metadata before
submission. Bioconda build/merge policy remains upstream's responsibility.

References: [Bioconda updating recipes](https://bioconda.github.io/contributor/updating.html),
[CI inventory](https://bioconda.github.io/developer/ci-inventory.html), and
[current autobump filters](https://github.com/bioconda/bioconda-utils/blob/master/bioconda_utils/autobump.py).
