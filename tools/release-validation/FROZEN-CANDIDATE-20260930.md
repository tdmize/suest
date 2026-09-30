# Release-validation candidate — 2026-09-30

Candidate version: 0.1.6.9000. Parent baseline: 29ae851.
This is a development candidate for platform checks, not a published release.

The four validation stages pass their declared numerical targets. The latest
Linux/R 4.5.2 check passes 1,680 source and installed-package assertions and
144 acceptance cases, with R CMD check Status: OK. The full reports specify
scope and retained numerical qualifications. Platform CI is still pending.

The freeze adds workflow_dispatch to R-CMD-check, pull_request to the existing
acceptance workflow, and a dependency preflight that requires every declared
Imports/Suggests package and records installed versions. The YAML and preflight
were checked locally. Executable statistical code, package tests, and built
help match the gate-4 checked source. No statistical suite rerun is claimed
for these workflow-only additions.

Run both workflows against this candidate commit:

1. R-CMD-check: macOS and Windows R release; Linux R release and R-devel.
2. paper-replications: the full acceptance suite on Linux.

Use a candidate branch and a pull request targeting main, or dispatch both
workflows manually with that branch selected. Push of this candidate branch
alone does not run the main-only push workflows; select the branch explicitly.
The pkgdown deployment should wait until the validated candidate reaches main.

For each job retain its commit SHA, R/dependency versions, package-check log,
and test totals. Require no errors, failed assertions, or unexplained skips.
Inspect every warning/note, even when Actions is green: R-CMD-check currently
fails on errors, while warnings still require review. Platform test counts
must be reconciled with the local 1,680 assertions; a count difference alone
needs investigation rather than an automatic pass/fail claim. Acceptance
should complete 144 cases with zero failures.

Preserve the FE average-level inference exclusion, native integration
conventions, and random-slope Stata GSEM qualifications, including the retained
failed absolute-tolerance diagnostic. Reviewed native fitting warnings in
audit scripts are documented separately from the clean package-check results.

Nothing has been pushed, merged, tagged, or published. GitHub push requires
Trent's explicit authorization. After platform results pass review, choose
the release version and verify that any release-metadata changes preserve the
tested source before tagging or website promotion.
