"""Helpers for asserting test output hashes match the reference JSON.

Every call to :func:`assert_matches_reference` records its computed values
into the module-level :data:`computed_results` dict, regardless of whether
the assertion passed or failed. ``conftest.py``'s ``pytest_sessionfinish``
hook then dumps the collected dict at the end of the run, so you can paste
it straight into ``results/results_<environ>.json`` to update references.
"""
import pytest

# Session-wide collector, populated by `assert_matches_reference` and read
# by the `pytest_sessionfinish` hook in conftest.py. Only tests that actually
# ran will appear here, so running a subset of the suite produces a partial
# dict rather than entries with missing values.
computed_results: dict[str, dict] = {}

# Set from the `--tolerate-hash-mismatch` command-line flag by
# `pytest_configure` in conftest.py. When true, a hash mismatch is reported as
# an expected failure (XFAIL) instead of a failure, so the run stays green.
# Every other way a test can fail -- a compile error, a ROMS crash, a missing
# output file -- is untouched and still fails the run.
#
# This exists because the reference hashes are not reproducible across
# platforms, so hash mismatches are chronic noise in CI while genuine
# regressions are not. Defaults to false: an ordinary local `pytest` run still
# treats a mismatch as a failure. See the CI workflow for where it is enabled.
tolerate_hash_mismatch: bool = False

# Names of the tests whose mismatches were tolerated this session, reported by
# the `pytest_sessionfinish` hook so a green run still says plainly what was
# let through.
tolerated_mismatches: list[str] = []


def assert_output_matches_reference(reference_results: dict, test_name: str, computed: dict):
    """Compare ``computed`` against ``reference_results[test_name]``.

    ``computed`` is a dict like ``{"physics": "<hash>"}`` or
    ``{"physics": "<hash>", "bgc": "<hash>"}``.

    The computed dict is recorded into :data:`computed_results` before the
    comparison is made, so the session-end dump captures every test that
    ran (not just the ones that failed).

    When :data:`tolerate_hash_mismatch` is set, a mismatch calls
    ``pytest.xfail`` rather than raising, which ends the test immediately and
    reports it as XFAIL. As with the ``AssertionError`` below, nothing after
    this call in the test body runs, so a test is only ever recorded once.
    """
    computed_results[test_name] = computed

    expected = reference_results.get(test_name)
    if expected != computed:
        message = (
            f"Hash mismatch for {test_name!r}.\n"
            f"  expected: {expected}\n"
            f"  computed: {computed}"
        )
        if tolerate_hash_mismatch:
            tolerated_mismatches.append(test_name)
            pytest.xfail(message)
        raise AssertionError(message)
