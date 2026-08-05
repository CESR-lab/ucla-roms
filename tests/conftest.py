"""Pytest configuration and shared fixtures for the ROMS test suite."""
import json
from pathlib import Path

import pytest

from . import _assertions
from ._assertions import computed_results, tolerated_mismatches
from ._roms_input_file_generation import create_roms_inputs


def pytest_addoption(parser):
    parser.addoption(
        "--environ",
        action="store",
        default="laptop",
        help="Reference-results environment name. "
             "The suite compares output hashes against tests/results/results_<environ>.json.",
    )
    parser.addoption(
        "--tolerate-hash-mismatch",
        action="store_true",
        default=False,
        help="Report reference-hash mismatches as XFAIL instead of failures, so they "
             "do not fail the run. Every other failure mode -- compile errors, ROMS "
             "crashes, missing output -- still fails. Used by CI, where the hashes are "
             "not reproducible across platforms; leave it off locally.",
    )


def pytest_configure(config):
    """Hand the --tolerate-hash-mismatch setting to the assertion helper.

    `assert_output_matches_reference` is called directly from test bodies with
    no access to the pytest config, so the flag is stashed on the module rather
    than threaded through all ten call sites.
    """
    _assertions.tolerate_hash_mismatch = config.getoption("--tolerate-hash-mismatch")


@pytest.fixture(scope="session")
def environ(request) -> str:
    return request.config.getoption("--environ")


@pytest.fixture(scope="session")
def reference_results(environ) -> dict:
    """Load the JSON of expected hashes for the chosen environment.

    Returns ``{}`` if the file does not exist, so a brand-new environment
    will fail loudly (with the computed hashes shown in the assertion message)
    rather than silently passing.
    """
    results_path = Path(__file__).parent / "results" / f"results_{environ}.json"
    if not results_path.exists():
        return {}
    with open(results_path) as f:
        return json.load(f)


@pytest.fixture(scope="session")
def input_dir(tmp_path_factory) -> Path:
    """One-time generation of every ROMS input file the suite needs.

    Built once per pytest session into a temporary directory and shared
    by every test via the shared inputs.
    """
    target = tmp_path_factory.mktemp("roms_inputs")
    create_roms_inputs(target)
    return target


def pytest_sessionfinish(session, exitstatus):
    """At the end of the session, dump the collected hashes as JSON.

    Only tests that actually ran appear in the dump, so a subset run
    produces a partial dict rather than entries with missing values.
    The dump goes to the terminal reporter so it isn't swallowed by
    pytest's output capture.
    """
    if not computed_results:
        return

    reporter = session.config.pluginmanager.get_plugin("terminalreporter")
    if reporter is None:
        return

    # Say plainly which mismatches were let through. On a run made green by
    # --tolerate-hash-mismatch this is the only unmissable record of it, since
    # nobody reads the per-test output of a passing job.
    if tolerated_mismatches:
        reporter.write_sep("=", "tolerated hash mismatches (reported XFAIL, did not fail the run)")
        for test_name in tolerated_mismatches:
            reporter.write_line(f"  {test_name}")

    current_environ = session.config.getoption("--environ")
    reporter.write_sep("=", f"computed results (paste into $ROMS_ROOT/results/results_{current_environ}.json to update if you understand the reason for - and expect - this discrepancy.)")
    reporter.write_line(json.dumps(computed_results, indent=2))
