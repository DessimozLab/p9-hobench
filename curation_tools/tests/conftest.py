import os
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
FIX = Path(__file__).parent / "fixtures"


def pytest_addoption(parser):
    parser.addoption("--network", action="store_true", default=False, help="run tests that call live web services")


def pytest_configure(config):
    config.addinivalue_line("markers", "network: test calls live web services (enable with --network or CURATION_NETWORK=1)")


def pytest_collection_modifyitems(config, items):
    if config.getoption("--network") or os.environ.get("CURATION_NETWORK") == "1":
        return
    skip = pytest.mark.skip(reason="network test (use --network)")
    for it in items:
        if "network" in it.keywords:
            it.add_marker(skip)


@pytest.fixture
def fix():
    return FIX
