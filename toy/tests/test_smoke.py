"""Phase 0 smoke test: the package and its dependencies import."""


def test_imports() -> None:
    import numpy  # noqa: F401
    import pygame  # noqa: F401
    import scipy  # noqa: F401

    import sandframe_toy

    assert sandframe_toy.__version__ == "0.0.0"
