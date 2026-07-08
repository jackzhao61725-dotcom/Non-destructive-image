from __future__ import annotations

from pathlib import Path

import pytest

np = pytest.importorskip("numpy")

from non_destructive_image import simulate_faraday_image


BASELINE_PATH = Path("regression/baseline/imaging/faraday_imaging_baseline_v1.npz")


def test_simulate_faraday_image_matches_faraday_baseline() -> None:
    with np.load(BASELINE_PATH) as baseline:
        result = simulate_faraday_image(
            baseline["theta_f_map_rad"],
            baseline["pupil"],
            return_intermediates=True,
        )

        for key in [
            "theta_f_map_rad",
            "sigma_plus_object_field",
            "sigma_minus_object_field",
            "sigma_plus_field",
            "sigma_minus_field",
            "output_ex_field",
            "output_ey_field",
            "dark_field_intensity",
            "dual_port_u_intensity",
            "dual_port_v_intensity",
            "dual_port_signal",
        ]:
            np.testing.assert_allclose(result[key], baseline[key], rtol=1e-10, atol=1e-12)
