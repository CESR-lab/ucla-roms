"""Realistic BGC tests (BEC and MARBL)."""
from ._assertions import assert_output_matches_reference
from ._helpers import (
    ADDITIONAL_FILES,
    BGC_REALISTIC_CPP_KEYS,
    MARBL_CPP_KEYS,
    OCEAN_PHYSICS_CPP_KEYS,
    REALISTIC_CPP_KEYS,
    UNIVERSAL_CPP_KEYS,
    ROMSConfiguration,
    create_test_namelist_dict,
    get_summary_value,
)


def test_bgc_real_bec(tmp_path, input_dir, reference_results):
    cpp_keys = (
        UNIVERSAL_CPP_KEYS + OCEAN_PHYSICS_CPP_KEYS
        + REALISTIC_CPP_KEYS + BGC_REALISTIC_CPP_KEYS
        + ["TIDES", "BIOLOGY_BEC2"]
    )
    conf = ROMSConfiguration(cpp_keys=cpp_keys, location=tmp_path / "build")
    conf.compile()

    nml = create_test_namelist_dict(input_dir)
    nml["TIME_STEPPING"]["dt"] = 20
    nml["TIME_STEPPING"]["ntimes"] = 10
    nml["BASIC_OUTPUT_SETTINGS"]["output_period_his"] = 100
    nml["BGC_SETTINGS"]["output_period_bgc_his"] = 400
    nml["PARAM_SETTINGS"]["nt_bgc"] = 26
    nml["TIDAL_FRC_SETTINGS"].update({"pot_tides": True, "ntides": 2})
    nml["FORCING_FILES"] = {
        "frcfiles": [
            str(input_dir / "example_input_boundary_forcing.nc"),
            str(input_dir / "example_input_surface_forcing.nc"),
            str(input_dir / "example_input_bgc_surface_forcing.nc"),
            str(input_dir / "example_input_co2_surface_forcing.nc"),
            str(input_dir / "example_input_river_forcing.nc"),
            str(input_dir / "example_input_bgc_boundary_forcing.nc"),
            str(input_dir / "example_input_tides.nc"),
        ]
    }

    conf.run(nml)

    pv = get_summary_value(conf.location, prefix="roms_his.20100101000000")
    bv = get_summary_value(conf.location, prefix="roms_bgc.20100101000000")
    assert_output_matches_reference(
        reference_results, "bgc_real_bec", {"physics": pv, "bgc": bv}
    )


def test_bgc_real_marbl(tmp_path, input_dir, reference_results):
    cpp_keys = (
        UNIVERSAL_CPP_KEYS + OCEAN_PHYSICS_CPP_KEYS
        + REALISTIC_CPP_KEYS + BGC_REALISTIC_CPP_KEYS
        + MARBL_CPP_KEYS + ["TIDES"] + ["SFLX_CORR", "CFLX_CORR"]
    )
    build_dir = tmp_path / "build"
    conf = ROMSConfiguration(cpp_keys=cpp_keys, location=build_dir)
    conf.compile()

    (build_dir / "marbl_in").symlink_to(ADDITIONAL_FILES / "marbl_in")

    nml = create_test_namelist_dict(input_dir)
    nml["BASIC_OUTPUT_SETTINGS"]["output_period_his"] = 100
    nml["BGC_SETTINGS"]["output_period_bgc_his"] = 100
    nml["MARBL_BIOGEOCHEMISTRY_SETTINGS"]["marbl_diagnostics_to_write"] = [
        "DCO2STAR_ALT_CO2", "graze_diat_zint", "FESEDFLUX",
    ]
    nml["MARBL_BIOGEOCHEMISTRY_SETTINGS"]["marbl_timestep"] = 20
    nml["TIME_STEPPING"].update({"dt": 20, "ntimes": 10})
    nml["TIDAL_FRC_SETTINGS"].update({"pot_tides": True, "ntides": 2})
    nml["PARAM_SETTINGS"]["nt_bgc"] = 32
    nml["SSS_CORRECTION"] = {"dSSSdt": 7.777}
    nml["DIC_ALK_CORRECTION"] = {"dCdt": 7.777}
    nml["FORCING_FILES"] = {
        "frcfiles": [
            str(input_dir / "example_input_boundary_forcing.nc"),
            str(input_dir / "example_input_surface_forcing.nc"),
            str(input_dir / "example_input_bgc_surface_forcing.nc"),
            str(input_dir / "example_input_restoring_sss.nc"),
            str(input_dir / "example_input_restoring_dic_alk.nc"),
            str(input_dir / "example_input_co2_surface_forcing.nc"),
            str(input_dir / "example_input_river_forcing.nc"),
            str(input_dir / "example_input_bgc_boundary_forcing.nc"),
            str(input_dir / "example_input_tides.nc"),
        ]
    }

    conf.run(nml)

    pv = get_summary_value(conf.location, prefix="roms_his.20100101000000")
    bv = get_summary_value(conf.location, prefix="roms_bgc.20100101000000")
    assert_output_matches_reference(
        reference_results, "bgc_real_marbl", {"physics": pv, "bgc": bv}
    )
