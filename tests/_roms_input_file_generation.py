"""ROMS input files built from synthetic upstream data via roms-tools.

The entry point ``create_roms_inputs`` first lays down the upstream synthetic
datasets in ``target_dir`` (via :mod:`upstream_data`), then uses ``roms-tools``
to produce all the inputs the test suite needs: grid, initial conditions,
boundary forcing (physics + BGC), surface forcing (physics + BGC), tides,
CDR forcing in three flavours, river/pipe forcing, and a flux-form surface
forcing file derived from the bulk surface forcing.
"""
import datetime as dt
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import roms_tools as rt
import xarray as xr

from ._roms_tools_input_file_generation import create_roms_tools_inputs


def create_roms_inputs(target_dir: Path):
    """Populate ``target_dir`` with every ROMS input file the tests need."""
    create_roms_tools_inputs(target_dir)

    grid = create_roms_grid()
    grid.save(target_dir / "example_input_grid.nc")

    create_roms_physical_boundary_forcing(grid, target_dir)
    create_roms_bgc_boundary_forcing(grid, target_dir)

    create_roms_physical_surface_forcing(grid, target_dir)
    create_roms_bgc_surface_forcing(grid, target_dir)
    create_roms_surface_forcing_restoring_sss(grid, target_dir)
    create_roms_surface_forcing_restoring_dic_alk(grid, target_dir)

    create_roms_co2_surface_forcing(grid, target_dir)

    create_roms_initial_conditions(grid, target_dir)

    create_roms_tides(grid, target_dir)

    create_roms_cdr_forcing_3d(
        grid,
        irange=range(8, 12),
        jrange=range(18, 22),
        alk_pert=1e-3,
        target_dir=target_dir,
    )

    create_roms_cdr_forcing_dp(
        grid,
        release_lon=-120.64,
        release_lat=34.5349,
        alk_pert=1e-3,
        target_dir=target_dir,
    )

    create_roms_cdr_forcing_parm(
        grid,
        release_lon=-120.64,
        release_lat=34.5349,
        release_dep=10,
        alk_pert=1e-3,
        target_dir=target_dir,
    )

    river_forcing_ds = create_river_or_pipe_data(
        input_type="river",
        forcing_i=[10, 11, 12, 13, 14],
        forcing_j=[20, 19, 19, 18, 18],
        grid=grid,
    )
    river_forcing_ds.to_netcdf(target_dir / "example_input_river_forcing.nc")

    pipe_forcing_ds = create_river_or_pipe_data(
        input_type="pipe",
        forcing_i=[3, 4, 5],
        forcing_j=[20, 19, 18],
        grid=grid,
    )
    pipe_forcing_ds.to_netcdf(target_dir / "example_input_pipe_forcing.nc")

    create_flux_forcing_from_bulk(target_dir / "example_input_surface_forcing.nc")

    # Pre-processing: split each input file across the 2x2 MPI decomposition.
    subprocess.run(
        "for F in example_input_*;do partit 2 2 ${F};done",
        text=True, shell=True, cwd=target_dir,
    )
    subprocess.run(
        "partit 2 2 cdr_forcing_3d.nc",
        text=True, shell=True, cwd=target_dir,
    )


def create_roms_grid():
    return rt.Grid(
        theta_s=6.0,
        theta_b=6.0,
        hc=25.0,
        N=10,
        nx=39, ny=19,
        center_lon=-120.6475,
        center_lat=34.54,
        rot=-37.4,
        size_x=25,
        size_y=12,
    )


def create_roms_physical_boundary_forcing(grid, target_dir: Path):
    rpf = rt.BoundaryForcing(
        grid=grid,
        start_time=dt.datetime(2010, 1, 1),
        end_time=dt.datetime(2010, 1, 2),
        source={"name": "GLORYS", "path": target_dir / "fake_phys_3d_data.nc"},
        apply_2d_horizontal_fill=True,
    )
    rpf.save(target_dir / "example_input_boundary_forcing.nc", group=False)


def create_roms_bgc_boundary_forcing(grid, target_dir: Path):
    rbf = rt.BoundaryForcing(
        grid=grid,
        type="bgc",
        start_time=dt.datetime(2010, 1, 1),
        end_time=dt.datetime(2010, 1, 2),
        source={"name": "CESM_REGRIDDED",
                "path": target_dir / "fake_bgc_3d_data.nc"},
        apply_2d_horizontal_fill=True,
    )
    # Add DOFE for BEC:
    for bry in ["west", "north", "south", "east"]:
        rbf.ds[f"DOFE_{bry}"] = rbf.ds[f"DOP_{bry}"]
        rbf.ds[f"DOFE_{bry}"][:] = 0.
        rbf.ds[f"DOFE_{bry}"].attrs = {
            "long_name": f"{bry}ern boundary dissolved organic iron",
            "units": "mMol Fe m-3",
        }

    rbf.save(target_dir / "example_input_bgc_boundary_forcing.nc", group=False)


def create_roms_physical_surface_forcing(grid, target_dir: Path):
    sf = rt.SurfaceForcing(
        grid=grid,
        start_time=dt.datetime(2010, 1, 1),
        end_time=dt.datetime(2010, 1, 2),
        source={"name": "ERA5", "path": target_dir / "fake_phys_surf_data.nc"},
    )
    sf.save(target_dir / "example_input_surface_forcing.nc", group=False)


def create_roms_bgc_surface_forcing(grid, target_dir: Path):
    rsf = rt.SurfaceForcing(
        grid=grid,
        type="bgc",
        start_time=dt.datetime(2010, 1, 1),
        end_time=dt.datetime(2010, 1, 2),
        source={"name": "CESM_REGRIDDED",
                "path": target_dir / "fake_bgc_surf_data.nc"},
    )
    rsf.save(target_dir / "example_input_bgc_surface_forcing.nc", group=False)


def create_roms_surface_forcing_restoring_sss(grid, target_dir: Path):
    rsf = rt.SurfaceForcing(
        grid=grid,
        type="restoring",
        restoring_forces=['sss'],
        start_time=dt.datetime(2010, 1, 1),
        end_time=dt.datetime(2010, 1, 2),
        source={"name": "WOA",
                "path": target_dir / "fake_restore_sss_surf_data.nc",
                "climatology": True},
    )
    rsf.save(target_dir / "example_input_restoring_sss.nc", group=False)


def create_roms_surface_forcing_restoring_dic_alk(grid, target_dir: Path):
    rsf = rt.SurfaceForcing(
        grid=grid,
        type="restoring",
        restoring_forces=['sDIC', 'sALK'],
        start_time=dt.datetime(2010, 1, 1),
        end_time=dt.datetime(2010, 1, 2),
        source={"name": "SODA",
                "path": target_dir / "fake_restore_dic_alk_surf_data.nc"},
    )
    rsf.save(target_dir / "example_input_restoring_dic_alk.nc", group=False)


def create_roms_co2_surface_forcing(grid, target_dir: Path):
    rsf = rt.SurfaceForcing(
        grid=grid,
        type="bgc",
        start_time=dt.datetime(2010, 1, 1),
        end_time=dt.datetime(2010, 1, 2),
        source={"name": "MBL_co2"}
    )
    rsf.save(target_dir / "example_input_co2_surface_forcing.nc", group=False)


def create_roms_initial_conditions(grid, target_dir: Path):
    ic = rt.InitialConditions(
        grid=grid,
        ini_time=dt.datetime(2010, 1, 1),
        source={"name": "GLORYS", "path": target_dir / "fake_phys_3d_data.nc"},
        bgc_source={"name": "CESM_REGRIDDED", "path": target_dir / "fake_bgc_3d_data.nc"},
    )
    # Add DOFE for BEC
    ic.ds["DOFE"] = ic.ds["DOP"]
    ic.ds["DOFE"][:] = 0.
    ic.ds["DOFE"].attrs = {
        "long_name": "Dissolved organic iron",
        "units": "mMol Fe m-3",
    }

    ic.save(target_dir / "example_input_bgc_initial_conditions.nc")


def create_roms_tides(grid, target_dir: Path):
    tpxo_dict = {
        "grid": target_dir / "fake_tides_data_g.nc",
        "h":    target_dir / "fake_tides_data_h.nc",
        "u":    target_dir / "fake_tides_data_u.nc",
    }
    tidal_forcing = rt.TidalForcing(
        grid=grid,
        source={"name": "TPXO", "path": tpxo_dict},
        ntides=2,
        model_reference_date=dt.datetime(2000, 1, 1),
        use_dask=True,
    )
    tidal_forcing.save(target_dir / "example_input_tides.nc")


def create_roms_cdr_forcing_3d(grid, irange, jrange, alk_pert, target_dir: Path):
    grid_ds = grid.ds

    cdr_time = xr.DataArray(
        np.array([3653.0, 3655.0]),
        dims=("cdr_time",),
        attrs={"long_name": "Time for CDR release", "units": "days"},
    )
    coords = {
        "cdr_time": cdr_time,
        "s_rho":    grid_ds.s_rho.astype("int32"),
        "xi_rho":   grid_ds.xi_rho.astype("int32"),
        "eta_rho":  grid_ds.eta_rho.astype("int32"),
    }
    shape_3d = (
        len(cdr_time),
        len(grid_ds.s_rho),
        len(grid_ds.eta_rho),
        len(grid_ds.xi_rho),
    )

    alk = xr.DataArray(
        np.zeros(shape_3d, dtype="f8"),
        dims=("cdr_time", "s_rho", "eta_rho", "xi_rho"),
        attrs={"long_name": "tracer flux [mmol/s]", "units": "mmol/s"},
    )
    alk[:, :, irange, jrange] = alk_pert

    dic = xr.DataArray(
        np.zeros(shape_3d, dtype="f8"),
        dims=("cdr_time", "s_rho", "eta_rho", "xi_rho"),
        attrs={"long_name": "tracer flux [mmol/s]", "units": "mmol/s"},
    )
    ds_3d = xr.Dataset(
        data_vars={"cdr_trcflx_3d_ALK": alk, "cdr_trcflx_3d_DIC": dic},
        coords=coords,
    )
    ds_3d.to_netcdf(target_dir / "cdr_forcing_3d.nc")


def create_roms_cdr_forcing_dp(grid, release_lon, release_lat, alk_pert, target_dir: Path):
    nt = 34
    ncdr = 1

    grid_ds = grid.ds
    cdr_time = xr.DataArray(
        np.array([3653.0, 3655.0]),
        dims=("cdr_time",),
        attrs={"long_name": "Time for CDR release", "units": "days"},
    )

    coords = {
        "s_rho": grid_ds.s_rho.astype("int32"),
        "cdr_time": cdr_time,
        "ncdr_prof": ("ncdr_prof", np.arange(1).astype("int32")),
        "nt": ("nt", np.arange(nt).astype("int32")),
        "one": ("one", np.arange(1).astype("int32")),
        "two": ("two", np.arange(2).astype("int32")),
    }

    depth_profiles = xr.DataArray(
        np.array([b"T"]),
        dims=("one",),
        attrs={"long_name": "depth profiles (T) or Gaussian (F)", "units": "nondim"},
    )
    profile = xr.DataArray(
        np.zeros((len(cdr_time), len(grid_ds.s_rho), 2, ncdr), dtype="f8"),
        dims=("cdr_time", "s_rho", "two", "ncdr_prof"),
        attrs={"long_name": "tracer flux [mmol/s]", "units": "mmol/s"},
    )
    profile[:, :, 0, :] = alk_pert

    cdr_lon = xr.DataArray(
        np.full((ncdr,), release_lon, dtype="f8"),
        dims=("ncdr_prof",),
        attrs={"long_name": "longitude of CDR release [degrees East]", "units": "deg E"},
    )
    cdr_lat = xr.DataArray(
        np.full((ncdr,), release_lat, dtype="f8"),
        dims=("ncdr_prof",),
        attrs={"long_name": "latitude of CDR release [degrees North]", "units": "deg N"},
    )
    cdr_layer_thickness = xr.DataArray(
        np.full((len(cdr_time), len(grid_ds.s_rho), ncdr), 1.0, dtype="f8"),
        dims=("cdr_time", "s_rho", "ncdr_prof"),
        attrs={
            "long_name": "layer thicknesses of CDR release given in a vertical profile [m]",
            "units": "m",
        },
    )

    ds_dp = xr.Dataset(
        {
            "depth_profiles": depth_profiles,
            "cdr_trcflx_profile": profile,
            "cdr_lon": cdr_lon,
            "cdr_lat": cdr_lat,
            "cdr_layer_thickness": cdr_layer_thickness,
        },
        coords=coords,
    )
    ds_dp["cdr_time"].attrs["long_name"] = "Time for CDR release"
    ds_dp["cdr_time"].attrs["units"] = "days"
    ds_dp.to_netcdf(target_dir / "cdr_forcing_dp.nc")


def create_roms_cdr_forcing_parm(grid, release_lon, release_lat, release_dep, alk_pert, target_dir: Path):
    ncdr = 1
    nt = 34
    cdr_time = xr.DataArray(
        np.array([3653.0, 3655.0]),
        dims=("cdr_time",),
        attrs={"long_name": "Time for CDR release", "units": "days"},
    )

    coords = {
        "cdr_time": cdr_time,
        "ncdr_parm": ("ncdr_parm", np.arange(1).astype("int32")),
        "nt": ("nt", np.arange(nt).astype("int32")),
    }

    trcflx_parm = xr.DataArray(
        np.zeros((len(cdr_time), nt, ncdr), dtype="f8"),
        dims=("cdr_time", "nt", "ncdr_parm"),
        attrs={"long_name": "tracer flux [mmol/s]", "units": "mmol/s"},
    )
    trcflx_parm[:, 11, :] = alk_pert

    cdr_lon = xr.DataArray(
        np.full((ncdr,), release_lon, dtype="f8"),
        dims=("ncdr_prof",),
        attrs={"long_name": "longitude of CDR release [degrees East]", "units": "deg E"},
    )
    cdr_lat = xr.DataArray(
        np.full((ncdr,), release_lat, dtype="f8"),
        dims=("ncdr_prof",),
        attrs={"long_name": "latitude of CDR release [degrees North]", "units": "deg N"},
    )
    cdr_dep = xr.DataArray(
        np.full((ncdr,), release_dep, dtype="f8"),
        dims=("ncdr_parm",),
        attrs={"long_name": "depth of CDR release [m]", "units": "m"},
    )
    cdr_hsc = xr.DataArray(
        np.full((ncdr,), 360.0, dtype="f8"),
        dims=("ncdr_parm",),
        attrs={"long_name": "horizontal scale of CDR release", "units": "m"},
    )
    cdr_vsc = xr.DataArray(
        np.full((ncdr,), 5.0, dtype="f8"),
        dims=("ncdr_parm",),
        attrs={"long_name": "vertical scale of CDR release", "units": "m"},
    )

    ds_parm = xr.Dataset(
        {
            "cdr_trcflx": trcflx_parm,
            "cdr_lon":    cdr_lon,
            "cdr_lat":    cdr_lat,
            "cdr_dep":    cdr_dep,
            "cdr_hsc":    cdr_hsc,
            "cdr_vsc":    cdr_vsc,
        },
        coords=coords,
    )

    ds_parm.to_netcdf(target_dir / "cdr_forcing_parm.nc")


def create_river_or_pipe_data(input_type, forcing_i, forcing_j, grid):
    grid_ds = grid.ds
    n_inputs = 1
    ntracers = 2
    forcing_time = np.arange(0.5, 360.5)

    ref_time = pd.Timestamp("2000-01-01 00:00:00")
    abs_time = ref_time + pd.to_timedelta(forcing_time, unit="D")

    forcing_time_coords = {
        f"{input_type}_time": forcing_time,
        "abs_time": (f"{input_type}_time", abs_time),
        f"n{input_type}": np.arange(1, n_inputs + 1),
        f"{input_type}_name": (
            f"n{input_type}",
            [f"Example{input_type.capitalize()}"],
        ),
    }
    tracer_coords = {
        "ntracers": np.arange(1, ntracers + 1),
        "tracer_name": ("ntracers", ["temp", "salt"]),
        "tracer_unit": ("ntracers", ["degrees Celsius", "PSU"]),
        "tracer_long_name": ("ntracers", ["potential temperature", "salinity"]),
    }

    constant_value = 2000.  # m3/s
    forcing_volume_data = np.full((len(forcing_time), n_inputs), constant_value, dtype=np.float32)
    forcing_volume = xr.DataArray(
        forcing_volume_data,
        dims=(f"{input_type}_time", f"n{input_type}"),
        coords=forcing_time_coords,
        attrs={
            "long_name": f"{input_type.capitalize()} volume flux",
            "units": "m^3/s",
        },
    )
    forcing_volume[f"{input_type}_time"].attrs.update({
        "long_name": "relative time: days since 2000-01-01 00:00:00",
        "units": "yearday",
        "cycle_length": 360,
    })

    forcing_tracer_data = np.zeros((len(forcing_time), ntracers, n_inputs), dtype=np.float32)
    forcing_tracer_data[:, 0, :] = 17.0  # Temperature
    forcing_tracer_data[:, 1, :] = 1.0   # Salinity
    forcing_tracer = xr.DataArray(
        forcing_tracer_data,
        dims=(f"{input_type}_time", "ntracers", f"n{input_type}"),
        coords={**forcing_time_coords, **tracer_coords},
        attrs={"long_name": f"{input_type.capitalize()} tracer data"},
    )
    forcing_tracer[f"{input_type}_time"].attrs.update(
        forcing_volume[f"{input_type}_time"].attrs
    )

    forcing_index_data = np.zeros((len(grid_ds.eta_rho), len(grid_ds.xi_rho)))
    forcing_index_data[forcing_i, forcing_j] = 1.
    forcing_index = xr.DataArray(
        forcing_index_data,
        dims=("eta_rho", "xi_rho"),
        attrs={"long_name": f"{input_type.capitalize()} ID", "units": "none"},
    )

    forcing_fraction_data = np.zeros((len(grid_ds.eta_rho), len(grid_ds.xi_rho)))
    forcing_fraction_data[forcing_i, forcing_j] = 1. / len(forcing_i)
    forcing_fraction = xr.DataArray(
        forcing_fraction_data,
        dims=("eta_rho", "xi_rho"),
        attrs={"long_name": f"{input_type.capitalize()} volume fraction", "units": "none"},
    )

    return xr.Dataset(
        data_vars={
            f"{input_type}_volume":   forcing_volume,
            f"{input_type}_tracer":   forcing_tracer,
            f"{input_type}_index":    forcing_index,
            f"{input_type}_fraction": forcing_fraction,
        },
    )


def create_flux_forcing_from_bulk(bulk_file: Path):
    """Derive a flux-form surface forcing file from the bulk forcing."""
    bulk_ds = xr.open_dataset(bulk_file, decode_times=False)

    def rho_to_u(A):
        return (0.5 * (A[..., :-1] + A[..., 1:])).rename({"xi_rho": "xi_u"})

    def rho_to_v(A):
        return (0.5 * (A[:, :-1, :] + A[:, 1:, :])).rename({"eta_rho": "eta_v"})

    rho_air = 1.22
    Cd = 1.3e-3

    U = bulk_ds["uwnd"]
    V = bulk_ds["vwnd"]
    W = np.hypot(U, V)

    tau_x = rho_air * Cd * W * U
    tau_y = rho_air * Cd * W * V

    sustr = rho_to_u(tau_x).astype("float32")
    svstr = rho_to_v(tau_y).astype("float32")

    shflux = (bulk_ds["swrad"] + bulk_ds["lwrad"]).astype("float32")
    swrad  = bulk_ds["swrad"].astype("float32")

    # freshwater flux: E - P (placeholder)
    swflux = (-bulk_ds["rain"]).astype("float32")

    template = bulk_ds["Tair"]
    SST    = xr.full_like(template, 20.0).astype("float32")
    SSS    = xr.full_like(template, 35.0).astype("float32")
    dQdSST = xr.full_like(template, -30.0).astype("float32")

    time = bulk_ds["time"].astype("float32")

    forcing = xr.Dataset(
        data_vars=dict(
            sustr  = (("sms_time", "eta_rho", "xi_u"),   sustr.data),
            svstr  = (("sms_time", "eta_v",   "xi_rho"), svstr.data),
            shflux = (("shf_time", "eta_rho", "xi_rho"), shflux.data),
            swflux = (("swf_time", "eta_rho", "xi_rho"), swflux.data),
            SST    = (("sst_time", "eta_rho", "xi_rho"), SST.data),
            SSS    = (("sss_time", "eta_rho", "xi_rho"), SSS.data),
            dQdSST = (("sst_time", "eta_rho", "xi_rho"), dQdSST.data),
            swrad  = (("srf_time", "eta_rho", "xi_rho"), swrad.data),
        ),
        coords=dict(
            sms_time=("sms_time", time.data),
            shf_time=("shf_time", time.data),
            swf_time=("swf_time", time.data),
            sst_time=("sst_time", time.data),
            sss_time=("sss_time", time.data),
            srf_time=("srf_time", time.data),
        ),
    )

    forcing.to_netcdf(bulk_file.parent / "example_input_surface_flux_forcing.nc")
