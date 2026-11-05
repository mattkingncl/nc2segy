#!/usr/bin/env python3
"""
Convert airborne ice-penetrating-radar profiles from NetCDF to SEG-Y rev 1.

Each radar trace (one position along the flight line) becomes one SEG-Y trace,
and the fast-time axis (two-way travel time) becomes the SEG-Y sample axis.
Samples are written as 32-bit IEEE floats (format code 5).

Two file layouts are recognised automatically:
  bas   BAS / UK Polar Data Centre PASIN files (e.g. WISE_ISODYN_W15.nc)
  utig  UTIG / ICECAP HiCARS L1B files (e.g. ER2HI1B_*.nc)
Anything else can be converted by naming the variables on the command line.

Time convention
---------------
SEG-Y stores the sample interval as a whole number of microseconds, which is
too coarse for radar. As in the BAS SEG-Y files, the default here stores the
interval in NANOSECONDS, so seismic software that labels its time axis "ms"
is really showing microseconds of two-way travel time. If the interval is not
a whole number of nanoseconds (PASIN samples at 22 MHz, i.e. 45.45 ns), the
traces are by default linearly resampled onto the nearest whole-nanosecond
interval, so travel times stay correct. "--inexact round" instead keeps the
samples and rounds the header value, which is what the BAS files do; that
stretches the time axis by the rounding error (about 1 % for 45.45 -> 45 ns).

Examples
--------
    python nc2segy.py WISE_ISODYN_W15.nc --outdir segy/
    python nc2segy.py ER2HI1B_2016009_ASB1_JKB2n_R44a_000.nc --outdir segy/
    python nc2segy.py W15.nc --data chirp_data --date 2005-12-01   # BAS has no flight date

Requires: numpy, pandas, xarray, netCDF4, segyio; pyproj when projecting lat/lon.
"""

import argparse
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import segyio
import xarray as xr

# ---------------------------------------------------------------- constants --

# Seconds in one unit of each fast-time unit string we expect to meet.
TIME_UNITS_TO_S = {
    "s": 1.0, "seconds": 1.0, "ms": 1e-3, "milliseconds": 1e-3,
    "us": 1e-6, "µs": 1e-6, "microseconds": 1e-6, "microsecond": 1e-6,
    "ns": 1e-9, "nanoseconds": 1e-9, "ps": 1e-12, "picoseconds": 1e-12,
}

# What one stored header unit means for each choice of --time-units.
HEADER_UNIT_S = {"us": 1e-6, "ns": 1e-9, "ps": 1e-12}

MAX_INT16 = 32767          # SEG-Y rev 1 sample interval and count are 16-bit
MISSING = -9999.0          # BAS marks missing navigation values this way

# Variable names used by the two known layouts. "time_of_day" is seconds after
# UTC midnight (BAS has no flight date in the file); "time" is a full timestamp.
PRESETS = {
    "bas": dict(
        detect="PriNumber_layerData",
        data=["chirp_data", "pulse_data"], trace_dim="traces",
        sample_dim="fast_time", sample_var="fast_time",
        x="x_coordinates", y="y_coordinates", crs="EPSG:3031",
        elev="aircraft_altitude_layerData", surface="surface_altitude_layerData",
        time_of_day="UTC_time_layerData", record="PriNumber_layerData",
        trace_no="traces",
        # The fast_time axis in these files spans 0 to N*dt over N samples, so its
        # spacing is slightly too large; the 22 MHz digitiser (45.4545 ns) is the
        # interval BAS's own ice thicknesses are consistent with.
        prefer_metadata_dt=True),
    "utig": dict(
        detect="amplitude_high_gain",
        data=["amplitude_high_gain", "amplitude_low_gain"], trace_dim="time",
        sample_dim="fasttime_dim", sample_var="fasttime",
        lat="lat", lon="lon", crs="EPSG:3031",
        elev="altitude", time="time"),
}

# ------------------------------------------------------------------ helpers --


def detect_preset(ds):
    """Return the name of the matching layout, or None."""
    for name, p in PRESETS.items():
        if p["detect"] in ds.variables:
            return name
    return None


def clean(values):
    """Return a float copy with BAS -9999 flags and non-finite values as NaN."""
    v = np.asarray(values, dtype=float).copy()
    v[(v == MISSING) | ~np.isfinite(v)] = np.nan
    return v


def fast_time_seconds(ds, sample_var, sample_dim):
    """Fast-time axis in seconds, taken from the named 1-D variable."""
    var = ds[sample_var]
    if var.dims != (sample_dim,):
        raise ValueError(f"'{sample_var}' should lie along '{sample_dim}', not {var.dims}.")
    units = var.attrs.get("units", "").strip().lower()
    if units not in TIME_UNITS_TO_S:
        raise ValueError(f"Unknown fast-time units '{units}'; use --dt-ns.")
    return var.values.astype(float) * TIME_UNITS_TO_S[units]


def metadata_interval_s(ds):
    """
    Look for a sampling interval or rate quoted in the global attributes,
    so it can be checked against the fast-time axis. Returns seconds or None.
    """
    text = " ".join(str(v) for v in ds.attrs.values())
    m = re.search(r"sampling interval of\s*([\d.]+)\s*ns", text, re.I)
    if m:
        return float(m.group(1)) * 1e-9
    m = re.search(r"sampl\w* (?:rate|frequency)[^\d]{0,20}([\d.]+)\s*MHz", text, re.I)
    if m:
        return 1.0 / (float(m.group(1)) * 1e6)
    return None


def regular_axis(t):
    """
    Check the fast-time axis is evenly spaced and return (t0, dt) in seconds.
    Small jitter from float32 storage is tolerated; real gaps are not.
    """
    n = t.size
    dt = (t[-1] - t[0]) / (n - 1)
    resid = np.abs(t - (t[0] + dt * np.arange(n)))
    if resid.max() > 0.01 * dt:
        raise ValueError(f"Fast time is unevenly spaced (worst sample is "
                         f"{resid.max() / dt:.1%} of an interval off a regular grid).")
    return float(t[0]), float(dt)


def resample_weights(t_old, t_new):
    """Indices and weights for linear interpolation from t_old onto t_new."""
    i = np.clip(np.searchsorted(t_old, t_new, side="right") - 1, 0, t_old.size - 2)
    w = (t_new - t_old[i]) / (t_old[i + 1] - t_old[i])
    return i, np.clip(w, 0.0, 1.0)


def text_header(lines):
    """Build the 40-line, 80-column textual header (segyio writes it as EBCDIC)."""
    lines = [l.upper()[:76] for l in lines][:40]
    lines += [""] * (40 - len(lines))
    return "".join(f"C{i + 1:2d} {l:<76}" for i, l in enumerate(lines))


# --------------------------------------------------------------- converter --


def convert(ds, src, data_name, out_path, cfg, args):
    """Write one radar variable of an opened dataset to one SEG-Y file."""
    tdim, sdim = cfg["trace_dim"], cfg["sample_dim"]
    da = ds[data_name]
    if set(da.dims) != {tdim, sdim}:
        raise ValueError(f"'{data_name}' has dims {da.dims}; expected {tdim} and {sdim}.")
    da = da.transpose(tdim, sdim)
    n_tr = da.sizes[tdim]

    # ---- fast-time axis: what the file says, and what we will store --------
    t_file = fast_time_seconds(ds, cfg["sample_var"], sdim)
    t0, dt = regular_axis(t_file)
    if args.dt_ns:                                   # user override
        dt = args.dt_ns * 1e-9
    t_old = t0 + dt * np.arange(t_file.size)         # regularised true times

    meta_dt = metadata_interval_s(ds)
    if meta_dt and not args.dt_ns and cfg.get("prefer_metadata_dt"):
        if abs(meta_dt - dt) / meta_dt > 2e-4:
            print(f"    using metadata interval {meta_dt * 1e9:.4f} ns instead of the "
                  f"axis spacing {dt * 1e9:.4f} ns")
        dt = meta_dt
        t_old = t0 + dt * np.arange(t_file.size)
    elif meta_dt and not args.dt_ns and abs(meta_dt - dt) / meta_dt > 2e-4:
        warnings.warn(
            f"Fast-time axis spacing is {dt * 1e9:.4f} ns, but the metadata quote "
            f"{meta_dt * 1e9:.4f} ns ({(dt / meta_dt - 1):+.3%}). If the metadata are "
            f"right, rerun with --dt-ns {meta_dt * 1e9:.4f}.")

    unit = HEADER_UNIT_S[args.time_units]
    exact = abs(dt / unit - round(dt / unit)) < 1e-4 * dt / unit
    resample = (not exact) and args.inexact == "resample"
    if not exact and args.inexact == "fail":
        raise ValueError(f"Interval {dt * 1e9:.4f} ns is not a whole number of {args.time_units}.")
    dt_stored = max(1, int(round(dt / unit)))
    if resample:
        dt_new = dt_stored * unit
        n_new = int(np.floor((t_old[-1] - t0) / dt_new + 1e-9)) + 1
        t_new = t0 + dt_new * np.arange(n_new)
        idx, wgt = resample_weights(t_old, t_new)
        n_smp = n_new
        note = (f"RESAMPLED LINEARLY FROM {dt * 1e9:.4f} NS TO {dt_stored} {args.time_units}")
    else:
        n_smp = t_file.size
        stretch = dt_stored * unit / dt - 1
        if not exact:
            warnings.warn(f"Header interval rounded from {dt * 1e9:.4f} ns to "
                          f"{dt_stored} {args.time_units}: every travel time read from "
                          f"the SEG-Y is off by {stretch:+.2%} "
                          f"({stretch * t_old[-1] * 1e6:+.3f} us at the last sample).")
        note = f"NATIVE SAMPLES; HEADER INTERVAL ERROR {stretch:+.3%}"
    if dt_stored > MAX_INT16 or n_smp > MAX_INT16:
        raise ValueError("Interval or sample count does not fit SEG-Y rev 1 16-bit fields; "
                         "choose coarser --time-units or crop fast time.")

    # ---- per-trace navigation ----------------------------------------------
    def var(key):
        name = cfg.get(key)
        if not name or name not in ds.variables:
            return None
        v = ds[name]
        if v.dims != (tdim,):
            raise ValueError(f"'{name}' should lie along '{tdim}', not {v.dims}.")
        return v.values

    crs = cfg.get("crs")
    if var("x") is not None and var("y") is not None:
        xs, ys = clean(var("x")), clean(var("y"))
    elif var("lon") is not None and var("lat") is not None:
        from pyproj import Transformer   # only needed for lat/lon input
        lon, lat = clean(var("lon")), clean(var("lat"))
        xs, ys = Transformer.from_crs("EPSG:4326", crs, always_xy=True).transform(lon, lat)
    else:
        xs = ys = None
    elev, surf = (clean(var(k)) if var(k) is not None else None for k in ("elev", "surface"))
    record = var("record")
    trace_no = var("trace_no")

    # Acquisition time: full timestamps, or seconds of day plus --date.
    stamps = None
    if var("time") is not None and np.issubdtype(var("time").dtype, np.datetime64):
        stamps = pd.to_datetime(var("time"))
    elif var("time_of_day") is not None:
        sod = clean(var("time_of_day"))
        sod = sod + 86400 * np.cumsum(np.r_[0, np.diff(sod) < -43200])  # past midnight
        day = pd.Timestamp(args.date) if args.date else pd.Timestamp("1970-01-01")
        stamps = day + pd.to_timedelta(sod, unit="s")
        if not args.date:
            warnings.warn("File gives time of day but no date; year and day-of-year "
                          "are left at zero. Use --date YYYY-MM-DD to fill them.")

    cs, es = args.coord_scalar, -100            # elevations always to 1 cm
    def scaled(v, s):
        # SEG-Y scalar convention: negative s divides, positive s multiplies.
        return int(round(v * (-s if s < 0 else 1 / s)))

    # ---- textual and binary headers ----------------------------------------
    units = da.attrs.get("units", "unknown")
    lines = [
        f"RADAR PROFILE CONVERTED FROM NETCDF: {Path(src).name}",
        f"VARIABLE {data_name} UNITS {units}",
        f"INSTITUTION {str(ds.attrs.get('institution', '')).splitlines()[0] if ds.attrs.get('institution') else 'UNKNOWN'}",
        f"FLIGHT {ds.attrs.get('flight', '')} CAMPAIGN {ds.attrs.get('campaign', '')}",
        f"SAMPLE INTERVAL HEADER UNIT = 1 {args.time_units.upper()} (NOT US)",
        "SO A TIME AXIS LABELLED MS SHOWS US OF TWO-WAY TRAVEL TIME" if args.time_units == "ns" else "",
        note,
        f"FIRST SAMPLE {t0 * 1e6:.4f} US OF RECORD TIME",
        f"COORDS {crs} METRES, SCALAR {cs}; ELEVATIONS M WGS84 SCALAR {es}",
        "BYTE 41 ANTENNA ALTITUDE; BYTE 45 ICE SURFACE ELEVATION IF GIVEN",
        "BYTE 9 FIELD RECORD = PRI NUMBER IF GIVEN; BYTE 17 = TRACE NUMBER",
        "TIMES UTC" + ("" if stamps is None or args.date or var("time") is not None
                      else "; DATE UNKNOWN"),
        "SEG-Y REV 1, IEEE FLOAT32, BIG-ENDIAN",
    ]
    if "TX-record_offset" in ds.attrs:
        lines.append(str(ds.attrs["TX-record_offset"]).strip())

    spec = segyio.spec()
    spec.format, spec.sorting, spec.tracecount = 5, 0, n_tr
    spec.samples = np.arange(n_smp) * dt_stored

    n_bad = n_nopos = 0
    with segyio.create(out_path, spec) as f:
        f.text[0] = text_header(lines)
        f.bin.update({
            segyio.BinField.Interval: dt_stored, segyio.BinField.Samples: n_smp,
            segyio.BinField.Format: 5, segyio.BinField.MeasurementSystem: 1,
            segyio.BinField.SEGYRevision: 1, segyio.BinField.SEGYRevisionMinor: 0,
            segyio.BinField.TraceFlag: 1,
        })
        for start in range(0, n_tr, args.block):
            stop = min(start + args.block, n_tr)
            block = da.isel({tdim: slice(start, stop)}).values.astype(np.float64)
            if start == 0 and "db" in units.lower() and not args.to_db:
                lo, hi = np.nanmin(block), np.nanmax(block)
                if lo >= 0 and hi > 1e4 * max(lo, 1e-30):
                    warnings.warn(f"Units say '{units}', but values span {lo:.3g} to "
                                  f"{hi:.3g}, which looks like linear power, not dB. "
                                  "Use --to-db to write 10*log10 values.")
            if args.to_db:
                block = 10 * np.log10(np.maximum(block, np.finfo(np.float32).tiny))
            if resample:
                block = block[:, idx] * (1 - wgt) + block[:, idx + 1] * wgt
            block = block.astype(np.float32)
            bad = ~np.isfinite(block)
            n_bad += int(bad.sum())
            block[bad] = 0.0
            for j, i in enumerate(range(start, stop)):
                tn = int(trace_no[i]) if trace_no is not None else i + 1
                h = {
                    segyio.TraceField.TRACE_SEQUENCE_LINE: i + 1,
                    segyio.TraceField.TRACE_SEQUENCE_FILE: i + 1,
                    segyio.TraceField.FieldRecord: int(record[i]) if record is not None else i + 1,
                    segyio.TraceField.TraceNumber: 1,
                    segyio.TraceField.EnergySourcePoint: tn,
                    segyio.TraceField.CDP: tn,
                    segyio.TraceField.TraceIdentificationCode: 1,
                    segyio.TraceField.NStackedTraces: 1,
                    segyio.TraceField.TRACE_SAMPLE_COUNT: n_smp,
                    segyio.TraceField.TRACE_SAMPLE_INTERVAL: dt_stored,
                    segyio.TraceField.SourceGroupScalar: cs,
                    segyio.TraceField.ElevationScalar: es,
                    segyio.TraceField.CoordinateUnits: 1,
                }
                if xs is not None and np.isfinite(xs[i]) and np.isfinite(ys[i]):
                    x, y = scaled(xs[i], cs), scaled(ys[i], cs)
                    for fx, fy in ((segyio.TraceField.SourceX, segyio.TraceField.SourceY),
                                   (segyio.TraceField.GroupX, segyio.TraceField.GroupY),
                                   (segyio.TraceField.CDP_X, segyio.TraceField.CDP_Y)):
                        h[fx], h[fy] = x, y
                else:
                    n_nopos += 1
                if elev is not None and np.isfinite(elev[i]):
                    h[segyio.TraceField.ReceiverGroupElevation] = scaled(elev[i], es)
                if surf is not None and np.isfinite(surf[i]):
                    h[segyio.TraceField.SourceSurfaceElevation] = scaled(surf[i], es)
                if stamps is not None and not pd.isna(stamps[i]):
                    s = stamps[i]
                    if args.date or var("time") is not None:
                        h[segyio.TraceField.YearDataRecorded] = s.year
                        h[segyio.TraceField.DayOfYear] = s.dayofyear
                    h[segyio.TraceField.HourOfDay] = s.hour
                    h[segyio.TraceField.MinuteOfHour] = s.minute
                    h[segyio.TraceField.SecondOfMinute] = s.second
                    h[segyio.TraceField.TimeBaseCode] = 4          # UTC
                f.header[i] = h
                f.trace[i] = np.ascontiguousarray(block[j])

    print(f"  {data_name}: {n_tr} traces x {n_smp} samples, interval {dt_stored} "
          f"{args.time_units} -> {out_path}")
    if n_bad:
        print(f"    {n_bad} non-finite samples written as 0")
    if n_nopos:
        print(f"    {n_nopos} traces without position (coordinates left at 0)")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("inputs", nargs="+", help="NetCDF file(s)")
    p.add_argument("--outdir", default=".", help="output folder (default: current)")
    p.add_argument("--preset", choices=list(PRESETS), help="file layout (default: detect)")
    p.add_argument("--data", action="append", help="radar variable(s); default: all in preset")
    p.add_argument("--trace-dim"), p.add_argument("--sample-dim")
    p.add_argument("--sample-var", help="1-D fast-time variable with time units")
    p.add_argument("--x"), p.add_argument("--y")
    p.add_argument("--lat"), p.add_argument("--lon")
    p.add_argument("--crs", help="projected CRS for coordinates (default EPSG:3031)")
    p.add_argument("--elev", help="antenna altitude variable (m)")
    p.add_argument("--surface", help="ice-surface elevation variable (m)")
    p.add_argument("--time", help="per-trace timestamp variable")
    p.add_argument("--date", help="flight date YYYY-MM-DD when the file gives only time of day")
    p.add_argument("--dt-ns", type=float, help="override the sample interval (ns)")
    p.add_argument("--time-units", choices=list(HEADER_UNIT_S), default="ns",
                   help="meaning of one header time unit (default ns, as BAS)")
    p.add_argument("--inexact", choices=["resample", "round", "fail"], default="resample",
                   help="what to do if the interval is not a whole header unit")
    p.add_argument("--to-db", action="store_true", help="write 10*log10 of linear power")
    p.add_argument("--coord-scalar", type=int, default=-100,
                   help="SEG-Y coordinate scalar (default -100 = cm; 1 = whole metres)")
    p.add_argument("--block", type=int, default=2000, help="traces read at a time")
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    for src in args.inputs:
        ds = xr.open_dataset(src)
        name = args.preset or detect_preset(ds)
        cfg = dict(PRESETS[name]) if name else {"crs": "EPSG:3031"}
        for key in ("trace_dim", "sample_dim", "sample_var", "x", "y", "lat", "lon",
                    "crs", "elev", "surface", "time"):
            if getattr(args, key):
                cfg[key] = getattr(args, key)
        data = args.data or cfg.get("data")
        if not data or not cfg.get("trace_dim") or not cfg.get("sample_dim"):
            raise SystemExit(f"{src}: layout not recognised; give --data, --trace-dim, "
                             "--sample-dim and --sample-var.")
        cfg.setdefault("sample_var", cfg["sample_dim"])
        print(f"{src} (layout: {name or 'custom'})")
        for d in data:
            convert(ds, src, d, str(outdir / f"{Path(src).stem}_{d}.sgy"), cfg, args)
        ds.close()


if __name__ == "__main__":
    main()
