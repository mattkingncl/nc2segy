#!/usr/bin/env python3
"""
Convert airborne ice-penetrating-radar profiles to SEG-Y rev 1, from NetCDF or
from existing radar SEG-Y files whose time axis needs correcting.

Each radar trace (one position along the flight line) becomes one SEG-Y trace,
and the fast-time axis (two-way travel time) becomes the SEG-Y sample axis.
Samples are written as 32-bit IEEE floats (format code 5).

Recognised inputs:
  bas   BAS / UK Polar Data Centre PASIN NetCDF (e.g. WISE_ISODYN_W15.nc)
  utig  UTIG / ICECAP HiCARS L1B NetCDF (e.g. ER2HI1B_*.nc)
  .segy / .sgy files, including the BAS-supplied PASIN SEG-Y (e.g. W15_chirp.segy)
Other NetCDF layouts can be converted by naming the variables.

Time axis (--axis)
------------------
SEG-Y stores the sample interval as a whole number of microseconds, which is
too coarse for radar, so the header unit is redefined. --axis says what a time
axis labelled "ms" in interpretation software will show:
  us  (default) header unit = 1 ns, so "ms" reads as microseconds, as in the
      BAS SEG-Y files. A 40 us bed echo appears at 40 "ms".
  ns  header unit = 1 ps, so "ms" reads as nanoseconds (40 us -> 40000 "ms").
      Intervals above 32.767 ns do not fit the 16-bit field, so such traces
      are upsampled by a whole factor; every original sample is kept.
  ms  true SEG-Y units. Refused for radar, because whole microseconds would
      discard most of the samples.
If the interval is not a whole number of header units (PASIN: 45.45 ns), the
traces are by default linearly resampled so travel times stay correct;
"--inexact round" keeps the samples and rounds the header value instead,
which is what the BAS files do and makes every travel time about 1 % short.

Examples
--------
    python nc2segy.py WISE_ISODYN_W15.nc --outdir segy/
    python nc2segy.py ER2HI1B_2016009_ASB1_JKB2n_R44a_000.nc --axis ns --outdir segy/
    python nc2segy.py W15_chirp.segy --outdir segy/          # repair BAS SEG-Y
    python nc2segy.py WISE_ISODYN_W15.nc W15_chirp.segy --kingdom --outdir kingdom/
    python nc2segy.py W15.nc --data chirp_data --date YYYY-MM-DD   # BAS files give no flight date

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

# Show warnings as one plain line instead of Python's file-and-line format.
warnings.formatwarning = lambda message, *args, **kwargs: f"  WARNING: {message}\n"

# ---------------------------------------------------------------- constants --

# Seconds in one unit of each fast-time unit string we expect to meet.
TIME_UNITS_TO_S = {
    "s": 1.0, "seconds": 1.0, "ms": 1e-3, "milliseconds": 1e-3,
    "us": 1e-6, "µs": 1e-6, "microseconds": 1e-6, "microsecond": 1e-6,
    "ns": 1e-9, "nanoseconds": 1e-9, "ps": 1e-12, "picoseconds": 1e-12,
}

# What a time axis labelled "ms" in interpretation software should show, and
# the length of one stored header unit that achieves it. Software reads the
# header interval as microseconds, so storing nanoseconds makes "ms" read as us.
AXIS_UNIT_S = {"ms": 1e-6, "us": 1e-9, "ns": 1e-12}
AXIS_UNIT_NAME = {"ms": "us", "us": "ns", "ns": "ps"}

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
        linear_power=True,           # values are linear power despite "dBm" label
        trace_no="traces",
        # The fast_time axis in these files spans 0 to N*dt over N samples, so its
        # spacing is slightly too large. Use the 22 MHz digitiser interval
        # (45.4545 ns; the metadata round it to 45.4546), which BAS's own ice
        # thicknesses are consistent with. The same value is used for BAS SEG-Y.
        dt_s=1 / 22e6),
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


def plan_time_axis(t0, dt, n, args):
    """
    Decide how the fast-time axis will be stored, given the true first-sample
    time t0 and interval dt (seconds) and the sample count n.

    SEG-Y holds the interval as a signed 16-bit whole number of header units
    (1 to 32767). Three things can force the samples to be recomputed:
      * the interval is too long for 32767 units (e.g. 45 ns in picoseconds),
        so the trace is upsampled by a whole factor k; every original sample
        is kept, and linear interpolation fills the points in between;
      * the interval is not a whole number of units (45.45 ns in nanoseconds),
        so with --inexact resample the trace is interpolated onto a grid that is;
      * the units are so coarse that whole units would discard samples, which
        is refused.
    Returns a dict with the stored interval, sample count, interpolation
    indices and weights (None if the samples are written unchanged) and a note.
    """
    unit = AXIS_UNIT_S[args.axis]
    uname = AXIS_UNIT_NAME[args.axis]

    # Upsampling factor needed for the interval to fit in 16 bits.
    k = int(np.ceil(dt / unit / MAX_INT16))
    dt_k = dt / k
    stored = max(1, int(round(dt_k / unit)))
    if stored * unit > 1.05 * dt_k:
        raise ValueError(
            f"A {dt * 1e9:.4f} ns interval cannot be written in whole {uname}: it would "
            f"need resampling to {stored} {uname}, discarding about "
            f"{1 - dt_k / (stored * unit):.0%} of the samples. Choose a finer --axis.")
    exact = abs(dt_k / unit - stored) < 1e-4 * dt_k / unit
    if not exact and args.inexact == "fail":
        raise ValueError(f"Interval {dt * 1e9:.4f} ns is not a whole number of {uname}.")

    # Interval of the grid actually written.
    grid_dt = stored * unit if (not exact and args.inexact == "resample") else dt_k
    t_old = t0 + dt * np.arange(n)
    if grid_dt == dt:
        n_new, idx, wgt = n, None, None
    else:
        n_new = int(np.floor((t_old[-1] - t0) / grid_dt + 1e-9)) + 1
        idx, wgt = resample_weights(t_old, t0 + grid_dt * np.arange(n_new))

    # Relative error of the stored header interval against the written grid.
    stretch = stored * unit / grid_dt - 1
    if abs(stretch) > 1e-4:
        warnings.warn(f"Header interval rounded from {grid_dt / unit:.4f} to {stored} {uname}: "
                      f"every travel time read from the SEG-Y is off by {stretch:+.2%} "
                      f"({stretch * (n_new - 1) * grid_dt * 1e6:+.3f} us at the last sample).")
    if n_new > MAX_INT16:
        raise ValueError(f"{n_new} samples per trace exceeds the SEG-Y rev 1 limit of "
                         f"{MAX_INT16}; choose a coarser --axis or crop fast time.")

    if idx is None:
        note = "ORIGINAL SAMPLES, UNCHANGED"
    elif grid_dt == dt_k:
        note = f"UPSAMPLED X{k} BY LINEAR INTERPOLATION; ORIGINAL SAMPLES KEPT"
    else:
        note = f"RESAMPLED LINEARLY FROM {dt * 1e9:.4f} NS TO {grid_dt * 1e9:.4f} NS"
    err = (f"HEADER INTERVAL ERROR {stretch:+.4%} "
           f"({stretch * (n_new - 1) * grid_dt * 1e9:+.2f} NS AT LAST SAMPLE)")
    return dict(stored=stored, n=n_new, idx=idx, wgt=wgt, note=note, err=err, uname=uname)


def axis_text_lines(args, plan):
    """Text-header lines explaining the time convention."""
    shows = {"ms": "MS", "us": "US", "ns": "NS"}[args.axis]
    return [f"SAMPLE INTERVAL HEADER UNIT = 1 {plan['uname'].upper()}",
            f"SO A TIME AXIS LABELLED MS SHOWS {shows} OF TWO-WAY TRAVEL TIME",
            plan["note"], plan["err"]]


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

    if cfg.get("dt_s") and not args.dt_ns:          # known digitiser interval
        if abs(cfg["dt_s"] - dt) / dt > 2e-4:
            print(f"    using digitiser interval {cfg['dt_s'] * 1e9:.4f} ns instead of "
                  f"the axis spacing {dt * 1e9:.4f} ns")
        dt = cfg["dt_s"]
    meta_dt = metadata_interval_s(ds)
    if meta_dt and not args.dt_ns and abs(meta_dt - dt) / meta_dt > 2e-4:
        warnings.warn(
            f"Fast-time axis spacing is {dt * 1e9:.4f} ns, but the metadata quote "
            f"{meta_dt * 1e9:.4f} ns ({(dt / meta_dt - 1):+.3%}). If the metadata are "
            f"right, rerun with --dt-ns {meta_dt * 1e9:.4f}.")

    plan = plan_time_axis(t0, dt, t_file.size, args)
    dt_stored, n_smp, note = plan["stored"], plan["n"], plan["note"]

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
    # --kingdom converts known linear-power data (BAS) to dB for display.
    to_db = args.to_db or (args.kingdom and cfg.get("linear_power", False))
    def scaled(v, s):
        # SEG-Y scalar convention: negative s divides, positive s multiplies.
        return int(round(v * (-s if s < 0 else 1 / s)))

    # ---- textual and binary headers ----------------------------------------
    units = da.attrs.get("units", "unknown")
    lines = [
        f"RADAR PROFILE CONVERTED FROM NETCDF: {Path(src).name}",
        f"VARIABLE {data_name} UNITS {units}",
        "DATA WRITTEN AS 10*LOG10 OF ORIGINAL VALUES (DB)" if to_db else "DATA WRITTEN AS ORIGINAL VALUES",
        f"INSTITUTION {str(ds.attrs.get('institution', '')).splitlines()[0] if ds.attrs.get('institution') else 'UNKNOWN'}",
        f"FLIGHT {ds.attrs.get('flight', '')} CAMPAIGN {ds.attrs.get('campaign', '')}",
        *axis_text_lines(args, plan),
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
            if start == 0 and "db" in units.lower() and not to_db:
                lo, hi = np.nanmin(block), np.nanmax(block)
                if lo >= 0 and hi > 1e4 * max(lo, 1e-30):
                    warnings.warn(f"Units say '{units}', but values span {lo:.3g} to "
                                  f"{hi:.3g}, which looks like linear power, not dB. "
                                  "Use --to-db to write 10*log10 values.")
            if to_db:
                block = 10 * np.log10(np.maximum(block, np.finfo(np.float32).tiny))
            if plan["idx"] is not None:
                i0, w = plan["idx"], plan["wgt"]
                block = block[:, i0] * (1 - w) + block[:, i0 + 1] * w
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
          f"{plan['uname']} -> {out_path}")
    if n_bad:
        print(f"    {n_bad} non-finite samples written as 0")
    if n_nopos:
        print(f"    {n_nopos} traces without position (coordinates left at 0)")


def convert_segy(src, out_path, args):
    """
    Rewrite an existing radar SEG-Y file with a correct time axis.

    Made for the BAS-supplied PASIN SEG-Y files, which store a rounded 45 ns
    interval in the trace headers, 1000 in the binary header, and file-creation
    dates in the acquisition-time fields. For BAS files the true interval is
    taken as 45.4546 ns (the 22 MHz PASIN digitiser); other SEG-Y files need
    --dt-ns. All other trace-header values are copied unchanged, except that
    empty CDP number and CDP x/y fields are filled from the source fields.
    """
    T = segyio.TraceField
    with segyio.open(src, ignore_geometry=True) as f:
        n_tr, n = f.tracecount, f.header[0][T.TRACE_SAMPLE_COUNT]
        trace_dt, bin_dt = f.header[0][T.TRACE_SAMPLE_INTERVAL], f.bin[segyio.BinField.Interval]
        is_bas = trace_dt == 45 and bin_dt == 1000
        if args.dt_ns:
            dt = args.dt_ns * 1e-9
        elif is_bas:
            dt = PRESETS["bas"]["dt_s"]
            print(f"  BAS PASIN SEG-Y recognised (trace interval 45, binary 1000); "
                  f"using {dt * 1e9:.4f} ns (22 MHz)")
        else:
            raise ValueError("true sample interval unknown; give --dt-ns.")
        t0 = 0.0                                       # delays not used by BAS
        to_db = args.to_db or (args.kingdom and is_bas)
        plan = plan_time_axis(t0, dt, n, args)

        lines = [
            f"RADAR SEG-Y REWRITTEN FROM {Path(src).name}",
            f"TRUE SAMPLE INTERVAL {dt * 1e9:.4f} NS; ORIGINAL HEADERS: TRACE "
            f"{trace_dt}, BINARY {bin_dt}",
            *axis_text_lines(args, plan),
            "DATA " + ("10*LOG10 OF ORIGINAL VALUES (DB)" if to_db else "ORIGINAL VALUES"),
            "OTHER TRACE HEADERS COPIED FROM ORIGINAL; CDP FILLED FROM SOURCE",
        ]
        if is_bas:
            lines += ["BAS VALUES ARE LINEAR POWER DESPITE DBM LABEL IN NETCDF",
                      "ORIGINAL YEAR/DAY/TIME WERE FILE-CREATION VALUES; SET TO ZERO",
                      "COORDS EPSG:3031 METRES; BYTE 9 = PRI NUMBER"]
        lines.append("SEG-Y REV 1, IEEE FLOAT32, BIG-ENDIAN")

        spec = segyio.spec()
        spec.format, spec.sorting, spec.tracecount = 5, 0, n_tr
        spec.samples = np.arange(plan["n"]) * plan["stored"]
        with segyio.create(out_path, spec) as g:
            g.text[0] = text_header(lines)
            g.bin.update({
                segyio.BinField.Interval: plan["stored"], segyio.BinField.Samples: plan["n"],
                segyio.BinField.Format: 5, segyio.BinField.MeasurementSystem: 1,
                segyio.BinField.SEGYRevision: 1, segyio.BinField.SEGYRevisionMinor: 0,
                segyio.BinField.TraceFlag: 1,
            })
            for start in range(0, n_tr, args.block):
                stop = min(start + args.block, n_tr)
                block = np.asarray(f.trace.raw[start:stop], dtype=np.float64)
                if to_db:
                    block = 10 * np.log10(np.maximum(block, np.finfo(np.float32).tiny))
                if plan["idx"] is not None:
                    i0, w = plan["idx"], plan["wgt"]
                    block = block[:, i0] * (1 - w) + block[:, i0 + 1] * w
                block = block.astype(np.float32)
                for j, i in enumerate(range(start, stop)):
                    h = dict(f.header[i])
                    h[T.TRACE_SAMPLE_COUNT] = plan["n"]
                    h[T.TRACE_SAMPLE_INTERVAL] = plan["stored"]
                    if not h[T.CDP]:
                        h[T.CDP] = h[T.EnergySourcePoint] or i + 1
                    if not h[T.CDP_X] and not h[T.CDP_Y]:
                        h[T.CDP_X], h[T.CDP_Y] = h[T.SourceX], h[T.SourceY]
                    if is_bas:
                        for fld in (T.YearDataRecorded, T.DayOfYear, T.HourOfDay,
                                    T.MinuteOfHour, T.SecondOfMinute):
                            h[fld] = 0
                    g.header[i] = h
                    g.trace[i] = np.ascontiguousarray(block[j])
    print(f"  {n_tr} traces x {plan['n']} samples, interval {plan['stored']} "
          f"{plan['uname']} -> {out_path}")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("inputs", nargs="+", help="NetCDF or SEG-Y file(s)")
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
    p.add_argument("--axis", choices=list(AXIS_UNIT_S), default=None,
                   help="what a time axis labelled ms shows: us (default, as BAS), "
                        "ns, or true ms (radar intervals are too short for ms)")
    p.add_argument("--inexact", choices=["resample", "round", "fail"], default="resample",
                   help="what to do if the interval is not a whole header unit")
    p.add_argument("--to-db", action="store_true", help="write 10*log10 of linear power")
    p.add_argument("--coord-scalar", type=int, default=None,
                   help="SEG-Y coordinate scalar (default -100 = cm; 1 = whole metres)")
    p.add_argument("--kingdom", action="store_true",
                   help="settings for Kingdom: --axis us, --coord-scalar 1, and dB for "
                        "BAS linear power (explicit options still win)")
    p.add_argument("--block", type=int, default=2000, help="traces read at a time")
    args = p.parse_args()

    # Fill options the user did not set, from --kingdom or the general defaults.
    if args.axis is None:
        args.axis = "us"
    if args.coord_scalar is None:
        args.coord_scalar = 1 if args.kingdom else -100
    if args.kingdom:
        print(f"Kingdom settings: axis {args.axis}, coordinate scalar {args.coord_scalar}, "
              f"BAS linear power written in dB" + ("" if not args.to_db else "; --to-db for all"))

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    failed = []
    for src in args.inputs:
        # A problem with one file is reported and the rest are still converted.
        try:
            convert_file(src, outdir, args)
        except (ValueError, KeyError, OSError) as err:
            msg = err.args[0] if err.args else repr(err)
            print(f"  ERROR: {msg}")
            failed.append(src)
    if failed:
        raise SystemExit(f"{len(failed)} of {len(args.inputs)} file(s) not converted: "
                         + ", ".join(failed))


def convert_file(src, outdir, args):
    """Convert one input file (NetCDF or SEG-Y) into one or more SEG-Y files."""
    if Path(src).suffix.lower() in (".sgy", ".segy"):
        print(f"{src} (SEG-Y input)")
        out = outdir / f"{Path(src).stem}_{args.axis}axis.sgy"
        if out.resolve() == Path(src).resolve():
            raise ValueError("Output would overwrite the input; choose another --outdir.")
        convert_segy(src, str(out), args)
        return
    with xr.open_dataset(src) as ds:
        name = args.preset or detect_preset(ds)
        cfg = dict(PRESETS[name]) if name else {"crs": "EPSG:3031"}
        for key in ("trace_dim", "sample_dim", "sample_var", "x", "y", "lat", "lon",
                    "crs", "elev", "surface", "time"):
            if getattr(args, key):
                cfg[key] = getattr(args, key)
        data = args.data or cfg.get("data")
        if not data or not cfg.get("trace_dim") or not cfg.get("sample_dim"):
            raise ValueError("layout not recognised; give --data, --trace-dim, "
                             "--sample-dim and --sample-var.")
        cfg.setdefault("sample_var", cfg["sample_dim"])
        print(f"{src} (layout: {name or 'custom'})")
        for d in data:
            convert(ds, src, d, str(outdir / f"{Path(src).stem}_{d}.sgy"), cfg, args)


if __name__ == "__main__":
    main()
