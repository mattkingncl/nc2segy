# nc2segy — radar NetCDF to SEG-Y

`nc2segy.py` converts airborne ice-penetrating-radar profiles from NetCDF to SEG-Y revision 1, so they can be loaded into seismic-interpretation software. Each radar trace (one position along the flight line) becomes one SEG-Y trace. The fast-time axis (two-way travel time) becomes the SEG-Y sample axis. Samples are written as 32-bit IEEE floats (format code 5).

Two file layouts are recognised automatically:

| Layout | Example file | Radar variables written |
|---|---|---|
| `bas` | BAS / UK Polar Data Centre PASIN, e.g. `WISE_ISODYN_W15.nc` | `chirp_data`, `pulse_data` |
| `utig` | UTIG / ICECAP HiCARS L1B, e.g. `ER2HI1B_2016009_ASB1_JKB2n_R44a_000.nc` | `amplitude_high_gain`, `amplitude_low_gain` |

Other layouts can be converted by naming the variables on the command line (see [Other file layouts](#other-file-layouts)).

## Requirements

Python 3 with `numpy`, `pandas`, `xarray`, `netCDF4` and `segyio`. You also need `pyproj` when positions are given only as latitude and longitude, as in the UTIG files.

```
pip install numpy pandas xarray netCDF4 segyio pyproj
```

## Quick start

```
python nc2segy.py WISE_ISODYN_W15.nc ER2HI1B_2016009_ASB1_JKB2n_R44a_000.nc --outdir segy/
```

Each radar variable goes to its own file, named `<input stem>_<variable>.sgy`.

## Example runs

These are real runs on the two example files. The outputs were checked as described under [Verification](#verification).

### 1. Default conversion of both files

```
$ python nc2segy.py WISE_ISODYN_W15.nc ER2HI1B_2016009_ASB1_JKB2n_R44a_000.nc --outdir segy
WISE_ISODYN_W15.nc (layout: bas)
    using metadata interval 45.4546 ns instead of the axis spacing 45.4870 ns
UserWarning: File gives time of day but no date; year and day-of-year are left at zero. Use --date YYYY-MM-DD to fill them.
UserWarning: Units say 'Power in decibel-milliwatts (dBm)', but values span 0.413 to 4.61e+06, which looks like linear power, not dB. Use --to-db to write 10*log10 values.
  chirp_data: 12278 traces x 1414 samples, interval 45 ns -> segy/WISE_ISODYN_W15_chirp_data.sgy
    using metadata interval 45.4546 ns instead of the axis spacing 45.4870 ns
UserWarning: Units say 'Power in decibel-milliwatts (dBm)', but values span 1.11e-05 to 2.21, which looks like linear power, not dB. Use --to-db to write 10*log10 values.
  pulse_data: 12278 traces x 1414 samples, interval 45 ns -> segy/WISE_ISODYN_W15_pulse_data.sgy
ER2HI1B_2016009_ASB1_JKB2n_R44a_000.nc (layout: utig)
  amplitude_high_gain: 4000 traces x 3200 samples, interval 20 ns -> segy/ER2HI1B_2016009_ASB1_JKB2n_R44a_000_amplitude_high_gain.sgy
  amplitude_low_gain: 4000 traces x 3200 samples, interval 20 ns -> segy/ER2HI1B_2016009_ASB1_JKB2n_R44a_000_amplitude_low_gain.sgy
```

The BAS traces grow from 1400 to 1414 samples because they are resampled from 45.45 ns to 45 ns over the same 63.6 µs record (see [Time convention](#time-convention)). The UTIG interval of 20 ns is stored exactly, so its samples are unchanged.

### 2. Reproducing the BAS-supplied SEG-Y

```
$ python nc2segy.py WISE_ISODYN_W15.nc --data chirp_data --inexact round --coord-scalar 1 --outdir bas_style
WISE_ISODYN_W15.nc (layout: bas)
    using metadata interval 45.4546 ns instead of the axis spacing 45.4870 ns
UserWarning: Header interval rounded from 45.4546 ns to 45 ns: every travel time read from the SEG-Y is off by -1.00% (-0.636 us at the last sample).
...
  chirp_data: 12278 traces x 1400 samples, interval 45 ns -> bas_style/WISE_ISODYN_W15_chirp_data.sgy
```

This keeps the native samples and whole-metre coordinates, like `W15_chirp.segy`. The trace data are bit-identical to that file. The warning states the timing error that this choice carries.

### 3. Other common variations

```
# Write BAS power in decibels, and fill in the flight date (the file has none)
python nc2segy.py WISE_ISODYN_W15.nc --to-db --date YYYY-MM-DD --outdir segy

# Only the high-gain UTIG channel
python nc2segy.py ER2HI1B_2016009_ASB1_JKB2n_R44a_000.nc --data amplitude_high_gain --outdir segy

# Refuse to write anything if the interval cannot be stored exactly
python nc2segy.py WISE_ISODYN_W15.nc --inexact fail
```

Replace `YYYY-MM-DD` with the real flight date. I have not been able to verify the W15 flight date.

## Time convention

SEG-Y stores the sample interval as a whole number of microseconds, which is far too coarse for radar. Following the BAS SEG-Y files, the converter stores the interval in **nanoseconds** by default (`--time-units ns`). Seismic software reads the header as microseconds, so:

> A time axis labelled "ms" in your interpretation software shows **µs** of two-way travel time.

For example, a bed echo at 40 µs appears at "40 ms". Any velocity entered in the software must allow for this factor of 1000.

If the true interval is not a whole number of nanoseconds, `--inexact` decides what happens. The PASIN 22 MHz digitiser gives 45.4545 ns, so this applies to BAS files.

| `--inexact` | What it does | Travel times | Amplitudes |
|---|---|---|---|
| `resample` (default) | Linearly interpolates onto the nearest whole-ns interval | Correct | Slightly smoothed (see below) |
| `round` | Keeps samples, rounds the header value (as BAS did) | Off by the rounding error, −1.00 % for PASIN | Unchanged |
| `fail` | Stops with an error | — | — |

With `round`, the −1.00 % error is about 0.40 µs at a 40 µs bed. At 168 m/µs in ice, that is roughly 34 m of apparent ice thickness.

## Trace header layout

| Bytes | Field | Contents |
|---|---|---|
| 1, 5 | Trace sequence (line, file) | 1, 2, 3 … |
| 9 | Field record | PRI number (BAS); otherwise trace sequence |
| 13 | Trace number | 1 (as in the BAS files) |
| 17, 21 | Energy source point, CDP | Trace number from the file (BAS `traces`), otherwise sequence |
| 29 | Trace identification code | 1 |
| 41 | Receiver group elevation | Antenna (aircraft) altitude, m WGS84 |
| 45 | Surface elevation at source | Ice-surface elevation, m WGS84 (BAS only) |
| 69 | Elevation scalar | −100 (values in cm) |
| 71 | Coordinate scalar | −100 by default (values in cm); set with `--coord-scalar` |
| 73–88, 181–188 | Source, group and CDP x/y | EPSG:3031 polar stereographic, m |
| 89 | Coordinate units | 1 (length) |
| 115, 117 | Samples, interval | As in the binary header, which now agrees with the traces |
| 157–167 | Year, day, hour, minute, second | UTC; year and day are 0 for BAS unless `--date` is given |
| 167 | Time basis code | 4 (UTC) |

BAS files supply projected x/y directly. UTIG latitude and longitude are projected to EPSG:3031 with pyproj. BAS missing-value flags (−9999) and non-finite values are left at 0 in the headers, and the script reports how many traces were affected. Non-finite samples are written as 0 and counted.

The textual header (EBCDIC, 40 lines) records the source file, variable, units, time convention, any resampling, the coordinate system and this byte layout.

## Verification

All of these checks were run on the two example files.

* **BAS, `--inexact round --coord-scalar 1`.** The trace data are bit-identical to `W15_chirp.segy`. PRI numbers and trace numbers match exactly, and x/y agree within 1 m. The 1 m differences come from rounding, because the reference appears to round differently.
* **BAS, default.** Resampled samples match an independent `numpy.interp` to a relative error of 3 × 10⁻⁶ (float32 precision). Coordinates match the NetCDF exactly, elevations within 0.5 cm, and time of day within the 1 s truncation that SEG-Y imposes.
* **Effect of resampling on bed echoes.** Across 9394 picked bed echoes, the peak power changes by a median of −0.04 dB. The 5th percentile is −0.52 dB, the 1st percentile is −1.0 dB, and the worst single trace is −3.1 dB. The cause is that linear interpolation lowers sharp peaks sampled barely above the Nyquist rate. Use `--inexact round` if amplitudes matter more than timing.
* **UTIG.** The samples are bit-identical to the NetCDF. x/y agree within 5 mm of an independent pyproj projection, and elevations within 0.5 cm. The UTC date and time are written in full.
* **Headers.** The binary and trace-header sample intervals agree, and the revision is set to 1.

You can inspect any output with segyio:

```python
import segyio

# Open without assuming a 3-D survey geometry
with segyio.open("segy/WISE_ISODYN_W15_chirp_data.sgy", ignore_geometry=True) as f:
    print(segyio.tools.wrap(f.text[0]))           # textual header
    print(f.bin[segyio.BinField.Interval])        # sample interval (ns)
    print(f.header[0])                            # first trace header
    data = f.trace.raw[:]                         # traces x samples array
```

## Notes on the source data

The converter detects or works around each of the following problems. These are my findings from the example files, not statements from BAS or UTIG.

1. **The reference BAS SEG-Y has the wrong time scale.** `W15_chirp.segy` stores a 45 ns interval for 45.45 ns data, so travel times read from it are 1.0 % short. Its binary header also gives an interval of 1000 while the traces give 45, so different programs may read different intervals. Its textual header is unreadable.
2. **The BAS `fast_time` axis is slightly stretched.** It runs from 0 to N·dt over N samples, which gives a spacing of 45.487 ns. The metadata quote 45.4546 ns (22 MHz). A fit of BAS's own ice thicknesses against the pick separation gives 3.818 m per sample, i.e. 84 m/µs × 45.454 ns, with a residual standard deviation of 0.14 m. For BAS files the converter therefore uses the metadata interval. `--dt-ns` overrides this.
3. **BAS radar data labelled "dBm" are linear power.** Chirp values span 0.41 to 4.6 × 10⁶. Their 10·log10 values fit the 20–80 dB grey scale in `WISE_ISODYN_W15.pdf`. The script warns about this, and `--to-db` converts.
4. **BAS picks labelled "microseconds" are sample numbers.** Bed picks reach 871 within a 63.6 µs record. In the PDF, 871 samples × 45.45 ns ≈ 39.6 µs matches the deepest bed.
5. **The BAS thickness offset is unexplained.** The fit in note 2 has a 23.8 m intercept, whereas the metadata state a 10 m firn correction. I can't account for the difference; pick conventions may play a part.
6. **The reference SEG-Y's dates are not acquisition dates.** They read 2021, day 200, but the survey was flown in 2005/06. Their times of day don't match `UTC_time_layerData` either. The BAS NetCDF gives time of day only, so year and day stay 0 unless you give `--date`.
7. **The BAS map-projection text disagrees with the coordinates.** The metadata list a false northing of 2 082 760 m. However, `x_coordinates` and `y_coordinates` agree with standard EPSG:3031 (false northing 0) to within 5 mm, so EPSG:3031 is used.
8. **UTIG fast time starts at the start of the record, not at transmission.** The metadata give a transmit delay of 2.33 µs, so zero travel time sits at 2.33 µs on the SEG-Y time axis. This is recorded in the textual header, but the data are not shifted.
9. **UTIG amplitudes are already in dB** ("counts in dB"), so do not use `--to-db` on them.

## Other file layouts

Name the radar variable, its two dimensions and the fast-time variable. Add whichever navigation variables exist:

```
python nc2segy.py other.nc --data echo --trace-dim trace --sample-dim sample \
    --sample-var twtt --lat latitude --lon longitude --elev height --time time --outdir segy
```

The fast-time variable needs a `units` attribute in s, ms, µs (`microseconds`), ns or ps, and must be evenly spaced. If it has no usable units, give `--dt-ns`.

## Options

```
--outdir OUTDIR       output folder (default: current)
--preset {bas,utig}   file layout (default: detect)
--data DATA           radar variable(s); repeat for several; default: all in preset
--trace-dim, --sample-dim, --sample-var
                      dimension and variable names for other layouts
--x, --y              projected coordinates (m)
--lat, --lon          geographic coordinates, projected to --crs
--crs CRS             projected CRS (default EPSG:3031)
--elev ELEV           antenna altitude variable (m)
--surface SURFACE     ice-surface elevation variable (m)
--time TIME           per-trace timestamp variable
--date DATE           flight date YYYY-MM-DD when the file gives only time of day
--dt-ns DT_NS         override the sample interval (ns)
--time-units {us,ns,ps}
                      meaning of one header time unit (default ns, as BAS)
--inexact {resample,round,fail}
                      what to do if the interval is not a whole header unit
--to-db               write 10*log10 of linear power
--coord-scalar N      SEG-Y coordinate scalar (default -100 = cm; 1 = whole metres)
--block BLOCK         traces read at a time (memory use)
```

## Limitations

* **Size limits.** SEG-Y rev 1 allows at most 32 767 samples per trace and an interval of at most 32 767 header units.
* **Regular sampling only.** Fast time must be evenly spaced; residual jitter must be under 1 % of an interval.
* **Two dimensions only.** Data with extra dimensions, such as channels, need converting to separate 2-D variables first.
* **No picks in the output.** Surface and bed picks have no standard SEG-Y field and are not written. Load them into your software as horizons from the NetCDF.
* **Untested in commercial software.** I have not tested the output in Petrel, Kingdom or OpendTect. Check one line in your software before converting a whole survey.
