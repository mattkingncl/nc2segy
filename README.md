# nc2segy — radar NetCDF to SEG-Y

`nc2segy.py` converts airborne ice-penetrating-radar profiles to SEG-Y revision 1, so they can be loaded into seismic-interpretation software. It reads NetCDF files, and it can also repair existing radar SEG-Y files whose time axis is wrong, such as those supplied by BAS. Each radar trace (one position along the flight line) becomes one SEG-Y trace. The fast-time axis (two-way travel time) becomes the SEG-Y sample axis. Samples are written as 32-bit IEEE floats (format code 5).

These inputs are recognised automatically:

| Input | Example file | What is written |
|---|---|---|
| BAS / UK Polar Data Centre PASIN NetCDF | `WISE_ISODYN_W15.nc` | `chirp_data`, `pulse_data` |
| UTIG / ICECAP HiCARS L1B NetCDF | `ER2HI1B_2016009_ASB1_JKB2n_R44a_000.nc` | `amplitude_high_gain`, `amplitude_low_gain` |
| BAS-supplied PASIN SEG-Y | `W15_chirp.segy` | The same traces with a corrected time axis |

Other NetCDF layouts can be converted by naming the variables on the command line (see [Other file layouts](#other-file-layouts)).

## Requirements

Python 3 with `numpy`, `pandas`, `xarray`, `netCDF4` and `segyio`. You also need `pyproj` when positions are given only as latitude and longitude, as in the UTIG files.

```
pip install -r requirements.txt
```

## Quick start

```
python nc2segy.py WISE_ISODYN_W15.nc ER2HI1B_2016009_ASB1_JKB2n_R44a_000.nc --outdir segy/
python nc2segy.py W15_chirp.segy --outdir segy/
```

Each NetCDF radar variable goes to its own file, named `<input stem>_<variable>.sgy`. A SEG-Y input is written as `<input stem>_<axis>axis.sgy`, for example `W15_chirp_usaxis.sgy`. Input files are never overwritten.

## Time axis

SEG-Y stores the sample interval as a whole number of microseconds, which is far too coarse for radar. The converter therefore redefines the header unit. `--axis` sets what a time axis labelled "ms" in your interpretation software will show.

| `--axis` | Header unit | A time axis labelled "ms" shows | A 40 µs bed echo appears at |
|---|---|---|---|
| `us` (default, as BAS) | 1 ns | µs | 40 "ms" |
| `ns` | 1 ps | ns | 40 000 "ms" |
| `ms` | 1 µs (true SEG-Y) | ms | Refused for radar |

Any velocity entered in the software must allow for the chosen factor.

**Why `ms` is refused.** The PASIN interval of 45.45 ns would have to become 1 µs, discarding about 95 % of the samples. The script stops with an error rather than do that.

**Upsampling for `ns`.** The interval field holds at most 32 767 units, which is 32.767 ns when the unit is 1 ps. The PASIN interval of 45.45 ns is too long for that, so the trace is upsampled by a whole factor (×2, to 22.727 ns). Every original sample is kept unchanged at every second position, and linear interpolation fills the points in between. The stored 22 727 ps interval is 0.0012 % short, which is 0.76 ns at the last sample. The UTIG interval of 20 ns fits as 20 000 ps, so UTIG traces are written unchanged.

**Intervals that aren't a whole number of header units.** With `--axis us`, the PASIN 22 MHz interval (45.4545 ns) can't be stored exactly, and `--inexact` decides what happens:

| `--inexact` | What it does | Travel times | Amplitudes |
|---|---|---|---|
| `resample` (default) | Linearly interpolates onto the nearest whole-ns interval (45 ns) | Correct | Slightly smoothed (see [Verification](#verification)) |
| `round` | Keeps samples, rounds the header value (as BAS did) | 1.00 % short for PASIN | Unchanged |
| `fail` | Stops with an error | — | — |

With `round`, the −1.00 % error is about 0.40 µs at a 40 µs bed. At 168 m/µs in ice, that is roughly 34 m of apparent ice thickness. If you need both correct timing and untouched amplitudes, use `--axis ns`, which keeps every original sample.

## Example runs

These are real runs on the example files; long warning lines are shortened to the message. The outputs were checked as described under [Verification](#verification).

### 1. Default conversion of the NetCDF files

```
$ python nc2segy.py WISE_ISODYN_W15.nc ER2HI1B_2016009_ASB1_JKB2n_R44a_000.nc --outdir segy
WISE_ISODYN_W15.nc (layout: bas)
    using digitiser interval 45.4545 ns instead of the axis spacing 45.4870 ns
UserWarning: File gives time of day but no date; year and day-of-year are left at zero. Use --date YYYY-MM-DD to fill them.
UserWarning: Units say 'Power in decibel-milliwatts (dBm)', but values span 0.413 to 4.61e+06, which looks like linear power, not dB. Use --to-db to write 10*log10 values.
  chirp_data: 12278 traces x 1414 samples, interval 45 ns -> segy/WISE_ISODYN_W15_chirp_data.sgy
    using digitiser interval 45.4545 ns instead of the axis spacing 45.4870 ns
UserWarning: Units say 'Power in decibel-milliwatts (dBm)', but values span 1.11e-05 to 2.21, which looks like linear power, not dB. Use --to-db to write 10*log10 values.
  pulse_data: 12278 traces x 1414 samples, interval 45 ns -> segy/WISE_ISODYN_W15_pulse_data.sgy
ER2HI1B_2016009_ASB1_JKB2n_R44a_000.nc (layout: utig)
  amplitude_high_gain: 4000 traces x 3200 samples, interval 20 ns -> segy/ER2HI1B_2016009_ASB1_JKB2n_R44a_000_amplitude_high_gain.sgy
  amplitude_low_gain: 4000 traces x 3200 samples, interval 20 ns -> segy/ER2HI1B_2016009_ASB1_JKB2n_R44a_000_amplitude_low_gain.sgy
```

The BAS traces grow from 1400 to 1414 samples because they are resampled from 45.45 ns to 45 ns over the same 63.6 µs record. The UTIG interval of 20 ns is stored exactly, so its samples are unchanged.

### 2. Repairing a BAS-supplied SEG-Y file

```
$ python nc2segy.py W15_chirp.segy --outdir segy
W15_chirp.segy (SEG-Y input)
  BAS PASIN SEG-Y recognised (trace interval 45, binary 1000); using 45.4545 ns (22 MHz)
  12278 traces x 1414 samples, interval 45 ns -> segy/W15_chirp_usaxis.sgy

$ python nc2segy.py W15_chirp.segy --axis ns --outdir segy
W15_chirp.segy (SEG-Y input)
  BAS PASIN SEG-Y recognised (trace interval 45, binary 1000); using 45.4545 ns (22 MHz)
  12278 traces x 2799 samples, interval 22727 ps -> segy/W15_chirp_nsaxis.sgy
```

A BAS PASIN SEG-Y file is recognised by its two conflicting header intervals: 45 in the trace headers and 1000 in the binary header. Its true interval is then taken as 45.4545 ns. The repair:

* resamples or upsamples the traces as set by `--axis` and `--inexact`;
* writes the same interval in the binary and trace headers;
* copies every other trace header (coordinates, PRI number, trace number) unchanged;
* fills the empty CDP number and CDP x/y from the source fields;
* sets the year, day and time fields to zero, because BAS filled them with the file-creation date;
* writes a readable textual header describing all of this.

For any other SEG-Y file the true interval can't be inferred, so give it with `--dt-ns`.

### 3. Reproducing the BAS-supplied SEG-Y from NetCDF

```
$ python nc2segy.py WISE_ISODYN_W15.nc --data chirp_data --inexact round --coord-scalar 1 --outdir bas_style
WISE_ISODYN_W15.nc (layout: bas)
    using digitiser interval 45.4545 ns instead of the axis spacing 45.4870 ns
UserWarning: Header interval rounded from 45.4545 to 45 ns: every travel time read from the SEG-Y is off by -1.00% (-0.636 us at the last sample).
...
  chirp_data: 12278 traces x 1400 samples, interval 45 ns -> bas_style/WISE_ISODYN_W15_chirp_data.sgy
```

This keeps the native samples and whole-metre coordinates, like `W15_chirp.segy`, and the trace data are bit-identical to that file. The warning states the timing error that this choice carries.

### 4. Other common variations

```
# Time axis in nanoseconds
python nc2segy.py ER2HI1B_2016009_ASB1_JKB2n_R44a_000.nc --axis ns --outdir segy

# Write BAS power in decibels, and fill in the flight date (the file has none)
python nc2segy.py WISE_ISODYN_W15.nc --to-db --date YYYY-MM-DD --outdir segy

# Repair a whole folder of BAS SEG-Y files
python nc2segy.py bas_segy/*.segy --outdir segy

# Refuse to write anything if the interval cannot be stored exactly
python nc2segy.py WISE_ISODYN_W15.nc --inexact fail
```

Replace `YYYY-MM-DD` with the real flight date; I have not been able to verify the W15 flight date. In Windows PowerShell, wildcards such as `*.segy` are not expanded for Python. Run `python nc2segy.py (Get-ChildItem bas_segy\*.segy).FullName --outdir segy` instead.

## Trace header layout (NetCDF input)

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
| 115, 117 | Samples, interval | As in the binary header |
| 157–165 | Year, day, hour, minute, second | UTC; year and day are 0 for BAS unless `--date` is given |
| 167 | Time basis code | 4 (UTC) |

BAS files supply projected x/y directly. UTIG latitude and longitude are projected to EPSG:3031 with pyproj. BAS missing-value flags (−9999) and non-finite values are left at 0 in the headers, and the script reports how many traces were affected. Non-finite samples are written as 0 and counted.

The textual header (EBCDIC, 40 lines) records the source file, variable, units, time axis, any resampling and its header error, the coordinate system and this byte layout.

## Verification

All of these checks were run on the example files.

* **BAS NetCDF, `--inexact round --coord-scalar 1`.** The trace data are bit-identical to `W15_chirp.segy`. PRI numbers and trace numbers match exactly, and x/y agree within 1 m. The 1 m differences come from rounding, because the reference appears to round differently.
* **BAS NetCDF, default.** Resampled samples match an independent `numpy.interp` to a relative error of 3 × 10⁻⁶ (float32 precision). Coordinates match the NetCDF exactly, elevations within 0.5 cm, and time of day within the 1 s truncation that SEG-Y imposes.
* **Effect of resampling on bed echoes (`--axis us`).** Across 9394 picked bed echoes, the peak power changes by a median of −0.04 dB. The 5th percentile is −0.52 dB, the 1st percentile is −1.0 dB, and the worst single trace is −3.1 dB. The cause is that linear interpolation lowers sharp peaks sampled barely above the Nyquist rate. Use `--axis ns` or `--inexact round` if amplitudes matter.
* **`--axis ns`.** For BAS, every second output sample is bit-identical to the original samples. For UTIG, the output is identical to the NetCDF.
* **BAS SEG-Y repair.** For both `--axis us` and `--axis ns`, the repaired `W15_chirp.segy` traces are bit-identical to those converted from `WISE_ISODYN_W15.nc`. Coordinates, PRI numbers and trace numbers are copied unchanged.
* **UTIG.** The samples are bit-identical to the NetCDF. x/y agree within 5 mm of an independent pyproj projection, and elevations within 0.5 cm. The UTC date and time are written in full.
* **Headers.** The binary and trace-header sample intervals agree, and the revision is set to 1. segyio reads the time axis as intended (step 0.045 for `us`, 22.727 for `ns`).

You can inspect any output with segyio:

```python
import segyio

# Open without assuming a 3-D survey geometry
with segyio.open("segy/W15_chirp_usaxis.sgy", ignore_geometry=True) as f:
    print(segyio.tools.wrap(f.text[0]))           # textual header
    print(f.bin[segyio.BinField.Interval])        # stored sample interval
    print(f.samples[:5])                          # time axis as software shows it
    print(f.header[0])                            # first trace header
    data = f.trace.raw[:]                         # traces x samples array
```

## Notes on the source data

The converter detects or works around each of the following problems. These are my findings from the example files, not statements from BAS or UTIG.

1. **The BAS-supplied SEG-Y has the wrong time scale.** `W15_chirp.segy` stores a 45 ns interval for 45.45 ns data, so travel times read from it are 1.0 % short. Its binary header gives an interval of 1000 while the traces give 45, so different programs read different intervals; segyio, for example, falls back to 4 ms. Its textual header is unreadable. Repair such files as shown in example 2.
2. **The BAS `fast_time` axis is slightly stretched.** It runs from 0 to N·dt over N samples, which gives a spacing of 45.487 ns. The digitiser runs at 22 MHz (45.4545 ns; the metadata round this to 45.4546 ns). A fit of BAS's own ice thicknesses against the pick separation gives 3.818 m per sample, i.e. 84 m/µs × 45.454 ns, with a residual standard deviation of 0.14 m. For BAS files the converter therefore uses 1/22 MHz. `--dt-ns` overrides this.
3. **BAS radar data labelled "dBm" are linear power.** Chirp values span 0.41 to 4.6 × 10⁶. Their 10·log10 values fit the 20–80 dB grey scale in `WISE_ISODYN_W15.pdf`. The script warns about this, and `--to-db` converts.
4. **BAS picks labelled "microseconds" are sample numbers.** Bed picks reach 871 within a 63.6 µs record. In the PDF, 871 samples × 45.45 ns ≈ 39.6 µs matches the deepest bed.
5. **The BAS thickness offset is unexplained.** The fit in note 2 has a 23.8 m intercept, whereas the metadata state a 10 m firn correction. I can't account for the difference; pick conventions may play a part.
6. **The BAS SEG-Y dates are not acquisition dates.** They read 2021, day 200, but the survey was flown in 2005/06. Their times of day don't match `UTC_time_layerData` either, so the repair sets them to zero. The BAS NetCDF gives time of day only, so year and day stay 0 unless you give `--date`.
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
--preset {bas,utig}   NetCDF layout (default: detect)
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
--dt-ns DT_NS         true sample interval (ns); needed for non-BAS SEG-Y input
--axis {ms,us,ns}     what a time axis labelled "ms" shows (default us, as BAS)
--inexact {resample,round,fail}
                      what to do if the interval is not a whole header unit
--to-db               write 10*log10 of linear power
--coord-scalar N      SEG-Y coordinate scalar (default -100 = cm; 1 = whole metres)
--block BLOCK         traces read at a time (memory use)
```

## Limitations

* **Size limits.** SEG-Y rev 1 allows at most 32 767 samples per trace and an interval of at most 32 767 header units. segyio reads larger values as negative, so they are never written.
* **Regular sampling only.** Fast time must be evenly spaced; residual jitter must be under 1 % of an interval.
* **Two dimensions only.** Data with extra dimensions, such as channels, need converting to separate 2-D variables first.
* **No picks in the output.** Surface and bed picks have no standard SEG-Y field and are not written. Load them into your software as horizons from the NetCDF.
* **SEG-Y repair assumes BAS conventions.** It is tested only on `W15_chirp.segy`. Delay-recording times in other SEG-Y files are ignored.
* **Untested in commercial software.** I have not tested the output in Petrel, Kingdom or OpendTect. Check one line in your software before converting a whole survey.
