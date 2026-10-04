"""
plume_area_timeseries_figure.py
===============================

Makes ONE figure: plume area through time (2022-2025) for Sentinel-1 and
Sentinel-2, in two stacked panels.

    - each dot          = one satellite scene with a plume
    - the solid line    = 10-day rolling mean of plume area
    - grey bands        = wet-weather days: rain >= 0.1 inch plus the 3 days
                          after it, or high river flow

    plume_areas_s1_s2.csv (step 07) ──► scene dots + 10-day rolling mean ─┐
    IBWC daily flow  ─┐                                                    ├─► figure
    TJ NERR rainfall ─┴──────────────► wet-weather bands ─────────────────┘

No statistics linking plume area to flow or rain are calculated here; flow and
rain are only used to decide where the grey bands go.

Inputs
------
1. plume_areas_s1_s2.csv   plume area of every mask (written by 07_measure_plumes.py)
2. IBWC daily river flow   (CSV, Million Gallons per day)
3. TJ NERR met station     (CSV, 15-minute precipitation)

Output
------
fig_plume_area_timeseries.png, saved next to this script (in figures/)

Needs only: numpy, pandas, matplotlib
Run with:   python plume_area_timeseries_figure.py   (after 07_measure_plumes.py)
"""

import os
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")                     # write to file, no window
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


# =============================================================================
# 1. FILE LOCATIONS
# =============================================================================
# TIDE server if /home/jovyan/s2 exists, otherwise the Mac external drive
# (set the TJ_ROOT environment variable if the drive is mounted elsewhere).
TIDE = Path("/home/jovyan/s2")
if TIDE.exists():
    AREAS_CSV  = TIDE / "07_figures" / "plume_areas_s1_s2.csv"
    FLOW_CSV   = TIDE / "ancillary" / "Total Flow.Preliminary Daily.US Million Gallons.csv"
    PRECIP_CSV = TIDE / "ancillary" / "348805.csv"
else:
    TJ = Path(os.environ.get("TJ_ROOT", "/Volumes/External/TJ"))
    AREAS_CSV  = TJ / "010_optical" / "07_figures" / "plume_areas_s1_s2.csv"
    FLOW_CSV   = TJ / "011_environmentalData" / "IBWC" / "Total Flow.Preliminary Daily.US Million Gallons.csv"
    PRECIP_CSV = TJ / "011_environmentalData" / "TJNERR_metData" / "348805.csv"
OUT_PNG    = Path(__file__).resolve().parent / "fig_plume_area_timeseries.png"


# =============================================================================
# 2. SETTINGS (everything you might want to change is here)
# =============================================================================
FIRST_YEAR, LAST_YEAR = 2022, 2025

# Rolling average line
WINDOW = "10D"           # window length, centred on each scene (5 days either side)
MIN_SCENES = 2           # scenes needed inside the window; fewer -> the line breaks

# A day counts as "wet weather" if EITHER is true:
#   (a) it rained >= WET_RAIN_IN that day, or within the previous RAIN_FOLLOWING_DAYS days
#       (the wet-weather definition used in southern California beach bacteria TMDLs)
#   (b) daily river flow was above WET_FLOW_MGD (about the 90th percentile, 2022-2025)
WET_RAIN_IN = 0.1                # inches of rain in one day
RAIN_FOLLOWING_DAYS = 3          # days after a rain day that still count as wet
WET_FLOW_MGD = 100.0             # Million Gallons/day
MM_PER_INCH = 25.4               # the rain data are in mm

# Precipitation quality control (as in the original step 07 figures)
GOOD_QC_FLAGS = {"0", "1", "4", "5"}     # SWMP flags kept: passed, suspect, historical, corrected
MIN_DAY_COVERAGE = 0.8                   # a day needs >= 80% of its 96 15-min records

# Colours (blue/orange is a colour-blind-safe pair)
S1_COLOR = "#2a6fb0"
S2_COLOR = "#e07b00"
WET_COLOR = "0.86"       # light grey


# =============================================================================
# 3. LOAD THE DATA
# =============================================================================
def load_plume_areas():
    """Plume area (km²) per scene. Keeps scenes with a plume (area > 0)."""
    df = pd.read_csv(AREAS_CSV)
    df["Timestamp"] = pd.to_datetime(df["Timestamp"], utc=True).dt.tz_localize(None)
    df = df[df["Area_km2"] > 0]
    df = df[df["Timestamp"].dt.year.between(FIRST_YEAR, LAST_YEAR)]
    return df


def load_daily_flow():
    """Daily river flow (Million Gallons/day), one value per date."""
    df = pd.read_csv(FLOW_CSV)
    time_col = next(c for c in df.columns if c.strip().startswith("Timestamp"))
    value_col = next(c for c in df.columns if c.strip().startswith("Value"))
    dates = pd.to_datetime(df[time_col], format="%m/%d/%y %H:%M").dt.floor("D")
    flow = pd.Series(df[value_col].values, index=dates)
    return flow.groupby(level=0).mean()


def load_daily_rain():
    """Daily rainfall total (mm) from 15-minute records, after quality control."""
    df = pd.read_csv(PRECIP_CSV, usecols=["DateTimeStamp", "TotPrcp", "F_TotPrcp"])

    # keep only records whose QC flag (e.g. "<0>") is in GOOD_QC_FLAGS
    flag = df["F_TotPrcp"].astype(str).str.extract(r"<(-?\d+)>")[0]
    df = df[flag.isin(GOOD_QC_FLAGS)]

    # sum each day; blank out days with too few records
    dates = pd.to_datetime(df["DateTimeStamp"], format="%m/%d/%Y %H:%M").dt.floor("D")
    daily = df.groupby(dates)["TotPrcp"].agg(["sum", "count"])
    daily.loc[daily["count"] < MIN_DAY_COVERAGE * 96, "sum"] = np.nan
    return daily["sum"]


# =============================================================================
# 4. FIND WET-WEATHER PERIODS
# =============================================================================
def find_wet_periods(flow, rain):
    """Return a list of (start_date, end_date) for runs of consecutive wet days."""
    all_days = pd.date_range(f"{FIRST_YEAR}-01-01", f"{LAST_YEAR}-12-31", freq="D")
    # (a) rain days, then extend each one forward by RAIN_FOLLOWING_DAYS days
    rain_day = (rain.reindex(all_days) >= WET_RAIN_IN * MM_PER_INCH).astype(int)
    rain_wet = rain_day.rolling(RAIN_FOLLOWING_DAYS + 1, min_periods=1).max() == 1

    # (b) high-flow days
    flow_wet = flow.reindex(all_days) > WET_FLOW_MGD

    is_wet = rain_wet | flow_wet

    periods = []
    start = None
    for day, wet in zip(all_days, is_wet):
        if wet and start is None:              # a wet run begins
            start = day
        elif not wet and start is not None:    # the wet run ended yesterday
            periods.append((start, day))
            start = None
    if start is not None:                      # still wet on the last day
        periods.append((start, all_days[-1] + pd.Timedelta(days=1)))
    return periods


# =============================================================================
# 5. ROLLING AVERAGE
# =============================================================================
def rolling_mean(areas_one_sensor):
    """10-day rolling mean of plume area, evaluated at each scene's date.

    For each scene, average all scenes of the same sensor within +/- 5 days.
    If fewer than MIN_SCENES scenes fall in that window the value is left
    blank (NaN), which shows up as a break in the line.
    """
    series = areas_one_sensor.set_index("Timestamp")["Area_km2"].sort_index()
    return series.rolling(WINDOW, center=True, min_periods=MIN_SCENES).mean()


# =============================================================================
# 6. PLOT
# =============================================================================
def draw_panel(ax, scenes, line, color, marker, title, wet_periods):
    """One panel: grey wet-weather bands, faint scene dots, rolling-mean line."""
    for start, end in wet_periods:
        ax.axvspan(start, end, color=WET_COLOR, lw=0.6, zorder=0)   # lw keeps 1-day bands visible

    ax.scatter(scenes["Timestamp"], scenes["Area_km2"], s=12, marker=marker,
               color=color, alpha=0.3, lw=0, zorder=2)
    ax.plot(line.index, line.values, color=color, lw=2, zorder=3)

    ax.set_ylim(0, scenes["Area_km2"].max() * 1.08)
    ax.set_ylabel("Plume area (km²)")
    ax.text(0.005, 0.97, f"{title} (n = {len(scenes)} scenes)", transform=ax.transAxes,
            ha="left", va="top", fontweight="bold", color=color)

    # thin vertical line at each new year
    for year in range(FIRST_YEAR + 1, LAST_YEAR + 1):
        ax.axvline(pd.Timestamp(f"{year}-01-01"), color="0.55", lw=0.8, zorder=1)

    ax.grid(axis="y", color="0.92", lw=0.8, zorder=0)
    ax.spines[["top", "right"]].set_visible(False)


def main():
    # --- load ---
    areas = load_plume_areas()
    s1 = areas[areas["Sensor"] == "S1"]
    s2 = areas[areas["Sensor"] == "S2"]
    wet_periods = find_wet_periods(load_daily_flow(), load_daily_rain())

    # --- figure: two panels sharing the time axis ---
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(13, 6), sharex=True, layout="constrained")
    draw_panel(ax1, s1, rolling_mean(s1), S1_COLOR, "o", "Sentinel-1", wet_periods)
    draw_panel(ax2, s2, rolling_mean(s2), S2_COLOR, "s", "Sentinel-2, cloud-free", wet_periods)

    # --- time axis: month ticks, year labels above the top panel ---
    ax2.set_xlim(pd.Timestamp(f"{FIRST_YEAR}-01-01"), pd.Timestamp(f"{LAST_YEAR + 1}-01-01"))
    ax2.xaxis.set_major_locator(mdates.MonthLocator(bymonth=(1, 4, 7, 10)))
    ax2.xaxis.set_major_formatter(mdates.DateFormatter("%b"))
    for year in range(FIRST_YEAR, LAST_YEAR + 1):
        ax1.text(pd.Timestamp(f"{year}-07-02"), ax1.get_ylim()[1], str(year),
                 ha="center", va="bottom", fontweight="bold", color="0.25")

    # --- legend ---
    legend_items = [
        Line2D([], [], ls="", marker="o", ms=5, color="0.4", alpha=0.4, label="Individual scenes"),
        Line2D([], [], color="0.3", lw=2, label=f"{WINDOW[:-1]}-day rolling mean"),
        Patch(color=WET_COLOR, label=f"Wet weather (rain ≥ {WET_RAIN_IN} in + {RAIN_FOLLOWING_DAYS} days after, "
                                     f"or flow > {WET_FLOW_MGD:.0f} MGD)"),
    ]
    fig.legend(handles=legend_items, loc="outside upper center", ncol=3, frameon=False)

    # --- save ---
    fig.savefig(OUT_PNG, dpi=300)
    plt.close(fig)
    print(f"S1: {len(s1)} scenes, S2: {len(s2)} scenes, {len(wet_periods)} wet-weather periods")
    print(f"Saved {OUT_PNG}")


if __name__ == "__main__":
    main()
