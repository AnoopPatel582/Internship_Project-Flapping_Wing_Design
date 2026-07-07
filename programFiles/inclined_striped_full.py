import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import Polygon, Patch
from matplotlib.lines import Line2D
from scipy.interpolate import make_interp_spline

# HARD-CODED PARAMETERS
MEM_HEIGHT      = 6.0    # mm
SPAR_HEIGHT     = 1.0    # mm
SPAR_BOTTOM     = 6.0    # mm
CHORD           = 7.0    # mm
THETA           = np.radians(35.0)
STIFF_WIDTH     = 1.0    # mm
NUM_SUBSTRIPS   = 10
NUM_PAIRS       = 3

# Strip dimensions
TYPE_A_WIDTH    = 2.0    # mm

# Derived geometry
STIFF_X_BOTTOM  = MEM_HEIGHT / np.tan(THETA)
STIFF_TOP_X     = STIFF_WIDTH / np.sin(THETA)
TYPE_B_WIDTH    = STIFF_X_BOTTOM + STIFF_TOP_X
WINGSPAN        = NUM_PAIRS * (TYPE_A_WIDTH + TYPE_B_WIDTH)

# Material thicknesses
SPAR_THICKNESS  = 0.5
MEM_THICKNESS   = 0.05
STIFF_THICKNESS = 0.5

# Material densities
SPAR_DENSITY    = 1.6
MEM_DENSITY     = 1.3
STIFF_DENSITY   = 1.6

# Young's Modulus
E_CARBON        = 80.0
E_MEMBRANE      = 5.0

# HELPER FUNCTIONS
def x_left(y):
    return STIFF_X_BOTTOM * (1 - y / MEM_HEIGHT)

def x_right(y):
    return STIFF_TOP_X + STIFF_X_BOTTOM * (1 - y / MEM_HEIGHT)

def clamp(val, lo=0.0, hi=MEM_HEIGHT):
    return max(lo, min(hi, val))

def y_transition(x_val, boundary='left'):
    if boundary == 'left':
        y = MEM_HEIGHT * (1 - x_val / STIFF_X_BOTTOM)
    else:
        y = MEM_HEIGHT * (1 - (x_val - STIFF_TOP_X) / STIFF_X_BOTTOM)
    return clamp(y)

# REGION CALCULATOR
def calc_region(x_s, x_e, y_b, y_top, region_type, global_x_offset=0.0):
    h      = y_top - y_b
    w      = x_e - x_s
    cx_mid = global_x_offset + x_s + w / 2.0
    cy_mid = y_b + h / 2.0

    if region_type in ('mem_rect', 'right_mem_rect'):
        area = w * h
        mass = area * MEM_THICKNESS * MEM_DENSITY
        EA   = E_MEMBRANE * MEM_THICKNESS * w
        EAy  = EA * cy_mid
        return mass, cx_mid, cy_mid, EA, EAy

    elif region_type == 'stiff_rect':
        area = w * h
        mass = area * STIFF_THICKNESS * STIFF_DENSITY
        EA   = E_CARBON * STIFF_THICKNESS * w
        EAy  = EA * cy_mid
        return mass, cx_mid, cy_mid, EA, EAy

    elif region_type == 'mem_tri_left':
        xl_top    = x_left(y_top)
        xl_bot    = x_left(y_b)
        mw_top    = max(0, xl_top - x_s)
        mw_bot    = max(0, xl_bot - x_s)
        sw_top    = w - mw_top
        sw_bot    = w - mw_bot
        mem_area  = 0.5 * (mw_top + mw_bot) * h
        stiff_area = w * h - mem_area

        if (mw_top + mw_bot) > 1e-10:
            mem_cg_y   = y_b + h*(2*mw_top + mw_bot) / (3*(mw_top + mw_bot))
        else:
            mem_cg_y   = cy_mid

        if (sw_top + sw_bot) > 1e-10:
            stiff_cg_y = y_b + h*(2*sw_top + sw_bot) / (3*(sw_top + sw_bot))
        else:
            stiff_cg_y = cy_mid

        mem_avg_w   = mem_area   / h if h > 0 else 0
        stiff_avg_w = stiff_area / h if h > 0 else 0
        mem_mass    = mem_area   * MEM_THICKNESS   * MEM_DENSITY
        stiff_mass  = stiff_area * STIFF_THICKNESS * STIFF_DENSITY
        mem_EA      = E_MEMBRANE * MEM_THICKNESS   * mem_avg_w
        stiff_EA    = E_CARBON   * STIFF_THICKNESS * stiff_avg_w
        total_mass  = mem_mass  + stiff_mass
        total_EA    = mem_EA    + stiff_EA
        total_EAy   = mem_EA*mem_cg_y + stiff_EA*stiff_cg_y
        cg_y = (mem_mass*mem_cg_y + stiff_mass*stiff_cg_y) / total_mass if total_mass > 0 else cy_mid
        return total_mass, cx_mid, cg_y, total_EA, total_EAy

    elif region_type == 'stiff_tri_left':
        xr_top     = x_right(y_top)
        xr_bot     = x_right(y_b)
        sw_top     = max(0, min(xr_top, x_e) - x_s)
        sw_bot     = max(0, min(xr_bot, x_e) - x_s)
        mw_top     = w - sw_top
        mw_bot     = w - sw_bot
        stiff_area = 0.5 * (sw_top + sw_bot) * h
        mem_area   = w * h - stiff_area

        if (sw_top + sw_bot) > 1e-10:
            stiff_cg_y = y_b + h*(2*sw_top + sw_bot) / (3*(sw_top + sw_bot))
        else:
            stiff_cg_y = cy_mid

        if (mw_top + mw_bot) > 1e-10:
            mem_cg_y   = y_b + h*(2*mw_top + mw_bot) / (3*(mw_top + mw_bot))
        else:
            mem_cg_y   = cy_mid

        stiff_avg_w = stiff_area / h if h > 0 else 0
        mem_avg_w   = mem_area   / h if h > 0 else 0
        stiff_mass  = stiff_area * STIFF_THICKNESS * STIFF_DENSITY
        mem_mass    = mem_area   * MEM_THICKNESS   * MEM_DENSITY
        stiff_EA    = E_CARBON   * STIFF_THICKNESS * stiff_avg_w
        mem_EA      = E_MEMBRANE * MEM_THICKNESS   * mem_avg_w
        total_mass  = mem_mass  + stiff_mass
        total_EA    = mem_EA    + stiff_EA
        total_EAy   = mem_EA*mem_cg_y + stiff_EA*stiff_cg_y
        cg_y = (mem_mass*mem_cg_y + stiff_mass*stiff_cg_y) / total_mass if total_mass > 0 else cy_mid
        return total_mass, cx_mid, cg_y, total_EA, total_EAy

    return 0, cx_mid, cy_mid, 0, 0

# TYPE B SUB-STRIP CALCULATOR
def calc_type_b_substrip(x_s, x_e, global_x_offset=0.0):
    w = x_e - x_s
    yt_xl_xs = y_transition(x_s, 'left')
    yt_xl_xe = y_transition(x_e, 'left')
    yt_xr_xs = y_transition(x_s, 'right')
    yt_xr_xe = y_transition(x_e, 'right')
    y_bounds = sorted(set([0.0, MEM_HEIGHT,
                            yt_xl_xs, yt_xl_xe,
                            yt_xr_xs, yt_xr_xe]))
    y_bounds = [y for y in y_bounds if 0.0 <= y <= MEM_HEIGHT]
    y_bounds = sorted(set(y_bounds))

    total_mass = 0.0
    total_EA   = 0.0
    total_EAy  = 0.0
    sum_mcgx   = 0.0
    sum_mcgy   = 0.0

    for j in range(len(y_bounds) - 1):
        y_b    = y_bounds[j]
        y_top  = y_bounds[j + 1]
        y_mid  = (y_b + y_top) / 2.0
        height = y_top - y_b

        if height < 1e-10:
            continue

        xl_mid = x_left(y_mid)
        xr_mid = x_right(y_mid)

        in_left_mem   = xl_mid > x_e
        in_right_mem  = xr_mid < x_s
        in_stiff_full = xl_mid <= x_s and xr_mid >= x_e

        if in_left_mem:
            region_type = 'mem_rect'
        elif in_right_mem:
            region_type = 'right_mem_rect'
        elif in_stiff_full:
            region_type = 'stiff_rect'
        elif xl_mid > x_s and xl_mid <= x_e and xr_mid >= x_e:
            region_type = 'mem_tri_left'
        elif xr_mid >= x_s and xr_mid < x_e and xl_mid <= x_s:
            region_type = 'stiff_tri_left'
        else:
            region_type = 'mem_rect'

        mass, cg_x, cg_y, EA, EAy = calc_region(
            x_s, x_e, y_b, y_top, region_type, global_x_offset)

        total_mass += mass
        total_EA   += EA
        total_EAy  += EAy
        sum_mcgx   += mass * cg_x
        sum_mcgy   += mass * cg_y

    # Spar contribution
    spar_mass = w * SPAR_HEIGHT * SPAR_THICKNESS * SPAR_DENSITY
    spar_cg_x = global_x_offset + x_s + w / 2.0
    spar_cg_y = SPAR_BOTTOM + SPAR_HEIGHT / 2.0
    spar_EA   = E_CARBON * w * SPAR_THICKNESS
    spar_EAy  = spar_EA * spar_cg_y

    total_mass += spar_mass
    total_EA   += spar_EA
    total_EAy  += spar_EAy
    sum_mcgx   += spar_mass * spar_cg_x
    sum_mcgy   += spar_mass * spar_cg_y

    cg_x_final = sum_mcgx / total_mass
    cg_y_final = sum_mcgy / total_mass
    na_y_final = total_EAy / total_EA

    return cg_x_final, cg_y_final, na_y_final, total_mass

# TYPE A STRIP CALCULATOR
def calc_type_a_strip(x_start):
    w          = TYPE_A_WIDTH
    spar_mass  = w * SPAR_HEIGHT * SPAR_THICKNESS * SPAR_DENSITY
    spar_cg_x  = x_start + w / 2.0
    spar_cg_y  = SPAR_BOTTOM + SPAR_HEIGHT / 2.0
    spar_EA    = E_CARBON * w * SPAR_THICKNESS
    spar_EAy   = spar_EA * spar_cg_y

    mem_mass   = w * MEM_HEIGHT * MEM_THICKNESS * MEM_DENSITY
    mem_cg_x   = x_start + w / 2.0
    mem_cg_y   = MEM_HEIGHT / 2.0
    mem_EA     = E_MEMBRANE * w * MEM_THICKNESS
    mem_EAy    = mem_EA * mem_cg_y

    total_mass = spar_mass + mem_mass
    total_EA   = spar_EA + mem_EA
    cg_x = (spar_mass * spar_cg_x + mem_mass * mem_cg_x) / total_mass
    cg_y = (spar_mass * spar_cg_y + mem_mass * mem_cg_y) / total_mass
    na_y = (spar_EAy + mem_EAy) / total_EA

    return cg_x, cg_y, na_y, total_mass

# FULL WING CALCULATION
all_cg_x   = []
all_cg_y   = []
all_na_y   = []
all_masses = []
all_types  = []

sub_width  = TYPE_B_WIDTH / NUM_SUBSTRIPS
x_current  = 0.0

print(f"\n{'Strip':<10} {'CG_x':>10} {'CG_y':>10} {'NA_y':>10} {'Mass':>10}")
print("-" * 55)

for i in range(NUM_PAIRS):

    # ── TYPE A ──
    cg_x, cg_y, na_y, mass = calc_type_a_strip(x_current)
    all_cg_x.append(cg_x)
    all_cg_y.append(cg_y)
    all_na_y.append(na_y)
    all_masses.append(mass)
    all_types.append('A')
    print(f"  A{i+1:<7}  {cg_x:>10.4f} {cg_y:>10.4f} {na_y:>10.4f} {mass:>10.4f}")
    x_current += TYPE_A_WIDTH

    # ── TYPE B — 10 sub-strips ──
    for k in range(NUM_SUBSTRIPS):
        x_s = k * sub_width
        x_e = x_s + sub_width
        cg_x, cg_y, na_y, mass = calc_type_b_substrip(x_s, x_e,
                                                        global_x_offset=x_current)
        all_cg_x.append(cg_x)
        all_cg_y.append(cg_y)
        all_na_y.append(na_y)
        all_masses.append(mass)
        all_types.append('B')
        print(f"  B{i+1}-S{k+1:<4}  {cg_x:>10.4f} {cg_y:>10.4f} {na_y:>10.4f} {mass:>10.4f}")

    x_current += TYPE_B_WIDTH

# Overall CG
total_mass  = sum(all_masses)
overall_cgx = sum(m * x for m, x in zip(all_masses, all_cg_x)) / total_mass
overall_cgy = sum(m * y for m, y in zip(all_masses, all_cg_y)) / total_mass

print("-" * 55)
print(f"\n  Total Mass   : {total_mass:.4f} mg")
print(f"  Overall CG   : x = {overall_cgx:.4f} mm, y = {overall_cgy:.4f} mm")

# HELPER: smooth spline
def smooth_spline(x_pts, y_pts, k=3):
    x_arr    = np.array(x_pts)
    y_arr    = np.array(y_pts)
    x_smooth = np.linspace(x_arr.min(), x_arr.max(), 500)
    spline   = make_interp_spline(x_arr, y_arr, k=k,
                                   bc_type=([(1, 0.0)], [(1, 0.0)]))
    return x_smooth, spline(x_smooth)

# PLOTTING
fig, ax = plt.subplots(1, 1, figsize=(16, 7))
fig.patch.set_facecolor('#0f1117')
ax.set_facecolor('#1a1d27')

# ── Draw wing layout ──
x_current = 0.0
for i in range(NUM_PAIRS):

    # Type A
    ax.add_patch(patches.Rectangle((x_current, 0), TYPE_A_WIDTH, MEM_HEIGHT,
                                    linewidth=0.5, edgecolor='white',
                                    facecolor='#2a9d8f', alpha=0.4))
    ax.add_patch(patches.Rectangle((x_current, SPAR_BOTTOM), TYPE_A_WIDTH, SPAR_HEIGHT,
                                    linewidth=0.5, edgecolor='white',
                                    facecolor='#e9c46a', alpha=0.7))
    x_current += TYPE_A_WIDTH

    xs = x_current

    # Left triangle (membrane)
    ax.add_patch(Polygon([(xs, 0), (xs, MEM_HEIGHT),
                           (xs + STIFF_X_BOTTOM, 0)],
                          closed=True, linewidth=0.5, edgecolor='white',
                          facecolor='#2a9d8f', alpha=0.4))

    # Inclined stiffener (parallelogram)
    ax.add_patch(Polygon([(xs, MEM_HEIGHT), (xs + STIFF_TOP_X, MEM_HEIGHT),
                           (xs + TYPE_B_WIDTH, 0), (xs + STIFF_X_BOTTOM, 0)],
                          closed=True, linewidth=0.5, edgecolor='white',
                          facecolor='#e9c46a', alpha=0.7))

    # Right triangle (membrane)
    ax.add_patch(Polygon([(xs + STIFF_TOP_X, MEM_HEIGHT),
                           (xs + TYPE_B_WIDTH, MEM_HEIGHT),
                           (xs + TYPE_B_WIDTH, 0)],
                          closed=True, linewidth=0.5, edgecolor='white',
                          facecolor='#2a9d8f', alpha=0.4))

    # Leading edge spar
    ax.add_patch(patches.Rectangle((xs, SPAR_BOTTOM), TYPE_B_WIDTH, SPAR_HEIGHT,
                                    linewidth=0.5, edgecolor='white',
                                    facecolor='#e9c46a', alpha=0.7))

    # Sub-strip boundaries for Type B
    for k in range(NUM_SUBSTRIPS + 1):
        x = xs + k * sub_width
        ax.axvline(x=x, color='white', linewidth=0.3, alpha=0.3, linestyle=':')

    x_current += TYPE_B_WIDTH

# ── Plot CG points ──
for idx, (cx, cy, stype) in enumerate(zip(all_cg_x, all_cg_y, all_types)):
    color = '#ff9f1c' if stype == 'A' else '#e76f51'
    ax.scatter(cx, cy, color=color, s=60, zorder=5,
               edgecolors='white', linewidths=0.6)

# ── CG smooth curve ──
x_cg_smooth, y_cg_smooth = smooth_spline(all_cg_x, all_cg_y)
ax.plot(x_cg_smooth, y_cg_smooth, color='#e76f51', linewidth=2.0,
        linestyle='--', alpha=0.9, label='CG Line')

# ── Plot NA points ──
for idx, (cx, ny, stype) in enumerate(zip(all_cg_x, all_na_y, all_types)):
    color = '#00b4d8' if stype == 'A' else '#90e0ef'
    ax.scatter(cx, ny, color=color, s=60, zorder=5,
               marker='D', edgecolors='white', linewidths=0.6)

# ── NA smooth curve ──
x_na_smooth, y_na_smooth = smooth_spline(all_cg_x, all_na_y)
ax.plot(x_na_smooth, y_na_smooth, color='#00b4d8', linewidth=2.0,
        linestyle='-', alpha=0.9, label='Neutral Axis')

# ── Overall CG ──
ax.scatter([overall_cgx], [overall_cgy], color='#ff4d6d', s=200,
           marker='*', zorder=6, edgecolors='white', linewidths=1)

# ── Legend ──
legend_elements = [
    Patch(facecolor='#2a9d8f', alpha=0.5,  label='Membrane (Polyimide)'),
    Patch(facecolor='#e9c46a', alpha=0.8,  label='Spar & Stiffener (Carbon Fibre)'),
    Line2D([0],[0], color='#e76f51', linewidth=2, linestyle='--', label='CG Line'),
    Line2D([0],[0], color='#00b4d8', linewidth=2, linestyle='-',  label='Neutral Axis'),
    Line2D([0],[0], marker='o', color='w', markerfacecolor='#ff9f1c',
           markersize=8, label='Strip CG (Type A)'),
    Line2D([0],[0], marker='o', color='w', markerfacecolor='#e76f51',
           markersize=8, label='Sub-strip CG (Type B)'),
    Line2D([0],[0], marker='D', color='w', markerfacecolor='#00b4d8',
           markersize=8, label='Strip NA (Type A)'),
    Line2D([0],[0], marker='D', color='w', markerfacecolor='#90e0ef',
           markersize=8, label='Sub-strip NA (Type B)'),
    Line2D([0],[0], marker='*', color='w', markerfacecolor='#ff4d6d',
           markersize=12, label=f'Overall CG ({overall_cgx:.2f}, {overall_cgy:.2f})'),
]

ax.set_xlim(-0.5, WINGSPAN + 0.5)
ax.set_ylim(-0.5, CHORD + 1.0)
ax.set_xlabel('x — Wingspan (mm)', color='white', fontsize=11)
ax.set_ylabel('y — Chord (mm)', color='white', fontsize=11)
ax.set_title(f'Full Wing — CG and Neutral Axis (θ = {np.degrees(THETA):.0f}°)',
             color='white', fontsize=12, fontweight='bold')
ax.tick_params(colors='white')
ax.spines[:].set_color('#444')
ax.legend(handles=legend_elements, loc='upper right', fontsize=8,
          facecolor='#1a1d27', edgecolor='#444', labelcolor='white', ncol=2)

plt.suptitle(f'Rectangular Wing with Inclined Stiffeners (θ = {np.degrees(THETA):.0f}°)\n'
             f'CG Line vs Neutral Axis — 3 Pairs (Type A + Type B with 10 Sub-strips)',
             color='white', fontsize=13, fontweight='bold')

plt.tight_layout()
plt.show()
