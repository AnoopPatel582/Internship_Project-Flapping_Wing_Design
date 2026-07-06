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
THETA           = np.radians(35.0)
STIFF_WIDTH     = 1.0    # mm (perpendicular width)
NUM_SUBSTRIPS   = 10

# Derived geometry
STIFF_X_BOTTOM  = MEM_HEIGHT / np.tan(THETA)
STIFF_TOP_X     = STIFF_WIDTH / np.sin(THETA)
TYPE_B_WIDTH    = STIFF_X_BOTTOM + STIFF_TOP_X

# Material thicknesses
SPAR_THICKNESS  = 0.5    # mm
MEM_THICKNESS   = 0.05   # mm
STIFF_THICKNESS = 0.5    # mm

# Material densities
SPAR_DENSITY    = 1.6    # mg/mm³
MEM_DENSITY     = 1.3    # mg/mm³
STIFF_DENSITY   = 1.6    # mg/mm³

# Young's Modulus
E_CARBON        = 80.0   # GPa
E_MEMBRANE      = 5.0    # GPa

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
    h       = y_top - y_b
    w       = x_e - x_s
    cx_mid  = global_x_offset + x_s + w / 2.0
    cy_mid  = y_b + h / 2.0

    if region_type == 'mem_rect':
        area    = w * h
        mass    = area * MEM_THICKNESS * MEM_DENSITY
        EA      = E_MEMBRANE * MEM_THICKNESS * w
        EAy     = EA * cy_mid
        cg_x    = cx_mid
        cg_y    = cy_mid

    elif region_type == 'stiff_rect':
        area    = w * h
        mass    = area * STIFF_THICKNESS * STIFF_DENSITY
        EA      = E_CARBON * STIFF_THICKNESS * w
        EAy     = EA * cy_mid
        cg_x    = cx_mid
        cg_y    = cy_mid

    elif region_type == 'mem_tri_left':
        # Left membrane triangle + right stiffener triangle
        # Diagonal goes from (x_left(y_top), y_top) to (x_left(y_b), y_b)
        xl_top  = x_left(y_top)   # x_left at top of region
        xl_bot  = x_left(y_b)     # x_left at bottom of region

        # Membrane triangle (left of diagonal)
        # Vertices: (x_s, y_b), (x_s, y_top), (xl_top, y_top) approximately
        mem_base   = xl_top - x_s   # width at top
        mem_area   = 0.5 * (xl_bot - x_s + xl_top - x_s) * h  # trapezoid approx
        mem_cg_x   = global_x_offset + x_s + (xl_bot - x_s + xl_top - x_s) / 3.0
        mem_cg_y   = y_b + h * (2*(xl_top-x_s) + (xl_bot-x_s)) / (3*((xl_top-x_s) + (xl_bot-x_s))) if (xl_top-x_s + xl_bot-x_s) > 0 else cy_mid
        mem_avg_w  = mem_area / h if h > 0 else 0
        mem_mass   = mem_area * MEM_THICKNESS * MEM_DENSITY
        mem_EA     = E_MEMBRANE * MEM_THICKNESS * mem_avg_w
        mem_EAy    = mem_EA * mem_cg_y

        # Stiffener triangle (right of diagonal)
        stiff_area = w * h - mem_area
        stiff_cg_x = global_x_offset + x_s + (xl_bot - x_s + xl_top - x_s) / 3.0 + (w - (xl_bot-x_s+xl_top-x_s)/2) * 0.5
        stiff_cg_y = y_b + h * ((xl_top-x_s) + 2*(xl_bot-x_s)) / (3*((xl_top-x_s) + (xl_bot-x_s))) if (xl_top-x_s + xl_bot-x_s) > 0 else cy_mid
        stiff_avg_w = stiff_area / h if h > 0 else 0
        stiff_mass  = stiff_area * STIFF_THICKNESS * STIFF_DENSITY
        stiff_EA    = E_CARBON * STIFF_THICKNESS * stiff_avg_w
        stiff_EAy   = stiff_EA * stiff_cg_y

        total_mass  = mem_mass + stiff_mass
        total_EA    = mem_EA + stiff_EA
        total_EAy   = mem_EAy + stiff_EAy
        cg_x        = (mem_mass * mem_cg_x + stiff_mass * stiff_cg_x) / total_mass if total_mass > 0 else cx_mid
        cg_y        = (mem_mass * mem_cg_y + stiff_mass * stiff_cg_y) / total_mass if total_mass > 0 else cy_mid

        return total_mass, cg_x, cg_y, total_EA, total_EAy

    elif region_type == 'stiff_tri_left':
        # Left stiffener triangle + right membrane triangle
        # Diagonal is x_right(y)
        xr_top  = x_right(y_top)
        xr_bot  = x_right(y_b)

        # Stiffener triangle (left of x_right diagonal)
        stiff_w_top  = xr_top - x_s
        stiff_w_bot  = xr_bot - x_s
        stiff_area   = 0.5 * (stiff_w_top + stiff_w_bot) * h
        stiff_cg_y   = y_b + h * (2*stiff_w_top + stiff_w_bot) / (3*(stiff_w_top + stiff_w_bot)) if (stiff_w_top + stiff_w_bot) > 0 else cy_mid
        stiff_avg_w  = stiff_area / h if h > 0 else 0
        stiff_mass   = stiff_area * STIFF_THICKNESS * STIFF_DENSITY
        stiff_EA     = E_CARBON * STIFF_THICKNESS * stiff_avg_w
        stiff_EAy    = stiff_EA * stiff_cg_y
        stiff_cg_x   = global_x_offset + x_s + (stiff_w_top + stiff_w_bot) / 3.0

        # Membrane triangle (right of x_right diagonal)
        mem_area    = w * h - stiff_area
        mem_cg_y    = y_b + h * (stiff_w_top + 2*stiff_w_bot) / (3*(stiff_w_top + stiff_w_bot)) if (stiff_w_top + stiff_w_bot) > 0 else cy_mid
        mem_avg_w   = mem_area / h if h > 0 else 0
        mem_mass    = mem_area * MEM_THICKNESS * MEM_DENSITY
        mem_EA      = E_MEMBRANE * MEM_THICKNESS * mem_avg_w
        mem_EAy     = mem_EA * mem_cg_y
        mem_cg_x    = global_x_offset + x_s + (stiff_w_top + stiff_w_bot) / 3.0 + (w - (stiff_w_top+stiff_w_bot)/2) * 0.5

        total_mass  = mem_mass + stiff_mass
        total_EA    = mem_EA + stiff_EA
        total_EAy   = mem_EAy + stiff_EAy
        cg_x        = (mem_mass * mem_cg_x + stiff_mass * stiff_cg_x) / total_mass if total_mass > 0 else cx_mid
        cg_y        = (mem_mass * mem_cg_y + stiff_mass * stiff_cg_y) / total_mass if total_mass > 0 else cy_mid

        return total_mass, cg_x, cg_y, total_EA, total_EAy

    # Common return for rect types
    return mass, cg_x, cg_y, EA, EAy

# SUB-STRIP CALCULATOR FOR TYPE B
def calc_type_b_substrip(x_s, x_e, global_x_offset=0.0):
    h        = MEM_HEIGHT
    w        = x_e - x_s

    # Transition points (clamped to [0, 6])
    yt_xl_xs = y_transition(x_s, 'left')    # x_left = x_s
    yt_xl_xe = y_transition(x_e, 'left')    # x_left = x_e
    yt_xr_xs = y_transition(x_s, 'right')   # x_right = x_s
    yt_xr_xe = y_transition(x_e, 'right')   # x_right = x_e

    # Collect all unique y boundaries and sort descending
    y_bounds = sorted(set([0.0, MEM_HEIGHT,
                            yt_xl_xs, yt_xl_xe,
                            yt_xr_xs, yt_xr_xe]), reverse=False)
    y_bounds = [y for y in y_bounds if 0.0 <= y <= MEM_HEIGHT]
    y_bounds = sorted(set(y_bounds))

    total_mass = 0.0
    total_EA   = 0.0
    total_EAy  = 0.0
    sum_mcgx   = 0.0
    sum_mcgy   = 0.0

    for j in range(len(y_bounds) - 1):
        y_b   = y_bounds[j]
        y_top = y_bounds[j + 1]
        y_mid = (y_b + y_top) / 2.0
        height = y_top - y_b

        if height < 1e-10:
            continue

        # Determine region type at y_mid
        xl_mid = x_left(y_mid)
        xr_mid = x_right(y_mid)

        # Check what's in this sub-strip at y_mid
        in_left_mem   = xl_mid > x_e          # full left membrane
        in_right_mem  = xr_mid < x_s          # full right membrane
        in_stiff_full = xl_mid <= x_s and xr_mid >= x_e  # full stiffener

        if in_left_mem:
            region_type = 'mem_rect'
        elif in_right_mem:
            region_type = 'mem_rect'
        elif in_stiff_full:
            region_type = 'stiff_rect'
        elif xl_mid > x_s and xl_mid < x_e and xr_mid >= x_e:
            # x_left crosses sub-strip: left mem + right stiffener
            region_type = 'mem_tri_left'
        elif xr_mid > x_s and xr_mid < x_e and xl_mid <= x_s:
            # x_right crosses sub-strip: left stiffener + right mem
            region_type = 'stiff_tri_left'
        else:
            region_type = 'mem_rect'   # fallback

        mass, cg_x, cg_y, EA, EAy = calc_region(
            x_s, x_e, y_b, y_top, region_type, global_x_offset)

        total_mass += mass
        total_EA   += EA
        total_EAy  += EAy
        sum_mcgx   += mass * cg_x
        sum_mcgy   += mass * cg_y

    # Add spar contribution
    spar_mass  = w * SPAR_HEIGHT * SPAR_THICKNESS * SPAR_DENSITY
    spar_cg_x  = global_x_offset + x_s + w / 2.0
    spar_cg_y  = SPAR_BOTTOM + SPAR_HEIGHT / 2.0
    spar_A     = w * SPAR_THICKNESS
    spar_EA    = E_CARBON * spar_A
    spar_EAy   = spar_EA * spar_cg_y

    total_mass += spar_mass
    total_EA   += spar_EA
    total_EAy  += spar_EAy
    sum_mcgx   += spar_mass * spar_cg_x
    sum_mcgy   += spar_mass * spar_cg_y

    cg_x_final = sum_mcgx / total_mass
    cg_y_final = sum_mcgy / total_mass
    na_y_final = total_EAy / total_EA

    return cg_x_final, cg_y_final, na_y_final, total_mass

# CALCULATE 10 SUB-STRIPS OF TYPE B
sub_width   = TYPE_B_WIDTH / NUM_SUBSTRIPS
sub_cg_x    = []
sub_cg_y    = []
sub_na_y    = []
sub_masses  = []

print(f"\n{'Sub':<5} {'x_s':>8} {'x_e':>8} {'CG_x':>10} {'CG_y':>10} {'NA_y':>10}")
print("-" * 55)

for k in range(NUM_SUBSTRIPS):
    x_s = k * sub_width
    x_e = x_s + sub_width
    cg_x, cg_y, na_y, mass = calc_type_b_substrip(x_s, x_e, global_x_offset=0.0)
    sub_cg_x.append(cg_x + x_s)   # global x position
    sub_cg_y.append(cg_y)
    sub_na_y.append(na_y)
    sub_masses.append(mass)
    print(f"  {k+1:<3} {x_s:>8.4f} {x_e:>8.4f} {cg_x:>10.4f} {cg_y:>10.4f} {na_y:>10.4f}")

print("=" * 65)

# PLOTTING — TYPE B STRIP ONLY
fig, ax = plt.subplots(1, 1, figsize=(12, 7))
fig.patch.set_facecolor('#0f1117')
ax.set_facecolor('#1a1d27')

# Draw Type B strip
# Left triangle (membrane)
ax.add_patch(Polygon([(0, 0), (0, MEM_HEIGHT), (STIFF_X_BOTTOM, 0)],
                      closed=True, linewidth=1.0, edgecolor='white',
                      facecolor='#2a9d8f', alpha=0.5))

# Inclined stiffener (parallelogram)
ax.add_patch(Polygon([(0, MEM_HEIGHT), (STIFF_TOP_X, MEM_HEIGHT),
                       (TYPE_B_WIDTH, 0), (STIFF_X_BOTTOM, 0)],
                      closed=True, linewidth=1.0, edgecolor='white',
                      facecolor='#e9c46a', alpha=0.7))

# Right triangle (membrane)
ax.add_patch(Polygon([(STIFF_TOP_X, MEM_HEIGHT), (TYPE_B_WIDTH, MEM_HEIGHT),
                       (TYPE_B_WIDTH, 0)],
                      closed=True, linewidth=1.0, edgecolor='white',
                      facecolor='#2a9d8f', alpha=0.5))

# Leading edge spar
ax.add_patch(patches.Rectangle((0, SPAR_BOTTOM), TYPE_B_WIDTH, SPAR_HEIGHT,
                                 linewidth=1.0, edgecolor='white',
                                 facecolor='#e9c46a', alpha=0.9))

# Sub-strip boundaries
for k in range(NUM_SUBSTRIPS + 1):
    x = k * sub_width
    ax.axvline(x=x, color='white', linewidth=0.5, alpha=0.4, linestyle=':')

# Plot CG points
for k, (cx, cy) in enumerate(zip(sub_cg_x, sub_cg_y)):
    ax.scatter(cx, cy, color='#e76f51', s=80, zorder=5,
               edgecolors='white', linewidths=0.8)
    ax.annotate(f'{k+1}', (cx, cy), textcoords="offset points",
                xytext=(0, 8), ha='center', fontsize=7,
                color='white', fontweight='bold')

# CG smooth curve
x_arr    = np.array(sub_cg_x)
y_arr    = np.array(sub_cg_y)
x_smooth = np.linspace(x_arr.min(), x_arr.max(), 500)
spline   = make_interp_spline(x_arr, y_arr, k=3,
                               bc_type=([(1, 0.0)], [(1, 0.0)]))
ax.plot(x_smooth, spline(x_smooth), color='#e76f51', linewidth=2.0,
        linestyle='--', alpha=0.9, label='CG Line')

# Plot NA points
for k, (cx, ny) in enumerate(zip(sub_cg_x, sub_na_y)):
    ax.scatter(cx, ny, color='#00b4d8', s=80, zorder=5,
               marker='D', edgecolors='white', linewidths=0.8)

# NA smooth curve
y_na_arr = np.array(sub_na_y)
spline_na = make_interp_spline(x_arr, y_na_arr, k=3,
                                bc_type=([(1, 0.0)], [(1, 0.0)]))
ax.plot(x_smooth, spline_na(x_smooth), color='#00b4d8', linewidth=2.0,
        linestyle='-', alpha=0.9, label='Neutral Axis')

# Legend
legend_elements = [
    Patch(facecolor='#2a9d8f', alpha=0.5, label='Membrane (Polyimide)'),
    Patch(facecolor='#e9c46a', alpha=0.8, label='Spar & Stiffener (Carbon Fibre)'),
    Line2D([0],[0], color='#e76f51', linewidth=2, linestyle='--', label='CG Line'),
    Line2D([0],[0], color='#00b4d8', linewidth=2, linestyle='-',  label='Neutral Axis'),
    Line2D([0],[0], marker='o', color='w', markerfacecolor='#e76f51',
           markersize=8, label='Sub-strip CG'),
    Line2D([0],[0], marker='D', color='w', markerfacecolor='#00b4d8',
           markersize=8, label='Sub-strip NA'),
]

ax.set_xlim(-0.3, TYPE_B_WIDTH + 0.3)
ax.set_ylim(-0.5, 8.0)
ax.set_xlabel('x — Wingspan (mm)', color='white', fontsize=11)
ax.set_ylabel('y — Chord (mm)', color='white', fontsize=11)
ax.set_title(f'Type B Strip — CG and Neutral Axis (θ = {np.degrees(THETA):.0f}°, 10 Sub-strips)',
             color='white', fontsize=12, fontweight='bold')
ax.tick_params(colors='white')
ax.spines[:].set_color('#444')
ax.legend(handles=legend_elements, loc='upper right', fontsize=9,
          facecolor='#1a1d27', edgecolor='#444', labelcolor='white')

plt.suptitle('Inclined Stiffener Wing — Type B Strip Analysis',
             color='white', fontsize=13, fontweight='bold')

plt.tight_layout()
plt.show()