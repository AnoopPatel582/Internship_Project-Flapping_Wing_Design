import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import Polygon, Patch
from matplotlib.lines import Line2D
from scipy.interpolate import make_interp_spline

#varying parameters
THETA_DEG    = 45.0   # degrees (0 to 90)
SPAR_HEIGHT  = 1.0    # mm (spar height in y-direction)
STIFF_WIDTH  = 1.0    # mm (perpendicular width of stiffener)

NUM_A_SUBSTRIPS = 5
NUM_B_SUBSTRIPS = 10  

#fixed parameters
TOTAL_SPAN      = 50.0
CHORD           = 7.0
STIFF_STARTS    = [5.0, 20.0, 35.0]

SPAR_THICKNESS  = 0.5
MEM_THICKNESS   = 0.05
STIFF_THICKNESS = 0.5

SPAR_DENSITY    = 1.6
MEM_DENSITY     = 1.3
STIFF_DENSITY   = 1.6

E_CARBON        = 80.0
E_MEMBRANE      = 5.0
G_CARBON        = 27.0
G_MEMBRANE      = 1.5

#Auto driven geometry
MEM_HEIGHT   = CHORD - SPAR_HEIGHT
SPAR_BOTTOM  = MEM_HEIGHT
THETA        = np.radians(THETA_DEG)
NO_STIFFENER = (THETA_DEG == 0.0)

if NO_STIFFENER:
    STIFF_X_BOTTOM = 0.0
    STIFF_TOP_X    = 0.0
    TYPE_B_WIDTH   = 0.0
elif THETA_DEG == 90.0:
    STIFF_X_BOTTOM = 0.0
    STIFF_TOP_X    = STIFF_WIDTH / np.sin(THETA)
    TYPE_B_WIDTH   = STIFF_TOP_X
else:
    STIFF_X_BOTTOM = MEM_HEIGHT / np.tan(THETA)
    STIFF_TOP_X    = STIFF_WIDTH / np.sin(THETA)
    TYPE_B_WIDTH   = STIFF_X_BOTTOM + STIFF_TOP_X

if NO_STIFFENER:
    strips = [('A', 0.0, TOTAL_SPAN, NUM_A_SUBSTRIPS)]
else:
    strips = []
    x_prev = 0.0
    for x_start in STIFF_STARTS:
        if x_start > x_prev:
            strips.append(('A', x_prev, x_start, NUM_A_SUBSTRIPS))
        strips.append(('B', x_start, x_start + TYPE_B_WIDTH, NUM_B_SUBSTRIPS))
        x_prev = x_start + TYPE_B_WIDTH
    if x_prev < TOTAL_SPAN:
        strips.append(('A', x_prev, TOTAL_SPAN, NUM_A_SUBSTRIPS))

print("  WING PARAMETERS")
print(f"  THETA            : {THETA_DEG:.1f} degrees")
print(f"  SPAR_HEIGHT      : {SPAR_HEIGHT:.4f} mm")
print(f"  STIFF_WIDTH      : {STIFF_WIDTH:.4f} mm")
print(f"  MEM_HEIGHT       : {MEM_HEIGHT:.4f} mm")
print(f"  STIFF_X_BOTTOM   : {STIFF_X_BOTTOM:.4f} mm")
print(f"  STIFF_TOP_X      : {STIFF_TOP_X:.4f} mm")
print(f"  TYPE_B_WIDTH     : {TYPE_B_WIDTH:.4f} mm")
print(f"\n  Strip layout:")
for s in strips:
    print(f"    Type {s[0]}: x={s[1]:.3f} to x={s[2]:.3f} "
          f"(width={s[2]-s[1]:.3f}mm, {s[3]} sub-strips)")

# HELPER FUNCTIONS
def x_left(y):
    if STIFF_X_BOTTOM == 0:
        return 0.0
    return STIFF_X_BOTTOM * (1 - y / MEM_HEIGHT)

def x_right(y):
    return STIFF_TOP_X + STIFF_X_BOTTOM * (1 - y / MEM_HEIGHT)

def clamp(val, lo=0.0, hi=MEM_HEIGHT):
    return max(lo, min(hi, val))

def y_transition(x_val, boundary='left'):
    if boundary == 'left':
        if STIFF_X_BOTTOM == 0:
            return MEM_HEIGHT
        y = MEM_HEIGHT * (1 - x_val / STIFF_X_BOTTOM)
    else:
        if STIFF_X_BOTTOM == 0:
            y = MEM_HEIGHT * (1 - (x_val - STIFF_TOP_X) / 1e-10)
        else:
            y = MEM_HEIGHT * (1 - (x_val - STIFF_TOP_X) / STIFF_X_BOTTOM)
    return clamp(y)

# SECTION PROPERTY FUNCTIONS
def calc_I(width, thickness):
    """Base formulation for I (will be integrated/multiplied by h later)"""
    if width < 1e-10:
        return 0.0
    return (1.0 / 12.0) * width * thickness**3

def calc_J(width, thickness):
    """Continuous integration torsional constant: J = w*t^3/3 (No finite-edge correction)"""
    if width < 1e-10:
        return 0.0
    return (width * thickness**3 / 3.0)

# REGION CALCULATOR
def calc_region(x_s, x_e, y_b, y_top, region_type, global_x_offset=0.0):
    h      = y_top - y_b
    w      = x_e - x_s
    cx_mid = global_x_offset + x_s + w / 2.0
    cy_mid = y_b + h / 2.0

    if region_type in ('mem_rect', 'right_mem_rect'):
        mass = w * h * MEM_THICKNESS * MEM_DENSITY
        I    = calc_I(w, MEM_THICKNESS) * h
        EI   = E_MEMBRANE * I
        EIy  = EI * cy_mid
        J    = calc_J(w, MEM_THICKNESS) * h
        GJ   = G_MEMBRANE * J
        GJy  = GJ * cy_mid
        return mass, cx_mid, cy_mid, EI, EIy, GJ, GJy

    elif region_type == 'stiff_rect':
        mass = w * h * STIFF_THICKNESS * STIFF_DENSITY
        I    = calc_I(w, STIFF_THICKNESS) * h
        EI   = E_CARBON * I
        EIy  = EI * cy_mid
        J    = calc_J(w, STIFF_THICKNESS) * h
        GJ   = G_CARBON * J
        GJy  = GJ * cy_mid
        return mass, cx_mid, cy_mid, EI, EIy, GJ, GJy

    elif region_type == 'mem_tri_left':
        xl_top     = x_left(y_top)
        xl_bot     = x_left(y_b)
        mw_top     = max(0, xl_top - x_s)
        mw_bot     = max(0, xl_bot - x_s)
        sw_top     = w - mw_top
        sw_bot     = w - mw_bot
        mem_area   = 0.5 * (mw_top + mw_bot) * h
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
        
        # Multiply by height h for proper calculus integration weighting
        mem_I       = calc_I(mem_avg_w,   MEM_THICKNESS) * h
        stiff_I     = calc_I(stiff_avg_w, STIFF_THICKNESS) * h
        mem_EI      = E_MEMBRANE * mem_I
        stiff_EI    = E_CARBON   * stiff_I
        
        mem_J       = calc_J(mem_avg_w,   MEM_THICKNESS) * h
        stiff_J     = calc_J(stiff_avg_w, STIFF_THICKNESS) * h
        mem_GJ      = G_MEMBRANE * mem_J
        stiff_GJ    = G_CARBON   * stiff_J

        total_mass  = mem_mass  + stiff_mass
        total_EI    = mem_EI    + stiff_EI
        total_EIy   = mem_EI*mem_cg_y   + stiff_EI*stiff_cg_y
        total_GJ    = mem_GJ    + stiff_GJ
        total_GJy   = mem_GJ*mem_cg_y   + stiff_GJ*stiff_cg_y
        cg_y = (mem_mass*mem_cg_y + stiff_mass*stiff_cg_y) / total_mass if total_mass > 0 else cy_mid
        return total_mass, cx_mid, cg_y, total_EI, total_EIy, total_GJ, total_GJy

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
        
        # Multiply by height h for proper calculus integration weighting
        stiff_I     = calc_I(stiff_avg_w, STIFF_THICKNESS) * h
        mem_I       = calc_I(mem_avg_w,   MEM_THICKNESS) * h
        stiff_EI    = E_CARBON   * stiff_I
        mem_EI      = E_MEMBRANE * mem_I
        
        stiff_J     = calc_J(stiff_avg_w, STIFF_THICKNESS) * h
        mem_J       = calc_J(mem_avg_w,   MEM_THICKNESS) * h
        stiff_GJ    = G_CARBON   * stiff_J
        mem_GJ      = G_MEMBRANE * mem_J

        total_mass  = mem_mass  + stiff_mass
        total_EI    = mem_EI    + stiff_EI
        total_EIy   = mem_EI*mem_cg_y   + stiff_EI*stiff_cg_y
        total_GJ    = mem_GJ    + stiff_GJ
        total_GJy   = mem_GJ*mem_cg_y   + stiff_GJ*stiff_cg_y
        cg_y = (mem_mass*mem_cg_y + stiff_mass*stiff_cg_y) / total_mass if total_mass > 0 else cy_mid
        return total_mass, cx_mid, cg_y, total_EI, total_EIy, total_GJ, total_GJy

    return 0, cx_mid, cy_mid, 0, 0, 0, 0

# CORE SUB-STRIP CALCULATOR
def calc_substrip_core(x_s, x_e, global_x_offset=0.0):
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

    total_mass = total_EI = total_EIy = 0.0
    total_GJ   = total_GJy = 0.0
    sum_mcgx   = sum_mcgy  = 0.0

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

        mass, cg_x, cg_y, EI, EIy, GJ, GJy = calc_region(
            x_s, x_e, y_b, y_top, region_type, global_x_offset)

        total_mass += mass
        total_EI   += EI
        total_EIy  += EIy
        total_GJ   += GJ
        total_GJy  += GJy
        sum_mcgx   += mass * cg_x
        sum_mcgy   += mass * cg_y

    # Spar contribution
    spar_mass = w * SPAR_HEIGHT * SPAR_THICKNESS * SPAR_DENSITY
    spar_cg_x = global_x_offset + x_s + w / 2.0
    spar_cg_y = SPAR_BOTTOM + SPAR_HEIGHT / 2.0
    spar_I    = calc_I(w, SPAR_THICKNESS) * SPAR_HEIGHT
    spar_EI   = E_CARBON * spar_I
    spar_EIy  = spar_EI * spar_cg_y
    spar_J    = calc_J(w, SPAR_THICKNESS) * SPAR_HEIGHT
    spar_GJ   = G_CARBON * spar_J
    spar_GJy  = spar_GJ * spar_cg_y

    total_mass += spar_mass
    total_EI   += spar_EI
    total_EIy  += spar_EIy
    total_GJ   += spar_GJ
    total_GJy  += spar_GJy
    sum_mcgx   += spar_mass * spar_cg_x
    sum_mcgy   += spar_mass * spar_cg_y

    return total_mass, sum_mcgx, sum_mcgy, total_EI, total_EIy, total_GJ, total_GJy


# TYPE B SUB-STRIP CALCULATOR
def calc_type_b_substrip(x_s, x_e, global_x_offset=0.0):
    x_div_left  = STIFF_X_BOTTOM
    x_div_right = TYPE_B_WIDTH

    left_in  = x_s < x_div_left  < x_e
    right_in = x_s < x_div_right < x_e

    if left_in and right_in:
        p1 = calc_substrip_core(x_s,         x_div_left,  global_x_offset)
        p2 = calc_substrip_core(x_div_left,  x_div_right, global_x_offset)
        p3 = calc_substrip_core(x_div_right, x_e,         global_x_offset)
        parts = [p1, p2, p3]
    elif left_in:
        p1 = calc_substrip_core(x_s,        x_div_left, global_x_offset)
        p2 = calc_substrip_core(x_div_left, x_e,        global_x_offset)
        parts = [p1, p2]
    elif right_in:
        p1 = calc_substrip_core(x_s,         x_div_right, global_x_offset)
        p2 = calc_substrip_core(x_div_right, x_e,         global_x_offset)
        parts = [p1, p2]
    else:
        parts = [calc_substrip_core(x_s, x_e, global_x_offset)]

    total_mass = sum(p[0] for p in parts)
    sum_mcgx   = sum(p[1] for p in parts)
    sum_mcgy   = sum(p[2] for p in parts)
    total_EI   = sum(p[3] for p in parts)
    total_EIy  = sum(p[4] for p in parts)
    total_GJ   = sum(p[5] for p in parts)
    total_GJy  = sum(p[6] for p in parts)

    cg_x = sum_mcgx / total_mass
    cg_y = sum_mcgy / total_mass
    na_y = total_EIy / total_EI if total_EI > 0 else 0
    ta_y = total_GJy / total_GJ if total_GJ > 0 else 0

    return cg_x, cg_y, na_y, ta_y, total_mass

# TYPE A SUB-STRIP CALCULATOR
def calc_type_a_substrip(x_s, x_e, global_x_offset=0.0):
    w         = x_e - x_s
    cx_mid    = global_x_offset + x_s + w / 2.0

    # Membrane
    mem_mass  = w * MEM_HEIGHT * MEM_THICKNESS * MEM_DENSITY
    mem_cg_y  = MEM_HEIGHT / 2.0
    mem_I     = calc_I(w, MEM_THICKNESS) * MEM_HEIGHT
    mem_EI    = E_MEMBRANE * mem_I
    mem_EIy   = mem_EI * mem_cg_y
    mem_J     = calc_J(w, MEM_THICKNESS) * MEM_HEIGHT
    mem_GJ    = G_MEMBRANE * mem_J
    mem_GJy   = mem_GJ * mem_cg_y

    # Spar
    spar_mass = w * SPAR_HEIGHT * SPAR_THICKNESS * SPAR_DENSITY
    spar_cg_y = SPAR_BOTTOM + SPAR_HEIGHT / 2.0
    spar_I    = calc_I(w, SPAR_THICKNESS) * SPAR_HEIGHT
    spar_EI   = E_CARBON * spar_I
    spar_EIy  = spar_EI * spar_cg_y
    spar_J    = calc_J(w, SPAR_THICKNESS) * SPAR_HEIGHT
    spar_GJ   = G_CARBON * spar_J
    spar_GJy  = spar_GJ * spar_cg_y

    total_mass = mem_mass  + spar_mass
    total_EI   = mem_EI    + spar_EI
    total_EIy  = mem_EIy   + spar_EIy
    total_GJ   = mem_GJ    + spar_GJ
    total_GJy  = mem_GJy   + spar_GJy

    cg_y = (mem_mass*mem_cg_y + spar_mass*spar_cg_y) / total_mass
    na_y = total_EIy / total_EI
    ta_y = total_GJy / total_GJ

    return cx_mid, cg_y, na_y, ta_y, total_mass

# ─────────────────────────────────────────────
# FULL WING CALCULATION
# ─────────────────────────────────────────────

all_cg_x   = []
all_cg_y   = []
all_na_y   = []
all_ta_y   = []
all_masses = []
all_types  = []

print(f"\n{'Strip':<12} {'CG_x':>10} {'CG_y':>10} {'NA_y':>10} {'TA_y':>10} {'Mass':>10}")
print("-" * 65)

for strip_type, x_start, x_end, num_sub in strips:
    strip_width = x_end - x_start
    sub_width   = strip_width / num_sub

    for k in range(num_sub):
        x_s = k * sub_width
        x_e = x_s + sub_width

        if strip_type == 'A':
            cg_x, cg_y, na_y, ta_y, mass = calc_type_a_substrip(
                x_s, x_e, global_x_offset=x_start)
        else:
            cg_x, cg_y, na_y, ta_y, mass = calc_type_b_substrip(
                x_s, x_e, global_x_offset=x_start)

        all_cg_x.append(cg_x)
        all_cg_y.append(cg_y)
        all_na_y.append(na_y)
        all_ta_y.append(ta_y)
        all_masses.append(mass)
        all_types.append(strip_type)

        label = f"{strip_type}-S{k+1}"
        print(f"  {label:<10} {cg_x:>10.4f} {cg_y:>10.4f} {na_y:>10.4f} {ta_y:>10.4f} {mass:>10.4f}")

total_mass  = sum(all_masses)
overall_cgx = sum(m*x for m,x in zip(all_masses, all_cg_x)) / total_mass
overall_cgy = sum(m*y for m,y in zip(all_masses, all_cg_y)) / total_mass

print("-" * 65)
print(f"\n  Total Mass   : {total_mass:.4f} mg")
print(f"  Overall CG   : x = {overall_cgx:.4f} mm, y = {overall_cgy:.4f} mm")

# SMOOTH SPLINE
def smooth_spline(x_pts, y_pts):
    x_arr    = np.array(x_pts)
    y_arr    = np.array(y_pts)
    x_smooth = np.linspace(x_arr.min(), x_arr.max(), 500)
    spline   = make_interp_spline(x_arr, y_arr, k=3,
                                  bc_type=([(1, 0.0)], [(1, 0.0)]))
    return x_smooth, spline(x_smooth)

# PLOTTING
fig, ax = plt.subplots(1, 1, figsize=(18, 7))
fig.patch.set_facecolor('#0f1117')
ax.set_facecolor('#1a1d27')

# Full spar
ax.add_patch(patches.Rectangle((0, SPAR_BOTTOM), TOTAL_SPAN, SPAR_HEIGHT,
                                 linewidth=0.5, edgecolor='white',
                                 facecolor='#e9c46a', alpha=0.7))

for strip_type, x_start, x_end, num_sub in strips:
    if strip_type == 'A':
        ax.add_patch(patches.Rectangle((x_start, 0), x_end-x_start, MEM_HEIGHT,
                                        linewidth=0.5, edgecolor='white',
                                        facecolor='#2a9d8f', alpha=0.4))
        sub_w = (x_end - x_start) / num_sub
        for k in range(num_sub + 1):
            ax.axvline(x=x_start + k*sub_w, color='white',
                       linewidth=0.3, alpha=0.3, linestyle=':')
    else:
        xs = x_start
        if THETA_DEG == 90.0:
            ax.add_patch(patches.Rectangle((xs, 0), TYPE_B_WIDTH, MEM_HEIGHT,
                                            linewidth=0.5, edgecolor='white',
                                            facecolor='#e9c46a', alpha=0.7))
        else:
            ax.add_patch(Polygon([(xs,0),(xs,MEM_HEIGHT),(xs+STIFF_X_BOTTOM,0)],
                                  closed=True, linewidth=0.5, edgecolor='white',
                                  facecolor='#2a9d8f', alpha=0.4))
            ax.add_patch(Polygon([(xs,MEM_HEIGHT),(xs+STIFF_TOP_X,MEM_HEIGHT),
                                   (xs+TYPE_B_WIDTH,0),(xs+STIFF_X_BOTTOM,0)],
                                  closed=True, linewidth=0.5, edgecolor='white',
                                  facecolor='#e9c46a', alpha=0.7))
            ax.add_patch(Polygon([(xs+STIFF_TOP_X,MEM_HEIGHT),
                                   (xs+TYPE_B_WIDTH,MEM_HEIGHT),
                                   (xs+TYPE_B_WIDTH,0)],
                                  closed=True, linewidth=0.5, edgecolor='white',
                                  facecolor='#2a9d8f', alpha=0.4))
        sub_w = TYPE_B_WIDTH / num_sub
        for k in range(num_sub + 1):
            ax.axvline(x=xs + k*sub_w, color='white',
                       linewidth=0.3, alpha=0.3, linestyle=':')

# CG points and curve
for cx, cy, stype in zip(all_cg_x, all_cg_y, all_types):
    color = '#ff9f1c' if stype == 'A' else '#e76f51'
    ax.scatter(cx, cy, color=color, s=50, zorder=5,
               edgecolors='white', linewidths=0.5)
x_cg_s, y_cg_s = smooth_spline(all_cg_x, all_cg_y)
ax.plot(x_cg_s, y_cg_s, color='#e76f51', linewidth=2.0,
        linestyle='--', alpha=0.9, label='CG Line')

# NA points and curve
for cx, ny, stype in zip(all_cg_x, all_na_y, all_types):
    color = '#00b4d8' if stype == 'A' else '#90e0ef'
    ax.scatter(cx, ny, color=color, s=50, zorder=5,
               marker='D', edgecolors='white', linewidths=0.5)
x_na_s, y_na_s = smooth_spline(all_cg_x, all_na_y)
ax.plot(x_na_s, y_na_s, color='#00b4d8', linewidth=2.0,
        linestyle='-', alpha=0.9, label='Neutral Axis')

# TA points and curve
for cx, ty, stype in zip(all_cg_x, all_ta_y, all_types):
    color = '#06d6a0' if stype == 'A' else '#74c69d'
    ax.scatter(cx, ty, color=color, s=50, zorder=5,
               marker='^', edgecolors='white', linewidths=0.5)
x_ta_s, y_ta_s = smooth_spline(all_cg_x, all_ta_y)
ax.plot(x_ta_s, y_ta_s, color='#06d6a0', linewidth=2.0,
        linestyle=':', alpha=0.9, label='Torsional Axis')

# Overall CG
ax.scatter([overall_cgx], [overall_cgy], color='#ff4d6d', s=200,
           marker='*', zorder=6, edgecolors='white', linewidths=1)

legend_elements = [
    Patch(facecolor='#2a9d8f', alpha=0.5, label='Membrane (Polyimide)'),
    Patch(facecolor='#e9c46a', alpha=0.8, label='Spar & Stiffener (Carbon Fibre)'),
    Line2D([0],[0], color='#e76f51', linewidth=2, linestyle='--', label='CG Line'),
    Line2D([0],[0], color='#00b4d8', linewidth=2, linestyle='-',  label='Neutral Axis'),
    Line2D([0],[0], color='#06d6a0', linewidth=2, linestyle=':',  label='Torsional Axis'),
    Line2D([0],[0], marker='o', color='w', markerfacecolor='#ff9f1c',
           markersize=8, label='CG (Type A)'),
    Line2D([0],[0], marker='o', color='w', markerfacecolor='#e76f51',
           markersize=8, label='CG (Type B)'),
    Line2D([0],[0], marker='D', color='w', markerfacecolor='#00b4d8',
           markersize=8, label='NA (Type A)'),
    Line2D([0],[0], marker='D', color='w', markerfacecolor='#90e0ef',
           markersize=8, label='NA (Type B)'),
    Line2D([0],[0], marker='^', color='w', markerfacecolor='#06d6a0',
           markersize=8, label='TA (Type A)'),
    Line2D([0],[0], marker='^', color='w', markerfacecolor='#74c69d',
           markersize=8, label='TA (Type B)'),
    Line2D([0],[0], marker='*', color='w', markerfacecolor='#ff4d6d',
           markersize=12, label=f'Overall CG ({overall_cgx:.2f}, {overall_cgy:.2f})'),
]

ax.set_xlim(-0.5, TOTAL_SPAN + 0.5)
ax.set_ylim(-0.5, CHORD + 1.0)
ax.set_xlabel('x — Wingspan (mm)', color='white', fontsize=11)
ax.set_ylabel('y — Chord (mm)',    color='white', fontsize=11)
ax.set_title(f'Wing Layout — CG, Neutral Axis and Torsional Axis\n'
             f'(THETA={THETA_DEG}°, SPAR_HEIGHT={SPAR_HEIGHT}mm, STIFF_WIDTH={STIFF_WIDTH}mm)',
             color='white', fontsize=12, fontweight='bold')
ax.tick_params(colors='white')
ax.spines[:].set_color('#444')
ax.legend(handles=legend_elements, loc='upper right', fontsize=8,
          facecolor='#1a1d27', edgecolor='#444', labelcolor='white', ncol=2)

plt.suptitle(f'Rectangular Wing with Inclined Stiffeners\n'
             f'TYPE_B_WIDTH={TYPE_B_WIDTH:.3f}mm  |  WINGSPAN={TOTAL_SPAN}mm',
             color='white', fontsize=13, fontweight='bold')

plt.tight_layout()
plt.show()