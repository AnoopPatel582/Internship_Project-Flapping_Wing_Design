
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import Polygon, Patch
from matplotlib.lines import Line2D
from scipy.interpolate import make_interp_spline

# Wing geometry
CHORD           = 7.0    # mm (y: 0 to 7)
MEM_HEIGHT      = 6.0    # mm (y: 0 to 6)
SPAR_HEIGHT     = 1.0    # mm (y: 6 to 7)
SPAR_BOTTOM     = 6.0    # mm

# Stiffener angle
THETA           = np.radians(50.0)   # 50 degrees in radians

# Strip dimensions
TYPE_A_WIDTH    = 2.0    # mm
STIFF_WIDTH     = 2.0    # mm (perpendicular width)

# Derived geometry
STIFF_X_BOTTOM  = MEM_HEIGHT / np.tan(THETA)              # = 5.035 mm
STIFF_TOP_X     = STIFF_WIDTH / np.sin(THETA)             # = 1.3054 mm
TYPE_B_WIDTH    = STIFF_X_BOTTOM + STIFF_TOP_X            # = 6.3054 mm

# Total wingspan
NUM_PAIRS       = 5
WINGSPAN        = NUM_PAIRS * (TYPE_A_WIDTH + TYPE_B_WIDTH)

# Material thicknesses (z-direction)
SPAR_THICKNESS  = 0.5    # mm
MEM_THICKNESS   = 0.05   # mm
STIFF_THICKNESS = 0.5    # mm

# Material densities
SPAR_DENSITY    = 1.6    # mg/mm³
MEM_DENSITY     = 1.3    # mg/mm³
STIFF_DENSITY   = 1.6    # mg/mm³

# Young's Modulus (GPa)
E_CARBON        = 80.0   # GPa
E_MEMBRANE      = 5.0    # GPa

h = MEM_HEIGHT   # = 6

# Left triangle integrals
left_int_wy  = STIFF_X_BOTTOM * h / 2.0
left_int_wyy = STIFF_X_BOTTOM * h**2 / 6.0

# Right triangle integrals
right_int_wy  = STIFF_X_BOTTOM * h / 2.0
right_int_wyy = STIFF_X_BOTTOM * h**2 / 3.0

# EA contributions from triangles
left_EA      = E_MEMBRANE * MEM_THICKNESS * left_int_wy
left_EAy     = E_MEMBRANE * MEM_THICKNESS * left_int_wyy
right_EA     = E_MEMBRANE * MEM_THICKNESS * right_int_wy
right_EAy    = E_MEMBRANE * MEM_THICKNESS * right_int_wyy

# STRIP CG CALCULATIONS

strip_cg_x   = []
strip_cg_y   = []
strip_na_y   = []
strip_masses = []
strip_types  = []

print(f"\n{'Strip':<6} {'Type':<6} {'x_start':>8} {'x_end':>8} {'CG_x':>10} {'CG_y':>10} {'NA_y':>10} {'Mass(mg)':>10}")
print("-" * 80)

x_current = 0.0

for i in range(NUM_PAIRS):

    # ── TYPE A STRIP (2mm, no stiffener) ──
    x_start = x_current
    x_end   = x_start + TYPE_A_WIDTH

    # Component 1: Leading edge spar
    spar_mass  = TYPE_A_WIDTH * SPAR_HEIGHT * SPAR_THICKNESS * SPAR_DENSITY
    spar_cg_x  = x_start + TYPE_A_WIDTH / 2.0
    spar_cg_y  = SPAR_BOTTOM + SPAR_HEIGHT / 2.0
    spar_A     = TYPE_A_WIDTH * SPAR_THICKNESS
    spar_EA    = E_CARBON * spar_A

    # Component 2: Membrane (rectangle)
    mem_mass   = TYPE_A_WIDTH * MEM_HEIGHT * MEM_THICKNESS * MEM_DENSITY
    mem_cg_x   = x_start + TYPE_A_WIDTH / 2.0
    mem_cg_y   = MEM_HEIGHT / 2.0
    mem_A      = TYPE_A_WIDTH * MEM_THICKNESS
    mem_EA     = E_MEMBRANE * mem_A

    # ── Combined CG ──
    total_mass = spar_mass + mem_mass
    cg_x = (spar_mass * spar_cg_x + mem_mass * mem_cg_x) / total_mass
    cg_y = (spar_mass * spar_cg_y + mem_mass * mem_cg_y) / total_mass

    # ── Combined NA ──
    total_EA = spar_EA + mem_EA
    na_y     = (spar_EA * spar_cg_y + mem_EA * mem_cg_y) / total_EA

    strip_cg_x.append(cg_x)
    strip_cg_y.append(cg_y)
    strip_na_y.append(na_y)
    strip_masses.append(total_mass)
    strip_types.append('A')

    print(f"  {2*i+1:<4} {'A':<6} {x_start:>8.2f} {x_end:>8.2f} {cg_x:>10.4f} {cg_y:>10.4f} {na_y:>10.4f} {total_mass:>10.6f}")

    x_current = x_end

    # ── TYPE B STRIP (with inclined stiffener) ──
    x_start = x_current
    x_end   = x_start + TYPE_B_WIDTH

    # Component 1: Left Triangle (membrane)
    # Vertices: (0,0), (0,6), (STIFF_X_BOTTOM, 0)
    left_tri_area  = 0.5 * STIFF_X_BOTTOM * MEM_HEIGHT
    left_tri_mass  = left_tri_area * MEM_THICKNESS * MEM_DENSITY
    left_tri_cg_x  = x_start + STIFF_X_BOTTOM / 3.0
    left_tri_cg_y  = MEM_HEIGHT / 3.0               # = 2.0 mm

    # Component 2: Right Triangle (membrane)
    # Vertices: (STIFF_TOP_X, 6), (TYPE_B_WIDTH, 6), (TYPE_B_WIDTH, 0)
    right_tri_mass = left_tri_mass                   # same area as left triangle
    right_tri_cg_x = (x_start + STIFF_TOP_X +
                      x_start + TYPE_B_WIDTH +
                      x_start + TYPE_B_WIDTH) / 3.0
    right_tri_cg_y = (MEM_HEIGHT + MEM_HEIGHT + 0.0) / 3.0   # = 4.0 mm

    # Component 3: Inclined Stiffener (parallelogram)
    # Vertices: (0,6), (STIFF_TOP_X,6), (TYPE_B_WIDTH,0), (STIFF_X_BOTTOM,0)
    stiff_length   = np.sqrt(MEM_HEIGHT**2 + STIFF_X_BOTTOM**2)
    stiff_area     = stiff_length * STIFF_WIDTH
    stiff_mass     = stiff_area * STIFF_THICKNESS * STIFF_DENSITY
    stiff_cg_x     = x_start + (0 + STIFF_TOP_X + TYPE_B_WIDTH + STIFF_X_BOTTOM) / 4.0
    stiff_cg_y     = MEM_HEIGHT / 2.0               # = 3.0 mm
    stiff_A        = STIFF_WIDTH * STIFF_THICKNESS
    stiff_EA       = E_CARBON * stiff_A

    # Component 4: Leading edge spar
    spar_mass      = TYPE_B_WIDTH * SPAR_HEIGHT * SPAR_THICKNESS * SPAR_DENSITY
    spar_cg_x      = x_start + TYPE_B_WIDTH / 2.0
    spar_cg_y      = SPAR_BOTTOM + SPAR_HEIGHT / 2.0
    spar_A         = TYPE_B_WIDTH * SPAR_THICKNESS
    spar_EA        = E_CARBON * spar_A

    # ── Combined CG ──
    total_mass = left_tri_mass + right_tri_mass + stiff_mass + spar_mass
    cg_x = (left_tri_mass  * left_tri_cg_x  +
             right_tri_mass * right_tri_cg_x +
             stiff_mass     * stiff_cg_x     +
             spar_mass      * spar_cg_x) / total_mass
    cg_y = (left_tri_mass  * left_tri_cg_y  +
             right_tri_mass * right_tri_cg_y +
             stiff_mass     * stiff_cg_y     +
             spar_mass      * spar_cg_y) / total_mass

    # ── Combined NA ──
    total_EA_b  = left_EA + right_EA + stiff_EA + spar_EA
    total_EAy_b = (left_EAy + right_EAy +
                   stiff_EA * stiff_cg_y +
                   spar_EA  * spar_cg_y)
    na_y = total_EAy_b / total_EA_b

    strip_cg_x.append(cg_x)
    strip_cg_y.append(cg_y)
    strip_na_y.append(na_y)
    strip_masses.append(total_mass)
    strip_types.append('B')

    print(f"  {2*i+2:<4} {'B':<6} {x_start:>8.2f} {x_end:>8.2f} {cg_x:>10.4f} {cg_y:>10.4f} {na_y:>10.4f} {total_mass:>10.6f}")

    x_current = x_end

# OVERALL WING CG AND NA

total_mass_wing = sum(strip_masses)
overall_cg_x    = sum(m * x for m, x in zip(strip_masses, strip_cg_x)) / total_mass_wing
overall_cg_y    = sum(m * y for m, y in zip(strip_masses, strip_cg_y)) / total_mass_wing
overall_na_y    = sum(strip_na_y) / len(strip_na_y)

print("-" * 80)
print(f"\n  Total Wing Mass  : {total_mass_wing:.6f} mg")
print(f"  Overall CG       : x = {overall_cg_x:.4f} mm,  y = {overall_cg_y:.4f} mm")
print(f"  Average NA_y     : y = {overall_na_y:.4f} mm")

# HELPER: smooth spline

def smooth_spline(x_pts, y_pts):
    x_arr    = np.array(x_pts)
    y_arr    = np.array(y_pts)
    x_smooth = np.linspace(x_arr.min(), x_arr.max(), 500)
    spline   = make_interp_spline(x_arr, y_arr, k=3,
                                   bc_type=([(1, 0.0)], [(1, 0.0)]))
    return x_smooth, spline(x_smooth)

# PLOTTING

fig, ax = plt.subplots(1, 1, figsize=(14, 7))
fig.patch.set_facecolor('#0f1117')
ax.set_facecolor('#1a1d27')

# ── Draw wing layout ──
x_current = 0.0
for i in range(NUM_PAIRS):

    # Type A: membrane rectangle + spar
    ax.add_patch(patches.Rectangle((x_current, 0), TYPE_A_WIDTH, MEM_HEIGHT,
                                    linewidth=0.5, edgecolor='white',
                                    facecolor='#2a9d8f', alpha=0.4))
    ax.add_patch(patches.Rectangle((x_current, SPAR_BOTTOM), TYPE_A_WIDTH, SPAR_HEIGHT,
                                    linewidth=0.5, edgecolor='white',
                                    facecolor='#e9c46a', alpha=0.7))
    x_current += TYPE_A_WIDTH

    xs = x_current

    # Left triangle (membrane)
    ax.add_patch(Polygon([(xs, 0), (xs, MEM_HEIGHT), (xs + STIFF_X_BOTTOM, 0)],
                          closed=True, linewidth=0.5, edgecolor='white',
                          facecolor='#2a9d8f', alpha=0.4))

    # Inclined stiffener (parallelogram)
    ax.add_patch(Polygon([(xs, MEM_HEIGHT),
                           (xs + STIFF_TOP_X, MEM_HEIGHT),
                           (xs + TYPE_B_WIDTH, 0),
                           (xs + STIFF_X_BOTTOM, 0)],
                          closed=True, linewidth=0.5, edgecolor='white',
                          facecolor='#e9c46a', alpha=0.7))

    # Right triangle (membrane)
    ax.add_patch(Polygon([(xs + STIFF_TOP_X, MEM_HEIGHT),
                           (xs + TYPE_B_WIDTH, MEM_HEIGHT),
                           (xs + TYPE_B_WIDTH, 0)],
                          closed=True, linewidth=0.5, edgecolor='white',
                          facecolor='#2a9d8f', alpha=0.4))

    # Leading edge spar (full TYPE_B_WIDTH)
    ax.add_patch(patches.Rectangle((xs, SPAR_BOTTOM), TYPE_B_WIDTH, SPAR_HEIGHT,
                                    linewidth=0.5, edgecolor='white',
                                    facecolor='#e9c46a', alpha=0.7))
    x_current += TYPE_B_WIDTH

# ── Plot CG points ──
for idx, (cx, cy, stype) in enumerate(zip(strip_cg_x, strip_cg_y, strip_types)):
    color = '#ff9f1c' if stype == 'A' else '#e76f51'
    ax.scatter(cx, cy, color=color, s=80, zorder=5,
               edgecolors='white', linewidths=0.8)
    ax.annotate(f'{idx+1}', (cx, cy), textcoords="offset points",
                xytext=(0, 8), ha='center', fontsize=7,
                color='white', fontweight='bold')

# ── CG smooth curve ──
x_cg_smooth, y_cg_smooth = smooth_spline(strip_cg_x, strip_cg_y)
ax.plot(x_cg_smooth, y_cg_smooth, color='#e76f51', linewidth=2.0,
        linestyle='--', alpha=0.9)

# ── Plot NA points ──
for idx, (cx, ny, stype) in enumerate(zip(strip_cg_x, strip_na_y, strip_types)):
    color = '#00b4d8' if stype == 'A' else '#90e0ef'
    ax.scatter(cx, ny, color=color, s=80, zorder=5,
               marker='D', edgecolors='white', linewidths=0.8)

# ── NA smooth curve ──
x_na_smooth, y_na_smooth = smooth_spline(strip_cg_x, strip_na_y)
ax.plot(x_na_smooth, y_na_smooth, color='#00b4d8', linewidth=2.0,
        linestyle='-', alpha=0.9)

# ── Overall CG ──
ax.scatter([overall_cg_x], [overall_cg_y], color='#ff4d6d', s=200,
           marker='*', zorder=6, edgecolors='white', linewidths=1)

# ── Legend ──
legend_elements = [
    Patch(facecolor='#2a9d8f', alpha=0.5,  label='Membrane (Polyimide)'),
    Patch(facecolor='#e9c46a', alpha=0.8,  label='Spar & Inclined Stiffener (Carbon Fibre)'),
    Line2D([0],[0], color='#e76f51', linewidth=2, linestyle='--', label='CG Line'),
    Line2D([0],[0], color='#00b4d8', linewidth=2, linestyle='-',  label='Neutral Axis'),
    Line2D([0],[0], marker='o', color='w', markerfacecolor='#ff9f1c',
           markersize=8, label='Strip CG (Type A)'),
    Line2D([0],[0], marker='o', color='w', markerfacecolor='#e76f51',
           markersize=8, label='Strip CG (Type B)'),
    Line2D([0],[0], marker='D', color='w', markerfacecolor='#00b4d8',
           markersize=8, label='Strip NA (Type A)'),
    Line2D([0],[0], marker='D', color='w', markerfacecolor='#90e0ef',
           markersize=8, label='Strip NA (Type B)'),
    Line2D([0],[0], marker='*', color='w', markerfacecolor='#ff4d6d',
           markersize=12, label=f'Overall CG ({overall_cg_x:.2f}, {overall_cg_y:.2f})'),
]

ax.set_xlim(-0.5, WINGSPAN + 0.5)
ax.set_ylim(-0.5, CHORD + 1.0)
ax.set_xlabel('x — Wingspan (mm)', color='white', fontsize=11)
ax.set_ylabel('y — Chord (mm)', color='white', fontsize=11)
ax.set_title(f'Wing Layout with CG Line and Neutral Axis (θ = {np.degrees(THETA):.0f}°)',
             color='white', fontsize=12, fontweight='bold')
ax.tick_params(colors='white')
ax.spines[:].set_color('#444')
ax.legend(handles=legend_elements, loc='upper right', fontsize=8,
          facecolor='#1a1d27', edgecolor='#444', labelcolor='white', ncol=2)

plt.suptitle(f'Rectangular Wing with Inclined Stiffeners (θ = {np.degrees(THETA):.0f}°)\nCG Line vs Neutral Axis (Y-Direction Bending)',
             color='white', fontsize=13, fontweight='bold')

plt.tight_layout()
plt.show()