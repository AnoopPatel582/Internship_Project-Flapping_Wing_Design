import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import Polygon


# WING PARAMETERS
WINGSPAN       = 15.0    # mm (x-direction)
SPAR_HEIGHT    = 1.0     # mm (y: 0 to 1)
K              = 0.044   # curve factor
ROOT_HEIGHT    = 5.0     # mm (height at x=0 and x=15)

# Material thicknesses (z-direction)
SPAR_THICKNESS = 0.50    # mm
MEM_THICKNESS  = 0.05    # mm

# Material densities
SPAR_DENSITY   = 1.6     # mg/mm³ (Carbon Fibre)
MEM_DENSITY    = 1.3     # mg/mm³ (Polyimide)

# Discretisation
NUM_STRIPS = 10

# ─────────────────────────────────────────────
# HELPER: chord height at any x
# ─────────────────────────────────────────────

def chord_height(x):
    """Total chord height at position x."""
    return ROOT_HEIGHT + K * x * (WINGSPAN - x)


# STRIP CG CALCULATIONS
strip_width = WINGSPAN / NUM_STRIPS   # 1.5 mm each strip

strip_cg_x   = []
strip_cg_y   = []
strip_masses = []

print(f"\n{'Strip':<6} {'x_start':>8} {'x_end':>8} {'CG_x':>10} {'CG_y':>10} {'Mass(mg)':>12}")
print("-" * 65)

for i in range(NUM_STRIPS):
    x_start = i * strip_width
    x_end   = x_start + strip_width

    # Chord heights at left and right edges
    h_left  = chord_height(x_start)
    h_right = chord_height(x_end)
    h_min   = min(h_left, h_right)
    delta_h = abs(h_right - h_left)

    # ── COMPONENT 1: SPAR (Rectangle, same for all strips) ──
    spar_mass  = strip_width * SPAR_HEIGHT * SPAR_THICKNESS * SPAR_DENSITY
    spar_cg_x  = x_start + strip_width / 2.0
    spar_cg_y  = SPAR_HEIGHT / 2.0                     # = 0.5 mm

    # ── COMPONENT 2: MEMBRANE RECTANGLE ──
    mem_rect_height = h_min - SPAR_HEIGHT               # height from y=1 to h_min
    mem_rect_mass   = strip_width * mem_rect_height * MEM_THICKNESS * MEM_DENSITY
    mem_rect_cg_x   = x_start + strip_width / 2.0
    mem_rect_cg_y   = SPAR_HEIGHT + mem_rect_height / 2.0

    # ── COMPONENT 3: MEMBRANE TRIANGLE ──
    mem_tri_mass  = 0.5 * strip_width * delta_h * MEM_THICKNESS * MEM_DENSITY

    if h_right >= h_left:
        # Triangle on RIGHT side (strips 1-7)
        mem_tri_cg_x = x_start + (2.0 / 3.0) * strip_width
    else:
        # Triangle on LEFT side (strips 8-10)
        mem_tri_cg_x = x_start + (1.0 / 3.0) * strip_width

    mem_tri_cg_y = SPAR_HEIGHT + mem_rect_height + (1.0 / 3.0) * delta_h

    # ── COMBINED CG OF THIS STRIP ──
    total_mass = spar_mass + mem_rect_mass + mem_tri_mass

    cg_x = (spar_mass     * spar_cg_x +
             mem_rect_mass * mem_rect_cg_x +
             mem_tri_mass  * mem_tri_cg_x) / total_mass

    cg_y = (spar_mass     * spar_cg_y +
             mem_rect_mass * mem_rect_cg_y +
             mem_tri_mass  * mem_tri_cg_y) / total_mass

    strip_cg_x.append(cg_x)
    strip_cg_y.append(cg_y)
    strip_masses.append(total_mass)

    print(f"  {i+1:<4} {x_start:>8.2f} {x_end:>8.2f} {cg_x:>10.4f} {cg_y:>10.4f} {total_mass:>12.6f}")

# OVERALL WING CG
total_mass_wing = sum(strip_masses)
overall_cg_x = sum(m * x for m, x in zip(strip_masses, strip_cg_x)) / total_mass_wing
overall_cg_y = sum(m * y for m, y in zip(strip_masses, strip_cg_y)) / total_mass_wing

print("-" * 65)
print(f"\n  Total Wing Mass : {total_mass_wing:.6f} mg")
print(f"  Overall CG      : x = {overall_cg_x:.4f} mm,  y = {overall_cg_y:.4f} mm")

# PLOTTING
fig, ax = plt.subplots(1, 1, figsize=(10, 8))
fig.patch.set_facecolor('#0f1117')
ax.set_facecolor('#1a1d27')

# Draw SPAR (rectangle)
spar_rect = patches.Rectangle((0, 0), WINGSPAN, SPAR_HEIGHT,
                                linewidth=0, facecolor='#e9c46a',
                                alpha=0.8, label='Spar (Carbon Fibre)')
ax.add_patch(spar_rect)

# Draw MEMBRANE (curved shape using polygon)
x_vals = np.linspace(0, WINGSPAN, 300)
y_vals = [chord_height(x) for x in x_vals]

# Build polygon: top edge (y=1) then curved bottom edge
mem_x = list(x_vals) + list(x_vals[::-1])
mem_y = [SPAR_HEIGHT] * len(x_vals) + list(y_vals[::-1])
mem_polygon = Polygon(list(zip(mem_x, mem_y)), closed=True,
                       linewidth=0, facecolor='#2a9d8f',
                       alpha=0.5, label='Membrane (Polyimide)')
ax.add_patch(mem_polygon)

# Draw curved bottom edge line
ax.plot(x_vals, y_vals, color='white', linewidth=1.5, alpha=0.6)

# Draw strip boundaries
for i in range(NUM_STRIPS + 1):
    x = i * strip_width
    ax.axvline(x=x, color='white', linewidth=0.5, alpha=0.3)

# Plot CG of each strip
ax.scatter(strip_cg_x, strip_cg_y, color='#e76f51', s=80, zorder=5,
           edgecolors='white', linewidths=0.8, label='Strip CG')

# Connect CG points
ax.plot(strip_cg_x, strip_cg_y, color='#e76f51', linewidth=1.5,
        linestyle='--', alpha=0.7)

# Label each strip CG
for i, (cx, cy) in enumerate(zip(strip_cg_x, strip_cg_y)):
    ax.annotate(f'{i+1}', (cx, cy), textcoords="offset points",
                xytext=(0, 8), ha='center', fontsize=7,
                color='white', fontweight='bold')

# Overall CG
ax.scatter([overall_cg_x], [overall_cg_y], color='#ff4d6d', s=200,
           marker='*', zorder=6, edgecolors='white', linewidths=1,
           label=f'Overall CG ({overall_cg_x:.2f}, {overall_cg_y:.2f})')

ax.set_xlim(-0.5, WINGSPAN + 0.5)
ax.set_ylim(-0.5, 8.5)
ax.set_xlabel('x — Wingspan (mm)', color='white', fontsize=11)
ax.set_ylabel('y — Chord (mm)', color='white', fontsize=11)
ax.set_title('Wing Layout with Strip CG Points', color='white', fontsize=12, fontweight='bold')
ax.tick_params(colors='white')
ax.spines[:].set_color('#444')
ax.invert_yaxis()
ax.legend(loc='lower left', fontsize=8, facecolor='#1a1d27',
          edgecolor='#444', labelcolor='white')

plt.suptitle('Curved Wing — Centre of Gravity Analysis\n(Discretisation Method, 10 Strips)',
             color='white', fontsize=13, fontweight='bold')

plt.tight_layout()
plt.show()