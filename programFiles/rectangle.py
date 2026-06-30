import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# Wing geometry
WINGSPAN        = 15.0      # mm  (x-direction)
MEMBRANE_HEIGHT = 6.0       # mm  (y: 0 to 6)
SPAR_HEIGHT     = 1.0       # mm  (y: 6 to 7)
TOTAL_CHORD     = 7.0       # mm  (y: 0 to 7)

# Material thicknesses (z-direction)
MEMBRANE_THICKNESS = 0.05   # mm
SPAR_THICKNESS     = 0.50   # mm

# Material densities
MEMBRANE_DENSITY = 1.3      # mg/mm³  (Polyimide)
SPAR_DENSITY     = 1.6      # mg/mm³  (Carbon Fibre)

# Discretisation
NUM_STRIPS = 10

# STRIP CG CALCULATIONS
strip_width = WINGSPAN / NUM_STRIPS   # 1.5 mm each strip

strip_cg_x = []
strip_cg_y = []
strip_masses = []

print(f"\n{'Strip':<6} {'x_start':>8} {'x_end':>8} {'CG_x':>8} {'CG_y':>8} {'Mass(mg)':>10}")
print("-" * 55)

for i in range(NUM_STRIPS):
    x_start = i * strip_width
    x_end   = x_start + strip_width
    x_mid   = (x_start + x_end) / 2.0   # CG x of strip (symmetric)

    # ── Membrane part of this strip ──
    mem_volume = strip_width * MEMBRANE_HEIGHT * MEMBRANE_THICKNESS
    mem_mass   = mem_volume * MEMBRANE_DENSITY
    mem_cg_y   = MEMBRANE_HEIGHT / 2.0   # = 3.0 mm

    # ── Spar part of this strip ──
    spar_volume = strip_width * SPAR_HEIGHT * SPAR_THICKNESS
    spar_mass   = spar_volume * SPAR_DENSITY
    spar_cg_y   = MEMBRANE_HEIGHT + SPAR_HEIGHT / 2.0   # = 6.5 mm

    # ── Combined CG of this strip ──
    total_mass = mem_mass + spar_mass
    cg_y = (mem_mass * mem_cg_y + spar_mass * spar_cg_y) / total_mass
    cg_x = x_mid   # strips are uniform along x, so CG is at midpoint

    strip_cg_x.append(cg_x)
    strip_cg_y.append(cg_y)
    strip_masses.append(total_mass)

    print(f"  {i+1:<4} {x_start:>8.2f} {x_end:>8.2f} {cg_x:>8.3f} {cg_y:>8.4f} {total_mass:>10.5f}")

# OVERALL WING CG
total_mass_wing = sum(strip_masses)
overall_cg_x = sum(m * x for m, x in zip(strip_masses, strip_cg_x)) / total_mass_wing
overall_cg_y = sum(m * y for m, y in zip(strip_masses, strip_cg_y)) / total_mass_wing

print("-" * 55)
print(f"\n  Total Wing Mass : {total_mass_wing:.5f} mg")
print(f"  Overall CG      : x = {overall_cg_x:.4f} mm,  y = {overall_cg_y:.4f} mm")

# PLOTTING
fig, ax1 = plt.subplots(1, 1, figsize=(8, 6))
fig.patch.set_facecolor('#0f1117')
ax1.set_facecolor('#1a1d27')

# ── Wing layout with strips and CG points ──

# Draw membrane
mem_rect = patches.Rectangle((0, 0), WINGSPAN, MEMBRANE_HEIGHT,
                               linewidth=0, facecolor='#2a9d8f', alpha=0.5, label='Membrane (Polyimide)')
ax1.add_patch(mem_rect)

# Draw spar
spar_rect = patches.Rectangle((0, MEMBRANE_HEIGHT), WINGSPAN, SPAR_HEIGHT,
                                linewidth=0, facecolor='#e9c46a', alpha=0.7, label='Spar (Carbon Fibre)')
ax1.add_patch(spar_rect)

# Draw strip boundaries
for i in range(NUM_STRIPS + 1):
    x = i * strip_width
    ax1.axvline(x=x, color='white', linewidth=0.5, alpha=0.3)

# Plot CG of each strip
ax1.scatter(strip_cg_x, strip_cg_y, color='#e76f51', s=80, zorder=5,
            edgecolors='white', linewidths=0.8, label='Strip CG')

# Connect CG points
ax1.plot(strip_cg_x, strip_cg_y, color='#e76f51', linewidth=1.5,
         linestyle='--', alpha=0.7)

# Label each strip CG
for i, (cx, cy) in enumerate(zip(strip_cg_x, strip_cg_y)):
    ax1.annotate(f'{i+1}', (cx, cy), textcoords="offset points",
                 xytext=(0, 8), ha='center', fontsize=7,
                 color='white', fontweight='bold')

# Overall CG
ax1.scatter([overall_cg_x], [overall_cg_y], color='#ff4d6d', s=200,
            marker='*', zorder=6, edgecolors='white', linewidths=1,
            label=f'Overall CG ({overall_cg_x:.2f}, {overall_cg_y:.2f})')

ax1.set_xlim(-0.5, WINGSPAN + 1.5)
ax1.set_ylim(-0.5, TOTAL_CHORD + 1.5)
ax1.set_xlabel('x — Wingspan (mm)', color='white', fontsize=11)
ax1.set_ylabel('y — Chord (mm)', color='white', fontsize=11)
ax1.set_title('Wing Layout with Strip CG Points', color='white', fontsize=11, fontweight='bold')
ax1.tick_params(colors='white')
ax1.spines[:].set_color('#444')
ax1.legend(loc='lower right', fontsize=8, facecolor='#1a1d27',
           edgecolor='#444', labelcolor='white')

# Y-axis labels
ax1.set_yticks([0, 3, 6, 6.5, 7])
ax1.set_yticklabels(['0', '3\n(mem CG)', '6', '6.5\n(spar CG)', '7'], color='white', fontsize=8)

plt.suptitle('Rectangular Wing — Centre of Gravity Analysis\n(Discretisation Method, 10 Strips)',
             color='white', fontsize=11, fontweight='bold', y=1.0)

plt.tight_layout()
plt.show()