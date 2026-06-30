"""
Presentation Diagram - Two Strip Types for Vertical Stiffener Wing (4th Profile)
==================================================================================
Shows Type A and Type B strips side by side with labels and dimensions.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import Polygon, Patch
from matplotlib.lines import Line2D

# ─────────────────────────────────────────────
# CONSTANTS
# ─────────────────────────────────────────────

MEM_HEIGHT      = 6.0    # mm (y: 0 to 6)
SPAR_HEIGHT     = 1.0    # mm (y: 6 to 7)
SPAR_BOTTOM     = 6.0    # mm
TYPE_A_WIDTH    = 2.0    # mm
TYPE_B_WIDTH    = 5.0    # mm
STIFF_WIDTH     = 1.0    # mm
MEM_SIDE_WIDTH  = (TYPE_B_WIDTH - STIFF_WIDTH) / 2.0   # = 2mm each side

# CG values (from main calculation)
A_CG_X  = 1.0
A_CG_Y  = 5.3529
B_CG_X  = 2.5
B_CG_Y  = 4.3514

# Colors
COLOR_MEM   = '#2a9d8f'
COLOR_SPAR  = '#e9c46a'
COLOR_STIFF = '#e9c46a'
COLOR_CG_A  = '#ff9f1c'
COLOR_CG_B  = '#e76f51'
COLOR_TEXT  = 'white'

# ─────────────────────────────────────────────
# FIGURE SETUP
# ─────────────────────────────────────────────

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 8),
                                 gridspec_kw={'width_ratios': [2, 5]})
fig.patch.set_facecolor('#0f1117')
for ax in [ax1, ax2]:
    ax.set_facecolor('#1a1d27')
    ax.tick_params(colors=COLOR_TEXT)
    ax.spines[:].set_color('#444')
    ax.set_aspect('equal')

# ─────────────────────────────────────────────
# TYPE A STRIP (ax1)
# ─────────────────────────────────────────────

# Membrane rectangle
ax1.add_patch(patches.Rectangle((0, 0), TYPE_A_WIDTH, MEM_HEIGHT,
                                  linewidth=1.5, edgecolor=COLOR_TEXT,
                                  facecolor=COLOR_MEM, alpha=0.6))

# Leading edge spar
ax1.add_patch(patches.Rectangle((0, SPAR_BOTTOM), TYPE_A_WIDTH, SPAR_HEIGHT,
                                  linewidth=1.5, edgecolor=COLOR_TEXT,
                                  facecolor=COLOR_SPAR, alpha=0.9))

# CG point
ax1.scatter([A_CG_X], [A_CG_Y], color=COLOR_CG_A, s=150, zorder=5,
            edgecolors='white', linewidths=1.0)
ax1.annotate(f'CG\n({A_CG_X:.2f}, {A_CG_Y:.2f})',
             (A_CG_X, A_CG_Y), textcoords="offset points",
             xytext=(30, 0), ha='left', fontsize=9,
             color=COLOR_CG_A, fontweight='bold',
             arrowprops=dict(arrowstyle='->', color=COLOR_CG_A, lw=1.2))

# ── Dimension arrows ──
# Width arrow (bottom)
ax1.annotate('', xy=(TYPE_A_WIDTH, -0.4), xytext=(0, -0.4),
             arrowprops=dict(arrowstyle='<->', color=COLOR_TEXT, lw=1.2))
ax1.text(TYPE_A_WIDTH/2, -0.65, '2 mm', ha='center', va='top',
         color=COLOR_TEXT, fontsize=9)

# Membrane height arrow (left)
ax1.annotate('', xy=(-0.4, MEM_HEIGHT), xytext=(-0.4, 0),
             arrowprops=dict(arrowstyle='<->', color=COLOR_TEXT, lw=1.2))
ax1.text(-0.5, MEM_HEIGHT/2, '6 mm', ha='right', va='center',
         color=COLOR_TEXT, fontsize=9)

# Spar height arrow (left)
ax1.annotate('', xy=(-0.4, SPAR_BOTTOM + SPAR_HEIGHT), xytext=(-0.4, SPAR_BOTTOM),
             arrowprops=dict(arrowstyle='<->', color=COLOR_TEXT, lw=1.2))
ax1.text(-0.5, SPAR_BOTTOM + SPAR_HEIGHT/2, '1 mm', ha='right', va='center',
         color=COLOR_TEXT, fontsize=8)

# ── Vertex coordinates ──
coords_A = [(0,0), (2,0), (0,6), (2,6), (0,7), (2,7)]
labels_A = ['(0,0)', '(2,0)', '(0,6)', '(2,6)', '(0,7)', '(2,7)']
offsets_A = [(-0.15,-0.3), (0.15,-0.3), (-0.15,0.2), (0.15,0.2),
             (-0.15,0.2), (0.15,0.2)]
for (x,y), label, (ox,oy) in zip(coords_A, labels_A, offsets_A):
    ax1.text(x+ox, y+oy, label, ha='center', va='center',
             fontsize=7.5, color='#cccccc')

# ── Component labels ──
ax1.text(TYPE_A_WIDTH/2, MEM_HEIGHT/2, 'MEMBRANE\n(Polyimide)',
         ha='center', va='center', fontsize=9,
         color=COLOR_TEXT, fontweight='bold')
ax1.text(TYPE_A_WIDTH/2, SPAR_BOTTOM + SPAR_HEIGHT/2, 'SPAR',
         ha='center', va='center', fontsize=8,
         color='#1a1d27', fontweight='bold')

# ── Axis settings ──
ax1.set_xlim(-1.5, 4.5)
ax1.set_ylim(-1.2, 8.5)
ax1.set_xlabel('x (mm)', color=COLOR_TEXT, fontsize=10)
ax1.set_ylabel('y (mm)', color=COLOR_TEXT, fontsize=10)
ax1.set_title('Type A Strip\n(2mm, No Stiffener)',
              color=COLOR_TEXT, fontsize=12, fontweight='bold', pad=10)
ax1.set_xticks([0, 1, 2])
ax1.set_yticks([0, 1, 2, 3, 4, 5, 6, 7])
ax1.tick_params(colors=COLOR_TEXT, labelsize=8)

# ─────────────────────────────────────────────
# TYPE B STRIP (ax2)
# ─────────────────────────────────────────────

x_stiff_start = MEM_SIDE_WIDTH          # = 2.0
x_stiff_end   = x_stiff_start + STIFF_WIDTH  # = 3.0

# Membrane left (rectangle: x=0 to 2, y=0 to 6)
ax2.add_patch(patches.Rectangle((0, 0), MEM_SIDE_WIDTH, MEM_HEIGHT,
                                  linewidth=1.5, edgecolor=COLOR_TEXT,
                                  facecolor=COLOR_MEM, alpha=0.6))

# Vertical stiffener (rectangle: x=2 to 3, y=0 to 6)
ax2.add_patch(patches.Rectangle((x_stiff_start, 0), STIFF_WIDTH, MEM_HEIGHT,
                                  linewidth=1.5, edgecolor=COLOR_TEXT,
                                  facecolor=COLOR_STIFF, alpha=0.9))

# Membrane right (rectangle: x=3 to 5, y=0 to 6)
ax2.add_patch(patches.Rectangle((x_stiff_end, 0), MEM_SIDE_WIDTH, MEM_HEIGHT,
                                  linewidth=1.5, edgecolor=COLOR_TEXT,
                                  facecolor=COLOR_MEM, alpha=0.6))

# Leading edge spar (full 5mm: x=0 to 5, y=6 to 7)
ax2.add_patch(patches.Rectangle((0, SPAR_BOTTOM), TYPE_B_WIDTH, SPAR_HEIGHT,
                                  linewidth=1.5, edgecolor=COLOR_TEXT,
                                  facecolor=COLOR_SPAR, alpha=0.9))

# CG point
ax2.scatter([B_CG_X], [B_CG_Y], color=COLOR_CG_B, s=150, zorder=5,
            edgecolors='white', linewidths=1.0)
ax2.annotate(f'CG\n({B_CG_X:.2f}, {B_CG_Y:.2f})',
             (B_CG_X, B_CG_Y), textcoords="offset points",
             xytext=(40, 10), ha='left', fontsize=9,
             color=COLOR_CG_B, fontweight='bold',
             arrowprops=dict(arrowstyle='->', color=COLOR_CG_B, lw=1.2))

# ── Dimension arrows ──
# Total width (bottom)
ax2.annotate('', xy=(TYPE_B_WIDTH, -0.4), xytext=(0, -0.4),
             arrowprops=dict(arrowstyle='<->', color=COLOR_TEXT, lw=1.2))
ax2.text(TYPE_B_WIDTH/2, -0.65, '5 mm', ha='center', va='top',
         color=COLOR_TEXT, fontsize=9)

# Membrane left width (bottom)
ax2.annotate('', xy=(MEM_SIDE_WIDTH, -0.9), xytext=(0, -0.9),
             arrowprops=dict(arrowstyle='<->', color='#aaaaaa', lw=1.0))
ax2.text(MEM_SIDE_WIDTH/2, -1.1, '2 mm', ha='center', va='top',
         color='#aaaaaa', fontsize=8)

# Stiffener width (bottom)
ax2.annotate('', xy=(x_stiff_end, -0.9), xytext=(x_stiff_start, -0.9),
             arrowprops=dict(arrowstyle='<->', color=COLOR_STIFF, lw=1.0))
ax2.text(x_stiff_start + STIFF_WIDTH/2, -1.1, '1 mm', ha='center', va='top',
         color=COLOR_STIFF, fontsize=8)

# Membrane right width (bottom)
ax2.annotate('', xy=(TYPE_B_WIDTH, -0.9), xytext=(x_stiff_end, -0.9),
             arrowprops=dict(arrowstyle='<->', color='#aaaaaa', lw=1.0))
ax2.text(x_stiff_end + MEM_SIDE_WIDTH/2, -1.1, '2 mm', ha='center', va='top',
         color='#aaaaaa', fontsize=8)

# Membrane height arrow (left)
ax2.annotate('', xy=(-0.4, MEM_HEIGHT), xytext=(-0.4, 0),
             arrowprops=dict(arrowstyle='<->', color=COLOR_TEXT, lw=1.2))
ax2.text(-0.5, MEM_HEIGHT/2, '6 mm', ha='right', va='center',
         color=COLOR_TEXT, fontsize=9)

# Spar height arrow (left)
ax2.annotate('', xy=(-0.4, SPAR_BOTTOM + SPAR_HEIGHT), xytext=(-0.4, SPAR_BOTTOM),
             arrowprops=dict(arrowstyle='<->', color=COLOR_TEXT, lw=1.2))
ax2.text(-0.5, SPAR_BOTTOM + SPAR_HEIGHT/2, '1 mm', ha='right', va='center',
         color=COLOR_TEXT, fontsize=8)

# ── Vertex coordinates ──
coords_B = [(0,0), (2,0), (3,0), (5,0),
            (0,6), (2,6), (3,6), (5,6),
            (0,7), (5,7)]
labels_B = ['(0,0)', '(2,0)', '(3,0)', '(5,0)',
            '(0,6)', '(2,6)', '(3,6)', '(5,6)',
            '(0,7)', '(5,7)']
offsets_B = [(-0.15,-0.3), (0.0,-0.3), (0.0,-0.3), (0.15,-0.3),
             (-0.15, 0.2), (-0.15, 0.2), (0.15, 0.2), (0.15, 0.2),
             (-0.15, 0.2), (0.15, 0.2)]
for (x,y), label, (ox,oy) in zip(coords_B, labels_B, offsets_B):
    ax2.text(x+ox, y+oy, label, ha='center', va='center',
             fontsize=7.5, color='#cccccc')

# ── Component labels ──
ax2.text(MEM_SIDE_WIDTH/2, MEM_HEIGHT/2, 'MEM\nLEFT',
         ha='center', va='center', fontsize=8,
         color=COLOR_TEXT, fontweight='bold')
ax2.text(x_stiff_start + STIFF_WIDTH/2, MEM_HEIGHT/2,
         'V\nS\nT\nI\nF\nF',
         ha='center', va='center', fontsize=7,
         color=COLOR_TEXT, fontweight='bold')
ax2.text(x_stiff_end + MEM_SIDE_WIDTH/2, MEM_HEIGHT/2, 'MEM\nRIGHT',
         ha='center', va='center', fontsize=8,
         color=COLOR_TEXT, fontweight='bold')
ax2.text(TYPE_B_WIDTH/2, SPAR_BOTTOM + SPAR_HEIGHT/2,
         'LEADING EDGE SPAR (Carbon Fibre)',
         ha='center', va='center', fontsize=8,
         color='#1a1d27', fontweight='bold')

# T-structure annotation
ax2.annotate('T - Structure', xy=(2.5, 6.0), xytext=(4.5, 7.5),
             fontsize=9, color='#ffffff', fontweight='bold',
             arrowprops=dict(arrowstyle='->', color='white', lw=1.2))

# ── Axis settings ──
ax2.set_xlim(-1.5, 8.0)
ax2.set_ylim(-1.5, 9.0)
ax2.set_xlabel('x (mm)', color=COLOR_TEXT, fontsize=10)
ax2.set_ylabel('y (mm)', color=COLOR_TEXT, fontsize=10)
ax2.set_title('Type B Strip\n(5mm, With Vertical Stiffener)',
              color=COLOR_TEXT, fontsize=12, fontweight='bold', pad=10)
ax2.set_xticks([0, 1, 2, 3, 4, 5])
ax2.set_yticks([0, 1, 2, 3, 4, 5, 6, 7])
ax2.tick_params(colors=COLOR_TEXT, labelsize=8)

# ─────────────────────────────────────────────
# SHARED LEGEND
# ─────────────────────────────────────────────

legend_elements = [
    Patch(facecolor=COLOR_MEM,   alpha=0.6, label='Membrane (Polyimide)'),
    Patch(facecolor=COLOR_SPAR,  alpha=0.9, label='Leading Edge Spar (Carbon Fibre)'),
    Patch(facecolor=COLOR_STIFF, alpha=0.9, label='Vertical Stiffener (Carbon Fibre)'),
    Line2D([0],[0], marker='o', color='w', markerfacecolor=COLOR_CG_A,
           markersize=10, label='CG — Type A Strip'),
    Line2D([0],[0], marker='o', color='w', markerfacecolor=COLOR_CG_B,
           markersize=10, label='CG — Type B Strip'),
]
fig.legend(handles=legend_elements, loc='lower center', ncol=5,
           fontsize=9, facecolor='#1a1d27', edgecolor='#444',
           labelcolor=COLOR_TEXT, bbox_to_anchor=(0.5, 0.0))

plt.suptitle('Wing Strip Types — Vertical Stiffener Wing (4th Profile)\nCentre of Gravity Calculation',
             color=COLOR_TEXT, fontsize=14, fontweight='bold', y=1.01)

plt.tight_layout(rect=[0, 0.06, 1, 1])
plt.show()