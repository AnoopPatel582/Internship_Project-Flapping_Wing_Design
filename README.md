# Flapping Wing Morphological and Aeroelastic Analysis

## Overview
This repository contains analytical models, calculation frameworks, and design methodologies for evaluating the aerodynamic morphology and structural dynamics of microrobotic flapping wings. 

The theoretical foundation of this project is heavily inspired by **Yufeng Chen's PhD dissertation** (Harvard University, 2017): *"Experimental and Computational Study of Flapping-Wing Dynamics and Locomotion in Aerial and Aquatic Environments"*, as well as Ellington's foundational studies on insect wing shape parameterization.

## Motivation
Researchers across the world are continuously trying to build a MAV(micro-aerial vehicle). A small insect-like flying robot that can fly like an insect by flapping its wings. Harvard researchers design a 80mg MAV. These robots are flying by power supply with tether means a very thin wire. 
These robots would be very useful in calamitous disasters like earthquakes, floods, and other accidents to save human lives. The main problem for these MAVs is that they are flying on power supply by the tether. We are trying to figure out the optimal wing design that can generate more lift and carry more payload capacity so that it can carry its own power source.

## Problem Statement
5 key parameters affect the aerodynamic performance of the wing. These are as follows.
1. Intertia
2. Aspect Ratio (AR)
3. First Area Moment ($\hat{r}_1$)
4. Leading Edge Sweep Ratio (LESR)
5. Size of the wing
    
We will study how these parameters affect the aerodynamic performance. We will try to design a wing that can carry more payload, generate more lift and dissipate less power. Our ultimate goal is to make a machine learning model and train it over thousands of data points so that it can predict wing performance when we change the wing parameters.

## Key Features & Capabilities

Our models process 2D/3D wing geometries and composite material properties to calculate critical parameters for flapping-wing flight:

### 1. Aerodynamic Morphology (Discretization Method)
By applying the strip/discretization method (slicing the wing into vertical rectangular strips), the project calculates:
* **Center of Gravity (CG)**: Determines the balance point of the wing.
* **First Area Moment (Spanwise Area Moment, $\hat{r}_1$)**: A purely geometric parameter dictating where the bulk of the wing's surface area is located along the span. The normalized version is critical for comparing wings of different scales (e.g., robotic vs. biological).
* **Second Moment of Area**: Crucial for understanding the distribution of area relative to the rotational axis.

### 2. Structural Aeroelasticity
Calculates structural axes to predict the passive aeroelastic response (bending and twisting) of the wing during the flapping stroke:
* **Elastic Neutral Axis**: Determined using stiffness-weighted calculations ($E \times I$) to find the axis where the wing undergoes zero stress during bending.
* **Torsional Neutral Axis**: Calculated using the Modulus of Rigidity ($G$) and the Polar Moment of Inertia ($J$) to find the spanwise line about which the wing twists without plunging.

### 3. Advanced Composite Wing Structures
The project supports modelling for multi-material, stiffened wing architectures, including:
* **Membrane and Spar Composites**: Calculations account for different material properties between the flexible membrane and the rigid leading edge.
* **Inclined Stiffeners & Chordwise Ribs**: Evaluates the shift from a 1D structural problem to advanced aeroelasticity by tracking both Local and Global Neutral Axes when carbon fiber ribs are introduced (e.g., ribs placed at specific spanwise coordinates like $x=5$ and $x=10$).

## Material Properties Used in Current Models
Our default structural calculations utilize the following material specifications (typical for carbon-fiber and polymer microrobots):

* **Membrane:**
    * Young's Modulus ($E_m$): 5 GPa
    * Modulus of Rigidity / Shear Modulus ($G_m$): 1.5 GPa
* **Structural Spars / Ribs (Carbon Fiber):**
    * Young's Modulus ($E_s$): 80 GPa
    * Modulus of Rigidity / Shear Modulus ($G_s$): 27 GPa

## Mathematical Background

**Torsional Axis Calculation:**
The torsional axis is found using the summation of torsional stiffness contributions:
$$X_{torsion} = \frac{\sum (G_i \cdot J_i \cdot x_i)}{\sum (G_i \cdot J_i)}$$
*(Where G is the modulus of rigidity and J is the polar moment of inertia for the cross-section components).*

**Spanwise Area Moment ($\hat{r}_1$):**
Calculated by holding the root and tip boundaries constant, mathematically controlling the trailing edge curve factor ($k$) to analyze its direct effect on aerodynamic performance.

## Future Work
* Integration of 2D and 3D Computational Fluid Dynamics (CFD) validations.
* Simulation of wing dynamics in both aerial and aquatic environments.
* Further optimization of inclined stiffener angles (e.g., 45-degree carbon fiber laminates) for targeted wing-pitching profiles.

## References
* Chen, Y. (2017). *Experimental and Computational Study of Flapping-Wing Dynamics and Locomotion in Aerial and Aquatic Environments*. [PhD Thesis, Harvard University].
* Ellington, C. P. (1984). *The Aerodynamics of Hovering Insect Flight*.
