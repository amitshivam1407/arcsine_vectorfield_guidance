# Arcsine Vector Field Guidance for UAV Path Following

This repository contains MATLAB implementations, simulation studies, and demonstration videos associated with the following publication:

> Amit Shivam and Ashwini Ratnoo,  
> “Arcsine Vector Field for Path Following Guidance,”  
> *Journal of Guidance, Control, and Dynamics (JGCD)*, 2023.

---

# Overview

This work proposes a geometry-driven vector field guidance method for autonomous path following of unmanned aerial vehicles (UAVs).

The proposed guidance law uses an arcsine shaping function of the path-following error to generate smooth commanded course angles with reduced curvature demand and lower control effort.

The framework is developed for:

- Straight-line path following
- Circular orbit following
- Wind-disturbed flight
- General curvature paths

The method is compared against the classical vector field guidance law of Nelson et al. and demonstrates:

- Significant reduction in maximum curvature
- Reduced control effort
- Smooth trajectory convergence
- Robustness in the presence of wind

---

# Key Contributions

- Direct geometry-driven vector field guidance
- Closed-form curvature analysis
- Reduced curvature and control effort
- Computationally efficient implementation
- Wind-resilient path following
- Extension to general curvature paths

---
## Representative Results

### Trajectory Demonstrations
<table>
  <tr>
    <td align="center">
      <b>Straight-Line Path Following</b><br>
      <img src="figures/straight_line.gif" width="430">
    </td>
    <td align="center">
      <b>Circular Orbit Following</b><br>
      <img src="figures/circular_orbit_demo.gif" width="430">
    </td>
  </tr>
</table>

# Repository Structure

```text
paper/      -> JGCD publication
figures/    -> simulation plots and illustrations
videos/     -> animated trajectory demonstrations
matlab/     -> MATLAB simulation codes
