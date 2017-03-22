# mav_trajectory_generation
This repository contains tools for polynomial trajectory generation and optimization based on methods described in [1].
These techniques are especially suitable for rotary-wing micro aerial vehicles (MAVs).
This README provides a brief overview of our trajectory generation utilities with some examples.

**Authors**: Markus Achtelik, Michael Burri, Helen Oleynikova, Rik Bähnemann, Marija Popović  
**Maintainer**: Rik Bähnemann, brik@ethz.ch  
**Affiliation**: Autonomous Systems Lab, ETH Zurich

## Basics
A **vertex** describes the properties of a support point of a **polynomial** path. Pairs of vertices are connected together to form **segments**.
Each vertex has a set of constraints: the values of position derivatives that must be matched during optimization procedures.
Typically, only the positions are specified for vertices along a path, while start and end vertices have all derivatives of position set to zero.
```
  x----------x-----------------x
vertex            segment
```
In the case of multi-dimensional vertices, the derivative constraints exist in all dimensions, with possibly different values.

## Bibliography
[1] C. Richter, A. Bry, and N. Roy, “Polynomial trajectory planning for quadrotor flight,” in International Conference on Robotics and Automation. Singapore: Springer, 2013.
