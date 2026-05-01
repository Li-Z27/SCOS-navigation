# SCOS Navigation

R code for: "Sequential Decision-Making with Bayesian Updates for Network Traversal under Spatially Correlated Uncertainty"

## Overview

This repository contains two-stage policy learning algorithms for navigation in grid environments under spatially correlated obstacles using Gaussian Random Field priors and Bayesian belief updates.

**Algorithms implemented:**
- **TS**: Two-stage policy with Monte Carlo value estimation
- **TS-DRL (Distributional RL)**: Two-stage policy with distributional reinforcement learning

**Environment:**
- Lattice grid graph with vertices at integer coordinates
- 8-connected grid (horizontal, vertical, and diagonal edges)
- Circular obstacles that block intersecting edges

## Repository Structure

```
SCOS-navigation/
├── README.md
└── R/
    ├── helpers_shared.R     # Shared helper functions for both MC & DRL
    ├── helpers_DRL.R        # Additional helper functions for DRL
    ├── main_TS.R            # Two-stage policy (MC)
    └── main_TS_DRL.R        # Two-stage policy (DRL)
```

## Dependencies

```r
install.packages(c("igraph", "MASS", "gtools", "data.table"))
```

## Usage

Each algorithm file needs to source the helper files:

```r
# Load helper functions
source("R/helpers_shared.R")
source("R/helpers_DRL.R")  # Only needed for Policy_TS_DRL

# Prepare obstacle data frame
# obstacle_df must have columns: x, y, disambiguation_cost, prob, status

# Example for Policy_TS (Monte Carlo)
result <- Policy_TS(
  obstacle_df,      # Data frame with obstacle information
  x,                # Grid width (e.g., 50, 100)
  y,                # Grid height (e.g., 25, 50)
  r,                # Obstacle radius (e.g., 3, 5)
  s,                # Start vertex index
  g,                # Goal vertex index
  status,           # True obstacle status vector
  method_decision,  # "greedy", "epsilon-greedy", or "softmax"
  sensing_range,    # Sensor detection radius
  sense_noise,      # Sensing precision parameter (lambda, e.g., 0.25, 0.75, 1.5, 2.5)
  noise_var,        # Observation noise variance
  l,                # GRF length scale parameter
  sigma_f,          # GRF marginal standard deviation
  alpha_ig,         # Information gain scaling factor (initial)
  div_IG            # Information gain scaling factor (divisor)
)

# Access results
result$Path_cost               # Total Euclidean path length + disambiguation cost
result$Disambiguation_cost     # Total disambiguation cost
result$Record_path             # Vertex sequence of traversed path
result$Record_obstacle         # Which obstacles were disambiguated

# Example for Policy_TS_DRL (Distributional RL)
result_drl <- Policy_TS_DRL(
  obstacle_df,      # Data frame with obstacle information
  x,                # Grid width (e.g., 50, 100)
  y,                # Grid height (e.g., 25, 50)
  r,                # Obstacle radius (e.g., 3, 5)
  s,                # Start vertex index
  g,                # Goal vertex index
  status,           # True obstacle status vector
  sensing_range,    # Sensor detection radius
  sense_noise,      # Sensing precision parameter (lambda, e.g., 0.25, 0.75, 1.5, 2.5)
  noise_var,        # Observation noise variance
  l,                # GRF length scale parameter
  sigma_f,          # GRF marginal standard deviation
  replace_condi,    # Support refinement threshold (delta)
  alpha_ig,         # Information gain scaling factor (initial)
  div_IG            # Information gain scaling factor (divisor)
)
```
