MATLAB implementations for constructing **continuous piecewise-linear (CPWL) δ-approximations** of univariate functions:

## Algorithm Source

- **Continuous Piecewise Linear δ-Approximations for Univariate Functions: Computing Minimal Breakpoint Systems**  
  *Steffen Rebennack, Josef Kallrath*  
  URL: https://doi.org/10.1007/s10957-014-0687-3
  
## Repository Layout
├── fbsd.m                 % FBSD implementation (solves minimal δ for fixed B using Gurobi)
├── alphaforward.m         % α-Forward Heuristic (greedy, solver-free approach)
├── runtime_experiments.m  % Experiments on sin(x)/x; plots time vs. B (FBSD) and time vs. δ (α-Forward)
└── README.md              % Project documentation

## Requirements

- **MATLAB** R2022a or newer
- **YALMIP** (for FBSD) — <https://yalmip.github.io/>
- **Optimzation Problem Solver** (for FBSD) — Gurobi

## How to Run

### Live Editor
Open `runtimeExperiments.m` in MATLAB Live Editor and click **Run**.  
This runs both algorithms on sample function and generates runtime scaling plots.

### Direct Script Calls

**Fixed Breakpoint System with Discretization**
```matlab
f = @(x) sin(x)
domain = [0, 2*pi];
B = 6; I = 100;
fixedBreakpointsDeltaApproximator(f, domain, B, I);
```
**α-Forward Heuristic Algorithm**
```matlab
f = @(x) sin(x)
domain = [0, 2*pi];
delta = 0.1; alpha = 0.99, D = 100;
alphaForwardHeuristic(f, domain(1),domain(2), delta, alpha, D);
