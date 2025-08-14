This work is based on the paper “Continuous Piecewise Linear δ-Approximations for Univariate Functions: Computing Minimal Breakpoint Systems” by Steffen Rebennack and Josef Kallrath. The paper addresses the problem of approximating a univariate function with a continuous piecewise-linear function that stays within a specified tube of radius δ. In my implementation, I focused on two methods described in the paper: the Fixed Breakpoint System with Discretization (FBSD), which for a fixed number of breakpoints B solves an optimization problem to find the minimal δ, and the α-Forward heuristic, which for a fixed δ greedily determines the number of breakpoints required. Both methods are implemented in MATLAB — FBSD using YALMIP with a solver such as Gurobi, and the α-Forward heuristic without solver dependency — along with runtime experiments to compare their scalability and performance.

Algorithm Source: https://doi.org/10.1007/s10957-014-0687-3
  
## Repository Layout
```
├── fixedBreakpointsDeltaApproximator.m                 % FBSD implementation (solves minimal δ for fixed B using Gurobi)
├── alphaForwardHeuristc.m         % α-Forward Heuristic (greedy, solver-free approach)
├── runtimeExperiments.m  % Experiments on sin(x)/x; plots time vs. B (FBSD) and time vs. δ (α-Forward)
└── README.md              % Project documentation
```
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
