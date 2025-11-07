# ActiveTrail-2D 

**ActiveTrail-2D** is a 2D simulation program written in modern Fortran
that models the collective behavior of *self-driven particles* (e.g.,
ants or active agents) moving on a discrete grid. Each particle performs
a random walk influenced by local traces (pheromone marks) left by other
particles. These traces evaporate over time, generating a self-organized
dynamics of aggregation and path formation.

## 🧠 Physical Model Overview

-   **Particles:** Represent self-propelled agents that move randomly
    while interacting indirectly through pheromone traces.
-   **Traces:** Binary field marking cells visited by particles.
-   **Evaporation:** Traces disappear probabilistically according to the
    evaporation probability `p`.
-   **Staying Rule:** When a particle moves into an unmarked cell, it
    decides whether to stay or leave according to a probability `q`.
-   **Motion:** Particle velocities and positions are integrated with a
    4th-order Runge--Kutta method including friction and random
    forcing.
-   **Boundary Conditions:** Reflective boundaries avoid particle loss
    from the simulation domain.

## ⚙️ Simulation Workflow

1.  Random initialization of particle positions.
2.  Mapping of particles into grid cells.
3.  Time evolution loop:
    -   Random forces generation.
    -   Velocity and position updates (RK4).
    -   Collision avoidance (one particle per cell).
    -   Decision to stay/leave depending on trace presence.
    -   Trace deposition and evaporation.
    -   Output of trajectories and traces.

## 📂 Output Files

| File | Description | Format |
|------|--------------|--------|
| `trajectories.csv` | Particle positions and velocities per time step | CSV (`t, id, x, y, u, v`) |
| `traces.csv` | Coordinates of marked cells (traces) | CSV (`t, x, y`) |
| `active_trail_trajectories.gif` | Animation of particle trajectories | GIF (via Gnuplot) |
| `active_trail_traces.gif` | Animation of trace points over time | GIF (via Gnuplot) |
 
## 🧩 Key Parameters

| Parameter | Meaning | Default |
|------------|----------|----------|
| `n, m` | Grid resolution | 50 × 50 |
| `Nf` | Number of particles | 10,000 |
| `dt` | Time step | 0.01 |
| `p` | Evaporation probability | 0.6 |
| `q` | Stay probability (if no trace) | 0.8 |
| `time` | Total simulated time | 20.0 |
| `out_stride` | Output frequency (steps between frames) | 10 |

## 🧮 Numerical Methods

-   **Integrator:** 4th-order Runge--Kutta (velocity and position).
-   **Random Forces:** Uniformly distributed in range `[-dx, dx]` and
    `[-dy, dy]`.
-   **One Particle per Cell Rule:** If a cell is already occupied,
    newcomers are slightly displaced according to their force vector.
-   **RNG Management:** Random seeds are initialized once at program
    start to ensure proper randomness.

## 🎞️ Visualization

Gnuplot scripts included (`postprocess.gnu`) can generate:
- Animated trajectories (`active_trail_trajectories.gif`)
- Trace evolution (`active_trail_traces.gif`)

Example usage:

``` bash
gnuplot postprocess.gnu
```


## 🧑‍💻 Authors

Developed and maintained at **LCEC--UNB** (Laboratório de Computação
Científica em Escoamentos Complexos).

Supervisor: Prof. Rafael Gabler Gontijo.

Refactored under the **ActiveTrail-2D** project (2025).

------------------------------------------------------------------------

© 2025 LCEC--UNB. Open research and educational use permitted with
citation.
