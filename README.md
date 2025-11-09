# ActiveTrail-2D

**ActiveTrail-2D** is a 2D simulation program written in modern Fortran that models the collective dynamics of *self-driven particles* (e.g., ants or active agents) moving on a discrete grid and interacting indirectly through *pheromone-like fields*. The model combines stochastic motion, local interactions, and probabilistic behavioral rules, leading to emergent collective organization and trail formation.

---

## 🧠 Physical and Mathematical Model

The system represents a population of $N_f$ active agents moving on a rectangular grid $(n \times m)$ of spacing $(\Delta x, \Delta y)$, confined between reflective boundaries.

### 1. **Particle Motion**

Each particle $i$ has position $(x_i, y_i)$ and velocity $(u_i, v_i)$ evolving according to a stochastic dynamics:

$$
m \frac{du_i}{dt} = -\gamma u_i + F_i^{\text{rand}} + F_i^{\text{rep}},
$$

where:
- $m$ is the particle mass;
- $\gamma$ is the friction coefficient;
- $F_i^{\text{rand}}$ is a stochastic force (Gaussian white noise in the Euler–Maruyama version);
- $F_i^{\text{rep}}$ is a short-range repulsive force preventing particle overlap.

In the v0.1 base version, random forces are uniform in $[-\Delta x, \Delta x]$, while velocities and positions are updated via a 4th-order Runge–Kutta scheme.  
In the extended versions (v0.2+), the update follows the **Euler–Maruyama** stochastic integrator with Gaussian noise.

### 2. **Repulsive Interaction**

Particles interact only when they are within a cutoff distance $d_0$. The deterministic repulsive force follows a linear penalty model:

$$
\mathbf{F}_{ij}^{\text{rep}} = k_{\text{rep}}(d_0 - r_{ij})\hat{\mathbf{r}}_{ij}, \quad r_{ij} < d_0,
$$

where $k_{\text{rep}}$ is a stiffness coefficient.  
The total deterministic force on particle $i$ is $\mathbf{F}_i^{\text{rep}} = \sum_j \mathbf{F}_{ij}^{\text{rep}}$.

### 3. **Pheromone (Trace) Field**

A binary field $F(x,y,t)$ marks the grid cells visited by particles:
- When a particle occupies a cell, that cell’s pheromone value becomes 1 (deposition rule).
- When a cell is empty, its pheromone value evaporates with probability $p$ per time step (evaporation rule).

This dynamic field provides indirect communication among agents — creating a memory of past trajectories that guides future motion statistically.

### 4. **Behavioral Rule: Stay or Leave**

When a particle enters a cell **without pheromone**, it must decide whether to *stay* or *move on*:
- With probability **q**, it stays (remains in that cell for the next step);
- With probability **1–q**, it receives a small random “kick” and moves away.

If the cell already contains pheromone, the particle always stays (reinforcement of trail-following). Here, we set $q < 0.5$ so the particle has a higher probability to stay in a cell that has a pheromone trace.

### 5. **Boundary Conditions**

Reflective conditions ensure confinement:
$$
x_i, y_i \in [x_{\min}, x_{\max}] \times [y_{\min}, y_{\max}].
$$

Whenever a particle attempts to exit the domain, its position is slightly displaced inward.

---

## ⚙️ Simulation Algorithm

1. Random initialization of particle positions.
2. Mapping of particles into discrete grid cells.
3. At each time step:
   - Generate stochastic forces (random or Gaussian).
   - Compute short-range repulsion (if enabled).
   - Update velocities and positions (RK4 or EM).
   - Handle collisions (one particle per cell).
   - Apply pheromone-based behavioral rules (stay/leave).
   - Evaporate pheromone field.
   - Record particle and field states to output.

---

## 📂 Output Files

| File | Description | Format |
|------|--------------|--------|
| `trajectories.csv` | Particle positions and velocities per time step | CSV (`t, id, x, y, u, v`) |
| `traces.csv` | Coordinates of marked cells (pheromone traces) | CSV (`t, x, y`) |
| `active_trail_trajectories.gif` | Animation of particle trajectories | GIF (via Gnuplot) |
| `active_trail_traces.gif` | Animation of trace points over time | GIF (via Gnuplot) |

---

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

---

## 🧮 Numerical and Algorithmic Details

- **Integrator:** 4th-order Runge–Kutta (v0.1); Euler–Maruyama (v0.2+).  
- **Random Forces:** Uniform (`[-dx, dx]`) or Gaussian (`σ√dt`).  
- **One-Particle-per-Cell Rule:** Prevents overlap by local displacement.  
- **RNG:** Seeded once for consistent random sequences.  
- **Boundary Handling:** Reflective displacement inward.  
- **Output:** CSV files for postprocessing and animation via Gnuplot.

---

## 🎞️ Visualization

The included `postprocess.gnu` script can generate:

- Animated trajectories (`active_trail_trajectories.gif`)
- Pheromone trace evolution (`active_trail_traces.gif`)

Usage:

```bash
gnuplot postprocess.gnu
```

---

## 🧑‍💻 Authors

Developed and maintained at **LCEC–UNB**  
(*Laboratório de Computação Científica em Escoamentos Complexos*).

Developer: **Prof. Rafael Gabler Gontijo**


---

© 2025 LCEC–UNB. Open research and educational use permitted with citation.
