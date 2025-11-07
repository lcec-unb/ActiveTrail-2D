# ActiveTrail-2D — Examples

This folder contains example visualizations generated from **ActiveTrail-2D** simulation outputs.

## 🔧 How it works

1. Copy the output files `trajectories.csv` and `traces.csv` from your simulation into this folder and place them here.
2. Run the script:

   ```bash
   ./Allrun.sh
   ```

This script automatically calls Gnuplot to execute `postprocess.gnu` and generate animation files.

### 🎞️ Output

After running the script, you should see:

`active_trail_trajectories.gif` — particle trajectories

`active_trail_traces.gif` — pheromone trace evolution

🧠 Note: Gnuplot must be installed on your system and available in the terminal path.
