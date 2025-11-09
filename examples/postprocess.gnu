### postprocess.gnu — ActiveTrail-2D (corrigido)

# >>> NÃO usar separador por vírgula: seus arquivos têm ESPAÇOS, não vírgulas
unset datafile separator
set datafile columnheaders   # ignora a primeira linha com cabeçalho

CSV_TRAJ   = "trajectories.csv"
CSV_TRACES = "traces.csv"

# Domínio (ajuste se necessário)
xmin = 0.0; xmax = 50.0
ymin = 0.0; ymax = 50.0

# Parâmetros de tempo (ajuste se mudou no código)
dt         = 0.01
out_stride = 5
t_step     = dt * out_stride
t_eps      = 1e-12
t_win      = t_step/2.0   # janela temporal para filtrar linhas

# Resolução e taxa
W=1000; H=900; FPS=20

# Grade do heatmap (ajuste conforme n,m)
NX=50; NY=50

# Limitar frames (0 = todos)
limit_frames = 0

# Descobre t_min/t_max a partir dos dados (coluna 1 = t)
stats CSV_TRAJ u 1 nooutput
t_min = STATS_min
t_max = STATS_max
frames = int( (t_max - t_min)/t_step + 0.5 ) + 1
if (limit_frames > 0 && limit_frames < frames) frames = limit_frames

print sprintf(">> Intervalo de tempo: [%.6f, %.6f], passo de saída: %.6f, frames: %d", t_min, t_max, t_step, frames)

# Estética base
set xrange [xmin:xmax]
set yrange [ymin:ymax]
set size ratio -1
set key off
set border lw 1.2
set grid back

########################################################
# 1) Trajetórias
########################################################
set term gif animate optimize size W,H delay int(100.0/FPS)
set output "active_trail_trajectories.gif"

do for [f=0:frames-1] {
    T = t_min + f * t_step
    set title sprintf("ActiveTrail-2D — Trajectories (t = %.3f)", T)
    # filtra linhas com t em [T - t_win, T + t_win)
    plot CSV_TRAJ u ( ($1>=T-t_win && $1<T+t_win) ? $3 : 1/0 ) : \
                    ( ($1>=T-t_win && $1<T+t_win) ? $4 : 1/0 ) \
         w p pt 6 ps 0.7 lc rgb "black"
}
unset output

########################################################
# 2) Pontos de feromônio
########################################################
set term gif animate optimize size W,H delay int(100.0/FPS)
set output "active_trail_traces.gif"

do for [f=0:frames-1] {
    T = t_min + f * t_step
    set title sprintf("ActiveTrail-2D — Pheromone traces (t = %.3f)", T)
    plot CSV_TRACES u ( ($1>=T-t_win && $1<T+t_win) ? $2 : 1/0 ) : \
                      ( ($1>=T-t_win && $1<T+t_win) ? $3 : 1/0 ) \
         w p pt 5 ps 0.7 lc rgb "#1f77b4"
}
unset output


print ">> The following gifs were generated:"
print "   - active_trail_trajectories.gif"
print "   - active_trail_traces.gif"
