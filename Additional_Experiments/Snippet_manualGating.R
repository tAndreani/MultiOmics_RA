# Introduction ------------------------------------------------------------

# This R snippet applies OpenCyto for manual gating 
# Configuration: >16 Cores, >64 RAM
# Data: GatingSet from CIA week8_20201117, after compensation and transformation
# Output: Saved fcs files with gated live cell events (with inverse transformation)
# Output:

install.packages("styler")

library(flowWorkspace)
library(openCyto)
library(ggcyto)

# Reload GatingSet
gs <- load_gs(file.path(tmp, "my_gs"))

# Extra: Manually add live cell gate --------------------------------------

gs_get_pop_paths(gs)
plot(gs[[1]])

# Add a single gate first
gs_add_gating_method(gs, alias = "singlets",
                     pop = "+",
                     parent = "nonDebris",
                     dims = "FSC-H,FSC-A",
                     gating_method = "singletGate",
                     gating_args = "wider_gate=FALSE,subsample_pct=0.1",
                     parallel_type = "multicore", mc.cores = 16)

# Add live cell gate
## gate_flowclust_2d: By default, the largest cluster is selected as the population of interest
gs_add_gating_method(gs,
                     alias = "liveCells",
                     pop = "-",
                     parent = "singlets",
                     dims = "DCM,FSC-A",
                     gating_method = "flowClust",
                     gating_args = "K=3,target=c(1e3,1.8e5), quantile = 0.75",
                     parallel_type = "multicore", mc.cores = 8
)

# Add CD45 gate
## gate_flowclust_2d: By default, the largest cluster is selected as the population of interest
gs_add_gating_method(gs, alias = "CD45Cells",
                     pop = "+",
                     parent = "liveCells",
                     dims = "CD45,SSC-A",
                     gating_method = "flowClust",
                     gating_args = "K=1,target=c(3e3,2.5e5), quantile = 0.5",
                     parallel_type = "multicore", mc.cores = 32
)

## tailgate: When there is only one major peak detected thus automatically disqualify the usage of mindensity
## tailgate: tol is to control how far the cut point should be placed away from the peak.
gs_add_gating_method(gs, alias = "CD45Cells",
                     pop = "+",
                     parent = "liveCells",
                     dims = "CD45",
                     gating_method = "tailgate",
                     gating_args = "side='left', tol=1e-4",
                     parallel_type = "multicore", mc.cores = 16)

## Now the gates are added to the gating tree but the actual data is not gated yet
## This is done by calling recompute method explictily
recompute(gs[[1]])

# Plot newly added gate ---------------------------------------------------

ggcyto(gs[[10]], aes(x = `FSC-H`, y = `FSC-A`), subset = "nonDebris") +
  geom_hex(bins = 6400) +
  geom_gate("singlets") +
  geom_stats(size = 6, color = "white", fill = "black", adjust = 0.3) +
  labs_cyto("marker") +
  ggcyto_par_set(
    limits = list(x = c(0, 2e6), y = c(0, 2e6)),
    hex_fill = scale_fill_viridis_c(option = "plasma", alpha = .8)
  )

ggcyto(gs[[2]], aes(x = DCM, y = `FSC-A`), subset = "singlets") +
  geom_hex(bins = 6400) +
  geom_gate("liveCells") +
  geom_stats(size = 6, color = "white", fill = "black", adjust = 0.3) +
  axis_x_inverse_trans() +
  labs_cyto("marker") +
  ggcyto_par_set(
    limits = list(x = c(0, 5e3), y = c(0, 5e5)),
    hex_fill = scale_fill_viridis_c(option = "plasma", alpha = .8)
  )

ggcyto(gs[[1]], aes(x = CD45, y = `SSC-A`), subset = "liveCells") +
  geom_hex(bins = 6400) +
  geom_gate("CD45Cells") +
  geom_stats(size = 6, color = "white", fill = "black", adjust = 0.3) +
  labs_cyto("marker") +
  ggcyto_par_set(
    limits = list(x = c(-4e3, 4e3), y = c(0, 1.5e6)),
    hex_fill = scale_fill_viridis_c(option = "plasma", alpha = .8)
  )

# Remove gates ------------------------------------------------------------

# Remove recently added new gate
gs_remove_gating_method(gs)

# Remove a specific gate
gs_pop_remove(gs[[1]], "CD45Cells")

# Get names of all nodes in GatingSet
gs_get_pop_paths(gs)
