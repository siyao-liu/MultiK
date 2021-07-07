library(splatter)
library(Matrix)
library(scater)

# set seeds
seed = 629

# set parameters
params <- newSplatParams(batchCells = 1000, nGenes=5000)
params

####################################################################
# Simulating 2 equal groups

sim.groups <- splatSimulate(params=params,
                            group.prob = c(0.5, 0.5), method = "groups",
                            verbose = FALSE)

sim.groups <- logNormCounts(sim.groups)
sim.groups <- runPCA(sim.groups)
plotPCA(sim.groups, colour_by = "Group")


#####################################################################
# Simulating 3 unequal groups
sim.groups <- splatSimulate(params=params,
                            group.prob = c(0.1,0.3,0.6), method = "groups",
                            verbose = FALSE)
dim(sim.groups)

sim.groups <- logNormCounts(sim.groups)
sim.groups <- runPCA(sim.groups)
plotPCA(sim.groups, colour_by = "Group")



######################################################################
# Simulating 5 equal groups
sim.groups <- splatSimulate(params=params,
                            group.prob = c(0.2, 0.2, 0.2, 0.2, 0.2), method = "groups",
                            verbose = FALSE)
dim(sim.groups)

sim.groups <- logNormCounts(sim.groups)
sim.groups <- runPCA(sim.groups)
plotPCA(sim.groups, colour_by = "Group")