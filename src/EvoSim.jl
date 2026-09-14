"""
    EvoSim

An agent-based evolution simulator: single-celled organisms with
gene-regulatory-network genomes in a 2D petri dish with diffusing,
cross-fed resource channels.

Quick start:

    using EvoSim
    r = run_sim(gens = 50_000, p = Params(grid = 256, capacity = 40_000))

Then chart it with `scripts/dashboard.jl \$(r.outdir)`.
"""
module EvoSim

using Random, Statistics, LinearAlgebra, Printf, Dates

include("genome.jl")       # defines NG, NRES, NPHEN — must come first
include("params.jl")
include("population.jl")
include("environment.jl")
include("sim.jl")
include("stats.jl")
include("run.jl")

export Params, Sim, Environment, Population, CladeTable,
       step!, run_sim, seed_population!,
       develop!, inherit!, mutate!, seed_genome!, seed_fast!, channel_cost,
       StatLog, record!, flush_stats!, read_stats, stat_index, STAT_NAMES,
       snapshot_schedule, write_snapshot, read_snapshot,
       replenish!, diffuse!, total_resource,
       claim!, kill!, place!, high_water, resolved_capacity, nclades, new_clade!, prune_clades!,
       phen_dist2,
       NG, NRES, NPHEN, NSTAT,
       P_UPTAKE, P_SIZE, P_MOTIL, P_DIVIDE, P_REPAIR, P_STRESS, P_LEAK

end # module
