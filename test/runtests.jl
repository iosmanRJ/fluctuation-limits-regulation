using EvoSim
using Test, Random, Statistics

@testset "EvoSim" begin

# =============================================================================
@testset "genome: development" begin
    rng = Xoshiro(1)
    W = zeros(Float32, NG, NG, 4); bias = zeros(Float32, NG, 4)
    for i in 1:4
        seed_genome!(W, bias, i, rng)
    end
    phen = zeros(Float32, NPHEN); e = zeros(Float32, NG); buf = zeros(Float32, NG)
    env = zeros(Float32, NRES)

    develop!(phen, e, buf, W, bias, 1, env)
    @test all(0.0f0 .<= phen .<= 1.0f0)
    @test all(isfinite, phen)

    # deterministic: same genome + same environment -> same phenotype
    p2 = copy(phen)
    develop!(phen, e, buf, W, bias, 1, env)
    @test phen == p2

    # excretion profile is a proportion and must sum to one
    @test sum(phen[P_EXCRETE]) ≈ 1.0f0 atol=1e-4

    # plasticity: the same genome in a different environment develops differently
    env2 = fill(2.0f0, NRES)
    phen_b = zeros(Float32, NPHEN)
    develop!(phen_b, e, buf, W, bias, 1, env2)
    @test phen_b != phen

    # distinct genomes give distinct phenotypes
    phen_c = zeros(Float32, NPHEN)
    develop!(phen_c, e, buf, W, bias, 2, env)
    @test phen_c != phen

    # development must not allocate: it runs on every birth
    develop!(phen, e, buf, W, bias, 1, env)
    @test @allocated(develop!(phen, e, buf, W, bias, 1, env)) == 0
end

# =============================================================================
@testset "genome: mutation and inheritance" begin
    rng = Xoshiro(7)
    p = Params(mut_rate = 0.01f0, p_silence = 0.0f0, p_duplicate = 0.0f0)
    W = zeros(Float32, NG, NG, 2); bias = zeros(Float32, NG, 2)
    seed_genome!(W, bias, 1, rng)

    nmut, _ = inherit!(W, bias, 2, 1, rng, p, 0.0f0)
    @test nmut >= 0
    # child is a copy of the parent apart from the mutated sites
    diffs = count(!=(0.0f0), W[:, :, 2] .- W[:, :, 1])
    @test diffs <= nmut

    # perfect repair suppresses the mutation rate relative to none
    p2 = Params(mut_rate = 0.02f0, p_silence = 0.0f0, p_duplicate = 0.0f0,
                repair_strength = 0.85f0)
    W2 = zeros(Float32, NG, NG, 2); b2 = zeros(Float32, NG, 2)
    seed_genome!(W2, b2, 1, rng)
    lo = sum(first(inherit!(W2, b2, 2, 1, rng, p2, 1.0f0)) for _ in 1:400)
    hi = sum(first(inherit!(W2, b2, 2, 1, rng, p2, 0.0f0)) for _ in 1:400)
    @test lo < hi

    # gene silencing zeroes every outgoing weight of exactly one node
    p3 = Params(mut_rate = 0.0f0, p_silence = 1.0f0, p_duplicate = 0.0f0)
    W3 = ones(Float32, NG, NG, 1); b3 = zeros(Float32, NG, 1)
    mutate!(W3, b3, 1, rng, p3, 0.0f0)
    silenced = [k for k in 1:NG if all(W3[k, :, 1] .== 0.0f0)]
    @test length(silenced) == 1
end

# =============================================================================
@testset "population: slot allocator" begin
    p = Params(capacity = 8)
    pop = Population(p)
    @test pop.n == 0

    slots = [claim!(pop) for _ in 1:8]
    @test sort(slots) == collect(1:8)
    @test pop.n == 8
    @test claim!(pop) == 0                 # dish is full

    kill!(pop, 3); kill!(pop, 5)
    @test pop.n == 6
    @test !pop.alive[3] && !pop.alive[5]

    reused = [claim!(pop), claim!(pop)]
    @test sort(reused) == [3, 5]           # freed slots come back
    @test pop.n == 8
    @test claim!(pop) == 0

    # allocation must be O(1), never a scan of capacity
    p2 = Params(capacity = 100_000)
    pop2 = Population(p2)
    for _ in 1:1000; claim!(pop2); end
    @test @allocated(claim!(pop2)) == 0
end

# =============================================================================
@testset "environment: diffusion" begin
    # With the cap out of reach, no-flux diffusion must conserve mass exactly.
    p = Params(grid = 32, res_cap = 1.0f9, diff_rate = 0.2f0)
    env = Environment(p)
    rng = Xoshiro(3)
    for r in 1:NRES, j in 1:32, i in 1:32
        env.res[r, i, j] = rand(rng, Float32) * 10
    end
    before = [total_resource(env, r) for r in 1:NRES]
    for _ in 1:50
        diffuse!(env, p)
    end
    after = [total_resource(env, r) for r in 1:NRES]
    for r in 1:NRES
        @test after[r] ≈ before[r] rtol=1e-3      # nothing leaks out of the dish
    end

    # diffusion must flatten a point source, not amplify it
    env2 = Environment(p)
    fill!(env2.res, 0.0f0)
    env2.res[1, 16, 16] = 100.0f0
    for _ in 1:30
        diffuse!(env2, p)
    end
    @test env2.res[1, 16, 16] < 100.0f0
    @test env2.res[1, 18, 16] > 0.0f0
    @test maximum(env2.res) <= 100.0f0

    # a corner source must not produce negative concentrations
    env3 = Environment(p)
    fill!(env3.res, 0.0f0)
    env3.res[1, 1, 1] = 50.0f0
    for _ in 1:30
        diffuse!(env3, p)
    end
    @test minimum(env3.res) >= -1.0f-5
end

# =============================================================================
@testset "environment: replenish and seasons" begin
    p = Params(grid = 16, inflow = 0.1f0, decay = 0.0f0, res_cap = 2.0f0,
               season_amp = 0.0f0)
    env = Environment(p)
    fill!(env.res, 0.0f0)
    replenish!(env, p, 1)
    @test all(env.res[1, :, :] .≈ 0.1f0)
    @test all(env.res[2, :, :] .== 0.0f0)   # waste channels get no abiotic inflow

    for _ in 1:100
        replenish!(env, p, 1)
    end
    @test maximum(env.res) <= p.res_cap + 1.0f-5   # cap is respected

    # famine suppresses inflow, then lifts
    p2 = Params(grid = 8, inflow = 1.0f0, decay = 0.0f0, res_cap = 1.0f9,
                season_amp = 0.0f0, famine_factor = 0.0f0)
    env2 = Environment(p2); fill!(env2.res, 0.0f0)
    env2.famine = 3
    for _ in 1:3
        replenish!(env2, p2, 1)
    end
    @test all(env2.res[1, :, :] .== 0.0f0)
    @test env2.famine == 0
    replenish!(env2, p2, 1)
    @test all(env2.res[1, :, :] .≈ 1.0f0)          # inflow resumes
end

# =============================================================================
@testset "sim: stepping" begin
    p = Params(grid = 48, capacity = 3000, shock_interval = 0, seed = 11)
    sim = Sim(p)
    @test sim.pop.n > 0

    n0 = sim.pop.n
    for _ in 1:200
        step!(sim)
    end
    @test sim.gen == 200
    @test sim.pop.n <= p.capacity
    @test all(isfinite, sim.env.res)
    @test minimum(sim.env.res) >= -1.0f-5

    # every live organism has a valid position and a developed phenotype
    for i in 1:high_water(sim.pop)
        sim.pop.alive[i] || continue
        @test 1 <= sim.pop.x[i] <= p.grid
        @test 1 <= sim.pop.y[i] <= p.grid
        @test all(isfinite, sim.pop.phen[:, i])
    end

    # same seed must reproduce exactly — otherwise no result is ever checkable
    s1 = Sim(Params(grid = 32, capacity = 1000, seed = 99))
    s2 = Sim(Params(grid = 32, capacity = 1000, seed = 99))
    for _ in 1:120
        step!(s1); step!(s2)
    end
    @test s1.pop.n == s2.pop.n
    @test s1.env.res == s2.env.res
    @test s1.pop.phen == s2.pop.phen

    # different seeds must diverge
    s3 = Sim(Params(grid = 32, capacity = 1000, seed = 100))
    for _ in 1:120
        step!(s3)
    end
    @test s3.env.res != s1.env.res
end

@testset "sim: the hot loop does not allocate" begin
    p = Params(grid = 48, capacity = 3000, shock_interval = 0, seed = 5)
    sim = Sim(p)
    for _ in 1:80                    # let the birth buffer reach steady size
        step!(sim)
    end
    a = @allocated step!(sim)
    # a handful of small allocations is tolerable; growth per generation is not
    @test a < 8_000
end

@testset "sim: extinction is reachable and handled" begin
    # brutal upkeep with no inflow: everything must starve, without erroring
    p = Params(grid = 32, capacity = 2000, inflow = 0.0f0, maint_base = 5.0f0,
               shock_interval = 0, seed = 3)
    sim = Sim(p)
    for _ in 1:50
        step!(sim)
    end
    @test sim.pop.n == 0
    step!(sim)                       # stepping an empty dish must be safe
    @test sim.pop.n == 0
end

@testset "sim: shocks kill" begin
    p = Params(grid = 48, capacity = 3000, shock_interval = 1,
               shock_kill = 0.9f0, seed = 21)
    sim = Sim(p)
    n0 = sim.pop.n
    step!(sim)
    @test sim.pop.n < n0             # a shock every generation must bite
    @test !isempty(sim.shock_log)
end

# =============================================================================
@testset "stats: schedule" begin
    s = snapshot_schedule(1_000_000, 12)
    @test first(s) == 1
    @test last(s) == 1_000_000
    @test issorted(s)
    @test allunique(s)
    @test length(s) < 200                      # log-spaced, not linear
    # roughly log-uniform: as many samples below 1000 as above
    @test count(<(1000), s) > 20

    s2 = snapshot_schedule(10, 4)
    @test first(s2) == 1 && last(s2) == 10 && issorted(s2)
end

@testset "stats: log roundtrip" begin
    dir = mktempdir()
    p = Params(grid = 32, capacity = 1000, shock_interval = 0, seed = 4)
    sim = Sim(p)
    log = StatLog(dir; chunk = 16)
    for _ in 1:40
        step!(sim)
        record!(log, sim)
    end
    close(log)

    M = read_stats(dir)
    @test size(M, 1) == NSTAT
    @test size(M, 2) == 40
    @test length(STAT_NAMES) == NSTAT
    @test M[stat_index("gen"), :] == Float32.(1:40)
    @test all(M[stat_index("pop"), :] .>= 0)
    @test all(isfinite, M)
    @test all(M[stat_index("clade_shannon"), :] .>= 0)    # entropy is non-negative
end

@testset "stats: snapshot roundtrip" begin
    dir = mktempdir()
    p = Params(grid = 24, capacity = 800, shock_interval = 0, seed = 6)
    sim = Sim(p)
    for _ in 1:30
        step!(sim)
    end
    write_snapshot(sim, dir; nsample = 50)
    f = only(filter(endswith(".bin"), readdir(dir; join = true)))
    snap = read_snapshot(f)

    @test snap.gen == sim.gen
    @test snap.G == 24
    @test size(snap.res) == (NRES, 24, 24)
    @test snap.res == sim.env.res
    @test size(snap.phen, 1) == NPHEN
    @test length(snap.x) <= 50
    @test all(1 .<= snap.x .<= 24)
    @test all(1 .<= snap.y .<= 24)
    @test !isempty(snap.W)
    @test size(first(snap.W)) == (NG, NG)
end

# =============================================================================
@testset "run_sim: end to end" begin
    dir = mktempdir()
    r = run_sim(gens = 300, outdir = dir, verbose = false,
                p = Params(grid = 32, capacity = 2000, shock_interval = 0, seed = 8))
    @test r.gens_run == 300
    @test isfile(joinpath(dir, "stats.bin"))
    @test isfile(joinpath(dir, "meta.txt"))
    @test isfile(joinpath(dir, "clades.tsv"))
    @test isfile(joinpath(dir, "shocks.tsv"))
    @test !isempty(filter(startswith("snap_"), readdir(dir)))

    M = read_stats(dir)
    @test size(M, 2) == 300
    @test M[stat_index("gen"), end] == 300f0
end

end # EvoSim testset
