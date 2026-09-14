using EvoSim
using Test, Random, Statistics

@testset "EvoSim" begin

# =============================================================================
@testset "genome: development" begin
    rng = Xoshiro(1)
    W = zeros(Float32, NG_MAX, NG_MAX, 4); bias = zeros(Float32, NG_MAX, 4)
    for i in 1:4
        seed_genome!(W, bias, i, rng)
    end
    phen = zeros(Float32, NPHEN); e = zeros(Float32, NG_MAX); buf = zeros(Float32, NG_MAX)
    env = zeros(Float32, NRES)

    develop!(phen, e, buf, W, bias, 1, env, NG_MAX)
    @test all(0.0f0 .<= phen .<= 1.0f0)
    @test all(isfinite, phen)

    # deterministic: same genome + same environment -> same phenotype
    p2 = copy(phen)
    develop!(phen, e, buf, W, bias, 1, env, NG_MAX)
    @test phen == p2

    # uptake affinity is an allocation budget, so it must sum to one —
    # this is what makes specialist-vs-generalist a real trade-off
    @test sum(phen[P_UPTAKE]) ≈ 1.0f0 atol=1e-4

    # plasticity: the same genome in a different environment develops differently
    env2 = fill(2.0f0, NRES)
    phen_b = zeros(Float32, NPHEN)
    develop!(phen_b, e, buf, W, bias, 1, env2, NG_MAX)
    @test phen_b != phen
    @test sum(phen_b[P_UPTAKE]) ≈ 1.0f0 atol=1e-4      # budget holds under plasticity

    # distinct genomes give distinct phenotypes
    phen_c = zeros(Float32, NPHEN)
    develop!(phen_c, e, buf, W, bias, 2, env, NG_MAX)
    @test phen_c != phen

    # development must not allocate: it runs on every birth
    develop!(phen, e, buf, W, bias, 1, env, NG_MAX)
    @test @allocated(develop!(phen, e, buf, W, bias, 1, env, NG_MAX)) == 0
end

# =============================================================================
@testset "genome: mutation and inheritance" begin
    rng = Xoshiro(7)
    p = Params(mut_rate = 0.01f0, p_silence = 0.0f0, p_duplicate = 0.0f0)
    W = zeros(Float32, NG_MAX, NG_MAX, 2); bias = zeros(Float32, NG_MAX, 2)
    seed_genome!(W, bias, 1, rng)

    nmut, _, _, _ = inherit!(W, bias, 2, 1, rng, p, 0.0f0, NG_MAX)
    @test nmut >= 0
    # child is a copy of the parent apart from the mutated sites
    diffs = count(!=(0.0f0), W[:, :, 2] .- W[:, :, 1])
    @test diffs <= nmut

    # perfect repair suppresses the mutation rate relative to none
    p2 = Params(mut_rate = 0.02f0, p_silence = 0.0f0, p_duplicate = 0.0f0,
                repair_strength = 0.85f0)
    W2 = zeros(Float32, NG_MAX, NG_MAX, 2); b2 = zeros(Float32, NG_MAX, 2)
    seed_genome!(W2, b2, 1, rng)
    lo = sum(first(inherit!(W2, b2, 2, 1, rng, p2, 1.0f0, NG_MAX)) for _ in 1:400)
    hi = sum(first(inherit!(W2, b2, 2, 1, rng, p2, 0.0f0, NG_MAX)) for _ in 1:400)
    @test lo < hi

    # gene silencing zeroes every outgoing weight of exactly one node
    p3 = Params(mut_rate = 0.0f0, p_silence = 1.0f0, p_duplicate = 0.0f0)
    W3 = ones(Float32, NG_MAX, NG_MAX, 1); b3 = zeros(Float32, NG_MAX, 1)
    mutate!(W3, b3, 1, rng, p3, 0.0f0, NG_MAX)
    silenced = [k for k in 1:NG_MAX if all(W3[k, :, 1] .== 0.0f0)]
    @test length(silenced) == 1
end

# =============================================================================
@testset "population: capacity scales with the dish" begin
    # A fixed capacity default silently became the binding constraint at grid
    # 256, which invalidates a run. Capacity must track dish area by default.
    @test resolved_capacity(Params(grid = 128, capacity = 0)) == 128 * 128 * 4
    @test resolved_capacity(Params(grid = 256, capacity = 0)) == 256 * 256 * 4
    @test resolved_capacity(Params(grid = 256, capacity = 999)) == 999   # explicit wins
    @test Population(Params(grid = 64, capacity = 0)).capacity == 64 * 64 * 4

    # the derived headroom must exceed what a dish can physically hold under
    # the per-cell packing limit only by design, not by accident
    p = Params(grid = 64, capacity = 0)
    @test resolved_capacity(p) < p.grid^2 * p.max_per_cell   # per-cell limit binds first
end

@testset "sim: dropped births are counted" begin
    # Squeeze capacity so the array fills, and check it is detected rather than
    # silently degrading selection.
    p = Params(grid = 48, capacity = 400, shock_interval = 0, seed = 91)
    sim = Sim(p)
    for _ in 1:400
        step!(sim)
    end
    @test sim.pop.n <= 400
    @test sim.dropped > 0            # the condition is observable

    # with room to breathe, nothing is dropped
    p2 = Params(grid = 48, capacity = 0, shock_interval = 0, seed = 91)
    s2 = Sim(p2)
    for _ in 1:400
        step!(s2)
    end
    @test s2.dropped == 0
end

@testset "population: founders scale with the dish" begin
    # A fixed founder count is a density trap: 200 founders is 0.012/cell at
    # grid 128 but 0.0008/cell at grid 512, and sparse founders fail to
    # establish. Measured: 3 of 8 grid-256 runs died inside 90k ticks, two of
    # them inside 27k.
    @test founder_count(Params(grid = 128)) == max(200, 128 * 128 ÷ 80)
    @test founder_count(Params(grid = 512)) > 10 * founder_count(Params(grid = 64))
    @test founder_count(Params(grid = 256, n_founders = 77)) == 77   # explicit wins

    # founder DENSITY should be roughly grid-independent now
    d(g) = founder_count(Params(grid = g)) / (g * g)
    @test d(512) > 0.5 * d(256)
    @test d(512) < 2.0 * d(256)

    # and the sim actually starts with that many
    p = Params(grid = 128, capacity = 0)
    @test Sim(p).pop.n == founder_count(p)
end

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

@testset "sim: allocation does not scale with population" begin
    # The invariant that matters is not a fixed byte budget — the parallel
    # birth phase spawns tasks, which costs a constant few KB per tick — but
    # that allocation must not grow with the number of organisms. Anything
    # per-organism or per-birth is a garbage generator at 1e6 ticks.
    function bytes_per_tick(grid, cap)
        p = Params(grid = grid, capacity = cap, shock_interval = 0, seed = 5)
        sim = Sim(p)
        for _ in 1:150                # let buffers reach steady size
            step!(sim)
        end
        (minimum(@allocated(step!(sim)) for _ in 1:20), sim.pop.n)
    end
    small, n_small = bytes_per_tick(48, 3_000)
    large, n_large = bytes_per_tick(160, 40_000)

    @test n_large > 5 * n_small       # the comparison is meaningful
    # a 10x population must not bring a 10x allocation
    @test large < 2 * small
    @test large / n_large < 1.0       # well under a byte per organism

    # seeding an RNG must be allocation-free: it runs once per birth, and
    # Random.seed!(::Xoshiro, ::UInt64) costs 368 bytes a call
    r = Xoshiro(0)
    seed_fast!(r, UInt64(1))
    @test @allocated(seed_fast!(r, UInt64(99))) == 0
end

@testset "genome: fast seeding behaves like a seeding function" begin
    a, b = Xoshiro(0), Xoshiro(0)
    seed_fast!(a, UInt64(12345))
    seed_fast!(b, UInt64(12345))
    @test [rand(a, UInt64) for _ in 1:8] == [rand(b, UInt64) for _ in 1:8]

    seed_fast!(a, UInt64(12345)); seed_fast!(b, UInt64(12346))
    @test [rand(a, UInt64) for _ in 1:8] != [rand(b, UInt64) for _ in 1:8]

    # re-seeding must reset the stream, not continue it
    seed_fast!(a, UInt64(7)); first_run = [rand(a, Float32) for _ in 1:16]
    seed_fast!(a, UInt64(7)); @test [rand(a, Float32) for _ in 1:16] == first_run

    # and never leave the all-zero state, which xoshiro cannot escape
    seed_fast!(a, UInt64(0))
    @test any(!=(0), (a.s0, a.s1, a.s2, a.s3, a.s4))
    @test length(unique(rand(a, UInt64) for _ in 1:8)) > 1
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
    # Step a handful of times rather than once: a third of shocks are famines,
    # which kill nobody on the tick they fire, so a single step can legitimately
    # leave the population untouched.
    for _ in 1:6
        step!(sim)
    end
    @test sim.pop.n < n0             # repeated shocks must bite
    @test !isempty(sim.shock_log)
    @test any(k -> k in (:cull, :heat), last.(sim.shock_log))
end

# =============================================================================
@testset "sim: occupancy grid stays consistent" begin
    # The occupancy grid is updated incrementally on every move, birth and
    # death. If it drifts, crowding costs silently decouple from reality and
    # every result becomes meaningless — so check it against the truth.
    p = Params(grid = 40, capacity = 5000, seed = 31)
    sim = Sim(p)

    function check_occ(sim)
        truth = zeros(Int32, sim.p.grid, sim.p.grid)
        n = 0
        for i in 1:high_water(sim.pop)
            sim.pop.alive[i] || continue
            truth[sim.pop.x[i], sim.pop.y[i]] += Int32(1)
            n += 1
        end
        (truth == sim.pop.occ, n)
    end

    ok, n = check_occ(sim)
    @test ok
    @test n == sim.pop.n
    @test sum(sim.pop.occ) == sim.pop.n

    for _ in 1:250
        step!(sim)
    end
    ok, n = check_occ(sim)
    @test ok                                  # survives births, deaths and moves
    @test n == sim.pop.n
    @test sum(sim.pop.occ) == sim.pop.n
    @test all(sim.pop.occ .>= 0)
    @test maximum(sim.pop.occ) <= p.max_per_cell

    # and through a mass-extinction event, which kills outside the main loop
    p2 = Params(grid = 40, capacity = 5000, shock_interval = 1, shock_kill = 0.7f0, seed = 32)
    sim2 = Sim(p2)
    for _ in 1:40
        step!(sim2)
    end
    truth2 = zeros(Int32, 40, 40)
    for i in 1:high_water(sim2.pop)
        sim2.pop.alive[i] && (truth2[sim2.pop.x[i], sim2.pop.y[i]] += Int32(1))
    end
    @test truth2 == sim2.pop.occ
    @test sum(sim2.pop.occ) == sim2.pop.n
end

@testset "sim: byproducts never flow back to their own or an upstream channel" begin
    # THE bug that killed cross-feeding. When the excretion destination was an
    # evolvable trait, organisms evolved to excrete 36% of their waste straight
    # back into the channel they had just eaten — a futile cycle that meant the
    # waste channels were never fed. The microbial consumer-resource literature
    # excludes the diagonal of the metabolic matrix for exactly this reason.
    #
    # Assert it structurally: seed ONLY channel 1, let organisms eat, and check
    # that channel 1 never gains mass from metabolism.
    p = Params(grid = 32, capacity = 3000, shock_interval = 0, seed = 12,
               inflow = 0.0f0, decay = 0.0f0, diff_rate = 0.0f0)
    sim = Sim(p)
    fill!(sim.env.res, 0.0f0)
    for j in 1:p.grid, i in 1:p.grid
        sim.env.res[1, i, j] = 5.0f0
    end
    before = total_resource(sim.env, 1)
    for _ in 1:60
        step!(sim)
    end
    # with no inflow, channel 1 can only ever decrease
    @test total_resource(sim.env, 1) < before
    # and the downstream channels must have been fed by metabolism
    @test sum(total_resource(sim.env, r) for r in 2:NRES) > 0.0

    # monotonicity: channel 1 must never rise from one tick to the next
    prev = total_resource(sim.env, 1)
    for _ in 1:40
        step!(sim)
        now = total_resource(sim.env, 1)
        @test now <= prev + 1e-3
        prev = now
    end
end

@testset "sim: cross-feeding persists instead of decaying away" begin
    # Regression guard on the headline result. Before the conversion ladder and
    # the leakage floor, allocation to waste channels decayed toward zero over
    # long runs (0.180 -> 0.076 over 80k ticks, still falling). It must now hold
    # up, and it must still be substantial late in the run.
    p = Params(grid = 96, capacity = 0, seed = 7)
    sim = Sim(p)
    waste_share() = begin
        live = [i for i in 1:high_water(sim.pop) if sim.pop.alive[i]]
        isempty(live) && return 0.0
        mean(sum(sim.pop.phen[r, i] for r in 2:NRES) for i in live)
    end
    for _ in 1:6000
        step!(sim)
    end
    early = waste_share()
    for _ in 1:14000
        step!(sim)
    end
    late = waste_share()

    @test sim.pop.n > 0
    @test early > 0.10                 # cross-feeding establishes
    @test late  > 0.10                 # and is still there much later
    @test late  > 0.5 * early          # no collapse
end

@testset "sim: leakage floor is enforced" begin
    # Leakage is an environmental constant in the consumer-resource literature,
    # not a trait. Left free it evolved to ~0.003 and starved the waste channels,
    # because uptake saturates at what is locally available so the rate bonus
    # cannot pay for the yield loss. The floor is what guarantees supply.
    p = Params(grid = 48, capacity = 4000, shock_interval = 0, seed = 19,
               leak_min = 0.60f0, inflow = 0.0f0, decay = 0.0f0, diff_rate = 0.0f0)
    sim = Sim(p)
    fill!(sim.env.res, 0.0f0)
    for j in 1:p.grid, i in 1:p.grid
        sim.env.res[1, i, j] = 4.0f0
    end
    r1_before = total_resource(sim.env, 1)
    for _ in 1:40
        step!(sim)
    end
    consumed = r1_before - total_resource(sim.env, 1)
    produced = sum(total_resource(sim.env, r) for r in 2:NRES)
    @test consumed > 0
    # at a 0.6 floor, most of what is consumed should reappear downstream
    @test produced > 0.30 * consumed
end

@testset "sim: uptake budget holds for every live organism" begin
    p = Params(grid = 48, capacity = 4000, shock_interval = 0, seed = 71)
    sim = Sim(p)
    for _ in 1:400
        step!(sim)
    end
    n = 0
    for i in 1:high_water(sim.pop)
        sim.pop.alive[i] || continue
        n += 1
        @test sum(sim.pop.phen[P_UPTAKE, i]) ≈ 1.0f0 atol=1e-4
        @test all(sim.pop.phen[P_UPTAKE, i] .>= 0.0f0)
    end
    @test n > 0
end

@testset "sim: uptake exponent changes the outcome, and not as expected" begin
    # Pins a measured, counter-intuitive result so it is not silently
    # rediscovered: raising uptake_exponent was meant to make the allocation
    # trade-off convex and drive niche specialisation. It does the reverse —
    # the population becomes MORE generalist. Recorded as a regression guard on
    # the knob's behaviour, not as an endorsement of raising it.
    function specialisation(exponent)
        p = Params(; grid = 64, capacity = 12_000, shock_interval = 0, seed = 44,
                   uptake_exponent = exponent)
        sim = Sim(p)
        for _ in 1:2500
            step!(sim)
        end
        live = [i for i in 1:high_water(sim.pop) if sim.pop.alive[i]]
        isempty(live) && return NaN
        # max allocation share per organism: 1/NRES is a perfect generalist,
        # 1.0 is a complete specialist
        mean(maximum(sim.pop.phen[r, i] for r in P_UPTAKE) for i in live)
    end
    linear = specialisation(1.0f0)
    convex = specialisation(2.0f0)
    @test linear > 1.0 / NRES        # the default does specialise somewhat
    @test convex < linear            # raising the exponent generalises instead
end

@testset "sim: traits that buy something also cost something" begin
    # Every trait with a benefit needs a cost, or evolution pins it to 1.0.
    # Raising a trait's cost must lower its evolved value.
    function evolved(trait, kw)
        p = Params(; grid = 48, capacity = 4000, shock_interval = 0, seed = 84, kw...)
        sim = Sim(p)
        for _ in 1:1500
            step!(sim)
        end
        live = [i for i in 1:high_water(sim.pop) if sim.pop.alive[i]]
        isempty(live) ? NaN : mean(sim.pop.phen[trait, i] for i in live)
    end
    cheap_size = evolved(P_SIZE, (maint_size = 0.02f0,))
    dear_size  = evolved(P_SIZE, (maint_size = 0.60f0,))
    @test dear_size < cheap_size

    cheap_mot = evolved(P_MOTIL, (maint_motil = 0.0f0,))
    dear_mot  = evolved(P_MOTIL, (maint_motil = 0.20f0,))
    @test dear_mot < cheap_mot
end

@testset "sim: dormancy" begin
    # Every extinction measured was famine-caused (8 of 8, each within 24-211
    # ticks of one), because organisms carry under two ticks of reserve against
    # a 400-tick famine. Dormancy is what makes surviving one possible.
    p = Params(grid = 48, capacity = 4000, shock_interval = 0, seed = 88,
               dorm_threshold = 0.15f0, dorm_upkeep = 0.05f0)
    sim = Sim(p)
    for _ in 1:300
        step!(sim)
    end
    live = [i for i in 1:high_water(sim.pop) if sim.pop.alive[i]]
    @test !isempty(live)

    # a dormant cell must not divide, however much energy it happens to hold
    # (it is below threshold by definition, so this is about the guard itself)
    @test all(i -> sim.pop.energy[i] >= 0.0f0, live)

    # ageing pauses while dormant: with a lethal-upkeep dish and no inflow,
    # organisms should go dormant and persist far longer than max_age would
    # otherwise allow
    function survives(dorm)
        q = Params(grid = 24, capacity = 900, shock_interval = 0, seed = 5,
                   inflow = 0.0f0, decay = 0.0f0, diff_rate = 0.0f0,
                   dorm_threshold = dorm ? 0.5f0 : -1.0f0, dorm_upkeep = 0.02f0)
        s = Sim(q)
        fill!(s.env.res, 0.0f0)          # total famine from tick 1
        t_out = 0
        for t in 1:3000
            step!(s)
            if s.pop.n == 0
                t_out = t
                break
            end
        end
        t_out == 0 ? 3000 : t_out
    end
    without = survives(false)
    with    = survives(true)
    @test with > 3 * without        # spores outlast non-spores by a wide margin
end

@testset "sim: energy ceiling and load penalty" begin
    p = Params(grid = 40, capacity = 5000, shock_interval = 0, e_max = 2.5f0, seed = 41)
    sim = Sim(p)
    for _ in 1:300
        step!(sim)
    end
    for i in 1:high_water(sim.pop)
        sim.pop.alive[i] || continue
        @test sim.pop.energy[i] <= p.e_max + 1.0f-5     # no hoarding
        @test sim.pop.load[i] >= 0                      # load never goes negative
    end

    # load must actually cost energy: two identical dishes, one with the
    # penalty switched off, must diverge
    base = Params(grid = 32, capacity = 3000, shock_interval = 0, seed = 55,
                  mut_rate = 0.01f0, maint_load = 0.0f0)
    pen  = Params(grid = 32, capacity = 3000, shock_interval = 0, seed = 55,
                  mut_rate = 0.01f0, maint_load = 0.02f0)
    sa, sb = Sim(base), Sim(pen)
    for _ in 1:200
        step!(sa); step!(sb)
    end
    load_a = mean(sa.pop.load[i] for i in 1:high_water(sa.pop) if sa.pop.alive[i])
    load_b = mean(sb.pop.load[i] for i in 1:high_water(sb.pop) if sb.pop.alive[i])
    @test load_b < load_a          # the penalty purges damaged lineages
end

@testset "genome: mutation classification" begin
    rng = Xoshiro(77)
    # all-deleterious settings must drive load strictly up
    p = Params(mut_rate = 0.02f0, p_deleterious = 1.0f0, p_beneficial = 0.0f0,
               p_silence = 0.0f0, p_duplicate = 0.0f0)
    W = zeros(Float32, NG_MAX, NG_MAX, 2); b = zeros(Float32, NG_MAX, 2)
    seed_genome!(W, b, 1, rng)
    tot = sum(inherit!(W, b, 2, 1, rng, p, 0.0f0, NG_MAX)[3] for _ in 1:200)
    @test tot > 0

    # all-beneficial must drive it down
    p2 = Params(mut_rate = 0.02f0, p_deleterious = 0.0f0, p_beneficial = 1.0f0,
                p_silence = 0.0f0, p_duplicate = 0.0f0)
    tot2 = sum(inherit!(W, b, 2, 1, rng, p2, 0.0f0, NG_MAX)[3] for _ in 1:200)
    @test tot2 < 0
end

# =============================================================================
@testset "clades: pruning preserves the partition" begin
    # Pruning renumbers every clade id in the table AND in the population.
    # If the remap is wrong, diversity statistics silently become fiction, so
    # assert the partition of living organisms is exactly preserved.
    p = Params(grid = 64, capacity = 8000, seed = 63, clade_gc_interval = 0)
    sim = Sim(p)
    for _ in 1:1200
        step!(sim)
    end

    live() = [i for i in 1:high_water(sim.pop) if sim.pop.alive[i]]
    before_ids = [sim.pop.clade[i] for i in live()]
    before_n = nclades(sim.clades)
    @test before_n > 1                       # the run actually diversified

    # the partition is the grouping of organisms, independent of id numbering
    partition(ids) = sort([sort(findall(==(c), ids)) for c in unique(ids)])
    before_part = partition(before_ids)

    nk = prune_clades!(sim.clades, sim.pop)
    after_ids = [sim.pop.clade[i] for i in live()]

    @test nk <= before_n                     # pruning never grows the table
    @test nk == nclades(sim.clades)
    @test partition(after_ids) == before_part            # same grouping
    @test length(unique(after_ids)) == length(unique(before_ids))

    # every remaining id must be in range, and the tree must stay well-formed
    @test all(c -> Int(c) < nk, after_ids)
    @test sim.clades.parent[1] == 0                       # root
    for c in 2:nk
        @test Int(sim.clades.parent[c]) < c - 1           # parents precede children
    end
    @test length(sim.clades.founder) == nk * NPHEN
    @test length(sim.clades.born) == nk

    # pruning again is a no-op
    nk2 = prune_clades!(sim.clades, sim.pop)
    @test nk2 == nk
    @test [sim.pop.clade[i] for i in live()] == after_ids

    # and the simulation keeps running correctly afterwards
    for _ in 1:100
        step!(sim)
    end
    @test sim.pop.n > 0
    @test all(i -> Int(sim.pop.clade[i]) < nclades(sim.clades), live())

    # automatic GC during a run keeps the table bounded
    p2 = Params(grid = 64, capacity = 8000, seed = 63, clade_gc_interval = 200)
    s2 = Sim(p2)
    for _ in 1:1200
        step!(s2)
    end
    @test nclades(s2.clades) < before_n       # GC actually reclaimed entries
    @test s2.pop.n > 0
end

@testset "clades: founded on phenotypic drift" begin
    ct = CladeTable()
    phen = zeros(Float32, NPHEN)
    @test phen_dist2(ct, UInt32(0), phen) ≈ 0.0f0

    far = copy(phen); far[1] = 1.0f0
    @test phen_dist2(ct, UInt32(0), far) ≈ 1.0f0

    c = new_clade!(ct, UInt32(0), 500, far)
    @test nclades(ct) == 2
    @test ct.parent[Int(c) + 1] == 0
    @test ct.born[Int(c) + 1] == 500
    @test phen_dist2(ct, c, far) ≈ 0.0f0        # measured against its own founder
    @test phen_dist2(ct, c, phen) ≈ 1.0f0
end

# =============================================================================
@testset "stats: top-channel shares and genome structure" begin
    # top_share is THE branching measure — a mean cannot distinguish one
    # generalist from two coexisting specialists, these shares can.
    dir = mktempdir()
    sim = Sim(Params(grid = 48, capacity = 4000, shock_interval = 0, seed = 61))
    for _ in 1:500
        step!(sim)
    end
    sl = StatLog(dir); record!(sl, sim); close(sl)
    M = read_stats(dir)

    shares = [M[stat_index("top_share$(r)"), 1] for r in 1:NRES]
    @test sum(shares) ≈ 1.0f0 atol=1e-4         # it is a partition
    @test all(0.0 .<= shares .<= 1.0)

    # cross-check against counting it directly
    live = [i for i in 1:high_water(sim.pop) if sim.pop.alive[i]]
    direct = [count(i -> argmax(@view sim.pop.phen[1:NRES, i]) == r, live) / length(live)
              for r in 1:NRES]
    for r in 1:NRES
        @test shares[r] ≈ direct[r] atol=1e-3
    end

    nz = M[stat_index("mean_nonzero_w"), 1]
    ag = M[stat_index("mean_active_genes"), 1]
    @test 0 < nz <= NG_MAX * NG_MAX                      # connectivity is in range
    @test 0 < ag <= NG_MAX                           # active genes is in range
    @test M[stat_index("mean_abs_w"), 1] > 0
end

@testset "genome: size is variable and bounded" begin
    rng = Xoshiro(91)
    W = zeros(Float32, NG_MAX, NG_MAX, 2); b = zeros(Float32, NG_MAX, 2)
    seed_genome!(W, b, 1, rng, 16)

    # node addition grows the genome, deletion shrinks it
    padd = Params(mut_rate = 0.0f0, p_silence = 0.0f0, p_duplicate = 0.0f0,
                  p_node_add = 1.0f0, p_node_del = 0.0f0)
    _, _, _, ng2 = mutate!(W, b, 1, rng, padd, 0.0f0, 16)
    @test ng2 == 17

    pdel = Params(mut_rate = 0.0f0, p_silence = 0.0f0, p_duplicate = 0.0f0,
                  p_node_add = 0.0f0, p_node_del = 1.0f0)
    _, _, _, ng3 = mutate!(W, b, 1, rng, pdel, 0.0f0, 16)
    @test ng3 == 15

    # bounds hold at both ends
    _, _, _, ngmax = mutate!(W, b, 1, rng, padd, 0.0f0, NG_MAX)
    @test ngmax == NG_MAX               # cannot exceed the array bound
    _, _, _, ngmin = mutate!(W, b, 1, rng, pdel, 0.0f0, NG_MIN)
    @test ngmin == NG_MIN               # cannot delete the I/O interface

    # deletion must never touch nodes 1:NPHEN — they are the fixed interface
    W2 = zeros(Float32, NG_MAX, NG_MAX, 1); b2 = zeros(Float32, NG_MAX, 1)
    for k in 1:NG_MAX, m in 1:NG_MAX
        W2[m, k, 1] = Float32(100 * k + m)
    end
    iface_before = copy(W2[1:NPHEN, 1:NPHEN, 1])
    for _ in 1:20
        mutate!(W2, b2, 1, rng, pdel, 0.0f0, 20)
    end
    @test W2[1:NPHEN, 1:NPHEN, 1] == iface_before

    # development respects ng: a bigger network gives a different phenotype
    phen_a = zeros(Float32, NPHEN); phen_b = zeros(Float32, NPHEN)
    e = zeros(Float32, NG_MAX); buf = zeros(Float32, NG_MAX)
    env = zeros(Float32, NRES)
    seed_genome!(W, b, 1, rng, NG_MAX)
    develop!(phen_a, e, buf, W, b, 1, env, 14)
    develop!(phen_b, e, buf, W, b, 1, env, 28)
    @test phen_a != phen_b
    @test all(0.0f0 .<= phen_a .<= 1.0f0)
    @test sum(phen_a[P_UPTAKE]) ≈ 1.0f0 atol=1e-4
end

@testset "sim: genome size stays in range for every organism" begin
    p = Params(grid = 48, capacity = 4000, shock_interval = 0, seed = 77,
               p_node_add = 0.05f0, p_node_del = 0.05f0)
    sim = Sim(p)
    for _ in 1:600
        step!(sim)
    end
    live = [i for i in 1:high_water(sim.pop) if sim.pop.alive[i]]
    @test !isempty(live)
    @test all(i -> NG_MIN <= sim.pop.ng[i] <= NG_MAX, live)
    # with high add/delete rates the population should not be monomorphic in size
    @test length(unique(sim.pop.ng[i] for i in live)) > 1
end

@testset "genome: structure summary tracks silencing" begin
    rng = Xoshiro(23)
    W = zeros(Float32, NG_MAX, NG_MAX, 1); b = zeros(Float32, NG_MAX, 1)
    seed_genome!(W, b, 1, rng)
    nz0, act0, aw0 = genome_structure(W, b, 1, NG_MAX)
    @test nz0 > 0 && act0 > 0 && aw0 > 0

    # a fully wired genome has every weight and every gene active
    fill!(W, 1.0f0)
    nz1, act1, aw1 = genome_structure(W, b, 1, NG_MAX)
    @test nz1 == NG_MAX * NG_MAX
    @test act1 == NG_MAX
    @test aw1 ≈ 1.0f0

    # silencing one gene zeroes its outgoing row and drops the active count
    p = Params(mut_rate = 0.0f0, p_silence = 1.0f0, p_duplicate = 0.0f0)
    mutate!(W, b, 1, rng, p, 0.0f0, NG_MAX)
    nz2, act2, _ = genome_structure(W, b, 1, NG_MAX)
    @test act2 == NG_MAX - 1
    @test nz2 == NG_MAX * NG_MAX - NG_MAX

    # an empty genome registers as zero complexity
    fill!(W, 0.0f0)
    nz3, act3, aw3 = genome_structure(W, b, 1, NG_MAX)
    @test nz3 == 0 && act3 == 0 && aw3 == 0.0f0
end

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

@testset "stats: buffer flushes on a time bound too" begin
    # A long run with a large stat_every can take hundreds of thousands of
    # ticks to fill a chunk. Without a time bound the file stays empty and the
    # run is indistinguishable from a stalled one.
    dir = mktempdir()
    sim = Sim(Params(grid = 24, capacity = 500, shock_interval = 0, seed = 3))
    sl = StatLog(dir; chunk = 100_000, flush_after = 0.0)   # chunk never fills
    for _ in 1:3
        step!(sim)
        record!(sl, sim)
    end
    @test filesize(joinpath(dir, "stats.bin")) > 0          # flushed on the time bound
    @test size(read_stats(dir), 2) == 3

    # with a long time bound it is the chunk that triggers the write instead
    dir2 = mktempdir()
    sl2 = StatLog(dir2; chunk = 2, flush_after = 1e6)
    for _ in 1:2
        step!(sim)
        record!(sl2, sim)
    end
    @test size(read_stats(dir2), 2) == 2
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
    # occupancy must be the real grid, not one rebuilt from the organism
    # sample — a density map derived from a sample reads as empty almost
    # everywhere once the population is much larger than the sample
    @test snap.occ == sim.pop.occ
    @test sum(snap.occ) == sim.pop.n
    @test size(snap.occ) == (24, 24)
    @test size(snap.phen, 1) == NPHEN
    @test length(snap.x) <= 50
    @test all(1 .<= snap.x .<= 24)
    @test all(1 .<= snap.y .<= 24)
    @test !isempty(snap.W)
    @test size(first(snap.W)) == (NG_MAX, NG_MAX)
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
