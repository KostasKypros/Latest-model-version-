### root data per voxel prints
## meta data prints 

# VOXEL SIZE SENSITIVITY ANALYSIS 2

# test runs mono crops


# ═══════════════════════════════════════════════════════════════════

using VirtualPlantLab, ColorTypes
using Base.Threads: @threads
using Plots
import Random
using FastGaussQuadrature
using Distributions
Random.seed!(123456789)
import GLMakie
using SkyDomes
import Parameters: @with_kw
using Printf
using Statistics
using DataFrames      # ← NEW
using CSV             # ← NEW
Random.seed!(123456789)

# ═══════════════════════════════════════════════════════════════════
# SECTION: PLANT STRUCTS
# ═══════════════════════════════════════════════════════════════════

module TreeTypes
    using VirtualPlantLab
    using Distributions

    Base.@kwdef mutable struct Meristem <: VirtualPlantLab.Node
        age::Int64 = 0
        age_previous::Float64 = 0.0
        ageD::Float64 = 0.0
        ageDprevious::Float64 = 0.0
    end

    struct Bud <: VirtualPlantLab.Node end
    struct Node <: VirtualPlantLab.Node end
    struct BudNode <: VirtualPlantLab.Node end

    Base.@kwdef mutable struct Internode <: VirtualPlantLab.Node
        age::Int64 = 0
        biomass::Float64 = 0.0
        length::Float64 = 0.0
        width::Float64 = 0.0
        material::Lambertian{1} = Lambertian(τ = 0.1, ρ = 0.05)
        sink_strength::Float64 = 0.0
        age_previous::Float64 = 0.0
        ageD::Float64 = 0.0
        ageDprevious::Float64 = 0.0

        active::Bool = true # first 4 internodes 0 biomass
    end

    Base.@kwdef mutable struct Leaf <: VirtualPlantLab.Node
        age::Int64 = 0
        biomass::Float64 = 0.0
        length::Float64 = 0.0
        width::Float64 = 0.0
        material::Lambertian{1} = Lambertian(τ = 0.1, ρ = 0.05)
        sink_strength::Float64 = 0.0
        age_previous::Float64 = 0.0
        ageD::Float64 = 0.0
        ageDprevious::Float64 = 0.0
        delay::Float64 = 0.0 #dd to wait before growth begins
        delay_remaining::Float64 = 0.0 #dd still left to wait
    end

    Base.@kwdef mutable struct Petiole <: VirtualPlantLab.Node
        age::Int64 = 0
        biomass::Float64 = 0.0
        length::Float64 = 0.0
        width::Float64 = 0.0
        material::Lambertian{1} = Lambertian(τ = 0.1, ρ = 0.05)
        sink_strength::Float64 = 0.0
        age_previous::Float64 = 0.0
        ageD::Float64 = 0.0
        ageDprevious::Float64 = 0.0
        delay::Float64 = 0.0 #dd to wait before growth begins
        delay_remaining::Float64 = 0.0 #dd still left to wait
    end

    Base.@kwdef mutable struct Rsystem <: VirtualPlantLab.Node
        biomass::Float64 = 0.0
        age::Int64 = 0
        age_previous::Float64 = 0.0
        ageD::Float64 = 0.0
        ageDprevious::Float64 = 0.0
        sink_strength::Float64 = 0.0
        pending_nodals::Int64 = 0 
    end

    Base.@kwdef mutable struct Rmeristem <: VirtualPlantLab.Node
        age::Int64 = 0
        age_previous::Float64 = 0.0
        ageD::Float64 = 0.0
        ageDprevious::Float64 = 0.0
        rootseg::Int64 = 0
        meristem_index::Int64 = 0
        local_cumu::Float64 = 0.0
    end

    Base.@kwdef mutable struct RmeristemLateral <: VirtualPlantLab.Node
        age::Int64 = 0
        age_previous::Float64 = 0.0
        ageD::Float64 = 0.0
        ageDprevious::Float64 = 0.0
        rootseg_lateral::Int64 = 0
        lateral_index::Int64 = 0
    end

    Base.@kwdef mutable struct Rnode <: VirtualPlantLab.Node
        cumu::Float64 = 0.0
        lateral_created::Bool = false
    end

    Base.@kwdef mutable struct Rbud <: VirtualPlantLab.Node
        cumu::Float64 = 0.0
        broken::Bool = false
        age::Int64 = 0
        owner_index::Int64 = 0
    end

    struct RnodeLateral <: VirtualPlantLab.Node end
    struct Rbudnode <: VirtualPlantLab.Node end

    Base.@kwdef mutable struct Rsegment <: VirtualPlantLab.Node
        age::Int64 = 0
        biomass::Float64 = 0.0
        length::Float64 = 0.0
        width::Float64 = 0.0
        material::Lambertian{1} = Lambertian(τ = 0.1, ρ = 0.05)
        x::Float64 = 0.0
        y::Float64 = 0.0
        z::Float64 = 0.0
        sink_strength::Float64 = 0.0
        age_previous::Float64 = 0.0
        ageD::Float64 = 0.0
        ageDprevious::Float64 = 0.0
        potNup::Float64 = 0.0
        actNup::Float64 = 0.0
    end

    Base.@kwdef mutable struct RsegmentLateral <: VirtualPlantLab.Node
        age::Int64 = 0
        biomass::Float64 = 0.0
        length::Float64 = 0.0
        width::Float64 = 0.0
        material::Lambertian{1} = Lambertian(τ = 0.1, ρ = 0.05)
        x::Float64 = 0.0
        y::Float64 = 0.0
        z::Float64 = 0.0
        sink_strength::Float64 = 0.0
        age_previous::Float64 = 0.0
        ageD::Float64 = 0.0
        ageDprevious::Float64 = 0.0
        potNup::Float64 = 0.0
        actNup::Float64 = 0.0
        finerootbiomass::Float64 = 0.0
        fineroot_sinkperday::Float64 = 0.0
    end

    Base.@kwdef mutable struct treeparams
        biomass::Float64
        PAR::Float64
        RUE::Float64
        ΔB_candidate::Float64 = 0.0
        ΔB_actual::Float64 = 0.0
        IB0::Float64
        SIL::Float64          
        max_internode_width::Float64 

        SIW::Float64
        IS::Float64
        
        LB0::Float64
        SLW::Float64
        LS::Float64
        budbreak::Float64
        plastochron::Float64
        phyllochron::Float64
        tb::Float64 = 2.0
        sowing_date::Int64 = 90
        leaf_expansion::Float64
        phyllotaxis::Float64
        leaf_angle::Float64
        branch_angle::Float64
        RTD::Float64
        species::Symbol = :unknown
        N_uptake::Float64 = 0.0
        N_uptake_cumu::Float64 = 0.0
        n_leaves::Int64 = 0
        max_leaves::Int64 = 20
        seg_width::Float64 = 0.0
        seg_max_width::Float64 = 0.0
        shapeCoeff::Float64 = 0.0
        leaf_potential_biomass::Float64 = 0.0
        leaf_growth_duration::Float64 = 0.0
        Leaf_growth_respiration::Float64 = 0.0
        internode_potential_biomass::Float64 = 0.0
        internode_growth_duration::Float64 = 0.0
        internode_growth_respiration::Float64 = 0.0
        petiole_potential_biomass::Float64 = 0.0
        petiole_growth_duration::Float64 = 0.0
        petiole_growth_respiration::Float64 = 0.0
       
        
        root_cumu_length::Float64 = 0.0
        root_local_cumu::Dict{Int,Float64} = Dict{Int,Float64}()
        
        root_initial_diameter::Float64 = 0.0
        root_lateral_scale::Float64 = 0.0
        root_system_potential_biomass::Float64 = 0.0
        root_system_growth_duration::Float64 = 0.0
        root_system_growth_respiration::Float64 = 0.0
        root_main_pool::Float64 = 0.0
        root_lateral_pool::Float64 = 0.0
        root_fineroot_pool::Float64 = 0.0
        EL_main_mm::Float64 = 23.0
        EL_lateral_mult::Float64 = 2.5
        root_alloc_per_meristem::Dict{Int,Float64} = Dict{Int,Float64}()
        root_alloc_per_lateral::Dict{Int,Float64} = Dict{Int,Float64}()
        root_next_lateral_id::Int64 = 0
        rank_root::Int64 = 0
        fineroot_density::Float64 = 8.0
        fineroot_diameter::Float64 = 1e-4
        Imax::Float64 = 20000.0
        Km::Float64 = 50.0
        K2::Float64 = 1.464
        Nmin::Float64 = 2.0
        EN_coeff::Float64 = 100.0
        EN_enabled::Bool = true
        IBD::Float64 = 0.00

        root_meristem_count::Int64 = 1
        base_tilt::Float64 = 53.3 #45.0 #60 # Oynagi et al.,1994 
        base_tilt_jitter::Float64 = 3.4 # variance Oynagi et al.,1994 
        seminal_tilt::Float64 = 0.0 # first root tilt
        seminal_tilt_jitter::Float64 = 10.0 

        ER::Float64 = 0.056 # thermal time emergence rate
        ER_accumulator::Float64 = 0.0 
        ER_juvenile_lag_days::Int64 = 0 #8 in jies lag as per jies model, set 0 to dissable
        nodal_count::Int64 = 0
        max_nodals::Int64 = 6 #30 # max axial roots Jie Liu et al 


        seed_reserve_days::Int64 = 10
        seed_daily_carbon::Float64 = 0.003635  # 36.35 mg / 10 days 
        seed_nitrogen::Float64 = 0.00015   # given every day for 10 days . irrelevant note 1mg = 0.001 g, 
        carbon_reserve::Float64 = 0.0
        root_carbon_reserve::Float64 = 0.0



        # ── Per-step allocation diagnostics (written by grow!, read by record!) ──
        _diag_carbon_supply::Float64 = 0.0
        _diag_total_demand::Float64 = 0.0
        _diag_leaf_demand::Float64 = 0.0
        _diag_petiole_demand::Float64 = 0.0
        _diag_internode_demand::Float64 = 0.0
        _diag_root_demand::Float64 = 0.0
        _diag_carbon_surplus::Float64 = 0.0
        _diag_root_alloc::Float64 = 0.0
        _diag_alloc_leaves::Float64 = 0.0
        _diag_alloc_petioles::Float64 = 0.0
        _diag_alloc_internodes::Float64 = 0.0
        _diag_alloc_roots::Float64 = 0.0

        # 
        soil_x_min::Float64 = 0.0
        soil_x_max::Float64 = 0.0
        soil_y_min::Float64 = 0.0
        soil_y_max::Float64 = 0.0

    end
end
import .TreeTypes

# ═══════════════════════════════════════════════════════════════════
# SECTION: SPECIES PARAMETERS
# ═══════════════════════════════════════════════════════════════════

module SpeciesParams
    export get_species_params, species_list
    const species_list = [:cereal1, :cereal2]

    function get_species_params(spec::Symbol)
        @assert spec in species_list "unknown species: $spec"
        tp = Main.TreeTypes.treeparams
        if spec === :cereal1
            return tp(
                biomass = 2e-3, PAR = 0.0, RUE = 3.85,
                IB0 = 0.0, SIL = 0.395, max_internode_width = 0.0065, SIW = 5e5, IS = 15.0,
                LB0 = 0.0, SLW = 49, LS = 19.38,
                budbreak = 1 / 0.5, plastochron = 27.0, phyllochron = 43.0, leaf_expansion = 245.0,
                phyllotaxis = 137.5, leaf_angle = 22.0, branch_angle = 30.0,
                n_leaves = 0, max_leaves = 10, RTD = 60000,
                species = :cereal, root_initial_diameter = 0.001,
                root_lateral_scale = 0.12, IBD = 0.003,
                root_meristem_count = 1, base_tilt = 53.3, base_tilt_jitter = 3.4
            )
        else
            return tp(
                biomass = 2e-3, PAR = 0.0, RUE = 3.85,
                IB0 = 0.0, SIL = 0.395, max_internode_width = 0.00325, SIW = 5e5, IS = 15.0,
                LB0 = 0.0, SLW = 49, LS = 19.38,
                budbreak = 1 / 0.5, plastochron = 43.0, phyllochron = 86.0, leaf_expansion = 245.0,
                phyllotaxis = 137.5, leaf_angle = 55.0, branch_angle = 48.0,
                n_leaves = 0, max_leaves = 10, RTD = 68800,
                species = :legume,
                leaf_potential_biomass = 0.217, leaf_growth_duration = 245.0,
                internode_potential_biomass = 0.218, internode_growth_duration = 200.0,
                petiole_potential_biomass = 0.000, petiole_growth_duration = 220.0,
                
                Leaf_growth_respiration = 0.18, internode_growth_respiration = 0.18,
                petiole_growth_respiration = 0.18,
                
                
                root_initial_diameter = 0.001, root_lateral_scale = 0.12,
                root_system_potential_biomass = 2.05,
                root_system_growth_duration = 800.0,
                root_system_growth_respiration = 0.18,
                IBD = 0.003, root_meristem_count = 1, base_tilt = 53.3, base_tilt_jitter = 3.4, #base_tilt = 53.3 
                seminal_tilt = 0.0, seminal_tilt_jitter = 10.0, ER = 0.056, ER_accumulator = 0.0, 
                ER_juvenile_lag_days = 0, nodal_count= 0, max_nodals = 6 # 30
            )
        end
    end
end

# ═══════════════════════════════════════════════════════════════════
# SECTION: GEOMETRY (feed! implementations)
# ═══════════════════════════════════════════════════════════════════

function VirtualPlantLab.feed!(turtle::Turtle, i::TreeTypes.Internode, vars)
    rh!(turtle, vars.phyllotaxis)
    HollowCylinder!(turtle, length = i.length, height = i.width,
                    width = i.width, move = true,
                    colors = RGB(0.5, 0.4, 0.0), materials = i.material)
    return nothing
end

function VirtualPlantLab.feed!(turtle::Turtle, l::TreeTypes.Leaf, vars::TreeTypes.treeparams)
    ra!(turtle, -vars.leaf_angle)
    seg = 24
    total_bend_deg = 50.0 * (1.0 / (1.0 + exp(-0.03 * (l.ageD - 0.5 * 0.7 * vars.leaf_expansion))))
    bend_per_seg = total_bend_deg / seg
    acc_rot = 0.0
    seg_len = l.length / seg
    gain = hasproperty(vars, :width_gain) ? vars.width_gain : 1.0
    min_frac = hasproperty(vars, :min_frac) ? vars.min_frac : 0.02
    max_frac = hasproperty(vars, :max_frac) ? vars.max_frac : 1.0
    Lm = hasproperty(vars, :Lm) ? clamp(vars.Lm, 0.01, 0.99) : 0.603
    C = hasproperty(vars, :C) ? max(vars.C, 1e-6) : 0.63
    cutoff = hasproperty(vars, :cutoff) ? vars.cutoff : 0.01
    soft = hasproperty(vars, :soft) ? vars.soft : 2.0
    tip_min = hasproperty(vars, :tip_min) ? vars.tip_min : 0.002
    for j in 1:seg
        Lnorm = 1.0 - (j - 1) / (seg - 1)
        inside = (Lnorm * (2.0 * Lm - Lnorm)) / (Lm * Lm)
        inside = max(inside, 0.0)
        Wnorm = inside^C
        Wnorm = clamp(Wnorm, 0.0, 1.0)
        dist_to_edge = min(Lnorm, 1.0 - Lnorm)
        edge_factor = clamp(tip_min + (1.0 - tip_min) * (0.5 * (1 + tanh((dist_to_edge - cutoff) * soft))), 0.0, 1.0)
        raw_val = Wnorm * gain * edge_factor
        value = clamp(raw_val, min_frac, max_frac)
        seg_w = 0.5 * max(l.width, 1e-9) * value
        Rectangle!(turtle, move = true, length = seg_len, width = 2 * seg_w,
                   colors = RGB(0.1, 0.5, 0.2), materials = l.material)
        ra!(turtle, -bend_per_seg)
        acc_rot += -bend_per_seg
    end
    ra!(turtle, -acc_rot)
    ra!(turtle, vars.leaf_angle)
    return nothing
end

function VirtualPlantLab.feed!(turtle::Turtle, i::TreeTypes.Petiole, vars)
    ra!(turtle, -vars.leaf_angle)
    Rectangle!(turtle, length = i.length, width = i.width,
               move = true, colors = RGB(1.0, 0.4, 0.0), materials = i.material)
    ra!(turtle, vars.leaf_angle)
    return nothing
end

function VirtualPlantLab.feed!(turtle::Turtle, b::TreeTypes.BudNode, vars)
    ra!(turtle, -vars.branch_angle)
end

function VirtualPlantLab.feed!(turtle::Turtle, r::TreeTypes.Rsegment, vars)
    rh!(turtle, vars.phyllotaxis)
    rv!(turtle, -0.05) # 0.1 gr 
    HollowCylinder!(turtle, length = r.length, height = r.width,
                    width = r.width, move = true,
                    colors = RGB(1, 0.2, 0.0), materials = r.material)
    t_pos = pos(turtle)
    t_head = head(turtle)
    center = t_pos .+ 0.5 * r.length * t_head
    r.x = center[1]; r.y = center[2]; r.z = center[3]

    # ── Periodic boundary: wrap turtle XY position for next segment ──
    if vars.soil_x_max > vars.soil_x_min
        wx = wrap_coord(t_pos[1], vars.soil_x_min, vars.soil_x_max)
        wy = wrap_coord(t_pos[2], vars.soil_y_min, vars.soil_y_max)
        if wx != t_pos[1] || wy != t_pos[2]
            t!(turtle; to = Vec(wx, wy, t_pos[3]))
        end
    end
    return nothing
end

function VirtualPlantLab.feed!(turtle::Turtle, rl::TreeTypes.RsegmentLateral, vars)
    rh!(turtle, vars.phyllotaxis)
    rv!(turtle, -0.05) # Jie
    HollowCylinder!(turtle, length = rl.length, height = rl.width,
                    width = rl.width, move = true,
                    colors = RGB(0.0, 0.2, 0.9), materials = rl.material)
    t_pos = pos(turtle)
    t_head = head(turtle)
    center = t_pos .+ 0.5 * rl.length * t_head
    rl.x = center[1]; rl.y = center[2]; rl.z = center[3]

    # ── Periodic boundary: wrap turtle XY position for next segment ──
    if vars.soil_x_max > vars.soil_x_min
        wx = wrap_coord(t_pos[1], vars.soil_x_min, vars.soil_x_max)
        wy = wrap_coord(t_pos[2], vars.soil_y_min, vars.soil_y_max)
        if wx != t_pos[1] || wy != t_pos[2]
            t!(turtle; to = Vec(wx, wy, t_pos[3]))
        end
    end
    return nothing
end

function VirtualPlantLab.feed!(turtle::Turtle, rb::TreeTypes.Rbudnode, vars)
    ra!(turtle, 60.0- rand() * 40.0)
    return nothing
end

# ═══════════════════════════════════════════════════════════════════
# SECTION: HELPERS (root costs, fine roots, wrapper for root cloner,  elongation)
# ═══════════════════════════════════════════════════════════════════

function segment_creation_cost(tvars; lateral = false)
    L = tvars.IBD 
    d = lateral ? tvars.root_initial_diameter * tvars.root_lateral_scale : tvars.root_initial_diameter
    vol = pi * (d / 2)^2 * L
    return vol * tvars.RTD
end

function finerootsink_per_lateral_segment(rseg_lateral, tvars)
    d = tvars.fineroot_diameter
    EL_main_m = tvars.EL_main_mm
    PER = EL_main_m * d
    A_fineroot = π * (d / 2)^2
    finerootsink_day = PER * A_fineroot * tvars.RTD
    L_segment = rseg_lateral.length
    maxfineroot = tvars.fineroot_density * L_segment * A_fineroot * tvars.RTD
    return finerootsink_day, maxfineroot
end

# ══════════════════════════════════════════════════════════════════
# SECTION: PERIODIC BOUNDARY — root coordinate wrapping
# ══════════════════════════════════════════════════════════════════

# max/min grid dimensions
@inline function wrap_coord(v::Float64, vmin::Float64, vmax::Float64)
    span = vmax - vmin
    span <= 0.0 && return v          # degenerate grid — do nothing
    return vmin + mod(v - vmin, span)
end


function compute_elongation!(tree, root_alloc_total)
    tvars = data(tree)
    empty!(tvars.root_alloc_per_meristem)
    empty!(tvars.root_alloc_per_lateral)
    root_alloc_total += tvars.root_carbon_reserve
    tvars.root_carbon_reserve = 0.0
    seg_length = tvars.IBD 
    d_main = tvars.root_initial_diameter
    d_lat = d_main * tvars.root_lateral_scale
    EL_main_m = tvars.EL_main_mm
    EL_lat_m = EL_main_m * tvars.EL_lateral_mult
    rmer_list = get_Rmeristems(tree)
    lrmer_list = get_RmeristemsLateral(tree)
    all_RsegmentsLateral = get_RsegmentsLateral(tree)

    # Potential costs: main
    pot_costs_main = Dict{Int,Float64}()
    total_pot_main = 0.0
    PER_main = EL_main_m * d_main
    for rmer in rmer_list
        merid = rmer.meristem_index
        area = pi * (d_main / 2)^2
        pot_cost = PER_main * area * tvars.RTD
        pot_costs_main[merid] = pot_cost
        total_pot_main += pot_cost
    end

    # Potential costs: lateral
    pot_costs_lat = Dict{Int,Float64}()
    total_pot_lat = 0.0
    PER_lat = EL_lat_m * d_lat
    for lrmer in lrmer_list
        latid = lrmer.lateral_index
        area = pi * (d_lat / 2)^2
        pot_cost = PER_lat * area * tvars.RTD
        pot_costs_lat[latid] = pot_cost
        total_pot_lat += pot_cost
    end

    # Potential costs: fine roots
    total_pot_fineroot = 0.0
    fineroot_pot_per_segment = Dict{TreeTypes.RsegmentLateral,Float64}()
    for rseg_lat in all_RsegmentsLateral
        fineroot_day, fineroot_Mpot = finerootsink_per_lateral_segment(rseg_lat, tvars)
        if rseg_lat.finerootbiomass < fineroot_Mpot
            fineroot_pot_per_segment[rseg_lat] = fineroot_day
            total_pot_fineroot += fineroot_day
        else
            fineroot_pot_per_segment[rseg_lat] = 0.0
        end
    end

    # 3-way proportional split
    denom = total_pot_main + total_pot_lat + total_pot_fineroot
    if denom <= 0
        for rmer in rmer_list; rmer.rootseg = 0; end
        for lrmer in lrmer_list; lrmer.rootseg_lateral = 0; end
        tvars.root_main_pool = 0.0; tvars.root_lateral_pool = 0.0
        return nothing
    end
    alloc_main_total = (total_pot_main / denom) * root_alloc_total
    alloc_lat_total = (total_pot_lat / denom) * root_alloc_total
    alloc_fineroot_total = (total_pot_fineroot / denom) * root_alloc_total
    tvars.root_main_pool = alloc_main_total
    tvars.root_lateral_pool = alloc_lat_total
    tvars.root_fineroot_pool = alloc_fineroot_total

    # Distribute to main meristems
    for rmer in rmer_list
        merid = rmer.meristem_index
        pot = get(pot_costs_main, merid, 0.0)
        alloc = (total_pot_main > 0) ? alloc_main_total * (pot / total_pot_main) : 0.0
        tvars.root_alloc_per_meristem[merid] = alloc
        frac = (pot > 0) ? min(1.0, alloc / pot) : 0.0
        actual_length = PER_main * frac
        rmer.rootseg = max(0, Int(floor(actual_length / seg_length)))
    end

    # Distribute to lateral meristems
    for lrmer in lrmer_list
        latid = lrmer.lateral_index
        pot = get(pot_costs_lat, latid, 0.0)
        alloc = (total_pot_lat > 0) ? alloc_lat_total * (pot / total_pot_lat) : 0.0
        tvars.root_alloc_per_lateral[latid] = alloc
        frac = (pot > 0) ? min(1.0, alloc / pot) : 0.0
        actual_length = PER_lat * frac
        lrmer.rootseg_lateral = max(0, Int(floor(actual_length / seg_length)))
    end

    # Distribute to fine roots
    for rseg_lat in all_RsegmentsLateral
        pot = get(fineroot_pot_per_segment, rseg_lat, 0.0)
        if pot > 0 && total_pot_fineroot > 0
            alloc = alloc_fineroot_total * (pot / total_pot_fineroot)
            rseg_lat.finerootbiomass += alloc
            _, fineroot_Mpot = finerootsink_per_lateral_segment(rseg_lat, tvars)
            if rseg_lat.finerootbiomass > fineroot_Mpot
                rseg_lat.finerootbiomass = fineroot_Mpot
            end
        end
    end
    return nothing
end

# ═══════════════════════════════════════════════════════════════════
# SECTION: DEVELOPMENT RULES
# ═══════════════════════════════════════════════════════════════════

# function create_meristem_rule(vleaf, vint, vpet)
#     Rule(TreeTypes.Meristem,
#          lhs = mer -> data(mer).ageD >= graph_data(mer).plastochron,
#          rhs = mer -> begin
#              tvars = graph_data(mer)
#              if tvars.n_leaves >= tvars.max_leaves
#                  return nothing
#              end
#              tvars.n_leaves += 1
#             #  int_biomass = tvars.n_leaves <= 4 ? 0.0 : vint.biomass
#             #  int_length  = tvars.n_leaves <= 4 ? 0.0 : vint.length
#             #  int_width   = tvars.n_leaves <= 4 ? 0.0 : vint.width
#              TreeTypes.Node() +
#              TreeTypes.Internode(biomass = vint.biomass, length = vint.length, width = vint.width) +
#              (TreeTypes.Bud(),
#               TreeTypes.Petiole(biomass = vpet.biomass, length = vpet.length, width = vpet.width) +
#               TreeTypes.Leaf(biomass = vleaf.biomass, length = vleaf.length, width = vleaf.width)) +
#              TreeTypes.Meristem(ageD = data(mer).ageD - graph_data(mer).plastochron)
#          end)
# end

# function create_meristem_rule(vleaf, vint, vpet)
#     Rule(TreeTypes.Meristem,
#          lhs = mer -> data(mer).ageD >= graph_data(mer).plastochron,
#          rhs = mer -> begin
#              tvars = graph_data(mer)
#              if tvars.n_leaves >= tvars.max_leaves
#                  return nothing
#              end
#              tvars.n_leaves += 1
#               int_biomass = tvars.n_leaves <= 4 ? 0.0 : vint.biomass
#               int_length  = tvars.n_leaves <= 4 ? 0.0 : vint.length
#               int_width   = tvars.n_leaves <= 4 ? 0.0 : vint.width
#              TreeTypes.Node() +
#              TreeTypes.Internode(biomass = vint.biomass, length = vint.length, width = vint.width) +
#              (TreeTypes.Bud(),
#               TreeTypes.Petiole(biomass = vpet.biomass, length = vpet.length, width = vpet.width) +
#               TreeTypes.Leaf(biomass = vleaf.biomass, length = vleaf.length, width = vleaf.width)) +
#              TreeTypes.Meristem(ageD = data(mer).ageD - graph_data(mer).plastochron)
#          end)
# end


function create_meristem_rule(vleaf, vint, vpet)
    Rule(TreeTypes.Meristem,
         lhs = mer -> data(mer).ageD >= graph_data(mer).plastochron,
         rhs = mer -> begin
             tvars = graph_data(mer)
             if tvars.n_leaves >= tvars.max_leaves
                 return nothing
             end
             tvars.n_leaves += 1

             # Fixed delay: phyllochron - plastochron
             organ_delay = tvars.phyllochron - tvars.plastochron

             if tvars.n_leaves <= 4
              int_bm = 0.0; int_l = 0.0; int_w = 0.0
              int_active = false
            else
              int_bm = vint.biomass; int_l = vint.length; int_w = vint.width
              int_active = true
end


TreeTypes.Node() +
TreeTypes.Internode(biomass = int_bm, length = int_l, width = int_w, active = int_active) +
             (TreeTypes.Bud(),
              TreeTypes.Petiole(biomass = vpet.biomass, length = vpet.length, width = vpet.width,
                                delay = organ_delay, delay_remaining = organ_delay) +
              TreeTypes.Leaf(biomass = vleaf.biomass, length = vleaf.length, width = vleaf.width,
                             delay = organ_delay, delay_remaining = organ_delay)) +
             TreeTypes.Meristem(ageD = data(mer).ageD - graph_data(mer).plastochron)
         end)
end


function fixed_break(_bud)
    return rand() < 0.000
end

function create_branch_rule(vint)
    Rule(TreeTypes.Bud,
         lhs = fixed_break,
         rhs = bud -> TreeTypes.BudNode() +
                      TreeTypes.Internode(biomass = vint.biomass, length = vint.length, width = vint.width) +
                      TreeTypes.Meristem())
end

function create_root_meristem_rule()
    Rule(TreeTypes.Rmeristem,
         lhs = rmer -> data(rmer).ageD >= graph_data(rmer).plastochron,
         rhs = rmer -> begin
             tvars = graph_data(rmer)
             nd = data(rmer)
             merid = nd.meristem_index
             cost = segment_creation_cost(tvars; lateral = false)
             if cost <= 0.0; return TreeTypes.Rmeristem(); end
             desired = Int(nd.rootseg)
             allocated = get(tvars.root_alloc_per_meristem, merid, 0.0)
             can_pay = Int(floor(allocated / cost))
             n_create = min(desired, can_pay)
             if n_create <= 0
                 return TreeTypes.Rmeristem(rootseg = nd.rootseg, meristem_index = merid,
                                            ageD = nd.ageD - graph_data(rmer).plastochron)
             end
             tvars.root_alloc_per_meristem[merid] = max(0.0, allocated - n_create * cost)
             tvars.root_carbon_reserve += max(0.0, allocated - n_create * cost)
             tvars.root_main_pool = max(0.0, tvars.root_main_pool - n_create * cost)
             created_units = [
                 begin
                     seg = TreeTypes.Rsegment(biomass = cost, length = tvars.IBD ,
                                              width = tvars.root_initial_diameter)
                     tvars.root_local_cumu[merid] += seg.length
                     local_cumu = tvars.root_local_cumu[merid]
                     node = TreeTypes.Rnode(cumu = local_cumu, lateral_created = false)
                     prev_cumu = local_cumu - seg.length
                     n_IBD_before = floor(prev_cumu / tvars.IBD)
                     n_IBD_after = floor(local_cumu / tvars.IBD)
                     if n_IBD_after > n_IBD_before
                         bud = TreeTypes.Rbud(cumu = local_cumu, broken = false,
                                              owner_index = merid)
                         seg + (node, bud)
                     else
                         seg + node
                     end
                 end
                 for i in 1:n_create
             ]
             return foldl((a, b) -> a + b, created_units) +
                    TreeTypes.Rmeristem(rootseg = nd.rootseg, meristem_index = merid,
                                        ageD = nd.ageD - graph_data(rmer).plastochron)
         end)
end

function create_rsystem_emergence_rule()
    Rule(TreeTypes.Rsystem,
        lhs = rs -> data(rs).pending_nodals > 0,
        rhs = rs -> begin
            tvars = graph_data(rs)
            n = data(rs).pending_nodals
            data(rs).pending_nodals = 0   # consume

            new_axes = Any[]
            for k in 1:n
                idx_local      = tvars.nodal_count - n + k
                meristem_index = idx_local + 1   # 1 = seminal

                az    = mod(idx_local * 137.50776405003785, 360.0) # golden angle for new roots
                tilt = tvars.base_tilt +
                       (rand() * 2 - 1) * tvars.base_tilt_jitter

                tvars.root_local_cumu[meristem_index] = 0.0

                axis = RH(az) + RA(-tilt) +
                       TreeTypes.Rmeristem(meristem_index = meristem_index)
                push!(new_axes, axis)
            end

            return TreeTypes.Rsystem(
                       biomass        = data(rs).biomass,
                       age            = data(rs).age,
                       age_previous   = data(rs).age_previous,
                       ageD           = data(rs).ageD,
                       ageDprevious   = data(rs).ageDprevious,
                       sink_strength  = data(rs).sink_strength,
                       pending_nodals = 0
                   ) + Tuple(a for a in new_axes)   # ← no RA(180.0) here
        end)
end

function create_root_node_to_lateral_rule(vrseg_lateral)
    Rule(TreeTypes.Rbud,
        lhs = bud -> begin
            tvars = graph_data(bud)
            bd = data(bud)
            owner = bd.owner_index
            progress = tvars.root_local_cumu[owner] - bd.cumu
            return (progress >= tvars.IBD) && (!bd.broken)
        end,
        rhs = bud -> begin
            bd = data(bud)
            tvars = graph_data(bud)
            owner = bd.owner_index
            bd.broken = true
            cost = segment_creation_cost(tvars; lateral = true)
            tvars.root_next_lateral_id += 1
            new_lat_id = tvars.root_next_lateral_id
            if tvars.root_lateral_pool >= cost && cost > 0.0
                tvars.root_lateral_pool -= cost
                segl = TreeTypes.RsegmentLateral(biomass = cost,
                    length = vrseg_lateral.length, width = vrseg_lateral.width)
                rmerlat = TreeTypes.RmeristemLateral(rootseg_lateral = 0, lateral_index = new_lat_id)
                tvars.root_alloc_per_lateral[new_lat_id] = get(tvars.root_alloc_per_lateral, new_lat_id, 0.0)
                return TreeTypes.Rbudnode() + segl + rmerlat
            else
                rmerlat = TreeTypes.RmeristemLateral(rootseg_lateral = 0, lateral_index = new_lat_id)
                tvars.root_alloc_per_lateral[new_lat_id] = get(tvars.root_alloc_per_lateral, new_lat_id, 0.0)
                return TreeTypes.Rbudnode() + rmerlat
            end
        end)
end

function create_lateral_root_meristem_rule()
    Rule(TreeTypes.RmeristemLateral,
         lhs = rmer -> data(rmer).ageD >= graph_data(rmer).plastochron,
         rhs = rmer -> begin
             tvars = graph_data(rmer)
             nd = data(rmer)
             cost = segment_creation_cost(tvars; lateral = true)
             if cost <= 0.0; return TreeTypes.RmeristemLateral(); end
             desired = Int(nd.rootseg_lateral)
             lat_id = nd.lateral_index
             allocated = get(tvars.root_alloc_per_lateral, lat_id, 0.0)
             can_pay = Int(floor(allocated / cost))
             n_create = min(desired, can_pay)
             if n_create <= 0
                 return TreeTypes.RmeristemLateral(rootseg_lateral = nd.rootseg_lateral,
                                                   lateral_index = lat_id,
                                                   ageD = nd.ageD - graph_data(rmer).plastochron)
             end
             tvars.root_alloc_per_lateral[lat_id] = max(0.0, allocated - n_create * cost)
             tvars.root_carbon_reserve += max(0.0, allocated - n_create * cost)
             tvars.root_lateral_pool = max(0.0, tvars.root_lateral_pool - n_create * cost)
             created = Any[]
             for i in 1:n_create
                 segl = TreeTypes.RsegmentLateral(biomass = cost, length = tvars.IBD,
                     width = tvars.root_initial_diameter * tvars.root_lateral_scale)
                 push!(created, segl)
                 push!(created, TreeTypes.RnodeLateral())
             end
             push!(created, TreeTypes.RmeristemLateral(rootseg_lateral = nd.rootseg_lateral,
                                                        lateral_index = lat_id,
                                                        ageD = nd.ageD - graph_data(rmer).plastochron))
             return foldl((a, b) -> a + b, created)
         end)
end

# ═══════════════════════════════════════════════════════════════════
# SECTION: GROWTH (dimensions, sinks, allocation)
# ═══════════════════════════════════════════════════════════════════

function leaf_dims(biomass, vars)
    leaf_area = biomass / vars.SLW
    leaf_length = sqrt(leaf_area * 4 * vars.LS / pi)
    leaf_width = leaf_length / vars.LS
    return leaf_length, leaf_width
end

function int_dims(biomass, vars, ageD)
    int_length = biomass * vars.SIL
    int_width  = vars.max_internode_width / (1 + exp(-0.02 * (ageD - 0.5 * vars.internode_growth_duration)))
    return int_length, int_width
end

function pet_dims(biomass, vars)
    pet_volume = biomass / vars.SIW
    pet_length = cbrt(pet_volume * 4 * vars.IS^2 / pi)
    pet_width = pet_length / vars.IS
    return pet_length, pet_width
end

function Rseg_dims(biomass, vars)
    return vars.IBD, vars.root_initial_diameter
end

function Rseg_lateral_dims(biomass, vars)
    return vars.IBD, vars.root_initial_diameter * vars.root_lateral_scale
end

function estimate_sink_strength!(organ, organ_type::Symbol, vars)
    if organ_type == :leaf
        wmax = vars.leaf_potential_biomass; te = vars.leaf_growth_duration; r = vars.Leaf_growth_respiration
    elseif organ_type == :petiole
        wmax = vars.petiole_potential_biomass; te = vars.petiole_growth_duration; r = vars.petiole_growth_respiration
    elseif organ_type == :internode
        wmax = vars.internode_potential_biomass; te = vars.internode_growth_duration; r = vars.internode_growth_respiration
    elseif organ_type == :root_system
        wmax = vars.root_system_potential_biomass; te = vars.root_system_growth_duration; r = vars.root_system_growth_respiration
    else
        error("Unknown organ type: $organ_type")
    end
    t = min(organ.ageD, te)
    tm = te / 2.0
    max_sink_strength = wmax * ((2 * te - tm) / (te * (te - tm))) * (tm / te)^(tm / (te - tm))
    sink_strength = max_sink_strength / (1 - r) * ((te - t) / (te - tm)) * (t / tm)^(tm / (te - tm))
    organ.sink_strength = sink_strength * (t - organ.ageDprevious)
    organ.ageDprevious = organ.ageD
    return organ.sink_strength
end

# Queries
get_leaves(tree) = apply(tree, Query(TreeTypes.Leaf))
get_petioles(tree) = apply(tree, Query(TreeTypes.Petiole))
get_internodes(tree) = apply(tree, Query(TreeTypes.Internode))
get_Rsegments(tree) = apply(tree, Query(TreeTypes.Rsegment))
get_meristems(tree) = apply(tree, Query(TreeTypes.Meristem))
get_Rmeristems(tree) = apply(tree, Query(TreeTypes.Rmeristem))
get_RsegmentsLateral(tree) = apply(tree, Query(TreeTypes.RsegmentLateral))
get_RmeristemsLateral(tree) = apply(tree, Query(TreeTypes.RmeristemLateral))
get_Rsystems(tree) = apply(tree, Query(TreeTypes.Rsystem))

function age!(all_leaves, all_petioles, all_internodes,
              all_meristems, all_Rmeristems, all_Rsegments,
              all_RsegmentsLateral, all_RmeristemsLateral, all_Rsystems,
              tvars, dayofyear)
    ΔTT = max(0.0, temperature(dayofyear) - tvars.tb)

    for leaf in all_leaves
        if leaf.delay_remaining > 0.0
            leaf.delay_remaining = max(0.0, leaf.delay_remaining - ΔTT)
        else
            leaf.age += 1
            leaf.ageDprevious = leaf.ageD
            leaf.ageD += ΔTT
        end
    end
    for pet in all_petioles
        if pet.delay_remaining > 0.0
            pet.delay_remaining = max(0.0, pet.delay_remaining - ΔTT)
        else
            pet.age += 1
            pet.ageDprevious = pet.ageD
            pet.ageD += ΔTT
        end
    end
    for int in all_internodes
        int.age += 1
        int.ageDprevious = int.ageD
        int.ageD += ΔTT
    end
    for mer in all_meristems
        mer.age += 1
        mer.ageDprevious = mer.ageD
        mer.ageD += ΔTT
    end
    for Rmer in all_Rmeristems
        Rmer.age += 1
        Rmer.ageDprevious = Rmer.ageD
        Rmer.ageD += ΔTT
    end
    for Rseg in all_Rsegments
        Rseg.age += 1
        Rseg.ageDprevious = Rseg.ageD
        Rseg.ageD += ΔTT
    end
    for RsegLat in all_RsegmentsLateral
        RsegLat.age += 1
        RsegLat.ageDprevious = RsegLat.ageD
        RsegLat.ageD += ΔTT
    end
    for RmerLat in all_RmeristemsLateral
        RmerLat.age += 1
        RmerLat.ageDprevious = RmerLat.ageD
        RmerLat.ageD += ΔTT
    end
    for Rsystem in all_Rsystems
        Rsystem.age += 1
        Rsystem.ageDprevious = Rsystem.ageD
        Rsystem.ageD += ΔTT
        # ── Nodal-root emergence (ER mechanism) ──
    # Cap reached? (1 = seminal already in axiom)
    current_axes = 1 + tvars.nodal_count
         if current_axes < tvars.max_nodals
             # Optional juvenile lag before the first nodal only
             past_lag = (tvars.nodal_count > 0) ||
                   (Rsystem.age >= tvars.ER_juvenile_lag_days)
            if past_lag
            tvars.ER_accumulator += ΔTT * tvars.ER
            while tvars.ER_accumulator >= 1.0 &&
                  (1 + tvars.nodal_count) < tvars.max_nodals
                tvars.ER_accumulator -= 1.0
                tvars.nodal_count    += 1
                Rsystem.pending_nodals += 1
            end
        end
    end
end
    #end
    return nothing
end

function total_above_biomass(tree)
    above = 0.0
    for l in get_leaves(tree); above += l.biomass; end
    for p in get_petioles(tree); above += p.biomass; end
    for i in get_internodes(tree); above += i.biomass; end
    return above
end

function total_below_biomass(tree)
    below = 0.0
    for r in get_Rsegments(tree); below += r.biomass; end
    for r in get_RsegmentsLateral(tree); below += r.biomass; end
    return below
end

# clamp CNcrit to a realistic maximum (e.g., 6% = 0.06):
function compute_fn(tree; a1 = 3.68, a2 = 0.28, clamp01 = true, eps = 1e-12)
    tvars = data(tree)
    above = total_above_biomass(tree)
    N_tot = max(tvars.N_uptake_cumu, eps)
    CN = N_tot / above
    safe_above = max(above, eps)
    CNcrit_raw = (a1 * safe_above^(-a2)) / 100.0

    # ── NEW: Clamp to physiological maximum (6% N = 0.06 g/g) ──
    CNcrit = min(CNcrit_raw, 0.06)

    safe_CNcrit = max(CNcrit, eps)
    fn = CN / safe_CNcrit
    if clamp01; fn = clamp(fn, 0.0, 1.0); end
    return fn, CN, CNcrit, above
end

# In grow!(), REPLACE the hardcoded lines:
#   seed_reserve_days = 20
#   seed_daily_carbon = 0.5
# WITH:

function grow!(tree, all_leaves, all_petioles, all_internodes, all_Rsegments, all_RsegmentsLateral)
    tvars = data(tree)
    fn, CN, CNcrit, above = compute_fn(tree; a1 = 3.68, a2 = 0.28, clamp01 = true)

    # ── Read seed params from treeparams (no longer hardcoded) ──
    seed_reserve_days = tvars.seed_reserve_days
    seed_daily_carbon = tvars.seed_daily_carbon

    all_Rsystems = apply(tree, Query(TreeTypes.Rsystem))
    plant_age = length(all_Rsystems) > 0 ? all_Rsystems[1].age : 0


   # ── Seed nitrogen: release daily over same window as carbon reserve ──
if plant_age >= 1 && plant_age <= seed_reserve_days
    tvars.N_uptake_cumu += tvars.seed_nitrogen / seed_reserve_days
end
    #
    if plant_age <= seed_reserve_days
        fn_floor = max(0.0, 1.0 - plant_age / seed_reserve_days)
        fn = max(fn, fn_floor)
    end

    ΔB_photo = tvars.RUE * (tvars.PAR / 1e6) * fn  
    
    #ΔB_photo = tvars.RUE * (3.0e4 / 1e6) * fn                ################################## boost db by 50


  if plant_age <= seed_reserve_days
    ΔB_seed = seed_daily_carbon
else
    ΔB_seed = 0.0
end

    ΔB = ΔB_photo + ΔB_seed
    tvars.ΔB_candidate = ΔB_photo
    tvars.ΔB_actual = ΔB
    tvars.biomass += ΔB
    ΔB += tvars.carbon_reserve
    tvars.carbon_reserve = 0.0

    # ... rest of grow! unchanged ...
    for leaf in all_leaves; estimate_sink_strength!(leaf, :leaf, tvars); end
    for pet in all_petioles; estimate_sink_strength!(pet, :petiole, tvars); end
    
    for int in all_internodes
    int.active || continue
    estimate_sink_strength!(int, :internode, tvars)
    end

    root_system = length(all_Rsystems) >= 1 ? all_Rsystems[1] : nothing
    if root_system !== nothing; estimate_sink_strength!(root_system, :root_system, tvars); end

    leaf_demands = Float64[]
    for leaf in all_leaves
        headroom = max(0.0, tvars.leaf_potential_biomass - leaf.biomass)
        gross_headroom = headroom / (1.0 - tvars.Leaf_growth_respiration)
        push!(leaf_demands, min(leaf.sink_strength, gross_headroom))
    end
    pet_demands = Float64[]
    for pet in all_petioles
        headroom = max(0.0, tvars.petiole_potential_biomass - pet.biomass)
        gross_headroom = headroom / (1.0 - tvars.petiole_growth_respiration)
        push!(pet_demands, min(pet.sink_strength, gross_headroom))
    end
    int_demands = Float64[]
    for int in all_internodes
        # if int.biomass == 0.0
        #     push!(int_demands, 0.0)
        # else
            headroom = max(0.0, tvars.internode_potential_biomass - int.biomass)
            gross_headroom = headroom / (1.0 - tvars.internode_growth_respiration)
            push!(int_demands, min(int.sink_strength, gross_headroom))
        #end
    end
    root_demand = (root_system !== nothing) ? root_system.sink_strength : 0.0
    total_demand = sum(leaf_demands) + sum(pet_demands) + sum(int_demands) + root_demand + 1e-12
    if total_demand <= 0
        tvars._diag_carbon_supply = ΔB
        tvars._diag_total_demand = 0.0
        tvars._diag_leaf_demand = 0.0
        tvars._diag_petiole_demand = 0.0
        tvars._diag_internode_demand = 0.0
        tvars._diag_root_demand = 0.0
        tvars._diag_carbon_surplus = ΔB
        tvars._diag_root_alloc = 0.0
        tvars._diag_alloc_leaves = 0.0
        tvars._diag_alloc_petioles = 0.0
        tvars._diag_alloc_internodes = 0.0
        tvars._diag_alloc_roots = 0.0
        tvars.carbon_reserve = ΔB
        return nothing
    end
    ΔB_used = min(ΔB, total_demand)

    total_leaf_alloc = 0.0
    for (i, leaf) in enumerate(all_leaves)
        alloc_gross = ΔB_used * leaf_demands[i] / total_demand
        leaf.biomass += alloc_gross * (1.0 - tvars.Leaf_growth_respiration)
        total_leaf_alloc += alloc_gross
    end
    total_pet_alloc = 0.0
    for (i, pet) in enumerate(all_petioles)
        alloc_gross = ΔB_used * pet_demands[i] / total_demand
        pet.biomass += alloc_gross * (1.0 - tvars.petiole_growth_respiration)
        total_pet_alloc += alloc_gross
    end
    total_int_alloc = 0.0
    for (i, int) in enumerate(all_internodes)
        alloc_gross = ΔB_used * int_demands[i] / total_demand
        int.biomass += alloc_gross * (1.0 - tvars.internode_growth_respiration)
        total_int_alloc += alloc_gross
    end
    root_alloc_total = ΔB_used * root_demand / total_demand * (1.0 - tvars.root_system_growth_respiration)
    tvars.carbon_reserve = max(0.0, ΔB - ΔB_used)
    tvars._diag_carbon_supply = ΔB
    tvars._diag_total_demand = total_demand
    tvars._diag_leaf_demand = sum(leaf_demands)
    tvars._diag_petiole_demand = sum(pet_demands)
    tvars._diag_internode_demand = sum(int_demands)
    tvars._diag_root_demand = root_demand
    tvars._diag_carbon_surplus = max(0.0, ΔB - ΔB_used)
    tvars._diag_root_alloc = root_alloc_total
    tvars._diag_alloc_leaves = total_leaf_alloc
    tvars._diag_alloc_petioles = total_pet_alloc
    tvars._diag_alloc_internodes = total_int_alloc
    tvars._diag_alloc_roots = ΔB_used * root_demand / total_demand
    compute_elongation!(tree, root_alloc_total)
    return nothing
end

function size_leaves!(all_leaves, tvars)
    for leaf in all_leaves; leaf.length, leaf.width = leaf_dims(leaf.biomass, tvars); end
end

function size_petioles!(all_petioles, tvars)
    for pet in all_petioles; pet.length, pet.width = pet_dims(pet.biomass, tvars); end
end

function size_internodes!(all_internodes, tvars)
    for (idx, int) in enumerate(all_internodes)
        if idx <= 4
            int.biomass = 0.0; int.length = 0.0; int.width = 0.0
        else
            int.length, int.width = int_dims(int.biomass, tvars, int.ageD)
        end
    end
end

function size_Rsegments!(all_Rsegments, tvars)
    for Rsegment in all_Rsegments
        Rsegment.length = tvars.IBD; Rsegment.width = tvars.root_initial_diameter
    end
end

function size_RsegmentsLateral!(all_RsegmentsLateral, tvars)
    for RsegmentLateral in all_RsegmentsLateral
        RsegmentLateral.length = tvars.IBD
        RsegmentLateral.width = tvars.root_initial_diameter * tvars.root_lateral_scale
    end
end

# ═══════════════════════════════════════════════════════════════════
# SECTION: SOIL
# ═══════════════════════════════════════════════════════════════════

module SoilTypes
    using VirtualPlantLab, ColorTypes
    import Parameters: @with_kw
    export Soilcell, get_soilcells

    @with_kw mutable struct Soilcell <: VirtualPlantLab.Node
        length::Float64 = 0.01
        width::Float64 = 0.01
        height::Float64 = 0.01
        x::Float64 = 0.0
        y::Float64 = 0.0
        z::Float64 = 0.0
        i::Int64 = 0
        j::Int64 = 0
        k::Int64 = 0
        N::Float64 = 2000.0 
        material::Lambertian{1} = Lambertian(τ = 0.1, ρ = 0.05)
    end

    function VirtualPlantLab.feed!(turtle::Turtle, ss::Soilcell, vars)
        corner = pos(turtle)
        ss.x = corner[1]
        ss.y = corner[2]
        ss.z = corner[3] + ss.height / 2.0
        HollowCube!(turtle,
            length = ss.length, width = ss.width, height = ss.height,
            colors = RGBA(0.0, 1.0, 0.0, 0.1), materials = ss.material)
        return nothing
    end

    get_soilcells(soil) = apply(soil, Query(Soilcell))
end
using .SoilTypes: Soilcell, get_soilcells

Base.@kwdef struct Soil_surface <: VirtualPlantLab.Node
    length::Float64
    width::Float64
end

function VirtualPlantLab.feed!(turtle::Turtle, s::Soil_surface, data)
    Rectangle!(turtle, length = s.length, width = s.width,
               colors = RGB(255 / 255, 236 / 255, 179 / 255),
               materials = Lambertian(τ = 0.0, ρ = 0.21))
end

# ═══════════════════════════════════════════════════════════════════
# SECTION: NITROGEN UPTAKE
# ═══════════════════════════════════════════════════════════════════

function total_soil_N(soil)
    sum(sc -> sc.N, get_soilcells(soil))
end

voxel_L(sc) = sc.length * sc.width * sc.height * 1000.0

function soil_conc(sc)
    V = voxel_L(sc)
    return V > 0 ? sc.N / V : 0.0
end

function soil_available_umol(sc, vars)
    V = voxel_L(sc)
    return max(0.0, sc.N - vars.Nmin * V)
end

function root_geom_area_mass(rs, vars)
    L = rs.length; D = rs.width
    A = π * D * L
    vol = π * (D / 2)^2 * L
    rootDW_g = vars.RTD * vol
    return A, rootDW_g
end

# ── Shared containment predicate (Feature 2, Step 1) ──
@inline function root_in_cell(rs, sc; tol = 1e-10)
    return abs(rs.x - sc.x) <= (sc.length / 2 + rs.width / 2 + tol) &&
           abs(rs.y - sc.y) <= (sc.width  / 2 + rs.width / 2 + tol) &&
           abs(rs.z - sc.z) <= (sc.height / 2 + rs.width / 2 + tol)
end

function nitrogen_uptake!(tree, soil)
    root_seg = get_Rsegments(tree)
    root_seg_lateral = get_RsegmentsLateral(tree)
    soil_cells = get_soilcells(soil)
    total_plant_N_umol = 0.0
    tol = 1e-10
    vars = data(tree)

    nsc = length(soil_cells)
    pot_per_cell = zeros(Float64, nsc)

    for rs in root_seg; rs.potNup = 0.0; rs.actNup = 0.0; end
    for rs in root_seg_lateral; rs.potNup = 0.0; rs.actNup = 0.0; end

   # EN downregulation factor
    if vars.EN_enabled
        above = total_above_biomass(tree)
        safe_above = max(above, 1e-12)
        # Use same basis as compute_fn: cumulative N / above-ground organ biomass
        Nc_plant = vars.N_uptake_cumu / safe_above
        CN_crit = (3.68 * safe_above^(-0.28)) / 100.0
        CN_crit = min(CN_crit, 0.06)  # same clamp as compute_fn

        if Nc_plant > CN_crit
            EN = clamp(exp(-vars.EN_coeff * (Nc_plant - CN_crit)), 0.0, 1.0)
        else
            EN = 1.0
        end
    else
        EN = 1.0
    end

    # Phase A: compute potentials
    for (si, sc) in enumerate(soil_cells)
        V = voxel_L(sc)
        C = V > 0 ? sc.N / V : 0.0
        Ceff = max(0.0, C - vars.Nmin)
        Ihat = (Ceff > 0.0) ? (vars.Imax * Ceff / (vars.Km + Ceff)) : 0.0
        Ilat = vars.K2 * Ceff

        for rs in root_seg
            if root_in_cell(rs, sc, tol=tol)
                A, rootDW_g = root_geom_area_mass(rs, vars)
                rs.potNup = EN * (Ihat * A + Ilat * rs.biomass)
                pot_per_cell[si] += rs.potNup
            end
        end

        for rs in root_seg_lateral
            if root_in_cell(rs, sc, tol=tol)
                A, rootDW_g = root_geom_area_mass(rs, vars)
                base_pot = EN * (Ihat * A + Ilat * rs.biomass)
                d_fine = vars.fineroot_diameter
                fr_surface_area = 4.0 * rs.finerootbiomass / (vars.RTD * d_fine + 1e-12)
                HAT_from_fine = Ihat * fr_surface_area
                LAT_from_fine = vars.K2 * rs.finerootbiomass
                rs.potNup = base_pot + EN * (HAT_from_fine + LAT_from_fine)
                pot_per_cell[si] += rs.potNup
            end
        end
    end

    # Phase B: allocate supply
    for (si, sc) in enumerate(soil_cells)
        available = soil_available_umol(sc, vars)
        total_pot = pot_per_cell[si]
        if total_pot <= 0.0 || available <= 0.0; continue; end

        for rs in root_seg
            if root_in_cell(rs, sc, tol=tol)
                fair = available * (rs.potNup / total_pot)
                actual = min(rs.potNup, fair)
                rs.actNup = actual; sc.N -= actual; total_plant_N_umol += actual
            end
        end
        for rs in root_seg_lateral
              if root_in_cell(rs, sc, tol=tol)      # ← fixed, was inline
                fair = available * (rs.potNup / total_pot)
                actual = min(rs.potNup, fair)
                rs.actNup = actual; sc.N -= actual; total_plant_N_umol += actual
            end
        end
        sc.N = max(0.0, sc.N)
    end

    total_plant_N_g = total_plant_N_umol * 14e-6
    vars.N_uptake = total_plant_N_g
    vars.N_uptake_cumu += total_plant_N_g
    return total_plant_N_g
end

# ═══════════════════════════════════════════════════════════════════
# SECTION: RAY TRACING & LIGHT
# ═══════════════════════════════════════════════════════════════════

# Global soil_surface placeholder (will be set by field setup )
soil_surface = nothing

function create_scene(forest)
    mesh = Mesh(vec(forest))
    if Main.soil_surface !== nothing
        add!(mesh, Main.soil_surface, materials = Lambertian(τ = 0.0, ρ = 0.21))
    end
    return mesh
end

function create_sky(; mesh, lat = 52.0 * π / 180.0, DOY = 182)
    fs = collect(0.1:0.1:0.9)
    dec = declination(DOY)
    DL = day_length(lat, dec) * 3600
    temp = [clear_sky(lat = lat, DOY = DOY, f = f) for f in fs]
    Idir = getindex.(temp, 2)
    f_dir = waveband_conversion(Itype = :direct, waveband = :PAR, mode = :power)
    f_dif = waveband_conversion(Itype = :diffuse, waveband = :PAR, mode = :power)
    Idir_PAR = f_dir .* Idir
    Idif_PAR = f_dif .* getindex.(temp, 3)
    dome = sky(mesh,
              Idir = 0.0, Idif = sum(Idir_PAR) / 10 * DL,
              nrays_dif = 1_000_000, sky_model = StandardSky,
              dome_method = equal_solid_angles, ntheta = 9, nphi = 12)
    for I in Idir_PAR
        push!(dome, sky(mesh, Idir = I / 10 * DL, nrays_dir = 100_000, Idif = 0.0)[1])
    end
    return dome
end

function create_raytracer(acc_mesh, sources, sl::Float64, sw::Float64)
    settings = RTSettings(pkill = 0.9, maxiter = 4, nx = 0, ny = 0, dx = sl, dy = sw+0.02, parallel = true) ##
    RayTracer(acc_mesh, sources, settings = settings)
end

# settings = RTSettings(pkill = 0.9, maxiter = 4, nx = 5, ny = 5,
#                        rdx = soil_length, rdy = soil_width, parallel = true)



function run_raytracer!(forest, sl::Float64, sw::Float64; DOY = 182) # sl::Float64, sw::Float64;
    mesh = create_scene(forest)
      if ntriangles(mesh) == 0
        return nothing
    end
    acc_mesh = accelerate(mesh, acceleration = BVH, rule = SAH{3}(5, 10))
    sources = create_sky(mesh = acc_mesh, DOY = DOY)
    rtobj = create_raytracer(acc_mesh, sources,sl,sw) 
    trace!(rtobj)
    return nothing
end

function reset_PAR!(forest)
    for tree in forest; data(tree).PAR = 0.0; end
    return nothing
end

function calculate_PAR!(forest, sl::Float64, sw::Float64; DOY = 182)
    reset_PAR!(forest)
    run_raytracer!(forest, sl, sw, DOY = DOY)
    @threads for tree in forest
        for l in get_leaves(tree)
            data(tree).PAR += power(l.material)[1]
        end
    end
    return nothing
end


# ═══════════════════════════════════════════════════════════════════
# SECTION: TEMPERATURE FUNCTION
# ═══════════════════════════════════════════════════════════════════

"""
    temperature(dayofyear)

Return the mean daily air temperature (°C) for a given day of year using a
sinusoidal approximation representative of a temperate climate:

    T(DOY) = 10.7 + 7.55 × sin(2π × (DOY − 111) / 365)

Parameters are calibrated for a Northern European / central climate zone
(mean annual temperature ≈ 10.7 °C, seasonal amplitude ≈ 7.55 °C, with the
temperature curve peaking around DOY 111 + 91 ≈ DOY 202, i.e., mid-July).
Adapt these constants for other climates.
"""
function temperature(dayofyear)
    tav = 10.7 + 7.55 * sin(2 * π * (dayofyear - 111) / 365)
    return tav
end

# # print the curve
# using Plots
# # Generate temperatures for all days of the year
# using Plots
# # Generate temperatures for all days of the year
# days = 1:365
# temps = [temperature(day) for day in days]
# # Plot the temperature curve
# plot(days, temps, label = "Temperature (°C)", xlabel = "Day of Year", ylabel = "Temperature (°C)", title = "Temperature Curve for a Year")


# ═══════════════════════════════════════════════════════════════════
# SECTION: DAILY STEP
# ═══════════════════════════════════════════════════════════════════

function daily_step!(forest, soil, soil_surface_arg, DOY, x_min, x_max, y_min, y_max, sl, sw)
    calculate_PAR!(forest, sl, sw, DOY = DOY)

    for tree in forest
        nitrogen_uptake!(tree, soil)
    end

    @threads for tree in forest
        all_leaves = get_leaves(tree)
        all_petioles = get_petioles(tree)
        all_internodes = get_internodes(tree)
        all_meristems = get_meristems(tree)
        all_Rsegments = get_Rsegments(tree)
        all_Rmeristems = get_Rmeristems(tree)
        all_RsegmentsLateral = get_RsegmentsLateral(tree)
        all_RmeristemsLateral = get_RmeristemsLateral(tree)
        all_Rsystems = get_Rsystems(tree)

        tvars = data(tree)
        age!(all_leaves, all_petioles, all_internodes,
             all_meristems, all_Rmeristems, all_Rsegments,
             all_RsegmentsLateral, all_RmeristemsLateral, all_Rsystems,
             tvars, DOY)

        grow!(tree, all_leaves, all_petioles, all_internodes, all_Rsegments, all_RsegmentsLateral)
        tvars = data(tree)
        size_leaves!(all_leaves, tvars)
        size_petioles!(all_petioles, tvars)
        size_internodes!(all_internodes, tvars)
        size_Rsegments!(all_Rsegments, tvars)
        size_RsegmentsLateral!(all_RsegmentsLateral, tvars)
        rewrite!(tree)
    end
    
end

# ═══════════════════════════════════════════════════════════════════
# SECTION: RENDERING
# ═══════════════════════════════════════════════════════════════════

function render_forest(forest, soil_surface_mesh, soil_mesh)
    tree_mesh = Mesh(vec(forest))
    all_meshes = [tree_mesh, soil_surface_mesh, soil_mesh]
    mesh = Mesh(all_meshes)
    fig = render(mesh)
    return fig
end
# ═══════════════════════════════════════════════════════════════════
# SECTION: PLANT CONSTRUCTOR
# ═══════════════════════════════════════════════════════════════════

function create_tree(origin, orientation, species::Symbol)
    vars = SpeciesParams.get_species_params(species)

    # Initialize organ dimensions based on initial biomass
    leaf_length, leaf_width = leaf_dims(vars.LB0, vars)
    vleaf = (biomass = vars.LB0, length = leaf_length, width = leaf_width)
    pet_length, pet_width = pet_dims(vars.IB0, vars)
    vpet = (biomass = vars.IB0, length = pet_length, width = pet_width)
    int_length, int_width = int_dims(vars.IB0, vars, 0.0)
    vint = (biomass = vars.IB0, length = int_length, width = int_width)
    Rseg_length, Rseg_width = Rseg_dims(vars.IB0, vars)
    vrseg = (biomass = vars.IB0, length = Rseg_length, width = vars.root_initial_diameter)

   vrseg_lateral = (biomass = vrseg.biomass * vars.root_lateral_scale,
                 length  = vrseg.length  * vars.root_lateral_scale,
                 width   = vars.root_initial_diameter * vars.root_lateral_scale)

    meristem_rule = create_meristem_rule(vleaf, vint, vpet)
    branch_rule = create_branch_rule(vint)
    root_rule = create_root_meristem_rule()
    root_node_to_lateral_rule = create_root_node_to_lateral_rule(vrseg_lateral)
    lateral_root_rule = create_lateral_root_meristem_rule()
    rsystem_emergence_rule = create_rsystem_emergence_rule()

    # NEW: only the seminal is in the axiom
# Build the seminal axis (just the first root)
seminal_az   = rand() * 360.0
seminal_tilt = vars.seminal_tilt + (rand() * 2 - 1) * vars.seminal_tilt_jitter
vars.root_local_cumu = Dict(1 => 0.0)
vars.nodal_count    = 0 

seminal_axis = RH(seminal_az) + RA(-seminal_tilt) +
               TreeTypes.Rmeristem(meristem_index = 1)

# Rsystem now LIVES on the root side, after the RA(180) flip,
# and is the parent of the seminal (and of all future nodals).
axiom = T(origin) + RH(orientation) + (
    TreeTypes.Meristem(),
    RA(180.0) + TreeTypes.Rsystem() + seminal_axis
)

tree = Graph(axiom = axiom,
             rules = (meristem_rule, branch_rule,
                      root_rule, root_node_to_lateral_rule, lateral_root_rule,
                      rsystem_emergence_rule),
             data = vars)
    data(tree).species = species
    return tree
end

# ═══════════════════════════════════════════════════════════════════
# SECTION: SOIL CELL CENTER FIX
# ═══════════════════════════════════════════════════════════════════

function fix_soilcell_centers!(soil)
    for sc in get_soilcells(soil)
        sc.z = sc.z + sc.height / 2.0
    end
    return nothing
end

function plant_N_percent_above(tree; eps = 1e-12)
    tvars = data(tree)
    above = max(total_above_biomass(tree), eps)
    return 100.0 * (tvars.N_uptake_cumu / above)
end

# ═══════════════════════════════════════════════════════════════════
# SECTION: UNIFIED DIAGNOSTICS COLLECTOR (NEW)
# ═══════════════════════════════════════════════════════════════════

"""
    DiagnosticTimeSeries

Single unified container for all model outputs needed for visualization.
Replaces all ad-hoc Dicts, record_daily!, and run_with_seed_tracking storage.
"""
Base.@kwdef mutable struct DiagnosticTimeSeries
    # ── Time axis ──
    day::Vector{Int}            = Int[]

    # ── Biomass dynamics ──
    total_biomass::Vector{Float64}    = Float64[]   # g
    above_biomass::Vector{Float64}    = Float64[]   # g
    below_biomass::Vector{Float64}    = Float64[]   # g
    daily_increment::Vector{Float64}  = Float64[]   # g/day (ΔB_actual from grow!)
    db_from_seed::Vector{Float64}     = Float64[]   # g/day seed reserve contribution
    db_from_photo::Vector{Float64}    = Float64[]   # g/day photosynthetic contribution

    # ── Nitrogen dynamics ──
    Cn::Vector{Float64}             = Float64[]   # plant N conc (g N / g above DW)
    Cn_crit::Vector{Float64}        = Float64[]   # critical N conc (g N / g above DW)
    N_uptake_daily::Vector{Float64} = Float64[]   # g N/day
    N_uptake_cumu::Vector{Float64}  = Float64[]   # g N cumulative
    soil_N_umol::Vector{Float64}    = Float64[]   # µmol total soil N

    # ── Stress indicator ──
    fn::Vector{Float64}             = Float64[]   # nitrogen stress factor [0,1]

    # ── Structural counts ──
    n_leaves::Vector{Int}           = Int[]
    n_internodes::Vector{Int}       = Int[]
    n_roots_main::Vector{Int}       = Int[]
    n_roots_lat::Vector{Int}        = Int[]
    leaf_area::Vector{Float64}      = Float64[]   # m²
    root_length_main::Vector{Float64} = Float64[] # m
    root_length_lat::Vector{Float64}  = Float64[] # m
    root_length_total::Vector{Float64}= Float64[] # m

    # ── Carbon budget (NEW) ──
    PAR_intercepted::Vector{Float64}       = Float64[]   # MJ/day (tvars.PAR)
    carbon_supply::Vector{Float64}         = Float64[]   # g/day  (ΔB available to allocate, including reserve)
    carbon_demand_total::Vector{Float64}   = Float64[]   # g/day (total_demand from grow!)
    carbon_demand_leaves::Vector{Float64}  = Float64[]   # g/day (sum of leaf demands)
    carbon_demand_petioles::Vector{Float64} = Float64[]  # g/day
    carbon_demand_internodes::Vector{Float64} = Float64[] # g/day
    carbon_demand_roots::Vector{Float64}   = Float64[]   # g/day (root system sink strength)
    carbon_surplus::Vector{Float64}        = Float64[]   # g/day (ΔB - ΔB_used, what goes to reserve)
    root_alloc_to_elongation::Vector{Float64} = Float64[] # g/day (root_alloc_total passed to compute_elongation!)
    carbon_alloc_leaves::Vector{Float64}    = Float64[]   # g/day gross allocation to leaves (before respiration)
    carbon_alloc_petioles::Vector{Float64}  = Float64[]   # g/day gross allocation to petioles (before respiration)
    carbon_alloc_internodes::Vector{Float64}= Float64[]   # g/day gross allocation to internodes (before respiration)
    carbon_alloc_roots::Vector{Float64}     = Float64[]   # g/day gross allocation to roots (before respiration)

    # ── Metadata ──
    label::String               = ""
    seed_reserve_days::Int      = 10
    seed_daily_carbon::Float64  = 0.003635
end

"""
    record!(diag, day, tree, soil)

Record all diagnostic variables for one simulation day.
Call AFTER grow!, size!, rewrite! have been applied for this day.
"""
function record!(diag::DiagnosticTimeSeries, day::Int, tree, soil)
    tvars = data(tree)
    eps = 1e-12

    # ── Organ queries ──
    all_l  = get_leaves(tree)
    all_p  = get_petioles(tree)
    all_i  = get_internodes(tree)
    all_rs = get_Rsegments(tree)
    all_rl = get_RsegmentsLateral(tree)

    # ── Biomass ──
    ab = sum(l.biomass for l in all_l; init=0.0) +
         sum(ii.biomass for ii in all_i; init=0.0) +
         sum(pp.biomass for pp in all_p; init=0.0)
    bl = sum(r.biomass for r in all_rs; init=0.0) +
         sum(r.biomass for r in all_rl; init=0.0)
    total_bm = ab + bl

    # ── Plant age for seed reserve calculation ──
    all_rsys = get_Rsystems(tree)
    plant_age = length(all_rsys) > 0 ? all_rsys[1].age : day

    # ── Seed vs photosynthetic carbon (mirrors grow! logic exactly) ──
    sd = diag.seed_reserve_days
    sc_daily = diag.seed_daily_carbon
    if plant_age <= sd
        seed_frac = max(0.0, 1.0 - plant_age / sd)
        db_seed = sc_daily * seed_frac
    else
        db_seed = 0.0
    end
    # ΔB_actual = ΔB_photo + ΔB_seed  →  ΔB_photo = ΔB_actual - ΔB_seed
    db_photo = max(0.0, tvars.ΔB_actual - db_seed)

    # ── Nitrogen concentrations ──
      safe_above = max(ab, eps)
    Cn_val = tvars.N_uptake_cumu / safe_above
    Cn_crit_raw = (3.68 * safe_above^(-0.28)) / 100.0
    Cn_crit_val = min(Cn_crit_raw, 0.06)   # ← clamp to 6% max
    fn_val = clamp(Cn_val / max(Cn_crit_val, eps), 0.0, 1.0)

    # ── Soil N ──
    soil_N_total = soil !== nothing ? total_soil_N(soil) : 0.0

    # ── Leaf area (m²) ──
    la = sum(l.length * l.width * π / 4.0 for l in all_l; init=0.0)

    # ── Root lengths (m) ──
    rl_main = sum(r.length for r in all_rs; init=0.0)
    rl_lat  = sum(r.length for r in all_rl; init=0.0)

    # ── Push ──
    push!(diag.day, day)
    push!(diag.total_biomass, total_bm)
    push!(diag.above_biomass, ab)
    push!(diag.below_biomass, bl)
    push!(diag.daily_increment, tvars.ΔB_actual)
    push!(diag.db_from_seed, db_seed)
    push!(diag.db_from_photo, db_photo)
    push!(diag.Cn, Cn_val)
    push!(diag.Cn_crit, Cn_crit_val)
    push!(diag.N_uptake_daily, tvars.N_uptake)
    push!(diag.N_uptake_cumu, tvars.N_uptake_cumu)
    push!(diag.soil_N_umol, soil_N_total)
    push!(diag.fn, fn_val)
    push!(diag.n_leaves, length(all_l))
    push!(diag.n_internodes, length(all_i))
    push!(diag.n_roots_main, length(all_rs))
    push!(diag.n_roots_lat, length(all_rl))
    push!(diag.leaf_area, la)
    push!(diag.root_length_main, rl_main)
    push!(diag.root_length_lat, rl_lat)
    push!(diag.root_length_total, rl_main + rl_lat)
    push!(diag.PAR_intercepted, tvars.PAR)
    push!(diag.carbon_supply, tvars._diag_carbon_supply)
    push!(diag.carbon_demand_total, tvars._diag_total_demand)
    push!(diag.carbon_demand_leaves, tvars._diag_leaf_demand)
    push!(diag.carbon_demand_petioles, tvars._diag_petiole_demand)
    push!(diag.carbon_demand_internodes, tvars._diag_internode_demand)
    push!(diag.carbon_demand_roots, tvars._diag_root_demand)
    push!(diag.carbon_surplus, tvars._diag_carbon_surplus)
    push!(diag.root_alloc_to_elongation, tvars._diag_root_alloc)
    push!(diag.carbon_alloc_leaves, tvars._diag_alloc_leaves)
    push!(diag.carbon_alloc_petioles, tvars._diag_alloc_petioles)
    push!(diag.carbon_alloc_internodes, tvars._diag_alloc_internodes)
    push!(diag.carbon_alloc_roots, tvars._diag_alloc_roots)

    return nothing
end

# ═══════════════════════════════════════════════════════════════════
# SECTION: VISUALIZATION PIPELINE (NEW)
# ═══════════════════════════════════════════════════════════════════

"""
    plot_diagnostics(diag; title_prefix="")

Standard 6-panel diagnostic figure for a single scenario:
  Row 1: Biomass accumulation  |  Growth sources (seed vs photo stacked area)
  Row 2: Cn vs Cn_crit         |  N uptake (daily + cumulative)
  Row 3: fn over time          |  Structural development
"""
function plot_diagnostics(diag::DiagnosticTimeSeries; title_prefix = "")
    d = diag.day
    pfx = isempty(title_prefix) ? "" : "$title_prefix — "
    seed_end = diag.seed_reserve_days

    # ── Panel 1: Biomass accumulation ──
    p1 = plot(d, diag.total_biomass .* 1000, lw = 2.5, label = "Total",
              xlabel = "Day", ylabel = "Biomass (mg)",
              title = "$(pfx)Biomass accumulation", legend = :topleft)
    plot!(p1, d, diag.above_biomass .* 1000, lw = 2, label = "Above", ls = :dash)
    plot!(p1, d, diag.below_biomass .* 1000, lw = 2, label = "Below", ls = :dash)
    if seed_end > 0 && seed_end <= d[end]
        vline!(p1, [seed_end], ls = :dot, color = :gray, lw = 1, alpha = 0.6,
               label = "Seed reserve ends")
    end

    # ── Panel 2: Growth sources (stacked area) ──
    seed_mg = diag.db_from_seed .* 1000
    photo_mg = diag.db_from_photo .* 1000
    total_mg = diag.daily_increment .* 1000

    p2 = plot(d, seed_mg, lw = 0, fillrange = 0, fillalpha = 0.5,
              color = :orange, label = "Seed reserve",
              xlabel = "Day", ylabel = "ΔB (mg/day)",
              title = "$(pfx)Growth sources", legend = :topright)
    plot!(p2, d, seed_mg .+ photo_mg, lw = 0,
          fillrange = seed_mg, fillalpha = 0.5,
          color = :green, label = "Photosynthesis")
    plot!(p2, d, total_mg, lw = 2, color = :black, label = "Total ΔB")
    if seed_end > 0 && seed_end <= d[end]
        vline!(p2, [seed_end], ls = :dot, color = :gray, lw = 1, alpha = 0.6, label = nothing)
    end

    # ── Panel 3: Cn vs Cn_crit ──
    p3 = plot(d, diag.Cn .* 100, lw = 2.5, color = :blue, label = "Cₙ (plant N%)",
              xlabel = "Day", ylabel = "N concentration (%)",
              title = "$(pfx)Nitrogen concentrations", legend = :topright, ylims = (0,6))

    plot!(p3, d, diag.Cn_crit .* 100, lw = 2.5, color = :red, ls = :dash,
          label = "Cₙ,crit (dilution curve)")
    # Mark where Cn crosses below Cn_crit
    deficit_idx = findfirst(i -> diag.Cn[i] < diag.Cn_crit[i], 1:length(d))
    if deficit_idx !== nothing
        vline!(p3, [d[deficit_idx]], ls = :dot, color = :red, lw = 1, alpha = 0.4, label = nothing)
        annotate!(p3, d[deficit_idx] + 1,
                  maximum(diag.Cn_crit .* 100) * 0.85,
                  text("← N stress", 7, :left, :red))
    end

    # ── Panel 4: N uptake ──
    p4 = plot(d, diag.N_uptake_daily .* 1e6, lw = 2, color = :teal,
              label = "Daily (µg/day)",
              xlabel = "Day", ylabel = "N uptake (µg/day)",
              title = "$(pfx)Nitrogen uptake", legend = :topleft)
    p4_twin = twinx(p4)
    plot!(p4_twin, d, diag.N_uptake_cumu .* 1000, lw = 2, color = :purple,
          label = "Cumul. (mg)", ylabel = "Cumul. N (mg)", ls = :dash,
          legend = :right)

    # ── Panel 5: fn ──
    p5 = plot(d, diag.fn, lw = 2.5, color = :darkred,
              xlabel = "Day", ylabel = "fₙ",
              title = "$(pfx)Nitrogen stress factor (fₙ)",
              ylims = (-0.05, 1.1), legend = false)
    hline!(p5, [1.0], ls = :dash, color = :green, lw = 1, alpha = 0.5)
    hline!(p5, [0.5], ls = :dot, color = :orange, lw = 1, alpha = 0.3)
    annotate!(p5, d[end] * 0.65, 1.05, text("No stress", 8, :green))
    annotate!(p5, d[end] * 0.65, 0.45, text("Moderate", 8, :orange))
    if seed_end > 0 && seed_end <= d[end]
        vline!(p5, [seed_end], ls = :dot, color = :gray, lw = 1, alpha = 0.6)
    end

    # ── Panel 6: Structural development ──
    p6 = plot(d, diag.n_leaves, lw = 2, label = "Leaves",
              xlabel = "Day", ylabel = "Count",
              title = "$(pfx)Structure", legend = :topleft)
    plot!(p6, d, diag.n_roots_main .+ diag.n_roots_lat, lw = 2, label = "Root segs")
    p6_twin = twinx(p6)
    plot!(p6_twin, d, diag.leaf_area .* 1e4, lw = 2, color = :green, ls = :dash,
          label = "Leaf area (cm²)", ylabel = "Leaf area (cm²)", legend = :right)

    fig = plot(p1, p2, p3, p4, p5, p6,
               layout = (3, 2), size = (1200, 1000), margin = 5Plots.mm)
    return fig
end

"""
    plot_diagnostics_comparison(diags; title="Comparison")

Overlay multiple DiagnosticTimeSeries on 6 shared panels for direct comparison.
"""
function plot_diagnostics_comparison(diags::Vector{DiagnosticTimeSeries};
                                     title = "Comparison")
    palette = [:blue, :orange, :red, :green, :purple, :brown, :cyan, :magenta]

    p1 = plot(xlabel = "Day", ylabel = "Biomass (mg)",
              title = "$title — Biomass", legend = :topleft)
    p2 = plot(xlabel = "Day", ylabel = "ΔB (mg/day)",
              title = "$title — Growth rate", legend = :topright)
    p3 = plot(xlabel = "Day", ylabel = "N conc. (%)",
              title = "$title — Cₙ vs Cₙ,crit", legend = :topright, ylims = (0,6))
    p4 = plot(xlabel = "Day", ylabel = "fₙ",
              title = "$title — Stress factor fₙ", ylims = (-0.05, 1.1))
    p5 = plot(xlabel = "Day", ylabel = "Cumul. N (mg)",
              title = "$title — N uptake", legend = :topleft)
    p6 = plot(xlabel = "Day", ylabel = "Soil N (µmol)",
              title = "$title — Soil N", legend = :topright)

    for (idx, diag) in enumerate(diags)
        c = palette[mod1(idx, length(palette))]
        lab = isempty(diag.label) ? "Run $idx" : diag.label
        d = diag.day

        # Panel 1: biomass
        plot!(p1, d, diag.total_biomass .* 1000, lw = 2.5, color = c, label = lab)

        # Panel 2: daily increment
        plot!(p2, d, diag.daily_increment .* 1000, lw = 2, color = c, label = lab)

        # Panel 3: Cn (solid) + Cn_crit (dashed, same color)
        plot!(p3, d, diag.Cn .* 100, lw = 2.5, color = c, label = "Cₙ $lab")
        plot!(p3, d, diag.Cn_crit .* 100, lw = 1.5, color = c, ls = :dash,
              label = "Cₙ,crit $lab")

        # Panel 4: fn
        plot!(p4, d, diag.fn, lw = 2.5, color = c, label = lab)

        # Panel 5: cumulative N uptake
        plot!(p5, d, diag.N_uptake_cumu .* 1000, lw = 2.5, color = c, label = lab)

        # Panel 6: soil N
        plot!(p6, d, diag.soil_N_umol, lw = 2.5, color = c, label = lab)
    end

    hline!(p4, [1.0], ls = :dash, color = :gray, lw = 0.5, label = nothing)

    fig = plot(p1, p2, p3, p4, p5, p6,
               layout = (3, 2), size = (1200, 1000), margin = 5Plots.mm)
    return fig
end

"""
    plot_carbon_budget(diag; title_prefix="")

10-panel diagnostic figure for carbon supply/demand analysis:
  Row 1: PAR intercepted              |  Carbon supply vs total demand
  Row 2: Organ demand breakdown       |  Actual allocation breakdown
  Row 3: Supply minus demand (surplus) |  Root alloc to elongation
  Row 4: Organ counts                 |  Leaf area
  Row 5: Demand fractions (%)         |  Allocation fractions (%)
"""
function plot_carbon_budget(diag::DiagnosticTimeSeries; title_prefix = "")
    d = diag.day
    pfx = isempty(title_prefix) ? "" : "$title_prefix — "

    # Panel 1: PAR intercepted
    p1 = plot(d, diag.PAR_intercepted, lw=2, label="PAR",
              xlabel="Day", ylabel="PAR (MJ/day?)",
              title="$(pfx)PAR Intercepted", legend=:topleft)

    # Panel 2: Supply vs Demand
    p2 = plot(d, diag.carbon_supply .* 1000, lw=2.5, label="Supply (ΔB)",
              xlabel="Day", ylabel="mg/day",
              title="$(pfx)Carbon Supply vs Demand", legend=:topleft)
    plot!(p2, d, diag.carbon_demand_total .* 1000, lw=2, ls=:dash, label="Total demand")

    # Panel 3: Organ demand breakdown
    p3 = plot(d, diag.carbon_demand_leaves .* 1000, lw=2, label="Leaves",
              xlabel="Day", ylabel="mg/day",
              title="$(pfx)Organ Demands", legend=:topleft)
    plot!(p3, d, diag.carbon_demand_petioles .* 1000, lw=2, label="Petioles")
    plot!(p3, d, diag.carbon_demand_internodes .* 1000, lw=2, label="Internodes")
    plot!(p3, d, diag.carbon_demand_roots .* 1000, lw=2, label="Root system")

    # Panel 4: Actual allocation breakdown
    p_alloc = plot(d, diag.carbon_alloc_leaves .* 1000, lw=2, label="Leaves",
                   xlabel="Day", ylabel="mg/day",
                   title="$(pfx)Actual Allocation", legend=:topleft)
    plot!(p_alloc, d, diag.carbon_alloc_petioles .* 1000, lw=2, label="Petioles")
    plot!(p_alloc, d, diag.carbon_alloc_internodes .* 1000, lw=2, label="Internodes")
    plot!(p_alloc, d, diag.carbon_alloc_roots .* 1000, lw=2, label="Roots")

    # Panel 5: Surplus (supply - demand)
    p4 = plot(d, diag.carbon_surplus .* 1000, lw=2, label="Surplus",
              xlabel="Day", ylabel="mg/day",
              title="$(pfx)Carbon Surplus (lost or reserved)", legend=:topleft,
              fillrange=0, fillalpha=0.3)

    # Panel 6: Root allocation that reaches elongation
    p5 = plot(d, diag.root_alloc_to_elongation .* 1000, lw=2, label="Root alloc → elongation",
              xlabel="Day", ylabel="mg/day",
              title="$(pfx)Root Alloc to Elongation", legend=:topleft)

    # Panel 7: Organ counts
    p6 = plot(d, diag.n_leaves, lw=2, label="Leaves",
              xlabel="Day", ylabel="Count",
              title="$(pfx)Organ Counts", legend=:topleft)
    plot!(p6, d, diag.n_roots_main, lw=2, label="Main root segs")
    plot!(p6, d, diag.n_roots_lat, lw=2, label="Lateral root segs")

    # Panel 8: Leaf area
    p7 = plot(d, diag.leaf_area .* 1e4, lw=2, label="Leaf area",
              xlabel="Day", ylabel="cm²",
              title="$(pfx)Leaf Area", legend=:topleft)

    # Panel 9: Demand fractions
    total_d = max.(diag.carbon_demand_total, 1e-12)
    p8 = plot(d, 100 .* diag.carbon_demand_leaves ./ total_d, lw=2, label="Leaves %",
              xlabel="Day", ylabel="% of total demand",
              title="$(pfx)Demand Fractions", legend=:topright, ylims=(0, 105))
    plot!(p8, d, 100 .* diag.carbon_demand_petioles ./ total_d, lw=2, label="Petioles %")
    plot!(p8, d, 100 .* diag.carbon_demand_internodes ./ total_d, lw=2, label="Internodes %")
    plot!(p8, d, 100 .* diag.carbon_demand_roots ./ total_d, lw=2, label="Roots %")

    # Panel 10: Allocation fractions
    total_alloc = diag.carbon_alloc_leaves .+ diag.carbon_alloc_petioles .+
                  diag.carbon_alloc_internodes .+ diag.carbon_alloc_roots
    total_alloc_safe = max.(total_alloc, 1e-12)
    p_alloc_frac = plot(d, 100 .* diag.carbon_alloc_leaves ./ total_alloc_safe, lw=2, label="Leaves %",
                        xlabel="Day", ylabel="% of total allocation",
                        title="$(pfx)Allocation Fractions", legend=:topright, ylims=(0, 105))
    plot!(p_alloc_frac, d, 100 .* diag.carbon_alloc_petioles ./ total_alloc_safe, lw=2, label="Petioles %")
    plot!(p_alloc_frac, d, 100 .* diag.carbon_alloc_internodes ./ total_alloc_safe, lw=2, label="Internodes %")
    plot!(p_alloc_frac, d, 100 .* diag.carbon_alloc_roots ./ total_alloc_safe, lw=2, label="Roots %")

    fig = plot(p1, p2, p3, p_alloc, p6, p7,
               layout=(3, 2), size=(1400, 1750), margin=5Plots.mm)
    return fig
end

# ═══════════════════════════════════════════════════════════════════
# SECTION: SIMPLIFIED RUNNER (single voxel, fixed PAR, no ray tracer)
# ═══════════════════════════════════════════════════════════════════

"""
    run_simplified(; kwargs...)

Run a single plant with fixed daily PAR and a single soil voxel.
Returns a DiagnosticTimeSeries ready for plotting.
"""
# In run_simplified, AFTER creating the tree, set the seed params on the plant:

# function run_simplified(; n_days::Int = 80,
#                           soil_N_init::Float64 = 10000.0,
#                           fixed_PAR::Float64 = 3.0e4,
#                           label::String = "",
#                           EN_enabled::Bool = false,
#                           seed_daily_carbon::Float64 = 0.003635,
#                           seed_reserve_days::Int = 10,
#                           seed_nitrogen::Float64 = 0.00015, ### 1 mg seed N
#                           sowing_date::Int = 90)
#     Random.seed!(123456789)

#     tree = create_tree(Vec(0.0, 0.0, 0.0), 0.0, :cereal2)
#     tvars = data(tree)
#     tvars.EN_enabled = EN_enabled
#     # ── NEW: Actually set seed params on the plant ──
#     tvars.seed_daily_carbon = seed_daily_carbon
#     tvars.seed_reserve_days = seed_reserve_days
#     tvars.seed_nitrogen = seed_nitrogen
#     tvars.sowing_date = sowing_date

#     # ... rest unchanged ...
#     soil = Graph(axiom =
#         RA(0.0) + T(Vec(0.0, 0.0, -0.5)) +
#         Soilcell(length = 0.5, width = 0.5, height = 0.5,
#                  i = 1, j = 1, k = 1, N = soil_N_init)
#     )
#     Mesh(soil)
#     fix_soilcell_centers!(soil)

# # Add this right after creating the tree and setting params in run_simplified:
# println("  DEBUG: EN_enabled = ", tvars.EN_enabled,
#         "  EN_coeff = ", tvars.EN_coeff)

#     diag = DiagnosticTimeSeries(label = label,
#                                 seed_reserve_days = seed_reserve_days,
#                                 seed_daily_carbon = seed_daily_carbon)

#     for day in 1:n_days
#         tvars = data(tree)
#         tvars.PAR = fixed_PAR

#         all_l   = get_leaves(tree);   all_p   = get_petioles(tree)
#         all_i   = get_internodes(tree); all_m = get_meristems(tree)
#         all_rs  = get_Rsegments(tree);  all_rm = get_Rmeristems(tree)
#         all_rl  = get_RsegmentsLateral(tree)
#         all_rlm = get_RmeristemsLateral(tree)
#         all_rsys = get_Rsystems(tree)

#         dayofyear = tvars.sowing_date + day
#         age!(all_l, all_p, all_i, all_m, all_rm, all_rs, all_rl, all_rlm, all_rsys, tvars, dayofyear)
#         Mesh(tree)
#         nitrogen_uptake!(tree, soil)

#         all_l  = get_leaves(tree);  all_p = get_petioles(tree)
#         all_i  = get_internodes(tree)
#         all_rs = get_Rsegments(tree); all_rl = get_RsegmentsLateral(tree)
#         grow!(tree, all_l, all_p, all_i, all_rs, all_rl)

#         tvars = data(tree)
#         size_leaves!(all_l, tvars)
#         size_petioles!(all_p, tvars)
#         size_internodes!(all_i, tvars)
#         size_Rsegments!(all_rs, tvars)
#         size_RsegmentsLateral!(all_rl, tvars)
#         rewrite!(tree)

#         record!(diag, day, tree, soil)
#     end

#     return diag, tree, soil
# end

# ═══════════════════════════════════════════════════════════════════
# SECTION: FULL MODEL RUNNER (ray tracer + multi-voxel soil)
# ═══════════════════════════════════════════════════════════════════

# ── Run meta-data writer ──
function write_run_metadata(path::String;
        label, sp1, sp2,
        n_plants_row_sp1, plant_spacing_sp1, n_rows_sp1, row_spacing_sp1,
        n_plants_row_sp2, plant_spacing_sp2, n_rows_sp2, row_spacing_sp2,
        soil_N_init, sv, n_i, n_j, n_k,
        n_plants_sp1, n_plants_sp2, EN_coeff)

    rue_sp1 = Main.SpeciesParams.get_species_params(sp1).RUE
    rue_sp2 = Main.SpeciesParams.get_species_params(sp2).RUE

    df = DataFrame(
        label              = label,
        sp1                = String(sp1),
        sp2                = String(sp2),
        n_plants_row_sp1   = n_plants_row_sp1,
        plant_spacing_sp1  = plant_spacing_sp1,
        n_rows_sp1         = n_rows_sp1,
        row_spacing_sp1    = row_spacing_sp1,
        n_plants_row_sp2   = n_plants_row_sp2,
        plant_spacing_sp2  = plant_spacing_sp2,
        n_rows_sp2         = n_rows_sp2,
        row_spacing_sp2    = row_spacing_sp2,
        n_plants_sp1       = n_plants_sp1,
        n_plants_sp2       = n_plants_sp2,
        n_plants_total     = n_plants_sp1 + n_plants_sp2,
        n_rows_total       = n_rows_sp1 + n_rows_sp2,
        soil_N_init        = soil_N_init,
        sv                 = sv,
        n_voxels_x         = n_i,
        n_voxels_y         = n_j,
        n_voxels_z         = n_k,
        RUE_sp1            = rue_sp1,
        RUE_sp2            = rue_sp2,
        EN_coeff           = EN_coeff,
    )
    CSV.write(path, df)
    return df
end

# collector
function root_voxel_occupancy(forest, soil, plant_ids, DOY::Int)

    soil_cells = get_soilcells(soil)

    # (cell_index, plant_id, seg_type) => (length, biomass)
    acc  = Dict{Tuple{Int,Int,Symbol}, NTuple{2,Float64}}()
    cellmeta = Dict{Int, NTuple{6,Float64}}()

    for (si, sc) in enumerate(soil_cells)
        cellmeta[si] = (Float64(sc.i), Float64(sc.j), Float64(sc.k), sc.x, sc.y, sc.z)
    end

    for (tree, pid) in zip(forest, plant_ids)
        vars = data(tree)

        # ── regular root segments ──
        for rs in get_Rsegments(tree)
            _, dw = root_geom_area_mass(rs, vars)
            for (si, sc) in enumerate(soil_cells)
                if root_in_cell(rs, sc)
                    L, B = get(acc, (si, pid, :normal), (0.0, 0.0))
                    acc[(si, pid, :normal)] = (L + rs.length, B + dw)
                    break
                end
            end
        end

        # ── lateral root segments ──
        for rs in get_RsegmentsLateral(tree)
            _, dw = root_geom_area_mass(rs, vars)
            for (si, sc) in enumerate(soil_cells)
                if root_in_cell(rs, sc)
                    L, B = get(acc, (si, pid, :lateral), (0.0, 0.0))
                    acc[(si, pid, :lateral)] = (L + rs.length, B + dw)
                    break
                end
            end
        end
    end

    # per-(voxel, plant, seg_type) long table
    pp = DataFrame(DOY=Int[], cell=Int[], i=Int[], j=Int[], k=Int[], x=Float64[], y=Float64[], z=Float64[],
                   plant_id=Int[], seg_type=Symbol[], root_length=Float64[], root_biomass=Float64[])
    for ((si, pid, stype), (L, B)) in acc
        (i,j,k,x,y,z) = cellmeta[si]
        push!(pp, (DOY, si, Int(i), Int(j), Int(k), x, y, z, pid, stype, L, B))
    end

    # per-voxel summary (still split by seg_type, aggregated across plants)
    vox = Dict{Tuple{Int,Symbol}, NTuple{2,Float64}}()
    contributors = Dict{Tuple{Int,Symbol}, Set{Int}}()
    for ((si, pid, stype), (L, B)) in acc
        key = (si, stype)
        Lt, Bt = get(vox, key, (0.0, 0.0)); vox[key] = (Lt + L, Bt + B)
        push!(get!(contributors, key, Set{Int}()), pid)
    end
    pv = DataFrame(DOY=Int[], cell=Int[], i=Int[], j=Int[], k=Int[], x=Float64[], y=Float64[], z=Float64[],
                   seg_type=Symbol[], total_root_length=Float64[], total_root_biomass=Float64[],
                   n_plants=Int[], plant_ids=String[])
    for ((si, stype), (Lt, Bt)) in vox
        (i,j,k,x,y,z) = cellmeta[si]
        ids = sort(collect(contributors[(si, stype)]))
        push!(pv, (DOY, si, Int(i), Int(j), Int(k), x, y, z, stype, Lt, Bt,
                   length(ids), join(ids, ";")))
    end

    return pp, pv
end


"""
    run_full_model(; kwargs...)

Run plant(s) with full ray tracing and multi-voxel soil.
Returns a DiagnosticTimeSeries for plant 1, plus the forest and soil.
"""
function run_full_model(; n_days::Int = 80,
                          start_DOY::Int = 180,
                          soil_N_init::Float64 = 30000.0, # umol per liter 
                          label::String = "Full model",
                          EN_enabled::Bool = true,
                          sv:: Float64 = 0.01,      # soil voxel dimension (m)
                          sowing_date::Int = 90,
                          #  Field design parameters 
                          # ── Species 1 ── #
                          sp1::Symbol = :cereal1,
                          n_plants_row_sp1::Int = 1,
                          plant_spacing_sp1::Float64 = 0.15,
                          n_rows_sp1::Int = 1,
                          row_spacing_sp1::Float64 = 0.25,
                          # ── Species 2 ──
                          sp2::Symbol = :cereal2,
                          n_plants_row_sp2::Int = 1,
                          plant_spacing_sp2::Float64 = 0.02,
                          n_rows_sp2::Int = 1,
                          row_spacing_sp2::Float64 = 0.12,
                          # ── Inter-species gap ──
                          inter_row_spacing::Float64 = 0.30)
    #Random.seed!(123456789)

    # ── Build positions /field designs  ──
    positions = Tuple{Vec{Float64}, Symbol, Int, Int}[]
    plant_id = 0
    y_cursor = 0.0

    # rows of sp1
    for row in 1:n_rows_sp1
        y = y_cursor + (row - 1) * row_spacing_sp1
        for col in 1:n_plants_row_sp1
            plant_id += 1
            x = (col - (n_plants_row_sp1 + 1) / 2) * plant_spacing_sp1
            push!(positions, (Vec(x, y, 0.0), sp1, 1, plant_id))
        end
    end

    # gap between species (only if both present)
    if n_rows_sp1 > 0 && n_rows_sp2 > 0
        y_cursor += (n_rows_sp1 - 1) * row_spacing_sp1 + inter_row_spacing
    elseif n_rows_sp1 > 0
        y_cursor += (n_rows_sp1 - 1) * row_spacing_sp1
    end

    # rows of sp2
    for row in 1:n_rows_sp2
        y = y_cursor + (row - 1) * row_spacing_sp2
        for col in 1:n_plants_row_sp2
            plant_id += 1
            x = (col - (n_plants_row_sp2 + 1) / 2) * plant_spacing_sp2
            push!(positions, (Vec(x, y, 0.0), sp2, 1, plant_id))
        end
    end
    #
    
    xs     = [Float64(p[1][1]) for p in positions]
    ys     = [Float64(p[1][2]) for p in positions]
    buffer = 0.0
# Example calls for different congigurations 
#     # Monoculture cereal2 (2 rows x 4 plants)
# run_full_model(n_rows_sp1 = 0, sp2 = :cereal2,
#                n_plants_row_sp2 = 4, plant_spacing_sp2 = 0.15,
#                n_rows_sp2 = 2, row_spacing_sp2 = 0.20)

# # Row intercropping (1 row each)
# run_full_model(sp1 = :cereal1, n_plants_row_sp1 = 4, n_rows_sp1 = 1,
#                sp2 = :cereal2, n_plants_row_sp2 = 4, n_rows_sp2 = 1,
#                inter_row_spacing = 0.30)

# # Strip intercropping (3 rows sp1 then 3 rows sp2)
# run_full_model(sp1 = :cereal1, n_plants_row_sp1 = 4, n_rows_sp1 = 3,
#                sp2 = :cereal2, n_plants_row_sp2 = 4, n_rows_sp2 = 3,
#                inter_row_spacing = 0.50)
    

    # RUE_dist = Normal(1.6, 0.10)     
    forest_local = [create_tree(pos, rand() * 360.0, sp) 
                    for (pos, sp, _, _) in positions]
    
    # forest_local = [create_tree(pos, 0.0, sp)                   #### instead of pos, rand() * 360.0, sp we try pos,0.0,sp to remove rotation
    #                 for (pos, sp, _, _) in positions]

    for tree in forest_local
    data(tree).EN_enabled = EN_enabled
    data(tree).sowing_date = sowing_date
end

# X: n_plants × spacing gives half-margin on each side automatically
sl_sp1 = n_plants_row_sp1 > 0 ? n_plants_row_sp1 * plant_spacing_sp1 : 0.0
sl_sp2 = n_plants_row_sp2 > 0 ? n_plants_row_sp2 * plant_spacing_sp2 : 0.0
sl = max(0.1, max(sl_sp1, sl_sp2))

# Y: total row extents + inter-species gap
sw_sp1 = n_rows_sp1 * row_spacing_sp1
sw_sp2 = n_rows_sp2 * row_spacing_sp2
gap    = (n_rows_sp1 > 0 && n_rows_sp2 > 0) ? inter_row_spacing : 0.0
sw = max(0.1, sw_sp1 + gap + sw_sp2)

    x_center = (maximum(xs) + minimum(xs)) / 2.0
    y_center  = (maximum(ys) + minimum(ys)) / 2.0

    x_min = x_center - sl / 2.0
    x_max = x_center + sl / 2.0
    y_min = y_center - sw / 2.0
    y_max = y_center + sw / 2.0

       for tree in forest_local
        data(tree).soil_x_min = x_min
        data(tree).soil_x_max = x_max
        data(tree).soil_y_min = y_min
        data(tree).soil_y_max = y_max
    end

    # sv = 0.1    #### cell dimension (m)
    n_j = max(1, Int(floor(sl / sv)))
    n_i = max(1, Int(floor(sw / sv)))
    n_k = round(0.5 / sv)  # 0.5 m soil depth

    # convert soil n from liters to voxel volume 

    voxel_volume_L = sv^3 * 1000.0        # m³ → litres  (e.g. 0.1³ × 1000 = 1.0 L)
    N_per_voxel    = soil_N_init * voxel_volume_L

    soil_axiom = RA(-90.0) + T(Vec(0.0, 0.0, 0.0)) + Tuple(
        RA(-90.0) + T(Vec(
            x_center + sv * (j - n_j / 2 - 0.5),
            y_center + sv * (i - n_i / 2 - 0.5),
            -sv * (k - 1.0)
        )) + Soilcell(length = sv, width = sv, height = sv,
                      i = i, j = j, k = k, N = N_per_voxel)
        for j in 1:n_j, i in 1:n_i, k in 1:n_k
    )

    soil_local = Graph(axiom = soil_axiom)
    soil_mesh_local = Mesh(soil_local)

    ss_graph = RA(-90.0) + T(Vec(x_min, y_center, 0.0001)) +
           Soil_surface(length = sl, width = sw)
    ss_local = Mesh(Graph(axiom = ss_graph))

    # ── One diag per plant ── (moved up, before it's used)
    diags = [DiagnosticTimeSeries(label = "$label — plant $p")
             for p in 1:length(forest_local)]

    # ── Feature 2 accumulators ──
    plant_ids = [p[4] for p in positions]
    pp_all = DataFrame[]
    pv_all = DataFrame[]

    # -- extract positions--

     # ── Export plant positions (for spatial plotting in R) ──
    pos_df = DataFrame(
        x         = [Float64(p[1][1]) for p in positions],
        y         = [Float64(p[1][2]) for p in positions],
        species   = [p[2] for p in positions],
        plant_id  = [p[4] for p in positions]
    )
    # ── Single daily loop (no duplication) ──
    for i in 1:n_days
        DOY = start_DOY + i - 1
        daily_step!(forest_local, soil_local, ss_local, DOY, x_min, x_max, y_min, y_max, sl, sw)
        for (p, tree) in enumerate(forest_local)
            record!(diags[p], i, tree, soil_local)
        end

        pp_day, pv_day = root_voxel_occupancy(forest_local, soil_local, plant_ids, DOY)
        push!(pp_all, pp_day)
        push!(pv_all, pv_day)
    end

    # ── Write Feature 2 CSVs (before return, not after!) ──
    pp_full = vcat(pp_all...)
    pv_full = vcat(pv_all...)
    CSV.write("$(label)_root_per_plant.csv", pp_full)
    CSV.write("$(label)_root_per_voxel.csv", pv_full)

    # plant positions
    CSV.write("$(label)_plant_positions.csv", pos_df)


    # ── Write run metadata (Feature 1) ──
    write_run_metadata(
        "run_metadata_$(label).csv";
        label = label, sp1 = sp1, sp2 = sp2,
        n_plants_row_sp1 = n_plants_row_sp1, plant_spacing_sp1 = plant_spacing_sp1,
        n_rows_sp1 = n_rows_sp1, row_spacing_sp1 = row_spacing_sp1,
        n_plants_row_sp2 = n_plants_row_sp2, plant_spacing_sp2 = plant_spacing_sp2,
        n_rows_sp2 = n_rows_sp2, row_spacing_sp2 = row_spacing_sp2,
        soil_N_init = soil_N_init, sv = sv,
        n_i = n_i, n_j = n_j, n_k = Int(n_k),
        n_plants_sp1 = n_plants_row_sp1 * n_rows_sp1,
        n_plants_sp2 = n_plants_row_sp2 * n_rows_sp2,
        EN_coeff = data(forest_local[1]).EN_coeff,
    )

    return diags, forest_local, soil_local, ss_local, soil_mesh_local
end

#end

# ═══════════════════════════════════════════════════════════════════
# SECTION: MULTI-PLANT EXPORT + MEAN-PLANT HELPERS (NEW)
# ═══════════════════════════════════════════════════════════════════

"""
    diag_to_dataframe(diag; N_level=missing, plant=missing)

Convert one plant's DiagnosticTimeSeries into a tidy DataFrame (one row per day).
Adds scenario tag columns (N_level, plant) so all runs can be stacked together.
"""
function diag_to_dataframe(diag::DiagnosticTimeSeries; N_level = missing, plant = missing)
    df = DataFrame()
    df.N_level = fill(N_level, length(diag.day))
    df.plant   = fill(plant,   length(diag.day))
    for f in fieldnames(DiagnosticTimeSeries)
        v = getfield(diag, f)
        if v isa AbstractVector && length(v) == length(diag.day)
            df[!, f] = v
        end
    end
    return df
end

"""
    mean_diag(diags)

Return a new DiagnosticTimeSeries whose vector fields are the per-day mean
across all plants in `diags`. Integer count fields are rounded back to Int.
Scalar fields (seed_reserve_days, label, …) are copied from plant 1.
"""
function mean_diag(diags::Vector{<:DiagnosticTimeSeries})
    @assert !isempty(diags) "no plants to average"
    out = DiagnosticTimeSeries()
    n   = length(diags)
    ref = diags[1]
    for f in fieldnames(DiagnosticTimeSeries)
        ref_val = getfield(ref, f)
        if ref_val isa AbstractVector
            L = length(ref_val)
            if eltype(ref_val) <: Integer
                acc = zeros(Float64, L)
                for d in diags; acc .+= getfield(d, f); end
                setfield!(out, f, round.(Int, acc ./ n))
            elseif eltype(ref_val) <: AbstractFloat
                acc = zeros(Float64, L)
                for d in diags; acc .+= getfield(d, f); end
                setfield!(out, f, acc ./ n)
            else
                setfield!(out, f, copy(ref_val))
            end
        else
            setfield!(out, f, ref_val)
        end
    end
    out.day = copy(ref.day)
    return out
end

# ═══════════════════════════════════════════════════════════════════
# SECTION: RUN EVERYTHING
# ═══════════════════════════════════════════════════════════════════

# ═══════════════════════════════════════════════════════════════════
# SECTION: RUN EVERYTHING
# ═══════════════════════════════════════════════════════════════════

### with new params for field config. 
# ═══════════════════════════════════════════════════════════════════
#  5 nitrogen scenarios · 2 rows × 6 plants · run to 55 days
#    - store ALL plants to one CSV
#    - plot the MEAN plant (representative) for each N level
#    - render + save each 3D scene
# ═══════════════════════════════════════════════════════════════════

# ═══════════════════════════════════════════════════════════════════
#  VOXEL SIZE SENSITIVITY
#  5 voxel sizes · 2 rows × N plants · run to 35 days
#    - store ALL plants to one CSV
#    - plot the MEAN plant for each sv level
#    - render + save each 3D scene
# ═══════════════════════════════════════════════════════════════════

sv_levels = [0.02]  # ← voxel sizes (m)

diags_all  = Vector{Any}(undef, length(sv_levels))
forest_all = Vector{Any}(undef, length(sv_levels))
ss_all     = Vector{Any}(undef, length(sv_levels))
sm_all     = Vector{Any}(undef, length(sv_levels))

master_df = DataFrame()   # every plant of every run goes here

for (k, sv) in enumerate(sv_levels)
    println("\n", "="^70)
    println("  RUN $k / $(length(sv_levels))  —  sv = $(sv) m  ($(sv*100) cm)")
    println("="^70)

    diag_full, forest_full, soil_full, ss_full, sm_full = run_full_model(
        n_days      = 35,
        start_DOY   = 180,
        soil_N_init = 7000.0,
        sv          = sv,           # ← the key change
        label       = "sv=$(sv)",
        # ── same field config as before ──
        n_rows_sp1        = 0,
        sp2               = :cereal2,
        n_plants_row_sp2  = 12,
        n_rows_sp2        = 2,
        plant_spacing_sp2 = 0.02,
        row_spacing_sp2   = 0.12,
    )

    diags_all[k]  = diag_full
    forest_all[k] = forest_full
    ss_all[k]     = ss_full
    sm_all[k]     = sm_full

    nplants = length(diag_full)
    diag_rep = mean_diag(diag_full)   # representative mean plant

    # ── record all plants into master DataFrame ──
    for (p, diag) in enumerate(diag_full)
        df_p = diag_to_dataframe(diag)
        df_p[!, :sv]       .= sv
        df_p[!, :plant_id] .= p
        # ── soil N context columns ──
        voxel_volume_L          = sv^3 * 1000.0          # m³ → litres
        df_p[!, :soil_N_init_umol]       .= 7000.0      # total µmol per voxel
        df_p[!, :voxel_volume_L]         .= voxel_volume_L
        df_p[!, :soil_N_init_umol_per_L] .= 7000.0 / voxel_volume_L

        append!(master_df, df_p; cols = :union)
    end

    # ── plots ──
    fig_full = plot_diagnostics(diag_rep; title_prefix = "sv=$(sv) m (mean of $nplants)")
    savefig(fig_full, "diag_sv$(sv)_mean.png")
    display(fig_full)

    fig_budget_full = plot_carbon_budget(diag_rep; title_prefix = "sv=$(sv) m (mean of $nplants)")
    savefig(fig_budget_full, "carbon_budget_sv$(sv)_mean.png")
    display(fig_budget_full)
end

# ── write full per-plant dataset to one CSV ──
CSV.write("all_plants_all_sv.csv", master_df)
println("\nSaved: all_plants_all_sv.csv  ($(nrow(master_df)) rows)")



# ── render + save each 3D scene ──
for (k, sv) in enumerate(sv_levels)
    fig = render_forest(forest_all[k], ss_all[k], sm_all[k])
    GLMakie.save("scene_sv$(sv).png", fig)
    println("  Saved: scene_sv$(sv).png")
end

k = 1
fig = render_forest(forest_all[k], ss_all[k], sm_all[k])


##### check umol convertion from L to per voxel volume
println(names(master_df))

select(
    combine(groupby(master_df, :sv),
        :voxel_volume_L        => first => :voxel_volume_L,
        :soil_N_init_umol      => first => :soil_N_init_umol,
        :soil_N_init_umol_per_L => first => :soil_N_init_umol_per_L
    ), :sv, :voxel_volume_L, :soil_N_init_umol, :soil_N_init_umol_per_L
)

select(
    combine(groupby(master_df, :sv),
        :voxel_volume_L   => first => :voxel_volume_L,
    ),
    :sv,
    :voxel_volume_L,
    [:sv, :voxel_volume_L] => ((sv, vol) -> 30000.0 .* vol) => :N_umol_per_voxel,
    [:sv, :voxel_volume_L] => ((sv, vol) -> fill(30000.0, length(sv))) => :concentration_umol_per_L
)

for sv in [0.005, 0.01, 0.02, 0.03, 0.04, 0.06, 0.08, 0.12]
    println("sv=$sv → n_k=$(round(0.5/sv))  exact=$(0.5/sv == round(0.5/sv))")
end

a
###
# check first sv run
forest = forest_all[1]
tree = forest[1]
tvars = data(tree)

println("soil_x_min = $(tvars.soil_x_min) m = $(tvars.soil_x_min*100) cm")
println("soil_x_max = $(tvars.soil_x_max) m = $(tvars.soil_x_max*100) cm")
println("soil_y_min = $(tvars.soil_y_min) m = $(tvars.soil_y_min*100) cm")
println("soil_y_max = $(tvars.soil_y_max) m = $(tvars.soil_y_max*100) cm")
println("X domain = $((tvars.soil_x_max - tvars.soil_x_min)*100) cm")
println("Y domain = $((tvars.soil_y_max - tvars.soil_y_min)*100) cm")
###

## more tests for layout and ray tracer
n_plants_row_sp2  = 12
plant_spacing_sp2 = 0.02
n_rows_sp2        = 2
row_spacing_sp2   = 0.12

# actual plant column positions (same formula as create_tree positions)
xs = [(col - (n_plants_row_sp2 + 1)/2) * plant_spacing_sp2 for col in 1:n_plants_row_sp2]
ys = [(row - 1) * row_spacing_sp2 for row in 1:n_rows_sp2]

@show xs
@show ys
@show extrema(xs)
@show extrema(ys)

real_x_span = maximum(xs) - minimum(xs)
real_y_span = maximum(ys) - minimum(ys)
@show real_x_span real_y_span



## cloner uses 

sl_sp2 = n_plants_row_sp2 * plant_spacing_sp2
sw_sp2 = n_rows_sp2 * row_spacing_sp2

sl = max(0.1, sl_sp2)
sw = max(0.1, sw_sp2)

@show sl sw

settings = RTSettings(pkill=0.9, maxiter=4, nx=5, ny=5, dx=sl, dy=sw, parallel=true)
@show settings.dx settings.dy settings.nx settings.ny

test
# 
x_center = 0.0   # since positions are already centered around 0 in this axis
y_center = (maximum(ys) + minimum(ys)) / 2.0

x_min = x_center - sl/2
x_max = x_center + sl/2
y_min = y_center - sw/2
y_max = y_center + sw/2

@show x_min x_max y_min y_max
@show (x_max - x_min) (y_max - y_min)   # should equal sl, sw
@show minimum(xs) maximum(xs)           # real plants should sit strictly inside x_min..x_max
@show minimum(ys) maximum(ys)           # real plants should sit strictly inside y_min..y_max


#####
# ── render + save each 3D scene with fixed camera angle ──
eye    = GLMakie.Vec3f(0.01326,  1.5750, 0.0016)
lookat = GLMakie.Vec3f(-0.0022, 0.1750, -0.2758)
up     = GLMakie.Vec3f(0.0012, -0.1750, 0.9845)

for (k, sv) in enumerate(sv_levels)
    fig = render_forest(forest_all[k], ss_all[k], sm_all[k])

    # Apply fixed camera angle
    lscene = fig.content[1]
    cam = GLMakie.cam3d!(lscene.scene)
    cam.eyeposition[] = eye
    cam.lookat[]      = lookat
    cam.upvector[]    = up

    GLMakie.save("scene_sv$(sv).png", fig)
    println("  Saved: scene_sv$(sv).png")
end

aaa

# 3 which treatments produce empty voxels? 
# i noticed that something is wrong between what i expect to see and what i actually see 
# I expect 11 plants with 0.01m distance to be 0.11m or 11 cm 
# I expct voxel width or length to be 1 cm so 
# for 1 cm voxel 11/1= 11 voxels but i see 22
# for 2cm voxels 11/2= 5.5 voxels but i see 11
# for 3 cm 11/4= 3.66 but i see         7





aaa

N_levels = [30000.0]   # ← edit 6 N levels

diags_all  = Vector{Any}(undef, length(N_levels))
forest_all = Vector{Any}(undef, length(N_levels))
ss_all     = Vector{Any}(undef, length(N_levels))
sm_all     = Vector{Any}(undef, length(N_levels))

master_df = DataFrame()   # every plant of every run goes here



for (k, Nlev) in enumerate(N_levels)
    println("\n", "="^70)
    println("  RUN $k / $(length(N_levels))  —  soil_N_init = $Nlev µmol")
    println("="^70)

    diag_full, forest_full, soil_full, ss_full, sm_full = run_full_model(
        n_days      = 35, ##
        start_DOY   = 180,
        soil_N_init = Nlev,
        label       = "N=$(Nlev)",
        # ── 2 rows × 6 plants, single species (cereal2) ──
        n_rows_sp1        = 0,
        sp2               = :cereal2,
        n_plants_row_sp2  = 20, # 
        n_rows_sp2        = 2,
        plant_spacing_sp2 = 0.01,
        row_spacing_sp2   = 0.12, #form 0.10
    )

    diags_all[k]  = diag_full
    forest_all[k] = forest_full
    ss_all[k]     = ss_full
    sm_all[k]     = sm_full

    nplants = length(diag_full)
    println("  This run has $nplants plants.")   

    # (a) store ALL plants in the master table
    for p in 1:nplants
        append!(master_df, diag_to_dataframe(diag_full[p]; N_level = Nlev, plant = p))
    end

    # (b) representative plant = MEAN of all plants
    diag_rep = mean_diag(diag_full)
    diag_rep.label = "N=$(Nlev) mean"
    # If plot_diagnostics errors on a scalar field, also copy it:
    # diag_rep.seed_reserve_days = diag_full[1].seed_reserve_days
    #
    # Prefer a MIDDLE plant instead of the mean? comment the lines above and use:
    # diag_rep = diag_full[4]

    # (c) the two combined graphs, for the representative plant
    fig_full = plot_diagnostics(diag_rep; title_prefix = "N=$(Nlev) (mean of $nplants)")
    savefig(fig_full, "diag_full_N$(Int(Nlev))_mean.png")
    display(fig_full)

    fig_budget_full = plot_carbon_budget(diag_rep; title_prefix = "N=$(Nlev) (mean of $nplants)")
    savefig(fig_budget_full, "carbon_budget_N$(Int(Nlev))_mean.png")
    display(fig_budget_full)
end

# ── write the full per-plant dataset to one CSV (opens directly in Excel) ──
CSV.write("all_plants_all_N.csv", master_df)
println("\nSaved: all_plants_all_N.csv  ($(nrow(master_df)) rows)")

# Native .xlsx instead? uncomment `using XLSX` (Change 1) and:
# XLSX.writetable("all_plants_all_N.xlsx", overwrite = true,
#     AllPlants = (collect(eachcol(master_df)), names(master_df)))
println(pwd())
filter(f -> occursin("plant_positions", f), readdir())
methods(run_full_model)
run_full_model(n_days=1, label="test_positions")
readdir() |> f -> filter(x -> occursin("test_positions", x), f)


# ── render + save each 3D scene (all 12 plants are in the forest) ──
for (k, Nlev) in enumerate(N_levels)
    fig = render_forest(forest_all[k], ss_all[k], sm_all[k])
    GLMakie.save("scene_N$(Int(Nlev)).png", fig)
    println("  Saved: scene_N$(Int(Nlev)).png")
end

k = 1   # 1=357, 2=714, 3=1071,4=1535.0, 4=2000, 5=2200

fig = render_forest(forest_all[k], ss_all[k], sm_all[k])
display(fig)


# ──────────────────────────────────────────────────────────────────
# Set the SAME camera angle for all three figures
# ──────────────────────────────────────────────────────────────────

eye    = GLMakie.Vec3f(0.01326,  1.5750, 0.0016)
lookat = GLMakie.Vec3f(-0.0022, 0.1750, -0.2758)
up     = GLMakie.Vec3f(0.0012, -0.1750, 0.9845)

for fig in [fig] #,fig_p2, fig_p3]
    lscene = fig.content[1]
    cam = GLMakie.cam3d!(lscene.scene)
    cam.eyeposition[] = eye
    cam.lookat[]      = lookat
    cam.upvector[]    = up
end

GLMakie.save("scene_horizontal.png", fig)

println("sl = $sl, sw = $sw, n_j = $n_j, n_i = $n_i")
