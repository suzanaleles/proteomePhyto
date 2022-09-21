#################################################################
# Equations defining the proteome model for a phytoplanktn cell #
#################################################################

# August 2022; by Suzana Goncalves Leles

using JuMP, Ipopt


# Create a JuMP model, using Ipopt as the solver.
model = Model(Ipopt.Optimizer)

# Set the environment: parameters that we can iterate over n times
@NLparameter(model, T   == 20.0 + 273.0)   # K,  temperature
@NLparameter(model, DIN == 20.0)           # uM, external dissolved inorganic nitrogen
@NLparameter(model, DIC == 100.0)          # uM, dissolved inorganic carbon
@NLparameter(model, I   == 120.0)          # μmol photon m-2 s-1, light

#Define the unknown variables that we will optimize for.
@variables(model, begin
        logμ    , (start = 0.000308)  # per day, growth rate
        β   ≥ 0, (start = 2.34)       # um, volume to surface area ratio
        ctr ≥ 0, (start = 411)        # protein: number of  transporters / um3
        cri ≥ 0, (start = 1.1e3)      # protein: number of ribosomes / um3
        clb ≥ 0, (start = 65.3)       # protein: number of lipid biosynthesis pathway proteins / um3
        cp  ≥ 0, (start = 1e3)        # protein: number of photosystems / um3
        cup ≥ 0, (start = 55.7)       # protein: number of damaged photosystems / um3
        cru ≥ 0, (start = 6e3)        # protein: number of rubisco / um3
        cld ≥ 0, (start = 0.0001)     # protein: number of lipid degradation pathway proteins / um3
        cgl ≥ 0, (start = 2.5e3)      # protein: number of glycolysis pathway prteins / um3
        cre ≥ 0, (start = 82.7)       # protein: number of repair proteins / um3
        ϕtr ≥ 0, (start = 0.002)      # relative proteome investment on: transporters
        ϕri ≥ 0, (start = 0.098)      # relative proteome investment on: ribosomes
        ϕlb ≥ 0, (start = 0.004)      # relative proteome investment on: lipid biosynthesis
        ϕp  ≥ 0, (start = 0.04)       # relative proteome investment on: photosystems
        ϕru ≥ 0, (start = 0.16)       # relative proteome investment on: rubisco
        ϕld ≥ 0, (start = 0.04)       # relative proteome investment on: lipid degradation
        ϕgl ≥ 0, (start = 0.18)       # relative proteome investment on: glycolysis
        ϕre ≥ 0, (start = 0.002)      # relative proteome investment on: repair
        cin ≥ 0, (start = 3e6)        # internal pool of: nitrogen
        clm ≥ 0, (start = 8e5)        # internal pool of: membrane lipids
        cic ≥ 0, (start = 3e6)        # internal pool of: carbon
        1 ≥ αlm ≥ 0,(start = 0.006)   # fraction of lipid biosynthesis flux used to fuel dark respiration
        1 ≥ αday ≥ 0,   (start = 0.0) # fraction of cgl proteins active during the day
        1 ≥ αnight ≥ 0, (start = 1.0) # fraction of cgl proteins active during the night
end)

# Set the objective: maximize the log growth rate (logμ).
@objective(model, Max, logμ)

# Define any intermediate calculations that will be called in the system of equations.
@NLexpressions(model, begin
        γTa,   exp(Ea*(1.0/(R*Tref) - 1.0/(R*T)))                                              # unitless
        γTd,   exp(Ed*(1.0/(R*Tu) - 1.0/(R*T)))                                                # unitless
        γTad,  exp(Etr*(1.0/(R*Tref) - 1.0/(R*T)))                                             # unitless
        γTm,   (Tmax - T)/(Tmax - Tmin)                                                        # unitless
        vp,  ( kref["p"] ) * cp * I / (K["p"] + I)                                             # moles of photons/μm3/min
        vru, ( kref["ru"] * γTa ) * cru * DIC / ( K["ru"] + DIC )                              # moles of carbon/μm3/min
        vtr, ( kref["tr"] * γTad ) * ctr * DIN/ ( K["tr"] + DIN)                               # moles of nitrogen/μm3/min
        vri, ( kref["ri"] * γTa ) * cri * (cic / ( K["ri"] + cic )) * (cin / ( K["ri"] + cin)) # moles of amino acids/μm3/min
        vlb, ( kref["lb"] * γTa ) * clb * cic / ( K["lb"] + cic )                              # moles of carbon/μm3/min
        vld, ( kref["ld"] * γTa ) * cld                                                        # moles of carbon/μm3/min
        vgl_day, ( kref["gl"] * γTa ) * (cgl * αday) * cic / ( K["gl"] + cic )                 # moles of carbon/μm3/min
        vgl_night, ( kref["gl"] * γTa ) * (cgl * αnight)                                       # moles of carbon/μm3/min
        vu,  ( kref["u"]  * γTd ) * cp                                                         # moles of photosystems/μm3/min
        vf,  ( kref["f"]  * γTa ) * cre * (cup/(cup + cp))                                     # moles of photosystems/μm3/min
        #vf, ( kref["f"] * γTa ) * cre * (cup/(cup + K["f"]))                                  # moles of photosystems/μm3/min; alternative approach
        vgl,   vgl_night + vgl_day                                                             # moles of carbon/μm3/min
        vres,  vld + vgl_night + vgl_day                                                       # moles of carbon/μm3/min
        cgu,  vgl_night / ηguc * dt                                                            # moles of glucose/μm3
        cli,  vld / ηlic * dt                                                                  # moles of lipid droplet (6C TAG)/μm3
        ηaan,  ηaac * ((ϕtr + ϕru + ϕlb + ϕld + ϕgl + ϕre)*Qpt + ϕri*Qri + ϕp*Qp)              # moles of nitrogen / mol of amino acid
        Qnc,    cin/(cic + clm*ηlic + cli*ηlic + cgu*ηguc)                                     # nitrogen to carbon ratio
        Qtot,   ctr*ηtr + cru*ηru + clb*ηlb + cld*ηld + cgl*ηgl + cre*ηre + cri*ηri +
                (cp + cup)*ηp + cin*ηin + cic*ηic + clm*ηli + cli*ηli + cgu*ηgu + cot*ηot
        Qcell,  ( (ctr*ηtr + cru*ηru + clb*ηlb + cld*ηld + cgl*ηgl + cre*ηre)*Qpt +            # nitrogen to carbon quota of the cell
                (cri*ηri)*Qri + (cp + cup)*ηp*Qp + cot*ηot*Qother +
                (cin*ηin + cic*ηic + clm*ηli + cli*ηli + cgu*ηgu)*Qnc ) / (Qtot)
        vol,     4/3*π*(β*3)^3                                                                 # cell volume in um3
        ϵ,       1.0 - (cp + cup) * Vp - cli * Vli                                             # unitless
        tot_p,   (cp + cup)*ηp + cru*ηru + ctr*ηtr + cri*ηri + clb*ηlb + cld*ηld + cgl*ηgl + cre*ηre + cot*ηot # Da/μm3, total protein content
        ecost, vtr*e_tr + vri*e_ri + vlb*e_lb + vu*e_u + vf*e_re                               # moles of energy/μm3/min
        D, ηru*cru + ηri*cri + ηlb*clb + ηld*cld +
           ηgl*cgl + ηre*cre + ηin*cin + ηic*cic + cgu*ηgu + cot*ηot                           # Da/μm3, density of the cell
end)

# Define the system of equations and constraints.
@constraints(model, begin
        ϕp + ϕru + ϕtr + ϕri + ϕlb + ϕld + ϕgl + ϕre ==  0.5
end)

@NLconstraints(model, begin
        ctr / clm ≤ Mmin + (Mmax * γTm)                                          # constraint on the protein to lipid ratio in the membrane
        β * (s_tr * ctr + s_lm * clm) == 1.0                                     # cell size constraint
        cup / (cup + cp) == ( kref["u"] * γTd * cp ) / ( kref["f"] * γTa * cre ) # unfolding (vu) and folding (vf) dynamics assuming vu = vf
        ϕp * vri * 1/(ηp/ηaa) - exp(logμ) * (cp + cup) == 0.0                    # protein: moles of photosystems/μm3
        # alternative approach to model folding/unfolding of photosystems:
        #ϕp * vri * 1/(ηp/ηaa) + vf - vu - exp(logμ) * cp == 0.0                  # protein: folded photosystems/μm3
        #vu - vf - exp(logμ) * cup == 0.0                                         # protein: unfolded photosystems/μm3
        ϕru * vri * 1/(ηru/ηaa) - exp(logμ) * cru == 0.0                         # protein: rubisco/μm3
        ϕtr  * vri * 1/(ηtr/ηaa)  - exp(logμ) * ctr  == 0.0                      # protein: transporters/μm3
        ϕri * vri * 1/(ηri/ηaa) - exp(logμ) * cri == 0.0                         # protein: ribosomes/μm3
        ϕlb * vri * 1/(ηlb/ηaa) - exp(logμ) * clb == 0.0                         # protein: clb/μm3
        ϕld * vri * 1/(ηld/ηaa) - exp(logμ) * cld == 0.0                         # protein: cld/μm3
        ϕgl * vri * 1/(ηgl/ηaa) - exp(logμ) * cgl == 0.0                         # protein: cgl/μm3
        ϕre * vri * 1/(ηre/ηaa) - exp(logμ) * cre == 0.0                         # protein: cre/μm3
        vtr - vri*ηaan - exp(logμ) * cin == 0.0                                  # internal pool of nitrogen/μm3
        vru - vri*ηaac - vlb - vgl_day - vgl_night - exp(logμ) * cic == 0.0      # internal pool of carbon/μm3
        (1 - αlm) * vlb/ηlic - exp(logμ) * clm == 0.0                            # membrane lipids/μm3
        vld == αlm * vlb                                                         # lipid degradation rate; moles of carbon/μm3/min
        # energetic requirements:
        vp * e_p ≥ vru * e_ru                                                    # rubisco costs are paid by photosystems; energy/μm3/min
        ecost ≤ (vp*e_p - vru*e_ru) + vgl_day*e_gl                               # total ATP cost; energy/μm3/min
        ecost*fdr ≤ vgl_night*e_gl + vld*e_ld                                    # dark respiration; energy/μm3/min
        1/ϵ * D ≤ Dmax                                                           # density constraint; Da/μm3
end)

# Solve the optimization problem.
status = optimize!(model)
