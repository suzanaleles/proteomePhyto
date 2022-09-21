##############################################################################
# This file contains the parameter values to run the proteome model (eqn.jl) #
##############################################################################

# August 2022; by Suzana Goncalves Leles

# Maximum turnover rates for all reactions (# molecules per protein per minute)
kref = Dict("tr" => 120.0, "ri" => 114.0, "lb" => 120.0, "p" => 7800.0,
            "ru" => 157.0, "ld" => 120.0, "gl" => 120.0, "f" => 1000.0, "u" => 1000.0)

# Half-saturation constants for all reactions
K = Dict("tr" => 1.0, "ri" => 10^4, "lb" => 10^4, "p" => 1,
          "ru" => 1.0, "ld" => 10^4, "gl" => 10^4, "f" => 10^4, "u" => 10^4)

# Temperature effects
Ttest = 273.0 .+ [3.0:1.0:45.0;] # K, temperature regime
Tmin = minimum(Ttest)            # K, minimum temperature
Tmax = maximum(Ttest)            # K, maximal temperature
Ea = 0.5                         # eV, activation energy  Maranon et al 2018; Barton et al 2020
Ed = 1.5                         # eV, deactivation energy
Etr  = Ea/2                      # eV, activation energy nutrient diffusion
R    = 8.62*10^-5                # eV/K, Boltzmann's constant
Tref = 20+273.0                  # K, reference temperature
Tu   = 43+273.0                  # K, temperature at which half of enzymes become non-functional

# Space and density constraints
s_tr = 1.26e-5                   # um2, specific surface area of membrane transporter molecules
s_lm = 0.5e-6                    # um2, specific surface area of membrane lipid molecules
Mmax = 0.00055                   # unitless, maximum membrane transporter to lipid ratio
Mmin = 0.00005                   # unitless, maximum membrane transporter to lipid ratio
Vp = 2e-6                        # um3, volume/photosystem
Vli = 1.5e-9                     # um3, volume/lipid droplet
Dmax = 270/1e15*(6.02e23)/9      # Da/um3, maximum density of the cell

# Molecular weights
ηaa = 110                        # Da / amino acid
ηp  = 9082*ηaa                   # Da / molecule of photosystem
ηtr = 1046*ηaa                   # Da / molecule of transpoter
ηru = 5933*ηaa                   # Da / molecule of rubisco
ηri = 19393*ηaa                  # Da / molecule of ribosome
ηlb = 13742*ηaa                  # Da / molecule of lipid synthesis pathway proteins
ηgl = 15220*ηaa                  # Da / molecule of glycolysis pathway proteins
ηld = 5128*ηaa                   # Da / molecule of lipid degradation pathway proteins
ηre = 6350*ηaa                   # Da / molecule of repair proteins
ηot = 9486*ηaa                   # Da / molecule of other proteome
ηin = 14                         # 14  Da
ηic = 12                         # 12  Da
ηgu = 180                        # 180 Da; 6C glucose
ηli = 800                        # 800 Da; 16C TAG
cot = Dmax/4.9/110/(ηot/ηaa)     # molecules of other proteome / um3

# Carbon and Nitrogen quotas
Qpt = 0.2                        # Nitrogen to Carbon ratio of proteins
Qp  = 0.10                       # Nitrogen to Carbon ratio of photosystems
Qri = 0.25                       # Nitrogen to Carbon ratio of ribosomes
Qother = 0.2                     # Nitrogen to Carbon ratio of other proteins not simulated here
ηlic = 16                        # moles of C/mol of lipid
ηguc = 6                         # moles of C/mol of glucose
ηaac = 5                         # moles of C/mol of amino acid

# Energy conversion factors
e_p  = 1.0                       # energy / photon absorbed
e_ru = 9.4                       # energy / carbon
e_tr = 1.0                       # energy / nitrogen
e_ri = 3+9*5                     # energy / amino acid
e_lb = 4.5                       # energy / carbon
e_u  = 4.5                       # energy / photosystem
e_re = 4.5                       # energy / photosystem
e_gl = 5.0                       # energy / carbon
e_ld = 6.6                       # energy / carbon

fdr = 0.25                       # fraction of respiration that happens in the dark
dt = 12*60                       # light period (12h) in which the cell stores energy that will be burned at night (unit: minutes)
