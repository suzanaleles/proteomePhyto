#----------------------------------------------------------------------------#
# This file contains the parameter values to run the proteome model (eqn.jl) #
#----------------------------------------------------------------------------#

# December 2022; by Suzana Goncalves Leles

# Maximum turnover rates for all reactions (# molecules per protein per minute)
kref = Dict("tr" => 120.0, "ri" => 114.0, "lb" => 120.0, "p" => 7800.0,
            "ru" => 157.0, "ld" => 120.0, "gl" => 120.0, "re" => 1000.0, "d" => 1000.0)
#map!(x -> x + x*0.5, values(kref))

# Half-saturation constants for all reactions
K = Dict("tr" => 1.0, "ri" => 10^4, "lb" => 10^4, "p" => 1,
          "ru" => 1.0, "ld" => 10^4, "gl" => 10^4, "re" => 10^4, "d" => 10^4)

# Space and density constraints
s_tr = 1.26e-5                   # um2, specific surface area of membrane transporter molecules
s_lm = 0.5e-6                    # um2, specific surface area of membrane lipid molecules
Mmax = 0.0005                    # unitless, maximum membrane transporter to lipid ratio
Mmin = 3e-5                      # unitless, minimum membrane transporter to lipid ratio
Vp = 3e-5                        # um3, volume / photosystem
Vli = 1.5e-9                     # um3, volume / lipid droplet
Dmax = 270/1e15*(6.02e23)/9      # Da/um3, maximum density of the cell
ctag_max = 2*10^7#4*10^7         # maximum number of lipid molecules / um3
ctag_min = 200                   # minimum number of lipid molecules / um3

# Temperature effects
Ttest = 273.0 .+ [3.0:1.0:45.0;] # K, temperature regime
Tmin = minimum(Ttest)            # K, minimum temperature
Tmax = maximum(Ttest)            # K, maximal temperature
Ea = 0.5                         # eV, activation energy  Maranon et al 2018; Barton et al 2020
Ed = 1.5                         # eV, deactivation energy
Etr  = Ea/2                      # eV, activation energy nutrient diffusion
R    = 8.62*10^-5                # eV/K, Boltzmann's constant
Tref = 20+273.0                  # K, reference temperature
Tu   = 43+273.0                  # K, critical high temperature

# Molecular weights
??aa = 110                        # Da / amino acid
??p  = 9082*??aa                   # Da / molecule of photosystem
??tr = 1046*??aa                   # Da / molecule of transpoter
??ru = 5933*??aa                   # Da / molecule of rubisco
??ri = 19393*??aa                  # Da / molecule of ribosome
??lb = 13742*??aa                  # Da / molecule of lipid synthesis pathway proteins
??gl = 15220*??aa                  # Da / molecule of glycolysis pathway proteins
??ld = 5128*??aa                   # Da / molecule of lipid degradation pathway proteins
??re = 6350*??aa                   # Da / molecule of repair proteins
??ot = 9486*??aa                   # Da / molecule of other proteome
??in = 14                         # 14  Da
??ic = 12                         # 12  Da
??gu = 180                        # 180 Da; 6C glucose
??li = 800                        # 800 Da; 16C TAG
cot = Dmax/4.9/110/(??ot/??aa)      # molecules of other proteome / um3

# Carbon and Nitrogen quotas
Qpt = 0.20                       # Nitrogen to Carbon ratio of proteins
Qp  = 0.10                       # Nitrogen to Carbon ratio of photosystems
Qri = 0.33                       # Nitrogen to Carbon ratio of ribosomes
Qother = 0.20                    # Nitrogen to Carbon ratio of other proteins not simulated here
??lic = 16                        # molecules of C/ molecule of lipid
??guc = 6                         # molecules of C/ molecule of glucose
??aac = 5                         # molecules of C/ molecule of amino acid

# Energy conversion factors
e_p  = 1.0                       # energy / photon absorbed
e_ru = 10                        # energy / carbon
e_tr = 1.0                       # energy / nitrogen
e_ri = 3+9*5                     # energy / amino acid
e_lb = 3.9                       # energy / carbon
e_u  = 9e7/ctag_max              # energy / photosystem
e_re = 4.5                       # energy / photosystem
e_gl = 5.0                       # energy / carbon
e_ld = 6.6                       # energy / carbon

fdr = 0.25                       # fraction of respiration that happens in the dark
td = 12*60                       # minutes, dark period (12 hours)
