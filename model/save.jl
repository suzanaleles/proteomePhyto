#########################################################################
# Iterate over different environmental conditions and save model output #
#########################################################################

# August 2022; by Suzana Goncalves Leles

# Optimal solutions looping over different environmental conditions.
env = "Temperature (C°)" # should be one of the following: Light, Temperature (C°) or [DIN]

if env == "[DIN]"
    # Run varying [DIN].
    change_DIN = 10.0 .^ [-7.0:0.18:log10(20.0);]
    change_T = 20.0 + 273.0
    change_I = 20.0
elseif env == "Temperature (C°)"
    # Run varying Temperature.
    change_DIN = 20.0 # 0.05 # for N-limited
    change_I = 20.0
    change_T  = 273.0 .+ [3.0:1.0:45.0;]
elseif env == "Light"
    # Run varying Light.
    change_DIN = 20.0
    change_T = 20.0 + 273.0
    change_I = 10.0 .^ [-7.0:0.18:log10(20.0);]
end

# Define empty vectors
s_μ     = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_β     = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_ctr   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_cri   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_clb   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_cup   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_cp    = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_cru   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_cld   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_cgl   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_cre   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_ϕtr   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_ϕri   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_ϕlb   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_ϕp    = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_ϕru   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_ϕld   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_ϕgl   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_ϕre   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_cin   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_clm   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_cic   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_Qcell = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_vru   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_vres  = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_vp    = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_vgl   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_vu    = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_vf    = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_vlb   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_vld   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_vri   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_vtr   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_vol   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_totp  = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_αlm   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_D     = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_γTa   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_γTd   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_γTad  = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_γTm   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_cli   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_cgu   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_ecost = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_ηaan  = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_ϵ     = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_αday  = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_αnight = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))

# Store model output in empty vectors
for i in 1:length(change_DIN)
        set_value(DIN, change_DIN[i])
        for j in 1:length(change_T)
                set_value(T, change_T[j])
                for k in 1:length(change_I)
                        set_value(I, change_I[k])
                        optimize!(model)
                        s_μ[j,i,k]   = exp.(JuMP.value.(logμ))
                        s_β[j,i,k]   = JuMP.value.(β)
                        s_ctr[j,i,k] = JuMP.value.(ctr)
                        s_cri[j,i,k] = JuMP.value.(cri)
                        s_clb[j,i,k] = JuMP.value.(clb)
                        s_cup[j,i,k] = JuMP.value.(cup)
                        s_cp[j,i,k]  = JuMP.value.(cp)
                        s_cru[j,i,k] = JuMP.value.(cru)
                        s_cld[j,i,k] = JuMP.value.(cld)
                        s_cgl[j,i,k] = JuMP.value.(cgl)
                        s_cre[j,i,k] = JuMP.value.(cre)
                        s_ϕtr[j,i,k] = JuMP.value.(ϕtr)
                        s_ϕri[j,i,k] = JuMP.value.(ϕri)
                        s_ϕlb[j,i,k] = JuMP.value.(ϕlb)
                        s_ϕp[j,i,k]  = JuMP.value.(ϕp)
                        s_ϕru[j,i,k] = JuMP.value.(ϕru)
                        s_ϕld[j,i,k] = JuMP.value.(ϕld)
                        s_ϕgl[j,i,k] = JuMP.value.(ϕgl)
                        s_ϕre[j,i,k] = JuMP.value.(ϕre)
                        s_cin[j,i,k] = JuMP.value.(cin)
                        s_clm[j,i,k] = JuMP.value.(clm)
                        s_cic[j,i,k]  = JuMP.value.(cic)
                        s_Qcell[j,i,k] = value(Qcell)
                        s_vru[j,i,k]  = value(vru)
                        s_vres[j,i,k] = value(vres)
                        s_vp[j,i,k]   = value(vp)
                        s_vgl[j,i,k]  = value(vgl)
                        s_vu[j,i,k]   = value(vu)
                        s_vf[j,i,k]   = value(vf)
                        s_vlb[j,i,k]  = value(vlb)
                        s_vld[j,i,k]  = value(vld)
                        s_vri[j,i,k]  = value(vri)
                        s_vtr[j,i,k]  = value(vtr)
                        s_vol[j,i,k]  = value(vol)
                        s_totp[j,i,k] = value(tot_p)
                        s_αlm[j,i,k]  = JuMP.value.(αlm)
                        s_D[j,i,k]    = value(D)
                        s_γTa[j,i,k]  = value(γTa)
                        s_γTd[j,i,k]  = value(γTd)
                        s_γTad[j,i,k] = value(γTad)
                        s_γTm[j,i,k]  = value(γTm)
                        s_cli[j,i,k]  = value(cli)
                        s_cgu[j,i,k]  = value(cgu)
                        s_ecost[j,i,k] = value(ecost)
                        s_ηaan[j,i,k]  = value(ηaan)
                        s_ϵ[j,i,k]     = value(ϵ)
                        s_αday[j,i,k]  = JuMP.value.(αday)
                        s_αnight[j,i,k]  = JuMP.value.(αnight)
                end
        end
end

using CSV
using DataFrames

output = DataFrame(temp = variable, mu = s_μ[a,b,c]*60*24, size = s_β[a,b,c]*3*2,
                   NC = s_Qcell[a,b,c], phot = s_vru[a,b,c], resp = s_vres[a,b,c],
                   totp = s_totp[a,b,c],
                   chlvol = s_cp[a,b,c] + s_cup[a,b,c],
                   ilpd = s_ϕld[a,b,c], icbd = s_ϕgl[a,b,c], ichl = s_ϕp[a,b,c], icha = s_ϕre[a,b,c],
                   itr = s_ϕtr[a,b,c], ilpb = s_ϕlb[a,b,c],
                   irib = s_ϕri[a,b,c], irub = s_ϕru[a,b,c], clip = s_clm[a,b,c], ctr = s_ctr[a,b,c],
                   cp = s_cp[a,b,c], cup = s_cup[a,b,c], cri = s_cri[a,b,c], cli = s_cli[a,b,c], cgu = s_cgu[a,b,c])
CSV.write("output.csv", output, header = true)
