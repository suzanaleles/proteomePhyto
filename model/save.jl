#-----------------------------------------------------------------------#
# Iterate over different environmental conditions and save model output #
#-----------------------------------------------------------------------#

# December 2022; by Suzana Goncalves Leles

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
s_ptr   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_pri   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_plb   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_pdp   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_pp    = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_pru   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_pld   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_pgl   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_pre   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
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
s_vd    = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_vre    = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_vlb   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_vld   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_vri   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_vtr   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_vol   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_totp  = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_αld   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_D     = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_γTa   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_γTd   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_γTad  = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_γTm   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_γTli  = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_cli   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_cgu   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_ecost = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_ηaan  = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_ϵ     = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_ctag = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_αlm  = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_αtag  = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))

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
                        s_ptr[j,i,k] = JuMP.value.(ptr)
                        s_pri[j,i,k] = JuMP.value.(pri)
                        s_plb[j,i,k] = JuMP.value.(plb)
                        s_pdp[j,i,k] = JuMP.value.(pdp)
                        s_pp[j,i,k]  = JuMP.value.(pp)
                        s_pru[j,i,k] = JuMP.value.(pru)
                        s_pld[j,i,k] = JuMP.value.(pld)
                        s_pgl[j,i,k] = JuMP.value.(pgl)
                        s_pre[j,i,k] = JuMP.value.(pre)
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
                        s_vd[j,i,k]   = value(vd)
                        s_vre[j,i,k]   = value(vre)
                        s_vlb[j,i,k]  = value(vlb)
                        s_vld[j,i,k]  = value(vld)
                        s_vri[j,i,k]  = value(vri)
                        s_vtr[j,i,k]  = value(vtr)
                        s_vol[j,i,k]  = value(vol)
                        s_totp[j,i,k] = value(tot_p)
                        s_αld[j,i,k]  = JuMP.value.(αld)
                        s_D[j,i,k]    = value(D)
                        s_γTa[j,i,k]  = value(γTa)
                        s_γTd[j,i,k]  = value(γTd)
                        s_γTad[j,i,k] = value(γTad)
                        s_γTm[j,i,k]  = value(γTm)
                        s_γTli[j,i,k] = value(γTli)
                        s_cli[j,i,k]  = value(cli)
                        s_cgu[j,i,k]  = value(cgu)
                        s_ecost[j,i,k] = value(ecost)
                        s_ηaan[j,i,k]  = value(ηaan)
                        s_ϵ[j,i,k]     = value(ϵ)
                        s_ctag[j,i,k]  = JuMP.value.(ctag)
                        s_αlm[j,i,k]  = JuMP.value.(αlm)
                        s_αtag[j,i,k]  = JuMP.value.(αtag)
                end
        end
end


using CSV
using DataFrames

output = DataFrame(temp = variable, mu = s_μ[a,b,c]*60*24, size = s_β[a,b,c]*3*2,
                   NC = s_Qcell[a,b,c], phot = s_vru[a,b,c], resp = s_vres[a,b,c],
                   totp = s_totp[a,b,c],
                   chlvol = s_pp[a,b,c] + s_pdp[a,b,c],
                   ilpd = s_ϕld[a,b,c], icbd = s_ϕgl[a,b,c], ichl = s_ϕp[a,b,c], icha = s_ϕre[a,b,c],
                   itr = s_ϕtr[a,b,c], ilpb = s_ϕlb[a,b,c],
                   irib = s_ϕri[a,b,c], irub = s_ϕru[a,b,c], clip = s_clm[a,b,c], ctr = s_ptr[a,b,c],
                   cp = s_pp[a,b,c], cup = s_pdp[a,b,c], cri = s_pri[a,b,c], 
                   cli = s_cli[a,b,c], cgu = s_cgu[a,b,c], ctag = s_ctag[a,b,c],
                   ctot = s_cli[a,b,c] * ηlic .+ s_cgu[a,b,c] * ηguc .+ s_ctag[a,b,c] * ηlic)
CSV.write("output_liang_gly_2.csv", output, header = true)
