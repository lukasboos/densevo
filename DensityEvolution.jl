using SpecialFunctions
using Plots

#-----------------variables--------------------------

#number of iterations
n = 1

d_v = 3

d_c = 6

x = -10:0.01:35 #-10.0:0.01:35.0

eb_n0 = 1.12

R(d_v, d_c) = 1 - d_v / d_c # Code Rate

r = R(d_v, d_c)

db_to_lin(x) = (10^(x/10)) # SNR from exponential to linear space

eb_n0_lin = db_to_lin(eb_n0) #1.2941958414499861

σ_ñ(eb_n0_lin, R) = 1/sqrt(2 * eb_n0_lin * R) # See Example 2.3.4.2

variance = σ_ñ(eb_n0_lin, r) #0.8790225168308843

μ = 2/(variance^2) #2.5883916828999722

σ = sqrt(4 / (variance^2)) #2.2752545716468617

# bits connected to considered check node (Annahme null codeword)
b = [0, 0, 0, 0, 0, 0]

#LLRs of b
LLR_b = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

#extrinsic LLRs for cn
LLR_ext_bn = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

#probability distribution of cn
p_L_μ = [0.0, 0.0, 0.0, 0.0]
p_L_σ = [0.0, 0.0, 0.0, 0.0]

#probability distribution from cn to vn after i-th iteration
p_L_c_v_μ = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
p_L_c_v_σ = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

#probability distribution of vn
p_L_dv_ext_μ = [0.0, 0.0, 0.0, 0.0]
p_L_dv_ext_σ = [0.0, 0.0, 0.0, 0.0]

#probability distribution from vn to cn after i-th iteration
p_L_v_c_μ = [0.0, 0.0, 0.0, 0.0]
p_L_v_c_σ = [0.0, 0.0, 0.0, 0.0]


#-----------------functions--------------------------

#p(L(0)vc)is identical to p(L(ch)vc) in initial decoding iteration
# also p_vc_init ist gleich der Normalverteilung
#p_vc_init(x,σ,μ) = 1 / (σ * sqrt(2 * pi)) * exp(-1/2 * ((x - μ) / σ)^2)
cdf(x, σ, μ) = 1/2 * (1 + erf((x - μ) / (σ * sqrt(2)))) # Siehe https://en.wikipedia.org/wiki/Normal_distribution
f(x, σ, μ) = 1 / (σ * sqrt(2 * pi)) * exp(-1/2 * ((x - μ) / σ)^2) # Gaussian Distribution (Normal Distribution)
LLR(σ, μ) = log((1-cdf(0, σ, μ))/cdf(0, σ, μ)) # Equation (2.50)

#Calculate Check node LLR by box-sum operation, see equation 2.51 & 2.52-------------------------------------

#LLRs_b sind in unserem Fall gleich LLR_ch
for j=1:length(b)
    LLR_b[j] = LLR(σ, μ)
    LLR_ext_bn[j] = LLR(σ, μ)
end
@show LLR_b

#Calculate LLR_ext_bn
for i=1:length(b)
    for j=1:length(b)
        if j==i
            continue
        else
            LLR_ext_bn[i]=log((1+exp(LLR_ext_bn[i]+LLR_b[j]))/(exp(LLR_ext_bn[i])+exp(LLR_b[j])))
        end
    end
end
@show LLR_ext_bn

#Calculate probability disribution of check node, see equation 2.59 (initial) & 2.60 & 2.61-------------------
#Fuer jeden checknode berechnen; Extrinsic beachten!

#initial see equation 2.59 (see equation on: https://ccrma.stanford.edu/~jos/sasp/Product_Two_Gaussian_PDFs.html)
p_L_μ[1]=(LLR_ext_bn[1]*σ^2+LLR_ext_bn[2]*σ^2)/(2*σ^2)
@show p_L_μ[1]
p_L_σ[1]=sqrt((σ^4)/2*σ^2)
@show p_L_σ[1]

#see equation 2.60 (see equation on: https://ccrma.stanford.edu/~jos/sasp/Product_Two_Gaussian_PDFs.html)
for i=2:4
    p_L_μ[i]=(p_L_μ[1]*σ^2+LLR_ext_bn[i+1]*σ^2)/(2*σ^2)
    p_L_σ[i]=sqrt((σ^4)/2*σ^2)
end

#see equation 2.60 (see equation on: https://ccrma.stanford.edu/~jos/sasp/Product_Two_Gaussian_PDFs.html)
p_L_c_v_μ[n]=(1/d_c)*(p_L_μ[1]+p_L_μ[2]+p_L_μ[3]+p_L_μ[4])
@show p_L_c_v_μ[1]
p_L_c_v_σ[n]=(1/d_c)*(p_L_σ[1]+p_L_σ[2]+p_L_σ[3]+p_L_σ[4])
@show p_L_c_v_σ[1]

#Calculate probability distribution of variable node, see equation 2.62--------------------------------
#Faltung bei Normalverteilung addiert sigma & my
p_L_dv_ext_μ[1]=p_L_μ[2]+p_L_μ[3]+p_L_μ[4]
p_L_dv_ext_σ[1]=p_L_σ[2]+p_L_σ[3]+p_L_σ[4]
p_L_dv_ext_μ[2]=p_L_μ[1]+p_L_μ[3]+p_L_μ[4]
p_L_dv_ext_σ[2]=p_L_σ[1]+p_L_σ[3]+p_L_σ[4]
p_L_dv_ext_μ[3]=p_L_μ[1]+p_L_μ[2]+p_L_μ[4]
p_L_dv_ext_σ[3]=p_L_σ[1]+p_L_σ[2]+p_L_σ[4]
p_L_dv_ext_μ[4]=p_L_μ[1]+p_L_μ[2]+p_L_μ[3]
p_L_dv_ext_σ[4]=p_L_σ[1]+p_L_σ[2]+p_L_σ[3]

#see equation 2.63
p_L_v_c_μ[n]=(1/d_v)*(p_L_dv_ext_μ[1]+p_L_dv_ext_μ[2]+p_L_dv_ext_μ[3])
p_L_v_c_σ[n]=(1/d_v)*(p_L_dv_ext_σ[1]+p_L_dv_ext_σ[2]+p_L_dv_ext_σ[3])
@show p_L_v_c_μ[n]
@show p_L_v_c_σ[n]