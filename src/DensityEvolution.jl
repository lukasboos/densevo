using SpecialFunctions
using Plots

struct Gaussian
    μ::AbstractFloat
    σ::AbstractFloat
end

############################################################################################
# Functions
############################################################################################
# Code Rate
R(d_v, d_c) = 1 - d_v / d_c

# SNR from exponential to linear space
db_to_lin(x) = (10^(x / 10))

# See Example 2.3.4.2
σ_ñ(eb_n0_lin, R) = 1 / sqrt(2 * eb_n0_lin * R)

# Gaussian Distribution (Normal Distribution)
f(x, σ, μ) = 1 / (σ * sqrt(2 * pi)) * exp(-1 / 2 * ((x - μ) / σ)^2)

# Cumulative density function
# Siehe https://en.wikipedia.org/wiki/Normal_distribution
cdf(x, σ, μ) = 1 / 2 * (1 + erf((x - μ) / (σ * sqrt(2))))

# Equation (2.50)
LLR(σ, μ) = log((1 - cdf(0, σ, μ)) / cdf(0, σ, μ))

#looks like result is similar to a multiplication
function box_plus(a::Vector, b::Vector, n)
    result = ones(Float64, n)
    for i in 1:n
        dividend = 1 + exp(a[i] + b[i])
        divisor = exp(a[i]) + exp(b[i])
        result[i] = log(dividend / divisor)
    end
    @show length(result)
    @show length(result)
    return result
end


function Base.:*(
    a::Gaussian,
    b::Gaussian
)::Gaussian
    dividend = a.μ * b.σ^2 + b.μ * a.σ^2
    divisor = b.σ^2 + a.σ^2
    μ = dividend / divisor

    dividend = a.σ^2 * b.σ^2
    divisor = a.σ^2 + b.σ^2
    σ = sqrt(dividend / divisor)

    return Gaussian(μ, σ)
end

function convolution(
    a::Gaussian,
    b::Gaussian
)::Gaussian
    μ = a.μ + b.μ
    σ = sqrt(a.σ^2 + b.σ^2)

    return Gaussian(μ, σ)
end

function Base.:+(
    a::Gaussian,
    b::Gaussian
)
    μ = a.μ + b.μ
    σ = sqrt(a.σ^2 + b.σ^2)

    return Gaussian(μ, σ)
end

function Base.:/(
    a::Gaussian,
    b::Real
)
    return Gaussian(a.μ / b, a.σ / b)
end

function Base.:*(
    a::Gaussian,
    b::Real
)
    return Gaussian(a.μ * b, a.σ * b)
end

function Base.:*(
    a::Real,
    b::Gaussian
)
    return b * a
end


############################################################################################
# Configuration parameters
############################################################################################
d_v = 3
d_c = 6
eb_n0 = 1.12


# Code rate
r = R(d_v, d_c)

# Lineariation of SNR in dB range
eb_n0_lin = db_to_lin(eb_n0)

# variance
variance = σ_ñ(eb_n0_lin, r)

# mean value
μ = 2 / (variance^2)

# sigma
σ = sqrt(4 / (variance^2))


λ_d_v = zeros(Float64, d_v)
λ_d_v[d_v] = 1.0
λ_d_c = zeros(Float64, d_c)
λ_d_c[d_c] = 1.0

# Number of samples
n = 2000

# Plot range
x = -10:0.01:35

# Sample vector LLRs per check node since every check node contains same LLRs a vector is sufficient. Otherwise we need one per cn
samples = ones(Float64, n) #Vector{Float64}(undef, n)
# Extrinsic LLRs per check node since every check node contains same LLRs a vector is sufficient. Otherwise we need one per cn
LLR_ext_cn = ones(Float64, n, d_c)#Matrix{Float64}(undef, n, d_c)

############################################################################################
# Sample Normal distribution
############################################################################################

# Calculate sample Values for each check node
function sample(x, n, σ, μ, n, σ, μ)
    # Calculate step size
    steps = (last(x)-x[1])/n
    x_n = x[1]
    for i in 1:n
        samples[i] = f(x_n, σ, μ)
        x_n += steps
    end
    return samples
end
samples = sample(x, n, σ, μ)
@show length(samples)

# TODO: plot function of samples
# Plot function (divided by 6 for normalization)
y = samples ./ d_c
x=1:2000
#y=rand(10)

#plot(x, y)

plot(x, y, title = "Plotted Samples")
savefig("myplot.png")


############################################################################################
# Calculate check node LLRs by boxplus operation
############################################################################################
# init
#LLR_ext_cn[1] = samples[1]

# calculate LLR_ext_cn for d_c checknodes
function cn_init(d_c, a::Vector, n)
    result = Matrix{Float64}(undef, n, d_c)
    for edge_out in 1:d_c
        first = true
        for edge_in in 1:d_c
            if edge_out==edge_in
                continue
            else
                if first
                    result[:, edge_out] = a
                    first = false
                    #@show result[:, 1]
                else
                    result[:, edge_out] = box_plus(a, result[:, edge_out], n)# .* (d_c-1)
            end
        end
    end 
    return result
end 
LLR_ext_cn = cn_init(d_c, samples, n)

# Plot LLR_ext_cn (divided by 6 for normalization)
y = LLR_ext_cn[:, 1] ./ d_c
x=1:2000
#y=rand(10)

#plot(x, y)

plot(x, y, title = "Plotted LLR_ext_cn")
savefig("LLR_ext_plot.png")

############################################################################################
# Calculate probability distribution see formula 2.50 (log of left part divided by right part) see formula 2.50 (log of left part divided by right part)
############################################################################################

#=
function prob()
    LLR_cn_1 = 0.0
    LLR_cn_sum = 0.0
    p_cn = 0.0
    # LLR sum for 1
    for i in 1001:2000
        LLR_cn_1 = LLRs[i] + LLR_cn_1
    end
    # LLRs sum total
    for i in 1:2000
        LLR_cn_sum = LLR_ext_cn[i] + LLR_cn_sum
    end
    # calculate probability distribution
    p_cn = LLR_cn_1/LLR_cn_sum
    @show p_cn
end
LLR_cn = LLR_cn = prob()
=#=#
############################################################################################
# Show initial results
############################################################################################
# LLR value for normal distribution
@show LLR(σ, μ)

# Plot function (divided by 6 for normalization)
#y(x) = f(x, σ, μ) / d_c

plot(x, y)




############################################################################################
# First Iteration check node
############################################################################################
function de_cn(λ_d_c::Vector{Float64}, p_L_vn_cn::Gaussian)::Gaussian
    d_c = length(λ_d_c)

    p_L = Vector{Gaussian}(undef, d_c - 2)
    # p_L[0+1] = (p_L_vn_cn / (d_c - 1)) * (p_L_vn_cn / (d_c - 1))
    p_L[0+1] = p_L_vn_cn * p_L_vn_cn
    for i in 1:d_c-3
        # p_L[i+1] = p_L[i-1+1] * (p_L_vn_cn / (d_c - 1))
        p_L[i+1] = p_L[i-1+1] * p_L_vn_cn
    end

    p_L_cn_vn = λ_d_c[3] * p_L[1]
    for i in 4:d_c
        p_L_cn_vn += λ_d_c[i] * p_L[i-3+1]
    end

    return p_L_cn_vn
end

############################################################################################
# First Iteration variable node
############################################################################################
function de_vn(λ_d_v::Vector{Float64}, p_L_chv::Gaussian, p_L_cn_vn::Gaussian)::Gaussian
    d_v = length(λ_d_v)

    p_L_ext = Vector{Gaussian}(undef, d_v)
    # p_L_ext[1] = p_L_chv / d_v
    p_L_ext[1] = p_L_chv
    for i in 0:d_v-2
        # p_L_ext[i+2] = convolution(p_L_ext[i+1], p_L_cn_vn / d_v)
        p_L_ext[i+2] = convolution(p_L_ext[i+1], p_L_cn_vn)
    end

    p_L_vn_cn = λ_d_v[1] * p_L_ext[1]
    for i in 2:d_v
        p_L_vn_cn += λ_d_v[i] * p_L_ext[i]
    end

    return p_L_vn_cn
end


############################################################################################
# Density Evolution
############################################################################################
function de(λ_d_c::Vector{Float64}, λ_d_v::Vector{Float64}, μ::AbstractFloat, σ::AbstractFloat)
    # Before first checknode
    p_L_chv = Gaussian(μ, σ)
    @show p_L_chv
    vn_cn = Gaussian(μ, σ)
    @show vn_cn

    for i in 1:2
        cn_vn = de_cn(λ_d_c, vn_cn)
        @show i, cn_vn
        vn_cn = de_vn(λ_d_v, p_L_chv, cn_vn)
        @show i, vn_cn
    end

    return vn_cn
end

de(λ_d_c, λ_d_v, μ, σ)







################################################
c = ones(Float64, 8, 4)
c0 = [-4.49, -4.49, -4.49, -4.49]
c1 = [0.00, 0.00, 0.00, 0.00]
c2 = [-4.49, -4.49, -4.49, -4.49]
c3 = [-4.49, -4.49, -4.49, -4.49]
c4 = [4.49, 4.49, 4.49, 4.49]
c5 = [4.49, 4.49, 4.49, 4.49]
c6 = [-4.49, -4.49, -4.49, -4.49]
c7 = [-4.49, -4.49, -4.49, -4.49]

c[1, :] = c0
c[2, :] = c1
c[3, :] = c2
c[4, :] = c3
c[5, :] = c4
c[6, :] = c5
c[7, :] = c6
c[8, :] = c7

@show(c)

L = box_plus(c0, c2, 4)