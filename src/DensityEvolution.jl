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


function box_plus(a::AbstractFloat, b::AbstractFloat)
    dividend = 1 + exp(a + b)
    divisor = exp(a) + exp(b)

    return log(dividend / divisor)
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

# Sample vector LLRs
LLRs = Vector{Float64}(undef, n)

############################################################################################
# Sample Normal distribution
############################################################################################

# Calculate step size
step = (last(x)-x[1])/n

# Calculate sample Values
function sample()
    x_n = x[1]
    for i in 1:n
        LLRs[i] = f(x_n, σ, μ)
        x_n += step
    end
end
sample()

############################################################################################
# Calculate check node LLRs by boxplus operation
############################################################################################
# init
LLR_cn[i] = LLRs[i]
for i in 2:n
    LLR_ext_cn[i] = boxplus(LLR_cn[i], LLRs[i])
end 


############################################################################################
# Show initial results
############################################################################################
# LLR value for normal distribution
@show LLR(σ, μ)

# Plot function (divided by 6 for normalization)
y(x) = f(x, σ, μ) / d_c

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
