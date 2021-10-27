using LinearAlgebra
using Distributions
using DelimitedFiles
using Distributed
using Random
Random.seed!(123);

cd("C:\\Users\\matth\\Documents\\JuliaWD")

function dust_ind(a, b, T, xi)
    theta = rand(Gamma(a,1/b))
    nu = rand(Uniform(0,T))
    lambda = (1/2)*T*theta
    particles_tot = rand(Poisson(lambda))
    times = zeros(1, particles_tot)
    p = zeroes(1, particles_tot)

    for i=1:particles_tot
        accept = 0
        while accept == 0
            x = rand(Uniform(0,T))
            y = rand(Uniform(0,theta))
            if x < nu
                if y < x*theta/nu
                    accept = 1
                    times[i] = x
            if x > nu
                if y < (nu-x)*theta/(T-nu) + theta
                    accept = 1
                    times[i] = x
            p[i] = exp(-(T-times[i])*xi)

    particles_surv = 0
    for j=1:particles_tot
        particle_surv = particles_surv + rand(Binomial(1,p[i]))

    return particles_surv
end
