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
    p = zeros(1, particles_tot)

    for i in 1:particles_tot
        accept = 0
        while accept == 0
            x = rand(Uniform(0,T))
            y = rand(Uniform(0,theta))
            if x < nu
                if y < x*theta/nu
                    accept = 1
                    times[i] = x
                end
            elseif x > nu
                if y < (nu-x)*theta/(T-nu) + theta
                    accept = 1
                    times[i] = x
                end
            end
            p[i] = exp((times[i]-T)*xi)
        end
    end
    particles_surv = zeros(1,particles_tot)
    for j in 1:particles_tot
        particles_surv[j] = rand(Binomial(1,p[j]))
    end
    return sum(particles_surv)
end

function dust_n(n, a, b, T, xi)
    total_dust = zeros(1, n)
    for i in 1:n
        total_dust[i] = dust_ind(a, b, T, xi)
    end
    return sum(total_dust)
end
