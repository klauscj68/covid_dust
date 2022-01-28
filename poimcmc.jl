using Random,CSV,DataFrames,SpecialFunctions,Distributions
VecVw = Union{Vector{Float64},
	      SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}
	      };

# data
"""
Return dictionary with default parameter values
 First niso entries of vectors correspond to isolation room data
 niso+1:niso+n are where number of infected people in building are given
"""
function data()
	prm = Dict{Symbol,Float64}();
	
	# cap number of infected people in building
	prm[:nmax] = 100.0; nmax = Int64(prm[:nmax]);

	# number of infected people in the building
	prm[:n] = 10.0;

	# individual infection times
	for i=1:nmax
		sym = Symbol("t"*string(i));
		prm[sym] = 0.0;
	end

	# λ-params
	# Γ-hyperparameters for shedding amplitude
	# Γ(α,β) ~ β^α*x^{α-1}exp(-βx)/Γ(α)
	#  mean: α/β
	#  var:  α/β^2
	prm[:Γα] = 0.015;
	prm[:Γβ] = 0.00875;

	#  shedding amplitude
	for i=1:nmax
		sym = Symbol("A"*string(i));
		prm[sym] = 1.7142857;
	end

	#  shedding amplitude position
	for i=1:nmax
		sym = Symbol("Aₓ"*string(i));
		prm[sym] = 3.5;
	end

	#  shedding duration
	for i=1:nmax
		sym = Symbol("L"*string(i));
		prm[sym] = 7.0;
	end

	# particle decay rate
	prm[:ξ] = log(2)/7;

	# Time after time 0 at which dust collected
	prm[:T] = 7.0;

	# dust measurement copies/mg dust
	prm[:Y] = 144.0;
	
	vkeys = [k for k in keys(prm)];
	return prm,vkeys
end

# mcmcrg
"""
Return dictionaries with bounding intervals for parameters and bools to 
say which parameters are varied
"""
function mcmcrg()
	prmrg = Dict{Symbol,Vector{Float64}}();
	prmvary = Dict{Symbol,Bool}();
	
	# max permitted number of infected people in building
	#  nmax should agree with what is in data and not be varied
	prmrg[:nmax] = [100.0,101.0]; nmax = Int64(prmrg[:nmax][2]);
	prmvary[:nmax] = false;

	# number of infected people in the building
	#  bound by nmax enforced in prior
	prmrg[:n] = [1.0,100.0]; # For rej stats have be 1 less than nmax
	prmvary[:n] = true;

	# individual infection times
	for i=1:nmax
		sym = Symbol("t"*string(i));
		prmrg[sym] = [-3.0,7.0]; # max t₀ should be less than the fixed T-value
					 # for efficient sampling
		prmvary[sym] = true;
	end

	# λ-params # maybe 50% pickup in dorms by vacuum
	# Γ-distribution hyperparameters for amplitude
	prmrg[:Γα] = [0.001725,0.1725];
	prmvary[:Γα] = false;

	prmrg[:Γβ] = [0.0002225,0.02225];
	prmvary[:Γβ] = false;

	#  shedding amplitude
	for i=1:nmax
		sym = Symbol("A"*string(i));
		prmrg[sym] = [0.0,Inf];
		prmvary[sym] = true;
	end

	#  shedding amplitude position
	#   presently prp! uses a Unif(Tri) for LxAx and could
	#   have problem if L is varied while Ax is not. Similarly
	#   min permitted value of Ax should be <= min permitted
	#   value of L
	for i=1:nmax
		sym = Symbol("Aₓ"*string(i));
		prmrg[sym] = [0.0,14.0];
		prmvary[sym] = true;
	end

	#  shedding duration 
	for i=1:nmax
		sym = Symbol("L"*string(i));
		prmrg[sym] = [7.0,14.0];
		prmvary[sym] = true;
	end

	# particle decay rate
	prmrg[:ξ] = log(2)./[7.0,14.0];
	prmvary[:ξ] = false;

	# Time after time 0 at which dust is collected
	prmrg[:T] = [5.0,10.0];
	prmvary[:T] = false; # Current code needs false because prp doesn't vary

	# dust measurement copies/mg dust
	prmrg[:Y] = [0.0,1000.0]; # with Delta numbers over 1000 can go up to 10,000
	prmvary[:Y] = false; # Current code needs false becase prp doesn't vary

	return prmrg,prmvary
end

# wrtprm
"""
Write the prm dictionary to a column vector for storing to csv's
Uses multiple dispatch
call with no args: returns the dimension of column vector needed to store 
                   and list of keys
call with dictionary etc: returns the column vector stored in order of keys(prm)
"""
function wrtprm()
	prm,vkeys = data();
	
	# Create a vector of aprp size
	nelm = length(vkeys);

	V = Vector{Float64}(undef,nelm);
	for i=1:nelm
		V[i] = prm[vkeys[i]];
	end
		       
	return prm,vkeys,V
end
function wrtprm!(prm::Dict{Symbol,Float64},vkeys::Vector{Symbol},
                   V::VecVw)

	for i=1:length(vkeys)
		V[i] = prm[vkeys[i]];
	end
	
end
function wrtprm!(prm1::Dict{Symbol,Float64},
		 prm2::Dict{Symbol,Float64})
	for key in keys(prm1)
		prm2[key] = prm1[key];
	end
end

# rdprm
"""
Read a column vector formatted like wrtprm into a dictionary for 
restarting runs assuming each parameter has size 1.
"""
function rdprm(V::Vector{Float64},vkeys::Vector{Symbol})
	prm=Dict{Symbol,Float64}();
	for i=1:length(vkeys)
		prm[vkeys[i]] = V[i];
	end

	return prm,vkeys
end

# shedλ
"""
Compute the shedding μp = ∫ᵀ₀exp[-ξ*(T-t)]λ(t-t₀;θ)dt as function of input 
parameters. Multiple dispatch for case of single param values and a dictionary.
Also include a mutating version of the dictionary case for mem alloc. Right now
routine can handle Aₓ=0 or L but not ξ=0.
"""
function shedλ(A::Float64,L::Float64,t₀::Float64,T::Float64,Aₓ::Float64,
	       ξ::Float64)
	a1 = 0.0 >= t₀ ? 0.0 : t₀; b1 = T <= Aₓ+t₀ ? T : Aₓ+t₀;
	a2 = 0.0 >= Aₓ+t₀ ? 0.0 : Aₓ+t₀; b2 = T <= L+t₀ ? T : L+t₀;
	lnA = log(A);

	# ∫_{[0,T]∩[t₀,Aₓ+t₀]}exp(-ξ(T-t))exp(A*(t-t₀)/Aₓ)dt
	I₁ = ( (a1>=b1)||(Aₓ==0.0) ? 
	      0.0 : exp(-ξ*T-lnA/Aₓ*t₀)*1/(ξ+lnA/Aₓ)*( exp((ξ+lnA/Aₓ)*b1)-exp((ξ+lnA/Aₓ)*a1) ) 
	     );

	# ∫_{[0,T]∩[Aₓ+t₀,L+t₀]}exp(-ξ(T-t))exp(A-A*(t-t₀-Aₓ)/(L-Aₓ))dt
	I₂ = ( (a2>=b2)||(L==Aₓ) ?
	      0.0 : exp(-ξ*T+lnA+lnA/(L-Aₓ)*(t₀+Aₓ))*1/(ξ-lnA/(L-Aₓ))*( exp((ξ-lnA/(L-Aₓ))*b2)-exp((ξ-lnA/(L-Aₓ))*a2)  )  
	     );

	return I₁+I₂
end
function shedλ!(prm::Dict{Symbol,Float64};
		λval::Vector{Float64}=Vector{Float64}(undef,Int64(prm[:nmax])));
	for i=1:Int64(prm[:nmax])
		λval[i] = shedλ(prm[Symbol(:A,i)],prm[Symbol(:L,i)],
				prm[Symbol(:t,i)],prm[:T],
				prm[Symbol(:Aₓ,i)],
				prm[:ξ]);
	end

end
function shedλ(prm::Dict{Symbol,Float64})
	λval = Vector{Float64}(undef,Int64(prm[:nmax]));
	shedλ!(prm;λval=λval);

	return λval
end

# logΓ
"""
Compute log density of a Γ-distribution
"""
function logΓ(α::Float64,β::Float64,x::Float64)
	val = x>0 ? α*log(β)+(α-1)*log(x) - β*log(x) - log(gamma(α)) : -Inf;
	
	return val
end

# lognrm
"""
Compute log density of a normal distribution
"""
function lognrm(μ::Float64,σ::Float64,x::Float64)
	val = -0.5*log(2*π) - log(σ) - 0.5*( (x-μ)/σ )^2;

	return val
end

# logπ!
"""
Evaluate the log unnormalized posterior density for given choice of model parameters
prm::Dict storing the parameters we are evaluating
prmrg:: Dict storing the ranges parameters must belong to
prmvary:: Dict storing which parameters are varied
flagλval:: Bool saying if λval has already been evaluated for this param 
           choice or not (to save computation)
"""
function logπ!(prm::Dict{Symbol,Float64},
	       prmrg::Dict{Symbol,Vector{Float64}},prmvary::Dict{Symbol,Bool};
	       λval::Vector{Float64}=Vector{Float64}(undef,Int64(prm[:nmax])),
	       flagλval::Bool=true)	
	
	n = Int64(floor(prm[:n])); nmax = Int64(floor(prm[:nmax]));
	# Indicator prior on all params save shedding amplitudes A
	for key in keys(prmvary)
		if prmvary[key]&&( (prm[key]<prmrg[key][1])||(prm[key]>prmrg[key][2]) )
			return -Inf
		end
	end
	
	#  Indicator for Aₓ<=L
	if prmvary[Symbol(:Aₓ1)]
		for i=1:nmax
			Ai = Symbol(:Aₓ,i);
			Li = Symbol(:L,i);
			if prm[Ai]>prm[Li]
				return -Inf
			end
		end
	end
	
	#  Indicator for t₀<=T
	if prmvary[Symbol(:t1)]
		for i=1:nmax
			ti = Symbol(:t,i);
			if prm[ti]>prm[:T]
				return -Inf
			end
		end
	end
	
	# Priors on shedding amplitude conditioned on Γα,Γβ
	val1 = 0.0;
	for i=1:nmax
		Ai = Symbol(:A,i);
		val1 += logΓ(prm[:Γα],prm[:Γβ],prm[Ai]);
	end

	# Likelihood
	if flagλval
		shedλ!(prm;λval=λval);
	end
	val2 = 0.0;
	for i=1:n
		val2 += λval[i];
	end
	val2 = -val2 + floor(prm[:Y])*log(val2);	
	
	return val1+val2
end

# prp!
"""
Metropolis proposal function
"""
function prp!(prm0::Dict{Symbol,Float64},prm::Dict{Symbol,Float64},
		  prmrg::Dict{Symbol,Vector{Float64}},prmvary::Dict{Symbol,Bool};
		  λval::Vector{Float64}=Vector{Float64}(undef,Int64(prm[:nmax])),
		  rng::MersenneTwister=MersenneTwister())
	n = Int64(floor(prm[:n])); nmax = Int64(floor(prm[:nmax]));
	# Joint prp density of L x Ax is a unif trapezoid, so propto indicator
	if prmvary[:L1]
		for i=1:nmax
			Li = Symbol(:L,i);
			prm[Li] = prmrg[Li][1]+rand(rng)*(prmrg[Li][2]-prmrg[Li][1]);
		end
	end

	if prmvary[:Aₓ1]
		for i=1:nmax
			Axi = Symbol(:Aₓ,i); Li = Symbol(:L,i);
			prm[Axi] = prmrg[Axi][1] + rand(rng)*(prm[Li]-prmrg[Axi][1]);
		end
	end
	
	# prp density on n unif
	if prmvary[:n]
		prm[:n] = prmrg[:n][1]+rand(rng)*(prmrg[:n][2]-prmrg[:n][1])
	end

	# prp density on t₀'s is unif
	if prmvary[:t1]
		for i=1:nmax
			ti = Symbol(:t,i);
			prm[ti] = prmrg[ti][1]+rand(rng)*(prmrg[ti][2]-prmrg[ti][1]);
		end
	end
	
	# prp density on ξ is uniform
	if prmvary[:ξ]
		prm[:ξ] = prmrg[:ξ][1]+rand(rng)*(prmrg[:ξ][2]-prmrg[:ξ][1]);
	end

	# prp density on Gamma hypers is random walk
	if prmvary[:Γα]
		ΔΓα = 0.02*(prmrg[:Γα][2]-prmrg[:Γα][1]);
		prm[:Γα] = prm0[:Γα] + Δα*randn(rng);
	end

	if prmvary[:Γβ]
                ΔΓβ = 0.02*(prmrg[:Γβ][2]-prmrg[:Γβ][1]);
                prm[:Γβ] = prm0[:Γβ] + Δβ*randn(rng);
        end

	# Patched MH rejection
	#  Point is if alpha and beta aren't in cnst reg, the chain automatically rejects prm
	#  and resets to prm0. Used Julia wont sample a Gamma with nonpositive hyperparams
	if (prm[:Γα]<=0)||(prm[:Γβ]<=0)
		for key in keys(prm0)	
			prm[key] = prm0[key];
		end
		return
	end

	# prp density of amplitudes conditioned on hypers is Gamma
	if prmvary[:A1]
		Γdistr = Gamma(prm[:Γα],1/prm[:Γβ]);
		for i=1:nmax
			Ai = Symbol(:A,i);
			prm[Ai] = rand(rng,Γdistr);
		end
	end
end

# logρ!
"""
Evaluate the log unnormalized proposal density ρ(y|x) needed in mhratio
for the subset of parameters being varied. Only Gamma's tracked bc the unif's
in proposal density take same values all x ind of y and the random walk only
depends on |x-y| which is same y|x and x|y
"""
function logρ!(prm0::Dict{Symbol,Float64},prm::Dict{Symbol,Float64},
	       prmrg::Dict{Symbol,Vector{Float64}},prmvary::Dict{Symbol,Bool};
	       λval::Vector{Float64}=Vector{Float64}(undef,Int64(prm[:nmax][1])),
	       flagλval::Bool=true)
	val = 0.0;

	n = Int64(floor(prm[:n])); nmax = Int64(floor(prm[:nmax]));
	if prmvary[:A1]
		for i=1:nmax
			Ai = Symbol(:A,i);
			val += logΓ(prm[:Γα],prm[:Γβ],prm[Ai]);
		end
	end

	return val
end

# init!
"""
Initialize the mcmc sampler by a uniform draw conditioned on the posterior being
nonzero
"""
function init!(prm::Dict{Symbol,Float64},
	       prmrg::Dict{Symbol,Vector{Float64}},prmvary::Dict{Symbol,Bool}; 
	       λval::Vector{Float64}=Vector{Float64}(undef,Int64(prm[:nmax])),
	       rng::MersenneTwister=MersenneTwister())
	flagfd = false;
	while !flagfd
		for key in keys(prmvary)
			if prmvary[key]
				prm[key] = prmrg[key][1] .+ rand(rng)*(
						prmrg[key][2]-prmrg[key][1] < Inf ? prmrg[key][2]-prmrg[key][1] : 1000.0
					                                     );
			end
		end
		if prmvary[:Aₓ1]
			for i=1:Int64(prm[:nmax])
				Axi = Symbol(:Aₓ,i); Li = Symbol(:L,i);
				prm[Axi] = prmrg[Axi][1]+rand(rng)*(prm[Li]-prmrg[Axi][1]);
			end
		end
		
		if logπ!(prm,prmrg,prmvary;λval=λval) != -Inf
			flagfd = true;
		end
	end
end

# logmh!
"""
Compute the log Metropolis-Hastings acceptance ratio
"""
function logmh!(prm0::Dict{Symbol,Float64},prm::Dict{Symbol,Float64},
		prmrg::Dict{Symbol,Vector{Float64}},prmvary::Dict{Symbol,Bool};
	        λval::Vector{Float64}=Vector{Float64}(undef,Int64(prm[:nmax])))
	
	# stationary distribution
	val  = logπ!(prm,prmrg,prmvary;λval=λval);
	val -= logπ!(prm0,prmrg,prmvary;λval=λval);

	# proposal distribution
	val += logρ!(prm,prm0,prmrg,prmvary;λval=λval);
	val -= logρ!(prm0,prm,prmrg,prmvary;λval=λval);

	return val
end

# mcmcsmp!
"""
Run a single sample of the mcmc kernel. Routine allows for optionally cycling
the MH kernel a fixed number of times. Store eventual output to prm0 while 
writing proposals to prm. Returns the number of rejections.
"""
function mcmcsmp!(prm0::Dict{Symbol,Float64},prm::Dict{Symbol,Float64},
	          prmrg::Dict{Symbol,Vector{Float64}},prmvary::Dict{Symbol,Bool};
		  λval::Vector{Float64}=Vector{Float64}(undef,Int64(prm0[:nmax])),
		  rng::MersenneTwister=MersenneTwister(),
		  ncyc::Int64=1)
	nrej = 0;
	for i=1:ncyc
		prp!(prm0,prm,prmrg,prmvary;λval=λval,rng=rng)

		coin = rand(rng) |> log;
		if coin <=logmh!(prm0,prm,prmrg,prmvary;λval=λval)
			# accept
			for key in keys(prm0)
				if prmvary[key]
					prm0[key] = prm[key];
				end
			end
		else
			# reject
			nrej += 1;
		end
	end
	return nrej
end

# mcmcrun
"""
Run Metropolis-Hastings MCMC on dust measurement
"""
function mcmcrun(nsmp::Int64;
		 rng::MersenneTwister=MersenneTwister(),
		 flagrst::Bool=false,
		 ncyc::Int64=1,
		 Δprg::Float64=0.05)

	prmrg,prmvary = mcmcrg();
	# Initialize based on whether restarting from previous mcmc run
	if !flagrst
		prm,vkeys,V=wrtprm();
		λval = Vector{Float64}(undef,Int64(prm[:nmax]));
		init!(prm,prmrg,prmvary;λval=λval,rng=rng);

		prm0 = deepcopy(prm);	
	else
		df0 = CSV.read("MCMCsmp.csv",DataFrame);
		vkeys = Symbol.(df0[:,1]); V = df0[:,end];
		prm0,_ = rdprm(V,vkeys); prm = deepcopy(prm0);
		rng = myloadrng();
      	        λval = Vector{Float64}(undef,Int64(prm[:nmax]));
	end

	# Create matrix for samples
	SMP = Matrix{Float64}(undef,length(V),nsmp);

	# Create a vector for storing rejection statistics
	mhrej = Vector{Int64}(undef,nsmp);

	# run mcmc
	prg = 0.0;
	for i=1:nsmp
		
		mhrej[i] = mcmcsmp!(prm0,prm,prmrg,prmvary;
			            λval=λval,rng=rng,ncyc=ncyc);
		smp = @view SMP[:,i];
		wrtprm!(prm0,vkeys,smp);

		while i/nsmp>=prg+Δprg
			println("Progress through mcmc: $prg/1 ...");
			CSV.write("MCMCsmp.csv",[DataFrame(:prm=>String.(vkeys)) DataFrame(SMP[:,1:i])], writeheader=false,append=false);
			CSV.write("rejstats.csv",DataFrame(:rejct=>mhrej));
			mysaverng(rng);
			prg+=Δprg;
		end
	end

	# Save final csv and report rejection rates
	CSV.write("MCMCsmp.csv",[DataFrame(:prm=>String.(vkeys)) DataFrame(SMP[:,1:nsmp])], writeheader=false,append=false);
	CSV.write("rejstats.csv",DataFrame(:rejct=>mhrej),append=false);
	mysaverng(rng);

	rejrt = sum(mhrej)/(ncyc*nsmp); aptrt = 1-rejrt; aptwt = 1/aptrt;
	println("Rejection rate: $rejrt");
	println("Average num proposals before an accept: $aptwt");
end;
