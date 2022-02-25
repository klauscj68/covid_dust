using Random,CSV,DataFrames,SpecialFunctions,Distributions
VecVw = Union{Vector{Float64},
	      SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}
	      };

# data
"""
Return dictionary with default parameter values
"""
function data()
	prm = Dict{Symbol,Float64}();
	
	# cap number of infected people across all buildings
	#  Fall 2020 Iso: 39 avg ppl from Oct 21-27 for dust collected Oct 28th
	#  		  38 avg ppl from Oct 28th to Nov 3rd for dust collected Nov 4th
	prm[:nmax] = 80.0; nmax::Int64 = prm[:nmax];

	# number of buildings (several buildings used in calibration)
	prm[:nbld] = 2.0; nbld::Int64 = prm[:nbld];

	# number of infected people in each building
	prm[:n1] = 50.0;
	prm[:n2] = 30.0;

	# individual infection times
	@inbounds for i=1:nmax
		sym = Symbol("t"*string(i));
		prm[sym] = 0.0;
	end

	# individual times taking up residence in building
	@inbounds for i=1:nmax
		sym = Symbol("te"*string(i));
		prm[sym] = 0.0;
	end

	# individual times exiting residence in building
	@inbounds for i=1:nmax
		sym = Symbol("tℓ"*string(i));
		prm[sym] = Inf;
	end

	# flag to say if t₀<=tₑ or t₀<=tℓ should be enforced. Useful bc
	# different calibrations and fittings require different scenarios.
	# Unlike other values, this is not used in prmvary and prmrg bc it
	# is not varied during mcmc
	#  "0.0"=>No relation enforced
	#  "-1.0"=> t₀<=tₑ
	#  "1.0"=> t₀<=tℓ
	prm[:flagt] = -1.0;	

	# λ-params
	#  Γ-hyperparameters for shedding amplitude
	#   Γ(α,β) ~ β^α*x^{α-1}exp(-βx)/Γ(α)
	#   mean: α/β
	#   var:  α/β^2
	prm[:Γα] = 1.0;
	prm[:Γβ] = 0.00875;

	#  Normal-hyperparameters for increment of pos peak rel inf time
	prm[:Aₓμ] = 4.2;
	prm[:Aₓσ] = 0.5;

	#  Normal-hyperparameters for increment of duration rel pos peak
	prm[:Lμ] = 7.3;
	prm[:Lσ] = 0.6;

	#  shedding amplitude
	@inbounds for i=1:nmax
		sym = Symbol("A"*string(i));
		prm[sym] = 1.7142857;
	end

	#  shedding amplitude position rel an inf at time 0
	@inbounds for i=1:nmax
		sym = Symbol("Aₓ"*string(i));
		prm[sym] = 3.5;
	end

	#  shedding duration rel an inf at time 0
	@inbounds for i=1:nmax
		sym = Symbol("L"*string(i));
		prm[sym] = 7.0;
	end

	# particle decay rate
	prm[:ξ] = log(2)/7;

	# Time after time 0 at which dust collected
	prm[:T] = 10.0;

	# dust measurement copies/mg dust in each building
	#  Fall 2020 Iso: bag 1 145.598
	#                 bag 2 199.85
	#                 bag 3 283.398
	#                 bag 4 3226.79
	prm[:Y1] = 172.724;
	prm[:Y2] = 150.0;

	# flag to say if any deterministic relationships exist
	# between parameters, used in some synthetic fitting
	# cases. Like flagt, this parameter is not varied or used
	# in mcmcrg
	# "0.0"=>No deterministic relationships
	# "1.0"=>Some deterministic relationships
	prm[:flagdet] = 0.0;
	
	vkeys = [k for k in keys(prm)];
	return prm,vkeys
end

# data!
""" 
Enfore any deterministic relations between parameters that match the synthetic
data generation case. Only used if flag below is set to true.
WARNING: Determined parameters should NOT be varied in corresponding mcmc, ie
         prmvary in mcmcrg should be false
WARNING: Current limitation is that proposed parameters should NOT have a conditional
         dependence on the determined parameters unless you have modified the proposal
	 and initialization codes explicitly to accommodate.
"""
function data!(prm::Dict{Symbol,Float64})
	if prm[:flagdet]!=1.0
		return
	end

	# Individuals should leave 10 days after entering iso
	nmax = Int64(prm[:nmax]);
	@inbounds for i=1:nmax
		tei = Symbol(:te,i);
		tℓi = Symbol(:tℓ,i);

		prm[tℓi] = prm[tei]+10.0;
	end
end

# mcmcrg
"""
Return dictionaries with bounding intervals for parameters and bools to 
say which parameters are varied
"""
function mcmcrg()
	prmrg = Dict{Symbol,Vector{Float64}}();
	prmvary = Dict{Symbol,Bool}();
	
	# max permitted number of infected people across all buildings
	#  nmax should agree with what is in data and not be varied
	prmrg[:nmax] = [80.0,80.0]; nmax::Int64 = prmrg[:nmax][2];
	prmvary[:nmax] = false;

	# max number of buildings
	#  nbld should agree with what is in data and not be varied
	prmrg[:nbld] = [2.0,2.0]; nbld::Int64 = prmrg[:nbld][2];
	prmvary[:nbld] = false;

	# number of infected people in each building
	#  bound by nmax enforced in prior
	@inbounds for i=1:nbld
		ni = Symbol(:n,i);
		prmrg[ni] = [1.0,200.0];
		prmvary[ni] = false;
	end

	# individual infection times
	@inbounds for i=1:nmax
		sym = Symbol("t"*string(i));
		prmrg[sym] = [-10.0,10.0]; # max t₀ should be less than the fixed T-value
					 # for efficient sampling
		prmvary[sym] = true;
	end

	# individuals taking up residence in building times
	@inbounds for i=1:nmax
		sym = Symbol("te"*string(i)); #shedλ will integrate curve from 
		                              #[max(0,te),T] and vacate integral
					      #if te>T
		prmrg[sym] = [-9.0,10.0];
		prmvary[sym] = false;
	end

	# individual times exiting residence in building
	@inbounds for i=1:nmax
		sym = Symbol("tℓ"*string(i));
		prmrg[sym] = [-9.0,10.0];
		prmvary[sym] = false;
	end

	# λ-params # maybe 50% pickup in dorms by vacuum
	#  Γ-distribution hyperparameters for amplitude
	prmrg[:Γα] = [0.0,25.0];
	prmvary[:Γα] = true;

	prmrg[:Γβ] = [0.0,12.5];
	prmvary[:Γβ] = true;

	#  Normal-distribution hyperparameters for Aₓ increment
	prmrg[:Aₓμ] = [3.0,5.0];
	prmvary[:Aₓμ] = false;

	prmrg[:Aₓσ] = [0.0,1.0];
	prmvary[:Aₓσ] = false;
	
	#  Normal-distribution hyperparameters for L increment
	prmrg[:Lμ] = [6.0,8.0];
	prmvary[:Lμ] = false;

	prmrg[:Lσ] = [0.0,0.25];
	prmvary[:Lσ] = false;


	#  shedding amplitude
	@inbounds for i=1:nmax
		sym = Symbol("A"*string(i));
		prmrg[sym] = [0.0,Inf];
		prmvary[sym] = true;
	end

	#  shedding amplitude position
	#   presently prp! uses a Unif(Tri) for LxAx and could
	#   have problem if L is varied while Ax is not. Similarly
	#   min permitted value of Ax should be <= min permitted
	#   value of L
	@inbounds for i=1:nmax
		sym = Symbol("Aₓ"*string(i));
		prmrg[sym] = [-Inf,Inf];
		prmvary[sym] = true;
	end

	#  shedding duration 
	@inbounds for i=1:nmax
		sym = Symbol("L"*string(i));
		prmrg[sym] = [-Inf,Inf];
		prmvary[sym] = true;
	end

	# particle decay rate
	prmrg[:ξ] = log(2)./[7.0,14.0];
	prmvary[:ξ] = false;

	# Time after time 0 at which dust is collected
	prmrg[:T] = [5.0,10.0];
	prmvary[:T] = false; # Current code needs false because prp doesn't vary

	# dust measurement copies/mg dust in each building
	@inbounds for i=1:nbld
		Yi = Symbol(:Y,i);
		prmrg[Yi] = [0.0,1000.0]; # with Delta numbers over 1000 can go up to 10,000
		prmvary[Yi] = false; # Current code needs false becase prp doesn't vary
	end

	@assert (!prmvary[:nmax])&&(!prmvary[:nbld])&&(!prmvary[:T])&&(!prmvary[:Y1]) "illegal parameter varied in mcmcrg"

	return prmrg,prmvary
end

# shedλI₁
""" 
Ancillary routine used by shedλ to compute integrals of the first type, whose
domain of integration is in rising phase of shedding function. The integral
is
∫ᵇₐexp(-ξ*(T-t))*(exp((t-t₀)/Aₓ*ln(A+1))-1)dt
"""
function shedλI₁(a::Float64,b::Float64,
		 A::Float64,L::Float64,t₀::Float64,T::Float64,Aₓ::Float64,
		 ξ::Float64,tₑ::Float64,tℓ::Float64)
	val = 0.0;
	if b<a
		return val
	end

	η = log(A+1)/Aₓ+ξ;
	val += exp(η*b)-exp(η*a);
	val *= exp(-log(A+1)/Aₓ*t₀)/η;

	val += -1/ξ*(exp(ξ*b)-exp(ξ*a));

	val *= exp(-ξ*T);

	return val
end

# shedλI₂
"""
Ancillary routine used by shedλ to compute integrals of the second type, whose
domain of integration is in the declining phase of shedding function. The integral
is
∫ᵇₐexp(-ξ*(T-t))*( exp(( 1-(t-t₀-Aₓ)/(L-Aₓ) )*ln(A+1)) -1 )dt
"""
function shedλI₂(a::Float64,b::Float64,
		 A::Float64,L::Float64,t₀::Float64,T::Float64,Aₓ::Float64,
		 ξ::Float64,tₑ::Float64,tℓ::Float64)
	val = 0.0;
	if b<a
		return val
	end

	η = ξ-log(A+1)/(L-Aₓ);
	val += exp(η*b)-exp(η*a);
	val *= exp(log(A+1)*(1+(t₀+Aₓ)/(L-Aₓ)))/η;
	
	val += -1/ξ*(exp(ξ*b)-exp(ξ*a));

	val *= exp(-ξ*T);

	return val
end

# shedλ
"""
Compute the shedding μp = ∫_D exp[-ξ*(T-t)]λ(t-t₀;θ)dt 
as function of input parameters, where dust is collected at time T. Multiple 
dispatch for case of single param values and a dictionary. Also include a 
mutating version of the dictionary case for mem alloc. DONT DO flagsynth bc
priors are blind to it
"""
function shedλ(A::Float64,L::Float64,t₀::Float64,T::Float64,Aₓ::Float64,
	       ξ::Float64,tₑ::Float64,tℓ::Float64)
	
	val = 0.0;
	if tₑ<tℓ
		# Type 1 integral over D = [max(t₀,tₑ),min(T,tℓ)]∩[t₀,t₀+Aₓ]
		a = tₑ>=t₀ ? tₑ : t₀; b = (T<tℓ ? T : tℓ)<t₀+Aₓ ? (T<tℓ ? T : tℓ) : t₀+Aₓ;  
		val += shedλI₁(a,b,
			       A,L,t₀,T,Aₓ,
			       ξ,tₑ,tℓ);

		# Type 2 integral over D = [max(t₀,tₑ),min(T,tℓ)]∩[t₀+Aₓ,t₀+L]
		a = tₑ>=t₀+Aₓ ? tₑ : t₀+Aₓ; b = (T<tℓ ? T : tℓ)<t₀+L ? (T<tℓ ? T : tℓ) : t₀+L;
		val += shedλI₂(a,b,
			       A,L,t₀,T,Aₓ,
			       ξ,tₑ,tℓ);
	elseif tℓ<tₑ
		# Type 1 integral over D = [t₀,min(T,tℓ)]∩[t₀,t₀+Aₓ]
		a = t₀; b = (T<=tℓ ? T : tℓ)<=t₀+Aₓ ? (T<=tℓ ? T : tℓ) : t₀+Aₓ;
		val += shedλI₁(a,b,
			       A,L,t₀,T,Aₓ,
			       ξ,tₑ,tℓ);

		# Type 2 integral over D = [t₀,min(T,tℓ)]∩[t₀+Aₓ,t₀+L]
		a = t₀+Aₓ; b = (T<=tℓ ? T : tℓ)<=t₀+L ? (T<=tℓ ? T : tℓ) : t₀+L;
		val += shedλI₂(a,b,
			       A,L,t₀,T,Aₓ,
			       ξ,tₑ,tℓ);

		# Type 1 integral over D = [max(tₑ,t₀),T]∩[t₀,t₀+Aₓ]
		a = tₑ>=t₀ ? tₑ : t₀; b = T<=t₀+Aₓ ? T : t₀+Aₓ;
		val += shedλI₁(a,b,
			       A,L,t₀,T,Aₓ,
			       ξ,tₑ,tℓ);

		# Type 2 integral over D = [max(tₑ,t₀),T]∩[t₀+Aₓ,t₀+L]
		a = tₑ>=t₀ ? tₑ : t₀; b = T<=t₀+L ? T : t₀+L;
		val += shedλI₂(a,b,
			       A,L,t₀,T,Aₓ,
			       ξ,tₑ,tℓ);
	end

	return val > 0 ? val : 1e-16
end
function shedλ!(prm::Dict{Symbol,Float64};
		λval::Vector{Float64}=Vector{Float64}(undef,Int64(prm[:nmax])));
	@inbounds for i=1:Int64(prm[:nmax])
		λval[i] = shedλ(prm[Symbol(:A,i)],prm[Symbol(:L,i)],
				prm[Symbol(:t,i)],prm[:T],
				prm[Symbol(:Aₓ,i)],
				prm[:ξ],
				prm[Symbol(:te,i)],
				prm[Symbol(:tℓ,i)]);
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
"""
function logπ!(prm::Dict{Symbol,Float64},
	       prmrg::Dict{Symbol,Vector{Float64}},prmvary::Dict{Symbol,Bool};
	       λval::Vector{Float64}=Vector{Float64}(undef,Int64(prm[:nmax])))	
	
	nbld::Int64 = floor(prm[:nbld]); nmax::Int64 = floor(prm[:nmax]);
	# Indicator prior on all params save Aₓ,L,A
	@inbounds for key in keys(prmvary)
		if prmvary[key]&&( (prm[key]<prmrg[key][1])||(prm[key]>prmrg[key][2]) )
			return -Inf
		end
	end
	
	#  Indicator for t₀<=T
	if prmvary[Symbol(:t1)]
		@inbounds for i=1:nmax
			ti = Symbol(:t,i);
			if prm[ti]>prm[:T]
				return -Inf
			end
		end
	end

	#  Indicator for t₀<=te if specified in flagt
	if (prm[:flagt]==-1.0)&&(prmvary[Symbol(:t1)]||prmvary[Symbol(:te1)])
		@inbounds for i=1:nmax
			ti = Symbol(:t,i);
			tei = Symbol(:te,i);
			if ti>tei
				return -Inf
			end
		end
	end

	# Indicator for t₀<=tℓ
	if (prm[:flagt]==1.0)&&(prmvary[Symbol(:t1)]||prmvary[Symbol(:tℓ1)])
		@inbounds for i=1:nmax
			ti = Symbol(:t,i);
			tℓi = Symbol(:tℓ,i);
			if ti>tℓi
				return -Inf
			end
		end
	end

	# Indicator for te<=tℓ
	if prmvary[Symbol(:te1)]||prmvary[Symbol(:tℓ1)]
		@inbounds for i=1:nmax
			tei = Symbol(:te,i);
			tℓi = Symbol(:tℓ,i);
			if prm[tei]>prm[tℓi]
				return -Inf
			end
		end
	end
	
	# Conditional priors on Ax, L, A
	val1 = 0.0;
	#  Ax|μ,σ is N(Aₓμ,Aₓσ)
	#@inbounds for i=1:nmax
	#	Axi = Symbol(:Aₓ,i);
	#	ti = Symbol(:t,i);
	#	val1 += lognrm(prm[:Aₓμ],prm[:Aₓσ],prm[Axi]);
	#end

	# L|Aₓ,μ,σ is N(Lμ+Aₓ,Lσ)
	#@inbounds for i=1:nmax
	#	Li = Symbol(:L,i);
	#	Axi = Symbol(:Aₓi);
	#	val1 += lognrm(prm[:Lμ]+prm[Axi],prm[:Lσ],prm[Li]);
	#end

	# A|Γα,Γβ is Γ(α,β)
	#@inbounds for i=1:nmax
	#	Ai = Symbol(:A,i);
	#	val1 += logΓ(prm[:Γα],prm[:Γβ],prm[Ai]);
	#end

	# Likelihood
	shedλ!(prm;λval=λval);
	val2 = 0.0; pos = 0;
	@inbounds for i=1:nbld
		ni = Symbol(:n,i); n::Int64 = floor(prm[ni]);
		Yi = Symbol(:Y,i); Y = prm[Yi];
		val = 0.0;
		@inbounds for j=pos+1:pos+n
			val += λval[j];
		end
		val2 += -val + floor(Y)*log(val);
		pos += n;
	end
	
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
	nbld::Int64 = floor(prm[:nbld]); nmax::Int64 = floor(prm[:nmax]);	
	# prp density on n unif
	if prmvary[:n1]
		@inbounds for i=1:nbld
			ni = Symbol(:n,i);
			prm[ni] = prmrg[ni][1]+rand(rng)*(prmrg[ni][2]-prmrg[ni][1])
		end
	end

	# prp density on ξ is uniform
	if prmvary[:ξ]
		prm[:ξ] = prmrg[:ξ][1]+rand(rng)*(prmrg[:ξ][2]-prmrg[:ξ][1]);
	end

	# prp density of tₑ's is unif
	if prmvary[:te1]
		@inbounds for i=1:nmax
			tei = Symbol(:te,i);
			tℓi = Symbol(:tℓ,i);
			if !prmvary[:tℓ1]
				temax = prm[tℓi] <= prmrg[tei][2] ? prm[tℓi] : prmrg[tei][2];
			else
				temax = prmrg[tei][2];
			end
			prm[tei] = prmrg[tei][1]+rand(rng)*(temax-prmrg[tei][1]);
		end
	end

	# prp density on tℓ's conditioned on tₑ is unif from tₑ to its upper rg value
	if prmvary[:tℓ1]
		@inbounds for i=1:nmax
			tei = Symbol(:te,i);
			tℓi = Symbol(:tℓ,i);
			tℓlow = prm[tei] >= prmrg[tℓi][1] ? prm[tei] : prmrg[tℓi][1]
			prm[tℓi] = tℓlow + rand(rng)*(prmrg[tℓi][2]-tℓlow)
		end
	end

	# prp density on t₀'s is unif potentially conditioned on (tₑ,tℓ)
	if prmvary[:t1]
		@inbounds for i=1:nmax
			ti = Symbol(:t,i);
			if prm0[:flagt]==-1.0
				tei = Symbol(:te,i);
				tupr = prm[tei]<=prmrg[ti][2] ? prm[tei] : prmrg[ti][2];
			elseif prm0[:flagt]==1.0
				tℓi = Symbol(:tℓ,i);
				tupr = prm[tℓi]<=prmrg[ti][2] ? prm[tℓi] : prmrg[ti][2];
			else
				tupr = prmrg[ti][2];
			end
			prm[ti] = prmrg[ti][1]+rand(rng)*(tupr-prmrg[ti][1]);
		end
	end

	# prp density on Gamma hypers is random walk
	if prmvary[:Γα]
		ΔΓα = 0.02*(prmrg[:Γα][2]-prmrg[:Γα][1]);
		prm[:Γα] = prm0[:Γα] + ΔΓα*randn(rng);
	end
	if prmvary[:Γβ]
                ΔΓβ = 0.02*(prmrg[:Γβ][2]-prmrg[:Γβ][1]);
                prm[:Γβ] = prm0[:Γβ] + ΔΓβ*randn(rng);
        end
	# Patched MH rejection
	#  Point is if alpha and beta aren't in cnst reg, the chain automatically rejects prm
	#  and resets to prm0. Used bc Julia wont sample a Gamma with nonpositive hyperparams
	if (prm[:Γα]<=0)||(prm[:Γβ]<=0)
		@inbounds for key in keys(prmvary)	
			prm[key] = prm0[key];
		end
		return
	end
	
	# prp density on Normal hypers is random walk
	if prmvary[:Aₓμ]
		ΔAₓμ = 0.02*(prmrg[:Aₓμ][2]-prmrg[:Aₓμ][1]);
		prm[:Aₓμ] = prm0[:Aₓμ] + ΔAₓμ*randn(rng);
	end
	if prmvary[:Aₓσ]
		ΔAₓσ = 0.02*(prmrg[:Aₓσ][2]-prmrg[:Aₓσ][1]);
		prm[:Aₓσ] = prm0[:Aₓσ]+ΔAₓσ*randn(rng);
	end

	if prmvary[:Lμ]
		ΔLμ = 0.02*(prmrg[:Lμ][2]-prmrg[:Lμ][1]);
		prm[:Lμ] = prm0[:Lμ] + ΔLμ*randn(rng);
	end
	if prmvary[:Lσ]
		ΔLσ = 0.02*(prmrg[:Lσ][2]-prmrg[:Lσ][1]);
		prm[:Lσ] = prm0[:Lσ]+ΔLσ*randn(rng);
	end

	# prp density of amplitudes conditioned on hypers is Gamma
	if prmvary[:A1]
		Γdistr = Gamma(prm[:Γα],1/prm[:Γβ]);
		@inbounds for i=1:nmax
			Ai = Symbol(:A,i);
			prm[Ai] = rand(rng,Γdistr);
		end
	end

	# prp density of pos of peak rel inf time is normal
	if prmvary[:Aₓ1]
		@inbounds for i=1:nmax
			ti = Symbol(:t,i);
			Aₓi = Symbol(:Aₓ,i);
			prm[Aₓi] = prm[:Aₓμ]+prm[:Aₓσ]*randn(rng);
		end
	end

	# prp density of length of shedding conditioned on Aₓ is normal
	if prmvary[:L1]
		@inbounds for i=1:nmax
			Aₓi = Symbol(:Aₓ,i);
			Li = Symbol(:L,i);
			prm[Li] = prm[Aₓi]+prm[:Lμ]+prm[:Lσ]*randn(rng);
		end
	end

	data!(prm);
end

# logρ!
"""
Evaluate the log unnormalized proposal density ρ(y|x) needed in mhratio
for the subset of parameters being varied. Only the ρ(A|L) normalizing
factor and Gamma's need to be tracked since the other factors are either
prop to indicators and cancel in the ratio or are the random walk which
is symmetric depending only on |x-y| which is same y|x or x|y.
"""
function logρ!(prm0::Dict{Symbol,Float64},prm::Dict{Symbol,Float64},
	       prmrg::Dict{Symbol,Vector{Float64}},prmvary::Dict{Symbol,Bool};
	       λval::Vector{Float64}=Vector{Float64}(undef,Int64(prm[:nmax][1])))
	val = 0.0;

	nmax::Int64 = floor(prm[:nmax]);

	if prmvary[:te1]&&prmvary[:tℓ1]
		@inbounds for i=1:nmax
			tei = Symbol(:te,i);
			tℓi = Symbol(:tℓ,i);
			tℓlow = prm[tei] >= prmrg[tℓi][1] ? prm[tei] : prmrg[tℓi][1];

			val += -log(prmrg[tℓi][2]-tℓlow);
		end
	end

	if (prm0[:flagt]!=0.0)&&(prmvary[:t1])
		@inbounds for i=1:nmax
			ti = Symbol(:t,i);
			teℓi = prm0[:flagt]==-1.0 ? Symbol(:te,i) : Symbol(:tℓ,i);
			tupr = prm[teℓi]<=prmrg[ti][2] ? prm[teℓi] : prmrg[ti][2];

			val += -log(tupr - prmrg[ti][1])
		end
	end
	
	#if prmvary[:A1]
	#	@inbounds for i=1:nmax
	#		Ai = Symbol(:A,i);
	#		val += logΓ(prm[:Γα],prm[:Γβ],prm[Ai]);
	#	end
	#end
	
	#  Ax|t₀ is N(Aₓμ+t₀,Aₓσ)
	#@inbounds for i=1:nmax
	#	Axi = Symbol(:Aₓ,i);
	#	ti = Symbol(:t,i);
	#	val1 += lognrm(prm[:Aₓμ],prm[:Aₓσ],prm[Axi]);
	#end

	# L|Aₓ is N(Lμ+Aₓ,Lσ)
	#@inbounds for i=1:nmax
	#	Li = Symbol(:L,i);
	#	Axi = Symbol(:Aₓi);
	#	val1 += lognrm(prm[:Lμ]+prm[Axi],prm[:Lσ],prm[Li]);
	#end

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
	flagfd = false; nmax::Int64 = prm[:nmax];
	while !flagfd
		@inbounds for key in keys(prmvary)
			if prmvary[key]
				vlow = prmrg[key][1]==-Inf ? -3.0 : prmrg[key][1];
				vhgh = prmrg[key][2]==Inf ? 3.0 : prmrg[key][2];
				prm[key] = vlow + rand(rng)*(vhgh-vlow);

			end
		end

		if prmvary[:tℓ1]
			nmax = Int64(prm[:nmax]);
			@inbounds for i=1:nmax
				tℓi = Symbol(:tℓ,i);
				tei = Symbol(:te,i);
				vhgh = prmrg[tℓi][2]==Inf ? prm[tei]+10.0 : prmrg[tℓi][2];
				prm[tℓi] = prm[tei] + rand(rng)*(vhgh-prm[tei])
			end
		end

		if (prm[:flagt]!=0.0)&&(prmvary[:t1])
			@inbounds for i=1:nmax
				ti = Symbol(:t,i);
				teℓi = prm[:flagt]==-1.0 ? Symbol(:te,i) : Symbol(:tℓ,i);
				tupr = prm[teℓi]<=prmrg[ti][2] ? prm[teℓi] : prmrg[ti][2];

				prm[ti] = prmrg[ti][1]+rand(rng)*(tupr-prmrg[ti][1]);
			end
		end

		data!(prm);
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
	@inbounds for i=1:ncyc
		prp!(prm0,prm,prmrg,prmvary;λval=λval,rng=rng)

		coin = rand(rng) |> log;
		if coin <=logmh!(prm0,prm,prmrg,prmvary;λval=λval)
			# accept
			@inbounds for key in keys(prmvary)
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
		df0 = CSV.read("MCMCsmp.csv",DataFrame,header=false);
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
	println("Progress through mcmc: $prg/1 ...");
	@inbounds for i=1:nsmp
		mhrej[i] = mcmcsmp!(prm0,prm,prmrg,prmvary;
			            λval=λval,rng=rng,ncyc=ncyc);
		smp = @view SMP[:,i];
		wrtprm!(prm0,vkeys,smp);

		while i/nsmp>=prg+Δprg
			prg+=Δprg;
			println("Progress through mcmc: $prg/1 ...");
			CSV.write("MCMCsmp.csv",[DataFrame(:prm=>String.(vkeys)) DataFrame(SMP[:,1:i])], writeheader=false,append=false);
			CSV.write("rejstats.csv",DataFrame(:rejct=>mhrej));
			mysaverng(rng);
		end
	end

	# Save final csv and report rejection rates
	println("Progress through mcmc: 1.0/1 ...");
	CSV.write("MCMCsmp.csv",[DataFrame(:prm=>String.(vkeys)) DataFrame(SMP[:,1:nsmp])], writeheader=false,append=false);
	CSV.write("rejstats.csv",DataFrame(:rejct=>mhrej),append=false);
	mysaverng(rng);

	rejrt = sum(mhrej)/(ncyc*nsmp); aptrt = 1-rejrt; aptwt = 1/aptrt;
	println("Rejection rate: $rejrt");
	println("Average num proposals before an accept: $aptwt");
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
	@inbounds for i=1:nelm
		V[i] = prm[vkeys[i]];
	end
		       
	return prm,vkeys,V
end
function wrtprm!(prm::Dict{Symbol,Float64},vkeys::Vector{Symbol},
                   V::VecVw)

	@inbounds for i=1:length(vkeys)
		V[i] = prm[vkeys[i]];
	end
	
end
function wrtprm!(prm1::Dict{Symbol,Float64},
		 prm2::Dict{Symbol,Float64})
	@inbounds for key in keys(prm1)
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
	@inbounds for i=1:length(vkeys)
		prm[vkeys[i]] = V[i];
	end

	return prm,vkeys
end
