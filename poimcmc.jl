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
	
	# cap of fall iso and number of infected people in building
	#  Should be strictly more than niso + n
	prm[:nmax] = 101.0; nmax = Int64(prm[:nmax]);

	# number of people in iso
	prm[:niso] = 1.0;

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
		prm[sym] = 1.0;
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
	prm[:Yiso] = 172.0;
	
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

	# number of people in iso
	prmrg[:niso] = [39.0,40.0];
	prmvary[:niso] = false;

	# number of infected people in the building
	#  bound by nmax enforced in prior
	prmrg[:n] = [1.0,100.0]; # For rej stats have be 1 less than nmax
	prmvary[:n] = true;

	# individual infection times
	for i=1:nmax
		sym = Symbol("t"*string(i));
		prmrg[sym] = [-14.0,0.0];
		prmvary[sym] = false;
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
		prmrg[sym] = [0.0,10000.0];
		prmvary[sym] = true;
	end

	#  shedding amplitude position
	for i=1:nmax
		sym = Symbol("Aₓ"*string(i));
		prmrg[sym] = [0.0,7.0];
		prmvary[sym] = true;
	end

	#  shedding duration 
	for i=1:nmax
		sym = Symbol("L"*string(i));
		prmrg[sym] = [7.0,14.0];
		prmvary[sym] = false;
	end

	# particle decay rate
	prmrg[:ξ] = log(2)./[7.0,14.0];
	prmvary[:ξ] = false;

	# Time after time 0 at which dust is collected
	prmrg[:T] = [5.0,10.0];
	prmvary[:T] = false;

	# dust measurement copies/mg dust
	prmrg[:Y] = [0.0,1000.0]; # with Delta numbers over 1000 can go up to 10,000
	prmvary[:Y] = false;

	prmrg[:Yiso] = [172.0,173.0];
	prmvary[:Yiso] = false;

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
Compute the shedding ∫ᵀ₀exp[-ξ*(T-t)]λ(t-t₀;θ)dt as function of input parameters
Multiple dispatch for case of single param values and a dictionary.
Also include a mutating version of the dictionary case for mem alloc
"""
function shedλ(A::Float64,L::Float64,t₀::Float64,T::Float64,Aₓ::Float64,
	       ξ::Float64)
	# exported from Maple
	cg = Aₓ; cg1 = L; cg3 = T; cg5 = t₀; xi = ξ;

	λval = -A * (-cg1 * ((-cg3 + cg + cg5 < 0 ? 0 : 1) - (cg + cg5 < 0 ? 0 : 1)) * exp((xi * (-cg3 + cg + cg5))) + cg * ((-cg3 + cg1 + cg5 < 0.0e0 ? 0 : 1) - (cg1 + cg5 < 0.0e0 ? 0 : 1)) * exp(xi * (-cg3 + cg1 + cg5)) - (-xi * cg - 1 + (cg3 - cg5) * xi) * cg1 * (-cg3 + cg + cg5 < 0 ? 0 : 1) - cg * (xi * cg1 + 0.1e1 + ((-cg3 + cg5) * xi)) * (-cg3 + cg1 + cg5 < 0.0e0 ? 0 : 1) + (-cg1 * (xi * cg + xi * cg5 + 1) * (cg + cg5 < 0 ? 0 : 1) + cg * (cg1 + cg5 < 0.0e0 ? 0 : 1) * (xi * cg1 + (xi * cg5) + 0.1e1)) * exp(-(xi * cg3))) / (xi ^ 2) / (cg1 - cg) / cg;


	return λval
end
function shedλ(prm::Dict{Symbol,Float64})
	λval = Vector{Float64}(undef,Int64(prm[:nmax]));
	for i=1:Int64(prm[:nmax])
		λval[i] = shedλ(prm[Symbol(:A,i)],prm[Symbol(:L,i)],
				prm[Symbol(:t,i)],prm[:T],
				prm[Symbol(:Aₓ,i)],
				prm[:ξ]);
	end

	return λval
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
	# Bounding box prior
	for key in keys(prmvary)
		if prmvary[key]&&( sum( (prm[key]<prmrg[key][1])+(prm[key]>prmrg[key][2]) ) !=0 )
			return -Inf
		end
	end

	# Likelihood
	niso = Int64(floor(prm[:niso])); n = Int64(floor(prm[:n])); nmax = Int64(floor(prm[:nmax]));
	val = 0.0;
	if flagλval
		shedλ!(prm;λval=λval);
	end
	for i=(niso+1):(niso+n)
		val += λval[i];
	end
	val = -val + floor(prm[:Y])*log(val);	
	
	# Prior calibrated from fall isolation data
#	val2 = 0.0;
#	for i=1:niso
#		val2 += λval[i];
#	end
#	val2 = -val2 + floor(prm[:Yiso])*log(val2);

#	val += val2;

	# Priors on shedding amplitude
	for i=1:nmax
		Ai = Symbol(:A,i);
		val += logΓ(prm[:Γα],prm[:Γβ],prm[Ai]);
	end
	
	return val
end

# prp!
"""
Metropolis proposal function
"""
function prp!(prm0::Dict{Symbol,Float64},prm::Dict{Symbol,Float64},
		  prmrg::Dict{Symbol,Vector{Float64}},prmvary::Dict{Symbol,Bool};
		  λval::Vector{Float64}=Vector{Float64}(undef,Int64(prm[:nmax])),
		  rng::MersenneTwister=MersenneTwister(),
		  key::Symbol=:ALL,
		  uenv::Float64=1.0)
	if prmvary[:Γα]
		ΔΓα = 0.02*(prmrg[:Γα][2]-prmrg[:Γα][1]);
		prm[:Γα] = prm0[:Γα] + Δα*randn(rng);
	end

	if prmvary[:Γβ]
                ΔΓβ = 0.02*(prmrg[:Γβ][2]-prmrg[:Γβ][1]);
                prm[:Γβ] = prm0[:Γβ] + Δβ*randn(rng);
        end

	# Double-check this patched-MH rejection is valid
	#  Point is if alpha and beta aren't in cnst reg, the chain automatically rejects it
	#  Used since nonpositive values cause Gamma sampling to fail
	if (prm[:Γα]<=0)||(prm[:Γβ]<=0)
		for key in keys(prm0)	
			prm[key] = prm0[key];
		end
		return
	end

	niso = Int64(floor(prm[:niso])); n = Int64(floor(prm[:n])); nmax = Int64(floor(prm[:nmax]));
	if prmvary[:A1]
		Γdistr = Gamma(prm[:Γα],1/prm[:Γβ]);
		for i=1:nmax
			Ai = Symbol(:A,i);
			prm[Ai] = rand(rng,Γdistr);
		end
	end

	if prmvary[:Aₓ1]
		for i=1:nmax
			Axi = Symbol(:Aₓ,i);
			prm[Axi] = prmrg[Axi][1]+rand(rng)*(prmrg[Axi][2]-prmrg[Axi][1]);
		end
	end

	if prmvary[:L1]
		for i=1:nmax
			Li = Symbol(:L,i);
			prm[Li] = prmrg[Li][1]+rand(rng)*(prmrg[Li][2]-prmrg[Li][1]);
		end
	end

	if prmvary[:n]
		prm[:n] = prmrg[:n][1]+rand(rng)*(prmrg[:n][2]-prmrg[:n][1])
	end

end

# logρ!
"""
Evaluate the log unnormalized proposal density ρ(y|x)
for the subset of parameters being varied
"""
function logρ!(prm0::Dict{Symbol,Float64},prm::Dict{Symbol,Float64},
	       prmrg::Dict{Symbol,Vector{Float64}},prmvary::Dict{Symbol,Bool};
	       λval::Vector{Float64}=Vector{Float64}(undef,Int64(prm[:nmax][1])),
	       flagλval::Bool=true)
	val = 0.0;

	if prmvary[:Γα]
		ΔΓα = 0.02*(prmrg[:Γα][2]-prmrg[:Γα][1]);
		val += lognrm(prm0[:Γα],ΔΓα,prm[:Γα]);
	end

	if prmvary[:Γβ]
		ΔΓβ = 0.02*(prmrg[:Γβ][2]-prmrg[:Γβ][1]);
		val += lognrm(prm0[:Γβ],ΔΓβ,prm[:Γβ]);
	end
	
	niso = Int64(floor(prm[:niso])); n = Int64(floor(prm[:n])); nmax = Int64(floor(prm[:nmax]));
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
					           prmrg[key][2]-prmrg[key][1]
					                                     );
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
	        flagcase::Symbol=:glbl,
	       λval::Vector{Float64}=Vector{Float64}(undef,Int64(prm[:nmax])))
	
	# stationary distribution
	val  = logπ!(prm,prmrg,prmvary;λval=λval);
	val -= logπ!(prm0,prmrg,prmvary;λval=λval);

	# proposal distribution
	val += logρ!(prm,prm0,prmrg,prmvary;λval=λval);
	val -= logρ!(prm0,prm,prmrg,prmvary;λval=λval);

	return val
end

# mcmcsmp
"""
mcmc sample the posterior distribution
"""
function mcmcsmp(nsmp::Int64;
		 rng::MersenneTwister=MersenneTwister(),
		 flagrst::Bool=false,
		 ncyc::Int64=1)
	prm,vkeys,V = wrtprm(); 
	prmrg,prmvary = mcmcrg();
	λval = Vector{Float64}(undef,Int64(prm[:nmax]));
	
	if flagrst
		# restart sampling from csv
		df0 = CSV.read("GibbsMCMC.csv",DataFrame,header=false);
		vkeys = Symbol.(df0[:,1]); V = df0[:,end];
		prm0,_ = rdprm(V,vkeys); prm = deepcopy(prm0); 
		rng = myloadrng();

		df0 = similar(df0,0);
	else
		init!(prm,prmrg,prmvary;λval=λval,rng=rng);
		prm0 = deepcopy(prm);
	end

	# Create matrix for samples
	SMP = Matrix{Float64}(undef,length(V),nsmp);

	# Create vectors for storing rejection statistics	
	mhcnt = Dict{Symbol,Vector{Int64}}(:pos=>[1],:rjt=>fill(0,nsmp),
					   :tot=>fill(0,nsmp));

	# Preallocate vectors for prm symb's that aren't used in likelihood
	# due to varying number of indiv's. Used in computing rej stats.
	Asymb = [Symbol(:A,i) for i=1:Int64(prm[:nmax])];
	Aₓsymb = [Symbol(:Aₓ,i) for i=1:Int64(prm[:nmax])];
	Lsymb = [Symbol(:L,i) for i=1:Int64(prm[:nmax])];
	psymb = [Symbol(:p,i) for i=1:Int64(prm[:nmax])];

	# run mcmc
	pos = 0.0; Δprg = 0.02;
	gen=[1,0];
	for k=1:nsmp*ncyc
		# Cycle generator
		if gen[2]!=ncyc
			gen[2]+=1;
		else
			gen[1]+=1; gen[2]=1;
		end
		i=gen[1]; j=gen[2];

		#  Set rejection counter to current position
		mhcnt[:pos][1] = i;

		#  propose
		prp!(prm0,prm,prmrg,prmvary;λval=λval,rng=rng);

		#  mh accept-reject
		if log(rand(rng)) < logmh!(prm0,prm,prmrg,prmvary;λval=λval)
			wrtprm!(prm,prm0);
		else
			wrtprm!(prm0,prm);
			mhcnt[:rjt][mhcnt[:pos][1]]+=1;
		end
			mhcnt[:tot][mhcnt[:pos][1]]+=1;
		
		if j==ncyc
			# Record the sample
			P0 = @view SMP[:,i];
			wrtprm!(prm0,vkeys,P0);

			# Save partial progress
			prg = i/nsmp;
			if prg >= pos + Δprg
				pos = floor(prg/Δprg)*Δprg;
				CSV.write("GibbsMCMC.csv",[DataFrame(:prm=>String.(vkeys)) DataFrame(SMP[:,1:i])], writeheader=false,append=false);
				dftemp = DataFrame(:mhrej=>mhcnt[:rjt][1:i],:mhtot=>mhcnt[:tot][1:i]);
				CSV.write("RejStats.csv",dftemp);
				mysaverng(rng);
				println("$pos" *"/1 complete with MCMC samples ...")
			end
		end
	end

	# Save chain parameter values to csv
	CSV.write("GibbsMCMC.csv", [DataFrame(:prm=>String.(vkeys)) DataFrame(SMP)], writeheader=false, append=false);
	
	# Save rejection statistics to csv
	dftemp = DataFrame(:mhrej=>mhcnt[:rjt],:mhtot=>mhcnt[:tot]);
	CSV.write("RejStats.csv",dftemp);

	# Save random number generator state
	mysaverng(rng);


	return SMP
end
