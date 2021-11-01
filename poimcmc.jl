using Random,CSV,DataFrames
VecVw = Union{Vector{Float64},
	      SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}
	      };

# data
"""
Return dictionary with default parameter values
"""
function data()
	prm = Dict{Symbol,Float64}();
	
	# max permitted number of infected people in building
	#  Should be at or more than 40 bc of Fall Iso using 40
	prm[:nmax] = 100.0; nmax = Int64(prm[:nmax]);

	# number of infected people in the building
	prm[:n] = 10.0;

	# individual infection times
	for i=1:nmax
		sym = Symbol("t"*string(i));
		prm[sym] = 0.0;
	end

	# λ-params
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

	# p-survival params
	for i=1:nmax
		sym = Symbol("p"*string(i));
		prm[sym] = 0.5;
	end

	# Time after time 0 at which dust collected
	prm[:T] = 7.0;

	# dust measurement copies/mg dust
	prm[:Y] = 5.0;
	
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
	prmrg[:nmax] = [0.0,100.0]; nmax = Int64(prmrg[:nmax][2]);
	prmvary[:nmax] = false;

	# number of infected people in the building
	#  bound by nmax enforced in prior
	prmrg[:n] = [0.0,100.0];
	prmvary[:n] = true;

	# individual infection times
	for i=1:nmax
		sym = Symbol("t"*string(i));
		prmrg[sym] = [-14.0,0.0];
		prmvary[sym] = false;
	end

	# λ-params # maybe 50% pickup in dorms by vacuum
	#  shedding amplitude
	for i=1:nmax
		sym = Symbol("A"*string(i));
		prmrg[sym] = [0.0,1.0];
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
		prmvary[sym] = true;
	end

	# p-survival params
	for i=1:nmax
		sym = Symbol("p"*string(i));
		prmrg[sym] = [0.0,1.0]; # lose an order of magnitude a week
		prmvary[sym] = true;
	end

	# Time after time 0 at which dust is collected
	prmrg[:T] = [5.0,10.0];
	prmvary[:T] = false;

	# dust measurement copies/mg dust
	prmrg[:Y] = [0.0,1000.0]; # with Delta numbers over 1000 can go up to 10,000
	prmvary[:Y] = false;

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
Compute the shedding ∫ᵀ₀λ(t-t₀;θ)dt as function of input parameters
Multiple dispatch for case of single param values and a dictionary.
Also include a mutating version of the dictionary case for mem alloc
"""
function shedλ(A::Float64,L::Float64,t0::Float64,T::Float64,Aₓ::Float64)
	# exported from Maple
	if T<=Aₓ+t0
		λval = A*T*(T-2*t0)/(2*Aₓ);
	else
		λval = A*(Aₓ^2*L-2*Aₓ*T*L+2*t0*L*Aₓ+Aₓ*T^2-2*Aₓ*T*t0+t0^2*L)/(2*Aₓ*(-L+Aₓ));
	end

	if Aₓ+t0<0
		λval += -A*L*(Aₓ+t0)^2/(2*Aₓ*(-L+Aₓ));
	end

	return λval
end
function shedλ(prm::Dict{Symbol,Float64})
	λval = Vector{Float64}(undef,Int64(prm[:nmax]));
	for i=1:Int64(prm[:nmax])
		λval[i] = shedλ(prm[Symbol(:A,i)],prm[Symbol(:L,i)],
				prm[Symbol(:t,i)],prm[:T],
				prm[Symbol(:Aₓ,i)]);
	end

	return λval
end
function shedλ!(prm::Dict{Symbol,Float64};
		λval::Vector{Float64}=Vector{Float64}(undef,Int64(prm[:nmax])));
	for i=1:Int64(prm[:nmax])
		λval[i] = shedλ(prm[Symbol(:A,i)],prm[Symbol(:L,i)],
				prm[Symbol(:t,i)],prm[:T],
				prm[Symbol(:Aₓ,i)]);
	end

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
	val = 0.0;
	if flagλval
		shedλ!(prm;λval=λval);
	end
	for i=1:Int64(floor(prm[:n]))
		val += prm[Symbol(:p,i)]*λval[i];
	end
	val = -val + prm[:Y]*log(val);	
	
	# Prior calibrated from fall isolation data
	#  Oct: 40 ppl at 172 copies/mg dust
	#  Nov: 40 ppl at 283 copies/mg dust
	val2 = 0.0;
	for i=1:40
		val2 += prm[Symbol(:p,i)]*λval[i];
	end
	val2 = -val2 + 172*log(val2);

	val += val2;
	
	return val
end

# logρ!
"""
Evaluate the log unnormalized proposal density used for global sampling
Density is prop to
1/[ (∑pᵢμᵢ-y)^2 + 1 ]
which has absolute bound 1
"""
function logρ!(prm::Dict{Symbol,Float64},
	       prmrg::Dict{Symbol,Vector{Float64}},prmvary::Dict{Symbol,Bool};
	       λval::Vector{Float64}=Vector{Float64}(undef,Int64(prm[:nmax][1])),
	       flagλval::Bool=true)
 
	# Bounding box support
        for key in keys(prmvary)
		if prmvary[key]&&( sum( (prm[key]<prmrg[key][1])+(prm[key]>prmrg[key][2]) ) !=0 )
			return -Inf
		end
	end

	# density
	val = 0.0;
        if flagλval
                shedλ!(prm;λval=λval);
	end
	for i=1:Int64(floor(prm[:n]))
		val += prm[Symbol(:p,i)]*λval[i];
	end
	
	val = -log( (val-prm[:Y])^2 + 1 );	

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

# acptrjt!
"""
Accept-reject propose from the global proposal density by uniformly
sampling the subgraph
"""
function acptrjt!(prm::Dict{Symbol,Float64},
		  prmrg::Dict{Symbol,Vector{Float64}},prmvary::Dict{Symbol,Bool};
		  λval::Vector{Float64}=Vector{Float64}(undef,Int64(prm[:nmax])),
		  rng::MersenneTwister=MersenneTwister(),
		  key::Symbol=:ALL,
		  uenv::Float64=1.0,
		  rjtcnt::Dict{Symbol,Vector{Int64}}=Dict{Symbol,Vector{Int64}}(
						      :pos=>[0],:rjt=>[0],:tot=>[0]) )
	if (key!=:ALL)&&(!prmvary[key])
		return
	end

	flagfd = false;
	while !flagfd
		if key!=:ALL
			prm[key] = prmrg[key][1] .+ rand(rng)*(
						       prmrg[key][2]-prmrg[key][1]
						                      );
		else
			for key0 in keys(prmvary)
				if prmvary[key0]
					prm[key0] = prmrg[key0][1] .+ rand(rng)*(
							     prmrg[key0][2]-prmrg[key0][1]
							                                );
				end
			end
		end

		gr = rand(rng)*uenv;
		if log(gr)<logρ!(prm,prmrg,prmvary;λval=λval)
			flagfd = true;
		elseif rjtcnt[:pos][1]!=0
			rjtcnt[:rjt][rjtcnt[:pos][1]] += 1;
		end
	end
end

# ranw!
"""
Random walk propose a new sample
 Note: Assumes that prm is a copy of prm0 except in the entries being varied
"""
function ranw!(prm0::Dict{Symbol,Float64},prm::Dict{Symbol,Float64},
	       prmrg::Dict{Symbol,Vector{Float64}},prmvary::Dict{Symbol,Bool};
	       rng::MersenneTwister=MersenneTwister(),
	       key::Symbol=:ALL,
	       relΔr::Float64=0.025)
	if (key!=:ALL)&&(!prmvary[key])
		return
	end

	if key!=:ALL	
		prm[key] = prm0[key] + randn(rng)*relΔr*(
				           prmrg[key][2]-prmrg[key][1]
				                                     );
	else
		for key0 in keys(prmvary)
			if prmvary[key0]
				prm[key0] = prm0[key0] + randn(rng)*relΔr*(
					                  prmrg[key0][2]-prmrg[key0][1]
					                                              );
			end
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
	if flagcase==:glbl
		val += logρ!(prm0,prmrg,prmvary;λval=λval);
		val -= logρ!(prm,prmrg,prmvary;λval=λval);
	end

	return val
end

# mcmcsmp
"""
mcmc sample the posterior distribution
"""
function mcmcsmp(nsmp::Int64;
		 rng::MersenneTwister=MersenneTwister(),
		 MHmix::Float64=0.1,MGrw::Float64=0.99,
		 flagrst::Bool=false)
	prm,vkeys,V = wrtprm(); 
	prmrg,prmvary = mcmcrg();
	λval = Vector{Float64}(undef,Int64(prm[:nmax]));
	
	if flagrst
		# restart sampling from csv
	else
		init!(prm,prmrg,prmvary;λval=λval,rng=rng);
		prm0 = deepcopy(prm);
	end

	# Create matrix for samples
	SMP = Matrix{Float64}(undef,length(V),nsmp);

	# Create vectors for storing rejection statistics
	rjtcnt = Dict{Symbol,Vector{Int64}}(:pos=>[1],:rjt=>fill(0,nsmp),
					    :tot=>fill(0,nsmp));
	mhcnt = Dict{Symbol,Vector{Int64}}(:pos=>[1],:rjt=>fill(0,nsmp),
					   :tot=>fill(0,nsmp));

	# run mcmc
	pos = 0.0; Δprg = 0.02; nMH = 0; nGibbs = 0;
	for i=1:nsmp
		if rand(rng) < MHmix
			# Global Metropolis-Hastings
			nMH += 1;
			#  glbl propose and record rej stats
			acptrjt!(prm,prmrg,prmvary;λval=λval,rng=rng,
				                   rjtcnt=rjtcnt);
			rjtcnt[:tot][rjtcnt[:pos][1]]+=1;

			#  mh accept-reject
			if log(rand(rng)) < logmh!(prm0,prm,prmrg,prmvary;λval=λval)
				wrtprm!(prm,prm0);
			end
		else
			# Metropolis-within-Gibbs by mixture of glbl prp and rw
			nGibbs += 1;
			for key in vkeys
				if rand(rng) >= MGrw
					# glbl propose and record rej stats
					acptrjt!(prm,prmrg,prmvary;λval=λval,rng=rng,key=key,
						                   rjtcnt=rjtcnt);
					rjtcnt[:tot][rjtcnt[:pos][1]]+=1;

					# mh accept-reject
					#  the resulting yk value stored to prm0 and prm 
					#  for next iteration of Gibbs
					if log(rand(rng)) < logmh!(prm0,prm,prmrg,prmvary;λval=λval)
						prm0[key] = prm[key];
					else
						prm[key] = prm0[key];
					end
				else 
					# rw
					ranw!(prm0,prm,prmrg,prmvary;key=key);

					# mh accept-reject and record rej stats
					#  the resulting yk value stored to prm0 and prm 
					#  for next iteration of Gibbs
					if log(rand(rng)) < logmh!(prm0,prm,prmrg,prmvary;
								  λval=λval,flagcase=:rw)
						prm0[key] = prm[key];
					else
						prm[key] = prm0[key];
						mhcnt[:rjt][mhcnt[:pos][1]]+=1;
					end
					mhcnt[:tot][mhcnt[:pos][1]]+=1;
				end
			end
		end

		# Record the sample
		P0 = @view SMP[:,i];
		wrtprm!(prm0,vkeys,P0); wrtprm!(prm0,prm);

		# Cycle rej stat tracker to next position
		rjtcnt[:pos][1] = i;
		mhcnt[:pos][1] = i;

		# Save partial progress
		prg = i/nsmp;
		if prg >= pos + Δprg
			pos = floor(prg/Δprg)*Δprg;
			CSV.write("GibbsMCMC.csv", DataFrame(SMP[:,1:i]), writeheader=false);
			dftemp = DataFrame(:rjtrej=>rjtcnt[:rjt][1:i],:rjttot=>rjtcnt[:tot][1:i],
					   :mhrej=>mhcnt[:rjt][1:i],:mhtot=>mhcnt[:tot][1:i]);
			CSV.write("RejStats.csv",dftemp);
			mysaverng(rng);
			println("$pos" *"/1 complete with MCMC samples ...")
		end
	end

	# Save chain parameter values to csv
	println("Number of MH samples: $nMH");
	println("Number of Gibbs samples: $nGibbs");
	CSV.write("GibbsMCMC.csv", [DataFrame(:prm=>vkeys) DataFrame(SMP)], writeheader=false);
	
	# Save rejection statistics to csv
	dftemp = DataFrame(:rjtrej=>rjtcnt[:rjt],:rjttot=>rjtcnt[:tot],
			   :mhrej=>mhcnt[:rjt],:mhtot=>mhcnt[:tot]);
	CSV.write("RejStats.csv",dftemp);

	# Save random number generator state
	mysaverng(rng);


	return SMP
end
