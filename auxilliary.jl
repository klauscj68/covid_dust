using Random,CSV,DataFrames
#%% mysaverng
"""
Given a MersenneTwister rng, save its internal states to csv files which may be reloaded for
restarting future runs.
"""
function mysaverng(rng::MersenneTwister)
	# Loop over fieldnames of rng and save to their own csv's
	#  One of the fields was a structure of its own and needed slightly different handling
	for key in fieldnames(typeof(rng))
		if key != :state
			df = DataFrame(string(key)=>getfield(rng,key));
		else
			ram = getfield(rng,key);
			df = DataFrame(string(key)=>ram.val);
		end
		CSV.write("RNG"*string(key)*".csv",df);
	end
	return true
end

#%% myloadrng
"""
Given a MersenneTwister saved into csv's like my mysaverng function, restore a rng with these settings.
fname is the prefix that the CSV files begin with.
"""
function myloadrng(fname::String="RNG")
	fields = ["seed","state","vals","ints","idxF","idxI","adv","adv_jump","adv_vals","adv_ints"];
	
	# Seed
	DF = CSV.read(fname*"seed.csv",DataFrame);
	myseed = convert(Vector{UInt32},DF[!,1]);

	# State
	DF = CSV.read(fname*"state.csv",DataFrame);
	mystate = Random.DSFMT.DSFMT_state(convert(Vector{Int32},DF[!,1]));

	# vals
	DF = CSV.read(fname*"vals.csv",DataFrame);
	myvals = convert(Vector{Float64},DF[!,1]);

	# ints
	DF = CSV.read(fname*"ints.csv",DataFrame);
	myints = convert(Vector{UInt128},DF[!,1]);

	# idxF
	DF = CSV.read(fname*"idxF.csv",DataFrame);
	myidxF = convert(Int64,DF[!,1][1]);

	# idxI
        DF = CSV.read(fname*"idxI.csv",DataFrame);
	myidxI = convert(Int64,DF[!,1][1]);

	# adv
	DF = CSV.read(fname*"adv.csv",DataFrame);
	myadv = convert(Int64,DF[!,1][1]);

	# adv_jump
	DF = CSV.read(fname*"adv_jump.csv",DataFrame);
	myadv_jump = convert(BigInt,DF[!,1][1]);

	# adv_vals
	DF = CSV.read(fname*"adv_vals.csv",DataFrame);
	myadv_vals = convert(Int64,DF[!,1][1]);

	# adv_ints
	DF = CSV.read(fname*"adv_ints.csv",DataFrame,type=Int64);
	myadv_ints = convert(Int64,DF[!,1][1]);

	return MersenneTwister(myseed,mystate,myvals,myints,myidxF,myidxI,myadv,myadv_jump,myadv_vals,myadv_ints)
end
