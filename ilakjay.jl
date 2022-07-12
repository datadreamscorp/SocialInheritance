### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ f7d02abc-8f79-11eb-03c9-310218be3fc4
begin
	using Pkg
	Pkg.activate( mktempdir() )						#activate temporary directory
	Pkg.add(["DrWatson", "Agents", "LightGraphs",		#add packages
			 "Distributions", "Random", "Pipe",
			 "Gadfly", "Colors", "GraphPlot"])

	using DrWatson, Agents, LightGraphs,
	Distributions, Random, Statistics, Pipe, 
	Gadfly, Colors, GraphPlot, DataFrames,
	CSV
end

# ╔═╡ c28efb9a-910d-11eb-18a5-61936d98cd1b
mutable struct Recruit <: AbstractAgent
		id::Int #agent id
		pos::Int #agent position
		pers::Float64 #generic trait
end

# ╔═╡ c952630e-910d-11eb-1c0b-fb550ac6db34
#initializing function
		
function social_inheritance(;
		p_b::Float64 = 1.0, #probability of linking with recruiter
		p_n::Float64 = 0.2, #probability of connection inheritance
		p_r::Float64 = 0.01, #probability of random connection at birth
		num_nodes::Int64 = 200, #number of nodes/agents
		init_dens::Float64 = 0.05, #density of initial erdos-renyi network
		seed::Int64 = 125, #seed value
		var::Float64 = 0.01, #variance of normal distribution for trait inheritance
		mu::Float64 = 0.05, #drastic mutation rate
		)
	
	#we calculate the expected mean degree as derived in the paper:
	
	expected_mean_degree = ((num_nodes - 1) * (1 + (num_nodes - 2) * p_r)) / 
						   (num_nodes - 1 - (num_nodes - 2) * (p_n - p_r))
	
	#define the model
	
	model = ABM(
	        Recruit,
			GraphSpace( erdos_renyi(num_nodes, init_dens) ); #we start with a
	        properties = Dict(								 #random graph
				:num_nodes => num_nodes,
				:p_b => p_b,
				:p_n => p_n,
				:p_r => p_r,
				:var => var,
				:mu => mu,
				:expected_mean_degree => expected_mean_degree,
				:mean_degree => 0.0,
				:clustering => 0.0,
				:assort => 0.0,
				:tick => 0,
	        	),
			rng = RandomDevice()
			)
			
	for position in 1:num_nodes #add the agents with random pers
		add_agent!(position, model, rand(model.rng))
	end
			
	return model
end

# ╔═╡ ec3071ae-910d-11eb-37fb-c9d49cc4984e
#agent birth function
		
function birth!(model, old_ids)
		
	recruiter = random_agent(model) #choose a random agent as the recruiter
	
	nbs = nearby_ids(recruiter.pos, model, 1) |> collect #degree 1 neighbors
															 #of recruiter
	
	new_node = add_node!(model) #add a new node and 																#save its position to a variable
	
	#if a large mutation happens, we choose a trait at random
	#else, we inherit a trait from recruiter:
	
	new_pers = (rand(model.rng) < model.mu) ? rand(model.rng) : 
	@pipe recruiter.pers |> Normal(_, sqrt(model.var)) |> rand(model.rng, _)
	
	#we add the agent with its position and trait value
	
	add_agent!( new_node, model, clamp(new_pers, 0, 1) )
	
	if rand(model.rng) < model.p_b #if I connect with recruiter
		add_edge!(model, new_node, recruiter.pos) #add edge with recruiter
		if length(nbs) != 0
			for nb in nbs
				if rand(model.rng) < model.p_n          #add the connections
					add_edge!(model, new_node, model[nb].pos) #we inherit
				end
			end
		end
	end
	
	non_nbs = setdiff(old_ids, nbs) #non_neighbors
	for nnb in non_nbs
		if rand(model.rng) < model.p_r 			#add the connections
			add_edge!(model, new_node, model[nnb].pos)	#we form at random
		end
	end
	
end

# ╔═╡ 27079ea6-910e-11eb-1c50-09381afd4195
#agent death function
		
function death!(model, old_ids)
	dies = rand(model.rng, old_ids) #id at random for agent that dies
	death_position = model[dies].pos #position of dying agent
	last_position = nv(model.space.graph) #position of "last" node in the graph
	last_id = collect(
			filter(id -> model[id].pos == last_position, allids(model))
			)[1] #collect id of model at last node
	rem_node!(model, death_position) #F	
	model[last_id].pos = death_position #manually change the position
										#of agent at last node
end

# ╔═╡ 3257b250-910e-11eb-209f-696a113c22f4
#function to calculate assortativity
	
function assort!(model)
	c = 0 #normalizing coefficient
	non_norm = 0 #non-normalized assortativity
	num_nodes = nv(model.space.graph)
	num_edges = ne(model.space.graph)
	adj = adjacency_matrix(model.space.graph)
	
	#we perform this calculation for every pair of nodes in our network,
	#summing up our results as we calculate them
	
	i = 1
	while i <= num_nodes
		j = 1
		while j <= num_nodes
			agt_i = (filter(id->model[id].pos==i, allids(model)) |> collect)[1] 
			#id of i
			deg_i = degree(model.space.graph, i) #degree of i
			trait_i = model[agt_i].pers #trait value of i
			
			agt_j = (filter(id->model[id].pos==j, allids(model)) |> collect)[1] 
			#id of j
			deg_j = degree(model.space.graph, j) #degree of j
			trait_j = model[agt_j].pers #trait value of j
			
			c += ( (i==j ? 1 : 0)*deg_i - deg_i*deg_j/(2*num_edges) ) * 
				 (trait_i*trait_j) #sd(i)*sd(j)
			
			non_norm += ( adj[i,j] - deg_i*deg_j/(2*num_edges) )*trait_i*trait_j
			#cov(i, j)
			j += 1
		end
		i += 1
	end
	
	model.assort = non_norm/c
	
end

# ╔═╡ 3acd6e0c-8ff1-11eb-268c-277fb0d8cc0a
#sub-step function, one life cycle

function lifecycle_step!(model)
		
	old_ids = allids(model) |> collect
	
	birth!(model, old_ids)
	death!(model, old_ids)
	
end

# ╔═╡ 1f2d34da-91d7-11eb-3566-3d8ea46c2564
#complex step

function step_social_inheritance!(model, n)
	t = 1
	while t <= n
		step!(model, dummystep, lifecycle_step!, 1) #every tick, we step model 1ce
		if t == n #if at end of process, calculate network properties
			assort!(model)
			model.mean_degree = degree(model.space.graph) |> collect |> mean
			model.clustering = local_clustering_coefficient(model.space.graph) |>
							   collect |> mean
		end
		model.tick = t #tiki
		t += 1 #toki
	end
end

# ╔═╡ ef75bb96-9282-11eb-32a6-8127bf200ed3
#parameter scanning function

function paramscan_social_inheritance(parameters, mdata, n, reps)
	df_list = [] #we will populate this list with dataframes
	combinations = dict_list(parameters) #list of key-value pairs; 
										 #all param combinations
	for combination in combinations
		i = 1
		while i <= reps #i batches of each parameter combination
			model = social_inheritance(; combination...) #init model
			df_model = init_model_dataframe(model, mdata) #init model dataframe
			step_social_inheritance!(model, n) #step the model
			collect_model_data!(df_model, model, mdata, n) #collect the data
			push!(df_list, df_model) #push the dataframe into our list
			i += 1
		end
	end
	return vcat(df_list...) #finally, vertically concatenate everything
end

# ╔═╡ 4614c074-9264-11eb-36b9-97b11899f3c2
#function that plots the model at a particular time step

function network_plots(model)
	assort = [a.pers for a in allagents(model)] |> collect #collect trait values
	assort_plot = plot(									   #and save plot
		x=assort,
		Geom.histogram,
		Guide.xlabel("trait value"),
		Theme(default_color=RGB(0.7,0.2,0.4))
	)
	#network plot
	net_plot = gplot(model.space.graph, nodefillc=[RGB(0.5,i,i) for i in assort])
	
	deg = degree(model.space.graph) |> collect #collect degrees and plot
	degree_plot= plot(
		x=deg,
		Geom.histogram,
		Guide.xlabel("degree"),
		Theme(default_color=RGB(0.5,0.7,0.1))	
	)
	
	clust = local_clustering_coefficient( model.space.graph, #collect clustering
											 1:nv(model.space.graph) ) |> collect
	clust_plot = plot( #plot clustering
		x=clust,
		Geom.histogram,
		Guide.xlabel("clustering coefficient")
	)
	
	return @pipe hstack(net_plot, assort_plot) |>
	vstack( _, hstack(degree_plot, clust_plot) )
end

# ╔═╡ b82183c2-9364-11eb-0cd9-770462a88c9b
function replicate()
	
	num_ticks = 2000
	num_reps = 500
	
	mdata = [ #model data to collect
		:num_nodes, :p_n, :p_r, 
		:expected_mean_degree, 
		:mean_degree, :assort
	]
	
	parameters = Dict( #parameter ranges we will use
		:p_n => collect(0.1:0.1:0.9),
		:p_r => 0.01,
		:num_nodes => 100
	)
	
	rep_df_raw = paramscan_social_inheritance( #compute it, sam
		parameters, 
		mdata, 
		num_ticks, 
		num_reps
	)
	
	funcs = [ #functions to aggregate our data
		mean, std,
		mean, std,
		mean
	]
	
	params = [ #parameters to aggregate over
		:assort, :assort, 
		:mean_degree, :mean_degree,
		:expected_mean_degree
	]
	
	rep_df = @pipe groupby(rep_df_raw, [:p_n, :p_r]) |> 
			 combine(_, params .=> funcs) #split-apply-combine
	#and we prepare some more columns that will be useful for plotting
	
	d = sqrt(num_reps)
	
	rep_df[:,"assort_min"] = rep_df[:,"assort_mean"] - 
							 rep_df[:,"assort_std"] ./ d
	
	rep_df[:,"assort_max"] = rep_df[:,"assort_mean"] + 
							 rep_df[:,"assort_std"] ./ d
	
	rep_df[:,"mean_min"] = rep_df[:,"mean_degree_mean"] + 
						   rep_df[:,"mean_degree_std"] ./ d
	
	rep_df[:,"mean_max"] = rep_df[:,"mean_degree_mean"] + 
						   rep_df[:,"mean_degree_std"] ./ d
	
	CSV.write("rep_df.csv", rep_df) #save to csv
	
	return rep_df
	
end

# ╔═╡ Cell order:
# ╠═f7d02abc-8f79-11eb-03c9-310218be3fc4
# ╠═c28efb9a-910d-11eb-18a5-61936d98cd1b
# ╠═c952630e-910d-11eb-1c0b-fb550ac6db34
# ╠═ec3071ae-910d-11eb-37fb-c9d49cc4984e
# ╠═27079ea6-910e-11eb-1c50-09381afd4195
# ╠═3257b250-910e-11eb-209f-696a113c22f4
# ╠═3acd6e0c-8ff1-11eb-268c-277fb0d8cc0a
# ╠═1f2d34da-91d7-11eb-3566-3d8ea46c2564
# ╠═ef75bb96-9282-11eb-32a6-8127bf200ed3
# ╠═4614c074-9264-11eb-36b9-97b11899f3c2
# ╠═b82183c2-9364-11eb-0cd9-770462a88c9b
