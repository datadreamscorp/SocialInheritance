### A Pluto.jl notebook ###
# v0.14.1

using Markdown
using InteractiveUtils

# ╔═╡ d71ed86c-9380-11eb-03b8-f9265659dd8a
using CSV, Pipe, Statistics, DataFrames

# ╔═╡ a56feb6a-936e-11eb-3f23-797dddf9cc7a
md"""
#### ilakjay.jl analysis
run the functions in this script to replicate main results of the model
"""

# ╔═╡ 6ca6445c-9371-11eb-033c-19ca15349274
module ilakcay
	include("ilakcay.jl")
end

# ╔═╡ 87b0261a-23db-4fa9-b6a5-525dac3bc5bf
"""
replicate assortativity results from the model
"""
function replicate_assortativity()
	
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
	
	rep_df_raw = ilakcay.paramscan_social_inheritance( #compute it, sam
		parameters, 
		mdata, 
		num_ticks, 
		num_reps
	)
	
	funcs = [ #functions to aggregate our data
		mean, std
	]
	
	agg_params = [ #parameters to aggregate over
		:assort, :assort
	]
	
	rep_df = @pipe groupby(rep_df_raw, [:p_n, :p_r]) |> 
			 combine(_, agg_params .=> funcs) #split-apply-combine
	#and we prepare some more columns that will be useful for plotting
	
	d = sqrt(num_reps)
	
	rep_df[:,"assort_min"] = rep_df[:,"assort_mean"] - 
							 rep_df[:,"assort_std"] ./ d
	
	rep_df[:,"assort_max"] = rep_df[:,"assort_mean"] + 
							 rep_df[:,"assort_std"] ./ d
	
	CSV.write("rep_assortativity.csv", rep_df) #save to csv
	
	return rep_df
	
end

# ╔═╡ ca5e6956-937d-11eb-228d-73a03f280142
"""
replicate degree and clustering results from the model 
"""
function replicate_degree_clustering()
	params = Dict(
		:p_n => collect(0:0.05:1), 
		:p_r => collect(0:0.05:1),
		:num_nodes => 100,
		)
	
	mdata = [
		:num_nodes, :p_n, :p_r, 
		:expected_mean_degree, 
		:mean_degree, :assort,
		:clustering
		]
	
	num_reps = 50
	
	rep_df_raw = ilakcay.paramscan_social_inheritance(params, mdata, 2000, num_reps)
	
	funcs = [ #functions to aggregate our data
		mean, std,
		mean, std,
		mean
	]
	
	agg_params = [ #parameters to aggregate over
		:clustering, :clustering, 
		:mean_degree, :mean_degree,
		:expected_mean_degree
	]
	
	rep_df = @pipe groupby(rep_df_raw, [:p_n, :p_r]) |> 
			 combine(_, agg_params .=> funcs) #split-apply-combine
	#and we prepare some more columns that will be useful for plotting
	
	rep_df[:, :clust_min] = rep_df[:, :clustering_mean] - 
							rep_df[:, :clustering_std]
	
	rep_df[:, :clust_max] = rep_df[:,:clustering_mean] + 
							rep_df[:,:clustering_std]
	
	rep_df[:, :degree_min] = rep_df[:, :mean_degree_mean] - 
						   	 rep_df[:, :mean_degree_std]
	
	rep_df[:, :degree_max] = rep_df[:, :mean_degree_mean] + 
						   	 rep_df[:, :mean_degree_std]
	
	CSV.write("rep_degree_clustering.csv", rep_df)
	return rep_df
end

# ╔═╡ 1064eb11-a8c9-4755-98fa-6f6aa3e22254
#replicate_degree_clustering()

# ╔═╡ b9785c91-25dd-4053-bd94-180bb50ac416
#replicate_assortativity()

# ╔═╡ Cell order:
# ╟─a56feb6a-936e-11eb-3f23-797dddf9cc7a
# ╠═6ca6445c-9371-11eb-033c-19ca15349274
# ╠═d71ed86c-9380-11eb-03b8-f9265659dd8a
# ╠═87b0261a-23db-4fa9-b6a5-525dac3bc5bf
# ╠═ca5e6956-937d-11eb-228d-73a03f280142
# ╠═1064eb11-a8c9-4755-98fa-6f6aa3e22254
# ╠═b9785c91-25dd-4053-bd94-180bb50ac416
