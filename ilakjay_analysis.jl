### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ d71ed86c-9380-11eb-03b8-f9265659dd8a
using CSV

# ╔═╡ a56feb6a-936e-11eb-3f23-797dddf9cc7a
md"""
#### ilakjay.jl analysis
run the first function in this script to replicate main results of the model
"""

# ╔═╡ 6ca6445c-9371-11eb-033c-19ca15349274
module ilakjay
	include("ilakjay.jl")
end

# ╔═╡ 4261fc12-9372-11eb-16d3-9f4a6b5770a7
"""
replicate some main results from the model (assortativity, mean degree, clustering)

running the function will return a DataFrame with the data, as well as save it as a csv file.
"""
function replicate()
	ilakjay.replicate()
end

# ╔═╡ e80851e0-937a-11eb-2c6d-47c787ba1af1
md"""
additionally, if you want to see the effects of varying $p_r$, you can generate and use a dataset. I have used the same number of ticks (2000) and repetitions (20) as before, as well as the same range of $p_n$ = $(0.1, 0.9)$ and 100 nodes, but this time $p_r \in {0.001, 0.002, 0.005, 0.01}$.
"""

# ╔═╡ ca5e6956-937d-11eb-228d-73a03f280142
"""
more complete data collection than `replicate()`: `p_r` is not constant anymore, but takes longer to run. 
"""
function run_varying_pr()
	params = Dict(
		:p_n => collect(0.1:0.1:0.9), 
		:p_r => [0.001, 0.002, 0.005, 0.01, 0.015, 0.02],
		:num_nodes => 100,
		)
	
	mdata = [
		:num_nodes, :p_n, :p_r, 
		:expected_mean_degree, 
		:mean_degree, :assort
		]
	
	df = ilakjay.paramscan_social_inheritance(params, mdata, 2000, 500)
	CSV.write("varying_pr.csv", df)
	return df
end

# ╔═╡ Cell order:
# ╠═a56feb6a-936e-11eb-3f23-797dddf9cc7a
# ╠═6ca6445c-9371-11eb-033c-19ca15349274
# ╠═d71ed86c-9380-11eb-03b8-f9265659dd8a
# ╠═4261fc12-9372-11eb-16d3-9f4a6b5770a7
# ╟─e80851e0-937a-11eb-2c6d-47c787ba1af1
# ╠═ca5e6956-937d-11eb-228d-73a03f280142
