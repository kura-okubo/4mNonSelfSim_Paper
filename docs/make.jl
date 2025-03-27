using Documenter

push!(LOAD_PATH,"../src/")

makedocs(
	# format = Documenter.HTML(prettyurls = false), # turn on when locally compiling
    sitename = "4mNonSelfSim Paper",
    pages = [
	"Index" => "index.md",
    "Recipe of figures" => "plot_figures_recipe.md"
    ])

deploydocs(
    repo = "github.com/kura-okubo/4mNonSelfSim_Paper.jl.git",
)
