using Documenter, Synapse, Setfield

const makepdf = false

if makepdf
	using DocumenterLaTeX
end

makedocs(doctest = false,
	sitename = "Model of excitatory synapse in Julia",
	format = makepdf ? DocumenterLaTeX.LaTeX(platform = "none") : Documenter.HTML(collapselevel = 1),
	pages = Any[
		"Home" => "index.md",
		"Simple example" => "tutorials.md",
		"Library" => "library.md"
	]
	)
