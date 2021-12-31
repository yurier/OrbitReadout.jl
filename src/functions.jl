

############ to read data inside diffeq
function last(full_previous,current)
	return findall(x->x <=current, full_previous)[end]
end


############ creating polytopes
function polytopes_create(data,space; plot_ = false)
	## find limits of polytope grid
	ymax=-Inf;xmax=-Inf;ymin=Inf;xmin=Inf
	for i in 1:size(data,1)
		Xmax=maximum(data[i][1])
		Ymax=maximum(data[i][2])
		Xmin=minimum(data[i][1])
		Ymin=minimum(data[i][2])
		if Xmax>xmax
			xmax = Xmax
		end
		if Ymax>ymax
			ymax = Ymax
		end
		if Xmin<xmin
			xmin = Xmin
		end
		if Ymin<ymin
			ymin = Ymin
		end
	end
	pointsy=collect(ymin : (ymax-ymin)/space: ymax)
	pointsx=collect(xmin : (xmax-xmin)/space: xmax)

	HULL=Array{Array{Float64,1},1}[]
	if plot_ == true
		c1 = colorant"red"
		c2 = colorant"blue"
		Plots.plot(layout=(1,3),windowsize=(2.5*300,200),grid=:none,border=:none)  |> display
	end
	k=0
	for j in -1:(space-2)
		for i in 0:(space-1)
			#@show j,[space-1-j,space-j]
			hull=[[i,j] for i in pointsx[[1+i,2+i]] for j in pointsy[[space-1-j,space-j]]]
			k = k + 1
			if plot_ == true
				Plots.plot!( VPolygon(hull), alpha=.25,label="",color=weighted_color_mean((k)/(space*space), c1, c2),subplot=1)
			end
			append!(HULL,[hull])
		end
		if plot_ == true
			Plots.plot!() |>display
		end
	end
	return HULL, xmin, xmax, ymin, ymax
end


############ average matrix of time spent per class
function average_classes(data,space,sol,nhull)
	#time spent matrix (the organisation of the matrix is the polytope list is to facilitate the plot)
	tspent_matrix_all=[zeros(space,space) for i in 1:size(data,1)]
	for j in 0:(size(data,1)-1)
		for i in 1:(space*space)
			tspent_matrix_all[j+1][i] = tspent_matrix_all[j+1][i] + sol.u[end][i+nhull*(j)]
		end
	end

	#sum up all matrices of a given the class
	classes=unique([data[i][5] for i in 1:size(data,1)])
	classes_matrices=[zeros(space,space) for i in 1:length(classes)]
	for j in 1:length(classes)
		n_class = 0
		aux = zeros(space,space)
		for i in 1:size(data,1)
			if data[i][5] == classes[j]
				n_class=n_class + 1
				aux = aux .+ tspent_matrix_all[i]
			end
		end
		#average of the class time spent
		classes_matrices[j]=aux/n_class
	end
	return classes_matrices , classes
end

############ find polytopes with intersection and usage of the most "busy polytopes"
function coincidence_detector(data,space,classes_matrices,degree)
	classes=unique([data[i][5] for i in 1:size(data,1)])
	# take the Int(round(space*space/degree)) first ranked positions which the class spend more time
	# for j in 1:length(classes)
	# 	push!(spent_time_positions,sortperm(reshape(classes_matrices[j],space*space),rev=true)) #[1:Int(round(space*space/degree))]
	# end

	#nonzero indexes within the time spent matrix per class
	nonzero_list = Vector{Int64}[]
	for j in 1:length(classes)
		aux=Int64[]
		for i in 1:(space*space)
			if reshape(classes_matrices[j],space*space)[i] > 0
				append!(aux,i)
			end
		end
		push!(nonzero_list,aux)
		#push!(spent_time_positions,nonzero_list[1:Int(round(length(nonzero_list)/degree))])
	end


	x=nonzero_list
	vectorize=[x[j] for x in x for j in eachindex(x)]

	d=[[i,count(x->x==i,vectorize)] for i in 1:(space*space)] #how many times a polytope is found within classes
	intersection_list=Int64[]
		for i in 1:size(d,1)  #list the polytopes which are not exclusive to a class
			if d[i][2] > 1#(rand(Bool)+1)
				append!(intersection_list,d[i][1])
			end
		end

	spent_time_positions = Vector{Int64}[]#zeros(Int64,Int(round(space*space/degree)),length(classes))
	for j in 1:length(classes)
		ranked_time_spent=sortperm(reshape(classes_matrices[j],space*space),rev=true)
		nointer=filter!(e->e∉ intersection_list,ranked_time_spent[1:Int(round((space*space)*degree))])
		push!(spent_time_positions, filter!(e->e ∈ nonzero_list[j],nointer))
	end


	return spent_time_positions, intersection_list
end

############ plot time spent for all orbits
function time_spent_plot(data,space,nhull,sol; plot_ = false)
		if plot_ == true
			c1 = colorant"red"
			c2 = colorant"blue"
			p=Plots.plot(layout=(space,space),windowsize=(6*300,6*200),axis=:none,grid=:none,border=:none,legend=:none)
		end
		score=[zeros(space,space) for i in 1:size(data,1)]
		for j in 0:(size(data,1)-1)
				for i in 1:(space*space)
					@show i
					if plot_ == true
						p=Plots.plot!(p,sol,vars=(i+nhull*(j)),subplot=(i),label="",xlabel="",w=3,color=weighted_color_mean((i/(space*space)), c1, c2))
					end
					score[j+1][i] = score[j+1][i] + sol.u[end][i+nhull*(j)]
				end
		end
		if plot_ == true
			plot!() |>display
		end
		return score
end
