

############ to read data inside diffeq
function last(full_previous::Vector{Float64},current::Float64)
	return findall(x->x <=current, full_previous)[end]
end

function minmax_x_y(data)
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
	return xmax,xmin,ymax,ymin
end

############ creating polytopes
function polytopes_create(data,space; plot_ = false)
	## find limits of polytope grid
	xmax,xmin,ymax,ymin = minmax_x_y(data)

	pointsy=collect(ymin : (ymax-ymin)/space: ymax)
	pointsx=collect(xmin : (xmax-xmin)/space: xmax)

	mesh=Array{Array{Float64,1},1}[]
	if plot_ == true
		c1 = colorant"red"
		c2 = colorant"blue"
		Plots.plot(layout=(1,3),windowsize=(2.5*300,200),grid=:none,border=:none)  |> display
	end
	k=0
	for j in -1:(space-2)
		for i in 0:(space-1)
			#@show j,[space-1-j,space-j]
			mesh_entry=[[i,j] for i in pointsx[[1+i,2+i]] for j in pointsy[[space-1-j,space-j]]]
			k = k + 1
			if plot_ == true
				Plots.plot!( VPolygon(mesh_entry), alpha=.25,label="",color=weighted_color_mean((k)/(space*space), c1, c2),subplot=1)
			end
			append!(mesh,[mesh_entry])
		end
	end
	if plot_ == true
		Plots.plot!() |>display
	end
	return mesh, xmin, xmax, ymin, ymax
end


############ average matrix of time spent per class
function average_classes(data,space,sol)
	#time spent matrix (the organisation of the matrix is the polytope list is to facilitate the plot)
	mesh_size = space*space
	tspent_matrix_all=[zeros(space,space) for i in 1:size(data,1)]
	for j in 0:(size(data,1)-1)
		for i in 1:(space*space)
			tspent_matrix_all[j+1][i] = tspent_matrix_all[j+1][i] + sol.u[end][i+mesh_size*(j)]
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


## Plot and get matrix of time spent plotted for every orbit (regardless of class)
############ plot time spent for all orbits
# Example: OrbitReadout.time_spent_plot(data,space,mesh_size,sol; plot_ = false)

function time_spent_plot(data,space,sol; plot_ = false)
		if plot_ == true
			c1 = colorant"red"
			c2 = colorant"blue"
			p=Plots.plot(layout=(space,space),windowsize=(6*300,6*200),axis=:none,grid=:none,border=:none,legend=:none)
		end
		mesh_size = space * space
		score=[zeros(space,space) for i in 1:size(data,1)]
		for j in 0:(size(data,1)-1)
				for i in 1:(space*space)
					@show i
					if plot_ == true
						p=Plots.plot!(p,sol,vars=(i+mesh_size*(j)),subplot=(i),label="",xlabel="",w=3,color=weighted_color_mean((i/(space*space)), c1, c2))
					end
					score[j+1][i] = score[j+1][i] + sol.u[end][i+mesh_size*(j)]
				end
		end
		if plot_ == true
			plot!() |>display
		end
		return score
end



## DYNAMICAL 2D SYSTEM - Time spent by each orbit in each polygones
function time_spent_mesh!(du,u,p,t)
	for j in 1:size(p[3],1)
		aux=[p[3][j][1][OrbitReadout.last(p[3][j][3],t)],p[3][j][2][OrbitReadout.last(p[3][j][3],t)]]
		for i in 1:p[1]
			du[i+p[1]*(j-1)] 	= aux ∈ OrbitReadout.VPolygon(p[2][i])
		end
	end
end

## GENERATE GEOMETRICAL READOUT

function generate_geometrical_readout(degree, neighbour_tiles, data, space, sol, mesh; plot_ = false)
	#'degree' - collect(.1: .1 : .7) #degree from 0.1 to 1 which is the ammount of interesection we would like to mantain
	#'neighbour_tiles' - 0(min) to 7(max), number of neighbours of one tile to compose the polygonal region.

	### Average time spent per class in the mesh
	classes_matrices, classes = OrbitReadout.average_classes(data,space,sol)
	### 'spent_time_positions' Mesh index of the matrix in which the classes spent most of its time; 'intersection_list' list of mesh index in which there is more than one class
	spent_time_positions, intersection_list = OrbitReadout.coincidence_detector(data,space,classes_matrices,degree)

	list_mesh = Vector{Int64}[]
	polytopes=Vector{Array{Float64,1}}[]

	for cls in 1:length(classes)
		list = copy(spent_time_positions[cls])
		fixed_list = copy(spent_time_positions[cls])

		#Is it an isolated region? Search around
		for i in sort(fixed_list)
			remove_isolate=0
			remove_isolate = remove_isolate + 1*(i+1 ∈ fixed_list)
			remove_isolate = remove_isolate + 1*(i-1 ∈ fixed_list)
			remove_isolate = remove_isolate + 1*(i+space ∈ fixed_list)
			remove_isolate = remove_isolate + 1*(i-space ∈ fixed_list)
			remove_isolate = remove_isolate + 1*(i+1+space ∈ fixed_list)
			remove_isolate = remove_isolate + 1*(i+1-space ∈ fixed_list)
			remove_isolate = remove_isolate + 1*(i-1+space ∈ fixed_list)
			remove_isolate = remove_isolate + 1*(i-1-space ∈ fixed_list)
			if remove_isolate < neighbour_tiles
				filter!(e->e ≠ i ,list) #remove from the list if there are not enough tiles around it
			end
		end
		A=[mesh[i] for i in list]
		Ax=[A[i][j][1] for i in 1:size(A,1) for j in 1:size(A[1],1)]
		Ay=[A[i][j][2] for i in 1:size(A,1) for j in 1:size(A[1],1)]
		A=[[Ax[i],Ay[i]] for i in 1:size(Ay,1)]
		polytopes=push!(polytopes,A)
		list_mesh=push!(list_mesh,list)
	end

	if plot_ == true
		Plots.plot(layout=(1,3),windowsize=(1.8*300,.7*200),grid=:none,border=:none)
		colls=	range(HSV(140,1,1), stop=HSV(360,1,1), length=3)
		for c in 1:length(classes)
			Plots.plot!( OrbitReadout.VPolygon(polytopes[c]), alpha=.1,label="",color=colls[c],subplot=2)
			Plots.plot!( OrbitReadout.VPolygon(polytopes[c]), alpha=.1,label="",color=colls[c],subplot=3)
			for i in list_mesh[c]
				if classes_matrices[c][i]>0
					Plots.plot!( OrbitReadout.VPolygon(mesh[i]), alpha=classes_matrices[c][i]/maximum(classes_matrices[c][list_mesh[c]]),label="",color=colls[c],subplot=1)
					Plots.plot!( OrbitReadout.VPolygon(mesh[i]), alpha=classes_matrices[c][i]/maximum(classes_matrices[c][list_mesh[c]]),label="",color=colls[c],subplot=2)
				end
			end
		end

		## Plot the training test against the regions
		for j in 1:length(classes)
			for i in 1:size(data,1)
				if data[i][5] == classes[j]
					plot!(data[i][1],data[i][2],color=colls[j],xlabel="\$x_1\$",ylabel="\$x_2\$",w=.5,subplot=2,label="" )
					plot!(data[i][1],data[i][2],color=colls[j],xlabel="\$x_1\$",ylabel="\$x_2\$",w=.5,subplot=3,label="" )
				end
			end
		end
		xmax,xmin,ymax,ymin = minmax_x_y(data)
		plot!(xlim=[xmin,xmax],ylim=[ymin,ymax]) |>display
	end
	return polytopes
end
