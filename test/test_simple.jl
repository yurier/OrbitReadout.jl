using Logging: global_logger

## ORBITS TO GENERATE THE REGIONS
pyplot()
	plot(layout = 4,windowsize=(700,300))
	function dynamics!(du,u,p,t)
		du[1] = -.00015*(u[1])*p[1]
		du[2] = -.001*(u[2])*p[1]
	end

	list=["A"; "B";  "C" ; "D"]
	list_weight=[0.1; -0.4; 0.75; 0.05]
	list_class=["N"; "D"; "P"; "N"]
	data=Vector{Array{Any,1}}[]
	pulses = 10
	after_pulses=1.2e4
	ipi=100.

for k in 1:4
	for j in 1:4
		u0 = [0,0]
		tspan = (0.0,ipi)
		prob = OrbitReadout.ODEProblem(dynamics!,u0,tspan,(j))
		sol = OrbitReadout.solve(prob)

		test=DataFrame(time=sol.t,
		          var1=sol[1,:],
		          var2=sol[2,:])

		for i in 1:pulses
		    u0=[sol[1,end] + j + 1.5*rand(), sol[2,end]+ j + 1.5*rand()]
		    tspan = (sol.t[end], sol.t[end] + ifelse(i==pulses,after_pulses, ipi))
		    prob = OrbitReadout.ODEProblem(dynamics!,u0,tspan,(j))
		    sol = OrbitReadout.solve(prob)
		    append!(test,DataFrame(time=sol.t,
		              var1=sol[1,:],
		              var2=sol[2,:]))
		end
			plot!(test[:,:time]/1000,test[:,:var1],subplot=1,label="",xlabel="time (s)",ylabel="var 1")
			plot!(test[:,:time]/1000,test[:,:var2],subplot=2,label="",xlabel="time (s)",ylabel="var 2")
			plot!(test[:,:var1],test[:,:var2],subplot=3,label="",xlabel="var2",ylabel="var1") |>display
			global data=[data;[[test[:,:var1],test[:,:var2],test[:,:time]/1000,list_weight[j],list_class[j]]]]
	end
end

## Create a mesh of polygones given the max and min value among all trajectories

space=15 #resolution
	HULL, xmin, xmax, ymin, ymax = OrbitReadout.polytopes_create(data,space; plot_ = true)
	dynmcl_pttrn = OrbitReadout.Dynamical_Pattern(nhull = nhull=size(HULL,1))


## DYNAMICAL 2D SYSTEM - Time spent by each orbit in each polygones
function jointt!(du,u,p,t)
	for j in 1:size(data,1)
		for i in 1:(p[1].nhull)
			du[i+p[1].nhull*(j-1)] 	= [p[3][j][1][OrbitReadout.last(p[3][j][3],t)],p[3][j][2][OrbitReadout.last(p[3][j][3],t)]] ∈ OrbitReadout.VPolygon(p[2][i])
		end
	end
end


u0 = zeros(size(data,1)*space*space)
    tspan = (0.0,after_pulses)
    prob = ODEProblem(jointt!,u0,tspan,(dynmcl_pttrn,HULL,data))
    sol = solve(prob,reltol=1e-5,saveat=0.75)

## Plot and get matrix of time spent plotted for every orbit (regardless of class)
# reduce plot density (To do)
OrbitReadout.time_spent_plot(data,space,nhull,sol; plot_ = false)


##
#average of matrices per class
for degree in collect(.1: .1 : .7) #degree from 0.1 to 1 which is the ammount of interesection we would like to mantain
	neighbour_tiles = 2 #0(min) to 7(max), number of neighbours of one tile to compose the polygonal region.
	Plots.plot(layout=(1,3),windowsize=(1.8*300,.7*200),grid=:none,border=:none)
	### Average time spent per class in the mesh
	classes_matrices, classes = OrbitReadout.average_classes(data,space,sol,nhull)
	### 'spent_time_positions' Mesh index of the matrix in which the classes spent most of its time; 'intersection_list' list of mesh index in which there is more than one class
	spent_time_positions, intersection_list = OrbitReadout.coincidence_detector(data,space,classes_matrices,degree)

	list_HULLs = Vector{Int64}[]
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
		A=[HULL[i] for i in list]
		Ax=[A[i][j][1] for i in 1:size(A,1) for j in 1:size(A[1],1)]
		Ay=[A[i][j][2] for i in 1:size(A,1) for j in 1:size(A[1],1)]
		A=[[Ax[i],Ay[i]] for i in 1:size(Ay,1)]
		polytopes=push!(polytopes,A)
		list_HULLs=push!(list_HULLs,list)
	end

	colls=	range(HSV(140,1,1), stop=HSV(360,1,1), length=3)
	for c in 1:length(classes)
		Plots.plot!( OrbitReadout.VPolygon(polytopes[c]), alpha=.1,label="",color=colls[c],subplot=2)
		Plots.plot!( OrbitReadout.VPolygon(polytopes[c]), alpha=.1,label="",color=colls[c],subplot=3)
		for i in list_HULLs[c]
			if classes_matrices[c][i]>0
				Plots.plot!( OrbitReadout.VPolygon(HULL[i]), alpha=classes_matrices[c][i]/maximum(classes_matrices[c][list_HULLs[c]]),label="",color=colls[c],subplot=2)
			end
		end
	end

	## Plot the training test against the regions
	for j in 1:length(classes)
		for i in 1:size(data,1)
			if data[i][5] == classes[j]
				plot!(data[i][1],data[i][2],color=colls[j],xlabel="\$x_1\$",ylabel="\$x_2\$",w=.5,subplot=3,label="" )
			end
		end
	end
	plot!(xlim=[xmin,xmax],ylim=[ymin,ymax]) |>display
end

## Test it
