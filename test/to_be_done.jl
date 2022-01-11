using OrbitReadout

##################################################
#VECTOR
#enzymes=[[CaMKII[i],CaN[i]] for i in 1:size(CaMKII,1)]

# list=["1Pre";
# 	"2Pre-10";
# 	"2Pre-50";
# 	"2Post-1Pre50";
# 	"2Post-1Pre20";
# 	"1Pre-2Post50";
# 	"1Pre-2Post10";
# 	"1Pre-1Post10"]

list=["1Pre1Post10"; "2Pre50";  "1Pre2Post10" ; "2Pre10"]


# list_weight=[0.;
# 	0.05;
# 	-0.4;
# 	0.;
# 	0.25;
# 	0.5;
# 	0.75;
# 	0.1]

list_weight=[0.1; -0.4; 0.75; 0.05]
#
# list_class=["N";
# 	"N";
# 	"D";
# 	"N";
# 	"P";
# 	"P";
# 	"P";
# 	"N"]

list_class=["N"; "D"; "P"; "N"]

data=Vector{Array{Float64,1}}[]
	for j in 1:size(list,1)
		for i in 1:2
		protocol=list[j]
		tt=OrbitReadout.load("test/data_test/$(i)_$(protocol).jld2","time")
			g=OrbitReadout.load("test/data_test/$(i)_$(protocol).jld2","CaN")
			h=OrbitReadout.load("test/data_test/$(i)_$(protocol).jld2","CaMKII")
			# idx=[]
			# for i in collect(0. : tt[end]/500: tt[end])
			# 	append!(idx, findall(x->x <=i, tt)[end])
			# 	@show i/ tt[end]
			# end
			data=[data;[[g,h,tt,list_weight[j],list_class[j]]]]
		end
	end

## DYNAMICAL SYSTEM
# nhull=size(HULL,1)

space=7 #resolution
HULL, xmin, xmax, ymin, ymax = OrbitReadout.polytopes_create(data,space; plot_ = true)
	dynmcl_pttrn = OrbitReadout.Dynamical_Pattern(
				   nhull = nhull=size(HULL,1))

# using LaTeXStrings, ColorSchemes
	# i=1;plot!(data[i][1],data[i][2],color=get(colorschemes[:viridis],.5),xlabel="\$x_1\$",ylabel="\$x_2\$",label="N",w=2,subplot=3)
	# i=2;plot!(data[i][1],data[i][2],color=get(colorschemes[:viridis],.5),xlabel="\$x_1\$",ylabel="\$x_2\$",label="N",w=2,subplot=3)
	# i=3;plot!(data[i][1],data[i][2],color=get(colorschemes[:viridis],0.),xlabel="\$x_1\$",ylabel="\$x_2\$",label="D",w=2,subplot=3)
	# i=4;plot!(data[i][1],data[i][2],color=get(colorschemes[:viridis],0.),xlabel="\$x_1\$",ylabel="\$x_2\$",label="D",w=2,subplot=3)
	# i=5;plot!(data[i][1],data[i][2],color=get(colorschemes[:viridis],1.),xlabel="\$x_1\$",ylabel="\$x_2\$",label="P",w=2,subplot=3)
	# i=6;plot!(data[i][1],data[i][2],color=get(colorschemes[:viridis],1.),xlabel="\$x_1\$",ylabel="\$x_2\$",label="P",w=2,subplot=3)
	# i=7;plot!(data[i][1],data[i][2],color=get(colorschemes[:viridis],.5),xlabel="\$x_1\$",ylabel="\$x_2\$",label="N",w=2,subplot=3)
	# i=8;plot!(data[i][1],data[i][2],color=get(colorschemes[:viridis],.5),xlabel="\$x_1\$",ylabel="\$x_2\$",label="N",w=2,subplot=3)
# 	#Plots.savefig("example5.svg")


# DYNAMICAL 2D SYSTEM
function jointt!(du,u,p,t)
	for j in 1:size(data,1)
		for i in 1:(p[1].nhull)
			#du[i+p[1].nhull*(j-1)] 	= [data[j][1][last(data[j][3],t)],data[j][2][last(data[j][3],t)]] ∈ VPolygon(p[2][i])
			du[i+p[1].nhull*(j-1)] 	= [p[3][j][1][OrbitReadout.last(p[3][j][3],t)],p[3][j][2][OrbitReadout.last(p[3][j][3],t)]] ∈ OrbitReadout.VPolygon(p[2][i])
		end
	end
end

#SIMULATION
u0 = zeros(size(data,1)*space*space)
    tend= 4e5
    tspan = (0.0,tend)
    prob = OrbitReadout.ODEProblem(jointt!,u0,tspan,(dynmcl_pttrn,HULL,data))
    @time  sol = OrbitReadout.solve(prob,reltol=1e-4,saveat=0.85)

#plot and get matrix plotted
time_spent_plot(data,space,nhull; plot_ = false)

#average of matrices per class
c1 = colorant"red"
	c2 = colorant"blue"
	for degree in 1
	@show space^2/degree
	Plots.plot(layout=(1,3),windowsize=(1.8*300,.7*200),grid=:none,border=:none)
	classes_matrices, classes = average_classes(data,space,sol,nhull)
	spent_time_positions, intersection_list = coincidence_detector(data,space,classes_matrices,degree)

	list_HULLs = Vector{Int64}[]
	polytopes=Vector{Array{Float64,1}}[]

	for cls in 1:3
		list = spent_time_positions[cls]
		#list = list[classes_matrices[cls][list] .> 0]

		for i in sort(list)
			remove_isolate=0
			remove_isolate = remove_isolate + 1*(i+1 ∈ list)
			remove_isolate = remove_isolate + 1*(i-1 ∈ list)
			remove_isolate = remove_isolate + 1*(i+space ∈ list)
			remove_isolate = remove_isolate + 1*(i-space ∈ list)
			remove_isolate = remove_isolate + 1*(i+1+space ∈ list)
			remove_isolate = remove_isolate + 1*(i+1-space ∈ list)
			remove_isolate = remove_isolate + 1*(i-1+space ∈ list)
			remove_isolate = remove_isolate + 1*(i-1-space ∈ list)
			if remove_isolate == 0
				filter!(e->e≠ i ,list)
			end
		end
		A=[HULL[i] for i in list]
		Ax=[A[i][j][1] for i in 1:size(A,1) for j in 1:size(A[1],1)]
		Ay=[A[i][j][2] for i in 1:size(A,1) for j in 1:size(A[1],1)]
		A=[[Ax[i],Ay[i]] for i in 1:size(Ay,1)]
		polytopes=push!(polytopes,A)
		list_HULLs=push!(list_HULLs,list)
	end


	Plots.plot!( VPolygon(polytopes[3]), alpha=.1,label="",color="red",subplot=2)
	Plots.plot!( VPolygon(polytopes[3]), alpha=.1,label="",color="red",subplot=3)
	for i in list_HULLs[3]
		if classes_matrices[3][i]>0
			Plots.plot!( VPolygon(HULL[i]), alpha=classes_matrices[3][i]/maximum(classes_matrices[3][list_HULLs[3]]),label="",color="red",subplot=2)
		end
	end

	Plots.plot!( VPolygon(polytopes[2]), alpha=.1,label="",color="blue",subplot=2)
	Plots.plot!( VPolygon(polytopes[2]), alpha=.1,label="",color="blue",subplot=3)
	for i in list_HULLs[2]
		if classes_matrices[2][i]>0
			Plots.plot!( VPolygon(HULL[i]), alpha=classes_matrices[2][i]/maximum(classes_matrices[2][list_HULLs[2]]),label="",color="blue",subplot=2)
		end
	end
	#
	# Plots.plot!( VPolygon(polytopes[1]), alpha=.1,label="",color="green",subplot=3)
	# for i in list_HULLs[1]
	# 	if classes_matrices[1][i]>0
	# 		Plots.plot!( VPolygon(HULL[i]), alpha=classes_matrices[1][i]/maximum(classes_matrices[1][list_HULLs[1]]),label="",color="green",subplot=2)
	# 	end
	# end
	plot!(xlim=[xmin,xmax],ylim=[ymin,ymax]) |>display
end

classes=unique([data[i][5] for i in 1:size(data,1)])
colls = ["green","blue","red"]
for j in 1:length(classes)
	for i in 1:size(data,1)
		if data[i][5] == classes[j]
			plot!(data[i][1],data[i][2],color=colls[j],xlabel="\$x_1\$",ylabel="\$x_2\$",w=.5,subplot=3,label="" ) |>display
		end
	end
end





# m1=(score[1]+score[2])/maximum(score[1]+score[2])
# plot(1:100,m1[sortperm(reshape((score[1]+score[2])/maximum(score[1]+score[2]),space*space))])
# m2=(score[1]+score[3])/maximum(score[1]+score[3])
# plot!(1:100,m2[sortperm(reshape((score[1]+score[3])/maximum(score[1]+score[3]),space*space))])
# m3=(score[2]+score[3])/maximum(score[2]+score[3])
# plot!(1:100,m3[sortperm(reshape((score[2]+score[3])/maximum(score[2]+score[3]),space*space))])
