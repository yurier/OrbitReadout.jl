using JLD2, FileIO, Plots
using LazySets, DifferentialEquations, ColorSchemes, Parameters

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
		for i in 1:10
		protocol=list[j]
		tt=load("$(i)_$(protocol).jld2","time")
			g=load("$(i)_$(protocol).jld2","CaN")
			h=load("$(i)_$(protocol).jld2","CaMKII")
			# idx=[]
			# for i in collect(0. : tt[end]/500: tt[end])
			# 	append!(idx, findall(x->x <=i, tt)[end])
			# 	@show i/ tt[end]
			# end
			global data=[data;[[g,h,tt,list_weight[j],list_class[j]]]]
		end
	end
