using OrbitReadout, DataFrames, Plots, Random, Profile, PProf
## Orbits to generate the regions


	function mollifier(t, duration; pw = 20)
		if abs(t/duration) > 10
			return 0.
		else
			return 1. / (1. + (t/duration)^pw)
		end
	end

	function input_train(t, bapTimes, duration, amp)
		res = 0.
		Δ = duration / 2
		for ts in bapTimes
			res += mollifier(t - (ts + Δ), Δ)
		end
		return res * amp
	end



println("Generating 2D orbits for class definition - 3 classes A, B, C (no B nor A) ...")

	function dynamics!(du,u,p,t)
		du[1] = -.00015*(u[1]) + input_train(t,p[2],20.,p[1]/5.)
		du[2] = -.002*(u[2]) + input_train(t,p[2],20.,p[1]/5.)
	end

	list_weight=[0.1; -0.4; 0.75; 0.05;0.1; -0.4; 0.75; 0.05]
	list_class=["N"; "D"; "P"; "N";"N"; "D"; "P"; "N"]
	data=Vector{Array{Any,1}}[]
	pulses = 8
	after_pulses=1.2e4
	Random.seed!(1234)

pyplot()
	plot(layout = 4,windowsize=(700,300))

	for repeat in 1:4
		k=0
		for j in [collect(2: .6 : 4) ;collect((2 - .2): .6 : (4 - .2))]
				k = k +1
				u0 = [0,0]
				sequence = [100.0*i for i in 1:pulses]
				tspan = (0.0,sequence[end]+after_pulses)
				prob = OrbitReadout.ODEProblem(dynamics!,u0,tspan,(j+rand()/2 - 0.25,sequence))
				sol = OrbitReadout.solve(prob,reltol=1e-3,OrbitReadout.BS3())
				aux=DataFrame(time=sol.t,
				          var1=sol[1,:],
				          var2=sol[2,:])

			  plot!(aux[:,:time]/1000,aux[:,:var1],subplot=1,label="",xlabel="time (s)",ylabel="var 1")
			  plot!(aux[:,:time]/1000,aux[:,:var2],subplot=2,label="",xlabel="time (s)",ylabel="var 2")
			  plot!(aux[:,:var1],aux[:,:var2],subplot=3,label="",xlabel="var2",ylabel="var1")
			  global data=[data;[[aux[:,:var1],aux[:,:var2],aux[:,:time]/1000,list_weight[k],list_class[k]]]]
		end
	end
	println("done")
## Create a mesh of polygones given the max and min value among all trajectories

println("Creating a 2D mesh to contain all orbits - 15x15 (space variable) ...")
	space=25 #resolution
		mesh, xmin, xmax, ymin, ymax = OrbitReadout.polytopes_create(data,space; plot_ = true)
		println("done")

println("Evaluating the time spent of each orbit in each mesh element ...")
	u0 = zeros(size(data,1)*space*space)
    tspan = (0.0,after_pulses)
    prob = OrbitReadout.ODEProblem(time_spent_mesh!,u0,tspan,(size(mesh,1), mesh, data))
    sol = OrbitReadout.solve(prob,reltol=1e-3,saveat=10)
	println("done")

## Generate geometrical readout

println("Generating regions ...")
	readout = OrbitReadout.generate_geometrical_readout(1., 1, data, space, sol, mesh; plot_ = true)
	println("done")

## Testing it
# 
# println("Generating 2D orbits for class definition - 3 classes A, B, C (no B nor A) ...")
# 	Random.seed!(4322)
# 	pyplot()
# 	plot(layout = 4,windowsize=(700,300))
# 	plot!(OrbitReadout.VPolygon(readout[2]),subplot=3)
# 	plot!(OrbitReadout.VPolygon(readout[3]),subplot=3)
# 	#plot!(OrbitReadout.VPolygon(readout[1]),subplot=3)
# 	list_class=["N"; "D"; "P"; "N"]
# 	data=Vector{Array{Any,1}}[]
# 	for repeat in 1:10
# 		k=0
# 		for j in collect(2: .6 : 4)
# 				k = k +1
# 				u0 = [0,0]
# 				sequence = [100.0*i for i in 1:pulses]
# 				tspan = (0.0,sequence[end]+after_pulses)
# 				prob = OrbitReadout.ODEProblem(dynamics!,u0,tspan,(j+rand()/2 ,sequence))
# 				sol = OrbitReadout.solve(prob,reltol=1e-4,OrbitReadout.BS3())
# 				aux=DataFrame(time=sol.t,
# 				          var1=sol[1,:],
# 				          var2=sol[2,:])
#
#
# 			  plot!(aux[:,:time]/1000,aux[:,:var1],subplot=1,label="",xlabel="time (s)",ylabel="var 1")
# 			  plot!(aux[:,:time]/1000,aux[:,:var2],subplot=2,label="",xlabel="time (s)",ylabel="var 2")
# 			  plot!(aux[:,:var1],aux[:,:var2],subplot=3,label="",xlabel="var2",ylabel="var1") |>display
# 			  global data=[data;[[aux[:,:var1],aux[:,:var2],aux[:,:time]/1000,list_weight[k],list_class[k]]]]
# 		end
# 	end
#
# 	println("done")
#
# println("Evaluating the time spent of each orbit in each mesh element ...")
# 	u0 = zeros(size(data,1)*(size(readout,1)-1))
#     tspan = (0.0,after_pulses)
#     prob = OrbitReadout.ODEProblem(time_spent_mesh!,u0,tspan,((size(readout,1)-1), readout[2:3], data))
#     sol = OrbitReadout.solve(prob,reltol=1e-5,saveat=5)
# 	println("done")
#
# 	colls=range(HSV(140+90,1,1), stop=HSV(360,1,1), length=2)
#
# outcome=Array{String,1}[]
# 	for k in 1:10
# 		for j in 1:4
# 			  if sol[1+2*(j+k-2),end] < sol[2+2*(j+k-2),end]
# 				  global outcome = [outcome;"P"]
# 				  plot!(sol.t/1000,sol[2+2*(j+k-2),:],subplot=4,label="",color="red")
#
# 			  elseif  sol[1+2*(j+k-2),end] > sol[2+2*(j+k-2),end]
# 				  global outcome = [outcome;"D"]
# 				  plot!(sol.t/1000,sol[1+2*(j+k-2),:],subplot=4,label="",color="blue")
#
# 			  else
# 				  global outcome = [outcome;"N"]
# 				  plot!(sol.t/1000,sol[1+2*(j+k-2),:],subplot=4,label="",color="green")
# 			  end
# 		end
# 	end
#
#
# truee=[data[i][5] for i in 1:40]
# outcome
#
# count([truee[i]==outcome[i] for i in 1:40])
