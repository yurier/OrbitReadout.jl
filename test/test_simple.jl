using OrbitReadout, DataFrames, Plots, Random, Profile, PProf
## Orbits to generate the regions

@time begin
	println("\n Generating 2D orbits for class definition - 3 classes A, B, C (no B nor A)")
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
				  global data=[data;[[aux[:,:var1],aux[:,:var2],aux[:,:time]/1000,list_weight[k],list_class[k]]]]
			end
		end

	printstyled(color=:green,"\n Done in: \n \n")
	end
## Create a mesh of polygones given the max and min value among all trajectories
@time begin
	println("\n Creating a 2D mesh to contain all orbits - 15x15 (space variable)")
	space=20 #resolution
	mesh, xmin, xmax, ymin, ymax = OrbitReadout.polytopes_create(data,space; plot_ = true)
	printstyled(color=:green,"\n Done in: \n \n")
end

@time begin
	println("\n Evaluating the time spent of each orbit in each mesh element")
	u0 = zeros(size(data,1)*space*space)
    tspan = (0.0,after_pulses)
    prob = OrbitReadout.ODEProblem(time_spent_mesh!,u0,tspan,(size(mesh,1), mesh, data))
    sol = OrbitReadout.solve(prob,reltol=1e-3,saveat=10)
	printstyled(color=:green,"\n Done in: \n \n")
end
## Generate geometrical readout

@time begin
	println("\n Generating regions")
	readout = OrbitReadout.generate_geometrical_readout(1., 1, data, space, sol, mesh; plot_ = true)
	printstyled(color=:green,"\n Done in: \n \n")
end
## Generating testing set
@time begin
	println("\n Generating testing set")
	list_class=["D" , "D" , "P" ,"P"]
	data_test=Vector{Array{Any,1}}[]
	pulses = 8
	after_pulses=1.2e4

	Random.seed!(65498732164)
		for repeat in 1:20
			k=0
			for j in [2.2, 2.4, 3.0, 3.2]
					k = k +1
					u0 = [0,0]
					sequence = [100.0*i for i in 1:pulses]
					tspan = (0.0,sequence[end]+after_pulses)
					prob = OrbitReadout.ODEProblem(dynamics!,u0,tspan,(j+rand()/2 - 0.25,sequence))
					sol = OrbitReadout.solve(prob,reltol=1e-3,OrbitReadout.BS3())
					aux=DataFrame(time=sol.t,
							  var1=sol[1,:],
							  var2=sol[2,:])

				  global data_test=[data_test;[[aux[:,:var1],aux[:,:var2],aux[:,:time]/1000,list_weight[k],list_class[k]]]]
			end
		end
	printstyled(color=:green,"\n Done in: \n \n")
end

## Generating testing set

@time begin
	println("\n Evaluating the time spent of each orbit in each mesh element")
	u0 = zeros(size(data_test,1)*(size(readout,1)-1))
    tspan = (0.0,after_pulses)
    prob = OrbitReadout.ODEProblem(time_spent_mesh!,u0,tspan,((size(readout,1)-1), readout[2:3], data_test))
    sol_test = OrbitReadout.solve(prob,reltol=1e-5,saveat=5)
	printstyled(color=:green,"\n Done in: \n \n")

	outcome=Array{String,1}[]
	for k in 1:(size(readout,1)-2)
		for j in 1:Int(size(sol_test[:,:],1)/2)
		#@show 1+2*(j+k-2), 2+2*(j+k-2), k, j
		  if sol_test[1+2*(j+k-2),end] < sol_test[2+2*(j+k-2),end]
			  global outcome = [outcome;"P"]
			 # plot!(sol_test.t/1000,sol_test[2+2*(j+k-2),:],subplot=4,label="",color="red")

		  elseif  sol_test[1+2*(j+k-2),end] > sol_test[2+2*(j+k-2),end]
			  global outcome = [outcome;"D"]
			  #plot!(sol_test.t/1000,sol_test[1+2*(j+k-2),:],subplot=4,label="",color="blue")
		  else
			  global outcome = [outcome;"N"] #or threshold
			  #plot!(sol_test.t/1000,sol_test[1+2*(j+k-2),:],subplot=4,label="",color="green")
		  end
		end
	end

	classes=([data_test[i][5] for i in 1:size(data_test,1)])
	println("\n Clasified the generated ", length(classes), " orbits with ", 100*sum(outcome .== classes)/length(classes), "% accuracy" )
	printstyled(color=:green,"\n Test finished \n \n")
end
