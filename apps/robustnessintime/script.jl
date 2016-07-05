#!/home/cdesantana/Softwares/julia/usr/bin/julia

infile=open("script.sh","a");
nreal=2
meangamma= -0.8 
meanlambda= 0.12
#lambdastdev = [0;logspace(-2,-1,100)];#20 values
#gammastdev = [0;logspace(-2,-1,100)];#20 values

lambdastdev = [0;collect(0.05:0.001:0.1)];#100 values
gammastdev = [0;collect(0.05:0.001:0.1)];#100 values


for sdgamma in gammastdev 
    for sdlambda in lambdastdev 
	println(infile,"./robustnessintime 2 $nreal 1 10 1 $meangamma $meanlambda 5 5 $sdgamma $sdlambda > log || true");
    end
end
