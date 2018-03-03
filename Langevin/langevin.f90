program Langevin

integer it, npart
real x(1000), array(1000)
real dt, xmin, xmax, dU, T


it = 10
npart = 10
dt=1e-3

do i=1, npart
	x(0) = 0
	do j=1, it
		x(j)= evol(x(j-1), dt)
	enddo
enddo	



end


real function evol(x, dt)
	evol = x-dt*dU(x)/g1 + sqrt(2*T(x)*dt/(g1*g2))*e
end

real function dU(x)
	dU = 
end
