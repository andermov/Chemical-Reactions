program Reaction_Simulation_Gillespie

!! Numerical simulation of the A+B->C chemical reaction network using the Gillespie algorithm.
!! There are 6 other reactions: A->0, 0->A, B->0, 0->B, C->0, 0->C.
!! The initial conditions of the concentrations are Poissonian.

real, dimension(3,1000000) :: X, correl, Xnew
real, dimension(3) :: InitialX, mean, correlmean
integer, dimension(7,3) :: p
real, dimension(7) :: a
real, dimension(1000000) :: t, tau
real :: r1, r2, taum, mu, A0, h, ti
real :: k1A, k2A, k1B, k2B, k1C, k2C, k3
integer :: i,j,k
integer :: i1,i2,i3
integer :: check

!! First column in X correspond to species A, second to B and third to C. 
!! Stores the copy numbers of each simulation in the 100,000 rows.
!! k1 creation rate of each species. k2 annihilation rate and k3 the rate of A+B->C.
k1A=8
k2A=3 
k1B=8
k2B=2
k1C=3
k2C=1.5
k3=0.5

!! Stoichiometric coefficient of each of the reactions.
p(1,:)=(/1,0,0/)
p(2,:)=(/0,1,0/)
p(3,:)=(/0,0,1/)
p(4,:)=(/-1,0,0/)
p(5,:)=(/0,-1,0/)
p(6,:)=(/0,0,-1/)
p(7,:)=(/-1,-1,1/)


!! Sampling of Poissonian initial conditions.

do k=1,1000000
	call PoissonSampling (k1A/k2A,k1B/k2B,k1C/k2C,i1,i2,i3)
	X(:,k)=(/i1,i2,i3/)
enddo

print *, "Initial conditions are :", X(:,1)

a(1)=k1A
a(2)=k1B
a(3)=k1C




open (unit=8,status="new", file="probpr", action="write")
open (unit=9,status="new", file="probcorrpr", action="write")

h=0.0005
t(:)=0
	mean(1)=(sum(X(1,1:1000000)))/1000000
	mean(2)=(sum(X(2,1:1000000)))/1000000
	mean(3)=(sum(X(3,1:1000000)))/1000000
write (unit=8, fmt=*) 0, mean(:)

call random_seed()
do k=1,1000000
!! First sample.
	a(4)=k2A*X(1,k)
	a(5)=k2B*X(2,k)
	a(6)=k2C*X(3,k)
	a(7)=k3*X(1,k)*X(2,k)
	A0=sum(a)
	call random_number(r1)
	tau(k)=1/A0*log(1/r1)
	t(k)=tau(k)+t(k)
	call random_number(r2)	
	mu=r2*A0
	do j=1,7
		if (mu .le. sum(a(1:j))) then
			Xnew(:,k)=X(:,k)+p(j,:)
			exit
		endif
	enddo
enddo


!! Time evolution.
!! i index corresponds to time evolution (then, i*h=real_time_of_the_dynamics), k for the number of parallel simulations running.
call random_seed()
do i=1,9999
	ti=i*h
	do k=1,1000000
		do while (t(k) .le. ti)
			!Reaction occurs
			X(:,k)=Xnew(:,k)
			!Sample
			!Computes the type of the next reaction.
			call random_number(r2)	
			a(4)=k2A*X(1,k)
			a(5)=k2B*X(2,k)
			a(6)=k2C*X(3,k)
			a(7)=k3*X(1,k)*X(2,k)
			A0=sum(a)		
			mu=r2*A0
			do j=1,7
				if (mu .le. sum(a(1:j))) then
					Xnew(:,k)=X(:,k)+p(j,:)
					exit
				endif	
			enddo
			!Computes the time of the next reaction.
			call random_number(r1)
			tau(k)=1/A0*log(1/r1)
			t(k)=tau(k)+t(k)
		enddo
		correl(1,k)=X(1,k)*X(2,k)
		correl(2,k)=X(1,k)*X(3,k)
		correl(3,k)=X(2,k)*X(3,k)
	enddo
	!Writing to files
	mean(1)=(sum(X(1,1:1000000)))/1000000
	mean(2)=(sum(X(2,1:1000000)))/1000000
	mean(3)=(sum(X(3,1:1000000)))/1000000
	write (unit=8, fmt=*) ti, mean(:)
	correlmean(1)=(sum(correl(1,1:1000000)))/1000000-mean(1)*mean(2)
	correlmean(2)=(sum(correl(2,1:1000000)))/1000000-mean(1)*mean(3)
	correlmean(3)=(sum(correl(3,1:1000000)))/1000000-mean(2)*mean(3)
	write (unit=9, fmt=*) ti, correlmean(:)
enddo


contains


!! Subroutine to sample from a Poissonian distribution given its mean.
subroutine PoissonSampling(x1,x2,x3,y1,y2,y3)
implicit none

real :: r1, r2, r3, mm, s2, s3, x1, x2, x3
integer ::  y1, y2, y3, i, j, k

call random_seed()
call random_number(r1)
call random_number(r2)
call random_number(r3)

mm=0
s2=0
s3=0
i=0
j=0
k=0
do while (1.gt.0) 
	mm=mm+((x1)**i)*exp(-x1)/fact(i)
	if (r1 .lt. mm) then
		exit
	else
		i=i+1
	endif
enddo
y1=i
do while (1.gt.0) 
	s2=s2+x2**j*exp(-x2)/fact(j)
	if (r2 .lt. s2) then
		exit
	else
		j=j+1
	endif
enddo
y2=j
do while (1.gt.0) 
	s3=s3+x3**k*exp(-x3)/fact(k)
	if (r3 .lt. s3) then
		exit
	else
		k=k+1
	endif
enddo
y3=k

return

end subroutine PoissonSampling

!! Function to compute the factorial of an integer.
function fact(i)
integer :: i, fact, j

fact=1
if (i==0) then
	fact=1
else
	do j=1, i
		fact=fact*j
	enddo
endif
return
end function fact


end program Reaction_Simulation_Gillespie

