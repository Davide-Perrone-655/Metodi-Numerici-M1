	program ising
	implicit none
	
	integer, parameter :: L=70						!lunghezza del reticolo
	integer, parameter :: meas=10000					!numero di misure
	integer, parameter :: t_term=100					!ogni quanto prendo un dato (deve essere maggiore della correlazione)
	integer :: seed
	
	integer :: sigma(L,L)							!matrice delle configurazioni	
	integer :: nn(L), pp(L)						!vettori con le condizioni al bordo
	real :: beta, h_ext	
	
	real :: exp1(5), exp2(2) 						!esponenziali inizializzati

	real :: temp1, temp2
	integer :: i, j, s, iter
	real :: eps, flag1


	seed=125			!imposta il seme (12345 brutto)
	eps= 0.003			!passo della temperatura
	beta=0.39			!da dove parto
	iter=30			!numero di iterazioni
	
	h_ext=0.			!campo esterno
	
	
	call boundary_geom(nn, pp, L)
	

	
	open(44, file='shape70.dat')
	write(44,*) '#numero interazioni, L_dati /rep, beta_0, incremento'
	write(44,*) iter, meas, beta, eps
	close(44)
	
	
	
	open(22, file='mag_70.dat')
	open(33, file='ene_70.dat')
	
	do s=1,iter
		
		
		call initialize(sigma,L,seed)
		call exp_init(exp1, exp2, beta, h_ext)				!inizializzo gli esponenziali 
		
		do i=1, t_term								!piccolo warm up
			call newchainstep(sigma, nn, pp, L, seed, exp1, exp2)
		end do
		
		do i=1, meas
			do j=1,t_term
				call newchainstep(sigma, nn, pp, L, seed, exp1, exp2)
			end do
			
			temp1=E_val(sigma, h_ext, nn, pp, L)
			temp2=M_val(sigma,L)
			
			write(33,*) temp1
			write(22,*) temp2
			
			
		end do
		
		
		beta = beta + eps			!temperatura
		seed=seed + 1
		
	end do
	
	
	close(22)
	close(33)
	

	
	
C=================================================
	
	CONTAINS



c=========================================================
	!INIZIALIZZAZIONI
C=========================================================	

	!Condizioni al bordo
	subroutine boundary_geom(nn, pp, len_row)
	integer :: len_row, nn(len_row), pp(len_row)
	integer :: i,j
	
	do i=1,len_row
		nn(i)=i-1
		pp(i)=i+1
	end do
	nn(1)=len_row
	pp(len_row)=1
	end subroutine boundary_geom
	
	
	
	
	!inizializzo gli esponenziali
	subroutine exp_init(f1, f2, beta, h_ext)
	real :: f1(5), f2(2), beta, h_ext
	
	if (h_ext==0) then
		f2(1)=1
		f2(2)=1
	else 
		f2(1)=exp(-h_ext)
		f2(2)=exp(h_ext)
	end if
	
	
	f1(1)= exp(-beta*8)
	f1(2)= exp(-beta*4)
	f1(3)= 1
	f1(4)= exp(beta*4)
	f1(5)= exp(beta*8)
	end subroutine exp_init
	
	
	
	
	

	!inizializza la matrice di spin in modo random
	subroutine initialize(matrix, len_row,seed)
	integer :: len_row, seed
	integer :: matrix(len_row, len_row)
	integer :: i, j
	real :: num
	
	do i=1,len_row
		do j=1,len_row
			num=ran2(seed)
			if (num.le.0.5) then
				matrix(i,j)=1
			else
				matrix(i,j)=-1
			end if
		end do
	end do
	end subroutine initialize
	
	
C================================================
	!CALCOLO ENERGIA	
c=================================================

	function E_val(matrix, ext_field, nn, pp, len_row) result(ene)
	integer, intent(in) :: len_row
	integer, intent(in) :: matrix(len_row, len_row)
	integer, intent(in) :: nn(len_row), pp(len_row)
	real, intent(in) :: ext_field
	
	integer :: i, j, jp, jm , ip, im, f0, f1, temp
	integer :: ene1
	real :: ene
	
	ene1=0
	do i=1,len_row
		im=nn(i)
		ip=pp(i)
		do j=1,len_row
			jp=pp(j)
			jm=nn(j)
			
			f1=matrix(i,j)
			f0=matrix(i,jp) + matrix(i,jm) + matrix(ip,j) + matrix(im,j)
			temp=f1*(0.5*f0 + ext_field)
			ene1 = ene1 - temp
		end do
	end do
	
	ene= real(ene1)/(len_row**2)
	
	end function E_val
	
	
	
c=======================================================
	!CALCOLO MAGNETIZZAZIONE		
C=======================================================
	function M_val(matrix, len_row) result(mag_res)
	
	integer, intent(in) :: len_row
	integer, intent(in) :: matrix(len_row, len_row)
	integer :: i_m, j_m 
	real :: F, G, temp, mag, mag_res
	
	F=len_row**2
	mag=0
	
	do i_m = 1, len_row
		do j_m = 1, len_row
			temp = real(matrix(i_m, j_m))/F
			mag = mag + temp
		end do
	end do
	mag_res=abs(mag)
	
	end function M_val
		
	
	
C==========================================================	
	!PASSO BASE METROPOLIS
C==========================================================
	subroutine chainstep(matrix, nn, pp, row, seed, f1, f2)
	
	integer :: matrix(row,row), nn(row), pp(row)
	integer :: row, seed
	!real :: ext_field, b
	real :: f1(5), f2(2)
	real :: z, r, dE
	integer :: i0, j0, ip, jm, im, jp, f0, T, temp, s
	
	do s=1,row*row
		i0=ceiling(ran2(seed)*row)				!prendo un elemento random
		j0=ceiling(ran2(seed)*row)				!della matrice
		z=ran2(seed)
		temp=matrix(i0,j0)
		
		!chiamo le funzioni per assicurare la corretta geometria sui bordi
		im=nn(i0)
		ip=pp(i0)
		jp=pp(j0)
		jm=nn(j0)
		
		
		!calcolo l'energia sommando i 4 più vicini ed il campo esterno 
		f0=matrix(i0,jp) + matrix(i0,jm) + matrix(ip,j0) + matrix(im,j0)
		T=(temp*f0)/2 
		f0= -T + 3 
		
		!write(*,*) (f0-3), f0
		!dE=2*temp*(f0 + ext_field)
		!r=exp(-b*dE)
		
		if (temp>0) then
			r=f1(f0)*f2(1)
		else
			r=f1(f0)*f2(2)
		end if 
		
		!calcolo la probabilità di accettare il valore

		if(z.le.r) matrix(i0,j0)= -temp
		
	end do
	
	end subroutine chainstep

C==========================================================	
	!PASSO BASE METROPOLIS NON RANDOMICO
C==========================================================
	subroutine newchainstep(matrix, nn, pp, row, seed, f1, f2)
	
	integer :: matrix(row,row), nn(row), pp(row)
	integer :: row, seed
	!real :: ext_field, b
	real :: f1(5), f2(2)
	real :: z, r, dE
	integer :: i0, j0, ip, jm, im, jp, f0, T, temp, s, i
	
	do s=1,row
		do i=1, row
			i0=i						!prendo un elemento
			j0=s						!della matrice in sequenza
			z=ran2(seed)
			temp=matrix(i0,j0)
		
			!chiamo le funzioni per assicurare la corretta geometria sui bordi
			im=nn(i0)
			ip=pp(i0)
			jp=pp(j0)
			jm=nn(j0)
		
		
			!calcolo l'energia sommando i 4 più vicini ed il campo esterno 
			f0=matrix(i0,jp) + matrix(i0,jm) + matrix(ip,j0) + matrix(im,j0)
			T=(temp*f0)/2 
			f0= -T + 3 
		
			if (temp>0) then
				r=f1(f0)*f2(1)
			else
				r=f1(f0)*f2(2)
			end if 
		
			!calcolo la probabilità di accettare il valore

			if(z.le.r) matrix(i0,j0)= -temp
		end do	
	end do
	
	end subroutine newchainstep	

C===========================================================================
C	FUNZIONI UTILI
C=========================================================================
	
	!stampa la matrice di spin
	subroutine printmatrix(matrix, len_row)
	integer, intent(in) :: len_row
	integer, intent(in) :: matrix(len_row, len_row)
	integer :: i, j
	
	do i=1,len_row
		do j=1,len_row
			if (matrix(i,j)==1) then
				write(*, '(A)', advance='NO') ' + '
			else
				write(*,'(A)', advance='NO') ' 0 '
			end if
		end do
		write(*,*) '|'
	end do
	write(*,*) ' '
	end subroutine printmatrix
	

c==============================================================
C	Random number generator: ran2 di numerical recipes
c==============================================================
	FUNCTION ran2(idum)
	INTEGER :: idum
	REAL :: ran2
	INTEGER, PARAMETER :: IM1=2147483563, IM2=2147483399, IR1=12211
	INTEGER, PARAMETER :: IR2=3791, NTAB=32,  IMM1=IM1-1
	INTEGER,PARAMETER :: IA1=40014,IA2=40692,IQ1=53668,IQ2=52774
	real, parameter :: AM=1./IM1, NDIV=1+IMM1/NTAB, EPS=1.2e-7, RNMX=1.-EPS
	
	INTEGER :: idum2,j,k,iv(NTAB),iy 
	SAVE iv,iy,idum2
	DATA idum2/123456789/, iv/NTAB*0/, iy/0/
	if (idum.le.0) then
		idum=max(-idum,1)
		idum2=idum
		do j=NTAB+8,1,-1
			k=idum/IQ1
			idum=IA1*(idum-k*IQ1)-k*IR1
			if (idum.lt.0) idum=idum+IM1
			if (j.le.NTAB) iv(j)=idum
		end do
		iy=iv(1)
	endif
	k=idum/IQ1
	idum=IA1*(idum-k*IQ1)-k*IR1
	if (idum.lt.0) idum=idum+IM1
	k=idum2/IQ2
	idum2=IA2*(idum2-k*IQ2)-k*IR2
	if (idum2.lt.0) idum2=idum2+IM2
	j=1+iy/NDIV
	iy=iv(j)-idum2
	iv(j)=idum
	if(iy.lt.1)iy=iy+IMM1
	ran2=min(AM*iy,RNMX)
	return
	end function
	
	
	end program ising
