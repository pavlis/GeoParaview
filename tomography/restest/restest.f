	program restest

	integer SIZEOFA,NMAX,MMAX
	parameter(NMAX=1000000)
	parameter(MMAX=100000)
	parameter(SIZEOFA=100000000)
        common /lsqrdat/ a(SIZEOFA),ja(SIZEOFA),
     +                   na(MMAX),b(MMAX),jndx(NMAX)
c--this vector holds the model that defines this resolution test
	real*4 model(NMAX)
	integer imod
	real*4 mtmp
	integer i
c--these are matrix size variables
        integer m,n
	character*128 filedts,fileres,outfile,modfile
	integer lunres,lundts
	integer lout
	integer lresout,lunmod
	parameter(lunres=10,lundts=11,lresout=12,lunmod=13)

c--This should work as standard out on any unix box
        lout=6

	write(*,*) 'program restest'
	write(*,*) 'Enter file name holding sparse matrix coeficients'
	read(*,1000) filedts
 1000	format(a80)
	write(*,*) 'Enter file name of data residuals to be simulated'
	read(*,1000) fileres
	open (unit=lundts, file=filedts, access='sequential',
     +		form='unformatted')
	open (unit=lunres, file=fileres, access='sequential',
     +		form='unformatted')

c--This is the same routine used in runlsqr - used here for convenience
	call makea(m,n,lundts,lunres,lout)
c--initialize the model vector to zero so we can insert new values 
	do 100 i=1,n
		model(i)=0.0
  100	continue
c--Now we load a model vector by column index
	write(*,*) 'Enter file with cells to be set ',
     $     'for the resolution test'
	read(*,1000) modfile
	open(lunmod,file=modfile,form='formatted',access='sequential')
  150	read(lunmod,*,end=200) imod,mtmp
  		if(imod.gt.n) then
			write(*,*)'illegal model index=',imod
			write(*,*)'Number of columns = ',n
			write(*,*)'Abort with no output'
		endif
		model(imod)=mtmp
	goto 150
  200	continue

c
c--clear the right hand side vector to 0 so when we call aprod 
c--we get G*m
c
	do 250 i=1,m
		b(i)=0.0
  250	continue
  	call aprod(1,m,n,model,b)

	write(*,*) 'Enter file name of fake residuals for ',
     $     'this resolution test'
	read(*,1000) outfile
	open (unit=lresout, file=outfile, access='sequential',
     $		form='unformatted',status='new')

	do 300 i=1,m
		write(lresout) b(i)
  300	continue
  	end


       subroutine aprod(mode,m,n,x,y)
c
       integer mode,m,n
       real x(n),y(m)
c
c  aprod performs the following functions:
c
c      if mode = 1, set y = y + a*x
c      if mode = 2, set x = x + a(transpose) * y
c
c  where a is a matrix stored by rows in the arrays a, ja, and na.
c  in this example, a, ja, na are stored  in common.
c
	integer SIZEOFA,NMAX,MMAX
	parameter(NMAX=1000000)
	parameter(MMAX=100000)
	parameter(SIZEOFA=100000000)
        common /lsqrdat/ a(SIZEOFA),ja(SIZEOFA),
     +                   na(MMAX),b(MMAX),jndx(NMAX)


	real a
	integer ja,na

c
       integer i,j,l,l1,l2
       real sum,yi,zero
c
       zero = 0.0
       l2 = 0
       if(mode.ne.1) go to 400
c
c--------------------------------------------------------------------------
c mode = 1  -- set y = y + a*x.
c--------------------------------------------------------------------------
c
	do 200 i = 1,m
	  if (na(i).gt.0) then
	    sum = zero
	    l1 = l2 + 1
	    l2 = l2 + na(i)
	    do 100 l = l1,l2
	      j = ja(l)
	     sum = sum + a(l) * x(j)
100	   continue
	   y(i) = y(i) + sum
	 end if
200	continue
	return
c
c--------------------------------------------------------------------------
c  mode = 2  --  set x = x + a(transpose) * y
c---------------------------------------------------------------------------
c
400	do 600 i = 1,m
	  if (na(i).gt.0) then
	    yi = y(i)
	    l1 = l2 + 1
	    l2 = l2 + na(i)
	    do 500 l = l1,l2
	      j = ja(l)
	      x(j) = x(j) + a(l) * yi
500	    continue
	  end if
600	continue
	return
c
c       end of aprod
	end
	subroutine makea(m,nbl,lundts,lunres,lout)
c
c---------------------------------------------------------------------------
c
c  a is the sparse matrix whose nonzero coefficients are stored by
c  rows.  let a have m rows, n columns, and nz nonzeroes.  we need
c  three arrays dimensioned as real a(nz), and integer ja(nz),
c  na(m) where:
c
c	a(l) is the lth nonzero of a, counting across row 1, then
c	 row 2, and so on.
c
c	ja(l) is the column in which the lth nonzero of a lies.
c
c	na(i) is the number of nonzero coefficients  in the ith row
c	 of a.
c
c  (see lsqr alg. paige and saunders p.19)
c  these re made available to aprod through common.  the actual array
c  dimensions are:
c
c	a(SIZEOFA),ja(SIZEOFA), and na(m)
c
c  SIZEOFA is the probable largest number of nonzero coefficients in a:
c  this is dependent upon details of how a is formed
c
c  other input dimensional parameters are:
c
c	m - the row dimension of calling program.  If the number
c	    of rows read in exceed this the program is aborted.
c	n - the number of columns in a.  The program is aborted if
c	    any elements of ja() exceed n
c
c  on output:
c
c	mm - counts the number of rows processed, used in error check
c	l - counts the total number of nonzero coefficients processed
c	    also used as an error check
c	ja - as described above
c	a - as described above
c	na - as described above
c	b - the rhs vector
c
c	lundts = velocity element matrix
c	lunres = data vector
c----------------------------------------------------------------------------
c
	real a
	integer ja,na
	integer SIZEOFA,NMAX,MMAX
	parameter(NMAX=1000000)
	parameter(MMAX=100000)
	parameter(SIZEOFA=100000000)
	common /lsqrdat/ a(SIZEOFA),ja(SIZEOFA),
     +                   na(MMAX),b(MMAX),jndx(NMAX)


	l = 0
	mm = 1
c
c  read in velocity information, along with the data vector, first
c
100	read(lundts,end = 195) namm
	if (namm.eq.-1) go to 150
	read(lunres,end = 196) b(mm)
c	write(*,*) 'namm, mm, b(mm) = ', namm, mm, b(mm)
c	read(*,*) pause
c
c  test to check for overflow of array
c
	l1 = l + 1
	l2 = l + namm
	l3 = l2

	if(l3.gt.SIZEOFA) then
	  write(*,*) 'FATAL ERROR (makea):  work array a is full'
	  write(*,*) 'Increase size parameter SIZEOFA in source code'
	  write(lout,*) 'FATAL ERROR (makea):  work array a is full'
	  write(lout,*) 'Increase size parameter SIZEOFA in source code'
	  stop
	else
	  read(lundts) (ja(i),a(i),i=l1,l2)
c	  write(*,*) ' l1, l2 = ', l1, l2
c	  write(*,*) ' ja(l1), a(l1) = ', ja(l1), a(l1)
c	  write(*,*) ' ja(l2), a(l2) = ', ja(l2), a(l2)
	  do 110 i = l1, l2
	    if(ja(i).gt.NMAX) then
	      write(*,*) 'FATAL ERROR (makea):  grid overflow'
	      write(*,*) 'Datum number (row) = ',mm
	      write(lout,*) 'FATAL ERROR (makea):  grid overflow'
	      write(lout,*) 'Datum number (row) = ',mm
	      stop
	    endif
c	    a(i) = a(i)/scale(ja(i))
110	  continue
	  na(mm) = namm

	  l = l2
	  mm = mm + 1
c
c  error exit for too many data values (too many rows)
c
	  if(mm.gt.MMAX) then
	    write(*,*) 'FATAL ERROR (makea):  Row overflow'
	    write(*,*) 'Maximum allowed number of rows = ',MMAX
	    write(lout,*) 'FATAL ERROR (makea):  Row overflow'
	    write(lout,*) 'Maximum allowed number of rows = ',MMAX
	    stop
	  endif
	end if
	go to 100

c------now read in the constraint matrix
150	write(*,*) mm-1, ' rows of vel, hyp, and data read in '
	write(lout,*) mm-1, ' rows of vel, hyp, and data read in '
c---read in the index array
	read(lundts) nbl
	read(lundts) (jndx(i),i=1,nbl)

	go to 200

c----Abnormal Endings
195	write(*,*) ' Error: Ran out of velocity info! '
	write(lout,*) ' Error: Ran out of velocity info! '
	stop
196	write(*,*) ' Error: Ran out of data! '
	write(lout,*) ' Error: Ran out of data! '
	stop
197	write(*,*) ' Error: Ran out of Hypo info '
	write(lout,*) ' Error: Ran out of Hypo info '
	stop


200	m = mm - 1
	write(*,*) ' A total of ',m,' rows read in '
	write(*,*) ' A total of ',l,' elements of A read in '
	write(*,*) ' Total number of variables =   ',nbl
	write(lout,*) ' A total of ',m,' rows read in '
	write(lout,*) ' A total of ',l,' elements of A read in '
	write(lout,*) ' Total number of variables =   ',nbl
	
	return
	end
