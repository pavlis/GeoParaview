c---Program to convert db tables into a standard arrival time
	include "db.i"

	character*6 sta
	character*1 phs
	real*8 time, dsec, otime
	integer iyr,jday,ihr,imn
	real*8 xlat,xlon,depth
	real*8 wgt,deltim
	integer evid,prefid,arrid
	integer db(4)
	integer dbev(4),dbo(4),dbassoc(4),dbarrival(4)
	integer dbs1(4),dbs2(4)
	character*80 dbnm
	character*80 ofname
	integer irs,ire
	integer lastevid

	write(*,*) "Database to pickfile converter"
	write(*,*) "Enter database name to convert"
	read(*,*) dbnm
	write(*,*) "Will read from database ",dbnm
	write(*,*) "Enter file to contain output"
	read(*,*) ofname
	write(*,*) "Writing to file ",ofname

	ierr = dbopen(dbnm,"r",db)
	if(ierr.ne.0) then
		write(*,*) "Open failure on database"
		stop
	endif
	open(12,file=ofname)
c-- get all the database pointes for event, orign, assoc, and arrival
	call dblookup(dbev,db,"","event","","")
	call dblookup(dbo,db,"","origin","","")
	call dblookup(dbassoc,db,"","assoc","","")
	call dblookup(dbarrival,db,"","arrival","","")
	call dbjoin(dbs1,dbev,dbo,0,0,0,0,"")
	call dbsubset(dbs2,dbs1,"orid==prefor","")
	call dbjoin(dbs1,dbs2,dbassoc,0,0,0,0,"")
	call dbjoin(db,dbs1,dbarrival,0,0,0,0,"")

	call dbquery(db,dbRECORD_COUNT,nrec)
	write(*,*) "Working database view has ",nrec," rows"
	irs=1
	do 500 i=1,nrec
		db(4)=i
		if(i.eq.1) then
			ierr=dbgetv(dbb,"","evid",lastevid,0)
		endif
		ierr=dbgetv(db,"","evid",evid,0)
		if((evid.ne.lastevid) .or. (i.eq.nrec) )then
			if(i.eq.nrec) then
				ire = i
			else
				ire = i - 1
			endif
			do 400 j=irs,ire
				db(4)=j
				ierr=dbgetv(db,"","arrival.time",
     1					"lat",xlat,"lon",xlon,
     2					"depth",depth,"origin.time",otime,
     3					"evid",evid,"prefor",prefid,
     4					"sta",sta,"arid",arrid,
     5					"jdate",jdate,"phs",phs,
     6					"deltim",deltim,0)
				wgt = 1.0/deltim
				if(j.eq.irs) then
					call etoh(otime, iyr, jday, ihr, imn, dsec)
					write(12,200) iyr, jday, ihr, imn, dsec,
     +						xlat, xlon, depth, evid,prefid
				endif
				call etoh(time, iyr, jday, ihr, imn, dsec)
				if (phs.ne.'P'.and.phs.ne.'S') then
					write(*,*)"Unknown phase = ",phs," at record ",db(4)
					write(*,*)"Datum skipped"
				else
					write(12,202) sta, iyr, jday, ihr, imn, dsec, 
     +						phs, rwt, arrid

				endif
  400			continue
		endif
		write(12,*) ''
  500	continue
	
200     format(2i4,2i3,f9.4,f10.5,1x,f10.5,f9.4,6x,2i10)
202     format(a6,2i4,2i3,f8.3,1x,a1,8x,f8.3,1x,i8)
	end
	  
	subroutine etoh(epoch, iyear, iday, ihour, imin, sec)

	integer iyear, iday, ihour, imin
	real*8 epoch, sec

	integer date
	integer year
	integer month
	character*3 mname
	integer day
	integer doy
	integer hour
	integer minute
	real second
	integer i,dim,leap

	integer dayinmon(13)
	character*3 monname(12)

	data dayinmon /31,28,31,30,31,30,31,31,
     +                 30,31,30,31,31/
	data monname /"Jan","Feb","Mar","Apr","May","Jun",
     +                "Jul","Aug","Sep","Oct","Nov","Dec"/

	integer diy
	real*8 secleft

	doy = epoch / 86400.
	secleft = mod(epoch,86400.0)
	hour = 0.
	minute = 0.
	second = 0
c---compute hours minutes seconds
	if (secleft.ne.0.) then
c---before 1970		
		if (secleft.lt.0) then	
c---subtract a day
			doy = doy - 1
c----add a day
			secleft = secleft + 86400
		endif
		hour = secleft/3600
		secleft = mod(secleft,3600.0)
		minute = secleft/60
		second = mod(secleft,60.0)
	endif

	if (doy.ge.0) then
		year = 1970
5	        isl = isleap (year)
	        if (isl.eq.1) then
		  diy = 366
		else
		  diy = 365
	        endif
		if (doy.lt.diy) go to 10
		doy = doy - diy
		year = year + 1
		go to 5
	else
		year = 1969
7	        isl = isleap (year)
	        if (isl.eq.1) then
		  diy = 366
		else
		  diy = 365
	        endif
		if (doy.ge.0) go to 10
		doy = doy + diy
		year = year - 1
		go to 7
	endif
10	doy = doy + 1
	date = year * 1000 + doy

	leap = isleap(year)
	day = doy
	do i = 1, 12
		dim = dayinmon(i)
		if( leap.eq.1.and.i.eq.1 ) dim = dim + 1
		if( day.le.dim )go to 12
		day = day - dim
	enddo
12	month = i + 1
	mname = monname(i)

	iyear = year
	iday  = doy
	ihour = hour
	imin  = minute
	sec   = second

	return
	end

c	function mod(a,b)

c	imod = a/b
c	amod = a - imod * b

c	return
c	end

	function isleap(year)
	integer year

	if (mod(year,4).eq.0.and.mod(year,100).ne.0.or.mod(year,400).eq.0) then
	  isleap = 1
	else
	  isleap = 0
	endif
c	return(year % 4 == 0 && year % 100 != 0 || year % 400 == 0)

	return
	end

	subroutine htoe (iyear, iday, ihour, imin, sec, epoch)

	integer iyear
	integer iday
	integer ihour
	integer imin
	real*8 sec
	real*8 epoch

	integer jdate

	jdate = 1000*iyear + iday
	call dtoepoch (jdate, epoch)
	epoch = epoch + 3600.0*ihour + 60.0*imin + sec

	return
	end



c----convert julian date to epoch time
	subroutine dtoepoch(date, time)
	integer date
	real*8 time
	integer i,year,day,day

	i    = 0
	year = 0
	day  = 0
	days = 0

	year = date / 1000
	day  = mod(date,1000)

	if( year .gt. 1970 ) then
	  do i = 1970, year-1
	    days = days + 365
	    if (isleap(i).eq.1) days = days + 1
	  enddo
	endif
	if( year .lt. 1970 ) then
	  do i = year, 1969
	    days = days - 365
	    if (isleap(i).eq.1) days = days - 1
	  enddo
	endif

	days = days + day - 1
	time = days * 86400.
	return
	end
