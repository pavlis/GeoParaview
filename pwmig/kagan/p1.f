      implicit none

      integer iel,inf
      parameter(iel=10*10*10,inf=11*11*11)
      integer nelem(iel,8)
      integer nx,ny,nz,nn,ix,iy,iz,ip,i,j
      double precision xcoord(inf,3),xx(3),func(inf),value
      double precision coord(8,3),func_element(8)
c
c     nx, ny ,nz elements in x,y,z directions, respectively. 
c
      nx=10
      ny=10
      nz=10
c
c     nn=number of nodes
c
      nn=(nx+1)*(ny+1)*(nz+1)
c
c     xcoord(i,j) -> i=1,nn , j=1,3
c
      CALL READNN(inf,nx,ny,nz,XCOORD) 
c
c     nelem(i,j) -> i=1,number of elements, j=1,8
c
      call READELR(NX,NY,NZ,iel,NELEM)
c
c     func(i) -> function to be interpolated  i=1,nn
c
      call simple_function(inf,nn,xcoord,func)
c
c     xx(i) -> i=1,3 coordinates
c
      xx(1)=2.
      xx(2)=3.5
      xx(3)=2.4
      
      ix=3
      iy=4
      iz=3
      ip=ix+(iy-1)*nx+(iz-1)*nx*ny
      do i=1,8
         func_element(i)=func(nelem(ip,i))
         do j=1,3
            coord(i,j)=xcoord(nelem(ip,i),j)
         enddo
      enddo
c
c     coord(i,j) i=1,8
c     j=1,3
c     func_element(i) i=1,8
c
      call interpolate(xx,func_element,coord,value)
      
      end
c xx - point to interpolate
c func_element = value of f at each node
c coord = array holding actual coordinates of 8 corner points
c value - result

      subroutine interpolate(xx,func_element,coord,value)
      implicit none
      integer i,itno,j
      double precision x(3),xx(3),xdum
      double precision func_element(8),value,coord(8,3)
      double precision der(3,8),fun(8),jac(3,3)
      double precision residual(3),solution(3),er,det,jac1(3,3),ff
      itno=0
      do i=1,3
         x(i)=0.
      enddo

      do i=1,8
         write(*,2) (coord(i,j),j=1,3),func_element(i)
      enddo
 2    format(4e15.7)

 1    itno=itno+1
c
c     given x, derivatives of the interpolation functions with
c     respect to k1,k2,k3
c     der(i,j) -> i=1,3 j=1,8
c     
      call fmlin3(der,3,fun,x)
c
c     matmul -> der x coord = jac
c
      CALL MATMUL(DER,3,COORD,8,JAC,3,3,8,3)   
c
c     x=coord^T*F
c
      do i=1,3
         xdum=0.
         do j=1,8
            xdum=xdum+coord(j,i)*fun(j)
         enddo
         residual(i)=xdum
      enddo
c
c
c
      er=0.
      do i=1,3
         residual(i)=-(residual(i)-xx(i))
         er=er+residual(i)*residual(i)
      enddo
c
c     treex3 -> Jac^-1=Jac1
c
      CALL TREEX3(JAC,3,JAC1,3,DET)     
c
c     mvmult -> matrix x vector
c     jac*residual = solution
c
      call MVMULT(JAC1,3,residual,3,3,solution)

      do i=1,3
         x(i)=x(i)+solution(i)
      enddo

      if(er.gt.1.0e-03)goto 1
      write(*,*) 'itno = ',itno
      call fmlin3(der,3,fun,x)
      CALL matmul(DER,3,COORD,8,JAC,3,3,8,3)   

      do i=1,3
         xdum=0.
         do j=1,8
            xdum=xdum+coord(j,i)*fun(j)
         enddo
         residual(i)=xdum
      enddo

      do i=1,3
         xdum=0.
         do j=1,8
            xdum=xdum+coord(j,i)*fun(j)
         enddo
         residual(i)=xdum
         write(*,*) i,residual(i),xx(i)
      enddo

      er=0.
      do i=1,3
         residual(i)=-(residual(i)-xx(i))
         er=er+residual(i)*residual(i)
      enddo
      value=0.
      do i=1,8
         value=value+func_element(i)*fun(i)
      enddo

      write(*,*) 'function      ',value
      write(*,*) 'interpolation ',ff(xx(1),xx(2),xx(3))
      return
      end

      SUBROUTINE GEO(IP,XCOORD,GG,COORD,ICOORD,NELEM,IEL,
     +NF,INF)
C
C     STEERING VECTOR AND COORDINATES
C 
      IMPLICIT NONE
      INTEGER INF,ICOORD,IEL
      DOUBLE PRECISION XCOORD(INF,*),COORD(ICOORD,*)
      INTEGER GG(24),NELEM(IEL,*),NF(*) 
      INTEGER INC,I,J,IP

      INC=0
      DO 1 I=1,8
      INC=INC+1
      GG(INC)=NF(NELEM(IP,I))
 1    CONTINUE
      DO 2 I=1,8
      DO 2 J=1,3
      COORD(I,J)=XCOORD(NELEM(IP,I),J)
 2    CONTINUE
 5    FORMAT(24I3)
      RETURN
      END

      SUBROUTINE FMLIN3(DER,IDER,FUN,xsi1)
C
C      FORMS THE SHAPE FUNCTIONS AND THEIR
C      DERIVATIVES FOR 8-NODED BRICK ELEMENTS
C
      IMPLICIT NONE
      INTEGER IDER
      DOUBLE PRECISION DER(IDER,8),FUN(8)
      DOUBLE PRECISION ETA,XI,ZETA
      DOUBLE PRECISION ETAM,XIM,ZETAM
      DOUBLE PRECISION ETAP,XIP,ZETAP,xsi1(3)
      XI=xsi1(1)
      ETA=xsi1(2)
      ZETA=xsi1(3)
      ETAM=1.-ETA
      XIM=1.-XI
      ZETAM=1.-ZETA
      ETAP=ETA+1.
      XIP=XI+1.
      ZETAP=ZETA+1.
      FUN(1)=.125*XIM*ETAM*ZETAM
      FUN(2)=.125*XIM*ETAM*ZETAP
      FUN(3)=.125*XIP*ETAM*ZETAP
      FUN(4)=.125*XIP*ETAM*ZETAM
      FUN(5)=.125*XIM*ETAP*ZETAM
      FUN(6)=.125*XIM*ETAP*ZETAP
      FUN(7)=.125*XIP*ETAP*ZETAP
      FUN(8)=.125*XIP*ETAP*ZETAM
      DER(1,1)=-.125*ETAM*ZETAM
      DER(1,2)=-.125*ETAM*ZETAP
      DER(1,3)=.125*ETAM*ZETAP
      DER(1,4)=.125*ETAM*ZETAM
      DER(1,5)=-.125*ETAP*ZETAM
      DER(1,6)=-.125*ETAP*ZETAP
      DER(1,7)=.125*ETAP*ZETAP
      DER(1,8)=.125*ETAP*ZETAM
      DER(2,1)=-.125*XIM*ZETAM
      DER(2,2)=-.125*XIM*ZETAP
      DER(2,3)=-.125*XIP*ZETAP
      DER(2,4)=-.125*XIP*ZETAM
      DER(2,5)=.125*XIM*ZETAM
      DER(2,6)=.125*XIM*ZETAP
      DER(2,7)=.125*XIP*ZETAP
      DER(2,8)=.125*XIP*ZETAM
      DER(3,1)=-.125*XIM*ETAM
      DER(3,2)=.125*XIM*ETAM
      DER(3,3)=.125*XIP*ETAM
      DER(3,4)=-.125*XIP*ETAM
      DER(3,5)=-.125*XIM*ETAP
      DER(3,6)=.125*XIM*ETAP
      DER(3,7)=.125*XIP*ETAP
      DER(3,8)=-.125*XIP*ETAP
      RETURN
      END

      SUBROUTINE TREEX3(JAC,IJAC,JAC1,IJAC1,DET)
C
C      FORMS THE INVERSE OF A 3 BY 3 MATRIX
C
      IMPLICIT NONE
      INTEGER IJAC,IJAC1
      DOUBLE PRECISION JAC(IJAC,*),JAC1(IJAC1,*)
      INTEGER K,L
      DOUBLE PRECISION DET

      DET=JAC(1,1)*(JAC(2,2)*JAC(3,3)-JAC(3,2)*JAC(2,3))
      DET=DET-JAC(1,2)*(JAC(2,1)*JAC(3,3)-JAC(3,1)*JAC(2,3))
      DET=DET+JAC(1,3)*(JAC(2,1)*JAC(3,2)-JAC(3,1)*JAC(2,2))
      JAC1(1,1)=JAC(2,2)*JAC(3,3)-JAC(3,2)*JAC(2,3)
      JAC1(2,1)=-JAC(2,1)*JAC(3,3)+JAC(3,1)*JAC(2,3)
      JAC1(3,1)=JAC(2,1)*JAC(3,2)-JAC(3,1)*JAC(2,2)
      JAC1(1,2)=-JAC(1,2)*JAC(3,3)+JAC(3,2)*JAC(1,3)
      JAC1(2,2)=JAC(1,1)*JAC(3,3)-JAC(3,1)*JAC(1,3)
      JAC1(3,2)=-JAC(1,1)*JAC(3,2)+JAC(3,1)*JAC(1,2)
      JAC1(1,3)=JAC(1,2)*JAC(2,3)-JAC(2,2)*JAC(1,3)
      JAC1(2,3)=-JAC(1,1)*JAC(2,3)+JAC(2,1)*JAC(1,3)
      JAC1(3,3)=JAC(1,1)*JAC(2,2)-JAC(2,1)*JAC(1,2)
      DO 1 K=1,3
      DO 1 L=1,3
      JAC1(K,L)=JAC1(K,L)/DET
    1 CONTINUE
      RETURN
      END      

      SUBROUTINE READELR(NX,NY,NZ,iel,NELEM)
C
C     FOR A DOMAIN WITH NX, NY, AND NZ IN X Y AND Z DIR.
C
      IMPLICIT NONE
      INTEGER IEL
      INTEGER N(8),NELEM(IEL,*)
      INTEGER NX,NY,NZ,nnn
      INTEGER II,I,K,J

      DO K=1,NZ
         DO J=1,NY
            DO  I=1,NX
               N(1)=(K-1)*(NX+1)*(NY+1)+(J-1)*(NX+1)+I
               N(2)=N(1)+(NX+1)*(NY+1)
               N(3)=N(2)+1
               N(4)=N(1)+1
               N(5)=N(1)+(NX+1)
               N(6)=N(5)+(NX+1)*(NY+1)
               N(7)=N(6)+1
               N(8)=N(5)+1
               NNN=NNN+1
               DO II=1,8
                  NELEM(NNN,II)=N(II)
               enddo
            enddo
         enddo
      enddo

      RETURN
      END
      SUBROUTINE READNN(inf,nx,ny,nz,XCOORD)
C
C     READS NODAL POINT DATA
C
      IMPLICIT NONE
      INTEGER INF,NX,NY,NZ
      DOUBLE PRECISION  XCOORD(INF,*),dx,dy,dz
      INTEGER i,j,k,nnn

      dx=1.
      dy=1.
      dz=1.
      nnn=0
      do k=1,nz+1
         do j=1,ny+1
            do i=1,nx+1
               nnn=nnn+1
               xcoord(nnn,1)=dx*(i-1)
               xcoord(nnn,2)=dy*(j-1)
               xcoord(nnn,3)=dz*(k-1)
            enddo
         enddo
      enddo
      return
      end

      subroutine simple_function(inf,nn,xcoord,func)
      implicit none
      integer inf,nn,i
      double precision xcoord(inf,3),func(*),x,y,z,ff

      do i=1,nn
         x=xcoord(i,1)
         y=xcoord(i,2)
         z=xcoord(i,3)
         func(i)=ff(x,y,z)
      enddo

      return
      end

      double precision function ff(x,y,z)
      double precision x,y,z
c
c     ff=xy-yxz^2+1+4x-x^2*y^2
c
      ff=x*y-y*x*z*z+1.+4*x-x*x*y*2.
      return
      end

      SUBROUTINE MATMUL(A,IA,B,IB,C,IC,L,M,N)
C
C      FORMS THE PRODUCT OF TWO MATRICES
C
      IMPLICIT NONE
      INTEGER IA,IB,IC
      DOUBLE PRECISION  A(IA,*),B(IB,*),C(IC,*)
      INTEGER L,M,N,I,J,K
      DOUBLE PRECISION X

      DO 1 I=1,L
      DO 1 J=1,N
      X=0.0
      DO 2 K=1,M
    2 X=X+A(I,K)*B(K,J)
      C(I,J)=X
    1 CONTINUE
      RETURN
      END

      SUBROUTINE MATRAN(A,IA,B,IB,M,N)
C
C      FORMS THE TRANSPOSE OF A MATRIX
C
      IMPLICIT NONE
      INTEGER IA,IB
      DOUBLE PRECISION A(IA,*),B(IB,*)
      INTEGER M,N,I,J

      DO 1 I=1,M
      DO 1 J=1,N
    1 A(J,I)=B(I,J)
      RETURN
      END

      SUBROUTINE MVMULT(M,IM,V,K,L,Y)
C
C      MULTIPLIES A MATRIX BY A VECTOR
C
      IMPLICIT NONE
      INTEGER IM
      DOUBLE PRECISION M(IM,*),V(*),Y(*)
      INTEGER K,L,I,J
      DOUBLE PRECISION X

      DO 1 I=1,K
      X=0.
      DO 2 J=1,L
    2 X=X+M(I,J)*V(J)
      Y(I)=X
    1 CONTINUE
      RETURN
      END
