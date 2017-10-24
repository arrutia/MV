! THESE ROUTINES ARE FOR A MOTION VECTOR COMPUTATION USING PYTHON-FORTRAN 
! COMPILE USING: f2py -c -m test /home/jruiz/PYTHON/motion_vectors.f90
! FOR PARALLEL FUNCTIONALITY f2py -c -lgomp --f90flags="-fopenmp -lgomp" -m motion_vectors motion_vectors_parallel.f90 
! FOR DEBUG  FUNCTIONALITY f2py -c -lgomp --f90flags="-g -traceback" -m motion_vectors motion_vectors_parallel.f90 

SUBROUTINE MOTION_VECTOR(field_t0,field_t1,dt,dx,box_size,sigma,nx,ny,desp_max,min_box_fraction,u_motion,v_motion,max_corr &
       &      ,trend,nref,aux_output,aux_inputi,aux_inputj,motion_vector_option,motion_vector_weigth,motion_vector_norm)

IMPLICIT NONE
INTEGER, PARAMETER :: sp = kind(1.0e0)
INTEGER, PARAMETER :: dp = kind(1.0d0)
INTEGER            :: ii , jj , iii , jjj , kkk , ini_i , ini_j , end_i , end_j , contador , n_points , n_points2
REAL(dp), PARAMETER :: UNDEF=-999.0d0 , min_dbz=0.0d0
INTEGER , PARAMETER :: METRIC_OPTION = 1 ! Metric option for similarity quantification (1 = correlation , ...)
REAL(dp), INTENT(IN) :: field_t0(nx,ny),field_t1(nx,ny),dt,dx,sigma,min_box_fraction
INTEGER , INTENT(IN) :: nx,ny,desp_max,box_size
INTEGER , INTENT(IN) :: aux_inputi , aux_inputj
REAL(dp), INTENT(OUT) :: u_motion(nx,ny) , v_motion(nx,ny) , max_corr(nx,ny) , trend(nx,ny)
INTEGER , INTENT(OUT) :: nref(nx,ny) !Number of points with ref > min_dbz within the box.
REAL(dp), INTENT(OUT) :: aux_output(2*desp_max+ 1,2*desp_max+1)
INTEGER  :: tmp_dim(1)
REAL(dp),ALLOCATABLE :: corr_matrix(:,:) , vector_1(:) , vector_2(:) 
REAL(dp),ALLOCATABLE :: correlation_weigths(:,:) , vector_cw(:)   !Gaussian weights for weighted correlation.
REAL(dp) :: box_fraction , dist , tmp_corr  ,  despi ,despj,weigth
INTEGER  :: maxi , maxj  , i_point
REAL(dp) , PARAMETER :: correlation_threshold = 0.3  ! $
INTEGER , INTENT(IN) :: MOTION_VECTOR_OPTION  !0 - Maximum corr , 1 Average sourrounding max corr, 2 - Average over entire displacements.
INTEGER , INTENT(IN) :: MOTION_VECTOR_WEIGTH  !0 - Uniform weigth , 1 - Gaussian weigths.
INTEGER , INTENT(IN) :: MOTION_VECTOR_NORM    !0 - Linear correlation, 1 CDF Correlation , 2 -MSE
INTEGER , PARAMETER :: NT=6
REAL(dp)  :: THRESHOLDS(NT)
REAL(dp), ALLOCATABLE  :: field2_t0(:,:) , field2_t1(:,:)
INTEGER  :: nx2 , ny2

! Zero padding 
nx2= nx + 2*( desp_max + box_size )
ny2= ny + 2*( desp_max + box_size )

ALLOCATE( field2_t0(nx2,ny2) , field2_t1(nx2,ny2) )

field2_t0=min_dbz
field2_t1=min_dbz

WRITE(*,*)field_t0(1,1)
WRITE(*,*)field_t0(nx,1)
WRITE(*,*)field_t0(1,ny)
WRITE(*,*)field_t0(nx,ny)

field2_t0(desp_max+box_size+1:desp_max+box_size+nx,desp_max+box_size+1:desp_max+box_size+ny)=field_t0
field2_t1(desp_max+box_size+1:desp_max+box_size+nx,desp_max+box_size+1:desp_max+box_size+ny)=field_t1

where( field2_t0 < min_dbz) 
   field2_t0=min_dbz
end where
where( field2_t1 < min_dbz)
   field2_t1=min_dbz
end where 


DO ii=1,NT
   THRESHOLDS(ii)=ii*10
ENDDO 

WRITE(*,*) " INSIDE W MOTION VECTOR ROUTINE "
WRITE(*,*) " FIRST REF FIELD AT T 0 = ",field_t0(1,1)
WRITE(*,*) " FIRST REF FIELD AT T 1 = ",field_t1(1,1)
WRITE(*,*) " DELTA T                  ",dt
WRITE(*,*) " DELTA X                  ",dx
WRITE(*,*) " SIGMA                    ",box_size
WRITE(*,*) " SIGMA_THRESHOLD          ",sigma
WRITE(*,*) " NX ,  NY                 ",nx,ny
WRITE(*,*) " DESP MAX                 ",desp_max
WRITE(*,*) " MIN_BOX_FRACTION         ",min_box_fraction
WRITE(*,*) " THRESHOLDS               ",THRESHOLDS

! First find the box size based on sigma and the threshold. In this implementation
! the box size is computed from sigma and sigma_threshold.

ALLOCATE(corr_matrix(desp_max*2+1,desp_max*2+1) , vector_1((box_size*2+1)**2) , vector_2((box_size*2+1)**2)  )

! Compute weigths for correlation.
ALLOCATE( correlation_weigths(box_size*2+1,box_size*2+1) , vector_cw((box_size*2+1)**2) )


SELECT CASE ( MOTION_VECTOR_WEIGTH )

CASE( 0 ) !Uniform weigths

   correlation_weigths=1.0d0

CASE( 1 ) !Gaussian weigths

   DO ii=-box_size,box_size
    DO jj=-box_size,box_size
      dist=( (ii*dx) ** 2 + (jj*dx) ** 2 )
      correlation_weigths( ii+box_size+1 , jj+box_size+1 ) = exp (-( dist / (sigma**2) ) )
    ENDDO
  ENDDO

END SELECT

!Initialize u_motion and v_motion with UNDEF values.
u_motion=undef
v_motion=undef
max_corr=-1.0d0
nref=0

n_points=(box_size * 2 + 1)**2
tmp_dim(1)=n_points

!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(ii,jj,ini_i,end_i,ini_j,end_j,iii,jjj,kkk,corr_matrix & 
!$OMP&           ,vector_1,vector_2,vector_cw,contador,box_fraction      &
!$OMP&           ,despi,despj,tmp_corr,weigth,maxi,maxj,n_points2) 
!$OMP DO 


DO ii = desp_max+box_size+1,desp_max+box_size+nx

  WRITE(*,*)"Processing row number ",ii
  DO jj = desp_max+box_size+1,desp_max+box_size+ny
 
        corr_matrix=-1.0d0

        ini_i=ii-box_size
        end_i=ii+box_size
        ini_j=jj-box_size
        end_j=jj+box_size

        vector_1=RESHAPE(field2_t0(ini_i:end_i,ini_j:end_j),tmp_dim)
        vector_cw=RESHAPE(correlation_weigths,tmp_dim)
        vector_cw=vector_cw/SUM(vector_cw)

        nref(ii-box_size-desp_max,jj-box_size-desp_max)=0.0d0
        DO iii=1,n_points
           if( vector_1(iii) > min_dbz)then
             nref(ii-box_size-desp_max,jj-box_size-desp_max)=nref(ii-box_size-desp_max,jj-box_size-desp_max)+1
           endif
        ENDDO

        if( nref(ii-box_size-desp_max,jj-box_size-desp_max) < min_box_fraction )CYCLE

        maxi=0
        maxj=0


        DO iii= -desp_max,desp_max
          DO jjj= -desp_max,desp_max

            vector_2=RESHAPE(field2_t1(ini_i+iii:end_i+iii,ini_j+jjj:end_j+jjj),tmp_dim)

            SELECT CASE( MOTION_VECTOR_NORM )
 
            CASE( 0 ) !Weigthed linear correlation. 
             CALL CORRELATION( vector_1,vector_2,vector_cw,n_points,tmp_corr)
            CASE( 1 ) !Weigthed CDF correlation
             CALL CDFCORRELATION(vector_1,vector_2,vector_cw,n_points,THRESHOLDS,NT,tmp_corr)
            CASE( 2 ) !Weighted MSE
             CALL MSE(vector_1,vector_2,vector_cw,n_points,tmp_corr)
             IF( tmp_corr .NE. UNDEF  )THEN
                tmp_corr = 1.0d0/(tmp_corr+1.0d0) !The index has to be possitively oriented.
             ENDIF
            END SELECT

            corr_matrix(iii+desp_max+1,jjj+desp_max+1)=tmp_corr
            IF( tmp_corr >= max_corr(ii-box_size-desp_max,jj-box_size-desp_max) )THEN
                max_corr(ii-box_size-desp_max,jj-box_size-desp_max)=tmp_corr
                maxj=jjj
                maxi=iii
                !Compute the mean reflectivity trend in the box.
                trend(ii-box_size-desp_max,jj-box_size-desp_max)=SUM( vector_2-vector_1, &
      &         vector_1>min_dbz .OR. vector_2 > min_dbz)/SUM(vector_1-vector_1+1.0d0,vector_1>min_dbz &
      &         .OR. vector_2 > min_dbz)
            ENDIF

          ENDDO
        ENDDO
 
     !!!!! FOR DEBUG
     IF( ii==aux_inputi+desp_max+box_size .AND. jj==aux_inputj+desp_max+box_size )THEN
         aux_output=corr_matrix
         
         !WRITE(*,*)maxj,maxi,max_corr(ii-box_size-desp_max,jj-box_size-desp_max)
         !vector_2=RESHAPE(field2_t1(ini_i+0:end_i+0,ini_j+1:end_j+1),tmp_dim)
         !DO iii=1,tmp_dim(1)
         !   WRITE(*,*)vector_2(iii),vector_1(iii)
         !ENDDO
       
     ENDIF 
     !!!!!!! $

     SELECT CASE ( MOTION_VECTOR_OPTION )

     CASE( 0 ) !Maximum correlation criteria. 

       despi=REAL(maxi,dp)
       despj=REAL(maxj,dp)
 
       IF( max_corr(ii-box_size-desp_max,jj-box_size-desp_max) > 0.0d0 )THEN
       u_motion(ii-box_size-desp_max,jj-box_size-desp_max)= despj*dx/dt
       v_motion(ii-box_size-desp_max,jj-box_size-desp_max)= despi*dx/dt
       ENDIF

     CASE( 1 ) !Mean correlation around the maximum value
         !Find the weighted location for the minimum.
         despi=0.0d0
         despj=0.0d0
         weigth=0.0d0
         DO iii=maxi-1,maxi+1
           DO jjj=maxj-1,maxj+1 
             IF( iii <= desp_max  .AND. iii >= -desp_max .AND. jjj <= desp_max .AND. jjj >= -desp_max )THEN
                IF( corr_matrix( iii+desp_max+1 , jjj+desp_max + 1) > 0)THEN
                    tmp_corr = corr_matrix( iii+desp_max+1 , jjj+desp_max + 1)**2
                    weigth= weigth + tmp_corr
                    despi = despi  + tmp_corr *iii
                    despj = despj  + tmp_corr *jjj
                ENDIF
             ENDIF                
           ENDDO
         ENDDO

         IF( weigth > 0.0)THEN
           u_motion(ii-box_size-desp_max,jj-box_size-desp_max)= despj*dx/(dt*weigth)
           v_motion(ii-box_size-desp_max,jj-box_size-desp_max)= despi*dx/(dt*weigth)
         ENDIF

     CASE( 2 ) !Mean correlation using the highest correlation points.
        weigth=0.0d0
        despi=0.0d0
        despj=0.0d0
        DO iii= -desp_max,desp_max 
           DO jjj= -desp_max,desp_max
             tmp_corr = corr_matrix( iii + desp_max + 1, jjj + desp_max + 1 )
             IF( max_corr(ii-box_size-desp_max,jj-box_size-desp_max) - tmp_corr < correlation_threshold .AND. tmp_corr > 0.0 )THEN
                 tmp_corr = tmp_corr ** 2
                 weigth=weigth+tmp_corr
                 despi=despi + iii * (tmp_corr)
                 despj=despj + jjj * (tmp_corr)
             ENDIF
           ENDDO
        ENDDO

         !Get the location of the maximum correlation within
         !corr_matrix
         IF( weigth > 0.0)THEN
           u_motion(ii-box_size-desp_max,jj-box_size-desp_max)= despj*dx/(dt*weigth)
           v_motion(ii-box_size-desp_max,jj-box_size-desp_max)= despi*dx/(dt*weigth)
         ENDIF

     END SELECT

  ENDDO
ENDDO

!$END OMP DO
!$OMP END PARALLEL

DEALLOCATE(corr_matrix , vector_1 , vector_2  )
DEALLOCATE( correlation_weigths , vector_cw )
DEALLOCATE( field2_t0 , field2_t1 )

WRITE(*,*)"FINISH MOTION VECTOR COMPUTATION"
END SUBROUTINE MOTION_VECTOR



SUBROUTINE CDFCORRELATION( A,B,W,L,THRESHOLDS,NT,CDFINDEX)
!This subroutine computes a correspondance between two vectors
!based on the comparison of the correspondent CDF of the two vectors.
!The idea is taking the probability of detection (POD) for different 
!thresholds and output a weigthted sum of the PODs corresponding to 
!different thresholds.
! CDFINDEX= (N1xPOD1 + ... NntxPODnt)/( N1 + ..... Nnt)
!Where nt is the number of thresholds (input), N is the number of points above
!the threshold in A and POD1 is the number of points which is above the threshold in A and B simultaneously.

INTEGER, PARAMETER :: sp = kind(1.0e0)
INTEGER, PARAMETER :: dp = kind(1.0d0)
INTEGER , INTENT(IN) :: L , NT  !Number of elements in inputs vectors and number of thresholds.
REAL(dp), INTENT(IN) :: A(L) , B(L) , THRESHOLDS(NT)  !Input vectors A is the reference vector.
REAL(dp), INTENT(IN) :: W(L) !Weigths.
REAL(dp), PARAMETER  :: UNDEF=-999.0d0
REAL(dp)             :: C(L) !Auxiliar arrays
REAL(dp), INTENT(OUT) :: CDFINDEX   !Output
INTEGER   :: ii,jj,N

    CDFINDEX=0.0d0
    N=0

DO ii=1,NT
       C=0.0d0
    WHERE( A >= THRESHOLDS(ii) .AND. B >= THRESHOLDS(ii) )
       C=1.0d0 * W
    ENDWHERE
    CDFINDEX= CDFINDEX + SUM(  C  )
       C=0.0d0
    WHERE( A >= THRESHOLDS(ii) )
       C=W
    ENDWHERE
    N=N + SUM( C )
ENDDO

!  IF( N .GT. 0)THEN
!  CDFINDEX=CDFINDEX/REAL(N,dp)
!  ELSE
!  CDFINDEX=0.0d0
!  ENDIF

END SUBROUTINE CDFCORRELATION

SUBROUTINE MSE( A,B,W,L,WMSE)
!This subroutine computes a weigthed MSE between two vectors.

INTEGER, PARAMETER :: sp = kind(1.0e0)
INTEGER, PARAMETER :: dp = kind(1.0d0)
INTEGER , INTENT(IN) :: L  !Number of elements in inputs vectors 
REAL(dp), INTENT(IN) :: A(L) , B(L)  !Input vectors A is the reference vector.
REAL(dp), INTENT(IN) :: W(L) !Weigths.
REAL(dp), PARAMETER  :: UNDEF=-999.0d0
REAL(dp)             :: TW       !Auxiliar arrays
REAL(dp), INTENT(OUT) :: WMSE   !Output
INTEGER   :: ii,jj,N

    WMSE=0.0d0
    TW=0.0d0

DO ii=1,L
    IF( A(ii) .NE. UNDEF .AND. B(ii) .NE. UNDEF )THEN
       WMSE = WMSE + W(ii)*( ( A(ii) - B(ii) )**2 )
       TW = TW + W(ii)
    ENDIF
ENDDO

IF( TW .GT. 0.0d0 )THEN
  WMSE= WMSE / TW
ELSE
  WMSE= UNDEF
ENDIF


END SUBROUTINE MSE

SUBROUTINE UNDEF_MEAN( field,nx,ny, mean)
INTEGER, PARAMETER :: sp = kind(1.0e0)
INTEGER, PARAMETER :: dp = kind(1.0d0)
INTEGER , INTENT(IN) :: nx,ny
REAL(dp), INTENT(IN) :: field(nx,ny)
REAL(dp), PARAMETER  :: UNDEF=-999.0d0
REAL(dp), INTENT(OUT) :: mean
INTEGER   :: ii,jj,n

mean=0.0d0
n=0

DO ii=1,nx
  DO jj=1,ny
    IF( field(ii,jj) .NE. UNDEF)THEN
      mean = mean + field(ii,jj)    
      n=n+1
    ENDIF
  ENDDO
ENDDO

IF( n > 0)THEN
  mean = mean / REAL(n,dp)

ELSE
 
  mean = UNDEF

endif

END SUBROUTINE UNDEF_MEAN


!-----------------------------------------------------------------------
! Correlation
!-----------------------------------------------------------------------
SUBROUTINE CORRELATION(A,B,W,N,COR)
IMPLICIT NONE
INTEGER, PARAMETER :: sp = kind(1.0e0)
INTEGER, PARAMETER :: dp = kind(1.0d0)


INTEGER, INTENT(IN) :: N
REAL(dp) , INTENT(IN) :: A(N),B(N),W(N)
REAL(dp) , INTENT(OUT) :: COR

!Weighted correlation, A and B are two variables and W are the weights that 
!corresponds to each element of A and B.

REAL(dp) :: STDA , STDB , COV

  CALL STDEV(A,W,N,STDA)
  CALL STDEV(B,W,N,STDB)
  CALL COVAR(A,B,W,N,COV)

  COR = COV/STDA/STDB

END SUBROUTINE CORRELATION

!-----------------------------------------------------------------------
! Covariance
!-----------------------------------------------------------------------
SUBROUTINE COVAR(A,B,W,N,COV)
  IMPLICIT NONE
  INTEGER, PARAMETER :: sp = kind(1.0e0)
  INTEGER, PARAMETER :: dp = kind(1.0d0)
  INTEGER,INTENT(IN) :: N
  REAL(dp),INTENT(IN) :: A(N)
  REAL(dp),INTENT(IN) :: B(N)
  REAL(dp),INTENT(IN) :: W(N)
  REAL(dp),INTENT(OUT) :: COV

  REAL(dp) :: MEANA,MEANB

  MEANA=SUM(W*A)/SUM(W)  !Weighted mean.
  MEANB=SUM(W*B)/SUM(W)  !Weighted mean

  COV = SUM( W*(A-MEANA)*(B-MEANB) )/SUM(W) !SUM( W*(A-MEANA)*(B-MEANB) )  !Weigthed covariance

  RETURN
END SUBROUTINE COVAR

!-----------------------------------------------------------------------
! Standard deviation
!-----------------------------------------------------------------------
SUBROUTINE STDEV(A,W,N,STD)
  IMPLICIT NONE
  INTEGER, PARAMETER :: sp = kind(1.0e0)
  INTEGER, PARAMETER :: dp = kind(1.0d0)
  INTEGER,INTENT(IN) :: N
  REAL(dp),INTENT(IN) :: A(N)
  REAL(dp),INTENT(IN) :: W(N)
  REAL(dp),INTENT(OUT) :: STD

  !STD = SQRT( SUM ((A - SUM(A)/REAL(N,dp))**2)/REAL(N,dp) )

  STD = SQRT( SUM(W*(A(:) - SUM(W*A)/SUM(W) )**2)/SUM(W) )

  RETURN
END SUBROUTINE STDEV


SUBROUTINE GAUSSIAN_FILTER(field0,field1,dx,sigma,nx,ny)
IMPLICIT NONE
INTEGER, PARAMETER :: sp = kind(1.0e0)
INTEGER, PARAMETER :: dp = kind(1.0d0)
INTEGER            :: ii , jj , iii , jjj 
REAL(dp), PARAMETER :: UNDEF=-999.0d0 
REAL(dp), INTENT(IN) :: field0(nx,ny),dx,sigma
REAL(dp), INTENT(OUT):: field1(nx,ny)
INTEGER , INTENT(IN) :: nx,ny
REAL(dp)             :: tmp , suma , sumaw , dist

WRITE(*,*) " INSIDE GAUSSIAN FILTER "
WRITE(*,*) " DELTA X                  ",dx
WRITE(*,*) " SIGMA                    ",sigma
WRITE(*,*) " NX ,  NY                 ",nx,ny

!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(ii,jj,iii,jjj,dist,tmp,suma,sumaw) 
!$OMP DO 
DO ii = 1,nx
  DO jj = 1,ny
        !COMPUTE AVERAGED VALUE.
        suma=0.0d0
        sumaw=0.0d0
        DO iii=1,nx
          DO jjj=1,ny
               dist=( ((ii-iii)*dx) ** 2 + ((jj-jjj)*dx) ** 2 )/( sigma**2)
               IF( dist <= 10)THEN
                 tmp=exp(-dist)
                 suma=suma+field0(iii,jjj)*tmp
                 sumaw=sumaw+tmp
               ENDIF
          ENDDO
        ENDDO

  IF( sumaw > 0.0d0)THEN
   field1(ii,jj)=suma/sumaw
  ELSE
   field1(ii,jj)=UNDEF
  ENDIF

  ENDDO
ENDDO

!$END OMP DO
!$OMP END PARALLEL


WRITE(*,*)"FINISH GAUSSIAN FILTER"
END SUBROUTINE GAUSSIAN_FILTER



SUBROUTINE FILTER_OUTLIER(field0,field1,nx,ny,threshold,box_size)
IMPLICIT NONE
INTEGER, PARAMETER :: sp = kind(1.0e0)
INTEGER, PARAMETER :: dp = kind(1.0d0)
INTEGER            :: ii , jj , iii , jjj 
REAL(dp), PARAMETER :: UNDEF=-999.0d0 
REAL(dp), INTENT(IN) :: field0(nx,ny),threshold
REAL(dp), INTENT(OUT):: field1(nx,ny)
INTEGER , INTENT(IN) :: nx,ny,box_size
REAL(dp)             :: suma,sumasq,sumaw
INTEGER              :: ini_i,end_i,ini_j,end_j,np

field1=field0

!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(ii,jj,iii,jjj,suma,sumasq,sumaw, &
!$OMP&                          ini_i,end_i,ini_j,end_j,np) 
!$OMP DO 
!DO ii = 1,nx
!  DO jj = 1,ny
        !COMPUTE LOCAL AVERAGE
        ini_i=ii-box_size
        end_i=ii+box_size
        ini_j=jj-box_size
        end_j=jj+box_size

        if(ini_i < 1)ini_i=1
        if(end_i > nx)end_i=nx
        if(ini_j < 1)ini_j=1
        if(end_j > ny)end_j=ny

        !Compute sample size
        np=(end_i-ini_i+1)*(end_j-ini_j+1)
        
        suma=0.0d0
        sumaw=0.0d0
        do iii=ini_i,end_i
          do jjj=ini_j,end_j
            if( field0(iii,jjj) .ne. UNDEF )then
            suma=suma+field0(iii,jjj)
            sumaw=sumaw+1.0d0
            endif
          enddo
        enddo
        if( sumaw == 0.0d0 )cycle

        suma=suma/real(sumaw)

        if( ABS(field0(ii,jj) - suma ) > threshold )then
            field1(ii,jj)=UNDEF
        endif


  ENDDO
ENDDO

!$END OMP DO
!$OMP END PARALLEL


WRITE(*,*)"FINISH FILTER OUTLIER"
END SUBROUTINE FILTER_OUTLIER










