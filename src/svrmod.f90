module svrmod

  implicit none

  ! SVR Penalization Parameter
  real(8), parameter :: C_SVR = 1.0D+10
  ! SVR Tolerance
  real(8), parameter :: EPS_SVR = 1.25D-15

  ! COMMON ARRAYS

  real(8), allocatable :: QQ_(:, :), v_(:)

  ! COMMON SCALARS

  logical :: OUTPUT = .true.

  private

  ! Public variables and subroutines

  public :: qsvm, finalize, svrToTRDF, fromScratch, svr_set_output

contains

  !----------------------------------------------------------!
  ! SUBROUTINE SVR_SET_OUTPUT                                !
  !                                                          !
  ! Sets the output ON or OFF.                               !
  !----------------------------------------------------------!

  subroutine svr_set_output(o)

    ! SCALAR ARGUMENTS
    logical, intent(in) :: o

    OUTPUT = o

  end subroutine svr_set_output

  !----------------------------------------------------------!
  ! SUBROUTINE FINALIZE                                      !
  !                                                          !
  ! Destroys the allocated common structure.                 !
  !----------------------------------------------------------!

  subroutine finalize(status)

    implicit none

    ! SCALAR ARGUMENTS

    integer, intent(out) :: status

    deallocate(QQ_, v_, STAT = status)

  end subroutine finalize

  !----------------------------------------------------------!
  ! SUBROUTINE INITIALIZE                                    !
  !                                                          !
  ! Initializes the global structure of the module, if not   !
  ! initialized. Returns the status of the initialization in !
  ! 'status'.                                                !
  !----------------------------------------------------------!

  subroutine initialize(n, npt, status)

    implicit none

    ! SCALAR ARGUMENTS

    integer :: n, npt, status

    intent(in ) :: n, npt
    intent(out) :: status

    ! Allocate QQ_

    status = 0

    if ( allocated(QQ_) .and. size(QQ_) .ne. (2 * npt) * (2 * npt) ) then

       deallocate(QQ_)

    end if

    if ( .not. allocated(QQ_)  ) &
         allocate(QQ_(2 * npt, 2 * npt), STAT = status)

    if ( status .ne. 0 ) return

    ! Allocate v_
       
    if ( allocated(v_) .and. size(v_) .ne. (2 * npt) ) then

       deallocate(v_)

    end if

    if ( .not. allocated(v_)  ) &
         allocate(v_(2 * npt), STAT = status)

  end subroutine initialize

  !----------------------------------------------------------!
  ! SUBROUTINE FROMSCRATCH                                   !
  !                                                          !
  ! Builds a well poised sampling set from scratch.          !
  !                                                          !
  !----------------------------------------------------------!

  SUBROUTINE fromScratch(N, X, NPT, DELTA, OBJF, H, Y, FF, FLAG)

    implicit none

    ! SCALAR ARGUMENTS
    integer :: n, npt, flag
    real(8) :: delta

    ! ARRAY ARGUMENTS
    real(8) :: FF(npt), H(npt + n + 1, npt + n + 1), x(n), &
               Y(npt, n)

    intent(in ) :: delta, n, npt, x
    intent(out) :: ff, flag, h, y

    ! EXTERNAL SUBROUTINES
    external :: OBJF

    ! LOCAL ARRAYS

    integer :: IP(NPT), IQ(NPT)
    real(8) :: E(N + 1,NPT), YY(N), OMEGA(NPT,NPT), &
         GAMA(N + 1,N + 1), Z(NPT,NPT-N-1)

    ! LOCAL SCALARS

    integer :: ACUMULADOR, i, ii, j, k

    ! NPT2N = 2*N+1
    DO I=1, N
       Y(1,I) = X(I) 
    END DO

    ! PRINT*, "Y(1,2) = ", Y(1,2)

    DO I = 1,N
       DO J = 1,N 
          Y(I + 1,    J) = X(J)
          Y(I + 1,    I) = X(I) + DELTA 
          Y(I + N + 1,J) = X(J)
          Y(I + N + 1,I) = X(I) - DELTA                 
       END DO
    END DO   ! Y até 2n + 1

    !  PRINT*, "Y(1,2) = ", Y(1,2)

    DO I = 1,2 * N + 1 
       DO J = 1, N
          YY(J) = Y(I,J) 
       END DO
       CALL OBJF(N, YY, FF(I), FLAG)
       !     PRINT*, FLAG
       IF ( FLAG .NE. 0 ) RETURN
    END DO

    ! NPT >= 2N+1       

    IF (NPT .GT. 2 * N + 1) THEN             
       ! SETTING THE POINTS M-2N+1
       !     IF (NPT .GT. 2*N+1) THEN
       DO J = 2 * N + 2, NPT
          IF ( J .LE. 3 * N + 1 )  IP(J) = J - 2 * N - 1
          IF ( J .GE. 3 * N + 2 )  IP(J) = IP(J - N) 
       END DO
       !     END IF 
       II = 1 
       DO J= 2 * N + 2, NPT 
          IF ( IP(J) + II .LE. N ) THEN
             IQ(J) = IP(J) + II
          ELSE 
             IQ(J) = IP(J) + II - N 
          END IF
          IF ( MOD(IP(J) ,N) .EQ. 0 ) II = II + 1 
       END DO

       !    OBTAIN THE POINTS Y OF 2N+1 TO NPT.
       DO I = 2 * N + 2, NPT
          DO J = 1, N 
             Y(I,J) = Y(IP(I) + 1, J) + Y(IQ(I) + 1, J) - Y(1,J)                
          END DO
       END DO

       DO I = 2 * N + 2, NPT
          DO J = 1, N
             YY(J) = Y(I,J) 
          END DO
          CALL OBJF(N, YY, FF(I), FLAG)   
          IF ( FLAG .NE. 0 ) RETURN
       END DO
    END IF

    !******************* FIRST INVERSE***************************

    DO I=1, NPT
       DO J=1, NPT
          OMEGA(I,J)=0D0
       END DO
    END DO

    DO I=1, N+1
       DO J=1, N+1
          GAMA(I,J)=0D0
       END DO
    END DO

    DO I=1, NPT  
       DO J=1, NPT-N -1
          Z(I,J)=0D0
       END DO
    END DO

    ! MATRIX E( N+1 X NPT)
    DO I=1, N+1
       DO J=1,NPT  
          E(I,J)=0.0D0
       END DO
    END DO

    E(1,1) = 1.0D0
    DO I=2, N+1
       E(I,I)= 1.0D0/(2.0D0*DELTA)
       E(I,I+N)= -1.0D0/(2.0D0*DELTA)
    END DO

    ! MATRIX Z(NPT X NPT-N-1)            
    DO I=1, N 
       Z(1,I)= -SQRT(2.0D0)/(DELTA**2)
       Z(I+1,I)= SQRT(2.0D0)/(2.0D0*DELTA**2)
       Z(N+I+1,I)= SQRT(2.0D0)/(2.0D0*DELTA**2)
    END DO

    ! THE NEW INVERSE FOR MORE OF 2N+1 POINTS IN Y
    IF (NPT .GT. 2*N+1) THEN             
       DO I=N+1, NPT-N-1 
          Z(1,I)= 1.0D0/(DELTA**2)
          Z(N+I+1, I) = Z(1,I) 
          Z(IP(N+I+1)+1,I) = -1.0D0/(DELTA**2)
          Z(IQ(N+I+1)+1,I) = -1.0D0/(DELTA**2)        
       END DO
    END IF

    ! MULTIPLYING ZZ^T FOR DETERMINE OMEGA          
    ACUMULADOR=0D0
    DO I=1, NPT   
       DO K=1, NPT              
          DO J=1, NPT-N-1
             ACUMULADOR =  ACUMULADOR + Z(I,J)*Z(K,J)
          END DO
          OMEGA(I,K) = ACUMULADOR               
          ACUMULADOR = 0D0              
       END DO
    END DO

    ! THE MATRIX INVERSE H     

    do i = 1, NPT + N + 1

       do j = 1, NPT + N + 1

          H(I,J) = 0.0D0

       end do

    end do
    
    DO I=1, NPT
       DO J=1, NPT
          H(I,J)=OMEGA(I,J)         
       END DO
    END DO

    ! THE N+1 LINES OF H                 
    DO I=NPT+1, NPT+N+1
       DO J= 1, NPT 
          H(I,J) = E(I-NPT,J)         
       END DO
    END DO

    ! THE N+1 COLUMNS OF H 

    DO I=1, NPT
       DO J= NPT+1, NPT+N+1
          H(I,J) = H(J,I)
       END DO
    END DO

    !  PRINT*, "Y(1,1) = ", Y(1,1)
    !  PRINT*, "Y(1,2) = ", Y(1,2)

  END SUBROUTINE fromScratch

  !----------------------------------------------------------!
  ! SUBROUTINE QSVM                                          !
  !                                                          !
  !----------------------------------------------------------!

  subroutine qsvm(Y, FF, N, NPT, c, eps, HQ, g, b, flag)

    implicit none

    ! A partir dos pontos Y(NPT,N) construir o modelo resolvendo

    !               MIN   0.5*z'*Q*z + v'*z
    !     (1)       S.A.  A*z = 0
    !                     0 <= z <= C_SVR

    ! Em que:

    ! Q = [M*M'  -M*M'; -M*M'   M*M'], com M*M' = Y*Y' + (Y*Y').^2

    ! v = [-f(Y) + EPS_SVR*ones(NPT,1); f(Y) - EPS_SVR*ones(NPT,1)]

    ! A = [-ones(NPT,1) ones(NPT,1)]

    ! INPUT

    ! Y(NPT,N)        = Sample set, NPT points in R^n
    ! FF(NPT)         = f evaluated at the sample set
    ! C_SVR           = SVR C constant
    ! EPS_SVR         = SVR epsilon parameter

    ! OUTPUT

    ! [HQ, g, b] = the quadratic model is
    !            m(x) = 0.5*x'*HQ*x + g'*x + b

    ! SCALAR ARGUMENTS

    integer :: N, NPT, flag
    real(8) :: b, c, eps

    ! ARRAY ARGUMENTS

    real(8) :: g(N), Y(NPT,N), FF(NPT), HQ(N,N)

    intent(in   ) :: eps
    intent(out  ) :: flag
    intent(inout) :: c

    ! ALGENCAN VARIABLES

    logical :: checkder
    integer :: m, nvparam, ENE, inform, jcnnzmax, hnnzmax

    real(8) :: cnorm,efacc,efstain,eoacc,epsfeas,epsopt,eostain, &
         f,nlpsupn,snorm

    ! LOCAL SCALARS

    integer :: I, J, K, CONT1, CONT2, status
    real(8) :: AUX, XHX, gX, TOL, SOMAB, absalpha, absgamma, wTw

    ! LOCAL ARRAYS

    real(8) :: RES(2 * NPT), M_VALUES(NPT), ALFA(NPT), &
         GAMA(NPT), XIS(N), HXIS(N), BAUX1(NPT), BAUX2(NPT),          &
         AUXO(NPT), l(2 * NPT), lambda(1), u(2 * NPT)
    character(80) :: specfnm, outputfnm, vparam(10)
    logical :: coded(11), equatn(1), linear(1)

    ! EXTERNAL SUBROUTINES

!    external :: svrevalf,svrevalg,svrevalh,svrevalc,svrevaljac,svrevalhc, &
!         svrevalfc,svrevalgjac,svrevalgjacp,svrevalhl,svrevalhlp

    !     1 - PROBLEM BUILDING      

    !      DO I = 1, NPT
    !         DO J = 1, N
    !            PRINT*, "Y(",I,J,") = ", Y(I,J)
    !         END DO
    !      END DO

    !      DO I = 1, NPT
    !         PRINT*, "FF(",I,") = ", FF(I)
    !      END DO   

    flag = 0

    if ( OUTPUT ) write(*,FMT=1000)

    call initialize(n, npt, status)

    if ( status .ne. 0 ) then

       write(*,*) 'Memory problems when initializing SVR structure.'

       flag = 3

       return

    end if

    ! do i = 1, NPT
    !    write(*,*) (Y(i,j), j = 1,n), FF(i)
    ! end do

    ! MATRIX Q

    DO I = 1,NPT

       AUX = 0.0D0

       DO K = 1,N
          AUX = AUX + Y(I,K) * Y(I,K)
       END DO

       QQ_(I,I)             = AUX + AUX**2
       QQ_(I + NPT,I + NPT) =   QQ_(I,I)
       QQ_(I,I + NPT)       = - QQ_(I,I)
       QQ_(I + NPT,I)       = - QQ_(I,I)

       DO J = I + 1,NPT

          AUX = 0.0D0

          DO K = 1,N
             AUX = AUX + Y(J,K) * Y(I,K)
          END DO

          QQ_(I,J)             = AUX + AUX**2
          QQ_(J,I)             =   QQ_(I,J)
          QQ_(I + NPT,J + NPT) =   QQ_(I,J)
          QQ_(J + NPT,I + NPT) =   QQ_(J,I)
          QQ_(I,J + NPT)       = - QQ_(I,J)
          QQ_(J,I + NPT)       = - QQ_(J,I)
          QQ_(I + NPT,J)       = - QQ_(I,J)
          QQ_(J + NPT,I)       = - QQ_(J,I)

       END DO

    END DO

         ! DO I = 1, 2*NPT
         !    DO J = 1, 2*NPT
         !       PRINT*, "Q(",I,J,")=", QQ_(I,J)
         !    END DO
         ! END DO

    ! VECTOR v       

    DO I = 1, NPT
       v_(I)       = - FF(I) + eps
       v_(I + NPT) =   FF(I) + eps
    END DO

    ! do i = 1, 2 * NPT
    !    write(*,*) 'v(',I,')=', v_(I)
    ! end do

    ! 2 - SOLVING THE QUADRATIC PROBLEM (USING ALGENCAN)

    ! Set lower bounds, upper bounds, and initial guess

    ENE = 2 * NPT

    ! Constraints

    m = 1

    equatn(1) = .true.

    lambda(1) = 1.0d0

    linear(1) = .true.

    ! Coded subroutines

    coded( 1) = .true.  ! fsub
    coded( 2) = .true.  ! gsub
    coded( 3) = .true.  ! hsub
    coded( 4) = .true.  ! csub
    coded( 5) = .true.  ! jacsub
    coded( 6) = .true.  ! hcsub
    coded( 7) = .false. ! fcsub
    coded( 8) = .false. ! gjacsub
    coded( 9) = .false. ! gjacpsub
    coded(10) = .false. ! hlsub
    coded(11) = .false. ! hlpsub

    ! Upper bounds on the number of sparse-matrices non-null elements

    jcnnzmax = ENE

    hnnzmax  = ENE * ENE

    ! Checking derivatives

    checkder = .false.

    ! Parameters setting

    epsfeas   = 1.0d-06
    epsopt    = 1.0d-06

    efstain   = sqrt( epsfeas )
    ! Disable early stopping criterium
    eostain   = - 1.0D0

    efacc     = sqrt( epsfeas )
    eoacc     = sqrt( epsopt )

    outputfnm = ''
    specfnm   = ''

    nvparam = 2
    vparam(1) = 'SAFEMODE'
    vparam(2) = 'LINEAR-SYSTEMS-SOLVER-IN-ACCELERATION-PROCESS MA57'
!    vparam(1) = 'OBJECTIVE-AND-CONSTRAINTS-SCALING-AVOIDED'
!    vparam(1) = 'OUTER-ITERATIONS-LIMIT 500'
!    vparam(2) = 'ITERATIONS-OUTPUT-DETAIL 56'

010    DO I = 1, ENE
       l(I)   = 0.0D0
       u(I)   = c
    END DO

    do i = 1, NPT
       ! TODO: Test if the previous solution can be used
       RES(I) = 1.0D0
       RES(NPT + i) = -1.0D0
    end do

    call algencan(svrevalf,svrevalg,svrevalh,svrevalc,svrevaljac,svrevalhc, &
    svrevalfc,svrevalgjac,svrevalgjacp,svrevalhl,svrevalhlp,jcnnzmax,      &
    hnnzmax,epsfeas,epsopt,efstain,eostain,efacc,eoacc,outputfnm,     &
    specfnm,nvparam,vparam,ENE,RES,l,u,m,lambda,equatn,linear,coded,  &
    checkder,f,cnorm,snorm,nlpsupn,inform)

    if ( OUTPUT ) write(*, FMT=1001) f, cnorm, nlpsupn

    if ( inform .ne. 0 ) then

       flag = 1

       write(*,*) 'Error in the solver'

       return

    end if
    
    !      PRINT*, "ALGENCAN = ", f, cnorm, snorm, nlpsupn, inform
    !      pause

    ! 3 - MODEL BUILDING 

    !      DO I = 1, 2*NPT
    !         PRINT*, "RES(", I, ")", RES(I)
    !      END DO

    DO I = 1, NPT
       ALFA(I) = RES(I)
       GAMA(I) = RES(NPT + I)
    END DO

    !      DO I = 1, NPT
    !         PRINT*, "ALFA(", I, ")", ALFA(I)
    !      END DO

    !      DO I = 1, NPT
    !         PRINT*, "GAMA(", I, ")", GAMA(I)
    !      END DO      

    !      g = Y'*(ALFA - GAMA)

    DO I = 1, NPT
       AUXO(I) = ALFA(I) - GAMA(I)
    END DO

         ! DO I = 1, NPT
         !    PRINT*, AUXO(I)
         ! END DO  

    DO I = 1, N
       g(I) = 0.0D0

       DO J = 1, NPT
          g(I) = g(I) + Y(J,I) * AUXO(J)
       END DO
    END DO

    !      DO I = 1, NPT
    !         DO J = 1, N
    !            PRINT*, "Y(", I, J, ")", Y(I,J)
    !         END DO
    !      END DO

    !      DO I = 1, N
    !         PRINT*, "g(", I, ")", g(I)
    !      END DO  

    !      H = 2*SUM_{i=1}^NPT (ALFA(i) - GAMA(i))*Y(i,:)*Y(i,:)'

    DO I = 1, N
       DO J = 1, N
          HQ(I,J) = 0.0D0

          DO K = 1, NPT
             HQ(I,J) = HQ(I,J) + (ALFA(K) - GAMA(K))*Y(K,I)*Y(K,J)
          END DO

          HQ(I, J) = 2.0D0 * HQ(I, J)
       END DO
    END DO

    ! Just for debugging purposes, evaluates w^T * w

    if ( OUTPUT ) then

       wTw = 0.0D0

       do j = 1, NPT

          wTw = wTw + auxo(j) * dot_product(auxo, QQ_(1:NPT,j))
          
       end do

       write(*, FMT = 9000) wTw

    end if

    ! Calculates 'b'

    ! TODO: Evaluates only when necessary
    DO J = 1, NPT

       DO I = 1, N
          XIS(I) = Y(J,I)
       END DO

       DO I = 1, N
          HXIS(I) = 0.0D0
          DO K = 1, N
             HXIS(I) = HXIS(I) + HQ(I,K) * XIS(K)
          END DO
       END DO
       !        CALL matvet(HQ,XIS,N,N,HXIS)
       !         PRINT*, HXIS(1)
       !         PRINT*, HXIS(2)
       ! CALL mvv(XIS,HXIS,N,XHX)
       ! CALL mvv(g,XIS,N,gX)

       XHX = dot_product(HXIS, XIS)

       gX  = dot_product(XIS, g)

       M_VALUES(J) = 0.5D0 * XHX + gX

    END DO

    !      DO I = 1, NPT
    !         PRINT*, "M_VALUES(",I,")", M_VALUES(I)
    !      END DO

    TOL = 1.0D-3

    CONT1 = 0

    absalpha = 0.0

    DO I = 1, NPT
       absalpha = max(absalpha, abs(alfa(i)))

       IF (ALFA(I) .GT. TOL) THEN
          IF (ALFA(I) .LT. (C_SVR - TOL)) THEN
             CONT1 = CONT1 + 1
             BAUX1(CONT1) = FF(I) - M_VALUES(I) - EPS_SVR
          END IF
       END IF
    END DO

    !      PRINT*, "CONT1 = ", CONT1

    !      DO I = 1, CONT1
    !         PRINT*, "BAUX1", BAUX1(I)
    !      END DO

    CONT2 = 0

    absgamma = 0.0D0

    DO I = 1, NPT
       absgamma = max(absgamma, abs(gama(i)))

       IF (GAMA(I) .GT. TOL) THEN
          IF (GAMA(I) .LT. (C_SVR - TOL)) THEN
             CONT2 = CONT2 + 1            
             BAUX2(CONT2) = FF(I) - M_VALUES(I) + EPS_SVR
          END IF
       END IF
    END DO

    if ( OUTPUT ) write(*, FMT=1002) c, eps, CONT1, CONT2, absalpha, &
                                     absgamma

    !      DO I = 1, CONT2
    !         PRINT*, "BAUX2", BAUX2(I)
    !      END DO

    !      PRINT*, "CONT2 = ", CONT2

    SOMAB = 0.0D0

    DO I = 1, CONT1
       SOMAB = SOMAB + BAUX1(I)
    END DO

    DO I = 1, CONT2
       SOMAB = SOMAB + BAUX2(I)
    END DO

    IF ((CONT1 + CONT2) .EQ. 0 ) THEN

       if ( absalpha .gt. TOL .or. absgamma .gt. TOL ) then
          
          c = c * 10.0D0

          goto 010

       else

          b = FF(1)

       end if

    ELSE
       b = SOMAB / (CONT1 + CONT2)
    END IF

!    write(*,*) b

    ! NON-EXECUTABLE STATEMENTS

1000 FORMAT(/,5X,'BUILDING MODEL by SVR')
1001 FORMAT(5X,3X,'Solved optimization problem (ALGENCAN)',/, &
          5X,3X,3X,'Objective value:',33X,1PD12.5,/, &
          5X,3X,3X,'Feasibility:',37X,1PD12.5,/,&
          5X,3X,3X,'Gradient norm:',35X,1PD12.5)
1002 FORMAT(5X,3X,'Penalization:',39X,1PD12.5,/, &
          5X,3X,'Eps:',48X,1PD12.5,/, &
          5X,3X,"Number of non-binding alpha's:",24X,I10,/, &
          5X,3X,"Number of non-binding gamma's:",24X,I10,/, &
          5X,3X,"||alpha||_inf:",38X,1PD12.5,/,&
          5X,3X,"||gamma||_inf:",38X,1PD12.5)
9000 FORMAT(5X,3X,'w^Tw:',47X,1PD12.5)
          

    

  END SUBROUTINE qsvm

  !----------------------------------------------------------!
  ! SUBROUTINE SVMTOTRDF                                     !
  !                                                          !
  ! Converts data of a model built by SVR to data used by    !
  ! the TRDF algorithm.                                      !
  !----------------------------------------------------------!

  subroutine svrToTRDF(n, HQ, g, b, TRDF_XBASE, TRDF_HQ, TRDF_GOPT, &
       TRDF_MOD, Q)

    implicit none

    ! SCALAR ARGUMENTS
    
    integer :: n
    real(8) :: b, TRDF_MOD

    ! ARRAY ARGUMENTS

    real(8) :: HQ(N,N), g(N), Q(1 + N + (N + 1) * N / 2), TRDF_GOPT(N), &
         TRDF_HQ(N * (N + 1) / 2), TRDF_XBASE(N)

    intent(in ) :: n, b, g, HQ
    intent(out) :: Q, TRDF_GOPT, TRDF_HQ, TRDF_MOD, TRDF_XBASE

    ! LOCAL SCALARS
    
    integer :: i, j, k

    ! Evaluate the model at the basis and initializes TRDF's common
    ! variables

    do i = 1, n

       TRDF_GOPT(i) = 0.0D0

    end do

    do j = 1, n

       do i = 1, n
          
          TRDF_GOPT(i) = TRDF_GOPT(i) + HQ(i,j) * TRDF_XBASE(j)

       end do

    end do

    TRDF_MOD = dot_product(TRDF_GOPT, TRDF_XBASE)

    TRDF_MOD = TRDF_MOD / 2.0D0

    TRDF_MOD = TRDF_MOD + dot_product(g, TRDF_XBASE) + b

    do i = 1, n

       TRDF_GOPT(i) = TRDF_GOPT(i) + g(i)

    end do

    ! Initializes TRDF's structure of Q

    Q(1) = TRDF_MOD

    do i = 1, n

       Q(i + 1) = TRDF_GOPT(i)

    end do

    ! It is important to follow the correct structure and
    ! insert the lower triangular part of HQ

    K = 1

    DO I = 1, N
       DO J = 1, I
          Q(K + N + 1)       = HQ(I,J)
          TRDF_HQ(k)         = HQ(i,j)
          K                  = K + 1
       END DO
    END DO

    ! do i = 1, 1 + N + N * (N + 1) / 2
    !    write(*,*) 'Qc(',i,')=',Q(i)
    ! end do

  end subroutine svrToTRDF

  !----------------------------------------------------------!
  ! SUBROUTINE SVREVALF                                      !
  !                                                          !
  !----------------------------------------------------------!

  SUBROUTINE svrevalf(n,x,f,flag)

    IMPLICIT NONE

    ! SCALAR ARGUMENTS
    integer :: flag, n
    real(8) :: f

    ! ARRAY ARGUMENTS
    real(8) :: x(n)

    ! LOCAL ARRAYS
    real(8) :: Qx(n)

    ! LOCAL SCALARS
    integer :: I, K
    real(8) :: xQx, vx

    flag = 0

    DO I = 1, n
       Qx(I) = 0.0D0
       DO K = 1, n
          Qx(I) = Qx(I) + QQ_(I,K) * x(K)
       END DO
!       write(*,*) 'Q(',i,')=',Qx(i)
    END DO

    xQx = dot_product(x, Qx)

    vx  = dot_product(x, v_)

    f = xQx / 2.0D0 + vx

!    write(*,*) '---->',f, xQx, vx

  END SUBROUTINE svrevalf

  ! ******************************************************************
  ! ******************************************************************

  SUBROUTINE svrevalg(n,x,g,flag)

    IMPLICIT NONE

    ! SCALAR ARGUMENTS
    integer :: flag, n

    ! ARRAY ARGUMENTS
    real(8) :: x(n), g(n)

    ! LOCAL SCALARS
    integer :: I, K

    flag = 0

    DO I = 1, n
       g(I) = v_(I)

       DO K = 1, n
          g(I) = g(I) + QQ_(I,K) * x(K)
       END DO
    END DO

  END SUBROUTINE svrevalg

  ! ******************************************************************
  ! ******************************************************************

  SUBROUTINE svrevalh(n,x,hrow,hcol,hval,hnnz,lim,lmem,flag)

    IMPLICIT NONE

    ! SCALAR ARGUMENTS
    logical :: lmem
    integer :: flag, n, hnnz, lim

    ! ARRAY ARGUMENTS
    integer :: hcol(lim), hrow(lim)
    real(8) :: hval(lim), x(n)

    ! LOCAL SCALARS
    integer :: I, J

    flag = 0

    lmem = .false.

    hnnz = 0

    DO J = 1, n
       DO I = J, n
          hnnz       = hnnz + 1
          hrow(hnnz) = I
          hcol(hnnz) = J
          hval(hnnz) = QQ_(I,J)
       END DO
    END DO

  END SUBROUTINE svrevalh

  ! ******************************************************************
  ! ******************************************************************
  subroutine svrevalc(n,x,ind,c,flag)

    IMPLICIT NONE

    ! SCALAR ARGUMENTS
    integer, intent(in) :: ind,n
    integer, intent(out) :: flag
    real(kind=8), intent(out) :: c

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: x(n)

    ! LOCAL SCALARS
    integer :: i

    flag = 0

    if ( ind .eq. 1 ) then

       c = 0.0D0

       do i = 1, n / 2

          c = c - x(i) + x(n / 2 + i)

       end do

    else

       flag = - 1

    end if

  end subroutine svrevalc

  ! ******************************************************************
  ! ******************************************************************

  subroutine svrevaljac(n,x,ind,jcvar,jcval,jcnnz,lim,lmem,flag)

    implicit none

    ! SCALAR ARGUMENTS
    logical, intent(out) :: lmem
    integer, intent(in) :: ind,lim,n
    integer, intent(out) :: flag,jcnnz

    ! ARRAY ARGUMENTS
    integer, intent(out) :: jcvar(lim)
    real(kind=8), intent(in) :: x(n)
    real(kind=8), intent(out) :: jcval(lim)

    ! LOCAL
    integer :: I

    flag = 0

    lmem = .false.

    if ( ind .eq. 1 ) then

       if ( n .gt. lim ) then

          lmem = .true.

          return

       end if

       jcnnz = n

       DO I = 1, n / 2
          jcvar(I)         = I
          jcval(I)         = - 1.0D0

          jcvar(n / 2 + I) = n / 2 + I
          jcval(n / 2 + I) = 1.0D0
       END DO

    else

       flag = - 1

    end if

  end subroutine svrevaljac

  ! ******************************************************************
  ! ******************************************************************

  subroutine svrevalhc(n,x,ind,hcrow,hccol,hcval,hcnnz,lim,lmem,flag)

    implicit none

    ! SCALAR ARGUMENTS
    logical, intent(out) :: lmem
    integer, intent(in) :: ind,lim,n
    integer, intent(out) :: flag,hcnnz

    ! ARRAY ARGUMENTS
    integer, intent(out) :: hccol(lim),hcrow(lim)
    real(kind=8), intent(in) :: x(n)
    real(kind=8), intent(out) :: hcval(lim)

    flag = 0

    lmem = .false.

    hcnnz = 0

  end subroutine svrevalhc

  ! ******************************************************************
  ! ******************************************************************

  subroutine svrevalfc(n,x,f,m,c,flag)

    implicit none

    ! SCALAR ARGUMENTS
    integer,      intent(in)  :: m,n
    integer,      intent(out) :: flag
    real(kind=8), intent(out) :: f

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in)  :: x(n)
    real(kind=8), intent(out) :: c(m)

    flag = - 1

  end subroutine svrevalfc

  ! ******************************************************************
  ! ******************************************************************

  subroutine svrevalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,lim,lmem,flag)

    implicit none

    ! SCALAR ARGUMENTS
    logical,      intent(out) :: lmem
    integer,      intent(in)  :: lim,m,n
    integer,      intent(out) :: flag,jcnnz

    ! ARRAY ARGUMENTS
    integer,      intent(out) :: jcfun(lim),jcvar(lim)
    real(kind=8), intent(in)  :: x(n)
    real(kind=8), intent(out) :: g(n),jcval(lim)

    flag = - 1

  end subroutine svrevalgjac

  ! ******************************************************************
  ! ******************************************************************

  subroutine svrevalgjacp(n,x,g,m,p,q,work,gotj,flag)

    implicit none

    ! SCALAR ARGUMENTS
    logical,   intent(inout) :: gotj
    integer,   intent(in)    :: m,n
    integer,   intent(out)   :: flag
    character, intent(in)    :: work

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in)    :: x(n)
    real(kind=8), intent(inout) :: p(m),q(n)
    real(kind=8), intent(out)   :: g(n)

    flag = - 1

  end subroutine svrevalgjacp

  ! ******************************************************************
  ! ******************************************************************

  subroutine svrevalhl(n,x,m,lambda,sf,sc,hlrow,hlcol,hlval,hlnnz,lim,lmem,flag)

    implicit none

    ! SCALAR ARGUMENTS
    logical,      intent(out) :: lmem
    integer,      intent(in)  :: lim,m,n
    integer,      intent(out) :: flag,hlnnz
    real(kind=8), intent(in)  :: sf

    ! ARRAY ARGUMENTS
    integer,      intent(out) :: hlcol(lim),hlrow(lim)
    real(kind=8), intent(in)  :: lambda(m),sc(m),x(n)
    real(kind=8), intent(out) :: hlval(lim)

    flag = - 1

  end subroutine svrevalhl

  ! ******************************************************************
  ! ******************************************************************

  subroutine svrevalhlp(n,x,m,lambda,sf,sc,p,hp,goth,flag)

    implicit none

    ! SCALAR ARGUMENTS
    logical,      intent(inout) :: goth
    integer,      intent(in)    :: m,n
    integer,      intent(out)   :: flag
    real(kind=8), intent(in)    :: sf

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in)  :: lambda(m),p(n),sc(m),x(n)
    real(kind=8), intent(out) :: hp(n)

    flag = - 1

  end subroutine svrevalhlp

end module svrmod
