! topopatric.f90 -- v1
!Ayana B. Martins - 07/Fev/2017

MODULE globals
!!Loop variables
LOGICAL keep_going
INTEGER(4), ALLOCATABLE :: id(:)

!!CPU time variables
REAL:: start, finish
INTEGER(4) hours, minutes,seconds

!!Input variables
INTEGER(4) :: ntime,nc,mnpm,nb,radius,disp_rad
REAL rg
INTEGER(4) nf,deltat,sampled_times,window,stable
INTEGER(4) iREAD,inis,inisbit
REAL diff,mut,critG
LOGICAL independent_loci, is_dmi

!!Genotype and phenotype variables
INTEGER(1), ALLOCATABLE :: g(:,:)
INTEGER(1), ALLOCATABLE :: e(:,:)
INTEGER(1), ALLOCATABLE :: markers(:,:)

!!Fitness variables
REAL, ALLOCATABLE :: fit(:)
REAL(8) e0, sigk !trait optimum, variance in k

!!Space and neighbor variables
INTEGER(4), ALLOCATABLE :: x(:),y(:)
INTEGER(4), ALLOCATABLE :: neig(:,:),neigsp(:,:),ispecies(:,:)
INTEGER(4), ALLOCATABLE, SAVE :: worg(:,:,:),nworg(:,:)
INTEGER(4), ALLOCATABLE, SAVE :: idqx(:),idqy(:)
INTEGER(4) idensity,idq,nsites,this_radius
INTEGER(4), ALLOCATABLE, SAVE :: ineighbor(:), ineighborg(:)

!!Species detection variables
INTEGER(4) igt
INTEGER(4), ALLOCATABLE :: ispv(:)
INTEGER(4), ALLOCATABLE :: previous_igt(:)
REAL, ALLOCATABLE:: previous_gdists(:)


END MODULE

 MODULE Dist
IMPLICIT NONE
      REAL(8), PARAMETER :: PI=3.141592653589793238462
      CONTAINS

!
! Random Sample from normal (Gaussian) distribution
!
      FUNCTION normal(x,mean,stdev) RESULT(c)
       REAL(8) :: mean,stdev,c,x,expo,st
        expo = exp(-1*((x - mean)**2.0d0)/(2.0d0*(stdev**2.0d0)))
        c = expo/(sqrt(2.0d0*PI)*stdev)
      END FUNCTION     
      
END MODULE

PROGRAM topopatric
USE globals
USE Dist
IMPLICIT NONE
INTEGER(4) i,j,k,l,o !loop variables
INTEGER(4) nsets,iparameter
INTEGER(4) replica,ireplica
INTEGER(4) kmother,kmate,kcross
INTEGER(4) ix,iy
REAL dista
INTEGER(4) itime,iitime,iptime
REAL aux, random
CHARACTER*30 simulationID
CHARACTER*50 filename
CHARACTER*20 simID,repID
REAL check_time
INTEGER(4), ALLOCATABLE :: replica_igt(:), replica_stable(:)
REAL(8) ecotrait, fitsum

!!Random seed variables
INTEGER, DIMENSION(:), ALLOCATABLE :: iseed
INTEGER :: nseed

LOGICAL get_replacement
REAL isum

INTEGER(4) :: shuffle,shuffled
INTEGER(4), ALLOCATABLE, save :: random_order(:)
INTEGER(4), ALLOCATABLE:: auxx(:), auxy(:), auxg(:,:), auxm(:,:)

!!Temporarily stores offspring genotype
INTEGER(1), ALLOCATABLE :: goff(:,:)
INTEGER(1), ALLOCATABLE :: moff(:,:)
INTEGER(1), ALLOCATABLE :: eoff(:,:)
CALL cpu_time(start)
! READ input data
!!Brief description of each parameter can be found in the input file
OPEN (unit=7,file='input.in',status='old',position='rewind')
    READ(7,*) ntime,nc,nf,deltat
    READ(7,*) mut,diff,disp_rad,mnpm
    READ(7,*) radius,critG,nb, is_dmi
    READ(7,*) e0,sigk
    READ(7,*) iREAD,inis,inisbit
    READ(7,*) independent_loci,window,nsets,replica
CLOSE(7)
OPEN (unit=50,file='seed.in',status='old')  

!Absolute threshold is calculated from relative input
rg=nint(nb*critG)

IF (mnpm /= 0) THEN
	WRITE(6,*) 'Maybe you forgot that this version cannot handle the mnpm rule?'
	WRITE(6,*) 'Duh!'
	STOP
END IF

IF (nsets > 1) OPEN (unit=16,file='par.in',status='old',position='rewind')


!Erase appendable output files
OPEN (unit=17,file='summary.dat',status='unknown')
CLOSE(17, status='delete')

OPEN (unit=18,file='summary_replica.dat',status='unknown')
CLOSE(18, status='delete')

simulation: DO iparameter=1,nsets!change parameters

    ALLOCATE (replica_igt(replica),replica_stable(replica))
    IF (nsets > 1) READ (16, *) 
    IF (iparameter > 99) THEN
        WRITE(simID,'(i3)') iparameter
    ELSEIF (iparameter > 9) THEN
        WRITE(simID,'(i2)') iparameter
        simID = '0'//simID
    ELSE
        WRITE(simID,'(i1)') iparameter
        simID = '00'//simID
    END IF
    DO ireplica=1,replica

        IF (ireplica > 99) THEN
            WRITE(repID,'(i3)') ireplica
        ELSEIF (ireplica > 9) THEN
            WRITE(repID,'(i2)') ireplica
            repID = '0'//repID
        ELSE
            WRITE(repID,'(i1)') ireplica
            repID = '00'//repID
        END IF
        simulationID = 'SIM'//trim(simID)//'.'//trim(repID)
        WRITE (6,*) '-----------------',trim(simulationID),'----------------------'

        ! initialize random number generator
        CALL random_seed(size=nseed)
        ALLOCATE(iseed(nseed))
        READ(50,*) iseed
        CALL random_seed(put=iseed)
        CALL random_number(aux)

        idensity = 100 !maximum individuals per site
        !idensity is not intended to be a carrying capacity and is used just to set
        !!the lenght of the matrix
        
        nsites = 3.5*radius**2 !!number of sites in a neighborhood
        
        !Start a new simulation
        ALLOCATE (id(nc))
        ALLOCATE (g(nc,nb),markers(nc,nb),e(nc,nb))
        ALLOCATE (goff(nc,nb), eoff(nc,nb),moff(nc,nb))
        ALLOCATE (fit(nc))
        ALLOCATE (x(nc),y(nc))
        ALLOCATE (neig(nc,nc),neigsp(nc,nc))
        ALLOCATE (ineighbor(nc), ineighborg(nc))
        ALLOCATE (worg(nf,nf,idensity),nworg(nf,nf))
        ALLOCATE (idqx(nsites),idqy(nsites))
        ALLOCATE (ispv(nc),ispecies(nc,nc))
        ALLOCATE (previous_igt(ntime),previous_gdists(ntime))
        ALLOCATE (random_order(nc))
        ALLOCATE (auxx(nc),auxy(nc),auxg(nc,nb),auxm(nc,nb))
        
        !Set initial conditions of genotypes
        g = 0
        markers = 0
        goff = 0
        moff = 0
        e = 0
        eoff = 0
        fit = 0.0
        IF(inisbit /= 0) THEN !random genomes
            DO k=1,nc                    
                DO i=1,nb
                    CALL random_number(aux)
                    IF(aux > 0.5) g(k,i) = 1
                    !markers are also random
                    CALL random_number(aux)
                    IF(aux > 0.5) markers(k,i) = 1
                END DO
            END DO
        END IF

        !Set initial location of individuals
        DO i=1,nc
            IF(inis == 0) THEN  !random position
                CALL random_number(aux)
                x(i) = int(aux*nf)+1
                CALL random_number(aux)
                y(i) = int(aux*nf)+1
             ELSE !localized at the center
                x(i) = nf/2   
                y(i) = nf/2
             END IF
        END DO

        !Initialize individual's ids
        DO k=1,nc
            id(k) = k
        END DO
        
        !Save the information of in which site individuals are located
        nworg =0 !number of individuals per site
        worg = 0 !identity of individuals per site
        DO i=1,nc     
            nworg(x(i),y(i)) = nworg(x(i),y(i)) + 1
            IF (nworg(x(i),y(i)) > idensity) THEN
                WRITE(6,*) 'There too many individuals here! ln196'
            END IF
            worg(x(i),y(i),nworg(x(i),y(i))) = i
        END DO
        
        ! Set up template for the locations that are part of the neighborhood for a radius
        idqx = 0
        idqy = 0
        idq = 0
        DO i=-radius,radius               
            DO j=-radius,radius      
                IF(i*i+j*j <= radius**2) THEN !check if site belongs to the neighborhood
                    idq = idq + 1
                    idqx(idq) = i
                    idqy(idq) = j
                END IF
            END DO
         END DO
             
        !Erase appendable output files       
        filename = 'dist'//trim(simulationID)//'.dat'
        OPEN (unit=22,file=filename,status='unknown')
        CLOSE(22, status='delete')
        
        ! Initialize time variable
        iitime = 0

        ! Time evolution: mating, mutation and diffusion
        keep_going=.true.
        itime=0
        sampled_times=0
        looptime: DO WHILE (keep_going)
            itime=itime+1
            DO k = 1,nc
                CALL FINDNEIG(k,radius)
                CALL FINDGNEIG(k)
                fit(k) = 0.37
                ecotrait = 0
                DO l = 1,nb
                    ecotrait = ecotrait + e(k,l) 
                END DO
                ecotrait = ecotrait/nb
                fit(k) = normal(ecotrait,e0,sigk)
            END DO
            !Mating
            DO k = 1,nc
                get_replacement = .false. !initially, we assume that this individual is going to reproduce
                kmother = k
                !CALL FINDNEIG(kmother,radius)
                !CALL FINDGNEIG(kmother)
                !There is a random chance (m) that k is going to die and be replaced by the offspring
                !! of a pair of individuals nearby
                DO i = 1, ineighbor(k)
                    fitsum = fitsum + fit(neigsp(k,i))
                END DO
                CALL random_number(random) 
                IF (random < fit(k)) THEN
                    get_replacement = .true.
                END IF
                !Additionally, if k has not enough potential partners, it will also
                !!!!! need to be replaced
                IF(ineighborg(kmother) == 0) THEN
                    get_replacement = .true.
                END IF
                !If necessary, find a candidate to replace k
                IF (get_replacement) CALL FINDSUB(k,kmother)
                IF(ineighborg(kmother) == 0) THEN
                    WRITE(6,*) 'ineighborg = 0!'
                    WRITE(6,*) 'This should not be possible here (ln408)'
                    EXIT simulation
                END IF
                    CALL random_number(aux)
                    kmate = neig(kmother,int(aux*ineighborg(kmother))+1)    ! get a kmate /= kmother
                    IF(independent_loci) THEN
                        DO l=1,nb
                            CALL random_number(aux)
                            IF (aux<0.5) THEN
                                goff(k,l)=g(kmother,l)
                            ELSE
                                goff(k,l)=g(kmate,l)
                            END IF
                        END DO
                        DO l=1,nb
                            CALL random_number(aux)
                            IF (aux < 0.5) THEN
                                moff(k,l)=markers(kmother,l)
                            ELSE
                                moff(k,l)=markers(kmate,l)
                            END IF
                        END DO
                        DO l=1,nb
                            CALL random_number(aux)
                            IF (aux < 0.5) THEN
                                eoff(k,l)=e(kmother,l)
                            ELSE
                                eoff(k,l)=e(kmate,l)
                            END IF
                        END DO                        
                    ELSE
                        CALL random_number(aux)
                        kcross = int(aux*(nb-1))+1         ! crossover point
                        CALL random_number(aux)
                        IF(aux < 0.5) THEN
                            DO l=1,kcross          ! copy first kcross bits of g(kmate) 
                                goff(k,l) = g(kmate,l)  ! and last bits of g(kmother) into g(k)
                            END DO
                            DO l=kcross+1,nb        
                                goff(k,l) = g(kmother,l)
                            END DO
                        ELSE
                            DO l=1,kcross           ! copy first kcross bits of g(kmother) 
                                goff(k,l) = g(kmother,l) ! and last bits of g(kmate) into g(k)
                            END DO
                            DO l=kcross+1,nb       
                                goff(k,l) = g(kmate,l)
                            END DO
                        END IF
                        CALL random_number(aux)
                        kcross = int(aux*(nb-1))+1         ! crossover point
                        CALL random_number(aux)
                        IF(aux < 0.5) THEN
                            DO l=1,kcross          ! copy first kcross bits of g(kmate) 
                                eoff(k,l) = e(kmate,l)  ! and last bits of g(kmother) into g(k)
                            END DO
                            DO l=kcross+1,nb        
                                eoff(k,l) = e(kmother,l)
                            END DO
                        ELSE
                            DO l=1,kcross           ! copy first kcross bits of g(kmother) 
                                eoff(k,l) = e(kmother,l) ! and last bits of g(kmate) into g(k)
                            END DO
                            DO l=kcross+1,nb       
                                eoff(k,l) = e(kmate,l)
                            END DO
                        END IF                        
                    END IF
            END DO  ! end matings of the season
            
            ! Generations are discrete. Parents are replaced and offspring
            !! mutants and disperse only after the end of the mating season 
            DO k=1,nc
                    DO l=1,nb
                        g(k,l) = goff(k,l)
                        CALL random_number(aux)
                        IF(aux < mut) THEN
                            g(k,l) = 1-g(k,l)
                        END IF
                    END DO
                    DO l=1,nb
                        markers(k,l) = moff(k,l)
                        CALL random_number(aux)
                        IF(aux < mut) THEN
                            markers(k,l) = 1-markers(k,l)
                        END IF
                     END DO
                    DO l=1,nb
                        e(k,l) = eoff(k,l)
                        CALL random_number(aux)
                        IF(aux < mut) THEN
                            e(k,l) = 1-e(k,l)
                        END IF
                     END DO                     
                     IF(diff /= 0.0) CALL DISPERSAL(k)
             END DO           
   
            ! calculate species and other metrics every deltat
            check_time = (float(itime)/float(deltat) - itime/deltat)
            IF (check_time == 0.0) THEN
                sampled_times=sampled_times+1
                CALL FINDSPECIES         
                IF (igt < 5) THEN
                    WRITE(6,*) itime,itime+iitime,nc,igt,(ispv(k),k=1,igt)
                ELSE
                    WRITE(6,*) itime,itime+iitime,nc,igt,(ispv(k),k=1,4),'...'
                END IF
                filename = 'dist'//trim(simulationID)//'.dat'
                OPEN (unit=22,file=filename,status='unknown', position='append')
                    CALL is_stable
                CLOSE(22)
                IF(itime==ntime) THEN
                    keep_going=.false.
                    stable = 0
                END IF
            ELSE
                !WRITE(6,*) itime,itime+iitime
            END IF
            IF (keep_going) THEN
                !Shuffle indexes so that individuals reproduce in random order
                DO k=nc,1,-1
                    CALL random_number(aux)
                    shuffle = int(aux*(nc-1))+1
                    shuffled = id(k)
                    id(k) = id(shuffle)
                    id(shuffle) = shuffled
                END DO
                DO k=1,nc
                    auxx(k) = x(k)
                    auxy(k) = y(k)
                    DO l=1,nb
                        auxg(k,l) = g(k,l)
                        auxm(k,l) = markers(k,l)
                    END DO
                END DO
                DO k=1,nc
                    x(k) = auxx(id(k))
                    y(k) = auxy(id(k))
                    DO l=1,nb
                        g(k,l) = auxg(id(k),l)
                        markers(k,l) = auxm(id(k),l)
                    END DO
                END DO
            END IF
            !Update the information of in which site individuals are located
            nworg =0 !number of individuals per site
            worg = 0 !identity of individuals per site
            DO k=1,nc     
                nworg(x(k),y(k)) = nworg(x(k),y(k)) + 1
                IF (nworg(x(k),y(k)) > idensity) THEN
                    WRITE(6,*) 'There too many individuals here! ln196'
                END IF
                worg(x(k),y(k),nworg(x(k),y(k))) = k
            END DO                   
        END DO looptime ! loop in time

        IF (stable /= 0) THEN
            WRITE(6,*) simulationID,'Stabilized after',stable,'generations'
        ELSE
            WRITE(6,*) 'Stability not reached'
        END IF

        902 FORMAT(a15,1x,i15,1x,i15,1x,f15.5,1x,f15.2,1x,i15,1x,i15,1x,f15.2,1x,i15,1x,f15.2,1x,f15.2,1x,i15)
        OPEN (unit=17,file='summary.dat',status='unknown', position='append')
            WRITE(17,902) simulationID,nc,nf,mut,diff,disp_rad,radius,rg,nb,e0,sigk,igt
        CLOSE(17)
        replica_igt(ireplica) = igt
        replica_stable(ireplica) = stable
        filename = 'speciesplot'//trim(simulationID)//'.dat'
        OPEN (unit=8,file=filename,status='unknown')
            DO i=1,igt
                j = ispv(i)
                DO k=1,j
                    l = ispecies(i,k)
                    WRITE(8,*) x(l),y(l),i, (g(l,o),o=1,nb), (markers(l,o),o=1,nb)
                END DO
            END DO
        CLOSE(8)
        iptime = itime
        iitime = iitime + iptime

        filename = 'geneticdata'//trim(simulationID)//'.dat'
        OPEN (unit=9,file=filename,status='unknown',position='rewind')
            DO i=1,nc
                WRITE (9,'(10000i2)') (g(i,j),j=1,nb)
            END DO
        CLOSE(9)

        901 FORMAT(i4,1x,i4,1x,1000i1)
        filename = 'pop-'//trim(simulationID)//'.dat'
        OPEN (unit=10,file=filename,status='unknown',position='rewind')
            WRITE(10,*) iitime,nc,nf,deltat
            WRITE(10,*) mut,diff,disp_rad,mnpm
            WRITE(10,*) radius,rg,nb, is_dmi
            WRITE(10,*) e0,sigk
            WRITE(10,*) iREAD,inis,inisbit
            WRITE(10,*) independent_loci,window,nsets,replica
            DO i=1,nc
                WRITE (10,901) x(i),y(i),(g(i,j),j=1,nb)
            END DO
            WRITE(10,*) iseed
        CLOSE(10)

        DEALLOCATE (id)
        DEALLOCATE (g,goff,markers,moff,e,eoff)
        DEALLOCATE (fit)
        DEALLOCATE (x,y,neig,neigsp)
        DEALLOCATE (ineighbor,ineighborg)
        DEALLOCATE (ispv,ispecies)
        DEALLOCATE (worg,nworg)
        DEALLOCATE (idqx,idqy)        
        DEALLOCATE (previous_igt,previous_gdists)
        DEALLOCATE (iseed)
        DEALLOCATE (random_order)
        DEALLOCATE (auxx,auxy,auxg,auxm)
    END DO !replica

    OPEN (unit=18,file='summary_replica.dat',status='unknown', position='append')
        WRITE(18,903) nc,nf,mut,diff,disp_rad,radius,rg,nb,e0,sigk,(replica_stable(i),i=1,replica),(replica_igt(j),j=1,replica)
    CLOSE(18)
    903 FORMAT(i15,1x,i15,1x,f15.5,1x,f15.2,1x,i15,1x,i15,1x,f15.2,1x,i15,1x,f15.2,1x,f15.2,1x,10i15,1x,10i15)
    DEALLOCATE (replica_igt, replica_stable)
END DO simulation!change parameters

CLOSE(16)
CLOSE(17)
CLOSE(50)

CALL cpu_time(finish)
WRITE(6,*) 'Time', finish-start
CALL ELLAPSED_TIME

END PROGRAM topopatric

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Find spatial neighbors for this_k                                           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FINDNEIG(this_k, this_radius)
USE globals, ONLY: nc,nb,ineighbor,neigsp,x,y,nf,idq,idqx,idqy,nworg,worg,radius
IMPLICIT NONE
INTEGER(4) i,j,l!loop variables
INTEGER(4), intent(in) :: this_k, this_radius
INTEGER(4) dista,ix,iy,thisx,thisy

DO i=1,ineighbor(this_k)
    neigsp(this_k,ineighbor) = 0
END DO
ineighbor(this_k) = 0
ix = x(this_k)
iy = y(this_k)
thisx = 0
thisy = 0

IF (this_radius > 0.25*nf .OR. this_radius > radius) THEN
    searchpop: DO i = 1,nc
        IF(this_k == i) CYCLE searchpop! the focal individual is not included as a neighbor of itself
        dista = sqrt(REAL(min(abs(ix-x(i)), nf - abs(ix-x(i)))**2 + min(abs(iy-y(i)), nf - abs(iy-y(i)))**2))
        IF (dista <= this_radius) THEN
            ineighbor(this_k) =  ineighbor(this_k) + 1
            neigsp(this_k,ineighbor(this_k)) = i
        END IF
    END DO searchpop
ELSE
    loop1: DO i = 1,idq 
            thisx = ix+idqx(i)
            thisy = iy+idqy(i)
            IF (thisx > nf) thisx = thisx - nf
            IF (thisx < 1) thisx = nf + thisx
            IF (thisy > nf) thisy = thisy - nf
            IF (thisy < 1) thisy = nf + thisy
            loop2: DO j = 1,nworg(thisx,thisy)
                IF(worg(thisx,thisy,j) /= this_k) THEN ! the focal individual is not included as a neighbor of itself
                    ineighbor(this_k) =  ineighbor(this_k) + 1
                    neigsp(this_k,ineighbor) = worg(thisx,thisy,j)
                END IF
            END DO loop2
    END DO loop1
END IF

END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Find genetic neighbors for this_ k                                           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FINDGNEIG(this_k)
USE globals, ONLY: nb,ineighbor,ineighborg,neig,neigsp,g,rg,is_dmi
IMPLICIT NONE
INTEGER(4) i,l !loop variables
INTEGER(4), intent(in) :: this_k
INTEGER(4) dista

DO i=1,ineighborg(this_k)
    neig(this_k,ineighborg) = 0
END DO
ineighborg(this_k) = 0

DO i=1,ineighbor(this_k)
    dista = 0
    DO l=1,nb
        IF(is_dmi) THEN
            IF (g(this_k,l) + g(neigsp(this_k,ineighbor(this_k)),l) > 0) dista = dista + 1
            IF(dista > rg) EXIT
        ELSE
            dista = dista + abs(g(this_k,l)-g(neigsp(this_k,ineighbor(this_k)),l))
            IF(dista > rg) EXIT
        END IF
    END DO
    IF(dista <= rg) THEN
        ineighborg(this_k) = ineighborg(this_k) + 1
        neig(this_k,ineighborg) = neigsp(this_k,i)
    END IF
END DO

END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Find substitute for this_ k                                           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FINDSUB(this_k,kmother)
USE globals, ONLY: nc,nb,ineighbor,ineighborg,neig,neigsp,g,rg,is_dmi,radius
IMPLICIT NONE
INTEGER(4) i,j,l !loop variables
INTEGER, INTENT(IN) :: this_k
INTEGER, INTENT(INOUT) :: kmother
INTEGER(4) ineighbork, this_radius
INTEGER(4), ALLOCATABLE :: neigspk(:)
LOGICAL get_replacement
INTEGER(4) :: shuffle,shuffled
INTEGER(4), ALLOCATABLE, save :: random_order(:)
INTEGER(4) candidate
REAL random
ALLOCATE(neigspk(nc),random_order(nc))

this_radius = radius
get_replacement = .true.
replace: DO WHILE (get_replacement)
            ineighbork = ineighbor(this_k)
                DO i=1,ineighbork
                    neigspk(i) = neigsp(this_k,i)                    
                END DO               
                ! k's neighbors are the first set of candidates
                !!! Arrange neighbors in a random order
                DO i=1,ineighbork
                    random_order(i) = neigspk(i)
                END DO
                DO i=ineighbork,1,-1
                    CALL random_number(random)
                    shuffle = int(random*(ineighbork-1))+1
                    shuffled = random_order(i)
                    random_order(i) = random_order(shuffle)
                    random_order(shuffle) = shuffled
                END DO
                !Sucessively check if each neighbor is a suitable candidate for replacement
                DO i=1,ineighbork
                    candidate = random_order(i)
                    !check if candidate is going to be able to reproduce
                    !CALL FINDNEIG(candidate, radius)
                    !CALL FINDGNEIG(candidate)
                    IF(ineighborg(candidate) > 0) THEN
                        kmother = candidate
                        get_replacement = .false.
                        EXIT replace ! If this candidate is accepted, the search stops
                    END IF
                END DO !end of neighbor loop
                this_radius = this_radius + 1
                CALL FINDNEIG(this_k,this_radius)
        END DO replace
        CALL FINDNEIG(this_k,radius)
DEALLOCATE(neigspk,random_order)
END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Dispersal                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE DISPERSAL(k)
USE globals, ONLY:nf,diff,x,y,disp_rad
IMPLICIT NONE
INTEGER(4), intent(in) :: k
INTEGER(4) iix,iiy
REAL random, jump, turn , pi

pi = 3.141593
CALL random_number(random)
IF(random < diff) THEN
    ! calculate new position
    CALL random_number(random)
    jump = disp_rad*random
    CALL random_number(random)
    turn = 2*random*pi
    iix = x(k)+nint(cos(turn)*jump)
    iiy = y(k)+nint(sin(turn)*jump)
    IF (iix < 1 ) iix = iix + nf
    IF (iiy < 1 ) iiy = iiy + nf
    IF (iix > nf ) iix = iix - nf
    IF (iiy > nf ) iiy = iiy - nf
    x(k) = iix
    y(k) = iiy
END IF

END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Find species                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FINDSPECIES
USE globals, ONLY: x,y,nc,rg,g,nb,igt,ispecies,ispv,is_dmi
IMPLICIT NONE
INTEGER(4) i,j,k,l !loop variables
INTEGER(4), ALLOCATABLE :: species(:), auxy1(:),auxy2(:)
INTEGER(4) itot,i1,i2,icr,isp,dista,itest,isp0,ispold,ii,ji,o

ALLOCATE (species(nc),auxy1(nc),auxy2(nc))

itot = 0  ! count total population in groups
igt = 0   ! count number of groups
i2 = nc
DO i=1,i2  ! initialize aux2 -- contains all individuals not yet classified
    auxy2(i) = i
END DO 

DO WHILE (itot < nc)
    icr = auxy2(1)   !take first individual and find its species
    isp = 0
    ispold = 1
    i1 = 0
    auxy1 = 0
    loop1: DO i=1,i2
        ii = auxy2(i)
        dista = 0
        DO l=1,nb
            IF(is_dmi) THEN
                IF (g(icr,l) + g(ii,l) > 0) dista = dista + 1
            ELSE
                dista = dista + abs(g(icr,l) - g(ii,l))
            END IF
            IF(dista > rg) THEN
                i1 = i1 + 1
                auxy1(i1) = ii      !put creatures with dist > rg into aux1
                CYCLE loop1
            END IF
        END DO
        isp = isp + 1
        species(isp) = ii   !collect individuals with dist <= rg from icr
    END DO loop1

    IF(isp /= 0) THEN
        !check IF individuals in aux1 have to be included; put the rest in aux2
        itest = 1
        DO WHILE(itest /= 0)
            i2 = 0
            auxy2 = 0
            itest = 0
            isp0 = isp
            IF(i1 /= 0) THEN
                loop2:	DO i=1,i1
                    DO ji=ispold+1,isp0  
                        dista = 0
                        DO l=1,nb
                            IF(is_dmi) THEN
                                IF (g(auxy1(i),l) + g(species(ji),l) > 0) dista = dista + 1
                                IF(dista > rg) EXIT
                            ELSE
                                dista = dista + abs(g(auxy1(i),l) - g(species(ji),l)) 
                                IF(dista > rg) EXIT
                            END IF
                        END DO
                        IF(dista <= rg) THEN
                            isp = isp + 1 
                            species(isp) = auxy1(i)   ! colect the aux1 individual
                            itest = 1                 ! indicates that the process has to be repeated
                            CYCLE loop2
                        END IF
                    END DO
                    i2 = i2 + 1
                    auxy2(i2) = auxy1(i)  ! put individual in aux2
                END DO loop2
            END IF
            auxy1 = auxy2   ! aux1 contains the creatures not in the species
            i1 = i2
            ispold = isp0
        END DO
    END IF
    !IF (isp == 0) THEN
    !isp = 1
    !species(isp) = auxy2(i)
    !END IF
    itot = itot + isp    !total number of individuals classified into species
    igt = igt + 1        !number of species
    ! save species info
    DO i=1,isp
        ispecies(igt,i) = species(i) 
    END DO
    ispv(igt) = isp          ! number of individuals in species
END DO !END DO WHILE

DEALLOCATE (species,auxy1,auxy2)

END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!
! Check if number of species is stable          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE is_stable
USE globals, ONLY:g,nc,nb,igt,previous_igt,previous_gdists,keep_going,window,stable,deltat,sampled_times
IMPLICIT NONE
INTEGER(4) i,j,k !loop variables
INTEGER(4) igt_sum,dist,idists
REAL distsum, gdists_sum,gdists_mean
REAL igt_mean,check_mean

distsum=0
idists=0

DO i=1,nc
    DO j=i+1,nc
        dist=0
        idists=idists+1
        DO k=1,nb
            dist = dist + abs(g(i,k)-g(j,k))
        END DO
        WRITE(22,*) dist
        distsum = distsum + dist
    END DO
END DO

previous_gdists(sampled_times)=distsum/idists

IF(sampled_times>=window) THEN
    gdists_sum = 0
    DO i=sampled_times-window+1,sampled_times
        gdists_sum = gdists_sum+previous_gdists(i)
    END DO
    gdists_mean = REAL(gdists_sum)/window
    check_mean = REAL(previous_gdists(sampled_times-window+1))/gdists_mean
    IF (abs(1-check_mean)<=0.05) THEN
        keep_going =  .false.
        stable = (sampled_times-window+1)*deltat
    END IF
END IF

END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Shows how long it took to run the PROGRAM     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ELLAPSED_TIME
USE globals, ONLY:start,finish,hours,minutes,seconds
REAL time

time=0
hours=0
minutes=0
seconds=0

time= finish-start

hours=time/3600
time=MOD(time,3600.0)
minutes=time/60
seconds=MOD(time,60.0)

WRITE(6,'(a9,1x,i5,1x,a5,1x,i5,1x,a11,1x,i5,1x,a7)') 'Duration:', hours, 'hours', minutes, 'minutes and',seconds, 'seconds'

END