! topopatric.f90 -- v0.1
!Ayana B. Martins - 13/Jun/2016

module globals
!!Loop variables
integer(4) i,j,k,l 
logical keep_going
integer(4), allocatable :: id(:)

!!CPU time variables
real:: start, finish
integer(4) hours, minutes,seconds

!!Input variables
integer(4) :: ntime,nc,mnpm,nb,radius
real rg
integer(4) nf,deltat,sampled_times,window,stable
integer(4) iread,inis,inisbit
real diff,mut,m,critG
logical independent_loci, is_dmi, discretegen

!!Genotype and phenotype variables
integer(1), allocatable :: g(:,:)

!!Space and neighbor variables
integer(4), allocatable :: x(:),y(:)
integer(4), allocatable :: neig(:),neigsp(:),ispecies(:,:)
integer(4)  ineighbor,ineighborg,this_radius

!Species detection variables
integer(4) igt
integer(4), allocatable :: ispv(:)
integer(4), allocatable :: previous_igt(:)
real, allocatable:: previous_gdists(:)

!Fitness
integer(4), allocatable :: noff(:)

real, allocatable :: Fst(:), FstNull(:)
end module

program topopatric
use globals
implicit none
integer(4) nsets,iparameter
integer(4) replica,ireplica
integer(4) kmother,kmate,kcross
integer(4) ix,iy
real dista
integer(4) itime,iitime,iptime
real aux, random
character*30 simulationID
character*50 filename
character*20 simID,repID
real check_time
integer(4), allocatable :: replica_igt(:), replica_stable(:)

!!Random seed variables
integer, dimension(:), allocatable :: iseed
integer :: nseed

integer(4) candidate,this_neighbor
logical get_replacement
real isum

integer(4) :: shuffle,shuffled
integer(4), allocatable, save :: random_order(:)
integer(4), allocatable:: auxx(:), auxy(:), auxg(:,:)

logical, allocatable :: dead(:)

!!Temporarily stores offspring genotype
integer(1), allocatable :: goff(:,:)
call cpu_time(start)
! Read input data
!!Brief description of each parameter can be found in the input file
open(unit=7,file='input.in',status='old',position='rewind')
    read(7,*) ntime,nc,nf,deltat
    read(7,*) mut,diff,m,mnpm
    read(7,*) radius,critG,nb, is_dmi
    read(7,*) iread,inis,inisbit, discretegen
    read(7,*) independent_loci,window,nsets,replica
close(7)
open(unit=50,file='seed.in',status='old')  

!Absolute threshold is calculated from relative input
rg=nint(nb*critG)

if (mnpm /= 0) then
	write(6,*) 'Maybe you forgot that this version cannot handle the mnpm rule?'
	write(6,*) 'Duh!'
	stop
end if

if (nsets > 1) open(unit=16,file='par.in',status='old',position='rewind')


!Erase appendable output files
open(unit=17,file='summary.dat',status='unknown')
close(17, status='delete')

open(unit=18,file='summary_replica.dat',status='unknown')
close(18, status='delete')

simulation: do iparameter=1,nsets!change parameters

    allocate (replica_igt(replica),replica_stable(replica))
    if (nsets > 1) read (16, *) 
    if (iparameter > 99) then
        write(simID,'(i3)') iparameter
    elseif (iparameter > 9) then
        write(simID,'(i2)') iparameter
        simID = '0'//simID
    else
        write(simID,'(i1)') iparameter
        simID = '00'//simID
    end if
    do ireplica=1,replica

        if (ireplica > 99) then
            write(repID,'(i3)') ireplica
        elseif (ireplica > 9) then
            write(repID,'(i2)') ireplica
            repID = '0'//repID
        else
            write(repID,'(i1)') ireplica
            repID = '00'//repID
        end if
        simulationID = 'SIM'//trim(simID)//'.'//trim(repID)
        write (6,*) '-----------------',trim(simulationID),'----------------------'

        ! initialize random number generator
        call random_seed(size=nseed)
        allocate(iseed(nseed))
        read(50,*) iseed
        call random_seed(put=iseed)
        call random_number(aux)

        !Start a new simulation
        allocate(id(nc))
        allocate (g(nc,nb),goff(nc,nb))
        allocate(x(nc),y(nc),neig(nc),neigsp(nc))
        allocate (ispv(nc),ispecies(nc,nc))
        allocate(previous_igt(ntime),previous_gdists(ntime))
        allocate(random_order(nc))
        allocate(noff(nc))
        allocate(Fst(nb),FstNull(nb))
        allocate(auxx(nc),auxy(nc),auxg(nc,nb))
        allocate(dead(nc))
        
        !Set initial conditions of genotypes
        g = 0
        if(inisbit /= 0) then !random genomes
            do k=1,nc                    
                do i=1,nb
                    call random_number(aux)
                    if(aux > 0.5) g(k,i) = 1
                end do
            end do
        end if

        !Set initial location of individuals
        do i=1,nc
            if(inis == 0) then  !random position
                call random_number(aux)
                x(i) = int(aux*nf)+1
                call random_number(aux)
                y(i) = int(aux*nf)+1
             else !localized at the center
                x(i) = nf/2   
                y(i) = nf/2
             end if
        end do

        !Initialize individual's ids
        do k=1,nc
            id(k) = k
        end do
        
        dead=.false.
        
        !Erase appendable output files
        filename = 'degree_distplot'//trim(simulationID)//'.dat'
        open(unit=19,file=filename,status='unknown')
        close(19, status='delete')
        
        filename = 'Fst'//trim(simulationID)//'.dat'
        open(unit=20,file=filename,status='unknown')
        close(20, status='delete')
        
        filename = 'FstNull'//trim(simulationID)//'.dat'
        open(unit=21,file=filename,status='unknown')
        close(21, status='delete')
        
        filename = 'dist'//trim(simulationID)//'.dat'
        open(unit=22,file=filename,status='unknown')
        close(22, status='delete')
        
        ! Initialize time variable
        iitime = 0

        ! Time evolution: mating, mutation and diffusion
        keep_going=.true.
        itime=0
        sampled_times=0
        looptime: do while (keep_going)
            itime=itime+1
            !Mating
            noff = 0 !keep track of number of offspring per individual per generation
            goff = 0
            do k=1,nc
                get_replacement = .false. !initially, we assume that this individual is going to reproduce
                kmother = k
                call findneig(kmother,radius)
                !There is a random chance (m) that k is going to die and be replaced by the offspring
                !! of a pair of individuals nearby
                call random_number(random) 
                if (random < m) then
                    get_replacement = .true.
                    dead(k) = .true.
                end if
                !Additionally, if k has not enough potential partners, it will also
                !!!!! need to be replaced
                if(ineighborg < mnpm .or. ineighborg == 0) then
                    get_replacement = .true.
                    dead(k) = .true.
                end if
                this_radius = radius
                this_neighbor = ineighbor
                !If necessary, find a candidate to replace k
                replace: do while (get_replacement)
                    ! k's neighbors are the first set of candidates
                    !!! Arrange neighbors in a random order
                    do i=1,ineighbor
                        random_order(i) = neigsp(i)
                    end do
                    do i=ineighbor,1,-1
                        call random_number(random)
                        shuffle = int(random*(ineighbor-1))+1
                        shuffled = random_order(i)
                        random_order(i) = random_order(shuffle)
                        random_order(shuffle) = shuffled
                    end do
                    !Sucessively check if each neighbor is a suitable candidate
                    do i=1,ineighbor
                        candidate = random_order(i)
                        !check_ineighborg = 0 !necessary if using mnpm
                        ix = x(candidate)
                        iy = y(candidate)
                        check: do j = 1,nc !check if candidate is going to be able to reproduce
                            if(candidate == j) cycle check
                            dista = sqrt(real(min(abs(ix-x(j)), nf - abs(ix-x(j)))**2 + min(abs(iy-y(j)), nf - abs(iy-y(j)))**2))
                            if (dista <= radius) then
                                dista = 0
                                do l=1,nb
                                    if(is_dmi) then
                                        if (g(candidate,l) + g(j,l) > 0) dista = dista + 1
                                        if(dista > rg) exit
                                    else
                                        dista = dista + abs(g(candidate,l)-g(j,l))
                                        if(dista > rg) exit
                                    end if
                                end do
                                if(dista <= rg) then
                                    kmother = candidate
                                    get_replacement = .false.
                                    exit replace ! if this candidate is accepted, the search stops
                                end if
                            end if
                        end do check
                    end do !end of neighbor loop
                                    
                    !If all none of the neighbors is a suitable candidate,
                    !! increase radius of spatial neighborhood by 1
                    if (k==kmother) then
                        this_radius = this_radius + 1
                        call findneig(k, this_radius)
                    end if
                end do replace
                ! If k has been replaced, find neighbors of the new kmother
                if (kmother /= k) call findneig(kmother,radius)
                if(ineighborg == 0) then
                    write(6,*) 'ineighborg = 0!'
                    write(6,*) 'This should not be possible here (ln408)'
                    exit simulation
                end if
                !Only replace the individual if it has died
                if (dead(k)) then
                    call random_number(aux)
                    kmate = neig(int(aux*ineighborg)+1)    ! get a kmate /= kmother
                    noff(kmother) = noff(kmother) + 1
                    noff(kmate) = noff(kmate) + 1
                    if(independent_loci) then
                        do l=1,nb
                            call random_number(aux)
                            if (aux<0.5) then
                                goff(k,l)=g(kmother,l)
                            else
                                goff(k,l)=g(kmate,l)
                            end if
                        end do
                    else
                        call random_number(aux)
                        kcross = int(aux*(nb-1))+1         ! crossover point
                        call random_number(aux)
                        if(aux < 0.5) then
                            do l=1,kcross          ! copy first kcross bits of g(kmate) 
                                goff(k,l) = g(kmate,l)  ! and last bits of g(kmother) into g(k)
                            end do
                            do l=kcross+1,nb        
                                goff(k,l) = g(kmother,l)
                            end do
                        else
                            do l=1,kcross           ! copy first kcross bits of g(kmother) 
                                goff(k,l) = g(kmother,l) ! and last bits of g(kmate) into g(k)
                            end do
                            do l=kcross+1,nb       
                                goff(k,l) = g(kmate,l)
                            end do
                        end if
                    end if

                    if (.not. discretegen) then
                        !Replacement of parental genome and mutation
                        do l=1,nb                    
                            g(k,l) = goff(k,l)
                            call random_number(aux)
                            if(aux < mut) then
                                g(k,l) = 1-g(k,l)
                            end if
                        end do
                        !Dispersal
                        if(diff /= 0.0) call dispersal(k)
                    end if
                end if
            end do  ! end matings of the season
            
            ! If generations are discrete, replace parents and disperse
            !! only after the end of the mating season 
            if (discretegen) then
                do k=1,nc
                    if(dead(k)) then
                        do l=1,nb
                            g(k,l) = goff(k,l)
                            call random_number(aux)
                            if(aux < mut) then
                                g(k,l) = 1-g(k,l)
                            end if
                        end do
                        if(diff /= 0.0) call dispersal(k)
                    end if
                end do
            end if
            ! calculate species and other metrics every deltat
            check_time = (float(itime)/float(deltat) - itime/deltat)
            if (check_time == 0.0) then
                sampled_times=sampled_times+1
                call findspecies
                if (igt < 5) then
                    write(6,*) itime,itime+iitime,nc,igt,(ispv(k),k=1,igt)
                else
                    write(6,*) itime,itime+iitime,nc,igt,(ispv(k),k=1,4),'...'
                end if
                filename = 'dist'//trim(simulationID)//'.dat'
                open(unit=22,file=filename,status='unknown', position='append')
                    call is_stable
                close(22)
                if(itime==ntime) then
                    keep_going=.false.
                    stable = 0
                end if
                !Calculate the number of potential partners for each individual
                filename = 'degree_distplot'//trim(simulationID)//'.dat'
                open(unit=19,file=filename,status='unknown', position='append')
                    call degree_dist
                close(19)
                !Calculate Fst and generare a list to calculate the 95% confidence interval
                if (igt > 1) then
                    filename = 'Fst'//trim(simulationID)//'.dat'
                    open(unit=20,file=filename,status='unknown',position='append')
                    filename = 'FstNull'//trim(simulationID)//'.dat'
                    open(unit=21,file=filename,status='unknown',position='append')
                        call calcfst
                    close(20)
                    close(21)
                end if
            end if
            !Shuffle indexes so that individuals reproduce in random order
            do k=nc,1,-1
                call random_number(aux)
                shuffle = int(aux*(nc-1))+1
                shuffled = id(k)
                id(k) = id(shuffle)
                id(shuffle) = shuffled
            end do
            do k=1,nc
                auxx(k) = x(k)
                auxy(k) = y(k)
                do l=1,nb
                    auxg(k,l) = g(k,l)
                end do
            end do
            do k=1,nc
                x(k) = auxx(id(k))
                y(k) = auxy(id(k))
                do l=1,nb
                    g(k,l) = auxg(id(k),l)
                end do
            end do           
        end do looptime ! loop in time

        if (stable /= 0) then
            write(6,*) simulationID,'Stabilized after',stable,'generations'
        else
            write(6,*) 'Stability not reached'
        end if

        902 format(a15,1x,i15,1x,i15,1x,f15.5,1x,f15.2,1x,f15.2,1x,i15,1x,f15.2,1x,i15,1x,i15,1x,i15)
        open(unit=17,file='summary.dat',status='unknown', position='append')
            write(17,902) simulationID,nc,nf,mut,diff,m,radius,rg,nb,stable,igt
        close(17)
        replica_igt(ireplica) = igt
        replica_stable(ireplica) = stable
        filename = 'speciesplot'//trim(simulationID)//'.dat'
        open(unit=8,file=filename,status='unknown')
            do i=1,igt
            j = ispv(i)
                do k=1,j
                    l = ispecies(i,k)
                    write(8,*) x(l),y(l),i
                end do
            end do
        close(8)
        iptime = itime
        iitime = iitime + iptime

        filename = 'geneticdata'//trim(simulationID)//'.dat'
        open(unit=9,file=filename,status='unknown',position='rewind')
            do i=1,nc
                write (9,'(10000i2)') (g(i,j),j=1,nb)
            end do
        close(9)

        901 format(i4,1x,i4,1x,1000i1)
        filename = 'pop-'//trim(simulationID)//'.dat'
        open(unit=10,file=filename,status='unknown',position='rewind')
            write(10,*) iitime,nc,nf,deltat
            write(10,*) mut,diff,m,mnpm
            write(10,*)  radius,rg,nb, is_dmi
            write(10,*) iread,inis,inisbit,discretegen
            write(10,*) independent_loci,window,nsets,replica
            do i=1,nc
                write (10,901) x(i),y(i),(g(i,j),j=1,nb)
            end do
            write(10,*) iseed
        close(10)

        deallocate(id)
        deallocate (g,goff)
        deallocate(x,y,neig,neigsp)
        deallocate (ispv,ispecies)
        deallocate(previous_igt,previous_gdists)
        deallocate(iseed)
        deallocate(random_order)
        deallocate(noff)
        deallocate(Fst, FstNull)
        deallocate(auxx,auxy,auxg)
        deallocate(dead)

    end do !replica

    open(unit=18,file='summary_replica.dat',status='unknown', position='append')
        write(18,903) nc,nf,mut,diff,m,radius,rg,nb,(replica_stable(i),i=1,replica),(replica_igt(j),j=1,replica)
    close(18)
    903 format(i15,1x,i15,1x,f15.5,1x,f15.2,1x,f15.2,1x,i15,1x,f15.2,1x,i15,1x,10i15,1x,10i15)
    deallocate (replica_igt, replica_stable)
end do simulation!change parameters

close(16)
close(17)
close(50)

call cpu_time(finish)
write(6,*) 'Time', finish-start
call ellapsed_time

end program topopatric

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Find spatial and genetic neighbors for kmother                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine findneig(this_k, this_radius)
use globals, only: i,nc,l,nb,ineighbor,ineighborg,neig,neigsp,x,y,g,rg,nf,is_dmi
implicit none
integer(4), intent(in) :: this_k, this_radius
integer(4) dista,ix,iy

ineighbor = 0
ineighborg = 0
neig = 0   ! genetic neighbors
neigsp = 0 ! spatial neighbors
ix = x(this_k)
iy = y(this_k)
searchpop: do i = 1,nc
    if(this_k == i) cycle searchpop
    dista = sqrt(real(min(abs(ix-x(i)), nf - abs(ix-x(i)))**2 + min(abs(iy-y(i)), nf - abs(iy-y(i)))**2))
    if (dista <= this_radius) then
        ineighbor = ineighbor + 1
        neigsp(ineighbor) = i
        dista = 0
        do l=1,nb
            if(is_dmi) then
                if (g(this_k,l) + g(neigsp(ineighbor),l) > 0) dista = dista + 1
                if(dista > rg) exit
            else
                dista = dista + abs(g(this_k,l)-g(neigsp(ineighbor),l))
                if(dista > rg) exit
            end if
        end do
        if(dista <= rg) then
            ineighborg = ineighborg + 1
            neig(ineighborg) = neigsp(ineighbor)
        end if
    end if
end do searchpop

end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Dispersal                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine dispersal(k)
use globals, only:nf,diff,x,y,radius
implicit none
integer(4), intent(in) :: k
integer(4) iix,iiy
real random, jump, turn , pi

pi = 3.141593
call random_number(random)
if(random < diff) then
    ! calculate new position
    call random_number(random)
    jump = radius*random
    call random_number(random)
    turn = 2*random*pi
    iix = x(k)+nint(cos(turn)*jump)
    iiy = y(k)+nint(sin(turn)*jump)
    if (iix < 1 ) iix = iix + nf
    if (iiy < 1 ) iiy = iiy + nf
    if (iix > nf ) iix = iix - nf
    if (iiy > nf ) iiy = iiy - nf
    x(k) = iix
    y(k) = iiy
end if

end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Find species                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine findspecies
use globals, only: i,j,k,l,x,y,nc,rg,g,nb,igt,ispecies,ispv,is_dmi
implicit none
integer(4), allocatable :: species(:), auxy1(:),auxy2(:)
integer(4) itot,i1,i2,icr,isp,dista,itest,isp0,ispold,ii,ji

allocate (species(nc),auxy1(nc),auxy2(nc))

itot = 0  ! count total population in groups
igt = 0   ! count number of groups
i2 = nc
do i=1,i2  ! initialize aux2 -- contains all individuals not yet classified
    auxy2(i) = i
end do 

do while (itot < nc)
    icr = auxy2(1)   !take first individual and find its species
    isp = 0
    ispold = 1
    i1 = 0
    auxy1 = 0
    loop1: do i=1,i2
        ii = auxy2(i)
        dista = 0
        do l=1,nb
            if(is_dmi) then
                if (g(icr,l) + g(ii,l) > 0) dista = dista + 1
            else
                dista = dista + abs(g(icr,l) - g(ii,l))
            end if
            if(dista > rg) then
                i1 = i1 + 1
                auxy1(i1) = ii      !put creatures with dist > rg into aux1
                cycle loop1
            end if
        end do
        isp = isp + 1
        species(isp) = ii   !collect individuals with dist <= rg from icr
    end do loop1

    if(isp /= 0) then
        !check if individuals in aux1 have to be included; put the rest in aux2
        itest = 1
        do while(itest /= 0)
            i2 = 0
            auxy2 = 0
            itest = 0
            isp0 = isp
            if(i1 /= 0) then
                loop2:	do i=1,i1
                    do ji=ispold+1,isp0  
                        dista = 0
                        do l=1,nb
                            if(is_dmi) then
                                if (g(auxy1(i),l) + g(species(ji),l) > 0) dista = dista + 1
                                if(dista > rg) exit
                            else
                                dista = dista + abs(g(auxy1(i),l) - g(species(ji),l)) 
                                if(dista > rg) exit
                            end if
                        end do
                        if(dista <= rg) then
                            isp = isp + 1 
                            species(isp) = auxy1(i)   ! colect the aux1 individual
                            itest = 1                 ! indicates that the process has to be repeated
                            cycle loop2
                        end if
                    end do
                    i2 = i2 + 1
                    auxy2(i2) = auxy1(i)  ! put individual in aux2
                end do loop2
            end if
            auxy1 = auxy2   ! aux1 contains the creatures not in the species
            i1 = i2
            ispold = isp0
        end do
    end if
    if (isp == 0) then
    isp = 1
    species(isp) = auxy2(i)
    end if
    itot = itot + isp    !total number of individuals classified into species
    igt = igt + 1        !number of species
    ! save species info
    do i=1,isp
        ispecies(igt,i) = species(i) 
    end do
    ispv(igt) = isp          ! number of individuals in species
end do !end do while

deallocate (species,auxy1,auxy2)

end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate dist of degrees                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine degree_dist
use globals, only:i,j,l,g,nc,nb,rg,is_dmi, noff
implicit none
integer(4) dista,partners

do i=1,nc
    partners= 0
    do j=1,nc
        if (i==j) cycle
        dista = 0
        do l=1,nb
            if(is_dmi) then
                if (g(i,l) + g(j,l) > 0) dista = dista + 1
            else
                dista = dista + abs(g(i,l) - g(j,l)) 
            end if
        end do
        if (dista <= rg) partners = partners + 1
    end do
    write(19,*) partners, noff(i)
end do

end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate Fst for n subpopulations            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calcfst
use globals, only:i,j,k,l,x,y,nf,g,nc,nb,Fst,igt,ispv,ispecies
implicit none
real, allocatable:: freqs(:,:),Hs(:,:)
integer(4) dimH, shuffle, shuffled, ini, fin
real meanHs
integer(4), allocatable:: randomized(:)
real random

dimH = igt + 1
allocate(freqs(dimH,nb),Hs(dimH,nb))
allocate(randomized(nc))

freqs = 0
Hs = 0
Fst=0

do i=1,igt
    do l=1,nb
        do j=1,ispv(i) !cycle through individuals of species j
            freqs(i,l) = freqs(i,l) + g(ispecies(i,j),l)
        end do
        freqs(i,l) = freqs(i,l)/real(ispv(i))
    end do
end do

do l=1,nb
    do i=1,igt
        freqs(dimH,l) = freqs(dimH,l) + freqs(i,l) 
    end do
    freqs(dimH,l) = freqs(dimH,l)/real(igt)
end do

do i=1,dimH
    do l=1,nb
        Hs(i,l) = 2*freqs(i,l)*(1-freqs(i,l))
    end do
end do

do l =1,nb
    meanHs = 0
    do i = 1, igt
        meanHs = meanHs + Hs(i,l)
    end do
    meanHs = meanHs/real(igt)
    Fst(l) = (Hs(dimH,l)- meanHs) / Hs(dimH,l)
end do

write(20,*) (Fst(l),l=1,nb)

do k =1,nc
    randomized(k) = k
end do

do k=nc,1,-1
    call random_number(random)
    shuffle = int(random*(nc-1))+1
    shuffled = randomized(k)
    randomized(k) = randomized(shuffle)
    randomized(shuffle) = shuffled
end do

freqs = 0
Hs = 0
Fst=0
ini = 1
fin = 0
do i=1,igt
    ini = ini + fin
    fin = fin + ispv(i)
    do l=1,nb
        do j=ini,fin !cycle through individuals of species j
            freqs(i,l) = freqs(i,l) + g(randomized(j),l)
        end do
        freqs(i,l) = freqs(i,l)/real(ispv(i))
    end do
end do

do l=1,nb
    do i=1,igt
        freqs(dimH,l) = freqs(dimH,l) + freqs(i,l) 
    end do
    freqs(dimH,l) = freqs(dimH,l)/real(igt)
end do

do i=1,dimH
    do l=1,nb
        Hs(i,l) = 2*freqs(i,l)*(1-freqs(i,l))
    end do
end do

do l =1,nb
    meanHs = 0
    do i = 1, igt
        meanHs = meanHs + Hs(i,l)
    end do
    meanHs = meanHs/real(igt)
    Fst(l) = (Hs(dimH,l)- meanHs) / Hs(dimH,l)
end do

write(21,*) (Fst(l),l=1,nb)

deallocate(freqs,Hs)
deallocate(randomized)

!! First implemented comparing groups in space
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! integer(4), allocatable:: nsub(:)
! integer(4), allocatable:: subpops(:,:)
! real, allocatable:: freqs(:,:),Hs(:,:)
! integer(4) division, limit
! 
! allocate(subpops(4,nc))
! allocate(nsub(4))
! allocate(freqs(5,nb),Hs(5,nb))
! division = 2
! limit = nf/division
! 
! freqs = 0
! nsub = 0
! do k=1,nc
!         if (x(k) > limit .and. y(k) > limit) then
!             nsub(1) = nsub(1) + 1
!             subpops(1,nsub(1)) = k
!         else if (x(k) <= limit .and. y(k) > limit) then
!             nsub(2) = nsub(2) + 1
!             subpops(2,nsub(2)) = k
!         else if (x(k) <= limit .and. y(k) <= limit) then
!             nsub(3) = nsub(3) + 1
!             subpops(3,nsub(3)) = k
!         else
!             nsub(4) = nsub(4) + 1
!             subpops(4,nsub(4)) = k
!         end if            
! end do
! 
! do j=1,4
!     do l=1,nb
!         do i=1,nsub(j)
!             freqs(j,l) = freqs(j,l) + g(subpops(j,i),l)
!         end do
!         freqs(j,l) = freqs(j,l)/real(nsub(j))
!     end do
! end do
! 
! do l=1,nb
!     do i=1,4
!         freqs(5,l) = freqs(5,l) + freqs(i,l) 
!     end do
!     freqs(5,l) = freqs(5,l)/4.0
! end do
! 
! do j=1,5
!     do l=1,nb
!         Hs(j,l) = 2*freqs(j,l)*(1-freqs(j,l))
!     end do
! end do
! 
! Fst=0
! do l =1,nb
!     Fst(l) = (Hs(5,l)-((Hs(1,l)+Hs(2,l)+Hs(3,l)+Hs(4,l))/4)) / Hs(5,l)
! end do
! 
! write(6,'(10f15.3)') (Fst(l),l=1,nb)
! !pause
! 
! deallocate(nsub,freqs,Hs)
! deallocate(subpops)
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!
! Check if number of species is stable          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine is_stable
use globals, only:i,j,k,g,nc,nb,igt,previous_igt,previous_gdists,keep_going,window,stable,deltat,sampled_times
implicit none
integer(4) igt_sum,dist,idists
real distsum, gdists_sum,gdists_mean
real igt_mean,check_mean

distsum=0
idists=0

do i=1,nc
    do j=i+1,nc
        dist=0
        idists=idists+1
        do k=1,nb
            dist = dist + abs(g(i,k)-g(j,k))
        end do
        write(22,*) dist
        distsum = distsum + dist
    end do
end do

previous_gdists(sampled_times)=distsum/idists

if(sampled_times>=window) then
    gdists_sum = 0
    do i=sampled_times-window+1,sampled_times
        gdists_sum = gdists_sum+previous_gdists(i)
    end do
    gdists_mean = real(gdists_sum)/window
    check_mean = real(previous_gdists(sampled_times-window+1))/gdists_mean
    if (abs(1-check_mean)<=0.05) then
        keep_going =  .false.
        stable = (sampled_times-window+1)*deltat
    end if
end if

end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Shows how long it took to run the program     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ellapsed_time
use globals, only:start,finish,hours,minutes,seconds
real time

time=0
hours=0
minutes=0
seconds=0

time= finish-start

hours=time/3600
time=MOD(time,3600.0)
minutes=time/60
seconds=MOD(time,60.0)

write(6,'(a9,1x,i5,1x,a5,1x,i5,1x,a11,1x,i5,1x,a7)') 'Duration:', hours, 'hours', minutes, 'minutes and',seconds, 'seconds'

end