! topopatric.f90 -- v0.2
! panmictic -- v0.1
!Ayana B. Martins - 20/Jun/2016

module globals
!!Loop variables
integer(4) i,j,k,l,o 
logical keep_going
integer(4), allocatable :: id(:)

!!CPU time variables
real:: start, finish,local_start,local_finish
integer(4) hours, minutes,seconds

!!Input variables
integer(4) :: ntime,nc,nb,ntrials
real rg
integer(4) deltat,sampled_times,window,stable
integer(4) iread,inis,inisbit
real mut,m,critG
logical independent_loci, is_dmi, discretegen

!!Genotype and phenotype variables
integer(1), allocatable :: g(:,:)
integer(1), allocatable :: markers(:,:)

!Species detection variables
integer(4) igt
integer(4), allocatable :: ispecies(:,:)
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

integer(4) candidate, mate, these_trials
logical get_replacement
real isum

integer(4) trials

integer(4) :: shuffle,shuffled
integer(4), allocatable, save :: random_order(:), random_order_mates(:)
integer(4), allocatable:: auxx(:), auxy(:), auxg(:,:), auxm(:,:)

logical, allocatable :: dead(:)

!!Temporarily stores offspring genotype
integer(1), allocatable :: goff(:,:)
integer(1), allocatable :: moff(:,:)
call cpu_time(start)
! Read input data
!!Brief description of each parameter can be found in the input file
open(unit=7,file='input.in',status='old',position='rewind')
    read(7,*) ntime,nc,ntrials
    read(7,*) mut,m,deltat
    read(7,*) critG,nb,window
    read(7,*) independent_loci,discretegen,is_dmi
    read(7,*) iread,inis,inisbit
    read(7,*) nsets,replica
close(7)
open(unit=50,file='seed.in',status='old')  

!Absolute threshold is calculated from relative input
rg=nint(nb*critG)

if (nsets > 1) open(unit=16,file='par.in',status='old',position='rewind')

!Erase appendable output files
open(unit=17,file='summary.dat',status='unknown')
close(17, status='delete')

open(unit=18,file='summary_replica.dat',status='unknown')
close(18, status='delete')

simulation: do iparameter=1,nsets!change parameters

    allocate (replica_igt(replica),replica_stable(replica))
    if (nsets > 1) read (16, *) ntrials
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
        allocate (g(nc,nb),markers(nc,nb))
        allocate (goff(nc,nb), moff(nc,nb))
        allocate (ispv(nc),ispecies(nc,nc))
        allocate(previous_igt(ntime),previous_gdists(ntime))
        allocate(random_order(nc),random_order_mates(nc))
        allocate(noff(nc))
        allocate(Fst(nb),FstNull(nb))
        allocate(auxg(nc,nb),auxm(nc,nb))
        allocate(dead(nc))
        
        !Set initial conditions of genotypes
        g = 0
        markers = 0
        if(inisbit /= 0) then !random genomes
            do k=1,nc                    
                do i=1,nb
                    call random_number(aux)
                    if(aux > 0.5) g(k,i) = 1
                    !markers are also random
                    call random_number(aux)
                    if(aux > 0.5) markers(k,i) = 1
                end do
            end do
        end if

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
            moff = 0
            do k=1,nc
                get_replacement = .false. !initially, we assume that this individual is going to reproduce
                kmother = k
                !There is a random chance (m) that k is going to die and be replaced by the offspring
                !! of a pair of individuals nearby
                call random_number(random) 
                if (random < m) then
                    get_replacement = .true.
                    dead(k) = .true.                                      
                end if
                !Only replace the individual if it has died
                if (dead(k)) then
                    candidate = 0
                    kmother = 0
                    mate = 0
                    kmate = 0
                    do i=1,nc
                        random_order(i) = i
                    end do
                    do i=nc,1,-1
                        call random_number(random)
                        shuffle = int(random*(nc-1))+1
                        shuffled = random_order(i)
                        random_order(i) = random_order(shuffle)
                        random_order(shuffle) = shuffled
                    end do
                    rep: do while (get_replacement)
                        candidate = candidate + 1
                        kmother = random_order(candidate)
!                       write(6,*) k, candidate, kmother                        
                        do i=1,nc
                            random_order_mates(i) = i
                        end do
                        do i=nc,1,-1
                            call random_number(random)
                            shuffle = int(random*(nc-1))+1
                            shuffled = random_order_mates(i)
                            random_order_mates(i) = random_order_mates(shuffle)
                            random_order_mates(shuffle) = shuffled
                        end do
                        trials = 0
                        these_trials = ntrials
                        searchpop: do while (trials < these_trials)
                            trials = trials + 1
                            if (trials > nc) cycle rep
                            mate = random_order_mates(trials)
                            if(mate == k) then
                                these_trials = these_trials + 1
                                cycle searchpop !don't mate with the dead D:
                            end if
                            if(mate == kmother) then
                                these_trials = these_trials + 1
                                cycle searchpop !no self-fertilization                            
                            end if
                            dista = 0
                            do l=1,nb
                                if(is_dmi) then
                                    if (g(kmother,l) + g(mate,l) > 0) dista = dista + 1
                                    if(dista > rg) exit
                                else
                                    dista = dista + abs(g(kmother,l)-g(mate,l))
                                    if(dista > rg) exit
                                end if
                            end do
                            if(dista <= rg) then
                                kmate = mate
                                get_replacement = .false.
                                exit searchpop
                            end if
                        end do searchpop
                        if (kmate == 0 .and. candidate == nc) then
                            write(6,*) 'There is no one to replace k!!'
                            exit
                        end if
                    end do rep!get replacement 

                    noff(kmother) = noff(kmother) + 1
                    noff(kmate) = noff(kmate) + 1
                    !Generate offspring                
                    if(independent_loci) then
                        do l=1,nb
                            call random_number(random)
                            if (random<0.5) then
                                goff(k,l)=g(kmother,l)
                            else
                                goff(k,l)=g(kmate,l)
                            end if
                        end do
                        do l=1,nb
                            call random_number(random)
                            if (random<0.5) then
                                moff(k,l)=markers(kmother,l)
                            else
                                moff(k,l)=markers(kmate,l)
                            end if
                        end do
                    else
                        call random_number(random)
                        kcross = int(random*(nb-1))+1         ! crossover point
                        call random_number(random)
                        if(random < 0.5) then
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
                        do l=1,nb                    
                            markers(k,l) = moff(k,l)
                            call random_number(aux)
                            if(aux < mut) then
                                markers(k,l) = 1-markers(k,l)
                            end if
                        end do
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
                        do l=1,nb
                            markers(k,l) = moff(k,l)
                            call random_number(aux)
                            if(aux < mut) then
                                markers(k,l) = 1-markers(k,l)
                            end if
                        end do
                    end if
                end do
            end if
            ! calculate species and other metrics every deltat
            check_time = (float(itime)/float(deltat) - itime/deltat)
            if (check_time == 0.0 .or. itime==ntime) then
                sampled_times=sampled_times+1
                call cpu_time(local_start)
                call findspecies
                call cpu_time(local_finish)
                write(6,*) 'Time: findspecies'
                call ellapsed_time(local_start,local_finish)
                if (igt < 5) then
                    write(6,*) itime,itime+iitime,nc,igt,(ispv(k),k=1,igt)
                else
                    write(6,*) itime,itime+iitime,nc,igt,(ispv(k),k=1,4),'...'
                end if
                filename = 'dist'//trim(simulationID)//'.dat'
                open(unit=22,file=filename,status='unknown', position='append')
                !Calculate the number of potential partners for each individual
                filename = 'degree_distplot'//trim(simulationID)//'.dat'
                open(unit=19,file=filename,status='unknown', position='append')
                    call cpu_time(local_start)
                    call is_stable
                    call cpu_time(local_finish)
                    write(6,*) 'Time: is_stable'
                    call ellapsed_time(local_start,local_finish)
                close(22)
                close(19)
                !Calculate Fst and generare a list to calculate the 95% confidence interval
                if (igt > 1) then
                    filename = 'Fst'//trim(simulationID)//'.dat'
                    open(unit=20,file=filename,status='unknown',position='append')
                    filename = 'FstNull'//trim(simulationID)//'.dat'
                    open(unit=21,file=filename,status='unknown',position='append')
                        call cpu_time(local_start)
                        call calcfst
                        call cpu_time(local_finish)
                        write(6,*) 'Time: calcfst'
                        call ellapsed_time(local_start,local_finish)                               
                    close(20)
                    close(21)
                end if
                if (itime==ntime) then
                    keep_going=.false.
                end if
            else
                write(6,*) itime,itime+iitime
            end if
            if (keep_going) then
                !Shuffle indexes so that individuals reproduce in random order
                do k=nc,1,-1
                    call random_number(aux)
                    shuffle = int(aux*(nc-1))+1
                    shuffled = id(k)
                    id(k) = id(shuffle)
                    id(shuffle) = shuffled
                end do
                do k=1,nc
                    do l=1,nb
                        auxg(k,l) = g(k,l)
                        auxm(k,l) = markers(k,l)
                    end do
                end do
                do k=1,nc
                    do l=1,nb
                        g(k,l) = auxg(id(k),l)
                        markers(k,l) = auxm(id(k),l)
                    end do
                end do
            end if
        end do looptime ! loop in time

        if (stable /= 0) then
            write(6,*) simulationID,'Stabilized after',stable,'generations'
        else
            write(6,*) 'Stability not reached'
        end if

        902 format(a15,1x,i15,1x,f15.5,1x,f15.2,1x,f15.2,1x,i15,1x,i15,1x,i15,1x,i15)
        open(unit=17,file='summary.dat',status='unknown', position='append')
            write(17,902) simulationID,nc,mut,m,critG,nb,ntrials,stable,igt
        close(17)
        replica_igt(ireplica) = igt
        replica_stable(ireplica) = stable
        filename = 'speciesplot'//trim(simulationID)//'.dat'
        open(unit=8,file=filename,status='unknown')
            do i=1,igt
                j = ispv(i)
                do k=1,j
                    l = ispecies(i,k)
                    write(8,*) i, (g(l,o),o=1,nb), (markers(l,o),o=1,nb)
                end do
            end do
        close(8)
        iptime = itime
        iitime = iitime + iptime

        filename = 'geneticdata'//trim(simulationID)//'.dat'
        open(unit=9,file=filename,status='unknown',position='rewind')
            do i=1,nc
                write (9,'(500000i1)') (g(i,j),j=1,nb)
            end do
        close(9)

        901 format(i4,1x,i4,1x,500000i1)
        filename = 'pop-'//trim(simulationID)//'.dat'
        open(unit=10,file=filename,status='unknown',position='rewind')
            write(10,*) ntime,nc,ntrials
            write(10,*) mut,m,deltat
            write(10,*) critG,nb,window
            write(10,*) independent_loci,discretegen,is_dmi
            write(10,*) iread,inis,inisbit
            write(10,*) nsets,replica,aux
            do i=1,nc
                write (10,901) (g(i,j),j=1,nb)
            end do
            write(10,*) iseed
        close(10)

        deallocate(id)
        deallocate (g,goff,markers,moff)
        deallocate (ispv,ispecies)
        deallocate(previous_igt,previous_gdists)
        deallocate(iseed)
        deallocate(random_order,random_order_mates)
        deallocate(noff)
        deallocate(Fst, FstNull)
        deallocate(auxg,auxm)
        deallocate(dead)

    end do !replica

    open(unit=18,file='summary_replica.dat',status='unknown', position='append')
        write(18,903) nc,mut,m,critG,nb,ntrials,(replica_stable(i),i=1,replica),(replica_igt(j),j=1,replica)
    close(18)
    903 format(i15,1x,f15.5,1x,f15.2,1x,f15.2,1x,i15,1x,i15,1x,10i15,1x,10i15)
    deallocate (replica_igt, replica_stable)
end do simulation!change parameters

close(16)
close(17)
close(50)

call cpu_time(finish)
write(6,*) 'Time', finish-start
call ellapsed_time(start,finish)

end program topopatric

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Find species                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine findspecies
use globals, only: i,j,k,l,nc,rg,g,nb,igt,ispecies,ispv,is_dmi
implicit none
integer(4), allocatable :: species(:), auxy1(:),auxy2(:)
integer(4) itot,i1,i2,icr,isp,dista,itest,isp0,ispold,ii,ji,o

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
    !if (isp == 0) then
    !isp = 1
    !species(isp) = auxy2(i)
    !end if
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

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! Calculate dist of degrees                     !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine degree_dist
! use globals, only:i,j,l,g,nc,nb,rg,is_dmi, noff
! implicit none
! integer(4) dista,partners
! 
! do i=1,nc
!     partners= 0
!     do j=1,nc
!         if (i==j) cycle
!         dista = 0
!         do l=1,nb
!             if(is_dmi) then
!                 if (g(i,l) + g(j,l) > 0) dista = dista + 1
!             else
!                 dista = dista + abs(g(i,l) - g(j,l)) 
!             end if
!         end do
!         if (dista <= rg) partners = partners + 1
!     end do
!     write(19,*) partners, noff(i)
! end do
! 
! end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate Fst for n subpopulations            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calcfst
use globals, only:i,j,k,l,g,nc,nb,Fst,igt,ispv,ispecies
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

end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!
! Check if number of species is stable          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine is_stable
use globals, only:i,j,l,g,nc,nb,rg,is_dmi,noff,igt,previous_igt,previous_gdists,keep_going,window,stable,deltat,sampled_times
implicit none
integer(4) igt_sum,dist,idists
real distsum, gdists_sum,gdists_mean
real igt_mean,check_mean
integer(4) dista,partners

distsum=0
idists=0

do i=1,nc
    partners = 0
    do j=1,nc
        if (i==j) cycle
        dist=0
        idists=idists+1
        do l=1,nb
            if(is_dmi) then
                if (g(i,l) + g(j,l) > 0) dist = dist + 1
            else
                dist = dist + abs(g(i,l) - g(j,l)) 
            end if
        end do
        if (dist <= rg) partners = partners + 1
        write(22,*) dist
        distsum = distsum + dist
    end do
    write(19,*) partners, noff(i)
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
subroutine ellapsed_time(this_start,this_finish)
use globals, only:start,finish,hours,minutes,seconds
real, intent(in) :: this_start,this_finish
real time

time=0
hours=0
minutes=0
seconds=0

time= this_finish-this_start

hours=time/3600
time=MOD(time,3600.0)
minutes=time/60
seconds=MOD(time,60.0)

write(6,'(a9,1x,i5,1x,a5,1x,i5,1x,a11,1x,i5,1x,a7)') 'Duration:', hours, 'hours', minutes, 'minutes and',seconds, 'seconds'

end