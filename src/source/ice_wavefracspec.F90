! Lettie Roach, NIWA/VUW, June 2016
! Some of this modified from Chris Horvat's Matlab original
! Cecilia Bitz changed it a lot in Oct 2018
! New module for CICE to generate waves and break the FSD
!
      module ice_wavefracspec

      use ice_kinds_mod
      use ice_constants
      use ice_domain_size, only: nfsd, ncat, max_ntrcr, max_blocks
      use ice_calendar, only: dt
      use ice_fileunits, only: flush_fileunit
 
      implicit none
      private
      public :: wave_frac_fsd 

      logical (kind=log_kind), public :: &
         wave_spec                ! namelist parameter that indicates using ww3 spectrum

      integer (kind=int_kind), parameter :: &
          nfreqww3    = 25        ! number of spectra in ww3

      real, parameter  :: &
                swh_minval = 0.01_dbl_kind, & ! minimum value of wave height (m)
                straincrit =0.00003,        & ! critical strain
                D=10000.,                   & ! domain size
                dx = c1,                    & ! domain spacing
                threshold=10                  ! peak-finding threshold 
                                              ! (points are defined to be extrema
                                              ! if they are a local maximum or minimum
                                              ! over a distance of 10m on both sides, 
                                              ! based on the observations of Toyota et al. 
                                              ! (2011) who find this to be the order of
                                              ! the smallest floe size affected by wave fracture 
      integer (kind=int_kind) :: &
          nx = 10000                          ! number of points in domain
 

!=======================================================================

      contains

!=======================================================================

     function get_damfstd_wave(amfstd_init, fracture_hist, frac, fracbreak) &
                              result(d_amfstd)

     real (kind=dbl_kind), dimension (nfsd), intent(in) :: &
         amfstd_init, fracture_hist

     real (kind=dbl_kind), dimension (nfsd,nfsd), intent(in) :: &
         frac

     real (kind=dbl_kind), intent(in) :: & 
         fracbreak

     real (kind=dbl_kind), dimension (nfsd) :: &
         d_amfstd, gain, omega, histtmp

     real (kind=dbl_kind) :: & 
         runsum

 
      integer (kind=int_kind) :: k

      histtmp=fracture_hist
      omega(1) = c0 ! never fracture the smallest floes
      runsum = c0
      do k = 2,nfsd
	   runsum=runsum+histtmp(k-1)
           omega(k) = amfstd_init(k)*runsum
      end do

      omega(:) = omega(:)*fracbreak/dt

      if (SUM(omega).gt.c1+puny) stop &
          'omega cannot be greater than 1, waves'
                        
      do k =1,nfsd-1
           gain(k) = SUM(omega*frac(:,k)) 
      end do
      gain(nfsd) = c0  ! largest floes never gain more by fracture

      d_amfstd(:) = gain(:) - omega(:)
 
      if (SUM(d_amfstd(:)).gt.puny) stop 'area not cons, waves'

      end  function get_damfstd_wave


!=======================================================================
!  Author: Lettie Roach, NIWA, 2018
! 
! Given fracture histogram computed from local wave spectrum, evolve FSD

     subroutine wave_frac_fsd

     use ice_communicate, only: my_task
     use ice_domain, only: nblocks, blocks_ice
     use ice_blocks, only: block, get_block, nx_block, ny_block
     use ice_state, only: nt_fsd, aice, aicen, vice, &
                          vicen, trcrn, nt_qice
     use ice_flux, only: wave_spectrum
     use ice_fsd, only: d_afsd_wave, d_amfstd_wave, check_mfstd
 
     ! local variables
      type (block) :: &
         this_block      ! block information for current block

      integer (kind=int_kind) :: &  
        iblk, n, k, i, j, t, &
        ilo, ihi, jlo, jhi    ! beginning and end of physical domain

      real (kind=dbl_kind), dimension (nfsd) :: &
        fracture_hist

      real (kind=dbl_kind) :: &
         fracbreak                  ! how much of cell experiences breaking

      real (kind=dbl_kind), dimension (nfsd, nfsd) :: &
        frac    

      real (kind=dbl_kind) :: &
        hbar, elapsed_t, subdt, normalizer

      real (kind=dbl_kind), dimension (nfsd) :: &
         amfstd_init, &     ! tracer array
         amfstd_tmp, d_amfstd_tmp, &
         check_dt


     !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block,n,k,hbar,fracbreak,frac,normalizer,fracture_hist,amfstd_init,amfstd_tmp,elapsed_t,d_amfstd_tmp,check_dt,subdt)
!!!     !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block,hbar,fracbreak,frac,normalizer,amfstd_init,amfstd_tmp,d_amfstd_tmp,n,k,check_dt,subdt,elapsed_t,trcrn,d_amfstd_wave,d_afsd_wave)

     do iblk = 1, nblocks

         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         d_afsd_wave(ilo:ihi,jlo:jhi,:,iblk) = c0

         if (ANY(trcrn(ilo:ihi,jlo:jhi,nt_fsd:nt_fsd+nfsd-1,:,iblk).lt.c0)) stop 'neg b4-wb'

         do j = jlo, jhi
         do i = ilo, ihi

            if ((aice(i,j,iblk).gt.0.0001).and.(.NOT.(ALL(wave_spectrum(i,j,:,iblk).lt.puny)))) then

                hbar = vice(i,j,iblk) / aice(i,j,iblk)

                call  wave_frac(hbar, wave_spectrum(i,j,:,iblk), fracture_hist, fracbreak)

!         write(*,10) 'fracbreak',hbar,fracbreak,wave_spectrum(i,j,:,iblk), fracture_hist, size(wave_spectrum(i,j,:,iblk)), size(fracture_hist)
!         write(*,*)'fracbreak*', fracbreak
!10       format(1X,A9,2x,F8.5,2X,F8.4,5x,25(F8.5,1x),5x,12(F8.5,1x),5x,I5,2x,I5)


                if (fracbreak.gt.puny) then

                  ! frac does not vary within subcycle and does not depend on cat
                  frac(:,:) = c0
                  do k = 2, nfsd
                    frac(k,:k-1) = fracture_hist(:k-1)
                    normalizer = SUM(frac(k,:))
                    if (normalizer.gt.c0) frac(k,:) = frac(k,:)/normalizer
                  end do

                  do n = 1, ncat
                    if ((aicen(i,j,n,iblk).gt.puny).and.&
                       (SUM(trcrn(i,j,nt_fsd+1:nt_fsd+nfsd-1,n,iblk)).gt.puny).and.&
                       (.NOT.(trcrn(i,j,nt_fsd+1,n,iblk).ge.c1))) then
                       
                        amfstd_init=trcrn(i,j,nt_fsd:nt_fsd+nfsd-1,n,iblk)

                        if (ABS(SUM(amfstd_init)-c1).gt.puny) stop &
                                'init mFSTD not norm, wave'
                        
                        ! protect against small numerical errors
                        WHERE (amfstd_init.lt.puny) amfstd_init = c0
                        amfstd_init = amfstd_init / SUM(amfstd_init)

                        amfstd_tmp =  amfstd_init

                        ! adaptive sub-timestep
                        elapsed_t = c0
                        DO WHILE (elapsed_t.lt.dt)

                             ! calculate d_amfstd using current afstd
                             d_amfstd_tmp = get_damfstd_wave(amfstd_tmp, fracture_hist, frac, fracbreak)
                             WHERE (ABS(d_amfstd_tmp).lt.puny) d_amfstd_tmp = c0 
 
                             ! ensure timestep does not allow afstd to have any negatives
			     check_dt(:) = bignum 
			     do k = 1, nfsd
			        if (d_amfstd_tmp(k).gt.puny) check_dt(k) = (c1-amfstd_tmp(k))/d_amfstd_tmp(k)
				if (d_amfstd_tmp(k).lt.-puny) check_dt(k) = -amfstd_tmp(k)/d_amfstd_tmp(k)
		             end do 
                             subdt = MINVAL(check_dt)
                             subdt = MIN(subdt, dt-elapsed_t)

                             ! update amfstd and elpased time
                             amfstd_tmp = amfstd_tmp + subdt * d_amfstd_tmp
!                write(*,"(15F7.4)") subdt, d_amfstd_tmp

                             ! check conservation and negatives
                             if (ANY(amfstd_tmp.lt.-puny)) stop &
                                     'wb, <0 in loop'

                             if (ANY(amfstd_tmp.gt.c1+puny)) stop &
                                     'wb, >1 in loop'

! am getting n=1 cat fill up with nan though
! adding this write statement eliminated the NaN but still have times when subdt is reduced
!          write(*,10)'fracbraak*', i,j,n,iblk,elapsed_t,subdt, fracbreak, amfstd_tmp, d_amfstd_tmp

                             ! in case of small numerical errors
                             WHERE (amfstd_tmp.lt.puny) amfstd_tmp = c0
                             amfstd_tmp = amfstd_tmp/SUM(amfstd_tmp)

                            ! update time
                             elapsed_t = elapsed_t + subdt 

                        END DO

!          write(*,10)'fracbreak*', i,j,n,iblk,elapsed_t,subdt, fracbreak, amfstd_tmp, amfstd_init
!10       format(1X,A9,2x,4(I2,1x),2x,2(F10.1,1x),2x,50(F8.5,1x))
!         call flush_fileunit(6)

                       ! was already normalized in loop
                        trcrn(i,j,nt_fsd:nt_fsd+nfsd-1,n,iblk) = amfstd_tmp
     
                        if (ABS(SUM(trcrn(i,j,nt_fsd:nt_fsd+nfsd-1,n,iblk))-c1).gt.puny) stop 'not 1 wb'

                        ! for diagnostics
                        d_amfstd_wave(i,j,:,n,iblk) = trcrn(i,j,nt_fsd:nt_fsd+nfsd-1,n,iblk) - amfstd_init
                        d_afsd_wave(i,j,:,iblk) = d_afsd_wave(i,j,:,iblk) + aicen(i,j,n,iblk)*d_amfstd_wave(i,j,:,n,iblk)

                    end if ! aicen>puny
                    end do ! n

                end if ! not all frac zero
            end if ! aice>p01

         end do ! i
         end do  !j

         if (ANY(trcrn(ilo:ihi,jlo:jhi,nt_fsd:nt_fsd+nfsd-1,:,iblk).lt.c0)) stop 'wavebreak neg wb'
         if (ANY(trcrn(ilo:ihi,jlo:jhi,nt_fsd:nt_fsd+nfsd-1,:,iblk).gt.c1)) stop 'wavebreak >1 wb'

     end do ! iblk
     !$OMP END PARALLEL DO

     end subroutine wave_frac_fsd


!=======================================================================
! Author: Lettie Roach, NIWA, 2018
!
! Based on MatLab code from Horvat & Tziperman (2015). Calculates functions 
! to describe the change in the FSD when waves fracture ice, given a wave
! spectrum (1D frequency, 25 frequency bins) in ice. 
! We calculate extrema and if these are successive maximum, minimum, maximum or
! vice versa, and have strain greater than a critical strain, break ice and 
! create new floes with lengths equal to these distances

     subroutine wave_frac(hbar, spec_efreq,frac_local,fracbreak)

     use ice_communicate, only: my_task, master_task
     use ice_fsd, only: floe_rad_l, floe_binwidth, floe_rad_c, floe_rad_h
     use ice_domain_size, only: nfsd
     use ice_flux, only: freq, dfreq
 
     real (kind=dbl_kind),  intent(in) :: &
         hbar   ! mean ice thickness

     real (kind=dbl_kind), dimension(nfreqww3), intent(in) :: &
         spec_efreq
 
     real (kind=dbl_kind),  intent(out) :: &
         fracbreak                  ! how much of cell experiences breaking

     real (kind=dbl_kind), dimension (nfsd), intent(out) :: &
         frac_local

     ! local variables

     integer (kind=int_kind) :: i, j, k

     real (kind=dbl_kind), dimension(nfreqww3) :: &
         reverse_spec_elambda, &
         reverse_lambda,       & ! wavelengths (m) in decreasing order
         reverse_dlambda,      &
         spec_coeff, &
         phi, rand_array, summand

     real (kind=dbl_kind) ::  normalizer

! this needs to be 2*nx/threshold
     real (kind=dbl_kind), dimension(20000) :: &
         fixfraclengths

! this needs to be nx/threshold
     real (kind=dbl_kind), dimension(nx) :: &
         fixfraclengths2

     real (kind=dbl_kind), dimension(nx) :: &
          Xsub, eta

     real (kind=dbl_kind), dimension(nfsd) :: &
          frachistogram 

     logical (kind=log_kind) :: &
           e_stop          ! if true, stop and return zero omega and fsdformed


!!!!!!!!
        ! local to old get_fraclengths subroutine
        integer (kind=int_kind) :: &
                spcing, &       ! distance in dx over which to search for extrema on each side of point
                first, last, &  ! indices over which to search for extrema
                j_neg, &        ! nearest extrema backwards
                j_pos, &        ! nearest extrema forwards
                nlengths        ! number of points where strain is above
                                ! critical strain
        

         real (kind=dbl_kind), dimension(nx) :: &
                strain          ! the strain between triplets of extrema

        logical (kind=log_kind), dimension(nx) :: &
                is_max, is_min, &       ! arrays to hold whether each point is a local max or min
                is_extremum, &          ! or extremum
                is_triplet              ! or triplet of extrema

         real (kind=dbl_kind) :: &
                delta, &        ! difference in x between current and prev extrema
                delta_pos       ! difference in x between next and current extrema

       integer (kind=int_kind), dimension(1) :: &
        maxj, minj  ! indices of local max and min
!!!!!!!!

     ! switch that aborts fracture calc if true
     e_stop=.false.
 
     ! dispersion relation
     reverse_lambda(:) = gravit/(c2*pi*freq(:)**c2)   ! CMB lambda in decreasing order
     reverse_dlambda(:) = gravit/(c2*pi*dfreq(:)**c2)

     ! convert to lambda spectrum
     reverse_spec_elambda(:) = spec_efreq(:) * &
                          (p5 * (gravit/(c2*pi*reverse_lambda(:)**c3) )**p5 )    
 
 
     ! spectral coefficients
     spec_coeff = sqrt(c2*reverse_spec_elambda*reverse_dlambda)
     if (ANY(ISNAN(spec_coeff))) stop 'NaN spec_coeff' 

         ! random phase for each Fourier component
         ! varies in each j loop
         call RANDOM_NUMBER(rand_array)
         phi = c2*pi*rand_array
!	 write(*,*) 'phi',phi
 
         do j=1,nx
             ! spatial domain
             Xsub(j)= j*dx
             !SSH field in space (sum over wavelengths, no attenuation)
             summand = spec_coeff*COS(2*pi*Xsub(j)/reverse_lambda+phi)
             eta(j)=SUM(summand)
         end do
         
    if ((.NOT.(ALL(eta.eq.c0))).and.(hbar.gt.puny)) then 

! CMB moved the subroutine here and removed allocating arrays, etc
!            call get_fraclengths(eta, fraclengths, hbar , e_stop)
!!!!!!!!!!

       ! -------equivalent of peakfinder2
       !given eta and spcing, compute extremelocs in ascending order
       spcing=nint(threshold/dx)

       is_max = .false.
       is_min = .false.
       is_extremum = .false.
       is_triplet = .false.
       strain = c0
       j_neg = 0
       j_pos = 0      
 
       ! search for local max and min within spacing of
       ! 10m on either side of each point

       nlengths = 0
       do j = 1, nx

                first = MAX(1,j-spcing)
                last = MIN(nx,j+spcing)


                maxj = MAXLOC(eta(first:last))
                minj = MINLOC(eta(first:last))

                if (COUNT(eta(first:last).eq.MAXVAL(eta(first:last))).gt.1) &
                        stop 'more than one max'
                if (COUNT(eta(first:last).eq.MINVAL(eta(first:last))).gt.1) &
                        stop 'more than one min'

                if (maxj(1)+first-1.eq.j) is_max(j) = .true.
                if (minj(1)+first-1.eq.j) is_min(j) = .true.

                if (is_min(j).and.is_max(j)) then
                         print *, 'X ',Xsub
                         print *, 'eta ',eta
                         print *, 'frst last' ,first, last
                         print *, 'maxj, minj ',maxj,minj
                         stop &
                        'error in extrema'
                end if
                if (is_min(j).or.is_max(j)) is_extremum(j) = .true.
        end do

        do j = 2, nx-1
                if (is_extremum(j)) then
                        if (j.eq.2) then
                                if (is_extremum(1)) j_neg = 1
                        else
                            do k = j-1, 1, -1
                                if (is_extremum(k)) then
                                        j_neg = k
                                        EXIT
                                end if
                             end do
                        end if
                       
                        do k = j+1, nx
                                if (is_extremum(k)) then
                                        j_pos = k
                                        EXIT
                                end if
                        end do
                       
                        if ((j_neg.gt.0).and.(j_pos.gt.0)) then 
                            if (is_max(j_neg).and.is_min(j).and.is_max(j_pos)) &
                                is_triplet(j) = .true.
                            if (is_min(j_neg).and.is_max(j).and.is_min(j_pos)) &
                                is_triplet(j) = .true.
                        end if

                        ! calculate strain
                        if (is_triplet(j)) then
                                delta_pos = Xsub(j_pos) - Xsub(j)
                                delta = Xsub(j) - Xsub(j_neg)
                                strain(j) = (hbar/c2) *(eta(j_neg)*delta_pos - &
                                        eta(j)*(delta_pos+delta) + eta(j)*delta ) &
                                        /(delta*delta_pos*(delta+delta_pos))


                                if (strain(j).gt.straincrit) then
					nlengths=nlengths+1
            				fixfraclengths(nlengths)=delta_pos
  					fixfraclengths2(nlengths)=delta

                                end if
                        end if

                end if

        end do

        fixfraclengths((nlengths+1):(2*nlengths)) = fixfraclengths2(1:nlengths)
	nlengths = 2*nlengths

     else
       nlengths = 1
       fixfraclengths(1) = c0
     end if

!!!!!!!!!!

     fracbreak=min(1.0,SUM(fixfraclengths(1:nlengths))/D)

     if (nlengths.eq.1) e_stop = .true.  ! do not bother to break
     if (fracbreak.le.puny) e_stop = .true.       ! do not bother to break 

!         write(*,20) 'subroutin',nlengths,sum(fixfraclengths(1:nlengths)),fracbreak
!20       format(1X,A9,2x,I5,2x,F12.1,5x,F8.3)


     if (.NOT. e_stop ) then

        frachistogram(:)=c0
        ! convert from diameter to radii
        fixfraclengths(:)=fixfraclengths(:)/c2
        
        ! bin into FS cats
        ! highest cat cannot be fractured into
        do j=1,nlengths
           do k=1,nfsd-1
              if ((fixfraclengths(j).ge.floe_rad_l(k)).and.(fixfraclengths(j).lt.floe_rad_l(k+1))) then
                 frachistogram(k)=frachistogram(k)+1
              end if
           end do
!         if (fixfraclengths(j).lt.(floe_rad_l(1)-puny)) stop 'fractures too small'
        end do
           
            !if (my_task.eq.32) print *, 'got fraclengths'

        frac_local(:)=floe_rad_c(:)*frachistogram(:)

        ! normalize to one (not D/2, as it is not needed when leaving off twice
        normalizer = SUM(frac_local) !*c2/D       ! Roach et al 2018 Eq 21
        if (normalizer.ne.c0) frac_local(:) = frac_local(:) / normalizer
     else
        frac_local(:)=c0
     end if


!         write(*,30) 'wave_frac',fracbreak,frac_local
!30       format(1X,A9,2x,F8.3,5x,12(F8.3,2x))

     end subroutine wave_frac
 
!=======================================================================
     
      end module ice_wavefracspec

!=======================================================================


