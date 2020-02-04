!
!  This module contains the subroutines required to fracture sea ice
!  by ocean surface waves
!
!  Theory based on:
!
!    Horvat, C., & Tziperman, E. (2015). A prognostic model of the sea-ice 
!    floe size and thickness distribution. The Cryosphere, 9(6), 2119–2134.
!    doi:10.5194/tc-9-2119-2015
!
!  and implementation described in:
!
!    Roach, L. A., Horvat, C., Dean, S. M., & Bitz, C. M. (2018). An emergent
!    sea ice floe size distribution in a global coupled ocean--sea ice model. 
!    Journal of Geophysical Research: Oceans, 123(6), 4322–4337. 
!    doi:10.1029/2017JC013692
!
!  now with some modifications to allow direct input of ocean surface wave spectrum.
!
!  We calculate the fractures that would occur if waves enter a fully ice-covered 
!  region defined in one dimension in the direction of propagation, and then apply 
!  the outcome proportionally to the ice-covered fraction in each grid cell. Assuming
!  that sea ice flexes with the sea surface height field, strains are computed on this
!  sub-grid-scale 1D domain. If the strain between successive extrema exceeds a critical
!  value new floes are formed with diameters equal to the distance between the extrema.
!
!  authors: 2016-8 Lettie Roach, NIWA/VUW
!
      module ice_wavefracspec

      use ice_kinds_mod
      use ice_constants
      use ice_flux, only: dfreq
      use ice_domain_size, only: nfsd, ncat, nfreq, &
                                 max_ntrcr, max_blocks
      use ice_fsd, only: icepack_cleanup_fsdn  
 
      implicit none
      private
      public :: wave_frac, icepack_wavefracfsd, get_subdt_fsd 



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
 
      logical (kind=log_kind), public :: &
         wave_spec      

!=======================================================================

      contains

!=======================================================================
!
!  Calculate the change in the FSD arising from wave fracture
!
!  authors: 2017 Lettie Roach, NIWA/VUW
!
     function get_damfstd_wave(amfstd_init, fracture_hist, frac) &
                              result(d_amfstd)

      real (kind=dbl_kind), dimension (nfsd), intent(in) :: &
         amfstd_init, fracture_hist

     real (kind=dbl_kind), dimension (nfsd,nfsd), intent(in) :: &
         frac

     real (kind=dbl_kind), dimension (nfsd) :: &
         d_amfstd, loss, gain, omega
 
      integer (kind=int_kind) :: k



      do k = 1,nfsd
           ! fracture_hist is already normalized
           omega(k) = amfstd_init(k)*SUM(fracture_hist(1:k-1)) 
      end do

      if (SUM(omega).gt.c1+puny) stop &
          'omega cannot be greater than 1, waves'
                        
      loss = omega

      do k =1,nfsd
           gain(k) = SUM(omega*frac(:,k)) 
      end do

      if (gain(nfsd).gt.puny) stop 'largest cat cannot gain, waves'
      if (loss(1).gt.puny) stop 'smallest cat cannot lose, waves'

      d_amfstd(:) = gain(:) - loss(:)
 
      if (ABS(SUM(d_amfstd(:))).gt.puny) stop 'area not cons, waves'

      WHERE (ABS(d_amfstd).lt.puny) d_amfstd = c0 

      end  function get_damfstd_wave

!=======================================================================
!
!  Adaptive timestepping for FSD, process-agnostic
!  See reference: Horvat & Tziperman (2017) JGR, Appendix A
!
!  authors: 2018 Lettie Roach, NIWA/VUW
!
!
     function get_subdt_fsd(ndim, amfstd_init, d_amfstd) &
                              result(subdt)

      integer (kind=int_kind), intent(in) :: ndim

      real (kind=dbl_kind), dimension (ndim), intent(in) :: &
         amfstd_init, d_amfstd

      real (kind=dbl_kind), dimension (ndim) :: &
         check_dt 
 
      integer (kind=int_kind) :: k

      real (kind=dbl_kind) :: subdt

      check_dt(:) = bignum 
      do k = 1, ndim
          if (d_amfstd(k).gt.puny) check_dt(k) = (c1-amfstd_init(k))/d_amfstd(k)
          if (d_amfstd(k).lt.-puny) check_dt(k) = amfstd_init(k)/ABS(d_amfstd(k))
      end do 
            
      subdt = MINVAL(check_dt)



      end function get_subdt_fsd

!=======================================================================
! Author: Lettie Roach, NIWA, 2019
! 
! Wrapper routine for wave_frac_fsd

     subroutine icepack_wavefracfsd

     use ice_domain, only: nblocks, blocks_ice
     use ice_blocks, only: block, get_block
     use ice_state, only: nt_fsd, aice, aicen, vice, &
                          trcrn
     use ice_flux, only: wave_spectrum, wave_hs_in_ice
     use ice_fsd, only: d_afsd_wave, d_amfstd_wave
 
     ! local variables
      type (block) :: &
         this_block      ! block information for current block

      integer (kind=int_kind) :: &  
        iblk, i, j, &
        ilo, ihi, jlo, jhi ! beginning and end of physical domain

     !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block)
     do iblk = 1, nblocks

         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
         do i = ilo, ihi

            d_afsd_wave(i,j,:,iblk) = c0
            ! LR this condition is FOR TESTING ONLY when using dummy wave spectrum
            ! do not use for actual runs!!
            !if (aice(i,j,iblk).lt.0.8_dbl_kind) &

            call wave_frac_fsd(aice(i,j,iblk),  vice(i,j,iblk),        & ! in 
                               aicen(i,j,:,iblk),                      & ! in
                               wave_spectrum(i,j,:,iblk),              & ! in
                               trcrn(i,j,nt_fsd:nt_fsd+nfsd-1,:,iblk), & ! inout
                               wave_hs_in_ice(i,j,iblk),               & ! out
                               d_afsd_wave(i,j,:,iblk),                & ! out
                               d_amfstd_wave(i,j,:,:,iblk)             ) ! out
            
         end do ! i
         end do  !j

     end do ! iblk
     !$OMP END PARALLEL DO


     end subroutine icepack_wavefracfsd



!=======================================================================
! 
!  Given fracture histogram computed from local wave spectrum, evolve 
!  the floe size distribution
!
!  authors: 2018 Lettie Roach, NIWA/VUW
!


     subroutine wave_frac_fsd(aice,  vice,        & ! in 
                              aicen,              & ! in
                              wave_spectrum,      & ! in
                              trcrn,              & ! inout
                              wave_hs_in_ice,     & ! out
                              d_afsd_wave,        & ! out
                              d_amfstd_wave       ) ! out
            


     use ice_calendar, only: dt

     real (kind=dbl_kind), intent(in) :: &
         aice, &
         vice

     real (kind=dbl_kind), dimension(ncat), intent(in) :: &
         aicen

     real (kind=dbl_kind), dimension(nfreq), intent(inout) :: &
         wave_spectrum

     real (kind=dbl_kind), dimension(nfsd,ncat), intent(inout) :: &
         trcrn

     real (kind=dbl_kind), intent(out) :: &
         wave_hs_in_ice

     real (kind=dbl_kind), dimension(nfsd), intent(out) :: &
         d_afsd_wave

      real (kind=dbl_kind), dimension(nfsd,ncat), intent(out) :: &
         d_amfstd_wave 

     ! local variables
      integer (kind=int_kind) :: &  
        n, k, ks, t, &
        nsubt ! number of subcycles 

      real (kind=dbl_kind), dimension (nfsd) :: &
        fracture_hist

      real (kind=dbl_kind), dimension (nfsd, nfsd) :: &
        frac    

      real (kind=dbl_kind) :: &
        hbar, elapsed_t, subdt, cons_error

      real (kind=dbl_kind), dimension (nfsd) :: &
           amfstd_init     ! tracer array

      real (kind=dbl_kind), dimension (nfsd) :: &
         amfstd_tmp, d_amfstd_tmp, &
         stability_test, &   ! check timestep stability
         omega, &
         check_dt, &
         gain, loss, &
         tempfracs

    !------------------------------------

    ! initialize 
    d_afsd_wave(:) = c0
    d_amfstd_wave(:,:) = c0
    wave_hs_in_ice = c0
    fracture_hist(:) = c0

    ! sanity check
    if (ANY(trcrn(:,:).lt.c0)) stop 'neg b4-wb'
   
    ! temporary fix
    if (ANY(wave_spectrum(:).gt.1.0e+10_dbl_kind)) wave_spectrum(:) = c0

    ! do not try to fracture for minimal ice concentration or zero wave spectrum
    if ((aice.gt.p01).and.(.NOT.(ALL(wave_spectrum(:).lt.puny)))) then

        ! save for diagnostics
        wave_hs_in_ice = c4*SQRT(SUM(wave_spectrum(:)*dfreq(:)))

        hbar = vice / aice

        ! calculate fracture histogram
        call  wave_frac(hbar, wave_spectrum(:), fracture_hist)

        ! if fracture occurs
        if (MAXVAL(fracture_hist).gt.puny) then

            do n = 1, ncat
               

                call icepack_cleanup_fsdn(nfsd, trcrn(:,n))

                if ((aicen(n).gt.puny).and.&
                    (SUM(trcrn(:,n)).gt.puny).and.&
                    (trcrn(1,n).lt.c1-puny)) then
                           
                    ! for diagnostics
                    amfstd_init(:) = trcrn(:,n)
                    amfstd_tmp =  amfstd_init

                    ! frac does not vary within subcycle
                    frac(:,:) = c0
                    do k = 2, nfsd
                            frac(k,:k-1) = fracture_hist(:k-1)
                            frac(k,k:) = c0 
                    end do
                    
                    do ks=1,nfsd
                        if (SUM(frac(ks,:)).gt.c0) frac(ks,:) = frac(ks,:)/SUM(frac(ks,:))
                    end do

                    ! adaptive sub-timestep
                    elapsed_t = c0
                    cons_error = c0
                    nsubt = 0
                    DO WHILE (elapsed_t.lt.dt)
                         nsubt = nsubt + 1

                         ! calculate d_amfstd using current afstd
                         d_amfstd_tmp = get_damfstd_wave(amfstd_tmp, fracture_hist, frac)

                         if (ALL(ABS(d_amfstd_tmp).lt.puny)) EXIT
                         if (amfstd_tmp(1).ge.c1-puny) EXIT

                         if (nsubt.gt.100) then
                             print *, 'wave frac not converging'
                             print *, 'd_amfstd_tmp',d_amfstd_tmp
                             print *, 'afsd ',amfstd_tmp
                             print *, 'frac ',frac
                             print *, 'fracture_hist ',fracture_hist
                             stop
                         end if


 
                         ! timestep required for this
                         subdt = get_subdt_fsd(nfsd, amfstd_tmp, d_amfstd_tmp)
                         subdt = MIN(subdt, dt) ! cannot be greater than CICE timestep

                         ! update amfstd and elpased time
                         amfstd_tmp = amfstd_tmp + subdt * d_amfstd_tmp(:)

                         ! check conservation and negatives
                         if (ANY(amfstd_tmp.lt.-puny)) then 
                               print *, 'afsd ',amfstd_tmp
                               print *, 'init afsd',amfstd_tmp - subdt * d_amfstd_tmp(:)
                               print *, 'subdt ',subdt
                               print *, subdt * d_amfstd_tmp(:)
                               stop &
                                 'wb, <0 in loop'
                         end if

                         if (ANY(amfstd_tmp.gt.c1+puny)) stop &
                                 'wb, >1 in loop'
                         ! update time
                         elapsed_t = elapsed_t + subdt 

                    END DO

                    ! In some cases---particularly for strong fracturing---the equation 
                    ! for wave fracture does not quite conserve area. With this test wave
                    ! forcing, the area conservation error is usually less than 10^-8.
                    ! Simply renormalizing may cause the first floe size category to reduce,
                    ! which is not physically allowed to happen. So as a rather blunt fix,
                    ! we adjust the largest floe size category possible to account for the
                    ! tiny extra area.
                    cons_error = SUM(amfstd_tmp) - c1
                    if (ABS(cons_error).gt.1.0e-8_dbl_kind) print *, 'Area conservation error, waves ',cons_error

                  ! area loss: add to first category
                  if (cons_error.lt.c0) then
                      amfstd_tmp(1) = amfstd_tmp(1) - cons_error
                  else
                  ! area gain: take it from the largest possible category 
 
                    do k = nfsd, 1, -1
                        if (amfstd_tmp(k).gt.cons_error) then
                            amfstd_tmp(k) = amfstd_tmp(k) - cons_error
                            EXIT
                        end if
                    end do
                  end if

                  ! update trcrn    
                  trcrn(:,n) = amfstd_tmp
                  call icepack_cleanup_fsdn(nfsd, trcrn(:,n))


                                        ! sanity checks
                  if (ANY(trcrn(:,n).lt.c0-puny)) stop 'neg wb'
                  if (ANY(trcrn(:,n).gt.c1+puny)) stop '>1 wb'

                  ! for diagnostics
                  d_amfstd_wave(:,n) = trcrn(:,n) - amfstd_init(:)  
                  d_afsd_wave(:) =  d_afsd_wave(:) + aicen(n)*d_amfstd_wave(:,n)


                end if ! aicen>puny
            end do ! n

        end if ! not all frac zero
    end if ! aice>p01
        
     end subroutine wave_frac_fsd


!=======================================================================
!
!  Calculates functions to describe the change in the FSD when waves 
!  fracture ice, given a wave spectrum (1D frequency, 25 frequency bins)
!  in ice. We calculate extrema and if these are successive maximum, 
!  minimum, maximum or vice versa, and have strain greater than a 
!  critical strain, break ice and create new floes with lengths equal
!  to these distances. Based on MatLab code written by Chris Horvat,
!  from Horvat & Tziperman (2015). 
!
!  Note that a realization of sea surface height requires a random phase.
!
!  authors: 2018 Lettie Roach, NIWA/VUW

     subroutine wave_frac(hbar, spec_efreq, frac_local)

     use ice_fsd, only: floe_rad_l, floe_rad_c
     use ice_domain_size, only: nfsd
     use ice_flux, only: freq, dfreq
     use ice_exit, only: abort_ice
 
     real (kind=dbl_kind),  intent(in) :: &
         hbar   ! mean ice thickness
! TEMP INOUT
     real (kind=dbl_kind), dimension(nfreq), intent(in) :: &
         spec_efreq
 
     real (kind=dbl_kind), dimension (nfsd), intent(out) :: &
         frac_local

     ! local variables

     integer (kind=int_kind) :: jj, k

     real (kind=dbl_kind), dimension(25) :: &
         lambda,                   &! wavelengths (m)
         spec_coeff, &
         phi, rand_array, summand

     real (kind=dbl_kind), dimension(nx) :: &
         fraclengths

      real (kind=dbl_kind), dimension(nx) :: &
         X,  &    ! spatial domain (m)
         eta      ! sea surface height field (m)

      real (kind=dbl_kind), dimension(nfsd) :: &
         frachistogram 

      ! initialize fracture lengths
      fraclengths(:) = c0
      frachistogram(:) = c0
      frac_local(:) = c0
 
      ! spatial domain
      do jj = 1, nx
         X(jj)= jj*dx
      end do

      ! dispersion relation
      lambda (:) = gravit/(c2*pi*freq (:)**2)

      ! spectral coefficients
      spec_coeff = sqrt(c2*spec_efreq*dfreq) 

      ! Fixed phase - not converging waves here 
      phi = pi

      do jj = 1, nx
         ! SSH field in space (sum over wavelengths, no attenuation)
         summand = spec_coeff*COS(2*pi*X(jj)/lambda+phi)
         eta(jj)  = SUM(summand)
      end do

      if ((SUM(ABS(eta)) > puny).and.(hbar > puny)) then 
         call get_fraclengths(X, eta, fraclengths, hbar)
      end if
  
      if (ANY(fraclengths.gt.floe_rad_l(1))) then

          ! convert from diameter to radii
          fraclengths(:) = fraclengths(:)/c2

          ! bin into FS cats
          do jj = 1, size(fraclengths)
          do k = 1, nfsd-1
             if ((fraclengths(jj) >= floe_rad_l(k)) .and. &
             (fraclengths(jj) < floe_rad_l(k+1))) then
             frachistogram(k) = frachistogram(k) + 1
             end if
          end do
          if (fraclengths(jj)>floe_rad_l(nfsd)) frachistogram(nfsd) = frachistogram(nfsd) + 1
          end do

          do k = 1, nfsd
             frac_local(k) = floe_rad_c(k)*frachistogram(k)
          end do

          ! normalize
          if (SUM(frac_local).gt.puny) frac_local(:) = frac_local(:) / SUM(frac_local(:))

      end if


      end subroutine wave_frac

!===========================================================================
!
!  Given the (attenuated) sea surface height, find the strain across triplets
!  of max, min, max or min, max, min (local extrema within 10m).
!  If this strain is greater than the  critical strain, ice can fracture
!  and new floes are formed with sizes equal to the distances between
!  extrema. Based on MatLab code written by Chris Horvat,
!  from Horvat & Tziperman (2015). 
!
!  authors: 2016 Lettie Roach, NIWA/VUW
!
      subroutine get_fraclengths(X, eta, fraclengths, hbar)

      real (kind=dbl_kind) :: &
         hbar             ! mean thickness (m)

      real (kind=dbl_kind), intent(in), dimension (nx) :: &
         X, &              ! spatial domain (m)
         eta               ! sea surface height field (m)

      real (kind=dbl_kind), intent(inout), dimension (nx) :: &
         fraclengths      ! The distances between fracture points
                          ! Size cannot be greater than nx.
                          ! In practice, will be least

      ! local variables
      integer (kind=int_kind) :: &
         spcing,        & ! distance over which to search for extrema on each side of point
         jj, k,          & ! indices to iterate over domain
         first, last,   & ! indices over which to search for extrema
         j_neg,         & ! nearest extrema backwards
         j_pos,         & ! nearest extrema forwards
         n_above          ! number of points where strain is above critical

      real (kind=dbl_kind), dimension(nx) :: &
         fracdistances, & ! distances in space where fracture has occurred 
         strain           ! the strain between triplets of extrema

      logical (kind=log_kind), dimension(nx) :: &
         is_max, is_min,& ! arrays to hold whether each point is a local max or min
         is_extremum,   & ! or extremum
         is_triplet       ! or triplet of extrema

      real (kind=dbl_kind) :: &
         delta,         & ! difference in x between current and prev extrema
         delta_pos        ! difference in x between next and current extrema

      integer (kind=int_kind), dimension(1) :: &
         maxj, minj       ! indices of local max and min

      ! ------- equivalent of peakfinder2
      ! given eta and spcing, compute extremelocs in ascending order
      spcing = nint(threshold/dx)

      is_max = .false.
      is_min = .false.
      is_extremum = .false.
      is_triplet = .false.
      strain = c0
      j_neg = 0
      j_pos = 0   
      n_above = 0   
      fraclengths(:) = c0
      fracdistances(:) = c0

      ! search for local max and min within spacing
      ! on either side of each point

      do jj = 1, nx

         ! indices within which to search for local max and min
         first = MAX(1,jj-spcing)
         last  = MIN(nx,jj+spcing)

         ! location of max and min within spacing
         maxj = MAXLOC(eta(first:last))
         minj = MINLOC(eta(first:last))

         ! current j is the max or the min, save it
         if (maxj(1)+first-1 == jj) is_max(jj) = .true.
         if (minj(1)+first-1 == jj) is_min(jj) = .true.

         ! save whether max or min in one array
         if (is_min(jj).or.is_max(jj)) is_extremum(jj) = .true.
      end do

      ! loop over points
      ! nothing can happen at the first or last
      do jj = 2, nx-1
         if (is_extremum(jj)) then
            if (jj == 2) then
               if (is_extremum(1)) j_neg = 1
            else
               do k = jj-1, 1, -1
                  if (is_extremum(k)) then
                     j_neg = k
                     EXIT
                  end if
               end do
            end if
                       
            do k = jj+1, nx
               if (is_extremum(k)) then
                  j_pos = k
                  EXIT
               end if
            end do

            ! find triplets of max and min                       
            if ((j_neg > 0).and.(j_pos > 0)) then 
               if (is_max(j_neg).and.is_min(jj).and.is_max(j_pos)) &
                  is_triplet(jj) = .true.
               if (is_min(j_neg).and.is_max(jj).and.is_min(j_pos)) &
                  is_triplet(jj) = .true.
            end if

            ! calculate strain
            if (is_triplet(jj)) then

               ! finite differences
               delta_pos = X(j_pos) - X(jj    )
               delta     = X(jj    ) - X(j_neg)


               ! NB differs from HT2015 - factor 2 in numerator
               ! and eta(j_pos)
               strain(jj) = ABS(hbar*(eta(j_neg)* delta_pos &
                                - eta(jj    )*(delta_pos+delta) &
                                + eta(j_pos)*           delta) &
                                / (delta*delta_pos*(delta+delta_pos)))

            end if ! is triplet
         end if ! is extremum

      end do ! loop over j

      n_above = COUNT(strain > straincrit)

      ! only do if we have some strains exceeding strain crit
      if (n_above>0) then
          
          k = 0
          do jj = 1, nx
            if (strain(jj) > straincrit) then
              k = k + 1
              fracdistances(k) = X(jj)
            end if
          end do

          do jj = 1, n_above
              fraclengths(jj) = fracdistances(jj+1) - fracdistances(jj)
          end do

          fraclengths(n_above) = c0 ! the last one will be 0 - a distance,
                                      ! so reset back to zero
          if (ANY(fraclengths.lt.c0)) then
              print *, fraclengths
              stop 'fl lt 0'
          end if
          if (ANY(fraclengths.ne.fraclengths)) then
              print *, fraclengths
              stop 'nan fl'
          end if



      end if ! n_above
      

      end subroutine get_fraclengths


!=======================================================================
     
      end module ice_wavefracspec

!=======================================================================


