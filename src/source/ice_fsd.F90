!=========================================================================
!
!  This module contains the subroutines required to define
!  a floe size distribution tracer for sea ice
!
!  Theory based on:
!
!    Horvat, C., & Tziperman, E. (2015). A prognostic model of the sea-ice
!    floe size and thickness distribution. The Cryosphere, 9(6), 2119-2134.
!    doi:10.5194/tc-9-2119-2015
!
!  and implementation described in:
!
!    Roach, L. A., Horvat, C., Dean, S. M., & Bitz, C. M. (2018). An emergent
!    sea ice floe size distribution in a global coupled ocean--sea ice model.
!    Journal of Geophysical Research: Oceans, 123(6), 4322-4337.
!    doi:10.1029/2017JC013692
!
!  with some modifications.
!
!  For floe welding parameter and tensile mode parameter values, see
!
!    Roach, L. A., Smith, M. M., & Dean, S. M. (2018). Quantifying
!    growth of pancake sea ice floes using images from drifting buoys.
!    Journal of Geophysical Research: Oceans, 123(4), 2851-2866.
!    doi: 10.1002/2017JC013693
!
!  Variable naming convention:
!  for k = 1, nfsd and n = 1, ncat
!    afsdn(k,n) = trcrn(:,:,nt_nfsd+k-1,n,:)
!    afsd (k) is per thickness category or averaged over n
!    afs is the associated scalar value for (k,n)
!
!  authors: Lettie Roach, VUW/NIWA
!           C. M. Bitz, UW
!  
!  2016: CMB started
!  2016-8: LR worked on most of it
!
!-----------------------------------------------------------------
 
     module ice_fsd

      use ice_domain_size, only: ncat, nfsd, max_blocks, nfreq
      use ice_blocks, only: nx_block, ny_block
      use ice_state, only: nt_fsd, tr_fsd 
      use ice_kinds_mod
      use ice_constants

      implicit none

      private
      public :: init_fsd_bounds, init_fsd,     &
         fsd_lateral_growth, fsd_add_new_ice, &
          write_restart_fsd, read_restart_fsd, &
          wave_dep_growth, partition_area, &
          icepack_weldfsd, &
          icepack_cleanup_fsd, icepack_cleanup_fsdn, &
          get_subdt_fsd 

      logical (kind=log_kind), public :: & 
         restart_fsd      ! if .true., read fsd tracer restart file
 
      real(kind=dbl_kind), dimension(nfsd), save, public ::  &
         floe_rad_l,    &  ! fsd size lower bound in m (radius)
         floe_rad_c,    &  ! fsd size bin centre in m (radius)
         floe_binwidth     ! fsd size bin width in m (radius)

      ! should this be here?
      real (kind=dbl_kind), parameter,public :: &
         floeshape = 0.66_dbl_kind  ! constant from Steele (unitless)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nfsd,max_blocks), public, save :: &
        d_afsd_latg, d_afsd_latm, d_afsd_newi, d_afsd_weld, d_afsd_wave

     real(kind=dbl_kind), dimension(nfsd) ::  &
         floe_rad_h,    &  ! fsd size higher bound in m (radius)
         floe_area_l,   &  ! fsd area at lower bound (m^2)
         floe_area_h,   &  ! fsd area at higher bound (m^2)
         floe_area_c,   &  ! fsd area at bin centre (m^2)
         floe_area_binwidth ! floe area bin width (m^2)

      integer(kind=int_kind), dimension(nfsd, nfsd) ::  &
         iweld            ! floe size categories that can combine
                          ! during welding (dimensionless)

      logical (kind=log_kind), public :: &
         wave_spec

      ! delete
      integer(kind=int_kind), public :: new_ice_fs


!=======================================================================

      contains

!=======================================================================
!
!  Initialize ice fsd bounds (call whether or not restarting)
!  Define the bounds, midpoints and widths of floe size
!  categories in area and radius
!
!  Note that radius widths cannot be larger than twice previous
!  category width or floe welding will not have an effect
!
!  Note also that the bound of the lowest floe size category is used
!  to define the lead region width and the domain spacing for wave fracture
!
!  authors: Lettie Roach, NIWA/VUW and C. M. Bitz, UW
!
       subroutine init_fsd_bounds

        use ice_domain_size, only: nfsd,ncat
        use ice_constants, only: puny, c2
        use ice_communicate, only: my_task, master_task
        use ice_calendar, only: dt

       integer (kind=int_kind) :: k, a, b, c, m, n

       real (kind = dbl_kind) :: test

                                              
      real (kind=dbl_kind), dimension(:), allocatable :: &
         lims


      if (nfsd.eq.24) then

         allocate(lims(24+1))

         lims = (/ 6.65000000e-02,   5.31030847e+00,   1.42865861e+01,   2.90576686e+01, &
                   5.24122136e+01,   8.78691405e+01,   1.39518470e+02,   2.11635752e+02, &
                   3.08037274e+02,   4.31203059e+02,   5.81277225e+02,   7.55141047e+02, &
                   9.45812834e+02,   1.34354446e+03,   1.82265364e+03,   2.47261361e+03, &
                   3.35434988e+03,   4.55051413e+03,   6.17323164e+03,   8.37461170e+03, &
                   1.13610059e+04,   1.54123510e+04,   2.09084095e+04,   2.83643675e+04, &
                   3.84791270e+04 /)
        
      elseif (nfsd.eq.16) then

         allocate(lims(16+1))

         lims = (/ 6.65000000e-02,   5.31030847e+00,   1.42865861e+01,   2.90576686e+01, &
                   5.24122136e+01,   8.78691405e+01,   1.39518470e+02,   2.11635752e+02, &
                   3.08037274e+02,   4.31203059e+02,   5.81277225e+02,   7.55141047e+02, &
                   9.45812834e+02,   1.34354446e+03,   1.82265364e+03,   2.47261361e+03, &
                   3.35434988e+03 /)
        
      else if (nfsd.eq.12) then

         allocate(lims(12+1))

         lims = (/ 6.65000000e-02,   5.31030847e+00,   1.42865861e+01,   2.90576686e+01, &
                   5.24122136e+01,   8.78691405e+01,   1.39518470e+02,   2.11635752e+02, &
                   3.08037274e+02,   4.31203059e+02,   5.81277225e+02,   7.55141047e+02, &
                   9.45812834e+02 /)
 
      else if (nfsd.eq.1) then ! default case

         allocate(lims(1+1))

         lims = (/ 6.65000000e-02,   3.0e+02 /)

      else
            stop 'floe size categories not defined for this nfsd'
      end if

        floe_rad_l = lims(:nfsd)
        floe_rad_h = lims(2:)
        floe_rad_c = (floe_rad_h+floe_rad_l)/c2

        floe_area_l = 4.*floeshape*floe_rad_l**2.
        floe_area_c = 4.*floeshape*floe_rad_c**2.
        floe_area_h = 4.*floeshape*floe_rad_h**2.

        floe_binwidth = floe_rad_h - floe_rad_l


        
        if (my_task == master_task) then
            write(*,*)&
            'floe size bin info: low, high, center, width, area_c'
            write(*,*) floe_rad_l
            write(*,*) floe_rad_h
            write(*,*) floe_rad_c
            write(*,*) floe_binwidth
            write(*,*) floe_area_l
            write(*,*) floe_area_h   
            write(*,*) floe_area_c
        end if

      
      ! floe size categories that can combine during welding
      iweld(:,:) = -999
      do n = 1, nfsd
      do m = 1, nfsd
         test = floe_area_c(n) + floe_area_c(m)
         do k = 1, nfsd-1
            if ((test >= floe_area_l(k)) .and. (test < floe_area_h(k))) &
                  iweld(n,m) = k
         end do
         if (test >= floe_area_l(nfsd)) iweld(n,m) = nfsd
      end do
      end do

     

      end subroutine init_fsd_bounds

!=======================================================================
!
!  Initialize the FSD 
!
!  When growing from no-ice conditions, initialize to zero.
!  This allows the FSD to emerge, as described in Roach, Horvat et al. (2018)
!
!  Otherwise initalize with a power law, following Perovich
!  & Jones (2014). The basin-wide applicability of such a 
!  prescribed power law has not yet been tested.
!
!  Perovich, D. K., & Jones, K. F. (2014). The seasonal evolution of 
!  sea ice floe size distribution. Journal of Geophysical Research: Oceans,
!  119(12), 8767–8777. doi:10.1002/2014JC010136
!
!  authors: Lettie Roach, NIWA/VUW


      subroutine init_fsd(ice_ic, nx_block, ny_block, iblk, ncat, nfsd, afsdn)
        
        use ice_state, only: nt_fsd
        use ice_domain_size, only: max_ntrcr

   character(len=char_len_long), intent(in) :: &
             ice_ic           ! method of ice cover initialization

   integer(kind=int_kind), intent(in) :: &
             nx_block , &
             ny_block , &
             iblk , &
             ncat, nfsd

   real (kind=dbl_kind), dimension (nx_block,ny_block,nfsd,ncat), &
   intent(inout) :: &
              afsdn     ! tracer array

   real (kind=dbl_kind), parameter :: alpha=2.1_dbl_kind
    
   real (kind=dbl_kind) :: &
              totfrac     ! total fraction of floes 

   integer (kind=int_kind) :: k

   real  (kind=dbl_kind), dimension (nfsd) :: &
             afsd,     &  ! areal distribution of floes
             num_fsd,  &  ! number distribution of floes
             frac_fsd     ! fraction distribution of floes

    do k = 1, nfsd
        afsdn(:,:,k,:) = c0
    enddo


    if (trim(ice_ic).ne.'none') then
   
        ! Perovich (2014)                        
        ! fraction of ice in each floe size and thickness category
        ! same for ALL cells (even where no ice) initially
        totfrac = c0
        do k = 1, nfsd
            num_fsd(k) = (c2*floe_rad_c(k))**(-alpha-c1)
            afsd(k) = num_fsd(k)*floe_area_c(k)*floe_binwidth(k)
            totfrac = totfrac + afsd(k)
        enddo

        do k = 1, nfsd 
            afsdn(:,:,k,:) = afsd(k)/totfrac
        enddo

    endif ! ice_ic

      end subroutine init_fsd

!=======================================================================
!
!  Clean up small values and renormalize
!
!  authors:  Elizabeth Hunke, LANL
!
      subroutine icepack_cleanup_fsd (afsdn)


      real (kind=dbl_kind), dimension(nfsd,ncat), intent(inout) :: &
         afsdn              ! floe size distribution tracer

      ! local variables

      integer (kind=int_kind) :: &
         n                  ! thickness category index



      if (tr_fsd) then

         do n = 1, ncat
            call icepack_cleanup_fsdn(afsdn(:,n))
         enddo

      endif ! tr_fsd

      end subroutine icepack_cleanup_fsd

!=======================================================================
!
!  Clean up small values and renormalize -- per category
!
!  authors:  Elizabeth Hunke, LANL
!

      subroutine icepack_cleanup_fsdn (afsd)


      real (kind=dbl_kind), dimension(nfsd), intent(inout) :: &
         afsd               ! floe size distribution tracer

      ! local variables

      integer (kind=int_kind) :: &
         k                  ! floe size category index

      real (kind=dbl_kind) :: &
         tot

      do k = 1, nfsd
         if (afsd(k) < puny) afsd(k) = c0
      enddo

      tot = sum(afsd(:))
      if (tot > puny) then
         do k = 1, nfsd
            afsd(k) = afsd(k) / tot ! normalize
         enddo
      else ! represents ice-free cell, so set to zero
         afsd(:) = c0
      endif
 
      ! remove
     if (any(afsd < -puny)) stop 'afsd<0 in cleanup'
     if (any(afsd > c1+puny)) stop 'afsd>1 in cleanup'
     if (any(afsd .ne. afsd)) stop 'afsd NaN in cleanup'



      end subroutine icepack_cleanup_fsdn

!=======================================================================
! 
!  Given the joint ice thickness and floe size distribution, calculate
!  the lead region and the total lateral surface area following Horvat
!  and Tziperman (2015).
!
!  authors: Lettie Roach, NIWA/VUW

      subroutine partition_area (        aice,     &
                               aicen,    vicen,    &
                               trcrn,    lead_area,&
                               latsurf_area )

      
      real (kind=dbl_kind), intent(in) :: &
         aice               ! ice concentration

      real (kind=dbl_kind), dimension(ncat), intent(in) :: &
         aicen          , & ! fractional area of ice 
         vicen              ! volume per unit area of ice

      real (kind=dbl_kind), dimension(nfsd,ncat), intent(in) :: &
         trcrn       ! tracer array

      real (kind=dbl_kind), intent(out) :: &
         lead_area      , & ! the fractional area of the lead region
         latsurf_area       ! the area covered by lateral surface of floes

      ! local variables

      integer (kind=int_kind) :: &
         i, j           , & ! horizontal indices
         n              , & ! thickness category index
         k                  ! floe size index

      real (kind=dbl_kind) :: &
         width_leadreg, &   ! width of lead region (m)
         thickness          ! actual thickness of ice in thickness cat (m)

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------
      lead_area    = c0
      latsurf_area = c0

      ! Set the width of the lead region to be the smallest
      ! floe size category, as per Horvat & Tziperman (2015)
      width_leadreg = floe_rad_c(1)

      ! Only calculate these areas where there is some ice
      if (aice > puny) then

         ! lead area = sum of areas of annuli around all floes
         do n = 1, ncat
            do k = 1, nfsd
               lead_area = lead_area + aicen(n) * trcrn(k,n) &
                         * (c2*width_leadreg    /floe_rad_c(k)     &
                             + width_leadreg**2/floe_rad_c(k)**2)
            enddo ! k
         enddo    ! n

         ! cannot be greater than the open water fraction
         lead_area=MIN(lead_area, c1-aice)
         if (lead_area < c0) then
            if (lead_area < -puny) then
!               stop 'lead_area lt0 in partition_area'
            else
               lead_area=MAX(lead_area,c0)
            end if
         end if

         ! Total fractional lateral surface area in each grid (per unit ocean area)
         do n = 1, ncat
            thickness = c0
            if (aicen(n) > c0) thickness = vicen(n)/aicen(n)

            do k = 1, nfsd
               latsurf_area = latsurf_area &
                   + trcrn(k,n) * aicen(n) * c2 * thickness/floe_rad_c(k)
            end do
         end do

         ! check
         if (latsurf_area < c0) stop 'negative latsurf_area'
         if (latsurf_area /= latsurf_area) stop 'latsurf_area NaN'

      end if ! aice

      end subroutine partition_area

!=======================================================================
!
!  Lateral growth at the edges of floes
!
!  Compute the portion of new ice growth that occurs at the edges of
!  floes. The remainder will grow as new ice frazil ice in open water
!  (away from floes).
!
!  See Horvat & Tziperman (2015) and Roach, Horvat et al. (2018).
!
!  authors: Lettie Roach, NIWA/VUW
!
      subroutine fsd_lateral_growth (dt,        aice,         &
                                     aicen,     vicen,        &
                                     vi0new,                  &
                                     frazil,                  &
                                     afsdn,                   &
                                     G_radial,  d_an_latg,    &
                                     tot_latg)


      real (kind=dbl_kind), intent(in) :: &
         dt             , & ! time step (s)
         aice               ! total concentration of ice

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         aicen          , & ! concentration of ice
         vicen              ! volume per unit area of ice          (m)

      real (kind=dbl_kind), dimension(:,:), intent(in) :: &
         afsdn              ! floe size distribution tracer

      real (kind=dbl_kind), intent(inout) :: &
         vi0new         , & ! volume of new ice added to cat 1 (m)
         frazil             ! frazil ice growth        (m/step-->cm/day)

      real (kind=dbl_kind), dimension(ncat), intent(out) :: &
         d_an_latg          ! change in aicen occuring due to lateral growth

      real (kind=dbl_kind), intent(out) :: &
         G_radial       , & ! lateral melt rate (m/s)
         tot_latg           ! total change in aice due to
                            ! lateral growth at the edges of floes

      ! local variables
      integer (kind=int_kind) :: &
         n, k               ! ice category indices

      real (kind=dbl_kind) :: &
         vi0new_lat         ! volume of new ice added laterally to FSD (m)

      real (kind=dbl_kind)  :: &
         lead_area      , & ! the fractional area of the lead region
         latsurf_area       ! the area covered by lateral surface of floes


      lead_area    = c0
      latsurf_area = c0
      G_radial     = c0
      tot_latg     = c0
      d_an_latg    = c0

      ! partition volume into lateral growth and frazil
      call partition_area (aice,      &
                           aicen,      vicen,     &
                           afsdn,      lead_area, &
                           latsurf_area)

      vi0new_lat = c0
      if (latsurf_area > puny) then
         vi0new_lat = vi0new * lead_area / (c1 + aice/latsurf_area)
      end if

      ! for history/diagnostics
      frazil = vi0new - vi0new_lat

      ! lateral growth increment
      if (vi0new_lat > puny) then
         G_radial = vi0new_lat/dt
         do n = 1, ncat

            do k = 1, nfsd
               d_an_latg(n) = d_an_latg(n) &
                            + c2*aicen(n)*afsdn(k,n)*G_radial*dt/floe_rad_c(k)
            end do
         end do ! n
         
         ! cannot expand ice laterally beyond lead region
         if (SUM(d_an_latg(:)).ge.lead_area) then
             d_an_latg(:) = d_an_latg(:)/SUM(d_an_latg(:))
             d_an_latg(:) = d_an_latg(:)*lead_area
         end if

      endif ! vi0new_lat > 0

      ! Use remaining ice volume as in standard model,
      ! but ice cannot grow into the area that has grown laterally
      vi0new = vi0new - vi0new_lat
      tot_latg = SUM(d_an_latg(:))

      end subroutine fsd_lateral_growth

!=======================================================================
!
!  Evolve the FSD subject to lateral growth and the growth of new ice
!  in thickness category 1.
!
!  If ocean surface wave forcing is provided, the floe size of new ice
!  (grown away from floe edges) can computed from the wave field
!  based on the tensile (stretching) stress limitation following
!  Shen et al. (2001). Otherwise, new floes all grow in the smallest
!  floe size category, representing pancake ice formation.
!
!  Shen, H., Ackley, S., & Hopkins, M. (2001). A conceptual model 
!  for pancake-ice formation in a wave field. 
!  Annals of Glaciology, 33, 361-367. doi:10.3189/172756401781818239
!
!  authors: Lettie Roach, NIWA/VUW
!
      subroutine fsd_add_new_ice (n,                         &
                                  dt,         ai0new,        &
                                  d_an_latg,  d_an_newi,     &
                                  G_radial,   area2,         &
                                  wave_sig_ht,               &
                                  wave_spectrum,             &
                                  d_afsd_latg,               &
                                  d_afsd_newi,               &
                                  afsdn,      aicen_init,    &
                                  aicen,      trcrn)

      integer (kind=int_kind), intent(in) :: &
         n              ! thickness category number

      real (kind=dbl_kind), intent(in) :: &
         dt         , & ! time step (s)
         ai0new     , & ! area of new ice added to cat 1
         G_radial   , & ! lateral growth rate (m/s)
         wave_sig_ht    ! wave significant height (everywhere) (m)

      real (kind=dbl_kind), dimension(nfreq), intent(in)  :: &
         wave_spectrum  ! ocean surface wave spectrum as a function of frequency
                        ! power spectral density of surface elevation, E(f) (units m^2 s)

      real (kind=dbl_kind), dimension(ncat), intent(in) :: &
         area2      , & ! area after lateral growth, before new ice formation
         d_an_latg  , & ! change in aicen due to lateral growth
         d_an_newi  , & ! change in aicen due to frazil ice formation
         aicen_init , & ! fractional area of ice
         aicen          ! after update

      real (kind=dbl_kind), dimension (nfsd,ncat), intent(in) :: &
         afsdn          ! floe size distribution tracer

      real (kind=dbl_kind), dimension (nfsd,ncat), intent(inout) :: &
         trcrn          ! ice tracers fsd

      real (kind=dbl_kind), dimension(nfsd), intent(inout) :: &
                        ! change in floe size distribution (area)
         d_afsd_latg, & ! due to fsd lateral growth
         d_afsd_newi    ! new ice formation

      integer (kind=int_kind) :: &
         k, &              ! floe size category index
         new_size, &       ! index for floe size of new ice
         nsubt             ! number of subtimesteps

      real (kind=dbl_kind) :: &
         tmp, elapsed_t, subdt  ! elapsed time, subtimestep (s)

      real (kind=dbl_kind), dimension (nfsd,ncat) :: &
         afsdn_latg     ! fsd after lateral growth

      real (kind=dbl_kind), dimension (nfsd) :: &
         dafsd_tmp,  &  ! tmp FSD
         df_flx     , & ! finite differences for fsd
         afsd_ni        ! fsd after new ice added

      real (kind=dbl_kind), dimension(nfsd+1) :: &
         f_flx          ! finite differences in floe size


      afsdn_latg(:,n) = afsdn(:,n)  ! default

      if (d_an_latg(n) > puny) then ! lateral growth

         ! adaptive timestep
         elapsed_t = c0
         nsubt = 0

         DO WHILE (elapsed_t.lt.dt)
        
             nsubt = nsubt + 1
 
             ! finite differences
             df_flx(:) = c0 ! NB could stay zero if all in largest FS cat
             f_flx (:) = c0
             do k = 2, nfsd
                f_flx(k) = G_radial * afsdn_latg(k-1,n) / floe_binwidth(k-1)
             end do
             do k = 1, nfsd
                df_flx(k) = f_flx(k+1) - f_flx(k)
             end do


             dafsd_tmp(:) = c0
             do k = 1, nfsd
                dafsd_tmp(k) = (-df_flx(k) + c2 * G_radial * afsdn_latg(k,n) &
                            * (c1/floe_rad_c(k) - SUM(afsdn_latg(:,n)/floe_rad_c(:))) )

             end do
             WHERE (abs(dafsd_tmp).lt.puny) dafsd_tmp = c0

            ! timestep required for this
            subdt = get_subdt_fsd(nfsd, afsdn_latg(:,n), dafsd_tmp(:)) 
            subdt = MIN(subdt, dt)


            ! update fsd and elapsed time
            afsdn_latg(:,n) = afsdn_latg(:,n) + subdt*dafsd_tmp(:)
            elapsed_t = elapsed_t + subdt

           if (nsubt.gt.100) print *, 'warning, latg not converging ',nsubt
            
         END DO

         call icepack_cleanup_fsdn (afsdn_latg(:,n))
         trcrn(:,n) = afsdn_latg(:,n)

      end if ! lat growth

      new_size = nfsd
      if (n == 1) then
         ! add new frazil ice to smallest thickness
         if (d_an_newi(n) > puny) then

            afsd_ni(:) = c0

            if (SUM(afsdn_latg(:,n)) > puny) then ! fsd exists

               if (wave_spec) then
                  if (wave_sig_ht > puny) then
                     call wave_dep_growth (wave_spectrum, &
                                           new_size)
                  end if

                  ! grow in new_size category
                  afsd_ni(new_size) = (afsdn_latg(new_size,n)*area2(n) + ai0new) &
                                                          / (area2(n) + ai0new)
                  do k = 1, new_size-1  ! diminish other floe cats accordingly
                     afsd_ni(k) = afsdn_latg(k,n)*area2(n) / (area2(n) + ai0new)
                  end do
                  do k = new_size+1, nfsd  ! diminish other floe cats accordingly
                     afsd_ni(k) = afsdn_latg(k,n)*area2(n) / (area2(n) + ai0new)
                  end do

               else ! grow in smallest floe size category
                  afsd_ni(1) = (afsdn_latg(1,n)*area2(n) + ai0new) &
                                             / (area2(n) + ai0new)
                  do k = 2, nfsd  ! diminish other floe cats accordingly
                     afsd_ni(k) = afsdn_latg(k,n)*area2(n) / (area2(n)+ai0new)
                  enddo
               end if ! wave spec

            else ! no fsd, so entirely new ice

               if (wave_spec) then
                  if (wave_sig_ht > puny) then
                     call wave_dep_growth (wave_spectrum, &
                                           new_size)
                  end if

                  afsd_ni(new_size) = c1
               else
                  afsd_ni(1) = c1
               endif      ! wave forcing

            endif ! entirely new ice or not

            trcrn(:,n) = afsd_ni(:)
            call icepack_cleanup_fsdn (trcrn(:,n))

         endif ! d_an_newi > puny
      endif    ! n = 1

      ! history/diagnostics
      do k = 1, nfsd
         ! sum over n
         d_afsd_latg(k) = d_afsd_latg(k) &
                + area2(n)*afsdn_latg(k,n) & ! after latg
                - aicen_init(n)*afsdn(k,n) ! at start
         d_afsd_newi(k) = d_afsd_newi(k) &
                + aicen(n)*trcrn(k,n) & ! after latg and newi
                - area2(n)*afsdn_latg(k,n) ! after latg
      enddo    ! k

      end subroutine fsd_add_new_ice

!=======================================================================
!
!  Given a wave spectrum, calculate size of new floes based on 
!  tensile failure, following Shen et al. (2001)
!
!  The tensile mode parameter is based on in-situ measurements
!  by Roach, Smith & Dean (2018).
!
!  authors: Lettie Roach, NIWA/VUW
!

      subroutine wave_dep_growth (local_wave_spec, & 
                                       new_size )

       use ice_flux, only: dwavefreq, wavefreq

      real (kind=dbl_kind), dimension(nfreq), intent(in) :: &
         local_wave_spec ! ocean surface wave spectrum as a function of frequency
                         ! power spectral density of surface elevation, E(f) (units m^2 s)
                         ! dimension set in ice_forcing


      integer (kind=int_kind), intent(out) :: &
         new_size        ! index of floe size category in which new floes will grow

      ! local variables
      real (kind=dbl_kind), parameter :: &
         tensile_param = 0.167_dbl_kind ! tensile mode parameter (kg m^-1 s^-2)
                                        ! value from Roach, Smith & Dean (2018)

      real (kind=dbl_kind)  :: &
         w_amp,       & ! wave amplitude (m)
         f_peak,      & ! peak frequency (s^-1)
         r_max          ! floe radius (m)

      integer (kind=int_kind) :: k

      w_amp = c2* SQRT(SUM(local_wave_spec*dwavefreq))   ! sig wave amplitude
      f_peak = wavefreq(MAXLOC(local_wave_spec, DIM=1))  ! peak frequency

      ! tensile failure
      if (w_amp > puny .and. f_peak > puny) then
         r_max = p5*SQRT(tensile_param*gravit/(pi**5*rhoi*w_amp*2))/f_peak**2
      else
         r_max = bignum
      end if

      new_size = nfsd
      do k = nfsd-1, 1, -1
         if (r_max <= floe_rad_h(k)) new_size = k
      end do

      end subroutine wave_dep_growth


!=======================================================================
!
! Wrapper for floe_weld_thermo
!
      subroutine icepack_weldfsd ( &
                               iblk, nx_block, ny_block, &
                               icells, indxi, indxj,     & 
                               dt, aicen, frzmlt,        &
                               trcrn,                    &
                               d_afsd_weld)

      integer (kind=int_kind), intent(in) :: &
         iblk, &
         nx_block, ny_block, & ! block dimensions
         icells                ! number of ice/ocean grid cells

      integer (kind=int_kind), dimension (nx_block*ny_block), &
         intent(in) :: &
         indxi,  indxj         ! compressed i/j indices

     real (kind=dbl_kind), intent(in) :: &
         dt        ! time step (s)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat), &
         intent(in) :: &
         aicen   ! concentration of ice
 
      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         frzmlt    ! freezing/melting potential (W/m^2)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nfsd,ncat), &
         intent(inout) :: &
         trcrn

      real (kind=dbl_kind), dimension (nx_block,ny_block,nfsd), &
         intent(inout) :: &
         d_afsd_weld

      ! local variables

      integer (kind=int_kind) :: &
        ij, i, j

      
    do ij = 1, icells
        
        i = indxi(ij)
        j = indxj(ij)

        call floe_weld_thermo(dt, & 
                              frzmlt(i,j),            & ! in
                              aicen(i,j,:),       & ! in
                              trcrn(i,j,:,:),         & ! inout
                              d_afsd_weld(i,j,:) )

    end do!ij 

    
            end subroutine icepack_weldfsd


!=======================================================================
!
!  Floes are perimitted to weld together in freezing conditions, according
!  to their geometric probability of overlap if placed randomly on the 
!  domain. The rate per unit area c_weld is the total number 
!  of floes that weld with another, per square meter, per unit time, in the 
!  case of a fully covered ice surface (aice=1), equal to twice the reduction
!  in total floe number. See Roach, Smith & Dean (2018).
!
!
!  authors: Lettie Roach, NIWA/VUW
!

      subroutine floe_weld_thermo (dt, frzmlt,  &
                                      aicen, trcrn, &
                                      d_afsd_weld)

                       

      real (kind=dbl_kind), intent(in) :: &
         dt             ! time step (s)

      real (kind=dbl_kind), dimension (ncat), intent(in) :: &
         aicen          ! ice concentration

      real (kind=dbl_kind), intent(in) :: &
         frzmlt         ! freezing/melting potential (W/m^2)

      real (kind=dbl_kind), dimension (nfsd,ncat), intent(inout) :: &
         trcrn          ! ice traceers

      real (kind=dbl_kind), dimension (nfsd), intent(inout) :: &
         d_afsd_weld ! change in fsd due to welding

      ! local variables
      real (kind=dbl_kind), parameter :: &
         aminweld = p1  ! minimum ice concentration likely to weld

      real (kind=dbl_kind), parameter :: &
         c_weld = 1.0e-8_dbl_kind    
                        ! constant of proportionality for welding
                        ! total number of floes that weld with another, per square meter,
                        ! per unit time, in the case of a fully covered ice surface
                        ! units m^-2 s^-1, see documentation for details

      integer (kind=int_kind) :: &
        nt          , & ! time step index
        n           , & ! thickness category index
        k, kx, ky, i, j ! floe size category indices

      real (kind=dbl_kind), dimension(nfsd,ncat) :: &
         afsdn          ! floe size distribution tracer

      real (kind=dbl_kind), dimension(nfsd,ncat) :: &
         d_afsdn_weld   ! change in afsdn due to welding

      real (kind=dbl_kind), dimension(nfsd) :: &
         stability  , & ! check for stability
         nfsd_tmp   , & ! number fsd
         afsd_init  , & ! initial values
         afsd_tmp   , & ! work array
         gain, loss     ! welding tendencies

      real(kind=dbl_kind) :: &
         prefac     , & ! multiplies kernel
         kern       , & ! kernel
         subdt      , & ! subcycling time step for stability (s)
         elapsed_t      ! elapsed subcycling time


      afsdn  (:,:) = c0
      afsd_init(:) = c0
      stability    = c0
      prefac       = p5

      do n = 1, ncat

         d_afsd_weld (:)   = c0
         d_afsdn_weld(:,n) = c0
         afsdn(:,n) = trcrn(:,n)
         call icepack_cleanup_fsdn (afsdn(:,n))

         ! If there is some ice in the lower (nfsd-1) categories
         ! and there is freezing potential
         if ((frzmlt > puny) .and. &               ! freezing potential
             (aicen(n) > aminweld) .and. &         ! low concentrations unlikely to weld
             (SUM(afsdn(1:nfsd-1,n)) > puny)) then ! some ice in nfsd-1 categories

            afsd_init(:) = afsdn(:,n)     ! save initial values
            afsd_tmp (:) = afsd_init(:)   ! work array
               

            ! adaptive sub-timestep
            elapsed_t = c0 
            DO WHILE (elapsed_t < dt) 

               ! calculate sub timestep
               nfsd_tmp = afsd_tmp/floe_area_c
               WHERE (afsd_tmp > puny) &
                  stability = nfsd_tmp/(c_weld*afsd_tmp*aicen(n))
               WHERE (stability < puny) stability = bignum
               subdt = MINVAL(stability)
               subdt = MIN(subdt,dt)

               loss(:) = c0
               gain(:) = c0

               do i = 1, nfsd ! consider loss from this category
               do j = 1, nfsd ! consider all interaction partners
                   k = iweld(i,j) ! product of i+j
                   if (k > i) then
                       kern = c_weld * floe_area_c(i) * aicen(n)
                       loss(i) = loss(i) + kern*afsd_tmp(i)*afsd_tmp(j)
                       if (i.eq.j) prefac = c1 ! otherwise 0.5
                       gain(k) = gain(k) + prefac*kern*afsd_tmp(i)*afsd_tmp(j)
                   end if
               end do
               end do

               ! update afsd   
               afsd_tmp(:) = afsd_tmp(:) + subdt*(gain(:) - loss(:))

               ! in case of minor numerical errors
               WHERE(afsd_tmp < puny) afsd_tmp = c0
               afsd_tmp = afsd_tmp/SUM(afsd_tmp)

               ! update time
               elapsed_t = elapsed_t + subdt

               ! stop if all in largest floe size cat
               if (afsd_tmp(nfsd) > (c1-puny)) exit

            END DO ! time

            call icepack_cleanup_fsdn (afsdn(:,n))
            
            do k = 1, nfsd
               afsdn(k,n) = afsd_tmp(k)
               trcrn(k,n) = afsdn(k,n)
               ! history/diagnostics
               d_afsdn_weld(k,n) = afsdn(k,n) - afsd_init(k)
            enddo
 
         endif ! try to weld
      enddo    ! ncat

      ! history/diagnostics
      do k = 1, nfsd
         d_afsd_weld(k) = c0
         do n = 1, ncat
            d_afsd_weld(k) = d_afsd_weld(k) + aicen(n)*d_afsdn_weld(k,n)
         end do ! n
      end do    ! k


      end subroutine floe_weld_thermo


!=======================================================================
!
!  Adaptive timestepping (process agnostic)
!  See reference: Horvat & Tziperman (2017) JGR, Appendix A
!
!  authors: 2018 Lettie Roach, NIWA/VUW
!
!
      function get_subdt_fsd(nfsd, afsd_init, d_afsd) &
                              result(subdt)

      integer (kind=int_kind), intent(in) :: &
         nfsd       ! number of floe size categories

      real (kind=dbl_kind), dimension (nfsd), intent(in) :: &
         afsd_init, d_afsd ! floe size distribution tracer 

      ! output
      real (kind=dbl_kind) :: &
         subdt ! subcycle timestep (s)

      ! local variables
      real (kind=dbl_kind), dimension (nfsd) :: &
         check_dt ! to compute subcycle timestep (s)

      integer (kind=int_kind) :: k

      check_dt(:) = bignum 
      do k = 1, nfsd
          if (d_afsd(k) >  puny) check_dt(k) = (1-afsd_init(k))/d_afsd(k)
          if (d_afsd(k) < -puny) check_dt(k) = afsd_init(k)/ABS(d_afsd(k))
      end do 

      subdt = MINVAL(check_dt)

      end function get_subdt_fsd

!=======================================================================

! Dumps all values needed for restarting
! author Elizabeth C. Hunke, LANL

      subroutine write_restart_fsd()

      use ice_domain_size, only: ncat,nfsd
      use ice_restart, only: write_restart_field
      use ice_state, only: trcrn, nt_fsd
      use ice_fileunits, only: nu_dump_fsd

      ! local variables

      character*2 ck
      logical (kind=log_kind) :: diag
      integer k

      diag = .true.

      !-----------------------------------------------------------------
      do k=1,nfsd
        write(ck,'(i2.2)') k
        call write_restart_field(nu_dump_fsd,0, trcrn(:,:,nt_fsd+k-1,:,:), &
                            'ruf8','fsd'//'_'//ck,ncat,diag)
      enddo

      end subroutine write_restart_fsd

!=======================================================================

! Reads all values needed for an ice fsd restart
! author Elizabeth C. Hunke, LANL

      subroutine read_restart_fsd()

      use ice_domain_size, only: ncat,nfsd
      use ice_restart, only: read_restart_field
      use ice_state, only: trcrn, nt_fsd
      use ice_fileunits, only: nu_restart_fsd

      ! local variables

      character*2 ck
      logical (kind=log_kind) :: diag
      integer k

      diag = .true.

      do k=1,nfsd
        write(ck,'(i2.2)') k
        call read_restart_field(nu_restart_fsd,0,trcrn(:,:,nt_fsd+k-1,:,:), &
                 'ruf8','fsd'//'_'//ck,ncat,diag, &
                 field_type=field_type_scalar,field_loc=field_loc_center)
      enddo

      end subroutine read_restart_fsd



!=======================================================================

      end module ice_fsd

!=======================================================================

