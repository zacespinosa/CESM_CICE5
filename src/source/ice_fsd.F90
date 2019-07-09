!
! This module contains the subroutines required to define
! a floe size distribution tracer for sea ice
!
! authors: liuxy
!          C. M. Bitz, UW
!          Lettie Roach, NIWA
!
! 2015: liuxy, modified from ice_fsd module
! 2016: CMB rewrote a lot of it
! 2016: LR made some modifications
 
      module ice_fsd

      use ice_domain_size, only: ncat, nfsd, max_blocks
      use ice_blocks, only: nx_block, ny_block
 
      use ice_kinds_mod
      use ice_constants

      implicit none

      private
      public :: init_fsd, init_fsd_bounds,     &
          write_restart_fsd, read_restart_fsd, &
          frzmlt_bottom_lateral_fsd, &
          renorm_mfstd, check_mfstd, wave_dep_growth 

      logical (kind=log_kind), public :: & 
         restart_fsd      ! if .true., read fsd tracer restart file

      logical (kind=log_kind), public :: & 
         write_diag_diff, &  ! if .true., calculate differences in mFSTD and a_n and save to history file
         write_diag_wave     ! if .true., write the lat/lons from find_wave to history file


      real(kind=dbl_kind), dimension(nfsd), save, public ::  &
         floe_rad_l,    &  ! fsd size lower bound in m (radius)
         floe_rad_h,    &  ! fsd size higher bound in m (radius)
         floe_rad_c,    &  ! fsd size bin centre in m (radius)
         floe_binwidth, &  ! fsd size bin width in m (radius)
         floe_area_l,   &  ! fsd area at lower bound (m^2)
         floe_area_h,   &  ! fsd area at higher bound (m^2)
         floe_area_c,   &  ! fsd area at bin centre (m^2)
         floe_area_binwidth, & ! floe area bin width (m^2)
         area_scaled_l, &  ! area bins scaled so they begin at zero
         area_scaled_h, &  ! and no binwidth is greater than 1
         area_scaled_c, &  ! (dimensionless)
         area_scaled_binwidth

      integer(kind=int_kind), dimension(nfsd, nfsd), save, public ::  &
         alpha_mrg

      real (kind=dbl_kind), parameter,public :: &
         floeshape = 0.66_dbl_kind  ! constant from Steele (unitless)

! LR
      real (kind=dbl_kind), dimension (nx_block,ny_block,nfsd,max_blocks), public, save :: &
        d_afsd_latg, d_afsd_latm, d_afsd_addnew, d_afsd_merge, d_afsd_wave

      real (kind=dbl_kind), dimension (nx_block,ny_block,nfsd,ncat,max_blocks), public, save :: &
        d_amfstd_latg, d_amfstd_latm, d_amfstd_addnew, d_amfstd_merge, d_amfstd_wave

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat,max_blocks), public, save :: &
        d_an_latg, d_an_latm, d_an_addnew

      real (kind=dbl_kind), public :: &
        c_mrg            ! constant of proportionality for merging
                         ! see documentation for details

      logical (kind=log_kind), public :: &
         rdc_frzmlt      ! if true, only (1-oo) of frzmlt can melt (lat and bot)

      integer(kind=int_kind), save, public ::  &
         nfreq           ! number of frequencies in wave spectrum   
! LR


!=======================================================================

      contains

!=======================================================================

!  Initialize ice fsd bounds (call whether or not restarting)
!  Define the bounds, midpoints and widths of floe size
!  categories in area and radius
!
!  Note that radius widths cannot be larger than twice previous
!  category width or floe merging will not have an effect
!
!  Note also that the bound of the lowest floe size category is used
!  to define the lead region width and the domain spacing for wave breaking
!
       subroutine init_fsd_bounds

        use ice_domain_size, only: nfsd,ncat
        use ice_constants, only: puny, c2
        use ice_communicate, only: my_task, master_task
        use ice_calendar, only: dt

       integer (kind=int_kind) :: k, a, b, c

       real (kind = dbl_kind) :: test

       real (kind = dbl_kind), dimension (nfsd+1) :: &
         lims, area_lims, area_lims_scaled ! local variables
                                              


      if (nfsd.eq.24) then

        lims =   (/  6.65000000e-02,   5.31030847e+00,   1.42865861e+01,   2.90576686e+01, &
                     5.24122136e+01,   8.78691405e+01,   1.39518470e+02,   2.11635752e+02, &
                     3.08037274e+02,   4.31203059e+02,   5.81277225e+02,   7.55141047e+02, &
                     9.45812834e+02,   1.34354446e+03,   1.82265364e+03,   2.47261361e+03, &
                     3.35434988e+03,   4.55051413e+03,   6.17323164e+03,   8.37461170e+03, &
                     1.13610059e+04,   1.54123510e+04,   2.09084095e+04,   2.83643675e+04, &
                     3.84791270e+04 /)
!                     9.45812834e+02 /) ! was the end of 12 cat bins

        floe_rad_l = lims(:nfsd)
        floe_rad_h = lims(2:)
        floe_rad_c = (floe_rad_h+floe_rad_l)/c2

        floe_area_l = 4.*floeshape*floe_rad_l**2.
        floe_area_c = 4.*floeshape*floe_rad_c**2.
        floe_area_h = 4.*floeshape*floe_rad_h**2.


       else

        do k=1,nfsd
                floe_area_l(k) = c2**k
                floe_area_h(k) = c2**(k+1)
        end do
     
        floe_area_c=(floe_area_l+floe_area_h)/2.

        floe_rad_l=(floe_area_l/(4*floeshape))**.5
        floe_rad_h=(floe_area_h/(4*floeshape))**.5
        floe_rad_c=(floe_area_c/(4*floeshape))**.5

       end if

        floe_binwidth=(floe_rad_h-floe_rad_l)
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


        ! scaled area for merging
        floe_area_binwidth = floe_area_h - floe_area_l
        area_lims(:nfsd) = floe_area_l
        area_lims(nfsd+1) = floe_area_h(nfsd)
        area_lims_scaled = (area_lims - area_lims(1))/MAXVAL(floe_area_binwidth)

        area_scaled_h = area_lims_scaled(2:)
        area_scaled_l = area_lims_scaled(:nfsd)
        area_scaled_c = (area_scaled_h + area_scaled_l) / c2
        area_scaled_binwidth = area_scaled_h - area_scaled_l

        ! which floe sizes can combine during merging
        alpha_mrg(:,:) = - 999
        do a  = 1, nfsd
                do b = 1, nfsd
                        test = area_scaled_h(a) - area_scaled_c(b)

                        do c = 1, nfsd
                                if ((test.ge.area_scaled_l(c)).and.(test.lt.area_scaled_h(c))) then
                                        alpha_mrg(a,b) = c + 1  
                                end if
                        end do
                end do
        end do


      end subroutine init_fsd_bounds

!=======================================================================
!
!  Initialize the FSD as zero everywhere
!  Commented out part - initalize with a power law, following Perovich (2014)
!  Alpha value from Perovich (2014
!

      subroutine init_fsd(nx_block, ny_block, iblk, ncat, nfsd, trcrn)
        
        use ice_state, only: nt_fsd
        use ice_domain_size, only: max_ntrcr

        integer(kind=int_kind), intent(in) :: &
             nx_block , &
             ny_block , &
             iblk , &
             ncat, nfsd

        real (kind=dbl_kind), dimension (nx_block,ny_block,max_ntrcr,ncat), &
         intent(inout) :: &
              trcrn     ! tracer array

	real (kind=dbl_kind) :: alpha, totfrac

        integer (kind=int_kind) :: k

	real  (kind=dbl_kind), dimension (nfsd) :: &
             num_fsd,  &  ! number distribution of floes
             frac_fsd     ! fraction distribution of floes

        do k=1,nfsd
           trcrn(:,:,nt_fsd+k-1,:)=c0
        end do
 
        !alpha=c2+p1
        !
        !do k=1,nfsd
        !   num_fsd(k)=(2*floe_rad_c(k))**(-alpha-c1) 
        !enddo                   
        
        ! total fraction of floes 
        !totfrac=c0

        !do k=1,nfsd
        !   frac_fsd(k)=num_fsd(k)*floe_area_c(k) ! convert to frac from num
        !   totfrac=totfrac+frac_fsd(k)
        !enddo

        ! normalize to one
        !frac_fsd=frac_fsd/totfrac  

        ! fraction of ice in each floe size and thickness category
        ! same for all thickness categories
        ! same for ALL cells (even where no ice) initially
        !do k=1,nfsd
        !   trcrn(:,:,nt_fsd+k-1,:)=frac_fsd(k)
        !end do
 
        !write(*,*)'init_fsd: initial number distribution of floes'
        !write(*,*) num_fsd(1:nfsd)
        !write(*,*)'init_fsd: initial fraction distribution of floes'
        !write(*,*) frac_fsd(1:nfsd)

      end subroutine init_fsd

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
! Modified from frzmlt_bottom_lateral 
! rside becomes a function of itd ncat because 
! floes are not assumed to be same size
! in fact rside is also a function of nfsd but 
! I sum it over the fsd here to get rside_itd
! that can later be used to melt apply to aicen, vicen, etc
! remember that the heat to the ocean is going to come out of variable ?
! in lateral_melt
! a lot of the work here is to simply ensure that no more ice 
! melts than is available
!
! Compute heat flux to bottom surface.
! Compute fraction of ice that melts laterally.
!
! authors C. M. Bitz, UW
!         William H. Lipscomb, LANL
!         Elizabeth C. Hunke, LANL
!         
! 2016: Lettie Roach modified slightly to allow fside and fbot to be 
!       diagnostic output

      subroutine frzmlt_bottom_lateral_fsd (nx_block, ny_block, &
                                        ilo, ihi, jlo, jhi, &
                                        ntrcr,    dt,       &
                                        aice,     aicen,    & 
                                        lead_area,          &
                                        frzmlt,             &
                                        vicen,    vsnon,    &
                                        trcrn,              &
                                        sst,      Tf,       &
                                        strocnxT, strocnyT, &
                                        Tbot,     fbot,     &
                                        fside,   Cdn_ocn    )

       use ice_domain_size, only: ncat,nfsd,nilyr,nslyr
       use ice_constants, only: c1
       use ice_therm_vertical, only: ustar_min, fbot_xfer_type
       use ice_state, only: nt_qice, nt_qsno, nt_fsd

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         ilo,ihi,jlo,jhi   , & ! beginning and end of physical domain
         ntrcr                 ! number of tracers

      real (kind=dbl_kind), intent(in) :: &
         dt                  ! time step

      real (kind=dbl_kind), dimension(nx_block,ny_block), intent(in) :: &
         lead_area,& ! area near ice
         aice    , & ! ice concentration
         frzmlt  , & ! freezing/melting potential (W/m^2)
         sst     , & ! sea surface temperature (C)
         Tf      , & ! freezing temperature (C)
         Cdn_ocn , & ! ocean-ice neutral drag coefficient
         strocnxT, & ! ice-ocean stress, x-direction
         strocnyT    ! ice-ocean stress, y-direction

      real (kind=dbl_kind), dimension(nx_block,ny_block,ncat), &
         intent(in) :: &
         aicen   , & ! ITD
         vicen   , & ! ice volume (m)
         vsnon       ! snow volume (m)

      real (kind=dbl_kind), dimension(nx_block,ny_block,ntrcr,ncat), &   ! needs to be ntrcr here to match call
         intent(in) :: &
         trcrn       ! tracer array

      real (kind=dbl_kind), dimension(nx_block,ny_block), &
         intent(out) :: &
         Tbot    , & ! ice bottom surface temperature (deg C)
         fbot    , & ! heat flux to ice bottom  (W/m^2)
         fside       ! lateral heat flux (W/m^2) !LR

      ! local variables
! LR made this a local variable
      real (kind=dbl_kind), dimension(nx_block,ny_block,ncat) :: &
        rside_itd       ! fraction of ice that melts laterally

      real (kind=dbl_kind), dimension(ncat) :: &
        delta_an       ! amount of ice that melts laterally
! LR

      integer (kind=int_kind) :: &
         i, j           , & ! horizontal indices
         n              , & ! thickness category index
         k              , & ! layer index
         ij             , & ! horizontal index, combines i and j loops
         imelt              ! number of cells with ice melting

      integer (kind=int_kind), dimension (nx_block*ny_block) :: &
         indxi, indxj     ! compressed indices for cells with ice melting

      real (kind=dbl_kind), dimension (:), allocatable :: &
         etot        ! total energy in column
! LR         fside       ! lateral heat flux (W/m^2)


      real (kind=dbl_kind) :: &
! LR
        smfloe_arealoss, & ! change in area due to floes melting out of the
                         ! smallest floe size category
        G_radial   , & ! lateral melt rate, equal to negative wlat
! LR
         deltaT    , & ! SST - Tbot >= 0
         ustar     , & ! skin friction velocity for fbot (m/s)
         xtmp          ! temporary variable

      ! Parameters for bottom melting

      real (kind=dbl_kind) :: &
         cpchr         ! -cp_ocn*rhow*exchange coefficient

      ! Parameters for lateral melting

      real (kind=dbl_kind), parameter :: &
! LR     floediam = 300.0_dbl_kind, & ! effective floe diameter (m)
! LR     floeshape = 0.66_dbl_kind , & ! constant from Steele (unitless)
         m1 = 1.6e-6_dbl_kind     , & ! constant from Maykut & Perovich
                                      ! (m/s/deg^(-m2))
         m2 = 1.36_dbl_kind           ! constant from Maykut & Perovich
                                      ! (unitless)

      real (kind=dbl_kind), dimension(nx_block,ny_block) :: &
         wlat        ! lateral melt rate (m/s)


      ! LR fbot and fside will be zero if aice<puny or frzmlt < 0


      Tbot = Tf
      fbot = c0
      fside = c0 ! LR
      rside_itd = c0

      !-----------------------------------------------------------------
      ! Identify grid cells where ice can melt.
      !-----------------------------------------------------------------

      imelt = 0
      do j = jlo, jhi
      do i = ilo, ihi
         if (aice(i,j) > puny .and. frzmlt(i,j) < c0) then ! ice can melt
            imelt = imelt + 1
            indxi(imelt) = i
            indxj(imelt) = j
         endif
      enddo                     ! i
      enddo                     ! j

      allocate(etot (imelt))
! LR      allocate(fside(imelt))

      do ij = 1, imelt  ! cells where ice can melt
         i = indxi(ij)
         j = indxj(ij)

! LR         fside(ij) = c0

      !-----------------------------------------------------------------
      ! Use boundary layer theory for fbot.
      ! See Maykut and McPhee (1995): JGR, 100, 24,691-24,703.
      !-----------------------------------------------------------------

         deltaT = max((sst(i,j)-Tbot(i,j)),c0)

         ! strocnx has units N/m^2 so strocnx/rho has units m^2/s^2
         ustar = sqrt (sqrt(strocnxT(i,j)**2+strocnyT(i,j)**2)/rhow)
         ustar = max (ustar,ustar_min)

         if (trim(fbot_xfer_type) == 'Cdn_ocn') then
            ! Note: Cdn_ocn has already been used for calculating ustar 
            ! (formdrag only) --- David Schroeder (CPOM)
            cpchr = -cp_ocn*rhow*Cdn_ocn(i,j)
         else ! fbot_xfer_type == 'constant'
            ! 0.006 = unitless param for basal heat flx ala McPhee and Maykut
            cpchr = -cp_ocn*rhow*0.006_dbl_kind
         endif

         fbot(i,j) = cpchr * deltaT * ustar ! < 0

         fbot(i,j) = max (fbot(i,j), frzmlt(i,j)) ! frzmlt < fbot < 0

!!! uncomment to use all frzmlt for standalone runs
!        fbot(i,j) = min (c0, frzmlt(i,j))  !liuxy: if uncommentted, ic is larger in the margin of ice cover in Summer

      !-----------------------------------------------------------------
      ! Compute rside.  See these references:
      !    Maykut and Perovich (1987): JGR, 92, 7032-7044
      !    Hovart and Tziperman (2015): TC ?
      !-----------------------------------------------------------------

         wlat(i,j) = m1 * deltaT**m2 ! Maykut & Perovich 
         ! choose that this equals -Gr in Hovart
         G_radial = - wlat(i,j)

         ! LR - Compute rside but only to compute the lateral heat flux, fside
         ! The ITD and mFSTD are only actually evolved in ice_fsd_thermo,
         ! lateral_melt_fsdtherm subroutine
         ! We need to compute fside and not just wlat, because fside may
         ! be reduced so that fside + bottom < frzmlt
         rside_itd(i,j,:) = 0
         do n = 1, ncat
                       
                smfloe_arealoss = - trcrn(i,j,nt_fsd+1-1,n) / floe_binwidth(1) * &
                                dt * G_radial * aicen(i,j,n)
                
                delta_an(n)=c0
                do k=1,nfsd
                        ! delta_an is negative
                        delta_an(n) = delta_an(n) + ((c2/floe_rad_c(k))*aicen(i,j,n)* &
                                               trcrn(i,j,nt_fsd+k-1,n)*G_radial*dt) 
                end do                                                 
          
                ! add negative area loss from fsd
                delta_an(n) = delta_an(n) - smfloe_arealoss

                if(delta_an(n).gt.c0) stop 'delta_an gt0'

                ! to give same definition as in orginal code
                if (aicen(i,j,n).gt.c0) rside_itd(i,j,n)=-delta_an(n)/aicen(i,j,n) 
                ! otherwise rside_itd remains zero

                if (rside_itd(i,j,n).lt.c0) stop 'rside lt0'

         enddo ! n



      enddo                     ! ij


      !-----------------------------------------------------------------
      ! Compute heat flux associated with this value of rside.
      !-----------------------------------------------------------------

      do n = 1, ncat

         do ij = 1, imelt
            etot(ij) = c0
         enddo

         ! melting energy/unit area in each column, etot < 0

         do k = 1, nslyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do ij = 1, imelt
               i = indxi(ij)
               j = indxj(ij)
               etot(ij) = etot(ij) + trcrn(i,j,nt_qsno+k-1,n) &
                                   * vsnon(i,j,n)/real(nslyr,kind=dbl_kind)
            enddo               ! ij
         enddo

         do k = 1, nilyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do ij = 1, imelt
               i = indxi(ij)
               j = indxj(ij)
               etot(ij) = etot(ij) + trcrn(i,j,nt_qice+k-1,n) &
                                   * vicen(i,j,n)/real(nilyr,kind=dbl_kind)
            enddo               ! ij
         enddo                  ! nilyr

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, imelt
            i = indxi(ij)
            j = indxj(ij)
            ! lateral heat flux
	    ! CMB note rside is unique for each itd category 
            fside(i,j) = fside(i,j) + rside_itd(i,j,n)*etot(ij)/dt ! fside < 0
         enddo                  ! ij

      enddo                     ! n

      !-----------------------------------------------------------------
      ! Limit bottom and lateral heat fluxes if necessary.
      !-----------------------------------------------------------------

      do ij = 1, imelt
         i = indxi(ij)
         j = indxj(ij)

         if (rdc_frzmlt) then
                if ((aice(i,j) + lead_area(i,j)).gt.(c1+puny)) & 
                  stop 'c+oo gt 1'

                xtmp = MIN(c1,(aice(i,j) + lead_area(i,j))) * &
                       frzmlt(i,j)/(fbot(i,j) + fside(i,j) + puny) 
         else
                xtmp = frzmlt(i,j)/(fbot(i,j) + fside(i,j) + puny) 
         end if

         xtmp = min(xtmp, c1)
         fbot (i,j) = fbot (i,j) * xtmp
! LR
	 fside (i,j) = fside (i,j) * xtmp  
! LR 
      enddo                     ! ij
        
      deallocate(etot)
! LR      deallocate(fside)

      end subroutine frzmlt_bottom_lateral_fsd

!=======================================================================
!
! Normalize the floe size distribution so it sums to one in cells with ice.
! The FSD is zero is cells with no ice
!
! Includes some sanity checks for negative numbers
!
      subroutine renorm_mfstd(nx_block,ny_block,ncat,nfsd,aicen,trcrn)

        use ice_constants, only: puny, c0, c1
        use ice_state, only: nt_fsd
        use ice_domain_size, only: max_ntrcr

        integer(kind=int_kind), intent(in) :: &
             nx_block , &
             ny_block , &
             ncat, nfsd

        real (kind=dbl_kind), dimension (nx_block,ny_block,max_ntrcr,ncat), &  ! needs to be max_ntrcr here
             intent(inout) :: &
             trcrn     ! tracer array

        real(kind=dbl_kind), dimension(nx_block,ny_block,  ncat), &
             intent(in) :: aicen

        ! local variables

        integer (kind=int_kind) :: &
             n, &          ! ice thickness category index
             i,j


        do n=1,ncat

                do j = 1,ny_block
                do i = 1,nx_block

                if (aicen(i,j,n).lt.c0) stop 'negative aice'
                if (ANY(trcrn(i,j,nt_fsd:nt_fsd+nfsd-1,n).lt.c0-1000*puny)) then
                        print *, 'mFSTD ',trcrn(i,j,nt_fsd:nt_fsd+nfsd-1,n)
                        print *, 'aicen ',aicen(i,j,n),i,j,n
                        print *, 'this is bad but going forward anyway'
!			stop 'negative mFSTD'
                end if
                if (ANY(trcrn(i,j,nt_fsd:nt_fsd+nfsd-1,n).lt.c0)) then
                        print *, 'slightly negative mFSTD, set to zero'
                        print *, 'mFSTD ',trcrn(i,j,nt_fsd:nt_fsd+nfsd-1,n)
                        print *, 'aicen ',aicen(i,j,n),i,j,n
                        WHERE(trcrn(i,j,nt_fsd:nt_fsd+nfsd-1,n).lt.c0) &
                                trcrn(i,j,nt_fsd:nt_fsd+nfsd-1,n) = c0
                end if
         

                ! mFSTD is zero when there is no ice
                if (aicen(i,j,n).le.puny) then 
                        trcrn(i,j,nt_fsd:nt_fsd+nfsd-1,n) = c0
                else
                        if (SUM(trcrn(i,j,nt_fsd:nt_fsd+nfsd-1,n)).lt.puny) then
                                print *, aicen(i,j,n)
                                print *, trcrn(i,j,nt_fsd:nt_fsd+nfsd-1,n)
                                print *, SUM(trcrn(i,j,nt_fsd:nt_fsd+nfsd-1,n))
                                stop 'mFSTD zero for non-zero aicen'
                        end if
                        if (ABS(SUM(trcrn(i,j,nt_fsd:nt_fsd+nfsd-1,n))-c1).gt.c0) then
                                !print *, 'renorm necessary (called renorm_mfstd)'
                                trcrn(i,j,nt_fsd:nt_fsd+nfsd-1,n) = trcrn(i,j,nt_fsd:nt_fsd+nfsd-1,n) / SUM(trcrn(i,j,nt_fsd:nt_fsd+nfsd-1,n))
                        end if
                end if

                enddo !i
                enddo !j
        enddo !n

      end subroutine renorm_mfstd

!=======================================================================
!=======================================================================
!
! Normalize the floe size distribution so it sums to one in cells with ice.
! The FSD is zero is cells with no ice
!
! Includes some sanity checks for negative numbers
!
      subroutine check_mfstd(nx_block,ny_block,ncat,nfsd,aicen,trcrn,mssg)

        use ice_constants, only: puny, c0, c1
        use ice_state, only: nt_fsd
        use ice_domain_size, only: max_ntrcr

        integer(kind=int_kind), intent(in) :: &
             nx_block , &
             ny_block , &
             ncat, nfsd

       character (len=*), intent(in) :: mssg

        real (kind=dbl_kind), dimension (nx_block,ny_block,max_ntrcr,ncat), &  ! needs to be max_ntrcr here
             intent(inout) :: &
             trcrn     ! tracer array

        real(kind=dbl_kind), dimension(nx_block,ny_block,  ncat), &
             intent(in) :: aicen

        ! local variables

        integer (kind=int_kind) :: &
             n, &          ! ice thickness category index
             i,j


        do n=1,ncat

                do j = 1,ny_block
                do i = 1,nx_block

                if (aicen(i,j,n).lt.c0) print*, 'check_mfstd: negative aicen',aicen(i,j,n),i,j,n
                if (ANY(trcrn(i,j,nt_fsd:nt_fsd+nfsd-1,n).lt.c0-1000*puny)) then
                     print*,'check mfstd: ',mssg
                        print *, 'mFSTD ',trcrn(i,j,nt_fsd:nt_fsd+nfsd-1,n)
                        print *, 'aicen ',aicen(i,j,n),i,j,n
!			stop 'negative mFSTD'
!                else if (ANY(trcrn(i,j,nt_fsd:nt_fsd+nfsd-1,n).lt.c0)) then
!                     print*,mssg
!                        print *, 'slightly negative mFSTD, NOT set to zero'
!                        print *, 'mFSTD ',trcrn(i,j,nt_fsd:nt_fsd+nfsd-1,n)
!                        print *, 'aicen ',aicen(i,j,n)
!                        WHERE(trcrn(i,j,nt_fsd:nt_fsd+nfsd-1,n).lt.c0) &
!                                trcrn(i,j,nt_fsd:nt_fsd+nfsd-1,n) = c0
                end if
         

                enddo !i
                enddo !j
        enddo !n

      end subroutine check_mfstd

!=======================================================================
! 
! Given a wave spectrum, calculate size of new floes based on tensile failire
! See Shen & Ackley (2004), Roach, Smith & Dean (2018) for further details
! Author: Lettie Roach (NIWA) 2018

      subroutine wave_dep_growth (local_wave_spec, &
                                       new_size )

       use ice_flux, only: dfreq, freq
 
      real (kind=dbl_kind), dimension(nfreq), intent(in) :: &
           local_wave_spec ! e(f), dimension set in ice_forcing or ice_wavebreaking

      integer (kind=int_kind), intent(out) :: &
           new_size ! index of floe size category in which new floes will growh

      ! local variables
      real (kind=dbl_kind), parameter :: &
          tensile_param = 0.167_dbl_kind

      real (kind=dbl_kind)  :: &
          mo,   &   ! zeroth moment of the spectrum (m)
          h_s,  &   ! significant wave height (m)
          w_a,  &   ! wave amplitude (m)
          f_p,  &   ! peak frequency (s^-1)
          w_l,  &   ! wavelength from peak freqency (m) 
          d_max, &  ! d_max from tensile failure mode
          r_max     ! radius

      integer (kind=int_kind) :: k


      ! zeroth moment
      mo = SUM(local_wave_spec*dfreq)

      ! sig wave height and amplitude
      h_s = c4*SQRT(mo)
      w_a = h_s/c2

      ! peak frequency
      f_p = freq(MAXLOC(local_wave_spec, DIM=1))

      ! wavelength from peak freq
      w_l = gravit / (c2*pi*f_p**c2)

      ! tensile failure
      if (w_a.gt.puny) then
          d_max = SQRT(c2*tensile_param*w_l**c2/(pi**c3*w_a*gravit*rhoi))
          r_max = d_max/c2
      else
          r_max = bignum
      end if

      new_size = nfsd
      do k = 1, nfsd - 1
          if (r_max.le.floe_rad_h(k)) then
              new_size = k
              EXIT
          end if
      end do


      end subroutine wave_dep_growth


!=======================================================================

      end module ice_fsd

!=======================================================================

