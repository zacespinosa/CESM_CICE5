

! Thermodynamic growth and melt for a floe size distribution
! Closely follows the method in Horvat & Tziperman (2015)
!
! Lettie Roach, 2016
! 
!
      module ice_fsd_thermo

      use ice_kinds_mod
      use ice_constants
      use ice_domain_size, only: ncat, nfsd
      use ice_state, only: nt_fsd
      use ice_fsd, only: write_diag_diff

      implicit none

      private
      public :: partition_area, & 
                lateral_melt_fsdtherm, &
                floe_merge_thermo, add_new_ice_lat

      integer (kind=int_kind), public :: &
        new_ice_fs     ! how does new ice grow?
                
!=======================================================================

      contains

!=======================================================================

!=======================================================================
! 
!  Given the joint ice thickness and floe size distribution, calculate
!  the lead region and the total lateral surface area following Horvat
!  and Tziperman (2015)
!
! author: Lettie Roach, NIWA/VUW

      subroutine partition_area (nx_block, ny_block, &
                                        ilo, ihi, jlo, jhi, &
                                        ntrcr,              &
                                        aice,               &
                                        aicen,    vicen,    &
                                        trcrn,    lead_area,&
                                        latsurf_area )

      use ice_fsd, only: floe_rad_l, floe_rad_h, floe_rad_c, &
                          floe_binwidth, floe_area_c, floeshape
      use ice_itd, only: hin_max_init

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         ilo,ihi,jlo,jhi,    & ! beginning and end of physical domain
         ntrcr                 ! number of tracers
       
      real (kind=dbl_kind), dimension(nx_block,ny_block), intent(in) :: &
         aice        ! ice concentration
        
      real (kind=dbl_kind), dimension(nx_block,ny_block,ncat), &
         intent(in) :: &
         aicen, &      ! fractional area of ice 
         vicen         ! volume per unit area of ice

      real (kind=dbl_kind), dimension(nx_block,ny_block,ntrcr,ncat), & 
         intent(in) :: &
         trcrn       ! tracer array

      real (kind=dbl_kind), dimension(nx_block,ny_block), &
         intent(out) :: &
          lead_area, &  ! the fractional area of the lead region
          latsurf_area  ! the area covered by lateral surface of floes
      
      ! local variables
      
      integer (kind=int_kind) :: &
         i, j           , & ! horizontal indices
         n              , & ! thickness category index
         k                  ! floe size index

      real (kind=dbl_kind) :: &
        width_leadreg, &   ! width of lead region
        thickness          ! actual thickness of ice in thickness cat

      ! -----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------
       lead_area=c0
       latsurf_area=c0

       ! Set the width of the lead region to be the smallest
       ! floe size category, as per Horvat & Tziperman (2015)
       width_leadreg=floe_rad_c(1)
      
      !-----------------------------------------------------------------
      ! Loop over all gridcells. Only calculate these areas where there
      ! is some ice
      !-----------------------------------------------------------------

       do j = jlo, jhi
       do i = ilo, ihi
      
          lead_area(i,j) = c0
          latsurf_area(i,j) = c0
   
          if (aice(i,j).gt.puny) then

                ! lead area = sum of areas of annuli around all floes
                do n=1,ncat       
                        do k=1,nfsd
                                lead_area(i,j) = lead_area(i,j) + &
                                                 aicen(i,j,n) * trcrn(i,j,nt_fsd+k-1,n) * &
                                                 ( c2*width_leadreg/floe_rad_c(k) + &
                                                 width_leadreg**c2/floe_rad_c(k)**2)
                        enddo !k
                enddo !n
       
                ! cannot be greater than the open water fraction
                lead_area(i,j)=MIN(lead_area(i,j),(c1-aice(i,j)))
      
                ! sanity checks
                if (lead_area(i,j).gt.c1) stop 'lead_area not frac!'
                if (lead_area(i,j).ne.lead_area(i,j)) stop 'lead_a NaN'
                if (lead_area(i,j).lt.c0) then
                        if (lead_area(i,j).lt.(c0-puny)) then
                                stop 'lead_area lt0 in partition_area'
                        else
                                lead_area(i,j)=MAX(lead_area(i,j),c0)
                        end if
                end if

                ! Total fractional lateral surface area in each grid (per unit ocean area)
                do n=1,ncat
                    thickness = c0
                    if (aicen(i,j,n).gt.c0) thickness = vicen(i,j,n)/aicen(i,j,n)

                        do k=1,nfsd
                             latsurf_area(i,j) = latsurf_area(i,j) + &
                                                    trcrn(i,j,nt_fsd+k-1,n) * aicen(i,j,n) * & ! FSD
                                                    c2 * thickness/floe_rad_c(k)
                        end do
                end do 

                ! check
                if (latsurf_area(i,j).lt.c0) stop &
                          'negative latsurf_ area'
                if (lead_area(i,j).ne.lead_area(i,j)) stop &
                          'latsurf_ area NaN'
         end if ! aice
       end do !i
       end do !j

              end subroutine partition_area

!=======================================================================

!=======================================================================
! Melts the ice thickness distribution and the floe size distribution, 
! given fside, which was computed in frzmlt_bottom_lateral_fsdtherm
! 
! author: C. M. Bitz, UW
! 2003:   Modified by William H. Lipscomb and Elizabeth C. Hunke, LANL
! 2016:   Modified by C. M. Bitz, UW, to add floe effects
!         Further modified by L. Roach, NIWA


      subroutine lateral_melt_fsdtherm & 
                              (nx_block,   ny_block,   &
                               ilo, ihi,   jlo, jhi,   &
                               dt,         fpond,      &
                               fresh,      fsalt,      &
                               fhocn,                  &
                               faero_ocn,              &
! LR for CESM
                               fiso_ocn,               &
! LR for CESM
                               meltl,      aice,       &
                               aicen,      vicen,      &
                               vsnon,      trcrn,      &
                               fside,      frzmlt,     &
                               sss,     dSin0_frazil,  &
                               phi_init,  salinz,      &
                               G_radial )

      use ice_state, only: nt_qice, nt_qsno, nt_aero, tr_aero, &
                           tr_pond_topo, nt_apnd, nt_hpnd, &
                           nt_iso, tr_iso ! LR for CESM
      use ice_domain_size,only: nilyr, nslyr, n_aero, &
                                max_aero, max_ntrcr, &
                                max_iso,  n_iso ! LR for CESM
      use ice_fsd, only: floe_rad_c, floe_binwidth
      use ice_therm_mushy, only: liquidus_temperature_mush, enthalpy_mush
      use ice_therm_shared, only: ktherm

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         ilo,ihi,jlo,jhi       ! beginning and end of physical domain

      real (kind=dbl_kind), intent(in) :: &
         dt        ! time step (s)

       real (kind=dbl_kind), intent(in) :: &
         phi_init     , & ! initial frazil liquid fraction
         dSin0_frazil     ! initial frazil bulk salinity reduction from sss

      real (kind=dbl_kind), dimension(nx_block,ny_block,nilyr+1), intent(in) :: &
         salinz           ! initial salinity profile 

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat), &
         intent(inout) :: & !inout
         aicen   , & ! concentration of ice
         vicen   , & ! volume per unit area of ice          (m)
         vsnon       ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_ntrcr,ncat), &  ! needs to be max_ntrcr here
         intent(inout) :: & !inout
         trcrn     ! tracer array

      real (kind=dbl_kind), dimension(nx_block,ny_block), intent(in) :: &
         fside, &     ! flux that goes to lateral melt (m/s)
         frzmlt, &    ! freezing/melting potential (Wm/m^2)
         sss, &       ! sea surface salinity (ppt)
         aice         ! ice concentration

      real (kind=dbl_kind), dimension(nx_block,ny_block), &
         intent(inout) :: & !inout
         G_radial  , & ! rate of lateral melt (m/s)
         fpond     , & ! fresh water flux to ponds (kg/m^2/s)
         fresh     , & ! fresh water flux to ocean (kg/m^2/s)
         fsalt     , & ! salt flux to ocean (kg/m^2/s)
         fhocn     , & ! net heat flux to ocean (W/m^2)
         meltl         ! lateral ice melt         (m/step-->cm/day)
  
      real (kind=dbl_kind), dimension(nx_block,ny_block,max_aero), &
         intent(inout) :: & !intout
         faero_ocn     ! aerosol flux to ocean (kg/m^2/s)

! LR for CESM
      real (kind=dbl_kind), dimension(nx_block,ny_block,max_iso), &
         intent(inout) :: &
         fiso_ocn     ! isotope flux to ocean (kg/m^2/s)
! LR for CESM

      ! local variables

      integer (kind=int_kind) :: &
         i, j        , & ! horizontal indices
         n           , & ! thickness category index
         k, kk       , & ! layer index
         ij          , & ! horizontal index, combines i and j loops
         icells          ! number of cells with aice > puny

      integer (kind=int_kind), dimension(nx_block*ny_block) :: &
         indxi, indxj    ! compressed indices for cells with aice > puny

      real (kind=dbl_kind) :: &
         totfrac , & ! used to renormalize floe size dist
         dfhocn  , & ! change in fhocn
         dfpond  , & ! change in fpond
         dfresh  , & ! change in fresh
         dfsalt  , & ! change in fsalt
         Ti          ! frazil temperature

      real (kind=dbl_kind), dimension (nx_block,ny_block) :: &
        vicen_init   ! volume per unit area of ice          (m)
      
      real (kind=dbl_kind), dimension (nfsd) :: &
         areal_mfstd_final, & ! modified areal FSTD (tilda) 
         fin_diff             ! finite difference for G_r * areal mFSTD tilda

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat) :: &
         aicen_init, &
         rside_itd         ! delta_an/aicen

      real (kind=dbl_kind), dimension (ncat) :: &
         delta_an        ! change in the ITD

     real (kind=dbl_kind), dimension (nx_block,ny_block,nfsd,ncat) :: & 
         areal_mfstd_init

      real (kind=dbl_kind), dimension (nx_block,ny_block) :: &
         qi0          , & ! frazil ice enthalpy
         Si0              ! frazil ice bulk salinity

      real (kind=dbl_kind) :: cat1_arealoss

      real (kind=dbl_kind), dimension(nfsd+1) :: &
        f_flx

     ! initialization
     rside_itd = c0
     aicen_init = aicen           
     areal_mfstd_init = trcrn(:,:,nt_fsd:nt_fsd+nfsd-1,:)
     !-----------------------------------------------------------------
     ! Identify grid cells with lateral melt
     !-----------------------------------------------------------------

     icells = 0
     do j = jlo, jhi
     do i = ilo, ihi
         if (fside(i,j) < c0) then
             icells = icells + 1
             indxi(icells) = i
             indxj(icells) = j
         endif
     enddo                  ! i
     enddo                  ! j
     !-----------------------------------------------------------------                                              
     ! In these cells, ice can melt:                                                      
     !-----------------------------------------------------------------                                              
                                                                               
 
     !DIR$ CONCURRENT !Cray                                                                                                
     !cdir nodep      !NEC
     !ocl novrec      !Fujitsu
     do ij = 1, icells
              i = indxi(ij)
              j = indxj(ij)
         
              !-----------------------------------------------------------------
              ! Compute average enthalpy of ice. (taken from add_new_ice)
              ! Not sure if this is transferrable to existing ice?!?!
              ! Sprofile is the salinity profile used when adding new ice to
              ! all categories, for ktherm/=2, and to category 1 for all ktherm.
              !
              ! NOTE:  POP assumes new ice is fresh!
              !-----------------------------------------------------------------

              if (ktherm == 2) then  ! mushy
        !DIR$ CONCURRENT !Cray
        !cdir nodep      !NEC
        !ocl novrec      !Fujitsu
                    if (sss(i,j) > c2 * dSin0_frazil) then
                       Si0(i,j) = sss(i,j) - dSin0_frazil
                    else
                       Si0(i,j) = sss(i,j)**2 / (c4*dSin0_frazil)
                    endif
                    Ti = min(liquidus_temperature_mush(Si0(i,j)/phi_init), -p1)
                    qi0(i,j) = enthalpy_mush(Ti, Si0(i,j))

              else

        !DIR$ CONCURRENT !Cray
        !cdir nodep      !NEC
        !ocl novrec      !Fujitsu
                   qi0(i,j) = -rhoi*Lfresh
                             
              endif    ! ktherm

              !-----------------------------------------------------------------
              ! G_r is the lateral melt rate                                               
              !-----------------------------------------------------------------
 
              G_radial(i,j) = -fside(i,j)/qi0(i,j) !negative
              if (G_radial(i,j).gt.c0) stop 'Gr pos for melt'                
              !-----------------------------------------------------------------
              ! Given G_r, compute change to the ITD
              !-----------------------------------------------------------------
              if (G_radial(i,j).lt.(c0-puny)) then
                do n = 1, ncat

                    if (aicen(i,j,n).gt.puny) then
                        if (ABS(SUM(areal_mfstd_init(i,j,:,n))-c1).gt.1.0e-9_dbl_kind) then
                                print *, ABS(SUM(areal_mfstd_init(i,j,:,n))-c1)
                                print *, SUM(areal_mfstd_init(i,j,:,n))
                                print *, SUM(trcrn(i,j,nt_fsd:nt_fsd+nfsd-1,n))
                                print *, &
                        'WARNING init mFSTD not normed, lm'
                        end if
                        areal_mfstd_init(i,j,:,n) = areal_mfstd_init(i,j,:,n)/SUM(areal_mfstd_init(i,j,:,n)) 
                    end if

                        cat1_arealoss = - trcrn(i,j,nt_fsd+1-1,n) / floe_binwidth(1) * &
                                        dt * G_radial(i,j)*aicen(i,j,n)

                        delta_an(n)=c0
                        do k=1,nfsd
                                ! delta_an is negative
                                delta_an(n) = delta_an(n) + ((c2/floe_rad_c(k))*aicen(i,j,n)* &
                                                       trcrn(i,j,nt_fsd+k-1,n)*G_radial(i,j)*dt) 
                        end do                                                 
                  
                        ! add negative area loss from fsd
                        delta_an(n) = delta_an(n) - cat1_arealoss

                        if(delta_an(n).gt.c0) stop 'delta_an gt0'
    
                        ! to give same definition as in orginal code
                        if (aicen(i,j,n).gt.c0) rside_itd(i,j,n)=-delta_an(n)/aicen(i,j,n) 
                        ! otherwise rside_itd remains zero

                        if (rside_itd(i,j,n).lt.c0) stop 'rside lt0'
       
                enddo ! n
              end if ! G_r

     enddo !ij

     !-----------------------------------------------------------------                                            
     ! Increment fluxes and melt the ice using rside as per 
     ! existing routine
     !-----------------------------------------------------------------                                            
     do n = 1, ncat

        !DIR$ CONCURRENT !Cray
        !cdir nodep      !NEC
        !ocl novrec      !Fujitsu
                 do ij = 1, icells
                    i = indxi(ij)
                    j = indxj(ij)

                    ! fluxes to coupler
                    ! dfresh > 0, dfsalt > 0, dfpond > 0

                    dfresh = (rhos*vsnon(i,j,n) + rhoi*vicen(i,j,n)) &
                           * rside_itd(i,j,n) / dt
                    dfsalt = rhoi*vicen(i,j,n)*ice_ref_salinity*p001 &
                           * rside_itd(i,j,n) / dt
                    fresh(i,j)      = fresh(i,j)      + dfresh
                    fsalt(i,j)      = fsalt(i,j)      + dfsalt

                    if (tr_pond_topo) then
                       dfpond = aicen(i,j,n) &
                              * trcrn(i,j,nt_apnd,n) & 
                              * trcrn(i,j,nt_hpnd,n) &
                              * rside_itd(i,j,n)
                       fpond(i,j) = fpond(i,j) - dfpond
                    endif

                    ! history diagnostics
                    meltl(i,j) = meltl(i,j) + vicen(i,j,n)*rside_itd(i,j,n)

                    ! state variables
                    vicen_init(i,j) = vicen(i,j,n)
                    aicen(i,j,n) = aicen(i,j,n) * (c1 - rside_itd(i,j,n))
                    vicen(i,j,n) = vicen(i,j,n) * (c1 - rside_itd(i,j,n))
                    vsnon(i,j,n) = vsnon(i,j,n) * (c1 - rside_itd(i,j,n))

 
                    !-----------------------------------------------------------------  
                    ! Now compute the change to the mFSTD
                    !-----------------------------------------------------------------                                            
                    ! FSD
                    if (rside_itd(i,j,n).gt.puny) then
                    if (aicen(i,j,n).gt.puny) then

                        fin_diff(:) = c0
                        f_flx(:) = c0
                        do k=  2, nfsd
                                f_flx(k) =  G_radial(i,j) * &
                                            areal_mfstd_init(i,j,k,n)/ &
                                            floe_binwidth(k)
 
                        end do

                        do k = 1, nfsd
                                fin_diff(k) = f_flx(k+1) - f_flx(k)
                        end do

                        if (ABS(SUM(fin_diff(:))).gt.puny) stop &
                                 'sum fnk diff not zero in lm'

                        do k = 1,nfsd
                                areal_mfstd_final(k) = &
                                areal_mfstd_init(i,j,k,n) +   &
                                dt * (  - fin_diff(k) + &
                                c2 * G_radial(i,j) * areal_mfstd_init(i,j,k,n) * &
                                (c1/floe_rad_c(k) - & 
                                SUM(areal_mfstd_init(i,j,:,n)/floe_rad_c(:))) )
                        end do
                       
 
                        if (ABS(SUM(areal_mfstd_final)-c1).gt.puny) then
                                print *, SUM(fin_diff)
                                print *, SUM(areal_mfstd_final)-c1
                                stop &
                                'mFSTD not normed, lm' 
                        end if

                        ! this fixes tiny (e-30) differences from 1
                        areal_mfstd_final = areal_mfstd_final/SUM(areal_mfstd_final)

                        if (ANY(areal_mfstd_final.lt.c0)) stop &
                                'neg mFSTD, lm'

                        trcrn(i,j,nt_fsd:nt_fsd+nfsd-1,n) = areal_mfstd_final
                   else
                        trcrn(i,j,nt_fsd:nt_fsd+nfsd-1,n) = c0
                   end if !aicen>0
                   end if ! rside>0, otherwise do nothing
            
                   ! remove?
                   if (ANY(trcrn(i,j,nt_fsd:nt_fsd+nfsd-1,n).gt.c1+puny)) stop  &
                        'mFSTD > 1 in lat melt'

                    if ((aicen(i,j,n).ne.aicen(i,j,n)).or.(vicen(i,j,n).ne.vicen(i,j,n))) stop &
                        'aicen or vicen NaN after update in latmelt'

                 enddo                  ! ij

                 do k = 1, nilyr
        !DIR$ CONCURRENT !Cray
        !cdir nodep      !NEC
        !ocl novrec      !Fujitsu
                    do ij = 1, icells
                       i = indxi(ij)
                       j = indxj(ij)

                       ! enthalpy tracers do not change (e/v constant)
                       ! heat flux to coupler for ice melt (dfhocn < 0)
                       dfhocn = trcrn(i,j,nt_qice+k-1,n)*rside_itd(i,j,n) / dt &
                              * vicen(i,j,n)/real(nilyr,kind=dbl_kind)
                       fhocn(i,j)      = fhocn(i,j)      + dfhocn
                    enddo               ! ij
                 enddo                  ! nilyr

                 do k = 1, nslyr
        !DIR$ CONCURRENT !Cray
        !cdir nodep      !NEC
        !ocl novrec      !Fujitsu
                    do ij = 1, icells
                       i = indxi(ij)
                       j = indxj(ij)

                       ! heat flux to coupler for snow melt (dfhocn < 0)

                       dfhocn = trcrn(i,j,nt_qsno+k-1,n)*rside_itd(i,j,n) / dt &
                              * vsnon(i,j,n)/real(nslyr,kind=dbl_kind)
                       fhocn(i,j)      = fhocn(i,j)      + dfhocn
                    enddo               ! ij
                 enddo                  ! nslyr

                 if (tr_aero) then
                    do k = 1, n_aero
        !DIR$ CONCURRENT !Cray
        !cdir nodep      !NEC
        !ocl novrec      !Fujitsu
                       do ij = 1, icells
                          i = indxi(ij)
                          j = indxj(ij)
                          faero_ocn(i,j,k) = faero_ocn(i,j,k) + (vsnon(i,j,n) &
                                           *(trcrn(i,j,nt_aero  +4*(k-1),n)   &
                                           + trcrn(i,j,nt_aero+1+4*(k-1),n))  &
                                                              +  vicen(i,j,n) &
                                           *(trcrn(i,j,nt_aero+2+4*(k-1),n)   &
                                           + trcrn(i,j,nt_aero+3+4*(k-1),n))) &
                                           * rside_itd(i,j,n) / dt
                       enddo
                    enddo
                 endif

		 if (tr_iso) then
		    do k = 1, n_iso
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
	               do ij = 1, icells
			  i = indxi(ij)
			  j = indxj(ij)
			  fiso_ocn(i,j,k) = fiso_ocn(i,j,k) + (vsnon(i,j,n) &
					   *(trcrn(i,j,nt_iso  +4*(k-1),n)   &
					   + trcrn(i,j,nt_iso+1+4*(k-1),n))  &
							      +  vicen(i,j,n) &
					   *(trcrn(i,j,nt_iso+2+4*(k-1),n)   &
					   + trcrn(i,j,nt_iso+3+4*(k-1),n))) &
					   * rside_itd(i,j,n) / dt 
		       enddo
		    enddo
		 endif

     enddo  ! n


                    end subroutine lateral_melt_fsdtherm

!=======================================================================
!
! Given the volume of new ice grown in open water, compute its area
! and thickness and add it to the appropriate category or categories.
!
! NOTE: Usually all the new ice is added to category 1.  An exception is
!       made if there is no open water or if the new ice is too thick
!       for category 1, in which case ice is distributed evenly over the
!       entire cell.  Subroutine rebin should be called in case the ice
!       thickness lies outside category bounds after new ice formation.
!
! When ice must be added to categories above category 1, the mushy
! formulation (ktherm=2) adds it only to the bottom of the ice.  When
! added to only category 1, all formulations combine the new ice and
! existing ice tracers as bulk quantities.
!
! authors William H. Lipscomb, LANL
!         Elizabeth C. Hunke, LANL
!         Adrian Turner, LANL
!
              subroutine add_new_ice_lat (nx_block,  ny_block,   &
                                      ntrcr,     icells,     &
                                      indxi,     indxj,      &
                                      dt,                    &
                                      lead_area, latsurf_area, & ! LR
                                      aicen,     trcrn,      &
                                      vicen,                 &
                                      aice0,     aice,       &
                                      frzmlt,    frazil,     &
                                      vlateral,              &
                                      frz_onset, yday,       &
                                      update_ocn_f,          &
                                      fresh,     fsalt,      &
                                      Tf,        sss,        &
                                      salinz,    phi_init,   &
                                      dSin0_frazil,          &
! LR for CESM
                                      frazil_diag,           &
				      fiso_ocn,              &
				      HDO_ocn,             &
				      H2_16O_ocn,          &
				      H2_18O_ocn,          & 
! LR for CESM
                                      nbtrcr,    flux_bio,   &
                                      ocean_bio, &
                                      l_stop,                &
                                      istop,     jstop      , &
                                      d_an_latg,        d_an_addnew,    &
                                      d_afsd_latg,      d_afsd_addnew,  &
                                      d_amfstd_latg,    d_amfstd_addnew,&
                                      G_radial, tarea, wave_spectrum )
         
              use ice_domain_size, only: nilyr, n_aero, & 
! LR
                                   max_iso, n_iso ! for CESM 
              use ice_fsd, only: floe_rad_c, floe_binwidth, &
                                  floe_area_l, nfreq, wave_dep_growth
! LR     
              use ice_itd, only: hin_max, column_sum, &
                                 column_conservation_check 
              use ice_state, only: nt_Tsfc, nt_iage, nt_FY, nt_alvl, nt_vlvl, nt_aero, &
                                   nt_sice, nt_qice, nt_fsd, &
                                   nt_apnd, tr_pond_cesm, tr_pond_lvl, tr_pond_topo, &
                                   tr_iage, tr_FY, tr_lvl, tr_aero, tr_brine, tr_fsd, & ! CMB
! LR for CESM
                                   nt_iso, tr_iso
              use ice_isotope, only: isoice_alpha, frac
! LR for CESM
              use ice_therm_itd, only: l_conservation_check,update_vertical_tracers
              use ice_therm_mushy, only: liquidus_temperature_mush, enthalpy_mush
              use ice_therm_shared, only: ktherm, hfrazilmin
              use ice_zbgc, only: add_new_ice_bgc
              use ice_zbgc_shared, only: skl_bgc

              integer (kind=int_kind), intent(in) :: &
                 nx_block, ny_block, & ! block dimensions
                 ntrcr             , & ! number of tracers in use
                 icells                ! number of ice/ocean grid cells

              integer (kind=int_kind), dimension (nx_block*ny_block), &
                 intent(in) :: &
                 indxi,  indxj         ! compressed i/j indices

              real (kind=dbl_kind), intent(in) :: &
                 dt        ! time step (s)

              real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
                 tarea, &  ! grid cell area (m^2)
                 lead_area, & ! fractional area of ice in lead region
                 latsurf_area, & ! fractional area of ice on sides of floes
                 aice  , & ! total concentration of ice
                 frzmlt, & ! freezing/melting potential (W/m^2)
                 Tf    , & ! freezing temperature (C)
                 sss       ! sea surface salinity (ppt)

              real (kind=dbl_kind), dimension (nx_block,ny_block,ncat), &
                 intent(inout) :: &
                 aicen , & ! concentration of ice
                 vicen     ! volume per unit area of ice          (m)

              real (kind=dbl_kind), dimension (nx_block,ny_block,ntrcr,ncat), &
                 intent(inout) :: &
                 trcrn     ! ice tracers
                           ! 1: surface temperature

              real (kind=dbl_kind), dimension (nx_block,ny_block), &
                 intent(inout) :: &
                 G_radial  , & ! lateral melt rate (m/s)
                 aice0     , & ! concentration of open water
                 frazil    , & ! frazil ice growth        (m/step-->cm/day)
! LR for CESM
                 frazil_diag, & ! frazil ice growth diagnostic (m/step-->cm/day)
! LR for CESM
                 vlateral    , & ! lateral ice growth        (m/step-->cm/day)
                 fresh     , & ! fresh water flux to ocean (kg/m^2/s)
                 fsalt         ! salt flux to ocean (kg/m^2/s)

              real (kind=dbl_kind), dimension (nx_block,ny_block), &
                 intent(inout), optional :: &
                 frz_onset ! day of year that freezing begins (congel or frazil)

              real (kind=dbl_kind), intent(in), optional :: &
                 yday      ! day of year

              real (kind=dbl_kind), dimension(nx_block,ny_block,nilyr+1), intent(in) :: &
                 salinz     ! initial salinity profile

              real (kind=dbl_kind), intent(in) :: &
                 phi_init     , & ! initial frazil liquid fraction
                 dSin0_frazil     ! initial frazil bulk salinity reduction from sss
! LR for CESM
	      real (kind=dbl_kind), dimension (nx_block,ny_block,nt_iso), &
		 intent(inout), optional :: &
		 fiso_ocn  ! kg/m^2/s

	      real (kind=dbl_kind), dimension (nx_block,ny_block), &
		 intent(inout), optional :: &
		 HDO_ocn , & !
		 H2_16O_ocn, & !
		 H2_18O_ocn
! LR for CESM

              logical (kind=log_kind), intent(in) :: &
                 update_ocn_f ! if true, update fresh water and salt fluxes

              logical (kind=log_kind), intent(out) :: &
                 l_stop    ! if true, abort on return

              integer (kind=int_kind), intent(out) :: &
                 istop, jstop    ! indices of grid cell where model aborts

              real (kind=dbl_kind), dimension(nx_block,ny_block,ncat), intent(out) :: &
                d_an_latg, d_an_addnew

              real (kind=dbl_kind), dimension(nx_block,ny_block,nfsd,ncat), intent(out) :: &
                d_amfstd_latg, d_amfstd_addnew

              real (kind=dbl_kind), dimension(nx_block,ny_block,nfsd), intent(out) :: &
                d_afsd_latg, d_afsd_addnew

              real (kind=dbl_kind), dimension(nx_block,ny_block, nfreq), intent(in)  :: &
                wave_spectrum


             ! BGC
              integer (kind=int_kind), intent(in) :: &
                 nbtrcr          ! number of biology tracers

              real (kind=dbl_kind), dimension (nx_block,ny_block,nbtrcr), &
                 intent(inout) :: &
                flux_bio   ! tracer flux to ocean from biology (mmol/m^2/s)
                
              real (kind=dbl_kind), dimension (nx_block,ny_block,nbtrcr), &
                 intent(in) :: &
                 ocean_bio   ! ocean concentration of biological tracer

              ! local variables

              integer (kind=int_kind) :: &
! LR
                 new_size     , & ! index for floe size of new ice 
! LR
                 i, j         , & ! horizontal indices
                 n            , & ! ice category index
                 k            , & ! ice layer index
                 it               ! aerosol tracer index

              real (kind=dbl_kind), dimension (icells) :: &
                 ai0new       , & ! area of new ice added to cat 1
                 vi0new       , & ! volume of new ice added to cat 1
                 vi0new_lat   , & ! LR
                 hsurp            ! thickness of new ice added to each cat

              real (kind=dbl_kind), dimension (icells) :: &
                 vice1        , & ! starting volume of existing ice
                 vice_init, vice_final, & ! ice volume summed over categories
                 eice_init, eice_final    ! ice energy summed over categories

              real (kind=dbl_kind) :: &
                 vsurp_lat, & ! LR
                 fnew         , & ! heat flx to open water for new ice (W/m^2)
                 hi0new       , & ! thickness of new ice
                 hi0max       , & ! max allowed thickness of new ice
                 vsurp        , & ! volume of new ice added to each cat
                 asurp        , &
                 vtmp         , & ! total volume of new and old ice
                 area1        , & ! starting fractional area of existing ice
                 alvl         , & ! starting level ice area
                 rnilyr       , & ! real(nilyr)
                 dfresh       , & ! change in fresh
                 dfsalt       , & ! change in fsalt
                 vi0tmp       , & ! frzmlt part of frazil, LR for CESM
                 Ti               ! frazil temperature
              
              real (kind=dbl_kind), dimension (icells) :: &
                 qi0new       , & ! frazil ice enthalpy
                 Si0new           ! frazil ice bulk salinity

              real (kind=dbl_kind), dimension (icells,nilyr) :: &
                 Sprofile         ! salinity profile used for new ice additions

              integer (kind=int_kind) :: &
                 jcells, kcells     , & ! grid cell counters
                 ij, m , &                 ! combined i/j horizontal indices
                 kk

              integer (kind=int_kind), dimension (icells) :: &
                 indxij2,  indxij3  , & ! compressed i/j indices
                 indxi2, indxj2     , &
                 indxi3, indxj3

              character (len=char_len) :: &
                 fieldid           ! field identifier

              ! BGC
              real (kind=dbl_kind), dimension (nx_block,ny_block,ncat) :: &
                 eicen, &     ! energy of melting for each ice layer (J/m^2)
                 aicen_init, &    ! fractional area of ice
                 vicen_init       ! volume per unit area of ice (m)

              real (kind=dbl_kind), dimension (icells) :: &
                 vi0_init         ! volume of new ice
! LR
              real (kind=dbl_kind) :: frazil_conc ! for CESM
       
              real (kind=dbl_kind), dimension (nx_block,ny_block,ncat) :: &
                 area2, d_an_tot       ! change in the ITD due to lateral growth and new ice

              real (kind = dbl_kind), dimension (nx_block, ny_block) :: &
                delta_a1_ani            ! change in the first ITD cat
                                        ! due to new ice growth
         
              real (kind=dbl_kind), dimension (icells,ncat) :: &
                 ain0new     , &  ! area of new ice added to any thickness cat
                 vin0new          ! volume of new ice added to any thickness cat
 
              real (kind=dbl_kind), dimension (nx_block,ny_block,nfsd,ncat) :: &
                areal_mfstd_latg, & ! areal mFSTD after lateral growth
                areal_mfstd_init    ! initial areal mFSTD (tilda)

              real (kind=dbl_kind), dimension (nfsd) :: &
                 fin_diff, &         ! finite differences for G_r*tilda(L)
                 areal_mfstd_ni      ! areal mFSTD after new ice added

              real (kind=dbl_kind) :: &
                 totfrac, &      ! for FSD normalization
                 amount_taken              

              real (kind=dbl_kind), dimension (nfsd+1) :: &
                f_flx

               ! initialize vars      
               areal_mfstd_init = trcrn(:,:,nt_fsd:nt_fsd+nfsd-1,:)
               areal_mfstd_latg = areal_mfstd_init
               area2 = aicen

               ! initialize these to zero
               d_an_latg = c0
               d_an_addnew = c0        ! compute these regardless of diags
               d_an_tot = c0
               hsurp(:)  = c0
               hi0new = c0
               ai0new(:) = c0
               ain0new(:,:) = c0
               vin0new(:,:) = c0
               vi0new_lat = c0

               ! diagnostics returned like this if no growth occurs
               if (write_diag_diff) then
                       d_amfstd_latg = c0
                       d_amfstd_addnew = c0        
                       d_afsd_latg = c0
                       d_afsd_addnew = c0       
               end if 
! LR
              !-----------------------------------------------------------------
              ! initialize
              !-----------------------------------------------------------------
              l_stop = .false.
              istop = 0
              jstop = 0

              jcells = 0
              kcells = 0

             rnilyr = real(nilyr,kind=dbl_kind)

              if (ncat > 1) then
                 hi0max = hin_max(1)*0.9_dbl_kind  ! not too close to boundary
              else
                 hi0max = bignum                   ! big number
              endif

              ! for bgc
              aicen_init(:,:,:) = aicen(:,:,:)
              vicen_init(:,:,:) = vicen(:,:,:)

              if (l_conservation_check) then

              ! initial ice volume and energy in each grid cell
              eicen(:,:,:) = c0
              do n = 1, ncat
              do k = 1, nilyr
              do ij = 1, icells
                 i = indxi(ij)
                 j = indxj(ij)
                 eicen(i,j,n) = eicen(i,j,n) + trcrn(i,j,nt_qice+k-1,n) &
                              * vicen(i,j,n)/real(nilyr,kind=dbl_kind)
              enddo
              enddo
              enddo

              call column_sum (nx_block, ny_block,       &
                               icells,   indxi,   indxj, &
                               ncat,                     &
                               vicen,    vice_init)

              call column_sum (nx_block, ny_block,       &
                               icells,   indxi,   indxj, &
                               ncat,                     &
                               eicen,    eice_init)

              endif ! l_conservation_check

              !-----------------------------------------------------------------
              ! Compute average enthalpy of new ice.
              ! Sprofile is the salinity profile used when adding new ice to
              ! all categories, for ktherm/=2, and to category 1 for all ktherm.
              !
              ! NOTE:  POP assumes new ice is fresh!
              !-----------------------------------------------------------------

              if (ktherm == 2) then  ! mushy
        !DIR$ CONCURRENT !Cray
        !cdir nodep      !NEC
        !ocl novrec      !Fujitsu
                 do ij = 1, icells
                    i = indxi(ij)
                    j = indxj(ij)
                    if (sss(i,j) > c2 * dSin0_frazil) then
                       Si0new(ij) = sss(i,j) - dSin0_frazil
                    else
                       Si0new(ij) = sss(i,j)**2 / (c4*dSin0_frazil)
                    endif
                    do k = 1, nilyr
                       Sprofile(ij,k) = Si0new(ij)
                    enddo
                    Ti = min(liquidus_temperature_mush(Si0new(ij)/phi_init), -p1)
                    qi0new(ij) = enthalpy_mush(Ti, Si0new(ij))
                 enddo ! ij

              else

        !DIR$ CONCURRENT !Cray
        !cdir nodep      !NEC
        !ocl novrec      !Fujitsu
                 do ij = 1, icells
                    i = indxi(ij)
                    j = indxj(ij)
                    do k = 1, nilyr
                       Sprofile(ij,k) = salinz(i,j,k)
                    enddo
                    qi0new(ij) = -rhoi*Lfresh
                 enddo ! ij
              endif    ! ktherm

              !-----------------------------------------------------------------
              ! Compute the volume, area, and thickness of new ice.
              !-----------------------------------------------------------------
              ! LR
              vlateral = c0
   
        !DIR$ CONCURRENT !Cray
        !cdir nodep      !NEC
        !ocl novrec      !Fujitsu
              do ij = 1, icells
                 i = indxi(ij)
                 j = indxj(ij)

                 fnew = max (frzmlt(i,j), c0)       ! fnew > 0 iff frzmlt > 0
                 vi0new(ij) = -fnew*dt / qi0new(ij) ! note sign convention, qi < 0
                 vi0_init(ij) = vi0new(ij)          ! for bgc

                ! increment ice volume and energy
                 if (l_conservation_check) then
                    vice_init(ij) = vice_init(ij) + vi0new(ij)
                    eice_init(ij) = eice_init(ij) + vi0new(ij)*qi0new(ij)
                 endif

                 ! history diagnostics
                 !frazil(i,j) = vi0new(ij)

                 if (present(frz_onset) .and. present(yday)) then
                    if (frazil(i,j) > puny .and. frz_onset(i,j) < puny) &
                         frz_onset(i,j) = yday
                 endif

              !-----------------------------------------------------------------
              ! Update fresh water and salt fluxes.
              !
              ! NOTE: POP assumes fresh water and salt flux due to frzmlt > 0
              !       is NOT included in fluxes fresh and fsalt.
              !-----------------------------------------------------------------

                 if (update_ocn_f) then
                    dfresh = -rhoi*vi0new(ij)/dt 
                    dfsalt = ice_ref_salinity*p001*dfresh

                    fresh(i,j)      = fresh(i,j)      + dfresh
                    fsalt(i,j)      = fsalt(i,j)      + dfsalt
                 endif

! LR
		  ! Need to return mushy-layer frazil to POP
		 if (ktherm == 2 .and. .not.update_ocn_f) then
		    vi0tmp = fnew*dt / (rhoi*Lfresh)
		    frazil_diag(i,j) = frazil(i,j) - vi0tmp
		    dfresh = -rhoi*(vi0new(ij)-vi0tmp)/dt 
		    dfsalt = ice_ref_salinity*p001*dfresh

		    fresh(i,j)      = fresh(i,j)      + dfresh
		    fsalt(i,j)      = fsalt(i,j)      + dfsalt
		 endif
! LR
               
              !-----------------------------------------------------------------
              ! Decide how to distribute the new ice.
              !-----------------------------------------------------------------
                
  
                 if (vi0new(ij) > c0) then

                        ! LR partition
                        if (latsurf_area(i,j).gt.puny) then ! otherwise remains zero
                                vi0new_lat(ij)=(vi0new(ij)*lead_area(i,j)) / &
                                               (c1+(aice(i,j)/latsurf_area(i,j)))
                        end if

                        if (vi0new_lat(ij).lt.c0) stop 'latlt0'
                                               
                        ! LR, for diagnostics
                        vlateral(i,j) = vi0new_lat(ij) 
                        frazil(i,j) = vi0new(ij) - vi0new_lat(ij)
                   
                        !-----------------------------------------------------------------
                        ! LR: changes for lateral growth
                        !-----------------------------------------------------------------

                        if (vi0new_lat(ij).gt.puny) then
                                
                                G_radial(i,j)=vi0new_lat(ij)/dt

                                ! compute change to ITD      
                                do n=1,ncat

                                        if (aicen(i,j,n).gt.puny) then
                                             if (ABS(SUM(areal_mfstd_init(i,j,:,n))-c1).gt.1.0e-9_dbl_kind) then
                                                print *, SUM(areal_mfstd_init(i,j,:,n)), ABS(SUM(areal_mfstd_init(i,j,:,n))-c1)
                                                print *, &
                                    'WARNING init not normed, ani'
                                             end if
                                             ! in case of 10e-10 errors
                                             areal_mfstd_init(i,j,:,n) = areal_mfstd_init(i,j,:,n)/SUM(areal_mfstd_init(i,j,:,n))
                                        end if
    
                                        d_an_latg(i,j,n) = c0
                                        
                                        do k=1,nfsd ! sum over k
                                                d_an_latg(i,j,n) = d_an_latg(i,j,n) + (c2/floe_rad_c(k))*aicen(i,j,n)* &
                                                       areal_mfstd_init(i,j,k,n)*G_radial(i,j)*dt
                                        end do

                                        if (d_an_latg(i,j,n).lt.c0) stop &
                                                'delta itd lt0, lg'
                                end do ! n 
                                
                                if (SUM(d_an_latg(i,j,:)).ge.lead_area(i,j)) stop &
                                         'Filled up lead region'


                        endif ! vi0new_lat > 0
                        ! otherwise d_an_latg stays zero

                        !-------------------------------------------------------------------------
                        ! Now use remaining ice volume as in standard
                        ! model, but ice cannot grow into the area that
                        ! has grown laterally

                        vi0new(ij) = vi0new(ij) - vi0new_lat(ij)
                        if (vi0new(ij).lt.c0) stop 'neg vol'
                        if (vi0new(ij).gt.vi0_init(ij)) stop &
                                'increased vol'


                        !-----------------------------------------------------------------
                        ! Decide how to distribute the new ice.
                        !-----------------------------------------------------------------
                        if (lead_area(i,j).gt.aice0(i,j)) stop &
                         'leadarewrong'

                        hsurp(ij)  = c0
                        ai0new(ij) = c0
                        amount_taken = SUM(d_an_latg(i,j,:)) 
 
                        if (vi0new(ij) > c0) then

                                ! new ice area and thickness
                                ! hin_max(0) < new ice thickness < hin_max(1)
                                if ((aice0(i,j)-amount_taken) > puny) then
                                        hi0new = max(vi0new(ij)/(aice0(i,j)-amount_taken), hfrazilmin)
                                        if (hi0new > hi0max .and. (aice0(i,j)-amount_taken)+puny < c1) then
                                                ! distribute excess volume over all categories (below)
                                                if (aice(i,j).eq.c0) stop 'aice 0'
                                                hi0new = hi0max
                                                ai0new(ij) = (aice0(i,j) - amount_taken)
                                                vsurp      = vi0new(ij) - ai0new(ij)*hi0new
                                                hsurp(ij)  = vsurp / aice(i,j)
                                                vi0new(ij) = ai0new(ij)*hi0new
                                        else
                                                ! put ice in a single category, with hsurp = 0
                                                ai0new(ij) = vi0new(ij)/hi0new
                                        endif
                                else                ! aice0 < puny
                                        hsurp(ij) = vi0new(ij)/aice(i,j) ! new thickness in each cat
                                        vi0new(ij) = c0
                                endif               ! aice0 > puny
                        end if
 
                       !--------------------------------------------------------------------------
                       ! Combine things
                       !--------------------------------------------------------------------
                        ! diagnostics
                        d_an_addnew(i,j,1) = ai0new(ij)

                        ! volume added to each from lateral growth only
                        vin0new(ij,:) = c0
                        do n=1,ncat
                                if (aicen(i,j,n).gt.c0) &
                                vin0new(ij,n) = d_an_latg(i,j,n) * vicen(i,j,n)/aicen(i,j,n)
                        end do

                        ! altogether
                        d_an_tot(i,j,2:ncat) = d_an_latg(i,j,2:ncat)
                        d_an_tot(i,j,1) = d_an_latg(i,j,1) + d_an_addnew(i,j,1)
                        vin0new(ij,1) = vin0new(ij,1) +ai0new(ij)*hi0new 
               
                       if (SUM(d_an_tot(i,j,:)).gt.(puny+aice0(i,j))) stop &
                        'too much d_an_tot'

                       if (ANY(d_an_tot(i,j,:).lt.-puny)) stop &
                        'neg d_an_tot'
                 
                       if ((SUM(vin0new(ij,:)).le.c0).and.(hsurp(ij).le.c0)) &
                                stop 'no ice growth'
                 endif  ! vi0new > puny

                 !-----------------------------------------------------------------
                 ! Identify grid cells receiving new ice.
                 !-----------------------------------------------------------------

                 i = indxi(ij)
                 j = indxj(ij)

                 if (SUM(vin0new(ij,:)) > c0) then  ! area growth in all categories 
                        jcells = jcells + 1
                        indxi2(jcells) = i
                        indxj2(jcells) = j
                        indxij2(jcells) = ij
                 endif

                if (hsurp(ij) > c0) then   ! add ice to all categories
                        kcells = kcells + 1
                        indxi3(kcells) = i
                        indxj3(kcells) = j
                        indxij3(kcells) = ij
                endif

              enddo                  ! ij

              !-----------------------------------------------------------------
              ! Distribute excess ice volume among ice categories by increasing
              ! ice thickness, leaving ice area unchanged.
              !
              ! NOTE: If new ice contains globally conserved tracers
              !       (e.g., isotopes from seawater), code must be added here.
              !
              ! The mushy formulation (ktherm=2) puts the new ice only at the
              ! bottom of existing ice and adjusts the layers accordingly.
              ! The other formulations distribute the new ice throughout the 
              ! existing ice column.
              !-----------------------------------------------------------------

              do n = 1, ncat

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
                 do ij = 1, kcells
                    i = indxi3(ij)
                    j = indxj3(ij)
                    m = indxij3(ij)

                    vsurp = hsurp(m) * aicen(i,j,n)

                    ! update ice age due to freezing (new ice age = dt)
                    vtmp = vicen(i,j,n) + vsurp
                    if (tr_iage .and. vtmp > puny) &
                    trcrn(i,j,nt_iage,n) = &
                    (trcrn(i,j,nt_iage,n)*vicen(i,j,n) + dt*vsurp) / vtmp

                    if (tr_lvl .and. vicen(i,j,n) > puny) then
                        trcrn(i,j,nt_vlvl,n) = &
                        (trcrn(i,j,nt_vlvl,n)*vicen(i,j,n) + &
                        trcrn(i,j,nt_alvl,n)*vsurp) / vtmp
                    endif

                    if (tr_aero .and. vtmp > puny) then
                        do it = 1, n_aero
                                trcrn(i,j,nt_aero+2+4*(it-1),n) = &
                                trcrn(i,j,nt_aero+2+4*(it-1),n)*vicen(i,j,n) / vtmp
                                trcrn(i,j,nt_aero+3+4*(it-1),n) = &
                                trcrn(i,j,nt_aero+3+4*(it-1),n)*vicen(i,j,n) / vtmp
                        enddo
                    endif
! LR for CESM
		    if (tr_iso) then
		     do it=1,n_iso
		       if (it==1)   &
			  frazil_conc = isoice_alpha(c0,'HDO',frac)      &
					*HDO_ocn(i,j)
		       if (it==2)   &
			  frazil_conc = isoice_alpha(c0,'H2_16O',frac)   &
					*H2_16O_ocn(i,j)
		       if (it==3)   &
			  frazil_conc = isoice_alpha(c0,'H2_18O',frac)   &
					*H2_18O_ocn(i,j)

		       ! dilution in the ssl 
		       trcrn(i,j,nt_iso+2+4*(it-1),n)  &
			   = (trcrn(i,j,nt_iso+2+4*(it-1),n)*vicen(i,j,n)) &
			   / vtmp
		       ! dilution and uptake in the int
		       trcrn(i,j,nt_iso+3+4*(it-1),n)  &
			   = (trcrn(i,j,nt_iso+3+4*(it-1),n)*vicen(i,j,n) &
			   + frazil_conc*rhoi*vsurp) &
			   / vtmp

		       fiso_ocn(i,j,it) = fiso_ocn(i,j,it) &
			   - frazil_conc*rhoi*vsurp/dt
		     enddo
		    endif
! LR for CESM

                    ! update category volumes
                    vicen(i,j,n) = vtmp

                 enddo                  ! ij

                 if (ktherm == 2) then

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
                    do ij = 1, kcells
                       i = indxi3(ij)
                       j = indxj3(ij)
                       m = indxij3(ij)
               
                       vsurp = hsurp(m) * aicen(i,j,n)  ! note - save this above?
                       vtmp = vicen(i,j,n) - vsurp      ! vicen is the new volume
                       if (vicen(i,j,n) > c0) then
                                call update_vertical_tracers(trcrn(i,j,nt_qice:nt_qice+nilyr-1,n), &
                                                                vtmp, vicen(i,j,n), qi0new(m))
                                call update_vertical_tracers(trcrn(i,j,nt_sice:nt_sice+nilyr-1,n), &
                                                                vtmp, vicen(i,j,n), Si0new(m))
                       endif
                    enddo               ! ij

                 else ! ktherm

                        do k = 1, nilyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
                                do ij = 1, kcells
                                        i = indxi3(ij)
                                        j = indxj3(ij)
                                        m = indxij3(ij)

                                        ! factor of nilyr cancels out
                                        vsurp = hsurp(m) * aicen(i,j,n)  ! note - save this above?
                                        vtmp = vicen(i,j,n) - vsurp      ! vicen is the new volume
                                        if (vicen(i,j,n) > c0) then
                                                ! enthalpy
                                                trcrn(i,j,nt_qice+k-1,n) = &
                                                (trcrn(i,j,nt_qice+k-1,n)*vtmp + qi0new(ij)*vsurp) / vicen(i,j,n)
                                                ! salinity
                                                trcrn(i,j,nt_sice+k-1,n) = &
                                                (trcrn(i,j,nt_sice+k-1,n)*vtmp + Sprofile(ij,k)*vsurp) / vicen(i,j,n) 
                                        endif
                                enddo               ! ij
                        enddo               ! k

                 endif                  ! ktherm
              enddo                     ! n

             !-----------------------------------------------------------------
             ! Combine new ice grown in open water with category n ice.
             ! Assume that vsnon and esnon are unchanged.
             ! The mushy formulation assumes salt from frazil is added uniformly
             ! to category 1, while the others use a salinity profile.
             !-----------------------------------------------------------------
!LR
             do n=1,ncat

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
                do ij = 1, jcells
                        i = indxi2(ij)
                        j = indxj2(ij)
                        m = indxij2(ij)
                       
                       if ((d_an_tot(i,j,n).gt.c0).and.(vin0new(m,n).gt.c0)) then

                               if ((aicen(i,j,n).eq.c0).and.(vicen(i,j,n).gt.puny)) stop &
                                'aicen zero, vicen nonzero pre-update'

                                area1        = aicen(i,j,n)   ! save
                                vice1(ij)    = vicen(i,j,n)   ! save
                                aicen(i,j,n) = aicen(i,j,n) + d_an_tot(i,j,n)
                                aice0(i,j)   = aice0(i,j)   - d_an_tot(i,j,n)
                                vicen(i,j,n) = vicen(i,j,n) + vin0new(m,n)

                                if ((aicen(i,j,n).eq.c0).and.(vicen(i,j,n).gt.puny)) stop &
                                'aicen zero, vicen nonzero after update'

                                if ((aicen(i,j,n).ne.aicen(i,j,n)).or.(vicen(i,j,n).ne.vicen(i,j,n))) stop &
                                'aicen or vicen NaN after update in ani'

                                if (aicen(i,j,n).gt.c0) trcrn(i,j,nt_Tsfc,n) = &
                                (trcrn(i,j,nt_Tsfc,n)*area1 + Tf(i,j)*d_an_tot(i,j,n))/aicen(i,j,n)
                                trcrn(i,j,nt_Tsfc,n) = min (trcrn(i,j,nt_Tsfc,n), c0)

                                if (tr_FY) then
                                        if (aicen(i,j,n).gt.c0) trcrn(i,j,nt_FY,n) = &
                                        (trcrn(i,j,nt_FY,n)*area1 + d_an_tot(i,j,n))/aicen(i,j,n)
                                        trcrn(i,j,nt_FY,n) = min(trcrn(i,j,nt_FY,n), c1)
                                endif

                              !-----------------------------------------------------------------
                              ! Evolve mFSTD according to lateral growth
                              ! and growth of new ice (in first category)
                              !-----------------------------------------------------------------
                                if (d_an_latg(i,j,n).gt.puny) then ! lateral growth

                                        ! area after lateral growth and
                                        ! before new ice formation
                                        area2(i,j,n) = aicen_init(i,j,n) + d_an_latg(i,j,n)

                                        fin_diff(:) = c0 ! NB could stay zero if all in largest FS cat
                                        f_flx(:) = c0
                                        do k = 2, nfsd!+1
                                                f_flx(k) = G_radial(i,j) * areal_mfstd_init(i,j,k-1,n) / &
                                                                        floe_binwidth(k-1)
                                        end do
                                        do k = 1, nfsd
                                                fin_diff(k) = f_flx(k+1) - f_flx(k)
                                        end do
                                   
                                        if (ABS(SUM(fin_diff(:))).gt.puny) stop &
                                           'sum fnk diff not zero in lg'

                                        areal_mfstd_latg(i,j,:,n) = c0       
                                        do k = 1,nfsd
                                                areal_mfstd_latg(i,j,k,n) = &
                                                areal_mfstd_init(i,j,k,n) +   &
                                                dt * (  - fin_diff(k) + &
                                                c2 * G_radial(i,j) * areal_mfstd_init(i,j,k,n) * &
                                                (c1/floe_rad_c(k) - & 
                                                SUM(areal_mfstd_init(i,j,:,n)/floe_rad_c(:))) )
                                        end do
                                        
                                        if (ABS(SUM(areal_mfstd_init(i,j,:,n))-c1).gt.puny) stop &
                                            'init mFSTD not normed, lg' 
         
                                        if (ABS(SUM(areal_mfstd_latg(i,j,:,n))-c1).gt.puny) stop &
                                                'mFSTD not normed, lg' 

                                        ! just in case (may be errors < 1e-11)
                                        areal_mfstd_latg(i,j,:,n) = areal_mfstd_latg(i,j,:,n)/SUM(areal_mfstd_latg(i,j,:,n))

                                        if (ANY(areal_mfstd_latg(i,j,:,n).lt.c0)) stop &
                                                'neg mFSTD, lg'

                                        trcrn(i,j,nt_fsd:nt_fsd+nfsd-1,n) = areal_mfstd_latg(i,j,:,n)
 
                                        if (write_diag_diff) d_amfstd_latg(i,j,:,n) = & 
                                         areal_mfstd_latg(i,j,:,n)-  areal_mfstd_init(i,j,:,n)

                                else
                                        areal_mfstd_latg(i,j,:,n) = areal_mfstd_init(i,j,:,n) 
                                end if
                                
                                if (n.eq.1) then

                                   ! add new frazil ice to smallest
                                   ! thickness and floe size
                                   if (d_an_addnew(i,j,n).gt.puny) then   
                                     
                                        if (d_an_addnew(i,j,n).gt.aicen(i,j,n)) stop &
                                           'area update neg somewhere'      

                                       areal_mfstd_ni(:) = c0
                                        if (SUM(areal_mfstd_latg(i,j,:,n)).gt.puny) then ! FSD exists

                                            if (new_ice_fs.eq.0) then
                                                ! grow in smallest
                                                areal_mfstd_ni(1) =  (area2(i,j,n) * areal_mfstd_latg(i,j,1,n) + &
                                                                         ai0new(m))/(area2(i,j,n)+ai0new(m))                                                        

                                                do k=2,nfsd  ! diminish other floe cats accordingly
                                                        areal_mfstd_ni(k) = areal_mfstd_latg(i,j,k,n) * &
                                                              area2(i,j,n)/(area2(i,j,n)+ai0new(m))
                                                enddo 

                                            else if (new_ice_fs.eq.1) then
                                                ! grow in largest
                                                areal_mfstd_ni(nfsd) =  (area2(i,j,n) * areal_mfstd_latg(i,j,nfsd,n) + &
                                                                         ai0new(m))/(area2(i,j,n)+ai0new(m))

                                                do k=1,nfsd-1  ! diminish other floe cats accordingly
                                                        areal_mfstd_ni(k) = areal_mfstd_latg(i,j,k,n) * &
                                                              area2(i,j,n)/(area2(i,j,n)+ai0new(m))
                                                end do

                                            else if (new_ice_fs.ge.2) then
                                                if (new_ice_fs.eq.2) then
                                                    ! check open water area 
                                                    new_size = nfsd      
                                                    do k = 1,nfsd
                                                      if (aicen(i,j,n)*tarea(i,j).lt.floe_area_l(k)) then
                                                          new_size = k - 1
                                                          EXIT
                                                      end if
                                                    end do
                                                    new_size = MAX(new_size,1)
                                             
                                                else if (new_ice_fs.eq.3) then
                                                    ! wave depedent size
                                                    call wave_dep_growth(wave_spectrum(i,j,:), new_size) 
                                                end if

                                                ! grow in new_size
                                                areal_mfstd_ni(new_size) =  (area2(i,j,n) * areal_mfstd_latg(i,j,new_size,n) + &
                                                                         ai0new(m))/(area2(i,j,n)+ai0new(m))

                                                do k=1,new_size-1  ! diminish other floe cats accordingly
                                                        areal_mfstd_ni(k) = areal_mfstd_latg(i,j,k,n) * &
                                                              area2(i,j,n)/(area2(i,j,n)+ai0new(m))
                                                end do

                                                do k=new_size+1,nfsd  ! diminish other floe cats accordingly
                                                        areal_mfstd_ni(k) = areal_mfstd_latg(i,j,k,n) * &
                                                              area2(i,j,n)/(area2(i,j,n)+ai0new(m))
                                                end do

                                            end if ! new_fs_option
                                        else ! entirely new ice
                                             if (new_ice_fs.eq.0) then
                                                 areal_mfstd_ni(1) = c1
                                             else if (new_ice_fs.eq.1) then
                                                 areal_mfstd_ni(nfsd) = c1
                                             else if (new_ice_fs.eq.2) then
                                                 ! check open water area
                                                 new_size = nfsd
                                                 do k = 1,nfsd
                                                      if (aicen(i,j,n)*tarea(i,j).lt.floe_area_l(k)) then
                                                          new_size = k
                                                          EXIT
                                                      end if
                                                 end do
                                                 new_size = MAX(new_size,1)
                                                 areal_mfstd_ni(new_size) = c1
                                             else if  (new_ice_fs.eq.3) then
                                                 ! wave depedent size
                                                 call wave_dep_growth(wave_spectrum(i,j,:), new_size)
                                                 areal_mfstd_ni(new_size) = c1
                                             end if      
                                        end if ! entirely new ice 

                                        if (ABS(SUM(areal_mfstd_ni)-c1).gt.puny) then
                                                 print *, 'areal_mfstd_ni',areal_mfstd_ni
                                                 print *, ABS(SUM(areal_mfstd_ni)-c1)
                                                 print *, 'mFSTD not normed, ni'
                                        end if 
                                        areal_mfstd_ni = areal_mfstd_ni / SUM(areal_mfstd_ni)

                                        if (ANY(areal_mfstd_ni.lt.c0)) stop &
                                                'neg mFSTD, ni'

                                        trcrn(i,j,nt_fsd:nt_fsd+nfsd-1,n) = areal_mfstd_ni
                                        if (SUM(trcrn(i,j,nt_fsd:nt_fsd+nfsd-1,n)).lt.puny) stop 'should not be punyy'  

                                        if (write_diag_diff) d_amfstd_addnew(i,j,:,n) = &
                                          trcrn(i,j,nt_fsd:nt_fsd+nfsd-1,n) -  areal_mfstd_latg(i,j,:,n)
                                       
                                        if (ANY(trcrn(i,j,nt_fsd:nt_fsd+nfsd-1,n).gt.c1+puny)) stop  &
                                         'mFSTD > 1 in ani'


                                   end if ! d_an_addnew > puny

                                endif ! n=1

                               if (vicen(i,j,n) > puny) then
                                        if (tr_iage) &
                                        trcrn(i,j,nt_iage,n) = &
                                        (trcrn(i,j,nt_iage,n)*vice1(ij) + dt*vin0new(m,n))/vicen(i,j,n)

                                        if (tr_aero) then
                                                do it = 1, n_aero
                                                        trcrn(i,j,nt_aero+2+4*(it-1),n) = &
                                                        trcrn(i,j,nt_aero+2+4*(it-1),n)*vice1(ij)/vicen(i,j,n)
                                                        trcrn(i,j,nt_aero+3+4*(it-1),n) = &
                                                        trcrn(i,j,nt_aero+3+4*(it-1),n)*vice1(ij)/vicen(i,j,n)
                                                enddo
                                        endif

                                        if (tr_iso) then
					      do it=1,n_iso
						 if (it==1)   &
						    frazil_conc = isoice_alpha(c0,'HDO',frac)      &
								  *HDO_ocn(i,j)
						 if (it==2)   &
						    frazil_conc = isoice_alpha(c0,'H2_16O',frac)   &
								  *H2_16O_ocn(i,j)
						 if (it==3)       &
						    frazil_conc = isoice_alpha(c0,'H2_18O',frac)   &
								*H2_18O_ocn(i,j)

						trcrn(i,j,nt_iso+2+4*(it-1),n) = &
						  (trcrn(i,j,nt_iso+2+4*(it-1),n)*vice1(ij))/vicen(i,j,n)
						trcrn(i,j,nt_iso+3+4*(it-1),n) = &
						  (trcrn(i,j,nt_iso+3+4*(it-1),n)*vice1(ij) &
						  +frazil_conc*rhoi*vin0new(m,n))/vicen(i,j,n)

						fiso_ocn(i,j,it) = fiso_ocn(i,j,it) &
						  - frazil_conc*rhoi*SUM(vin0new(m,:))/dt
					      enddo
					endif


                                        if (tr_lvl) then
                                                alvl = trcrn(i,j,nt_alvl,n)
                                                trcrn(i,j,nt_alvl,n) = &
                                                (trcrn(i,j,nt_alvl,n)*area1 + d_an_tot(i,j,n))/aicen(i,j,n)
                                                trcrn(i,j,nt_vlvl,n) = &
                                                (trcrn(i,j,nt_vlvl,n)*vice1(ij) + vin0new(m,n))/vicen(i,j,n)
                                        endif

                                        if (tr_pond_cesm .or. tr_pond_topo) then
                                                trcrn(i,j,nt_apnd,n) = &
                                                trcrn(i,j,nt_apnd,n)*area1/aicen(i,j,n)
                                        elseif (tr_pond_lvl) then
                                                if (trcrn(i,j,nt_alvl,n) > puny) then
                                                        trcrn(i,j,nt_apnd,n) = &
                                                        trcrn(i,j,nt_apnd,n) * alvl*area1 &
                                                        / (trcrn(i,j,nt_alvl,n)*aicen(i,j,n))
                                                endif
                                        endif
                                endif
                        endif
                enddo                     ! ij

                do k = 1, nilyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
                        do ij = 1, jcells
                                i = indxi2(ij)
                                j = indxj2(ij)
                                m = indxij2(ij)
         
                                     if ((vin0new(m,n).gt.c0).and.(d_an_tot(i,j,n).gt.c0) ) then
                                        if (vicen(i,j,n) > c0) then
                                                ! factor of nilyr cancels out
                                                ! enthalpy
                                                trcrn(i,j,nt_qice+k-1,n) = &
                                                (trcrn(i,j,nt_qice+k-1,n)*vice1(ij) &
                                                + qi0new(m)*vin0new(m,n))/vicen(i,j,n)
                                                ! salinity
                                                trcrn(i,j,nt_sice+k-1,n) = &
                                                (trcrn(i,j,nt_sice+k-1,n)*vice1(ij) &
                                                + Sprofile(m,k)*vin0new(m,n))/vicen(i,j,n)
                                        endif
                                    end if
                        enddo ! ij
                enddo ! k

             enddo !n

            ! sanity
            do n = 1,ncat
            do ij = 1, jcells
                i = indxi2(ij)
                j = indxj2(ij)
                if (aicen(i,j,n).gt.puny) then
                        if (ABS(SUM(trcrn(i,j,nt_fsd:nt_fsd+nfsd-1,n))-c1).gt.puny) then 
                                print *, SUM(trcrn(i,j,nt_fsd:nt_fsd+nfsd-1,n))
                                print *, ABS(SUM(trcrn(i,j,nt_fsd:nt_fsd+nfsd-1,n))-c1)
                                print *, trcrn(i,j,nt_fsd:nt_fsd+nfsd-1,n)
                                print *, 'n= ',n
                                print *, vi0new(ij), vi0new_lat(ij)
                                print *, d_amfstd_latg(i,j,:,n)
                                print *, d_amfstd_addnew(i,j,:,n)
                                print *, 'd an ',d_an_tot(i,j,n), d_an_latg(i,j,n), d_an_addnew(i,j,n)
                                print *, 'WARNING not norm after growth'
                        end if
                else
                        trcrn(i,j,nt_fsd:nt_fsd+nfsd-1,n) = c0
                end if
            end do
            end do

            ! diagnostics
            if (write_diag_diff) then
                    do k=1,nfsd
                    d_afsd_latg(:,:,k) = c0
                    d_afsd_addnew(:,:,k) = c0
                    do n=1,ncat
                        d_afsd_latg(:,:,k) = d_afsd_latg(:,:,k)  + &
                                (aicen_init(:,:,n)+d_an_latg(:,:,n))*areal_mfstd_latg(:,:,k,n) - &
                                aicen_init(:,:,n)*areal_mfstd_init(:,:,k,n)

                        d_afsd_addnew(:,:,k) = d_afsd_addnew(:,:,k)  + &
                                aicen(:,:,n)*trcrn(:,:,nt_fsd+k-1,n) - &
                                (aicen_init(:,:,n)+d_an_latg(:,:,n))*areal_mfstd_latg(:,:,k,n)
                    end do
                    end do
            end if
             !----------------------

             if (l_conservation_check) then

             ! initial ice volume in each grid cell
             eicen(:,:,:) = c0
             do n = 1, ncat
             do k = 1, nilyr
             do ij = 1, icells
                i = indxi(ij)
                j = indxj(ij)
                eicen(i,j,n) = eicen(i,j,n) + trcrn(i,j,nt_qice+k-1,n) &
                      * vicen(i,j,n)/real(nilyr,kind=dbl_kind)
             enddo
             enddo
             enddo

             call column_sum (nx_block, ny_block,       &
                       icells,   indxi,   indxj, &
                       ncat,                     &
                       vicen,    vice_final)

             call column_sum (nx_block, ny_block,       &
                       icells,   indxi,   indxj, &
                       ncat,                     &
                       eicen,    eice_final)

             fieldid = 'vice, add_new_ice'
             call column_conservation_check (nx_block,  ny_block,      &
                                      icells,   indxi,   indxj, &
                                      fieldid,                  &
                                      vice_init, vice_final,    &
                                      puny,      l_stop,        &
                                      istop,     jstop)

             fieldid = 'eice, add_new_ice'
             call column_conservation_check (nx_block,  ny_block,      &
                                      icells,   indxi,   indxj, &
                                      fieldid,                  &
                                      eice_init, eice_final,    &
                                      puny*Lfresh*rhoi, l_stop, &
                                      istop,     jstop)
             if (l_stop) return

             endif ! l_conservation_check

             !-----------------------------------------------------------------
             ! Biogeochemistry
             !-----------------------------------------------------------------     
             if (tr_brine .or. skl_bgc) &
                call add_new_ice_bgc (nx_block,  ny_block,   dt,       &
                           icells,     jcells,     kcells,   &
                           indxi,      indxj,                &
                           indxi2,     indxj2,     indxij2,  &
                           indxi3,     indxj3,     indxij3,  &
                           aicen_init, vicen_init, vi0_init, &
                           aicen,      vicen,      vi0new,   &
                           ntrcr,      trcrn,      nbtrcr,   &
                           sss,        ocean_bio,  flux_bio, &
                           hsurp,                            &
                           l_stop,     istop,      jstop)

      end subroutine add_new_ice_lat



!=======================================================================

!=======================================================================

        subroutine floe_merge_thermo(iblk,nx_block, ny_block, &
                                    ntrcr, icells, indxi, indxj, & 
                                    dt, &
                                    aice, aicen, frzmlt, areal_mfstd, &
                                    d_amfstd_merge, d_afsd_merge)

        use ice_fsd, only: area_scaled_h, &  ! and no binwidth is greater than 1
                           area_scaled_c, &  ! (dimensionless)
                           area_scaled_binwidth, &
                           alpha_mrg, &      ! defines floe combinations
                           c_mrg ! units (s^-1)
                        

      integer (kind=int_kind), intent(in) :: &
         iblk, &
         nx_block, ny_block, & ! block dimensions
         ntrcr             , & ! number of tracers in use
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
         aice  , & ! total concentration of ice
         frzmlt    ! freezing/melting potential (W/m^2)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nfsd,ncat), &
         intent(inout) :: &
         areal_mfstd, &
         d_amfstd_merge

      real (kind=dbl_kind), dimension (nx_block,ny_block,nfsd), &
         intent(inout) :: &
         d_afsd_merge

      ! local variables

      integer (kind=int_kind) :: &
        t, &
        i, j, n, k, ij, m, &
        kx, ky, kz, a

      real (kind=dbl_kind), dimension(nfsd) :: &
        amfstd_init, amfstd_tmp, coag_pos, coag_neg

      real(kind=dbl_kind) :: &
        subdt, &
        area_loss, &    !
        area_loss_mcat, &
        stability       ! what needs to be one to satisfy
                        ! stability condition for Smol. eqn.

      integer(kind=int_kind) :: &
        ndt_mrg         ! number of sub-timesteps required to satisfy
                        ! stability condition for Smol. eqn.


        d_afsd_merge = c0
        d_amfstd_merge = c0
       
        do n=1,ncat
                do ij = 1, icells
           
                        i = indxi(ij)
                        j = indxj(ij)
        
                      
                        !-----------------------------------------------------------------
                        ! If there is some ice in the lower (nfsd-1) categories
                        ! and there is freezing potential
                        !-----------------------------------------------------------------
                        if ((frzmlt(i,j).gt.puny).and.(aicen(i,j,n).gt.p1)) then
                       
                               ! time step limitations for merging
                                stability = dt * c_mrg * aicen(i,j,n) * area_scaled_h(nfsd)
                                ndt_mrg = NINT(stability+p5) ! add .5 to round up                        
                                subdt = dt/FLOAT(ndt_mrg)
 
                                amfstd_init(:) = areal_mfstd(i,j,:,n)
                                amfstd_tmp = amfstd_init

                                if (ABS(SUM(amfstd_init) - c1).gt.puny) stop 'not 1 b4 mrg'
                                if (ANY(amfstd_init.lt.c0-puny)) stop &
                                 'negative mFSTD b4 mrg'
                                if (ANY(amfstd_init.gt.c1+puny)) stop &
                                 'mFSTD>1 b4 mrg'

                                area_loss_mcat = c0
                                do t = 1, ndt_mrg
                                     do kx = 1, nfsd

                                         coag_pos(kx) = c0
                                         do ky = 1, kx
                                             a = alpha_mrg(kx,ky)
                                             coag_pos(kx) = coag_pos(kx) + &
                                                            area_scaled_c(ky) * amfstd_tmp(ky) * aicen(i,j,n) * ( &
                                                            SUM(amfstd_tmp(a:nfsd)) + &
                                                            (amfstd_tmp(a-1)/area_scaled_binwidth(a-1)) * ( &
                                                            area_scaled_h(a-1) - area_scaled_h(kx) + area_scaled_c(ky) ))
       
                                         end do
                                     end do

                                     coag_neg(1) = c0
                                     coag_neg(2:nfsd) = coag_pos(1:nfsd-1)

                                     amfstd_tmp = amfstd_tmp - subdt*c_mrg*(coag_pos - coag_neg)

                                     if (ANY(amfstd_tmp.lt.c0-puny)) then
                                            print *, 'amfstd_init ',amfstd_init
                                            print *, 'coag_pos',coag_pos
                                            print *, 'coag_neg',coag_neg
                                            print *, 'amfstd_tmp ',amfstd_tmp
                                            print *, &
                                      'WARNING negative mFSTD mrg, l'
                                     end if


                                     if (ANY(amfstd_tmp.lt.c0-puny)) &
                                            stop 'negative mFSTD mrg, l'

                                     if (ANY(amfstd_tmp.gt.c1+puny)) &
                                            stop ' mFSTD> 1 mrg, l'

                                     if (ANY(dt*c_mrg*coag_pos.lt.-puny)) &
                                         stop 'not positive'

                                     area_loss_mcat = area_loss_mcat + subdt*c_mrg*coag_pos(nfsd)

                                end do
                                
                                ! ignore loss in largest cat
                                amfstd_tmp(nfsd) = amfstd_tmp(nfsd) + area_loss_mcat

                                area_loss = SUM(amfstd_init) - SUM(amfstd_tmp)

                                if (area_loss.lt.-puny) &
                                    stop 'area gain'

                                if (ABS(area_loss).gt.puny) &
                                    stop 'area change after correction'
                                
                                ! in case of small numerical errors
                                areal_mfstd(i,j,:,n) = amfstd_tmp/SUM(amfstd_tmp)

                                if (ANY(areal_mfstd(i,j,:,n).lt.-puny)) stop 'neg, mrg'

                                WHERE(areal_mfstd(i,j,:,n).lt.c0) areal_mfstd(i,j,:,n) = c0
                                
                                if (areal_mfstd(i,j,1,n).gt.amfstd_init(1)+puny) & 
                                    stop 'gain in smallest cat'

                                if (write_diag_diff) &
                                        d_amfstd_merge(i,j,:,n) = areal_mfstd(i,j,:,n) - amfstd_init

                        end if
               end do !ij

        end do! n

        if (write_diag_diff) then
            do k=1,nfsd
                d_afsd_merge(:,:,k) = c0
                do n=1,ncat
                        d_afsd_merge(:,:,k) = d_afsd_merge(:,:,k)  + &
                        aicen(:,:,n)* d_amfstd_merge(:,:,k,n)
                end do
            end do
        end if

           end subroutine floe_merge_thermo

!=======================================================================
!

      end module ice_fsd_thermo

!=======================================================================
