! Lettie Roach, NIWA/VUW, June 2016
! Some of this modified from Chris Horvat's Matlab original
! New module for CICE to generate waves and break the FSD
!
      module ice_wavebreaking

      use ice_kinds_mod
      use ice_constants
      use ice_domain_size, only: nfsd,ncat,max_ntrcr, max_blocks
      use ice_blocks, only: nx_block, ny_block
      use ice_fsd, only: floe_rad_c, floe_binwidth

      implicit none
      private
      public :: wave_break, init_wave, get_fraclengths, find_wave

      integer (kind=int_kind), private :: & 
         iblk        , & ! block index 
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         n, k
     
      integer (kind=int_kind), save, public :: &
         nwcat =10! number of wave, floe property categories
        
      real (kind=dbl_kind), dimension(ncat), save, public :: &
         dh, h         ! width, midpoint of thickness categories (m)

      real (kind=dbl_kind), dimension(10), save, public :: &
         nfloe_c, hs_c, tz_c, & ! category midpoints
         ! limits calculated from u-aj881 in Jupyter/wave-cats, 21/2/2017
         hs_l = (/ 0., 1.55382597 , 2.01547122 , 2.37775707 , 2.69950223 , 2.97230601 , 3.33214235 , 3.78572774 ,4.30519629 , 5.19086981/), &
         hs_h = (/ 1.55382597 , 2.01547122 , 2.37775707 , 2.69950223 , 2.97230601 , 3.33214235 , 3.78572774 , 4.30519629 , 5.19086981 , 12.65346622/), &
         tz_l = (/ 0. , 7.24649715 , 7.99813461 , 8.49258804 , 8.90908432 , 9.27624321 , 9.6638813 , 10.08105183 , 10.59714985 , 11.28474522/), &
         tz_h = (/ 7.24649715 , 7.99813461 , 8.49258804 , 8.90908432 , 9.27624321 , 9.6638813 , 10.08105183 , 10.59714985 , 11.28474522 , 15.27169895/), &
         nfloe_l = (/-5.59100181e-06, 5.59100181e-06, 2.76497868e-03, 2.04329416e-01,   6.44558029e+01,   1.89945813e+03,   7.66608594e+03,   1.79483105e+04,   3.47952930e+04,   6.06928242e+04/), &
         nfloe_h = (/ 5.59100181e-06, 2.76497868e-03, 2.04329416e-01,   6.44558029e+01,   1.89945813e+03,   7.66608594e+03,   1.79483105e+04,   3.47952930e+04,   6.06928242e+04, 1.000000e+5/)

        integer, parameter  :: &
                Nl=800, &       ! number of wavelengths
                loopcts=1       ! number of times to loop SSH

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
 
       real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks), save, public :: &
         nfl

       character(char_len_long), public :: &         
         wave_fn_dir    ! directory containing look-up tables
 
       logical (kind=log_kind), public :: &
         calc_wave      ! if true, calculate look tables for omega and fsdformed
                        ! if false, read them in from text file
      
!=======================================================================

      contains

!=======================================================================
! Author: Lettie Roach, NIWA, 2016
!
! Define discrete gridcell properties such as concentration, floe radius
! categories. Call wave_change for each of 5 concentration, NCAT thickness, NFSD
! floe size and 10 cumulative distance from ice edge categories (3000
! combinations). Save results for omega and fsdformed (change in FSD functions)
! in large arrays to be called later.
!
       subroutine init_wave
 
       use ice_communicate, only: my_task, master_task                      
       use ice_flux, only: fracture_histogram, wave_tau
       use ice_fsd, only: floe_binwidth
       use ice_itd, only: hin_max_init
       use ice_broadcast, only: broadcast_array

       integer (kind=int_kind) :: &
         iblk, &
         i, ct,            &
         nf, nd, nh, nt,nn,   &
         nk           , & ! thickness category index
         kk, ks               ! floe size index
      
       real (kind=dbl_kind)  :: &
         dmn

 
      !-----------------------------------------------------------------
      ! Define discrete properties for look-up table
      !-----------------------------------------------------------------

      dh=c0
      h=c0

      do n=1,ncat
                dh(n) = hin_max_init(n) - hin_max_init(n-1)
                h(n) = hin_max_init(n) - (dh(n)/2) ! midpoint of cats
      enddo

      do i=1,nwcat
        nfloe_c(i) = nfloe_l(i) + (nfloe_h(i)- nfloe_l(i))/c2
        hs_c(i) = hs_l(i) + (hs_h(i)- hs_l(i))/c2
        tz_c(i) = tz_l(i) + (tz_h(i)- tz_l(i))/c2
      end do

      !-----------------------------------------------------------------
      ! Produce look up tables or read them in
      !---------------------------------------------------------------
     
      wave_tau = c0
      fracture_histogram = c0

      if (my_task==master_task) print *, 'calc_wave=',calc_wave

      if (calc_wave) then
 
           call wave_change(h(:),nfloe_c(:),hs_c(:),tz_c(:), &
                            fracture_histogram(:,:,:,:,:), &
                            wave_tau(:,:,:,:) )
              

           ! write text files
           open(4,file= trim(wave_fn_dir)//&
                                "hist_lookup_new.txt")
 
           open(5,file= trim(wave_fn_dir)//&
                                "tau_lookup_new.txt")

           do nn=1,ncat
           do nd=1,nwcat                                     
           do nt=1,nwcat
           do nh=1,nwcat

                write(4,*) ( fracture_histogram(nh,nt,nd,nn,k), k=1,nfsd )
                write(5,*) wave_tau(nh,nt,nd,nn)
           enddo
           enddo
           enddo
           enddo
        close(4)
        close(5)

           print *, 'wrote text files'



        else
           if (trim(wave_fn_dir).eq.'') then
               print *, 'NO WAVE FRACTURE WILL OCCUR'
               fracture_histogram = c0
               wave_tau = c0
           else


               open(4,file = trim(wave_fn_dir)//&
                                "hist_lookup_roach2017.txt")
               open(5,file= trim(wave_fn_dir)//&
                                "tau_lookup_roach2017.txt")

               do nn=1,ncat
               do nd=1,nwcat                                     
               do nt=1,nwcat
               do nh=1,nwcat

                    read(4,*) ( fracture_histogram(nh,nt,nd,nn,k), k=1,nfsd )
                    read(5,*) wave_tau(nh,nt,nd,nn)

               enddo
               enddo
               enddo
               enddo
               close(4)
               close(5)
           end if
        end if

      
          end subroutine init_wave

!=======================================================================
! Author: Lettie Roach, NIWA, 2016
!
! For each ice-covered or freezing grid cell, look for nearest north/south cell 
! with wave data. If cell is freezing, save spectrum. If cell is ice-covered, 
! count floes to wave data

        subroutine find_wave

       use ice_blocks, only: block, get_block, nx_block, ny_block
       use ice_broadcast, only: broadcast_array
       use ice_communicate, only: my_task, master_task
       use ice_domain, only: distrb_info, nblocks, blocks_ice, ns_boundary_type, ew_boundary_type
       use ice_domain_size, only: nx_global, ny_global, max_blocks, nilyr
       use ice_flux, only: wave_hs, wave_tz, &
                           ice_search_i, ice_search_j, &
                           wave_search_i, wave_search_j, & 
                           cml_nfloes, frzmlt, &
                           nearest_wave_hs, nearest_wave_tz, &
                           wave_spectrum, dfreq, freq
       use ice_fsd, only: floe_rad_c, floe_area_c, write_diag_wave
       use ice_gather_scatter, only: gather_global, scatter_global
       use ice_grid, only: dyt, & ! height of T cell (m)
                           tarea, & ! area of grid cells (m^2)
                           TLAT, TLON, &    ! latitude, longitude of centre of Tcell pts (radians)
                           hm ! land mask (1 for ocean, 0 for land)
       use ice_itd, only: hin_max_init
       use ice_state, only: nt_fsd, aice, aicen, vice, vicen, trcrn, nt_qice

     !------local variables---------------------------------
      type (block) :: &
         this_block      ! block information for current block

      integer (kind=int_kind) :: &
         i, j, ij, ii, jj, &
         i_g, j_g

      integer (kind=int_kind) :: &
         jcells          ! number of cells with aice>0.01 or frzmlt>0

      integer (kind=int_kind), dimension(nx_block*ny_block) :: &
         jndxi, jndxj    ! compressed indices for jcells

       real (kind=dbl_kind), dimension(nx_global,ny_global) :: &
         aice_g,  tlat_g, tlon_g, &
         wave_hs_g, wave_tz_g,    &
         nfl_g

       logical (kind=log_kind), dimension(nx_global, ny_global) :: &
        tmp_mask

       integer (kind=int_kind), dimension(nx_global,ny_global) :: &
         int_tlat_g, int_tlon_g, hm_g

       integer (kind=int_kind), dimension(nx_block,ny_block,max_blocks) :: &
         hm_int, my_tlat_int, my_tlon_int

      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks) :: &
         tlat_deg, tlon_deg

       logical (kind=log_kind) :: &
                found_wave          ! if true, calculate properties

       real (kind=dbl_kind) :: &
                meanfs, mean_diam

       real (kind=dbl_kind), dimension(nfsd) :: &
        fsdperice, inv_diam

       integer (kind=int_kind), dimension(2) :: &
        found

       ! for spectrum
       real (kind=dbl_kind), dimension(Nl) :: &
                lambda, T, S_B, alpha
      !-----------------------------------------------------------------
      ! Compute number of floes (and tlat/tlon in degrees)
      !-----------------------------------------------------------------
     if (trim(wave_fn_dir).ne.'') then 

     ! wavelengths span
     do j=1,Nl
        lambda(j) = c2*floe_rad_c(1) + j*floe_rad_c(1)
     end do

     ! periods from wavelengths	
     T=(c2*pi*lambda/gravit)**p5

     if (allocated(wave_spectrum)) deallocate(wave_spectrum) 
     allocate(wave_spectrum(nx_block, ny_block, Nl, max_blocks))

     if (allocated(dfreq)) deallocate(dfreq)
     allocate(dfreq(Nl))

     dfreq(:) = sqrt(gravit/(c2*pi*floe_rad_c(1)))

     if (allocated(freq)) deallocate(freq)
     allocate(freq(Nl))
     freq = 1/T

     !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block,jcells,jndxi,jndxj,mean_diam)
     do iblk = 1, nblocks

        ! convert rad to deg
        tlat_deg(:,:,iblk) = tlat(:,:,iblk) * rad_to_deg
        tlon_deg(:,:,iblk) = tlon(:,:,iblk) * rad_to_deg
       

        ! give nfl dummy value for non-ice covered cells
        nfl(:,:,iblk) = -555.0_dbl_kind

        this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

            jcells = 0
            do j = jlo, jhi
            do i = ilo, ihi
               if (aice(i,j,iblk) > p01) then
                  jcells = jcells + 1
                  jndxi(jcells) = i
                  jndxj(jcells) = j
               endif
        
            enddo               ! i
            enddo               ! j
       
           do ij=1,jcells
               i = jndxi(ij)
               j = jndxj(ij)

               mean_diam = c0
               do k=1,nfsd
                   do n=1,ncat
                       mean_diam = mean_diam + &
                                   2 * floe_rad_c(k) *  aicen(i,j,n,iblk) * trcrn(i,j,nt_fsd+k-1,n,iblk)
                   end do                
               end do

               ! c/2mean(r) is number of floes per unit distance
               ! dyt is length of grid cell north-south
               if (mean_diam.gt.c0) nfl(i,j,iblk) = dyt(i,j,iblk)*aice(i,j,iblk)/mean_diam
               
               if (nfl(i,j,iblk).lt.puny) then
                   !if (my_task.eq.master_task) then
                       print *, 'nfl ',nfl(i,j,iblk)
                       print *, dyt(i,j,iblk),aice(i,j,iblk),mean_diam
                       print *, 'fsd ',trcrn(i,j,nt_fsd:nt_fsd+nfsd-1,:,iblk)
                   !end if
                   stop 'neg nfl'
               end if
           end do
          
           ! give nfl another dummy value for cells with land 
           WHERE (hm(:,:,iblk).lt.p5) nfl(:,:,iblk) = -999.0_dbl_kind


      end do !iblk 
     !$OMP END PARALLEL DO

     !-----------------------------------------------------------------
     ! Make these global arrays
     !-----------------------------------------------------------------
 
     call gather_global(wave_hs_g, wave_hs, master_task, distrb_info)
     call gather_global(wave_tz_g, wave_tz, master_task, distrb_info)
     call gather_global(nfl_g, nfl, master_task, distrb_info)
     call gather_global(tlat_g, tlat*rad_to_deg, master_task, distrb_info)
     call gather_global(tlon_g, tlon*rad_to_deg, master_task, distrb_info)

     ! broadcast global arrays from master task to all other tasks
     call broadcast_array(nfl_g, master_task)
     call broadcast_array(wave_hs_g, master_task)
     call broadcast_array(wave_tz_g, master_task)
     call broadcast_array(tlat_g, master_task)
     call broadcast_array(tlon_g, master_task)

      ice_search_i(:,:,:) = c0
      ice_search_j(:,:,:) = c0

!$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block, found, i_g, j_g, tmp_mask, alpha)
        do iblk=1, nblocks 
         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
         do i = ilo, ihi

         nearest_wave_hs(i,j,iblk) = c0
         nearest_wave_tz(i,j,iblk) = c0 
         cml_nfloes(i,j,iblk) = c0

         if (write_diag_wave) then
                 wave_search_i(i,j,iblk) = -1
                 wave_search_j(i,j,iblk) = -1
                 ice_search_i(i,j,iblk) = -1
                 ice_search_j(i,j,iblk) = -1
         end if

         found(1) = 0
         found(2) = 0

         ! For each ice-covered cell in block
         if ((aice(i,j,iblk).gt.p01).or.(frzmlt(i,j,iblk).gt.puny)) then

           ! get global indices from local indices
           i_g = this_block%i_glob(i)
           j_g = this_block%j_glob(j)

               tmp_mask = ((tlon_g.ge.(tlon(i,j,iblk)*rad_to_deg)-1.5_dbl_kind).and.(tlon_g.lt.(tlon(i,j,iblk)*rad_to_deg)+1.5_dbl_kind))

               ! if in NH
               if (tlat(i,j,iblk).gt.c0) then
                       ! find max latitude for this longitude that is not ice-covered and is less than current latitude
                       found = MAXLOC(tlat_g, MASK = tmp_mask.and.(nfl_g.lt.c0).and.(tlat_g.lt.tlat(i,j,iblk)*rad_to_deg) )

                       if ((SUM(found).gt.1).and.(nfl_g(found(1),found(2)).gt.-900)) then
                               nearest_wave_hs(i,j,iblk) = wave_hs_g(found(1),found(2))
                               nearest_wave_tz(i,j,iblk) = wave_tz_g(found(1),found(2))

                               cml_nfloes(i,j,iblk) = SUM(nfl_g, MASK = ((tmp_mask).and.(nfl_g.gt.c0).and. & 
                               (tlat_g.lt.tlat(i,j,iblk)*rad_to_deg).and.(tlat_g.gt.tlat_g(found(1),found(2)))) ) 
                       end if

              ! if in SH
               else
                       ! find min latitude for this longitude that is not ice-covered and is greater than current latitude
                       found = MINLOC(tlat_g, MASK = tmp_mask.and.(nfl_g.lt.c0).and.(tlat_g.gt.tlat(i,j,iblk)*rad_to_deg) ) 

                       if ((SUM(found).gt.1).and.(nfl_g(found(1),found(2)).gt.-900)) then
                               nearest_wave_hs(i,j,iblk) = wave_hs_g(found(1),found(2))
                               nearest_wave_tz(i,j,iblk) = wave_tz_g(found(1),found(2))

                               cml_nfloes(i,j,iblk) = SUM(nfl_g, MASK = ((tmp_mask).and.(nfl_g.gt.c0).and. & 
                                (tlat_g.gt.tlat(i,j,iblk)*rad_to_deg).and.(tlat_g.lt.tlat_g(found(1),found(2)))) ) 
                       end if

               end if
               if ((write_diag_wave).and.(SUM(found).gt.1)) then
                        if (nfl_g(found(1),found(2)).gt.-900) then
                                ice_search_i(i,j,iblk) = i_g
                                ice_search_j(i,j,iblk) = j_g
                                wave_search_i(i,j,iblk) = found(1)
                                wave_search_j(i,j,iblk) = found(2)
                        end if
               end if

              if (nearest_wave_tz(i,j,iblk).gt.puny) then
                   ! Bretscheider spectrum 
                   wave_spectrum(i,j,:,iblk) = (1/(c4*gravit)) * (nearest_wave_hs(i,j,iblk)**c2/nearest_wave_tz(i,j,iblk)**c4) &
                                           * T**c2* exp((-c1/pi)*(T/nearest_wave_tz(i,j,iblk))**c4)

                   ! attenuation coeff - quadratic fit from Horvat & Tziperman (2015)
                   if (aice(i,j,iblk).gt.p01) then ! attenuate!
                       alpha = exp(-.3203 + 2.058*(vice(i,j,iblk)/aice(i,j,iblk)) - .9375*T - &
                       .4269*(vice(i,j,iblk)/aice(i,j,iblk))**2 + 0.1566*(vice(i,j,iblk)/aice(i,j,iblk))*T + 0.0006*T**2)
                       wave_spectrum(i,j,:,iblk) = EXP(-alpha*cml_nfloes(i,j,iblk)) * wave_spectrum(i,j,:,iblk)
                   end if 

                   ! convert to function of freq
                   wave_spectrum(i,j,:,iblk) = wave_spectrum(i,j,:,iblk) / (T**c2)
              end if ! TZ

         end if ! aice
        end do
        end do
        end do
!$OMP END PARALLEL DO

    end if

               end subroutine find_wave

!=======================================================================
! Author: Lettie Roach, NIWA, 2016
!
! Call routine to fracture floes
!

        subroutine wave_break


      use ice_domain, only: nblocks
      use ice_flux, only: cml_nfloes, nearest_wave_hs, nearest_wave_tz
      use ice_state, only: aice, aicen, vice, vicen, trcrn


      integer (kind=int_kind) :: &
         iblk     ! block index


!$OMP PARALLEL DO PRIVATE(iblk)
      do iblk = 1, nblocks

         call wave_break_blocks (iblk,     &
                              cml_nfloes(:,:,iblk),&
                              nearest_wave_hs(:,:,iblk), &  
                              nearest_wave_tz(:,:,iblk), &
                              aice(:,:,iblk), aicen(:,:,:,iblk), &
                              vice(:,:,iblk), vicen(:,:,:,iblk), &
                              trcrn(:,:,:,:,iblk))

      end do !iblk 
!$OMP END PARALLEL DO
 

      end subroutine wave_break

!=======================================================================
! Author: Lettie Roach, NIWA, 2016
!
! Given the large arrays for omega and the fsd formed by wavebreaking, for each
! gridcell calculate the indices of the various discrete properties, use these
! to select the indices of the large arrays for omega and the fsd formed, and compute
! the change to the FSD. 

        subroutine wave_break_blocks (iblk,     &
                              cml_nfloes,&
                              nearest_wave_hs, &
                              nearest_wave_tz, &
                              aice, aicen, &
                              vice, vicen, &
                              trcrn)

       use ice_calendar, only: dt, istep
       use ice_domain, only: blocks_ice 
       use ice_flux, only: wave_tau, fracture_histogram                      
       use ice_itd, only: hin_max_init
       use ice_fsd, only: floe_rad_l, floe_rad_h, floe_rad_c, floe_binwidth, &
                          floe_area_c, &
                          d_afsd_wave, d_amfstd_wave
       use ice_state, only: nt_fsd
       use ice_blocks, only: nx_block, ny_block, get_block, block
       use ice_communicate, only: my_task


      integer (kind=int_kind), intent(in) :: &
         iblk ! current block

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         cml_nfloes, aice, vice, &        ! volume of ice
         nearest_wave_hs, nearest_wave_tz

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat), &
         intent(in) :: &
         aicen   , & ! concentration of ice
         vicen	     ! volume per unit area of ice          (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_ntrcr,ncat), &
         intent(inout) :: &
         trcrn     ! tracer array

      ! local

      type (block) :: &
         this_block    ! block information for current block

      integer (kind=int_kind) :: &
         nsubt,         & ! number of sub cycles required for 
                          ! timestep stability
         t,             & ! time index
         i, j,ii, jj    , & ! horizontal indices
         ilo, ihi, jlo, jhi, &
         iitmp, jjtmp, &
         ct,g,           &
         nh, nt,&
         nc, nd, nn, nk, &
         n, jn, ik   , & ! thickness category index
         k, ks       , & ! floe size index
         ij	         ! horizontal index, combines i and j loops
      
       real (kind=dbl_kind)  :: &
         dt_sub,       & ! length of subtimestep
         thckness,rdius,dmn,distance_tmp,totfrac, intgrl
 
       real (kind=dbl_kind), dimension (nfsd) :: &
         amfstd_tmp, &
         stability_test, &   ! check timestep stability
         omega, &
         gain, loss, &
         tempfracs

       real (kind=dbl_kind), dimension (nfsd, nfsd) :: &
        frac    ! to test

       real (kind=dbl_kind), dimension (ncat) :: &
         dh

       real (kind=dbl_kind), dimension (nx_block,ny_block,nfsd,ncat) :: &
           afstd, &           ! joint floe size and ice thickness area distribution 
           amfstd_init, &     ! tracer array
           afstd_init         ! aicen*trcrn
       


        d_amfstd_wave(:,:,:,:,iblk) = c0
        d_afsd_wave(:,:,:,iblk) = c0

        afstd=c0

        dh=c0

        do n=1,ncat
                dh(n) = hin_max_init(n) - hin_max_init(n-1)
        enddo

      !-----------------------------------------------------------------
      ! Loop over ice-cells for which we found waves.
      ! In each cell, find the area and hence number
      ! FSTD. Also find the mean ice
      ! thickness (for attenuation coefficient calculation)
      !---------------------------------------------------------------
         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
         do i = ilo, ihi
 
         if ((nearest_wave_hs(i,j).gt.swh_minval).and.(aice(i,j).gt.p01)) then

            nd=0
            nn=0
            nh=0
            nt=0
             
            thckness=vice(i,j)/aice(i,j)
            do n=1,ncat
                if ((thckness.le.hin_max_init(n)).and.(thckness.gt.hin_max_init(n-1))) nn = n            
            enddo
            if (thckness.gt.hin_max_init(ncat)) nn=ncat
                            
            do n=1,nwcat
                if ((cml_nfloes(i,j).gt.nfloe_l(n)).and.(cml_nfloes(i,j).le.nfloe_h(n))) nd = n
                if ((nearest_wave_hs(i,j).gt.hs_l(n)).and.(nearest_wave_hs(i,j).le.hs_h(n))) nh=n
                if ((nearest_wave_tz(i,j).gt.tz_l(n)).and.(nearest_wave_tz(i,j).le.tz_h(n))) nt=n
            end do
            if (cml_nfloes(i,j).gt.nfloe_h(nwcat)) nd=nwcat          
            if (nearest_wave_hs(i,j).gt.hs_h(nwcat)) nh=nwcat          
            if (nearest_wave_tz(i,j).gt.tz_h(nwcat)) nt=nwcat 
            if (nearest_wave_tz(i,j).lt.puny) nt = 1         

            ! check assignment is correct            
            if (.not.((nd.gt.0).and.(nh.gt.0).and.(nt.gt.0).and.(nn.gt.c0))) then
                print *, 'indices=',nn, nd, nh, nt
                print *, 'nearest tz ',nearest_wave_tz(i,j)
                print *, 'nearest hs ',nearest_wave_hs(i,j)
                print *, 'cml nfloes ',cml_nfloes(i,j)
                print *, tz_h
                print *, tz_l
                stop
            end if


            if (.not. ALL(fracture_histogram(nh,nt,nd,nn,:).eq.c0)) then


               ! we expect some wave fracture to occur!!
               do n = 1, ncat
                if ((aicen(i,j,n).gt.puny).and.(SUM(trcrn(i,j,nt_fsd+1:nt_fsd+nfsd-1,n)).gt.puny)) then
                        amfstd_init(i,j,:,n)=trcrn(i,j,nt_fsd:nt_fsd+nfsd-1,n)

                        if (ABS(SUM(amfstd_init(i,j,:,n))-c1).gt.puny) stop &
                                'init mFSTD not norm, wave'
                        
                        ! protect against small numerical errors
                        WHERE (amfstd_init(i,j,:,n).lt.puny) amfstd_init(i,j,:,n) = c0
                        amfstd_init(i,j,:,n) = amfstd_init(i,j,:,n) / SUM(amfstd_init(i,j,:,n))

                        amfstd_tmp =  amfstd_init(i,j,:,n)

                        ! frac does not vary within subcycle
                        frac(:,:) = c0
                        do k = 2, nfsd
                                frac(k,:k-1) = fracture_histogram(nh,nt,nd,nn,:k-1)
                                frac(k,k:) = c0 
                        end do
                        
                        do ks=1,nfsd
                         if (SUM(frac(ks,:)).gt.c0) frac(ks,:) = frac(ks,:)/SUM(frac(ks,:))
                        end do

                        ! now subcycle
                        do t = 1, NINT(dt/3600.0_dbl_kind*c5)

                                ! omega changes each subcycle
                                do k = 1,nfsd
                                        omega(k) = amfstd_tmp(k)*SUM(fracture_histogram(nh,nt,nd,nn,1:k-1))/ (D/c2) &
                                                /wave_tau(nh,nt,nd,nn) 
                                end do

                                if (SUM(omega).gt.c1+puny) stop &
                                 'omega cannot be greater than 1'
                        
                                loss = omega

                                do k =1,nfsd
                                        gain(k) = SUM(omega*frac(:,k)) 
                                end do

                                if (gain(nfsd).gt.puny) stop 'largest cat cannot gain, waves'
                                if (loss(1).gt.puny) stop 'smallest cat cannot lose, waves'

                                d_amfstd_wave(i,j,:,n,iblk) = gain(:) - loss(:)
                                
                                if (SUM(d_amfstd_wave(i,j,:,n,iblk)).gt.puny) stop &
                                        'area not cons, waves'

                                ! update
                                amfstd_tmp = amfstd_tmp + (dt/c5) * d_amfstd_wave(i,j,:,n,iblk)

                                if (ANY(amfstd_tmp.lt.-puny)) stop &
                                        'timestep still too long'

                                WHERE(amfstd_tmp.lt.c0) amfstd_tmp = c0
                        end do

                        ! new value for mFSTD, renorm in case of small
                        ! numerical errors
                        trcrn(i,j,nt_fsd:nt_fsd+nfsd-1,n) = amfstd_tmp/SUM(amfstd_tmp)

                        if (ANY(trcrn(i,j,nt_fsd:nt_fsd+nfsd-1,n).lt.c0)) stop 'neg wb'

                        ! for diagnostics
                        d_amfstd_wave(i,j,:,n,iblk) = trcrn(i,j,nt_fsd:nt_fsd+nfsd-1,n) - amfstd_init(i,j,:,n)  
                        d_afsd_wave(i,j,:,iblk) = aicen(i,j,n)*(trcrn(i,j,nt_fsd:nt_fsd+nfsd-1,n) - amfstd_init(i,j,:,n))

                end if ! aicen>puny

               end do ! n

           end if  ! omega zero
   
          end if ! nearest wave > 0.01
          end do !i
          end do !j     
        

              end subroutine wave_break_blocks


!=======================================================================
! Author: Lettie Roach, NIWA, 2016
!
! Based on MatLab code from Horvat & Tziperman (2015). Calculates functions 
! to describe the change in the FSD when waves fracture ice. Waves generated
! from the Bretschneider spectrum with (at the moment) constant Hs and Tz. We 
! calculate extrema and if these are successive maximum, minimum, maximum or
! vice versa, and have strain greater than a critical strain, break ice and 
! create new floes with lengths equal to these distances
     
        subroutine wave_change( hbar, nfloes_cml,    &
                                Hs, Tz, &
                                frac, &
                                tau)
        
        use ice_constants,only: pi
        use ice_fsd, only: floe_rad_l, floe_binwidth, floe_rad_c, floe_rad_h
        use ice_domain_size, only: nfsd
        use ice_communicate, only: my_task, master_task

     real (kind=dbl_kind), intent(in), dimension (nwcat) :: &
        nfloes_cml,     & ! distance to nearest ocean cell (m)
        Hs ,  &           ! significant wave height in nearest ocean cell (m)
        Tz                ! zero crossing period in nearest ocean cell (s)

     real (kind=dbl_kind), intent(in), dimension (ncat) :: &
        hbar             ! mean ice thickness in this gridcell (m)
 
     real (kind=dbl_kind), dimension (nwcat,nwcat,nwcat,ncat,nfsd), &
         intent(out) :: &
         frac

      real (kind=dbl_kind), dimension (nwcat,nwcat,nwcat,ncat), intent(out) :: &
        tau

      !------local variables---------------------------------
 
        real (kind=dbl_kind) :: &
                lambdamin, dlambda

        integer (kind=int_kind) :: &
                spcing, first, last, maxct, minct, lctn, tmp, i, j, k, n, m, &
                nn, nf, nd, nt,nh

        real (kind=dbl_kind), dimension(:), allocatable :: &
               fraclengths

        real (kind=dbl_kind), dimension(nx) :: &
               X, eta

        real (kind=dbl_kind), dimension(nfsd) :: &
                frachistogram 

        real (kind=dbl_kind), dimension(Nl) :: &
                lambda,T, S_B, rand_array, phi, a_i, logalpha, alpha, &
                v_group, summand

        logical (kind=log_kind) :: &
                e_stop          ! if true, stop and return zero omega and fsdformed

     


     lambdamin=c2*floe_rad_c(1)
     dlambda = floe_rad_c(1)
     ! wavelengths span
     do j=1,Nl
        lambda(j)=lambdamin+j*dlambda
     end do

     ! periods from wavelengths	
     T=(c2*pi*lambda/gravit)**p5

     do j=1,nx
        X(j)= j*dx
     end do


     if (dx.gt.floe_binwidth(1)/c2 ) stop 'must reduce dx'

     do nh=1,nwcat  
     do nt=1,nwcat
     do nd=1,nwcat
     do nn=1,ncat

        frac(nh,nt,nd,nn,:)= c0
        tau(nh,nt,nd,nn) = c0

        if (my_task.eq.master_task) print *, 'now ',nh,nt,nd,nn

      !-----------------------------------------------------------------
      ! Setting up
      !-----------------------------------------------------------------
        e_stop=.false.
      
       
        ! Bretscheider spectrum
        S_B=(1/(c4*gravit)) * (Hs(nh)**c2/Tz(nt)**c4) * T**c2* exp((-c1/pi)*(T/Tz(nt))**c4)

        ! quadratic fit from Horvat & Tziperman (2015)
        logalpha = -.3203 + 2.058*hbar(nn) - .9375*T - .4269*hbar(nn)**2 + 0.1566*hbar(nn)*T + 0.0006*T**2
        alpha = exp(logalpha)

        ! spectral coefficients
        a_i=sqrt(c2*S_B*dlambda)


        allocate(fraclengths(1))
        fraclengths(1)=c0
        
        !----------------loop 'loopcts' times starts here
        do i=1,loopcts

                ! random phase for each Fourier component
                ! varies in each j loop
                call RANDOM_NUMBER(rand_array)
                phi = c2*pi*rand_array
                                
                do j=1,nx
                        ! exponential damping, ignore attenuation in
                        ! current grid cell
                        summand = a_i*EXP(-alpha*nfloes_cml(nd))*COS(2*pi*X(j)/lambda+phi)
                        eta(j)=SUM(summand)

                end do

                call get_fraclengths(X, eta, fraclengths, hbar(nn), e_stop)

       end do !i

        !---------------------end loop over loopcts
        ! bin into FS categories
        ! histogram values get bigger during loop

        if (SIZE(fraclengths).eq.1) e_stop = .true.
        if (ALL(fraclengths(:).eq.c0)) e_stop = .true.

        if (.not. e_stop) then

                ! convert from diameter to radii
                fraclengths(:)=fraclengths(:)/c2

                frachistogram(:)=c0

                ! highest cat cannot be fractured into
                do j=1,size(fraclengths)
                        do k=1,nfsd-1
                                if ((fraclengths(j).ge.floe_rad_l(k)).and.(fraclengths(j).lt.floe_rad_l(k+1))) then
                                        frachistogram(k)=frachistogram(k)+1
                                end if
                        end do

                        if (fraclengths(j).lt.(floe_rad_l(1)-puny)) stop 'fractures too small'
                
                end do
        end if
 
       if (SUM(frachistogram).eq.c0) e_stop = .true.

       if (.not. e_stop) then

                do k=1,nfsd
                        frac(nh,nt,nd,nn,k)=floe_rad_c(k)*frachistogram(k)
                end do
               
                frac(nh,nt,nd,nn,:) = (D/c2) * frac(nh,nt,nd,nn,:) / SUM(frac(nh,nt,nd,nn,:))


                ! group velocity of waves
                v_group=0.5*(gravit*lambda/(2*pi))**0.5

                ! time to transit the domain
                tau(nh,nt,nd,nn)=(16/Hs(nh)**2)*dlambda*D*SUM(S_B/v_group)

        end if! e_stop


        deallocate(fraclengths)

        if (my_task.eq.master_task) then
        print *, 'frac equals ',(frac(nh,nt,nd,nn,k),k=1,nfsd)
        print *, 'tau equals ',tau(nh,nt,nd,nn)
        end if
        
        end do
        end do
        end do
        end do

        end subroutine wave_change




!===========================================================================
!  Given the (attenuated) sea surface height, find the strain across triplets
!  of max, min, max or min, max, min (local extrema within 10m).
!  If this strain is greater than the  critical strain, ice can fracture
!  and new floes are formed with sizes equal to the distances between
!  extrema
!
        subroutine get_fraclengths(X, eta, fraclengths, hbar, e_stop)

        use ice_communicate, only: my_task

        real (kind=dbl_kind) :: &
                hbar

        real (kind=dbl_kind), intent(in), dimension (nx) :: &
                X, eta
 
        real (kind=dbl_kind), intent(inout), dimension (:), allocatable :: &
                fraclengths
 
        logical (kind=log_kind), intent(inout) :: &
                e_stop          ! if true, stop and return zero omega and fsdformed


        ! local
        integer (kind=int_kind) :: &
                spcing, &       ! distance in dx over which to search for extrema on each side of point
                j, k, &         ! indices to iterate over domain
                first, last, &  ! indices over which to search for extrema
                j_neg, &        ! nearest extrema backwards
                j_pos, &        ! nearest extrema forwards
                n_above         ! number of points where strain is above
                                ! critical strain
        

         real (kind=dbl_kind), dimension(nx) :: &
                strain, &        ! the strain between triplets of extrema
                frac_size_one, & !
                frac_size_two

        logical (kind=log_kind), dimension(nx) :: &
                is_max, is_min, &       ! arrays to hold whether each point is a local max or min
                is_extremum, &          ! or extremum
                is_triplet              ! or triplet of extrema

         real (kind=dbl_kind) :: &
                delta, &        ! difference in x between current and prev extrema
                delta_pos       ! difference in x between next and current extrema

       integer (kind=int_kind), dimension(1) :: &
        maxj, minj  ! indices of local max and min


       ! -------equivalent of peakfinder2
       !given eta and spcing, compute extremelocs in ascending order
       spcing=nint(threshold/dx)

       is_max = .false.
       is_min = .false.
       is_extremum = .false.
       is_triplet = .false.
       strain = c0
       frac_size_one = c0
       frac_size_two = c0
       j_neg = 0
       j_pos = 0      
 
       ! search for local max and min within spacing of
       ! 10m on either side of each point

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
                         print *, 'X ',X
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
                                delta_pos = X(j_pos) - X(j)
                                delta = X(j) - X(j_neg)
                                strain(j) = (hbar/c2) *(eta(j_neg)*delta_pos - &
                                        eta(j)*(delta_pos+delta) + eta(j)*delta ) &
                                        /(delta*delta_pos*(delta+delta_pos))


                                if (strain(j).gt.straincrit) then
                                        frac_size_one(j) = X(j_pos) - X(j)
                                        frac_size_two(j) = X(j) - X(j_neg)
                                end if
                        end if

                end if

        end do

        n_above = COUNT(strain.gt.straincrit)
        if (n_above.gt.0) then
                DEALLOCATE(fraclengths)
                ALLOCATE(fraclengths(2*n_above))
                fraclengths(1:n_above) = PACK(frac_size_one,(frac_size_one.gt.c0))
                fraclengths(n_above+1:2*n_above) = &
                        PACK(frac_size_two,(frac_size_two.gt.c0))

                e_stop = .false.
        else
                e_stop = .true.

        end if

        end subroutine get_fraclengths

!=======================================================================

     
   end module ice_wavebreaking

!=======================================================================



