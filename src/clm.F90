module clm_module
  use clm_type_module
  use tile_type_module
  use grid_type_module
  use clm1d_type_module
  
  implicit none

  private

  public:: clm_setup_unstructured, &
       clm_setup_structured, &
       clm_advance_unstructured, &
       clm_advance_structured

contains

  !
  ! Set up an unstructured, general use case of CLM
  !
  ! Assumes the drv data has been set already
  !------------------------------------------------------
  subroutine clm_setup_unstructured(clm)
    type(clm_type) :: clm
  end subroutine clm_setup_unstructured


  !
  ! Set up a structured, general use case of CLM
  !
  ! Assumes the basic data has been set already
  !------------------------------------------------------
  subroutine clm_setup_structured(clm, nx, ny, nz, ntiles, ntypes, nsteps, istep_host,
    nx_f, ny_f, nz_f, topo, dz_mult, soi_z)
    implicit none
    type(clm_type) clm
    integer,intent(in) :: soi_z ! reference temperature for vegetation
    integer,intent(in) :: istep_host, nsteps
    integer,intent(in) :: nx,ny,nz, ntiles, ntypes
    integer,intent(in) :: nx_f,ny_f,nz_f

    real(r8) :: topo((nx+2)*(ny+2)*(nz+2))         ! mask from host code 0 for inactive, 1 for active, on grid w/ ghost nodes for current proc
    real(r8) :: dz_mult((nx+2)*(ny+2)*(nz+2))   ! dz multiplier from ParFlow on PF grid w/ ghost nodes for current proc
    
    ! locals
    integer t,i,j,k,l
    integer j_incr, k_incr, counter
    real(r8) :: total
    integer ierr
    ierr = 0

    ! set up structure
    clm%grid_ncols = nx
    clm%drv%nr = nx
    clm%grid_ncols = ny
    clm%drv%nc = ny
    clm%ntiles = ntiles
    clm%drv%nt = ntypes    ! currently hard-coded, but could change

    clm%clm%soi_z = soi_z

    ! increments for getting around the structured arrays
    j_incr = nx_f
    k_incr = nx_f*ny_f

    
    call clm_setup_common(clm)

    !=== Initialize the CLM topography mask 
    !    This is two components: 
    !    1) a x-y mask of 0 o 1 for active inactive and 
    !    2) a z/k mask that takes three values 
    !      (1)= top of land surface/host code's domain 
    !      (2)= top-nlevsoi and 
    !      (3)= the bottom of the host code's domain.
    if (clm%io%log /= 0) write(clm%io%log,*) "Initialize the CLM topography mask"
    call ASSERT(nz >= nlevsoi, "Number of vertical cells must be at least 10 (nlevsoi).")
    
    do t=1,clm%drv%nch
       i=clm%tile(t)%col
       j=clm%tile(t)%row
       counter = 0
       clm%clm(t)%topo_mask(3) = 1

       do k = nz, 1, -1 ! loop over host nz
          l = 1+i + (nx+2)*(j) + (nx+2)*(ny+2)*(k)
          if (topo(l) > 0) then
             counter = counter + 1
             if (counter == 1) then 
                clm(t)%topo_mask(1) = k
                clm(t)%planar_mask = 1
             end if
          endif
          if (topo(l) == 0 .and. topo(l+k_incr) > 0) clm%clm(t)%topo_mask(3) = k+1
       enddo ! k
       clm%clm(t)%topo_mask(2) = clm%clm(t)%topo_mask(1)-nlevsoi
    enddo ! t

    !=== Initialize the variable dz over the root column
    !    -- Copy dz multipliers for root zone cells from PF grid to 1D array
    !    -- Then loop to recompute clm(t)%z(j), clm(t)%dz(j), clm(t)%zi(j) 
    !       (replaces values set in drv_clmini)
     do t = 1,clm%drv%nch
        i = clm%tile(t)%col
        j = clm%tile(t)%row

        ! BH: modification of the interfaces depths and layers thicknesses to match 
        clm%clm(t)%zi(0) = 0.   

        ! check if cell is active
        if (clm%clm(t)%planar_mask == 1) then
           ! reset node depths (clm%z) based on variable dz multiplier
           do k = 1, nlevsoi
              l                 = 1+i + j_incr*(j) + k_incr*(clm%clm(t)%topo_mask(1)-(k-1))
	      clm(t)%dz(k)	= clm%drv%dz * dz_mult(l) ! basile
              if (k==1) then
                 clm(t)%z(k)    = 0.5 * clm%drv%dz * dz_mult(l)
                 clm(t)%zi(k)	= clm%drv%dz * dz_mult(l) ! basile
              else
                 total          = 0.0
                 do k1 = 1, k-1
                    l1          = 1+i + j_incr*(j) + k_incr*(clm(t)%topo_mask(1)-(k1-1))
                    total       = total + (clm%drv%dz * dz_mult(l1))
                 enddo
                 clm%z(k)       = total + (0.5 * clm%drv%dz * dz_mult(l))
                 clm%zi(k)	= total + clm%drv%dz * dz_mult(l)! basile
              endif
           enddo

           !! Overwrite Rootfr disttribution: start
           !! the following overwrites the root fraction definition which is previously set up in drv_clmini.F90 
           !! but based on constant DZ, regardless of dz_mult.  Fix to use the dz_mult.
           do k = 1, nlevsoi-1
              clm%clm(t)%rootfr(k) = .5*( exp(-clm%tile(t)%roota*clm%clm(t)%zi(k-1))  &
                   + exp(-clm%tile(t)%rootb*clm%clm(t)%zi(k-1))  &
                   - exp(-clm%tile(t)%roota*clm%clm(t)%zi(k  ))  &
                   - exp(-clm%tile(t)%rootb*clm%clm(t)%zi(k  )) )
           enddo
           clm%clm(t)%rootfr(nlevsoi)=.5*( exp(-clm%tile(t)%roota*clm%clm(t)%zi(nlevsoi-1))&
                + exp(-clm%tile(t)%rootb*clm%clm(t)%zi(nlevsoi-1)))

           ! reset depth variables assigned by user in clmin file 
           do k = 1,nlevsoi
              if (clm%grid(clm%tile(t)%col,clm%tile(t)%row)%rootfr /= clm%drv%udef) &
                   clm%clm(t)%rootfr(bl) = clm%grid(clm%tile(t)%col,clm%tile(t)%row)%rootfr
           enddo
        endif ! active/inactive
     enddo !t 


     !=== Read restart file or set initial conditions
     call drv_restart(1,drv,tile,clm,rank,istep_pf)        ! (1=read,2=write)
     
     !=== Loop over CLM tile space to set keys/constants from PF
     !    (watsat, residual sat, irrigation keys)
     do t=1,clm%drv%nch  

        ! check if cell is active
        if (clm%clm(t)%planar_mask == 1) then

           ! for beta and veg stress formulations
           clm%clm(t)%beta_type          = beta_typepf
           clm%clm(t)%vegwaterstresstype = veg_water_stress_typepf
           clm%clm(t)%wilting_point      = wilting_pointpf
           clm%clm(t)%field_capacity     = field_capacitypf
           clm%clm(t)%res_sat            = res_satpf

           ! for irrigation
           clm%clm(t)%irr_type           = irr_typepf
           clm%clm(t)%irr_cycle          = irr_cyclepf
           clm%clm(t)%irr_rate           = irr_ratepf
           clm%clm(t)%irr_start          = irr_startpf
           clm%clm(t)%irr_stop           = irr_stoppf
           clm%clm(t)%irr_threshold      = irr_thresholdpf     
           clm%clm(t)%threshold_type     = irr_thresholdtypepf

           ! set clm watsat, tksatu from host porosity
           ! convert t to i,j index
           i=clm%tile(t)%col        
           j=clm%tile(t)%row
           do k = 1, nlevsoi ! loop over clm soil layers (1->nlevsoi)
              ! convert clm space to parflow space, note that PF space has ghost nodes
              l = 1+i + j_incr*(j) + k_incr*(clm(t)%topo_mask(1)-(k-1))
              clm%clm(t)%watsat(k)       = porosity(l)
              clm%clm(t)%tksatu(k)       = clm%clm(t)%tkmg(k)*0.57**clm%clm(t)%watsat(k)
           end do !k
        endif ! active/inactive
     end do !t

     !=== Read restart file or set initial conditions
     call drv_restart(1,clm%drv,clm%tile,clm%clm,clm%rank,istep_host)        ! (1=read,2=write)

  end subroutine clm_setup_structured

  !
  ! Set up a structured, general use case of CLM
  !
  ! Assumes the basic data has been set already
  !------------------------------------------------------
  subroutine clm_setup_unstructured(clm, ncells, ntiles, ntypes, nsteps, istep_host, soi_z)
    implicit none
    type(clm_type) clm
    integer,intent(in) :: soi_z ! reference temperature for vegetation
    integer,intent(in) :: istep_host
    integer,intent(in) :: nrows, ncols, ntiles, ntypes
    
    ! locals
    integer t,i,j,k
    integer ierr
    ierr = 0

    ! set up structure
    clm%grid_nrows = ncells
    clm%drv%nr = ncells
    clm%grid_ncols = 1
    clm%drv%nc = 1
    clm%ntiles = ntiles
    clm%drv%nt = ntypes    ! currently hard-coded, but could change

    clm%clm%soi_z = soi_z
    
    call clm_setup_common(clm)

  end subroutine clm_setup_unstructured
  

  
  !
  ! Set up a structured, general use case of CLM
  !
  ! Assumes the basic data has been set already
  !------------------------------------------------------
  subroutine clm_setup_common(clm)
    use clm1d_varpar, only : nlevsoi
    implicit none
    type(clm_type) clm
    if (clm%io%log /= 0) &
         write(clm%io%log,*) &
         "CLM (IDEAS version) Log file\n----------------------------------------\nUnstructured Setup"

    if (clm%io%log /= 0) &
         write(clm%io%log,*) &
         "  n_grid_rows = ", clm%grid_nrows, " n_grid_cols = ", clm%grid_ncols, " n_tiles = ", clm%ntiles

    clm%drv%nr = clm%grid_nrows
    clm%drv%nc = clm%grid_ncols
    clm%grid => grid_create(clm%grid_nrows, clm%grid_ncols, clm%drv%nt)

     !=== Read in the clm input (drv_clmin.dat)
     call drv_readclmin(clm%drv, clm%grid, clm%rank, clm%io%log)

     if (clm%io%log /= 0) then
        write(clm%io%log,*) "CLM startcode for date (1=restart, 2=defined):", clm%drv%startcode
        write(clm%io%log,*) "CLM IC (1=restart, 2=defined):", clm%drv%clm_ic

        !=== @RMM check for error in IC or starting time
        if (clm%drv%startcode == 0) stop
        if (clm%drv%clm_ic == 0) stop
     end if

     !=== Allocate memory for subgrid tile space
     !=== CURRENT =============================================================================================
     !=== Because we only use one tile per grid cell, we don't need to call readvegtf to determine actual nch
     !    (nch is just equal to number of cells (nr*nc))
     clm%drv%nch = clm%grid_nrows * clm%grid_ncols
     if (clm%io%log /= 0) write(clm%io%log,*) "Allocate arrays -- using num_tiles =", clm%drv%nch
     clm%tile => tile_create_n(clm%drv%nch)
     clm%clm => clm1d_create_n(clm%drv%nch)

     !=== Set clm diagnostic indices and allocate space
     clm%clm%surfind = clm%drv%surfind 
     clm%clm%soilind = clm%drv%soilind
     clm%clm%snowind = clm%drv%snowind

     do t=1,clm%drv%nch 
        allocate(clm%clm(t)%diagsurf(1:clm%drv%surfind),stat=ierr); call drv_astp(ierr) 
        allocate(clm%clm(t)%diagsoil(1:clm%drv%soilind,1:nlevsoi),stat=ierr); call drv_astp(ierr)
        allocate(clm%clm(t)%diagsnow(1:clm%drv%snowind,-nlevsno+1:0),stat=ierr); call drv_astp(ierr)
     end do

     !=== NBE: Define the reference layer for the seasonal soil
     if (clm%io%log /= 0) write(clm%io%log,*) "Check soi_z",clm%clm%soi_z

     !=== Initialize clm derived type components
     if (clm%io%log /= 0) write(clm%io%log,*) "Calling clm_typini"
     call clm_typini(clm%drv%nch,clm%clm,istep_host)

     !=== Read in vegetation data and set tile information accordingly
     if (clm%io%log /= 0) write(clm%io%log,*) "Read in vegetation data and set tile information accordingly"
     call drv_readvegtf (drv, grid, tile, clm, nx, ny, ix, iy, gnx, gny, rank)

     !=== Transfer grid variables to tile space 
     if (clm%io%log /= 0) write(clm%io%log,*) "Transfer grid variables to tile space ", clm%drv%nch
     do t = 1, clm%drv%nch
        call drv_g2clm(clm%drv%udef, clm%drv, clm%grid, clm%tile(t), clm%clm(t))   
     enddo

     !=== Read vegetation parameter data file for IGBP classification
     if (clm%io%log /= 0) write(clm%io%log,*) "Read vegetation parameter data file for IGBP classification"
     call drv_readvegpf(clm%drv, clm%grid, clm%tile, clm%clm)  

     !=== Initialize CLM and DIAG variables
     if (clm%io%log /= 0) write(clm%io%log,*) "Initialize CLM and DIAG variables"
     do t=1,clm%drv%nch 
        clm%clm%kpatch = t
        call drv_clmini(clm%drv, clm%grid, clm%tile(t), clm%clm(t), clm%istep_host) !Initialize CLM Variables
     enddo
   end subroutine clm_setup_common


    
   
  !
  ! Advances the timestep
  !------------------------------------------------------
  subroutine clm_advance_unstructured(clm)
    type(clm_type) clm
  end subroutine clm_advance_unstructured

  !
  ! Advances the timestep
  !------------------------------------------------------
  subroutine clm_advance_structured(clm, saturation, pressure,ix,iy, nx_f, ny_f, nz, ip)
    type(clm_type) clm

    ! local
    integer j_incr, k_incr
    
    ! increments for getting around the structured arrays
    j_incr = nx_f
    k_incr = nx_f*ny_f
    
    
    !=========================================================================
    !=== Time looping
    !=========================================================================
    
    !=== Call routine to copy PF variables to CLM space 
    !    (converts saturation to soil moisture)
    !    (converts pressure from m to mm)
    !    (converts soil moisture to mass of h2o)
    call pfreadout(clm%clm,clm%drv,clm%tile,saturation,pressure,clm%rank, &
         ix,iy,clm%grid_nx,clm%grid_ny,nz,j_incr,k_incr,ip)

    !=== Advance time (CLM calendar time keeping routine)
    clm%drv%endtime = 0
    call drv_tick(clm%drv)

    !RMM: writing a CLM.log.out file with basic information only from the master node (0 processor)
    !
    if (rank==0) then
       write(9919,*)
       write(9919,*) "CLM starting time =", time, "gmt =", drv%gmt,"istep_pf =",istep_pf 
       write(9919,*) "CLM day =", drv%da, "month =", drv%mo,"year =", drv%yr
    end if ! CLM log


    !=== Read in the atmospheric forcing for off-line run
    !    (values no longer read by drv_getforce, passed from PF)
    !    (drv_getforce is modified to convert arrays from PF input to CLM space)
    !call drv_getforce(drv,tile,clm,nx,ny,sw_pf,lw_pf,prcp_pf,tas_pf,u_pf,v_pf,patm_pf,qatm_pf,istep_pf)
    !BH: modification of drv_getforc to optionnaly force vegetation (LAI/SAI/Z0M/DISPLA): 
    !BH: this replaces values from clm_dynvegpar called previously from drv_clmini and 
    !BH: replaces values from drv_readvegpf
    call drv_getforce(drv,tile,clm,nx,ny,sw_pf,lw_pf,prcp_pf,tas_pf,u_pf,v_pf, &
         patm_pf,qatm_pf,lai_pf,sai_pf,z0m_pf,displa_pf,istep_pf,clm_forc_veg)
    !=== Actual time loop
    !    (loop over CLM tile space, call 1D CLM at each point)
    do t = 1, drv%nch     
       clm(t)%qflx_infl_old       = clm(t)%qflx_infl
       clm(t)%qflx_tran_veg_old   = clm(t)%qflx_tran_veg
       if (clm(t)%planar_mask == 1) then
          call clm_main (clm(t),drv%day,drv%gmt) 
       else
       endif ! Planar mask
    enddo ! End of the space vector loop

    !=== Write CLM Output (timeseries model results)
    if (clm_1d_out == 1) then 
       call drv_1dout (drv, tile,clm,clm_write_logs)
    endif


    !=== Call 2D output routine
    !     Only call for clm_dump_interval steps (not time units, integer units)
    !     Only call if write_CLM_binary is True
    if (mod(dble(istep_pf),clm_dump_interval)==0)  then
       if (write_CLM_binary==1) then

          ! Call subroutine to open (2D-) output files
          call open_files (clm,drv,rank,ix,iy,istep_pf,clm_output_dir,clm_output_dir_length,clm_bin_output_dir) 

          ! Call subroutine to write 2D output
          call drv_2dout  (drv,grid,clm)

          ! Call to subroutine to close (2D-) output files
          call close_files(clm,drv)

       end if ! write_CLM_binary
    end if ! mod of istep and dump_interval


    !=== Copy values from 2D CLM arrays to PF arrays for printing from PF (as Silo)
    do t=1,drv%nch
       i=tile(t)%col
       j=tile(t)%row
       l = 1+i + (nx+2)*(j) + (nx+2)*(ny+2) 
       if (clm(t)%planar_mask==1) then
          eflx_lh_pf(l)      = clm(t)%eflx_lh_tot
          eflx_lwrad_pf(l)   = clm(t)%eflx_lwrad_out
          eflx_sh_pf(l)      = clm(t)%eflx_sh_tot
          eflx_grnd_pf(l)    = clm(t)%eflx_soil_grnd
          qflx_tot_pf(l)     = clm(t)%qflx_evap_tot
          qflx_grnd_pf(l)    = clm(t)%qflx_evap_grnd
          qflx_soi_pf(l)     = clm(t)%qflx_evap_soi
          qflx_eveg_pf(l)    = clm(t)%qflx_evap_veg 
          qflx_tveg_pf(l)    = clm(t)%qflx_tran_veg
          qflx_in_pf(l)      = clm(t)%qflx_infl 
          swe_pf(l)          = clm(t)%h2osno 
          t_g_pf(l)          = clm(t)%t_grnd
          qirr_pf(l)         = clm(t)%qflx_qirr
          irr_flag_pf(l)     = clm(t)%irr_flag
       else
          eflx_lh_pf(l)      = -9999.0
          eflx_lwrad_pf(l)   = -9999.0
          eflx_sh_pf(l)      = -9999.0
          eflx_grnd_pf(l)    = -9999.0
          qflx_tot_pf(l)     = -9999.0
          qflx_grnd_pf(l)    = -9999.0
          qflx_soi_pf(l)     = -9999.0
          qflx_eveg_pf(l)    = -9999.0
          qflx_tveg_pf(l)    = -9999.0
          qflx_in_pf(l)      = -9999.0
          swe_pf(l)          = -9999.0
          t_g_pf(l)          = -9999.0
          qirr_pf(l)         = -9999.0
          irr_flag_pf(l)     = -9999.0
       endif
    enddo


    !=== Repeat for values from 3D CLM arrays
    do t=1,drv%nch            ! Loop over CLM tile space
       i=tile(t)%col
       j=tile(t)%row
       if (clm(t)%planar_mask==1) then
          do k = 1,nlevsoi       ! Loop from 1 -> number of soil layers (in CLM)
             l = 1+i + j_incr*(j) + k_incr*(nlevsoi-(k-1))
             t_soi_pf(l)     = clm(t)%t_soisno(k)
             qirr_inst_pf(l) = clm(t)%qflx_qirr_inst(k)
          enddo
       else
          do k = 1,nlevsoi
             l = 1+i + j_incr*(j) + k_incr*(nlevsoi-(k-1))
             t_soi_pf(l)     = -9999.0
             qirr_inst_pf(l) = -9999.0
          enddo
       endif
    enddo



    !=== Write Daily Restarts
    if (clm_write_logs==1) then
       write(999,*) "End of time advance:" 
       write(999,*) 'time =', time, 'gmt =', drv%gmt, 'endtime =', drv%endtime
    endif
    if (rank==0) then
       write(9919,*) "End of time advance:"
       write(9919,*) 'time =', time, 'gmt =', drv%gmt, 'endtime =', drv%endtime
    end if !! rank 0, write log info

    ! if ( (drv%gmt==0.0).or.(drv%endtime==1) ) call drv_restart(2,drv,tile,clm,rank,istep_pf)
    ! ----------------------------------
    ! NBE: Added more control over writing of the RST files
    if (clm_next == 1) then
       if (clm_last_rst==1) then
          d_stp=0
       else
          d_stp = istep_pf
       endif

       if (clm_daily_rst==1) then
          if ( (drv%gmt==0.0).or.(drv%endtime==1) ) then
             if (rank==0) write(9919,*) 'Writing restart time =', time, 'gmt =', drv%gmt, 'istep_pf =',istep_pf
             !! @RMM/LEC  add in a TCL file that sets an istep value to better automate restarts
             if (rank==0) then
                open(393, file="clm_restart.tcl",action="write")
                write(393,*) "set istep ",istep_pf
                close(393)
             end if  !  write istep corresponding to restart step

             call drv_restart(2,drv,tile,clm,rank,d_stp)
          end if
       else
          call drv_restart(2,drv,tile,clm,rank,d_stp)
       endif

    endif
    ! ---------------------------------

    !=== Call routine to calculate CLM flux passed to PF
    !    (i.e., routine that couples CLM and PF)
    call pf_couple(drv,clm,tile,evap_trans,saturation,pressure,porosity,nx,ny,nz,j_incr,k_incr,ip,d_stp)


    !=== LEGACY ===========================================================================================
    !    (no longer needed because current setup restricts one tile per grid cell) 
    !=== Return required surface fields to atmospheric model (return to grid space)
    !    (accumulates tile fluxes over grid space)
    ! call drv_clm2g (drv, grid, tile, clm)


    !=== Write spatially-averaged BC's and IC's to file for user
    if (clm_write_logs==1) then ! NBE
       if (istep_pf==1) call drv_pout(drv,tile,clm,rank)
    endif

    !=== If at end of simulation, close all files
    if (drv%endtime==1) then
       ! close(166)
       ! close(199)
       if (clm_write_logs==1) close(999)
       if (rank == 0) close (9919)
    end if

  end subroutine clm_advance_structured



end module clm_module



