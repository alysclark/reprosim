module abm_models

  implicit none
  !Module parameters

  !Module types

  !Module depvar

  !Interfaces
  private
  public initialise_abm
  public move_cells_force
  public get_current_t
  public get_current_dt
contains
!
!###################################################################################
!
subroutine initialise_abm(model_type,total_cells,num_forces,time_step,min_time_step)
    use arrays, only: dp, cell_type,cell_list, num_cells, num_cell_f,abm_control
    use diagnostics, only: enter_exit,get_diagnostics_level
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_EVALUATE_ABM" :: EVALUATE_ABM

    character(len=60), intent(in) :: model_type
    integer, intent(in) :: total_cells
    integer, intent(in) :: num_forces
    real(dp), intent(in) :: time_step
    real(dp), intent(in) :: min_time_step

    !local variables
    integer :: nseed,i,input_seed
    integer, allocatable :: seed(:)
    character(len=60) :: sub_name
    sub_name = 'initialise_abm'
    call enter_exit(sub_name,1)

    num_cells = total_cells
    num_cell_f = num_forces
    !Allocate memory
    call allocate_abm_memory

    !Define model control parameters
    abm_control%delta_t = time_step
    abm_control%delta_t_min = min_time_step
    !Intitialise time
    abm_control%current_time = 0.0_dp

    input_seed = 1 !This needs to be an input later

    !Initialise cell location
    if(model_type.eq.'trophoblast')then
      call place_cells_plug
    else
      print*, 'Warning: Cannot initialise cell location, model type not implemented'
    endif

    !Initialise lists of neighbouring cells
    call setup_nbrlists
    !Need to set up random number seeds
    call random_seed(size = nseed)
    allocate(seed(nseed))
    do i = 1,nseed
      seed(i) = input_seed+1-i
    enddo
    call random_seed(put = seed)



    !call deallocate_abm_memory
    call enter_exit(sub_name,2)
end subroutine initialise_abm
!
!###################################################################################
!

subroutine allocate_abm_memory
    use arrays, only: cell_type,cell_list, num_cells,num_cell_f,cell_field
    use diagnostics, only: enter_exit,get_diagnostics_level

    !local variables
    integer :: AllocateStatus
    character(len=60) :: sub_name

    sub_name = 'allocate_abm_memory'
    call enter_exit(sub_name,1)

    ! cells in the domain are stored in a list, while the
    ! occupancy array holds the indices of cells in the list.  When a cell
    ! leaves the domain or dies a gap is created in the list.
    ! The locations of such gaps are stored in the gaplist, the total number
    ! of such gaps is ngaps.  A cell entering the domain is allocated an index
    ! from the tail of this list, if ngaps > 0, or else it is added to the end of the cell list.

    allocate (cell_list(num_cells), STAT = AllocateStatus) !List of cells
    if (AllocateStatus /= 0) STOP "*** Not enough memory for cell_list array ***"
    allocate (cell_field(num_cells,num_cell_f,3), STAT = AllocateStatus) !List of cells
    if (AllocateStatus /= 0) STOP "*** Not enough memory for cell_list array ***"

    call enter_exit(sub_name,2)
end subroutine allocate_abm_memory

!
!###################################################################################
!
subroutine create_cell(kcell,rsite,ctype,gen,tag,region,dividing)
    use arrays, only: dp,cell_type,cell_list, cell_stat,plug_params,abm_control
    use diagnostics, only: enter_exit,get_diagnostics_level


    type(cell_type), pointer :: cp
    integer :: kcell, ctype, gen, tag, region
    real(dp) :: rsite(3)
    logical :: dividing
    real(dp) :: tnow
    integer :: kpar = 0
    character(len=60) :: sub_name
    sub_name = 'create_cell'
    call enter_exit(sub_name,1)

    tnow = plug_params%istep*abm_control%delta_t
    cp => cell_list(kcell)
    !cell%entrytime = tnow
    cp%birthtime = tnow
    !if (dividing) then  !???????????????
    !  cp%ID = 0
    !else
    !  lastID = lastID + 1
    !  cp%ID = lastID
    !endif
    cp%state = cell_stat%ALIVE
    cp%nspheres = 1
    cp%centre(:,1) = rsite
    cp%radius = plug_params%Raverage
    cp%ctype = ctype
    cp%tag = tag
    cp%step = 0
    !cell%lastdir = random_int(1,6,kpar)
    cp%dtotal = 0
    cp%mitosis = 0
    cp%Iphase = .true.



    call enter_exit(sub_name,2)
end subroutine create_cell
!
!###################################################################################
!

subroutine deallocate_abm_memory
    use arrays, only: cell_type,cell_list, num_cells,num_cell_f,cell_field
    use diagnostics, only: enter_exit,get_diagnostics_level

    !local variables
    integer :: AllocateStatus
    character(len=60) :: sub_name

    sub_name = 'deallocate_abm_memory'
    call enter_exit(sub_name,1)

    ! cells in the domain are stored in a list, while the
    ! occupancy array holds the indices of cells in the list.  When a cell
    ! leaves the domain or dies a gap is created in the list.
    ! The locations of such gaps are stored in the gaplist, the total number
    ! of such gaps is ngaps.  A cell entering the domain is allocated an index
    ! from the tail of this list, if ngaps > 0, or else it is added to the end of the cell list.
    if(allocated(cell_list))then
      deallocate (cell_list, STAT = AllocateStatus)
    endif
    if(allocated(cell_field))then
      deallocate (cell_field, STAT = AllocateStatus)
    endif

    call enter_exit(sub_name,2)
end subroutine deallocate_abm_memory

!
!###################################################################################
!

subroutine reallocate_abm_memory(current_cells)
    use arrays, only: dp,cell_type,cell_list, num_cells,num_cell_f,cell_field
    use diagnostics, only: enter_exit,get_diagnostics_level

    integer, intent(in) :: current_cells

    !local variables
    type(cell_type), allocatable :: cell_list_tmp(:) !List of cells
    real(dp),allocatable :: cell_field_tmp(:,:,:)

    integer :: AllocateStatus
    character(len=60) :: sub_name

    sub_name = 'reallocate_abm_memory'
    call enter_exit(sub_name,1)
    if(current_cells > num_cells) then
      allocate(cell_list_tmp(num_cells))
      cell_list_tmp = cell_list ! copy to temporary array
      deallocate(cell_list) !deallocate initially allocated memory
      allocate(cell_list(current_cells))
      cell_list(1:num_cells)=cell_list_tmp(1:num_cells)
      deallocate(cell_list_tmp) !deallocate the temporary array

      allocate(cell_field_tmp(num_cells,num_cell_f,3))
      cell_field_tmp = cell_field ! copy to temporary array
      deallocate(cell_field) !deallocate initially allocated memory
      allocate(cell_field(current_cells,num_cell_f,3))
      cell_field(1:num_cells,:,:)=cell_field_tmp(1:num_cells,:,:)
      deallocate(cell_field_tmp) !deallocate the temporary array
    else
      allocate(cell_list_tmp(current_cells))
      cell_list_tmp(1:current_cells) = cell_list(1:current_cells) ! copy to temporary array
      deallocate(cell_list) !deallocate initially allocated memory
      allocate(cell_list(current_cells))
      cell_list=cell_list_tmp
      deallocate(cell_list_tmp) !deallocate the temporary array

      allocate(cell_field_tmp(current_cells,num_cell_f,3))
      cell_field_tmp(1:current_cells,:,:) = cell_field(1:current_cells,:,:) ! copy to temporary array
      deallocate(cell_field) !deallocate initially allocated memory
      allocate(cell_field(current_cells,num_cell_f,3))
      cell_field=cell_field_tmp
      deallocate(cell_field_tmp) !deallocate the temporary array
    endif

    num_cells = current_cells

    call enter_exit(sub_name,2)
end subroutine reallocate_abm_memory

!
!##############################################################################
!

subroutine place_cells_plug

!-----------------------------------------------------------------------------------------
! Possible cell locations are in a close-packed rectangular block.  Locations outside the
! cylindrical vessel are eliminated, then cells in the assumed pre-existing gap are
! eliminated.
!
! The origin is at the centre of the tube, at the placenta end.
! The Z axis is along the tube centre-line.
!
! Tube dimensions:
!   radius = plug_params%tube_radius (um)
!   length = tube_length (um)
! Plug dimensions:
!   from   = plug_params%plug_zmin (um)
!   to     = plug_params%plug_zmax (um)
!   height = plug_params%plug_hmax (um)
!-----------------------------------------------------------------------------------------

    use arrays, only: dp,cell_type,cell_list, num_cells,plug_params, TROPHO_CELL
    use other_consts
    use diagnostics, only: enter_exit,get_diagnostics_level
    use cellpacker


    !local variables
    integer :: kcell, i, ix, iy, iz, nxx, nyy, nzz
    integer :: gen, region, ctype, tag
    real(dp) :: cell_radius, xrng, zrng, z0, dz, u, h
    real(dp) :: centre(3), rsite(3), block_centre(3), plug_centre(3), x, y, z, r
    integer, parameter :: out_unit=20
    integer :: diagnostics_level
    character(len=60) :: sub_name

    sub_name = 'place_cells_plug'
    call enter_exit(sub_name,1)
    call get_diagnostics_level(diagnostics_level)

    open(out_unit, file = 'cell_initialisation.txt',action='write')

    if (plug_params%use_packing .or. plug_params%use_makeRing) then
        z0 = (plug_params%plug_zmax + plug_params%plug_zmin)/2
        xrng = plug_params%tube_radius
        zrng = plug_params%plug_zmax - plug_params%plug_zmin
        dz = zrng/2
        plug_centre = [0.d0, 0.d0, z0]
        cell_radius = plug_params%Raverage
        kcell = 0
    else!if (plug_params%use_loosepack) then
        !these setting are needed for running the part consistent with MATLAB
        kcell = 0
        cell_radius = plug_params%Raverage
        dz=(plug_params%plug_zmax-plug_params%plug_zmin)
        xrng = plug_params%tube_radius
        nxx = (xrng/(2*cell_radius))*2
        nzz = dz/(2*cell_radius) + 1
    endif


    if (plug_params%use_packing) then
        nxx = 1.2*xrng/cell_radius
        nzz = 1.5*dz/cell_radius
        nyy = 1.3*nxx

        call SelectCellLocations(nxx, nyy, nzz, cell_radius, block_centre)
        if(diagnostics_level.gt.1)then
          write(*,*) 'nxx,nyy,nzz: ',nxx,nyy,nzz
          write(*,'(a,3f8.1)') 'plug_centre: ',plug_centre
          write(*,'(a,3f8.1)') 'block_centre: ',block_centre
        endif
        do i = 1,nxx*nyy*nzz
            if (cdist(i)+cell_radius > plug_params%tube_radius) cycle
            ix = xyz_lookup(i)%x
            iy = xyz_lookup(i)%y
            iz = xyz_lookup(i)%z
            centre = cloc(ix,iy,iz)%centre - block_centre + plug_centre
            if (centre(3) < plug_params%plug_zmin .or. centre(3) > plug_params%plug_zmax) cycle
            u = centre(3) - z0
            h = (1 + cos(u*PI/dz))*plug_params%plug_hmax/2
            kcell = kcell + 1
            rsite = centre
            if(diagnostics_level.gt.1)then
              write(*,'(2i6,3i4,3f8.1)') kcell,i,ix,iy,iz,rsite
            endif
            gen = 1
            region = 0
            ctype = TROPHO_CELL
            tag = 0
            if(kcell > num_cells)then
              call reallocate_abm_memory(kcell)
            endif
            call create_cell(kcell,rsite,ctype,gen,tag,region,.false.)
            if(diagnostics_level.gt.1)then
              write(out_unit,'(2i6,3i4,5f8.1)') kcell,i,ix,iy,iz,cell_list(kcell)%centre(:,1)
            endif
        enddo
        call FreeCellLocations
    elseif (plug_params%use_loosepack) then
        do ix = 0,nxx,2
            do iy = 0,nxx,2
                do iz = 0,nzz,2
                    x = -xrng+cell_radius+(ix*2+mod(REAL(iy+iz,dp),2.))*cell_radius!-plug_params%tube_radius+cell_radius+ix*2*cell_radius
                    y = -xrng+cell_radius+sqrt(3.)*(iy+1/3.*mod(REAL(iz,dp),2.))*cell_radius!-plug_params%tube_radius+cell_radius+iy*2*cell_radius
                    z = plug_params%plug_zmin+(iz*2*sqrt(6.)/3)*cell_radius!plug_params%plug_zmin+iz*2*cell_radius

                    r = sqrt(x*x + y*y)
                    if ((r+cell_radius<plug_params%tube_radius) .AND. (z>plug_params%plug_zmin .OR. z<plug_params%plug_zmax)) then
                        kcell = kcell + 1
                        rsite = [x, y, z]
                        gen = 1
                        region = 0
                        ctype = TROPHO_CELL
                        tag = 0
                        if(kcell > num_cells)then
                          call reallocate_abm_memory(kcell)
                        endif
                        call create_cell(kcell,rsite,ctype,gen,tag,region,.false.)
                    endif
                enddo
            enddo
        enddo
    elseif (plug_params%use_makeRing) then
        nxx = plug_params%tube_radius/(2*cell_radius) + 1
        nzz = INT(dz/(2*cell_radius) )+ 2 !+1
        do ix = -nxx,nxx
            do iy = -nxx,nxx
                x = ix*2*cell_radius    ! (x,y,z) is the offset from the centre of the plug
                y = iy*2*cell_radius
                r = sqrt(x*x + y*y)
                if (r + cell_radius > plug_params%tube_radius) cycle
                do iz = -nzz,nzz,2
                    z = iz*2*cell_radius
                    if (z < -dz .or. z > dz) cycle
                    h = (1 + cos(z*PI/dz))*plug_params%plug_hmax/2
                    if (r < plug_params%tube_radius - h) cycle
                    kcell = kcell + 1
                    rsite = [x, y, z0+z]
                    gen = 1
                    region = 0
                    ctype = TROPHO_CELL
                    tag = 0
                    if(kcell > num_cells)then
                      call reallocate_abm_memory(kcell)
                    endif
                    call create_cell(kcell,rsite,ctype,gen,tag,region,.false.)
                enddo
            enddo
        enddo
    else
        !!!!!!this code that generates packed cells is consistent with my MATLAB code!!!!!
        do ix = 0,nxx
            do iy = 0,nxx
                do iz = 0,nzz
                    x = -xrng+cell_radius+(ix*2+mod(REAL(iy+iz,dp),2.))*cell_radius!-plug_params%tube_radius+cell_radius+ix*2*cell_radius
                    y = -xrng+cell_radius+sqrt(3.)*(iy+1/3.*mod(REAL(iz,dp),2.))*cell_radius!-plug_params%tube_radius+cell_radius+iy*2*cell_radius
                    z = plug_params%plug_zmin+(iz*2*sqrt(6.)/3)*cell_radius!plug_params%plug_zmin+iz*2*cell_radius

                    r = sqrt(x*x + y*y)
                    if ((r+cell_radius<plug_params%tube_radius) .AND. (z>plug_params%plug_zmin .OR. z<plug_params%plug_zmax)) then
                        kcell = kcell + 1
                        rsite = [x, y, z]
                        gen = 1
                        region = 0
                        ctype = TROPHO_CELL
                        tag = 0
                        if(kcell > num_cells)then
                          call reallocate_abm_memory(kcell)
                        endif
                        call create_cell(kcell,rsite,ctype,gen,tag,region,.false.)
                    endif
                enddo
            enddo
        enddo

    endif
    close(out_unit)
    !if(kcell < num_cells)then
    !  call reallocate_abm_memory(kcell)
    !endif
    num_cells = kcell

    call enter_exit(sub_name,2)
end subroutine place_cells_plug
!
!##############################################################################
!
subroutine setup_nbrlists
    !-----------------------------------------------------------------------------------------
    ! Dumb version to start
    ! Note that now we seek neighbour only for cell_stat%ALIVE cells, but the neighbours can include
    ! WALL cells (which are immobile), because they also exert force.
    !-----------------------------------------------------------------------------------------
    use arrays, only: dp,cell_type,cell_list, num_cells,plug_params, cell_stat,neighbour_type,abm_control
    use other_consts
    use diagnostics, only: enter_exit,get_diagnostics_level
    type(cell_type), pointer :: cp1, cp2
    real(dp) :: R1, c1(3), R2, c2(3), v(3), d2, d, Rsum, dfactor, dmin,nearest_dist
    integer :: kcell, k2, k, isphere1, nspheres1, isphere2, nspheres2, nbrs
    type(neighbour_type) :: nbrlist(100)
    logical :: near, incontact, contact(2,2)
    logical :: dbug = .false.
    integer :: diagnostics_level
    character(len=60) :: sub_name

    sub_name = 'setup_nbrlists'
    call enter_exit(sub_name,1)
    call get_diagnostics_level(diagnostics_level)

    dmin = 1.0e10
    do kcell = 1,num_cells
        cp1 => cell_list(kcell)
        if (cp1%state /= cell_stat%ALIVE) cycle
        !if (cp1%Iphase) then
        !   nspheres1 = 1
        !else
        !   nspheres1 = 2
        !endif
        nspheres1 = cp1%nspheres
        nbrs = 0
        nearest_dist = 1.0e10_dp
        do k2 = 1,num_cells
            if (k2 == kcell) cycle
            cp2 => cell_list(k2)
            if (cp2%state > cell_stat%ALIVE) cycle
            nspheres2 = cp2%nspheres
            near = .false.
            incontact = .false.
            contact = .false.
            do isphere1 = 1,nspheres1
                R1 = cp1%radius(isphere1)
                c1 = cp1%centre(:,isphere1)
                do isphere2 = 1,nspheres2
                    R2 = cp2%radius(isphere2)
                    c2 = cp2%centre(:,isphere2)
                    Rsum = R1 + R2
                    v = c2 - c1
                    d2 = dot_product(v,v)
                    d = sqrt(d2)
                    dmin = min(d,dmin) !Overall minimum distance

                    !if(diagnostics_level.gt.1)then
                    !    write(*,'(a,9f8.1)') 'c1,R1,c2,R2,d: ',c1,R1,c2,R2,d
                    !endif
                    !if (use_hysteresis) then
                    !    if (d < plug_params%d_nbr_limit) then
                    !        if (dbug) write(*,'(5f8.2)') R1,R2,Rsum,d,plug_params%d_nbr_limit
                    !        dfactor = 1
                    !        do k = 1,cp1%nbrs
                    !            if (cp1%nbrlist(k)%indx == k2) then
                    !                if (cp1%nbrlist(k)%contact(isphere1,isphere2)) dfactor = k_detach
                    !                exit
                    !            endif
                    !        enddo
                    !        near = .true.
                    !        if (d < dfactor*Rsum) then
                    !            incontact = .true.
                    !            contact(isphere1,isphere2) = .true.
                    !        endif
                    !    endif
                    !else
                    if (d < plug_params%d_nbr_limit) then
                        near = .true.
                        exit
                    endif
                    !endif
                enddo
            enddo
            if (near) then
                nbrs = nbrs + 1
                if (nbrs > 100) then
                    write(*,*) 'Size of nbrlist exceeded:'
                    write(*,*)  plug_params%istep, kcell, 'ncells: ',plug_params%istep,kcell,num_cells,nbrs,100
                    write(*,'(10i6)') nbrlist(:)%indx
                    stop
                endif
                nbrlist(nbrs)%indx = k2
                nbrlist(nbrs)%contact = contact
                nbrlist(nbrs)%distance = d
                if(dmin.lt.nearest_dist)then
                  nearest_dist = d
                endif
            !           if (cp2%state == WALL) write(*,*) 'WALL nbr: ',kcell,k2
            endif
        enddo
        cp1%nbrs = nbrs
        cp1%nbrlist(1:nbrs) = nbrlist(1:nbrs)
        cp1%nearest_dist = nearest_dist

    enddo

    abm_control%delta_max = max(0.5*dmin-abm_control%delta_min, abm_control%delta_min) !This means cells can never'hit each other' and ensures sep

    call enter_exit(sub_name,2)
end subroutine setup_nbrlists

!
!###################################################################################
!

subroutine move_cells_force(cell_population,kdrag,input_dt)
    use arrays, only: dp,num_cells,cell_list,cell_stat,cell_field,abm_control,num_cell_f
    use math_utilities, only: vector_length
    use diagnostics, only: enter_exit,get_diagnostics_level
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_MOVE_CELLS_FORCE" :: CALC_MOVE_CELLS_FORCE

    integer, intent(in) :: cell_population
    real(dp), intent(in) :: kdrag
    real(dp), intent(in) :: input_dt

    integer :: kcell, kforce
    real(dp) :: max_force,force_mag, dt_move,max_dx
    integer :: k1, kpar, nd, nr, nc(0:8), kfrom(0:8), kto(0:8), tMnodes
    real(dp), allocatable :: force(:,:)
    real(dp) :: fmax, dx1(3), dx2(3)

    integer :: diagnostics_level
    character(len=60) :: sub_name

    sub_name = 'move_cells_force'
    call enter_exit(sub_name,1)
    call get_diagnostics_level(diagnostics_level)

    !Initialise timestep for move to current time step
    dt_move = input_dt !ndt*delta_tmove
     !!write(*,*) 'fmover: ',ndt,delta_tmove,dt_move
    allocate(force(num_cells,3))
    force = 0.0
    max_force = -1.0e6_dp
    do kcell = 1, num_cells
      do kforce = 1,num_cell_f
        force(kcell,:) = force(kcell,:) + cell_field(kcell,kforce,:)
      enddo
      force_mag =vector_length(force(kcell,:))
      !write(*,*) force_mag/kdrag
      if(force_mag > max_force)then
        max_force = force_mag
      endif
    enddo

    if(diagnostics_level>1)write(*,*) 'Max force', max_force
    write(*,*) 'Max force', max_force
    if (max_force > 0) then
      max_dx = dt_move*max_force/kdrag
      do while (max_dx > abm_control%delta_max)
        dt_move = dt_move/2.0_dp
        max_dx = dt_move*max_force/kdrag
      enddo
    endif
    if(dt_move < abm_control%delta_t_min)then
      dt_move = abm_control%delta_t_min
    endif

    do kcell = 1, num_cells
        if (cell_list(kcell)%state == cell_stat%GONE_BACK .or. cell_list(kcell)%state == cell_stat%GONE_THROUGH) cycle
        if (cell_list(kcell)%state == cell_stat%ALIVE)then
          if (cell_list(kcell)%Iphase) then
            dx1 = dt_move*force(kcell,:)/kdrag
            !write(*,*) dx1(:),force(kcell,:)/kdrag
            cell_list(kcell)%centre(:,1) = cell_list(kcell)%centre(:,1) + dx1
          else
          !need to implement so this works for a mitosing cell
          !dx1 = dt_move*force(:,k1,1)/kdrag
          !dx2 = dt_move*force(:,k1,2)/kdrag
          !cp%centre(:,1) = cp%centre(:,1) + dx1
          !cp%centre(:,2) = cp%centre(:,2) + dx2
          endif
        endif
!        if (cp%centre(3,1)>tube_length) then
!         !   call Tcell_death(kcell)
!            cp%state=GONE_BACK
!        elseif (cp%centre(3,1)<0) then
!            cp%state=GONE_THROUGH
!        endif
!    enddo
    enddo
!!omp end parallel do
    deallocate(force)

    abm_control%current_time = abm_control%current_time + dt_move
    abm_control%used_delta_t = dt_move
    if(diagnostics_level.gt.1)then
        write(*,*) 'current_time', abm_control%current_time
    endif
    call enter_exit(sub_name,2)
end subroutine move_cells_force

function get_current_t() result(time_params)
    use arrays, only: dp,abm_control
    use diagnostics, only: enter_exit,get_diagnostics_level
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_get_current_t" :: get_current_t

    real(dp) :: time_params
    character(len=60) :: sub_name

    sub_name = 'get_current_t'
    call enter_exit(sub_name,1)

    time_params = abm_control%current_time
    call enter_exit(sub_name,2)

end function get_current_t

function get_current_dt() result(time_params)
    use arrays, only: dp,abm_control
    use diagnostics, only: enter_exit,get_diagnostics_level
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_RETURN_CURRENT_DT" :: RETURN_CURRENT_DT

    real(dp) :: time_params
    character(len=60) :: sub_name

    sub_name = 'get_current_dt'
    call enter_exit(sub_name,1)

    time_params = abm_control%used_delta_t

    call enter_exit(sub_name,2)

end function get_current_dt


end module abm_models
