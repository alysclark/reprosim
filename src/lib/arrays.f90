module arrays
!*Description:* This module defines arrays.
!

  implicit none

  integer :: num_elems,num_nodes,num_units,maxgen,num_arterial_elems,num_cells,num_cell_f
  integer :: num_faces,num_inlet_faces,num_outlet_faces,num_wall_faces

  integer, parameter :: dp=kind(0.d0) !  for double precision

  integer,allocatable :: nodes(:) !allocated in define_node_geometry
  integer,allocatable :: elems(:) !allocated in define_1d_elements
  integer,allocatable :: elem_cnct(:,:,:)  !NXI(-ni:ni,1,ne)
  integer,allocatable :: elem_nodes(:,:)
  integer,allocatable :: elem_ordrs(:,:)
  integer,allocatable :: elem_symmetry(:)
  integer,allocatable :: elem_units_below(:)
  integer,allocatable :: elems_at_node(:,:)
  integer,allocatable :: units(:)
  integer, allocatable :: internal_faces(:,:)
  integer, allocatable :: inlet_faces(:,:)
  integer, allocatable :: outlet_faces(:,:)
  integer, allocatable :: wall_faces(:,:)
  integer, allocatable :: elem_3d(:,:)
  real(dp), allocatable :: node_3d(:,:)


  real(dp),allocatable :: cell_field(:,:,:) !properties of cells
  real(dp),allocatable :: elem_field(:,:) !properties of elements
  real(dp),allocatable :: elem_direction(:,:)
  real(dp),allocatable :: node_xyz(:,:)
  real(dp),allocatable :: unit_field(:,:) !properties of elastic units
  real(dp),allocatable :: node_field(:,:)

! temporary, for debugging:
  real(dp) :: unit_before

  type neighbour_type
    integer :: indx
!   logical :: incontact
    logical*1 :: contact(2,2)
    real(dp) :: distance
  end type

  type cell_type
    integer :: ID            ! cell ID
    integer :: state         !
    integer :: site(3)
    logical :: Iphase        !Relates to cell cycle - Iphase is where the cell grows and prepares for mitosis
    integer :: nspheres      ! =1 for Iphase, =2 for Mphase (Mphase is mitosis phase)
    real(dp) :: V            ! actual volume (um3)
    real(dp) :: dVdt         ! Change in cell volume in time
    real(dp) :: radius(2)    ! sphere radii (um) - could be two of them if in Mphase
    real(dp) :: centre(3,2)  ! sphere centre positions
    real(dp) :: centre_t_n(3,2)  ! sphere centre positions at last exported time
    real(dp) :: centre_t_0(3,2)  ! sphere centre positions at seeding
    real(dp) :: d            ! centre separation distance (um)
    real(dp) :: birthtime    ! time of birth
    real(dp) :: t_start_mitosis !Time mitosis started
    real(dp) :: mitosis      ! level of mitosis (0 - 1)
    real(dp) :: V_divide
    real(dp) :: d_divide     ! centre separation distance at the end of mitosis
    integer :: step
    integer :: tag
    integer(2) :: ctype
    integer(2) :: lastdir
    integer :: dtotal(3)
    integer :: nbrs
    real(dp) :: nearest_dist !distance from nearest neighbour
    real(dp) :: wall_distance !distance from wall
    real(dp) :: wall_dir(3) !wall direction
    type(neighbour_type) :: nbrlist(100)
  end type cell_type

  type(cell_type), allocatable, target :: cell_list(:) !List of cells

  type abm_control_params
    real(dp) :: delta_t_min         !minimum time step
    real(dp) :: delta_t             !default time strp
    real(dp) :: delta_min = 0.31_dp    !minimum allowed separation between cells (fraction of cell radius)
    real(dp) :: d_nbr_limit = 2.0_dp*2.0_dp*20.0_dp
    real(dp) :: delta_max = 0.10_dp    !maximum movement in a single timestep
    real(dp) :: current_time
    real(dp) :: used_delta_t
    logical :: Wall = .True. !Is there a wall in the model?
    !delta_tmove = dt_min
    !ndt = 5
  end type abm_control_params
  type(abm_control_params) :: abm_control

  integer, parameter :: NCTYPES = 1
  integer, parameter :: TROPHO_CELL = 1
  integer, parameter :: MAX_CELLTYPES = 1


  type cell_status
    integer :: WALL = 0
    integer :: ALIVE = 1
    integer :: DEAD = 2
    integer :: GONE_BACK = 3
    integer :: GONE_THROUGH = 4
    integer :: TAGGED_CELL   = 99
  end type cell_status

  type(cell_status) :: cell_stat

  type plug_model_properties
    real(dp) :: beta                         ! speed: 0 < beta < 1       (0.65)
    real(dp) :: rho                          ! persistence: 0 < rho < 1  (0.95)
    real(dp) :: days                         ! number of days to simulate
    real(dp) :: settle_hrs                   ! number of hours for settling
    real(dp) :: seed(2)                      ! seed vector for the random number generators
    !real(dp) :: ncpu_dummy                   ! # of processors - not used at the moment
    real(dp) :: tube_length = 1000.0_dp                 ! length of tube (um)
    real(dp) :: tube_radius = 200.0_dp                 ! radius of tube (um)
    real(dp) :: plug_zmin  = 100.0_dp                  ! plug z from (um)
    real(dp) :: plug_zmax  = 600.0_dp                  ! plug z to (um)
    real(dp) :: plug_hmax   = 500.0_dp                ! plug height (um)
    real(dp) :: Raverage = 20.0_dp                    ! average cell radius (um)
    real(dp) :: ntgui                        ! interval between GUI outputs (timesteps)
    real(dp) :: idelay                       ! simulation step delay (ms)
    real(dp) :: ichemo_1
    real(dp) :: grad_amp(2)                  ! chemokine gradient amplitude
    real(dp) :: coef1                        ! Z chemotaxis1 exp(coef1*)
    real(dp) :: coef2                        ! Radial chemotaxis1 exp(coef2*)
    real(dp) :: BG_flow_amp                  ! background velocity amplitude (um/min)
    real(dp) :: inletPressure                ! inlet pressure
    real(dp) :: a_separation
    real(dp) :: a_force
    real(dp) :: c_force
    real(dp) :: x0_force
    real(dp) :: x1_force
    real(dp) :: kdrag
    real(dp) :: frandom
    real(dp) :: n_cell_positions             ! number of cell positions to save each time step
    integer :: Nsteps
    integer :: nsteps_per_min
    integer :: istep
    integer :: lastID
    integer :: ncells
    integer :: nwallcells
    integer :: nlist
    integer :: ndt
    logical :: use_packing = .True.
    logical :: use_loosepack = .False.
    logical :: use_makeRing = .False.
  end type plug_model_properties

  type(plug_model_properties) :: plug_params

  private
  public cell_type, cell_list, set_node_field_value, cell_field, elem_field,internal_faces,&
    num_elems, elem_nodes, node_xyz, nodes, elems, num_faces, inlet_faces, num_inlet_faces, &
    num_outlet_faces,num_wall_faces,outlet_faces,wall_faces,elem_3d,node_3d,&
    num_cells, num_cell_f, num_nodes, units, num_units, unit_field, node_field, dp, elem_cnct, elem_ordrs, elem_direction, &
    elems_at_node, elem_symmetry, elem_units_below, maxgen, num_arterial_elems,plug_params
  public NCTYPES, TROPHO_CELL, MAX_CELLTYPES,cell_stat,neighbour_type,abm_control

contains
  subroutine set_node_field_value(row, col, value)  
  !*Description:* This subroutine sets the value of a node field
    implicit none
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_SET_NODE_FIELD_VALUE" :: SET_NODE_FIELD_VALUE
    integer, intent(in) :: row, col
    real(dp), intent(in) :: value
    
    node_field(row, col) = value
	
  end subroutine set_node_field_value


end module arrays
