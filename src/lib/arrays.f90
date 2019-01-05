module arrays
!*Description:* This module defines arrays.
!

  implicit none

  integer :: num_elems,num_nodes,num_units,maxgen,num_arterial_elems

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

  real(dp),allocatable :: elem_field(:,:) !properties of elements
  real(dp),allocatable :: elem_direction(:,:)
  real(dp),allocatable :: node_xyz(:,:)
  real(dp),allocatable :: unit_field(:,:) !properties of elastic units
  real(dp),allocatable :: node_field(:,:)

! temporary, for debugging:
  real(dp) :: unit_before

  type cell_type
    integer :: ID
    integer :: state
    integer :: site(3)
    logical :: Iphase
    integer :: nspheres             ! =1 for Iphase, =2 for Mphase
    real(dp) :: V            ! actual volume (um3)
    real(dp) :: dVdt
    real(dp) :: radius(2)    ! sphere radii (um)
    real(dp) :: centre(3,2)  ! sphere centre positions
    real(dp) :: d            ! centre separation distance (um)
    real(dp) :: birthtime
    real(dp) :: t_start_mitosis
    real(dp) :: mitosis      ! level of mitosis (0 - 1)
    real(dp) :: V_divide
    real(dp) :: d_divide     ! centre separation distance at the end of mitosis
    integer :: step
    integer :: tag
  !NULLIFY() instead.
    integer(2) :: ctype
    integer(2) :: lastdir
    integer :: dtotal(3)
    integer :: nbrs
    !type(neighbour_type) :: nbrlist(100)
end type

  private
  public set_node_field_value, elem_field, num_elems, elem_nodes, node_xyz, nodes, elems, &
    num_nodes, units, num_units, unit_field, node_field, dp, elem_cnct, elem_ordrs, elem_direction, &
    elems_at_node, elem_symmetry, elem_units_below, maxgen, num_arterial_elems

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
