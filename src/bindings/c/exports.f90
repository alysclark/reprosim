module exports_c
  implicit none

  private

contains
!!!################################################################

  subroutine export_1d_elem_field_c(ne_field, EXELEMFILE, filename_len, group_name, group_name_len, field_name, field_name_len) &
    bind(C, name="export_1d_elem_field_c")

    use iso_c_binding, only: c_ptr
    use utils_c, only: strncpy
    use exports, only: export_1d_elem_field
    use other_consts, only: MAX_STRING_LEN, MAX_FILENAME_LEN
    implicit none
    integer,intent(in) :: ne_field, filename_len, group_name_len, field_name_len
    type(c_ptr), value, intent(in) :: EXELEMFILE, group_name, field_name
    character(len=MAX_FILENAME_LEN) :: filename_f
    character(len=MAX_STRING_LEN) :: group_name_f, field_name_f

    call strncpy(filename_f, EXELEMFILE, filename_len)
    call strncpy(group_name_f, group_name, group_name_len)
    call strncpy(field_name_f, field_name, field_name_len)

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_export_1d_elem_field(ne_field, filename_f, group_name_f, field_name_f)
#else
    call export_1d_elem_field(ne_field, filename_f, group_name_f, field_name_f)
#endif

  end subroutine export_1d_elem_field_c

!!!############################################################################

  subroutine export_1d_elem_geometry_c(EXELEMFILE, filename_len, name, name_len) bind(C, name="export_1d_elem_geometry_c")

    use iso_c_binding, only: c_ptr
    use utils_c, only: strncpy
    use exports, only: export_1d_elem_geometry
    use other_consts, only: MAX_STRING_LEN, MAX_FILENAME_LEN
    implicit none
    integer,intent(in) :: filename_len, name_len
    type(c_ptr), value, intent(in) :: EXELEMFILE, name
    character(len=MAX_FILENAME_LEN) :: filename_f
    character(len=MAX_STRING_LEN) :: name_f

    call strncpy(filename_f, EXELEMFILE, filename_len)
    call strncpy(name_f, name, name_len)

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_export_1d_elem_geometry(filename_f, name_f)
#else
    call export_1d_elem_geometry(filename_f, name_f)
#endif

  end subroutine export_1d_elem_geometry_c


!!!##########################################################################

  subroutine export_node_geometry_c(EXNODEFILE, filename_len, name, name_len) bind(C, name="export_node_geometry_c")

    use iso_c_binding, only: c_ptr
    use utils_c, only: strncpy
    use other_consts, only: MAX_STRING_LEN, MAX_FILENAME_LEN
    use exports, only: export_node_geometry
    implicit none

    integer,intent(in) :: filename_len, name_len
    type(c_ptr), value, intent(in) :: EXNODEFILE, name
    character(len=MAX_FILENAME_LEN) :: filename_f
    character(len=MAX_STRING_LEN) :: name_f

    call strncpy(filename_f, EXNODEFILE, filename_len)
    call strncpy(name_f, name, name_len)

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_export_node_geometry(filename_f, name_f)
#else
    call export_node_geometry(filename_f, name_f)
#endif

  end subroutine export_node_geometry_c

  !!!########################################################################

  subroutine export_terminal_perfusion_c(EXNODEFILE, filename_len, name, name_len) bind(C, name="export_terminal_perfusion_c")

    use iso_c_binding, only: c_ptr
    use utils_c, only: strncpy
    use exports, only: export_terminal_perfusion
    use other_consts, only: MAX_STRING_LEN, MAX_FILENAME_LEN
    implicit none
    integer,intent(in) :: filename_len, name_len
    type(c_ptr), value, intent(in) :: EXNODEFILE, name
    character(len=MAX_FILENAME_LEN) :: filename_f
    character(len=MAX_STRING_LEN) :: name_f

    call strncpy(filename_f, EXNODEFILE, filename_len)
    call strncpy(name_f, name, name_len)

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_export_terminal_perfusion(filename_f, name_f)
#else
    call export_terminal_perfusion(filename_f, name_f)
#endif

  end subroutine export_terminal_perfusion_c


!!! #################################################################

  subroutine export_node_field_c(nj_field, EXNODEFIELD, filename_len, name, name_len, field_name, field_name_len) &
    bind(C, name="export_node_field_c")

    use iso_c_binding, only: c_ptr
    use utils_c, only: strncpy
    use exports, only: export_node_field
    use other_consts, only: MAX_STRING_LEN, MAX_FILENAME_LEN
    implicit none
    integer,intent(in) :: nj_field, filename_len, name_len, field_name_len
    type(c_ptr), value, intent(in) :: EXNODEFIELD, name, field_name
    character(len=MAX_FILENAME_LEN) :: filename_f
    character(len=MAX_STRING_LEN) :: name_f, field_name_f

    call strncpy(filename_f, EXNODEFIELD, filename_len)
    call strncpy(name_f, name, name_len)
    call strncpy(field_name_f, field_name, field_name_len)

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_export_node_field(nj_field, filename_f, name_f, field_name_f)
#else
    call export_node_field(nj_field, filename_f, name_f, field_name_f)
#endif

  end subroutine export_node_field_c

!!!##########################################################################

  subroutine export_cell_location_c(filename, filename_len, cell_population) bind(C, name="export_cell_location_c")

    use iso_c_binding, only: c_ptr
    use utils_c, only: strncpy
    use other_consts, only: MAX_STRING_LEN, MAX_FILENAME_LEN
    use exports, only: export_cell_location
    implicit none

    integer,intent(in) :: filename_len
    type(c_ptr), value, intent(in) :: filename
    integer, intent(in) :: cell_population
    character(len=MAX_FILENAME_LEN) :: filename_f

    call strncpy(filename_f, filename, filename_len)

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_export_cell_location(filename_f, cell_population)
#else
    call export_cell_location(filename_f, cell_population)
#endif

  end subroutine export_cell_location_c



!!!##########################################################################

  subroutine export_cell_plug_c(filename, filename_len, cell_population, cur_time,velocity,flow_rate) &
      bind(C, name="export_cell_plug_c")

    use iso_c_binding, only: c_ptr
    use utils_c, only: strncpy
    use other_consts, only: MAX_STRING_LEN, MAX_FILENAME_LEN
    use arrays, only: dp
    use exports, only: export_cell_plug
    implicit none

    integer,intent(in) :: filename_len
    type(c_ptr), value, intent(in) :: filename
    integer, intent(in) :: cell_population
    character(len=MAX_FILENAME_LEN) :: filename_f
    real(dp) :: cur_time,velocity, flow_rate

    call strncpy(filename_f, filename, filename_len)

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_export_cell_plug(filename_f, cell_population,cur_time,velocity,flow_rate)
#else
    call export_cell_plug(filename_f, cell_population,cur_time,velocity,flow_rate)
#endif

  end subroutine export_cell_plug_c

!!!##########################################################################

  subroutine export_cell_exnode_c(filename, filename_len, cell_population) bind(C, name="export_cell_exnode_c")

    use iso_c_binding, only: c_ptr
    use utils_c, only: strncpy
    use other_consts, only: MAX_STRING_LEN, MAX_FILENAME_LEN
    use exports, only: export_cell_exnode
    implicit none

    integer,intent(in) :: filename_len
    type(c_ptr), value, intent(in) :: filename
    integer, intent(in) :: cell_population
    character(len=MAX_FILENAME_LEN) :: filename_f

    call strncpy(filename_f, filename, filename_len)

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_export_cell_exnode(filename_f, cell_population)
#else
    call export_cell_exnode(filename_f, cell_population)
#endif

  end subroutine export_cell_exnode_c

end module exports_c
