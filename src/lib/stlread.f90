 module stlread

    use stla_io

    contains

    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    subroutine setup_mesh(filename)

        character*(*) :: filename
        integer ( kind = 4 ), allocatable, dimension(:,:) :: face_node
        real ( kind = 8 ), allocatable, dimension(:,:) :: face_normal
        integer ( kind = 4 ) face_num
        integer ( kind = 4 ) ierror
        integer ( kind = 4 ) node_num
        real ( kind = 8 ), allocatable, dimension(:,:) :: node_xyz
        integer ( kind = 4 ) solid_num
        integer ( kind = 4 ) text_num

        call stla_size(filename,solid_num, node_num, face_num, text_num)
        call stla_size_print ( filename, solid_num, node_num, face_num, & text_num )

        allocate ( face_node(3,face_num) )
        allocate ( face_normal(3,face_num) )
        allocate ( node_xyz(3,node_num) )

        call stla_read ( input_file_name, node_num, face_num, node_xyz, & face_node, face_normal, ierror )

        if ( ierror /= 0 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a,i8)' ) '  STLA_READ returned IERROR = ', ierror
            return
        end if

        call stla_size_print ( input_file_name, solid_num, node_num, face_num, & text_num )

        call stla_face_node_print ( face_num, face_node )
        call stla_face_normal_print ( face_num, face_normal )
        call stla_node_xyz_print ( node_num, node_xyz )

        deallocate ( face_node )
        deallocate ( face_normal )
        deallocate ( node_xyz )

    return

    end subroutine

end module
