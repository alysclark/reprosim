module mix_fem
  use arrays, only: dp
  implicit none
  !Module parameters

  !Module types

  !Module depvar

  !Interfaces
  private
  public assemble_sparse_matrices
  public read_b_matrix
  public read_e2face
  public create_sampling_grid
contains

!
!############################################################
!
subroutine assemble_sparse_matrices
    use arrays, only: num_elems,B_MATRIX,all_faces,num_all_faces,element2face
    use diagnostics, only: enter_exit,get_diagnostics_level
    use indices
    use other_consts, only: MAX_FILENAME_LEN
    implicit none
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_ASSEMBLE_SPARSE_MATRICES" :: READ_ASSEMBLE_SPARSE_MATRICES


!integer :: nel
!real(dp) :: B(6,6)
!real(dp) :: k(:), BB(6*6*nel-CommonFaceNo), C(6)
!real(dp) :: CC(6*nel)
!real(dp) :: A(8*6*nel-CommonFaceNo)
!integer :: IB(6*6*nel-CommonFaceNo), JB(6*6*nel-CommonFaceNo)!, ones(6)
!integer :: IC(6*nel), JC(6*nel),All_Surfaces(:,:), element2face(:,:)
!integer :: IA(8*6*nel-CommonFaceNo),JA(8*6*nel-CommonFaceNo)
!!real(dp), allocatable :: Ak(:,:)
    integer :: ii,j, I(6), signum(6), DiagSig(6,6)
    integer :: counter, Recount1, Ccounter, commonCount
!logical :: ok1
!character * ( 255 ) :: fileplace, Bmatrix_filename

    character(len=60) :: sub_name
    integer :: diagnostics_level

    sub_name = 'assemble_sparse_matrices'
    call enter_exit(sub_name,1)
    call get_diagnostics_level(diagnostics_level)


    Recount1=0
    counter=0
    Ccounter=0
    commonCount=0
    do j = 1,num_elems
      I(:) = element2face(j,:) ! diag matrix of face numbers of the element
                         ! gets face numbers in the element
      !print *, "I=", I(:)
      signum=(/1,1,1,1,1,1/)
      do ii=1,6
        if (j==all_faces(I(ii),6)) then
            signum(ii)=-1 ! -1 if edge is on T-
        endif
      enddo
    DiagSig(:,:) = diag(signum(:),size(signum(:)))
    !B(:,:) = blood_viscosity/k(j)*matmul(matmul(DiagSig, B_MATRIX(j)%Inside_B),DiagSig)



    enddo

!Each element has a b-matrix, etc,which Rojan has assembled as a number of files to read in.
!for now we will have to do this because we don't have the code that she wrote for this

  call enter_exit(sub_name,1)
end subroutine assemble_sparse_matrices

!
!############################################################
!

subroutine read_b_matrix(filename,filename_len)
   use arrays, only: num_elems,B_MATRIX
    use diagnostics, only: enter_exit,get_diagnostics_level
    use indices
    use other_consts, only: MAX_FILENAME_LEN
    implicit none
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_READ_B_MATRIX" :: READ_B_MATRIX

    character(len=MAX_FILENAME_LEN), intent(in) :: filename !Input nodefile
    integer, intent(in) :: filename_len

    character(len=MAX_FILENAME_LEN) :: appended_fn
    character(len = filename_len) ::trimfilename
    integer :: ne

    character(len=60) :: sub_name
    integer :: diagnostics_level

    sub_name = 'read_b_matrix'
    call enter_exit(sub_name,1)
    call get_diagnostics_level(diagnostics_level)

    trimfilename = filename(1:filename_len)
    allocate(B_MATRIX(num_elems))
    do ne = 1,num_elems
      write(appended_fn,'(A100,I0,".txt")')trimfilename,ne
      call read_real_array(adjustl(appended_fn),6,6,B_MATRIX(ne)%Inside_B)

      write(*,*) appended_fn,B_MATRIX(ne)%Inside_B
    enddo

    call enter_exit(sub_name,2)
end subroutine read_b_matrix


!
!############################################################
!

subroutine read_e2face(filename,filename_len)
   use arrays, only: num_elems,element2face
    use diagnostics, only: enter_exit,get_diagnostics_level
    use indices
    use other_consts, only: MAX_FILENAME_LEN
    implicit none
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_READ_B_MATRIX" :: READ_B_MATRIX

    character(len=MAX_FILENAME_LEN), intent(in) :: filename !Input nodefile
    integer, intent(in) :: filename_len

    character(len=MAX_FILENAME_LEN) :: appended_fn
    character(len = filename_len) ::trimfilename
    integer :: ne

    character(len=60) :: sub_name
    integer :: diagnostics_level

    sub_name = 'read_e2face'
    call enter_exit(sub_name,1)
    call get_diagnostics_level(diagnostics_level)

    trimfilename = filename(1:filename_len)
    allocate(element2face(num_elems,6))

    call read_integer_array(adjustl(trimfilename),num_elems,6,element2face)

    write(*,*) element2face(1,:)
    write(*,*) element2face(num_elems,:)


    call enter_exit(sub_name,2)
end subroutine read_e2face

!
!#########################################################
!
!subroutine find_cell_mesh(elZ, nel_z,N_Cell,Element,vertices,el,cell_elem)
!integer :: nel_z, N_Cell,i, m, k,e
!integer :: Element(:,:), el(:,:)
!integer :: cell_elem(:)
!integer :: base(6,4)
!real (dp):: vertices(:,:),elZ(:)
!real (dp):: normal(6,3),dP(6)
!
!do i=1,N_Cell
!    !if (cell_list(i)%state /= ALIVE) cycle
!    if (cell_list(i)%state == GONE_BACK .or. cell_list(i)%state == GONE_THROUGH) cycle
!    do k=2,nel_z+1
!        if ( cell_list(i)%centre(3,1)<=elZ(k) .AND. &
!            cell_list(i)%centre(3,1)>=elZ(k-1) ) then
!            !print *, "!!", el(k-1,:)
!            do m=1, size(el,2)
!                e=el(k-1,m)
!                base(1,1) = Element(e,1)
!                base(1,2) = Element(e,4)
!                base(1,3) = Element(e,8)
!                base(1,4) = Element(e,5)
!                normal(1,:) = cross(&
!                (vertices(base(1,2),:)-vertices(base(1,1),:)), &
!                (vertices(base(1,3),:)-vertices(base(1,1),:)))
!                dP(1) = dot_product(normal(1,:),(cell_list(i)%centre(:,1) &
!                -vertices(base(1,1),:)))/NORM2(normal(1,:))
!                !*******************************************************
!                !NOTE that norm2 here is not using the built-in fortran funmction
!                !instead it is using norm2 function defined in global
!                !because of that I needed to take sqrt of the norm2
!                !************************************************************
!
!                base(2,1) = Element(e,1)
!                base(2,2) = Element(e,5)
!                base(2,3) = Element(e,6)
!                base(2,4) = Element(e,2)
!                normal(2,:) = cross(&
!                (vertices(base(2,2),:)-vertices(base(2,1),:)), &
!                (vertices(base(2,3),:)-vertices(base(2,1),:)))
!                dP(2) = dot_product(normal(2,:),(cell_list(i)%centre(:,1) &
!                -vertices(base(2,1),:)))/NORM2(normal(2,:))
!
!                base(3,1) = Element(e,4)
!                base(3,2) = Element(e,3)
!                base(3,3) = Element(e,7)
!                base(3,4) = Element(e,8)
!                normal(3,:) = cross(&
!                (vertices(base(3,2),:)-vertices(base(3,1),:)), &
!                (vertices(base(3,3),:)-vertices(base(3,1),:)))
!                dP(3) = dot_product(normal(3,:),(cell_list(i)%centre(:,1) &
!                -vertices(base(3,1),:)))/NORM2(normal(3,:))
!
!                base(4,1) = Element(e,1)
!                base(4,2) = Element(e,2)
!                base(4,3) = Element(e,3)
!                base(4,4) = Element(e,4)
!                normal(4,:) = cross(&
!                (vertices(base(4,2),:)-vertices(base(4,1),:)), &
!                (vertices(base(4,3),:)-vertices(base(4,1),:)))
!                dP(4) = dot_product(normal(4,:),(cell_list(i)%centre(:,1) &
!                -vertices(base(4,1),:)))/NORM2(normal(4,:))
!
!                base(5,1) = Element(e,5)
!                base(5,2) = Element(e,8)
!                base(5,3) = Element(e,7)
!                base(5,4) = Element(e,6)
!                normal(5,:) = cross(&
!                (vertices(base(5,2),:)-vertices(base(5,1),:)), &
!                (vertices(base(5,3),:)-vertices(base(5,1),:)))
!                dP(5) = dot_product(normal(5,:),(cell_list(i)%centre(:,1) &
!                -vertices(base(5,1),:)))/NORM2(normal(5,:))
!
!                base(6,1) = Element(e,2)
!                base(6,2) = Element(e,6)
!                base(6,3) = Element(e,7)
!                base(6,4) = Element(e,3)
!                normal(6,:) = cross(&
!                (vertices(base(6,2),:)-vertices(base(6,1),:)), &
!                (vertices(base(6,3),:)-vertices(base(6,1),:)))
!                dP(6) = dot_product(normal(6,:),(cell_list(i)%centre(:,1) &
!                -vertices(base(6,1),:)))/NORM2(normal(6,:))
!
!                !if ((i==1 .or. i==2) .and. (e==101)) then
!                    !print *, "Element", Element(e,:)
!                    !print *, "base", base(1,:)
!                    !print *, "normal", normal(1,:)
!                    !print *, "!!", dp(1),dp(2),dp(3),dp(4),dp(5),dp(6)
!                !elseif (i==3 .and. e==111) then
!                !   print *, "!!", dp(1),dp(2),dp(3),dp(4),dp(5),dp(6)
!                !endif
!
!              if (NINT(dP(1))>=0 .AND. NINT(dP(2))>=0 .AND. NINT(dP(3))>=0 &
!                .AND. NINT(dP(4))>=0 .AND. NINT(dP(5))>=0 .AND.  NINT(dP(6))>=0) then
!                  cell_elem(i)=e
!              endif
!            enddo
!        endif
!    enddo
!enddo
!
!
!end subroutine find_cell_mesh


subroutine read_real_array(file_name, NoRow, Nocol,n2f)

character * ( * ) :: file_name
integer :: input, status1,size1
integer :: input_stat,I, ii,Nocol,NoRow
integer :: line_no
real (dp), allocatable :: n2f(:,:)

IF (ALLOCATED (n2f)) DEALLOCATE (n2f)

    !call get_unit ( input )
    input=1
    !write (*,*) "filename:", fileplace//file_name
    open ( unit = input, file = trim(file_name), status = 'old', &
     iostat = input_stat )
    !write (*,*) "iostat", input_stat
    if ( input_stat .ne. 0 ) then
        write ( *, '(a)' ) ''
        write ( *, '(a)' ) 'REAL_READ_READ - Fatal error!'
        write ( *, '(a)' ) '  Could not open input file:'// trim ( file_name )
        stop 1
    end if
    allocate(n2f(NoRow,Nocol))

    do ii=1,NoRow
        read(input, *) n2f(ii,:)
        !print *, "n2f", n2f(ii,:)
    enddo
    !print *, "size n2f", size(n2f,1), size(n2f,2)
close (input)
end subroutine read_real_array


subroutine read_integer_array(file_name, NoRow, Nocol,n2f)

character * ( * ) :: file_name
integer :: input, status1,size1
integer :: input_stat,I, ii,Nocol,NoRow
integer :: line_no
integer, allocatable :: n2f(:,:)

IF (ALLOCATED (n2f)) DEALLOCATE (n2f)

    !call get_unit ( input )
    input=1
    !write (*,*) "filename:", fileplace//file_name
    open ( unit = input, file = trim(file_name), status = 'old', &
     iostat = input_stat )
    !write (*,*) "iostat", input_stat
    if ( input_stat .ne. 0 ) then
        write ( *, '(a)' ) ''
        write ( *, '(a)' ) 'REAL_READ_READ - Fatal error!'
        write ( *, '(a)' ) '  Could not open input file:'// trim ( file_name )
        stop 1
    end if
    allocate(n2f(NoRow,Nocol))

    do ii=1,NoRow
        read(input, *) n2f(ii,:)
        !print *, "n2f", n2f(ii,:)
    enddo
    !print *, "size n2f", size(n2f,1), size(n2f,2)
close (input)
end subroutine read_integer_array
!
!##################3
!

subroutine create_sampling_grid()

  use arrays, only: dp,sampling_grid,sampling_nodes,sampling_elems

    integer :: i, j, k, counter
    real(dp) :: x, y, z, site(3)
    integer :: nel_x ,nel_y , nel_z
    integer :: nnp
    real (dp) :: x_min, y_min, z_min, x_max, y_max, z_max
    integer :: EN1, EN2, EN3, EN4, EN5, EN6, EN7, EN8
    integer :: NodeNo(8), ElCount


    sampling_grid%nsd  = 3                             ! number of space dimensions
    sampling_grid%nel_z = 10                       ! number of elements on z axis
    sampling_grid%nel_y = 4                       ! number of elements on y axis
    sampling_grid%nel_x = 4                      ! number of elements on x axis
    sampling_grid%nel  = sampling_grid%nel_y*sampling_grid%nel_x*sampling_grid%nel_z                 ! number of elements
    sampling_grid%nnp  = (sampling_grid%nel_z+1)*(sampling_grid%nel_y+1)*(sampling_grid%nel_x+1)     ! number of nodal nodes
    sampling_grid%nen  = 8                             ! number of element nodes
    sampling_grid%ndof = 1                            ! degrees-of-freedom per node
    sampling_grid%neq  = sampling_grid%nnp*sampling_grid%ndof                ! number of equations
    sampling_grid%x_min = -200.0_dp             !min x coordinates - in mm
    sampling_grid%x_max = 200.0_dp          !max x coordinates - in mm
    sampling_grid%y_min = -200.0_dp             !min y coordinates - in mm
    sampling_grid%y_max = 200.0_dp            !max y coordinates - in mm
    sampling_grid%z_min = 0             !min z coordinates - in mm
    sampling_grid%z_max = 1000            !max z coordinates - in mm
    sampling_grid%x_width = (sampling_grid%x_max-sampling_grid%x_min)/sampling_grid%nel_x
    sampling_grid%y_width = (sampling_grid%y_max-sampling_grid%y_min)/sampling_grid%nel_y
    sampling_grid%z_width = (sampling_grid%z_max-sampling_grid%z_min)/sampling_grid%nel_z
    sampling_grid%volume = sampling_grid%x_width*sampling_grid%y_width*sampling_grid%z_width

    counter=0
    allocate(sampling_nodes(sampling_grid%nnp))

    do k= 0,sampling_grid%nel_z
        z=(sampling_grid%z_max - sampling_grid%z_min)/sampling_grid%nel_z*k + sampling_grid%z_min
        do j= 0,sampling_grid%nel_y
            y=(sampling_grid%y_max - sampling_grid%y_min)/sampling_grid%nel_y*j + sampling_grid%y_min
            do i= 0,sampling_grid%nel_x
                x=(sampling_grid%x_max - sampling_grid%x_min)/sampling_grid%nel_x*i + sampling_grid%x_min
                counter=counter+1
                site=(/x,y,z/)
                sampling_nodes(counter)%ID = counter
                sampling_nodes(counter)%coordinates=site

                print *, "Coordinates!=", counter, sampling_grid%nnp, site
            enddo
        enddo
    enddo

    write(*,*) 'finished with the nodes'
    allocate(sampling_elems(sampling_grid%nel))
    ElCount=0
    do k =  1,sampling_grid%nel_z
        do j = 1,sampling_grid%nel_y
            do i = 1,sampling_grid%nel_x
                EN1 = i+(sampling_grid%nel_x+1)*(j-1)+(sampling_grid%nel_x+1)*(sampling_grid%nel_y+1)*(k-1)
                EN2 = EN1+1
                EN3 = EN1+sampling_grid%nel_x+1
                EN4 = EN3+1
                EN5 = EN1+(sampling_grid%nel_x+1)*(sampling_grid%nel_y+1)
                EN6 = EN2+(sampling_grid%nel_x+1)*(sampling_grid%nel_y+1)
                EN7 = EN3+(sampling_grid%nel_x+1)*(sampling_grid%nel_y+1)
                EN8 = EN4+(sampling_grid%nel_x+1)*(sampling_grid%nel_y+1)

                ElCount = ElCount+1
                NodeNo=(/EN1, EN2, EN3, EN4, EN5, EN6, EN7, EN8/)
                sampling_elems(ElCount)%ID = ElCount
                sampling_elems(ElCount)%nodes= NodeNo
                print *, "elements Node Numbers!=", sampling_elems(ElCount)%nodes
            enddo
        enddo
    enddo

    call mesh_cylinder_volume
    call cellcount
end subroutine create_sampling_grid

!
!!!!!!!!!!!!!!!!!!!!!!!!!
!
subroutine mesh_cylinder_volume()

    use arrays, only: dp,sampling_grid,sampling_nodes,sampling_elems

    use other_consts, only:PI

    integer :: nel_x ,nel_y , nel_z, total_elems
    integer :: indices(8)
    integer :: i , j
    real (dp) :: x(8),x1, x0,y(8), y1, y0,z(8), z1, z0
    real (dp) :: angle_between, segmentprop, triangle_height
    real (dp) :: chord_length, segment_area,corner_triangle_area
    real (dp) :: approx_circle_area, volume

    total_elems=sampling_grid%nel

    do i=1,total_elems
       sampling_elems(i)%cylinder_volume = 1.0_dp
       indices(:)=sampling_elems(i)%nodes
       do j=1, 8
            x(j)=abs(sampling_nodes(indices(j))%coordinates(1))
            y(j)=abs(sampling_nodes(indices(j))%coordinates(2))
            z(j)=abs(sampling_nodes(indices(j))%coordinates(3))
       enddo
       x1=max(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8))
       x0=min(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8))
       y1=max(y(1),y(2),y(3),y(4),y(5),y(6),y(7),y(8))
       y0=min(y(1),y(2),y(3),y(4),y(5),y(6),y(7),y(8))
       z1=max(z(1),z(2),z(3),z(4),z(5),z(6),z(7),z(8))
       z0=min(z(1),z(2),z(3),z(4),z(5),z(6),z(7),z(8))
       if ((sqrt(x1**2+y1**2)) < 200.0_dp) then
           sampling_elems(i)%cylinder_volume = 1.0_dp
       elseif ((sqrt(x1**2+y0**2)>=200.0_dp).AND.(sqrt(x0**2+y1**2)>=200.0_dp)) then
           angle_between=atan2(sqrt(200.0_dp**2-x0**2),x0)-atan2(y0,sqrt(200.0_dp**2-y0**2))
           segmentprop=abs((angle_between)/(2*Pi))
           triangle_height=200.0_dp*cos(0.5*(angle_between))
           chord_length=200.0_dp*2*sin(0.5*(angle_between))
           segment_area=Pi*200.0_dp**2*segmentprop-0.5_dp*triangle_height*chord_length
           corner_triangle_area=0.5_dp*((sqrt(200.0_dp**2-y0**2)-x0)*(sqrt(200.0_dp**2-x0**2)-y0))
        approx_circle_area=corner_triangle_area+segment_area

        sampling_elems(i)%cylinder_volume = approx_circle_area*(z1-z0)/sampling_grid%volume
        if(sampling_elems(i)%cylinder_volume.gt.1.0_dp)then
        write(*,*) 'here', x1,x0,y1,y0,z1,z0,corner_triangle_area,segment_area,sampling_grid%volume
        endif

       elseif ((sqrt(x1**2+y0**2)<200.0_dp).AND.(sqrt(x0**2+y0**2)<200.0_dp) ) then
            angle_between=atan2(sqrt(200.0_dp**2-x0**2),x0)-atan2(sqrt(200.0_dp**2-x1**2),x1)
           segmentprop=abs((angle_between)/(2*Pi))
           triangle_height=200.0_dp*cos(0.5*(angle_between))
           chord_length=200.0_dp*2*sin(0.5*(angle_between))
           segment_area=Pi*200.0_dp**2*segmentprop-0.5*triangle_height*chord_length
           approx_circle_area=segment_area+(min(sqrt(200.0_dp**2-x0**2), &
           sqrt(200.0_dp**2-x1**2))-y0)*(x1-x0)+abs(sqrt(200.0_dp**2-x0**2)-sqrt(200.0_dp**2-x1**2))*(x1-x0)*0.5
           sampling_elems(i)%cylinder_volume =approx_circle_area*(z1-z0)/sampling_grid%volume
       elseif ((sqrt(x0**2+y1**2)<=200.0_dp).AND.(sqrt(x0**2+y1**2)<=200.0_dp) ) then
            angle_between=atan2(sqrt(200.0_dp**2-y0**2),y0)-atan2(sqrt(200.0_dp**2-y1**2),y1)
           segmentprop=abs((angle_between)/(2*Pi))
           triangle_height=200.0_dp*cos(0.5*(angle_between))
           chord_length=200.0_dp*2*sin(0.5*(angle_between))
           segment_area=Pi*200.0_dp**2*segmentprop-0.5*triangle_height*chord_length
           approx_circle_area=segment_area+(min(sqrt(200.0_dp**2-y0**2), &
           sqrt(200.0_dp**2-y1**2))-x0)*(y1-y0)+abs(sqrt(200.0_dp**2-y0**2)-sqrt(200.0_dp**2-y1**2))*(y1-y0)*0.5
           sampling_elems(i)%cylinder_volume = approx_circle_area*(z1-z0)/sampling_grid%volume
       endif
       write(*,*) i,sampling_elems(i)%cylinder_volume
    enddo


end  subroutine


subroutine cellcount()
    use arrays, only: dp,sampling_grid,sampling_nodes,sampling_elems,num_cells, cell_list

    use other_consts, only:PI


!type(element_type) :: EPm
!real (REAL_KIND) :: EpmCoor(:,:), Cell(:,:)
!    real (dp) :: sum_PyVol,Pyramids_volume
    integer :: i
    integer :: xelem_num(4),yelem_num,zelem_num,nelem
    real(dp) :: xelem_real


    do i = 1,sampling_grid%nel
      sampling_elems(nelem)%cell_cnt = 0
    enddo


    do i=1,num_cells

      xelem_real = (cell_list(i)%centre(1,1)- sampling_grid%x_min) / sampling_grid%x_width
      xelem_num (1:4) = ceiling(xelem_real)

      write(*,*) xelem_num
    !    xelem_num = xelem_num - 1
      yelem_num = ceiling((cell_list(i)%centre(2,1) - sampling_grid%y_min) / sampling_grid%y_width)
      write(*,*) yelem_num
    !if yelem_num == int(nelem_y):
    !    yelem_num = yelem_num - 1
      zelem_num = ceiling((cell_list(i)%centre(3,1) - sampling_grid%z_min) / sampling_grid%z_width)
      write(*,*) zelem_num
    !if zelem_num == int(nelem_z):
    !    zelem_num = zelem_num - 1
      nelem = xelem_num(1) + yelem_num * sampling_grid%nel_x + zelem_num * (sampling_grid%nel_x * sampling_grid%nel_y)  ! this is the element where the point/node located
      write(*,*) nelem

      sampling_elems(nelem)%cell_cnt = sampling_elems(nelem)%cell_cnt + 1
   end do
end subroutine


!-----------------------------------------------
!-----------------------------------------------
function diag(matrix,rank) result(diagmatrix)
integer :: rank,i
integer :: matrix(rank)
real(dp) :: diagmatrix(rank,rank)


diagmatrix(:,:) = 0
!print *, "zeros=", diagmatrix(:,:)
do i = 1, rank
   diagmatrix(i,i) = matrix(i)
end do
end function

!-----------------------------------------
FUNCTION cross(a, b)
  real (dp), DIMENSION(3) :: cross
  real (dp), DIMENSION(3), INTENT(IN) :: a, b

  cross(1) = a(2) * b(3) - a(3) * b(2)
  cross(2) = a(3) * b(1) - a(1) * b(3)
  cross(3) = a(1) * b(2) - a(2) * b(1)
END FUNCTION cross


end module mix_fem
