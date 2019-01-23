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
  public read_face2e
  public create_sampling_grid
  public compute_body_forces
  public solve_fem
contains

!
!############################################################
!
subroutine assemble_sparse_matrices
    use arrays, only: num_elems,B_MATRIX,all_faces,num_all_faces,element2face,plug_params,&
      elem_field,num_common_face, feA, feIA, feJA
    use diagnostics, only: enter_exit,get_diagnostics_level
    use indices
    use other_consts, only: MAX_FILENAME_LEN
    implicit none
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_ASSEMBLE_SPARSE_MATRICES" :: READ_ASSEMBLE_SPARSE_MATRICES


!integer :: nel
    real(dp) :: B(6,6),C(6)
    real(dp) :: BB(6*6*num_elems-num_common_face)
    real(dp) :: CC(6*num_elems)
    integer :: IB(6*6*num_elems-num_common_face), JB(6*6*num_elems-num_common_face)!, ones(6)
    integer :: IC(6*num_elems), JC(6*num_elems)
!!real(dp), allocatable :: Ak(:,:)
    integer :: ii,j, jj, I(6), signum(6), DiagSig(6,6)
    integer :: counter, Recount1, Ccounter, commonCount
    logical :: ok1
!character * ( 255 ) :: fileplace, Bmatrix_filename

    character(len=60) :: sub_name
    integer :: diagnostics_level

    sub_name = 'assemble_sparse_matrices'
    call enter_exit(sub_name,1)
    call get_diagnostics_level(diagnostics_level)

    if(.not.allocated(feA)) allocate(feA(8*6*num_elems-num_common_face))
    if(.not.allocated(feIA)) allocate(feIA(8*6*num_elems-num_common_face))
    if(.not.allocated(feJA)) allocate(feJA(8*6*num_elems-num_common_face))

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
    B(:,:) = plug_params%blood_viscosity/elem_field(j,ne_cond)* &
        matmul(matmul(DiagSig, B_MATRIX(j)%Inside_B),DiagSig)
    C(:)=matmul(ones(6),diag(signum,size(signum(:))))
    do ii=1,6
        do jj =1,6
            if (B(ii,jj).ne.0.0_dp ) then
                if (counter > 6) then
                    ok1 = .false.
                    Recount1=findequal(IB(1:counter),JB(1:counter),I(ii),I(jj),ok1)
!
                    if ( ok1 ) then
                        BB(Recount1) = B(ii,jj)+BB(Recount1)
                    else
                        counter=counter+1
                        BB(counter) = B(ii,jj)
                        IB(counter) = I(ii)
                        if(IB(counter).eq.0)write(*,*)'else10',counter
                        JB(counter) = I(jj)
                        if(JB(counter).eq.0)write(*,*)'else10',counter

                    endif
                else
                        counter=counter+1
                        BB(counter) = B(ii,jj)
                        IB(counter) = I(ii)
                        if(IB(counter).eq.0)write(*,*)'else1',counter
                        JB(counter) = I(jj)
                        if(JB(counter).eq.0)write(*,*)'else20',counter

                endif
            endif
        enddo
        Ccounter = Ccounter+1
        CC(Ccounter)= C(ii)
        IC(Ccounter)=I(ii)
        if(IB(Ccounter).eq.0)write(*,*)'else2',Ccounter

        JC(Ccounter)=j
        if(JC(Ccounter).eq.0)write(*,*)'else20',Ccounter

      enddo
    enddo

    if(diagnostics_level.gt.1)write(*,*) 'Assembling global stiffness matrices'
    !Global stiffness matrix A
    feA(:) = [BB(:), CC(:), CC(:)]
    !Row identifiers
    feIA(:) = [IB(:), IC(:), JC(:) + num_all_faces*ones(size(JC(:)))]
    !column identifiers
    feJA(:) = [JB(:), JC(:) + num_all_faces*ones(size(JC(:))), IC(:)]

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

      if(diagnostics_level.gt.1)write(*,*) appended_fn,B_MATRIX(ne)%Inside_B
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

    if(diagnostics_level.gt.1)write(*,*) element2face(1,:)
    if(diagnostics_level.gt.1)write(*,*) element2face(num_elems,:)


    call enter_exit(sub_name,2)
end subroutine read_e2face


subroutine read_face2e(filename,filename_len)
   use arrays, only: num_elems,face2element,num_common_face,num_all_faces
    use diagnostics, only: enter_exit,get_diagnostics_level
    use indices
    use other_consts, only: MAX_FILENAME_LEN
    implicit none
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_READ_FACE2E" :: READ_FACE2E

    character(len=MAX_FILENAME_LEN), intent(in) :: filename !Input nodefile
    integer, intent(in) :: filename_len


    character(len=MAX_FILENAME_LEN) :: appended_fn
    character(len = filename_len) ::trimfilename
    integer :: ne,nf

    character(len=60) :: sub_name
    integer :: diagnostics_level

    sub_name = 'read_face2e'
    call enter_exit(sub_name,1)
    call get_diagnostics_level(diagnostics_level)

    trimfilename = filename(1:filename_len)
    allocate(face2element(num_all_faces,6))

    call read_integer_array(adjustl(trimfilename),num_all_faces,6,face2element)

    if(diagnostics_level.gt.1)write(*,*) face2element(1,:)
    if(diagnostics_level.gt.1)write(*,*) face2element(num_all_faces,:)
    num_common_face = 0
    do nf = 1,size(face2element,1)
      if (face2element(nf,6).ne.0) then
        num_common_face = num_common_face + 1
      endif
    enddo

    if(diagnostics_level.gt.1)write(*,*) 'num common face', num_common_face
    call enter_exit(sub_name,2)
end subroutine read_face2e

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
    use diagnostics, only: enter_exit,get_diagnostics_level


character * ( * ) :: file_name
integer :: input, status1,size1
integer :: input_stat,I, ii,Nocol,NoRow
integer :: line_no
real (dp), allocatable :: n2f(:,:)

    character(len=60) :: sub_name
    integer :: diagnostics_level

    sub_name = 'read_real_array'
    call enter_exit(sub_name,1)
    call get_diagnostics_level(diagnostics_level)

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
    use diagnostics, only: enter_exit,get_diagnostics_level


character * ( * ) :: file_name
integer :: input, status1,size1
integer :: input_stat,I, ii,Nocol,NoRow
integer :: line_no
integer, allocatable :: n2f(:,:)

    character(len=60) :: sub_name
    integer :: diagnostics_level

    sub_name = 'read_integer_array'
    call enter_exit(sub_name,1)
    call get_diagnostics_level(diagnostics_level)

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
  use diagnostics, only: enter_exit,get_diagnostics_level


    integer :: i, j, k, counter
    real(dp) :: x, y, z, site(3)
    integer :: nel_x ,nel_y , nel_z
    integer :: nnp
    real (dp) :: x_min, y_min, z_min, x_max, y_max, z_max
    integer :: EN1, EN2, EN3, EN4, EN5, EN6, EN7, EN8
    integer :: NodeNo(8), ElCount

    character(len=60) :: sub_name
    integer :: diagnostics_level

    sub_name = 'create_sampling_grid'
    call enter_exit(sub_name,1)
    call get_diagnostics_level(diagnostics_level)


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

    if(diagnostics_level.gt.1)write(*,*) 'finished with the nodes'
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
    call mesh_node2sample

end subroutine create_sampling_grid

!
!!!!!!!!!!!!!!!!!!!!!!!!!
!
subroutine mesh_cylinder_volume()

    use arrays, only: dp,sampling_grid,sampling_nodes,sampling_elems
    use diagnostics, only: enter_exit,get_diagnostics_level
    use other_consts, only:PI

    integer :: nel_x ,nel_y , nel_z, total_elems
    integer :: indices(8)
    integer :: i , j
    real (dp) :: x(8),x1, x0,y(8), y1, y0,z(8), z1, z0
    real (dp) :: angle_between, segmentprop, triangle_height
    real (dp) :: chord_length, segment_area,corner_triangle_area
    real (dp) :: approx_circle_area, volume

    character(len=60) :: sub_name
    integer :: diagnostics_level

    sub_name = 'mesh_cylinder_volume'
    call enter_exit(sub_name,1)
    call get_diagnostics_level(diagnostics_level)

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
        if(diagnostics_level.gt.1)write(*,*) 'here', x1,x0,y1,y0,z1,z0,corner_triangle_area,segment_area,sampling_grid%volume
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
       if(diagnostics_level.gt.1)write(*,*) i,sampling_elems(i)%cylinder_volume
    enddo


end  subroutine

!
!########################################
!

subroutine mesh_node2sample()
    use arrays, only: dp,sampling_grid,sampling_nodes,sampling_elems,num_nodes,node_xyz,&
      elem_3d,num_elems,elem_field
    use diagnostics, only: enter_exit,get_diagnostics_level
    use indices, only: ne_cond
    use other_consts, only:PI


    integer :: i,j,k,l,ne
    integer :: xelem_num,yelem_num,zelem_num,nelem
    real(dp), allocatable :: node_cond(:),node_cond_sm(:), nelem_node(:)

    character(len=60) :: sub_name
    integer :: diagnostics_level

    sub_name = 'mesh_node2sample'
    call enter_exit(sub_name,1)
    call get_diagnostics_level(diagnostics_level)

    do i = 1,sampling_grid%nel
      sampling_elems(i)%mesh_node_cnt = 0
    enddo
    allocate(node_cond(num_nodes))
    allocate(node_cond_sm(num_nodes))
    allocate(nelem_node(num_nodes))
    nelem_node = 0.0_dp
    node_cond_sm = 0.0_dp

    do i=1,num_nodes

      xelem_num = ceiling((node_xyz(i,1) -  sampling_grid%x_min) / sampling_grid%x_width)
      if(xelem_num.eq.0) xelem_num =1
      yelem_num = ceiling((node_xyz(i,2) - sampling_grid%y_min) / sampling_grid%y_width)
      if(yelem_num.eq.0) yelem_num =1
      zelem_num = ceiling((node_xyz(i,3) - sampling_grid%z_min) / sampling_grid%z_width)
      if(zelem_num.eq.0) zelem_num =1

      nelem = xelem_num + (yelem_num-1) * (sampling_grid%nel_x) &
              + (zelem_num-1) * (sampling_grid%nel_x * sampling_grid%nel_y)  ! this is the element where the point/node located
      sampling_elems(nelem)%mesh_node_cnt = sampling_elems(nelem)%mesh_node_cnt + 1

      node_cond(i) = sampling_elems(nelem)%k_conduct

      if(diagnostics_level.gt.1)write(*,*) node_cond(i)
   enddo

   do i = 1,sampling_grid%nel
      if(diagnostics_level.gt.1)write(*,*) i,sampling_elems(i)%mesh_node_cnt
   enddo

   do ne = 1,num_elems
     elem_field(ne,ne_cond) = 0.0_dp
     do i = 1,8
       elem_field(ne,ne_cond) = elem_field(ne,ne_cond) + node_cond(elem_3d(ne,i))
     enddo
     elem_field(ne,ne_cond)  = elem_field(ne,ne_cond) / 8.0_dp

     !smoothing step
     !Give back to nodes
     if(diagnostics_level.gt.1)write(*,*) ne, elem_field(ne,ne_cond)
     do i =1,8
       nelem_node(elem_3d(ne,i)) = nelem_node(elem_3d(ne,i)) + 1.0_dp
       node_cond_sm(elem_3d(ne,i)) = node_cond_sm(elem_3d(ne,i)) + elem_field(ne,ne_cond)
     enddo
   enddo

   do i = 1,num_nodes
     node_cond_sm(i) = node_cond_sm(i)/nelem_node(i)
   enddo
   !and back to elems
   do ne = 1,num_elems
     elem_field(ne,ne_cond) = 0.0_dp
     do i = 1,8
       elem_field(ne,ne_cond) = elem_field(ne,ne_cond) + node_cond_sm(elem_3d(ne,i))
     enddo
     elem_field(ne,ne_cond)  = elem_field(ne,ne_cond) / 8.0_dp
     if(diagnostics_level.gt.1)write(*,*) ne, elem_field(ne,ne_cond)
   enddo

   deallocate(node_cond)
   deallocate(nelem_node)
   deallocate(node_cond_sm)
end subroutine
!
!##################################
!

subroutine cellcount()
    use arrays, only: dp,sampling_grid,sampling_nodes,sampling_elems,num_cells, cell_list,plug_params
    use diagnostics, only: enter_exit,get_diagnostics_level
    use other_consts, only:PI


    integer :: i,j,k,l
    integer :: xelem_num(2),yelem_num(2),zelem_num(2),nelem
    real(dp) :: xelem_real,yelem_real,zelem_real

    character(len=60) :: sub_name
    integer :: diagnostics_level

    sub_name = 'cellcount'
    call enter_exit(sub_name,1)
    call get_diagnostics_level(diagnostics_level)


    do i = 1,sampling_grid%nel
      sampling_elems(i)%cell_cnt = 0.0_dp
    enddo


    do i=1,num_cells
      xelem_real = (cell_list(i)%centre(1,1)- sampling_grid%x_min) / sampling_grid%x_width
      if(xelem_real.eq.nint(xelem_real))then
        xelem_num (1) =ceiling(xelem_real)
        xelem_num(2) = xelem_num(1)-1
      else
        xelem_num (1:2) = ceiling(xelem_real)
      endif

      yelem_real = ((cell_list(i)%centre(2,1) - sampling_grid%y_min) / sampling_grid%y_width)
      if(yelem_real.eq.nint(yelem_real))then
        yelem_num (1) =ceiling(yelem_real)
        yelem_num(2) = yelem_num(1)-1
      else
        yelem_num (1:2) = ceiling(yelem_real)
      endif

      zelem_real = ((cell_list(i)%centre(3,1) - sampling_grid%z_min) / sampling_grid%z_width)
      if(zelem_real.eq.nint(zelem_real))then
        zelem_num (1) =ceiling(zelem_real)
        zelem_num(2) = zelem_num(1)-1
      else
        zelem_num (1:2) = ceiling(zelem_real)
      endif

      do l = 1,2
        do j = 1,2
          do k = 1,2
            nelem = xelem_num(l) + (yelem_num(j)-1) * sampling_grid%nel_x &
              + (zelem_num(k)-1) * (sampling_grid%nel_x * sampling_grid%nel_y)  ! this is the element where the point/node located
            sampling_elems(nelem)%cell_cnt = sampling_elems(nelem)%cell_cnt + 1.0_dp/8.0_dp
          enddo
        enddo
      enddo
   end do
   do i = 1,sampling_grid%nel
      sampling_elems(i)%volume_fraction = sampling_elems(i)%cell_cnt*4.0_dp*PI*plug_params%Raverage**3.0_dp&
        /(sampling_elems(i)%cylinder_volume*sampling_grid%volume*3.0_dp)
      if(sampling_elems(i)%volume_fraction.gt.1.0_dp)then
        sampling_elems(i)%volume_fraction = 0.98
      endif
      if(sampling_elems(i)%volume_fraction.eq.0.0_dp)then
        sampling_elems(i)%k_conduct = plug_params%k_empty
      else
       sampling_elems(i)%k_conduct = (2.0_dp*plug_params%Raverage)**2.0_dp&
          /180.0_dp*(1.0_dp-sampling_elems(i)%volume_fraction)**3.0_dp &
            /(sampling_elems(i)%volume_fraction**2.0_dp) ! in um^2
      endif
      if(sampling_elems(i)%k_conduct.gt.plug_params%k_empty)then
        sampling_elems(i)%k_conduct = plug_params%k_empty
      endif
      if(diagnostics_level.gt.1)then
        write(*,*) i,sampling_elems(i)%cell_cnt, sampling_elems(i)%volume_fraction,sampling_elems(i)%k_conduct
      endif
   enddo
end subroutine

!
!#######################################
!
subroutine compute_body_forces(inletPressure,outletPressure)
    use arrays, only: num_all_faces,num_elems,feAx,feBodyForce,feQ_P,feFreeFaces,num_wall_faces,&
       face_info, face_stat,feIA, feJA,feA
    use diagnostics, only: enter_exit,get_diagnostics_level
    use indices
    use other_consts, only: MAX_FILENAME_LEN
    implicit none
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_COMPUTE_BODY_FORCE" :: READ_COMPUTE_BODY_FORCE
    real(dp), intent(in) :: inletPressure, outletPressure

    integer :: i, tmp(num_all_faces+num_elems),cnt_wall
    character(len=60) :: sub_name
    integer :: diagnostics_level

    sub_name = 'compute_body_forces'
    call enter_exit(sub_name,1)
    call get_diagnostics_level(diagnostics_level)

    if(.not.allocated(feAx))allocate(feAx(num_all_faces+num_elems))
    if(.not.allocated(feBodyForce))allocate(feBodyForce(num_all_faces+num_elems))
    if(.not.allocated(feQ_P))allocate(feQ_P(num_all_faces+num_elems))
    if(.not.allocated(feFreeFaces))allocate(feFreeFaces(num_all_faces+num_elems-num_wall_faces))

    feBodyForce = 0.0_dp
    feQ_P = 0.0_dp
    tmp(:) = 0
    cnt_wall = 0
    do i= 1, num_all_faces
      if(face_info(i).eq.face_stat%INLET)then
        feBodyForce(i) = outletPressure!inletPressure
      elseif(face_info(i).eq.face_stat%OUTLET)then
        feBodyForce(i) = inletPressure!outletPressure
      elseif(face_info(i).eq.face_stat%WALL)then
        feQ_P(i) = 0.0_dp
        cnt_wall = cnt_wall+1
        tmp(i) = 1
      endif
    enddo

    !free faces
    write(*,*) size(tmp),size(feFreeFaces),cnt_wall,num_wall_faces
    feFreeFaces(:) = findNONzeros(tmp(:))

    call ax_st(size(feBodyForce),size(feA), feIA, feJA, &
    feA,feQ_P, feAx ) !computes A*x, A is a sparse matrix
    feBodyForce(:) = feBodyForce(:) - feAx(:)

    call enter_exit(sub_name,2)
end subroutine compute_body_forces

!
!#######################################
!
subroutine solve_fem
    use arrays, only: dp,feIA,feFreeFaces,feJA,feBodyForce,feQ_P,feA
    use diagnostics, only: enter_exit,get_diagnostics_level
    use indices
    use sparsekit
    implicit none
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SOLVE_FEM" :: READ_SOLVE_FEM

    integer :: count1,i,j,k,n,sizeFFace,rowcount,columncount,p
    real(dp) :: t1,t2


    integer, allocatable :: Ap_save(:), Ai_save(:), NewIA(:), NewJA(:),row_ind(:), col_ptr(:)
    integer, allocatable :: jao(:),iao(:)
    real(dp), allocatable :: Ax_save(:),newBForce(:), newQP(:),newA(:),ao(:),acsc(:),solution(:),residual(:)

    character(len=60) :: sub_name
    integer :: diagnostics_level

    sub_name = 'solve_fem'
    call enter_exit(sub_name,1)
    call get_diagnostics_level(diagnostics_level)

    count1=0
    do i = 1,size(feFreeFaces)
      do j= 1, size(feIA)
         if (feFreeFaces(i).eq.feIA(j)) then
            do k = 1, size(feFreeFaces)
                if (feFreeFaces(k).eq.feJA(j)) then
                count1 = count1+1
                endif
            enddo
         endif
      enddo
  enddo
  n = count1

  sizeFFace = size(feFreeFaces)
   write(*,*) n,sizeFFace
  allocate(Ap_save(sizeFFace+1))
  allocate(Ai_save(n))
  allocate(Ax_save(n))
  allocate(NewIA(n))
  allocate(newJA(n))
  allocate(newBForce(sizeFFace))
  allocate( newQP(sizeFFace))
  allocate(newA(n))
  allocate(ao(n))
  allocate(jao(n))
  allocate(iao(n))
  allocate(acsc(n))
  allocate(row_ind(n))
  allocate(col_ptr(n))
  allocate(solution(sizeFFace))
  allocate(residual(sizeFFace))

  count1 = 0
  rowcount=0
  do i = 1, sizeFFace
    rowcount=rowcount+1
    columncount=0
    do j= 1, size(feIA)
        if (feFreeFaces(i)==feIA(j)) then
            do k = 1, sizeFFace
                if (feJA(j)==feFreeFaces(k)) then
                    count1 = count1+1
                    NewA(count1) = feA(j)
                    NewIA(count1) = rowcount !IA(j)
                    NewJA(count1) = k !JA(j)
                    !write (out_unit,*) NewA(count1),NewIA(count1),NewJA(count1)
                endif
            enddo
        endif
    enddo
    newBForce(i) = feBodyForce(feFreeFaces(i))
    newQP(i) = feQ_P(feFreeFaces(i))
  enddo


    call coocsr (sizeFFace,size(NEWA(:)),NewA,NEWIA,NEWJA,ao,jao,iao )

    call rearrange_cr (sizeFFace, iao, jao, ao )

    call csrcsc (sizeFFace, 1, 1, ao, jao, iao, acsc, row_ind, col_ptr )

    !----------------------------------------------------------------
    ! convert from 1-based to 0-based
    !----------------------------------------------------------------
    do j = 1, sizeFFace+1
        col_ptr(j) = col_ptr(j) - 1
    enddo
    do p = 1, size(NEWA(:))
        row_ind(p) = row_ind(p) - 1
    enddo


     !========== start slu solver=============================
     solution= newBForce
     Ap_save = col_ptr
     Ai_save = row_ind
     Ax_save = acsc
     call cpu_time(t1)

    call c_fortran_cgssv(4, sizeFFace, size(NEWA(:)), acsc, row_ind, col_ptr, solution);

    write(*,'(a,f8.1)') 'Solve time (cpu_time/nprocs) (sec): ',(t2-t1)/4
!
!!! print the residual.
    write(*,*) 'After solve, residual:'
    call resid (sizeFFace, Ap_save, Ai_save, Ax_save, solution, newBForce, residual)
!!========== end slu solver=============================
!
    newQP = solution
    !map back to soln array

    do i = 1, sizeFFace
        feQ_P(feFreeFaces(i))= newQP(i)
    enddo

  deallocate(Ap_save)
  deallocate(Ai_save)
  deallocate(Ax_save)
  deallocate(NewIA)
  deallocate(newJA)
  deallocate(newBForce)
  deallocate( newQP)
  deallocate(newA)
  deallocate(ao)
  deallocate(jao)
  deallocate(iao)
  deallocate(acsc)
  deallocate(row_ind)
  deallocate(col_ptr)
  deallocate(solution)
  deallocate(residual)

    call enter_exit(sub_name,2)
end subroutine solve_fem



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

!-----------------------------------------------
!
!-----------------------------------------------
function ones(Msize)

integer :: Msize, i, ones(Msize)
    do i= 1, Msize
        ones(i) = 1
    enddo
!print *, "ones=", ones(:)
end function

!-----------------------------------------------
!
!-----------------------------------------------
function findequal(IB, JB,Row,Col,ok)

integer :: IB(:), JB(:), Row,Col , findequal, ii
logical :: ok
 ok = .false.
do ii=1,size(IB)
    if (IB(ii)==Row .AND. JB(ii)==Col) then
        findequal=ii
        !print *, "recount", findequal
        ok = .true.
        return
    endif
enddo

end function

!-----------------------------------------------
!
!-----------------------------------------------
function findNONzeros(TMP)
use arrays, only: num_all_faces,num_elems,num_wall_faces
integer :: TMP(num_all_faces+num_elems), findNONzeros(num_all_faces+num_elems-num_wall_faces)
integer :: i, count1

count1=0
do i = 1, size(TMP(:))
    if (TMP(i)==0) then
    count1=count1+1
        findNONzeros(count1) = i
    endif
enddo



end function

subroutine ax_st ( n, nz_num, ia, ja, a, x, w )

!*****************************************************************************80
!
!! AX_ST computes A*x for a matrix stored in sparset triplet form.
!
!  Discussion:
!
!    The matrix A is assumed to be sparse.  To save on storage, only
!    the nonzero entries of A are stored.  For instance, the K-th nonzero
!    entry in the matrix is stored by:
!
!      A(K) = value of entry,
!      IA(K) = row of entry,
!      JA(K) = column of entry.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2006
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
!    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
!    Charles Romine, Henk van der Vorst,
!    Templates for the Solution of Linear Systems:
!    Building Blocks for Iterative Methods,
!    SIAM, 1994.
!    ISBN: 0898714710,
!    LC: QA297.8.T45.
!
!    Tim Kelley,
!    Iterative Methods for Linear and Nonlinear Equations,
!    SIAM, 2004,
!    ISBN: 0898713528,
!    LC: QA297.8.K45.
!
!    Yousef Saad,
!    Iterative Methods for Sparse Linear Systems,
!    Second Edition,
!    SIAM, 2003,
!    ISBN: 0898715342,
!    LC: QA188.S17.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the system.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzeros.
!
!    Input, integer ( kind = 4 ) IA(NZ_NUM), JA(NZ_NUM), the row and column
!    indices of the matrix values.
!
!    Input, real ( kind = 8 ) A(NZ_NUM), the matrix values.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) W(N), the value of A*X.
!
  use arrays, only: dp
  use diagnostics, only: enter_exit,get_diagnostics_level

  implicit none

  integer :: n
  integer :: nz_num

  real (dp) :: a(nz_num)
  integer :: i
  integer :: ia(nz_num)
  integer :: j
  integer :: ja(nz_num)
  integer :: k
  real (dp) :: w(n)
  real (dp) :: x(n)

  w(1:n) = 0.0D+00
  write(*,*) nz_num
  do k = 1, nz_num
    i = ia(k)
    j = ja(k)
    w(i) = w(i) + a(k) * x(j)
  end do

  return
end subroutine

!-------------------------------------------------------
!*****************************************************************************
 subroutine rearrange_cr ( n, ia, ja, a )
!!! REARRANGE_CR sorts a sparse compressed row matrix.
    !    This routine guarantees that the entries in the CR matrix
    !    are properly sorted.
    !
    !    After the sorting, the entries of the matrix are rearranged in such
    !    a way that the entries of each column are listed in ascending order
    !    of their column values.
    !    Input, integer :: N, the order of the system.
    !
    !    Input, integer :: NZ_NUM, the number of nonzeros.
    !
    !    Input, integer :: IA(N+1), the compressed row indices.
    !
    !    Input/output, integer :: JA(NZ_NUM), the column indices.
    !    On output, these may have been rearranged by the sorting.
    !
    !    Input/output, real (REAL_KIND) :: A(NZ_NUM), the matrix values.  On output,
    !    the matrix values may have been moved somewhat because of the sorting.
    !
    use arrays, only: dp
    implicit none

    integer :: n
    integer :: ia(*) !ia(n+1)
    integer :: ja(*) !ja(nz_num)
    real (dp) :: a(*) !a(nz_num)

    integer :: i
    integer :: i4temp
    integer :: k
    integer :: l
    real (dp) :: r8temp

    do i = 1, n

       do k = ia(i), ia(i+1) - 2
          do l = k + 1, ia(i+1) - 1

             if ( ja(l) < ja(k) ) then
                i4temp = ja(l)
                ja(l)  = ja(k)
                ja(k)  = i4temp

                r8temp = a(l)
                a(l)   = a(k)
                a(k)   = r8temp
             end if

          end do
       end do

    end do

    return
 end subroutine rearrange_cr

 !== resid ==============================================================
! Compute the residual, r = Ax-b, its max-norm, and print the max-norm
! Note that A is zero-based.
!=======================================================================
subroutine resid (n, Ap, Ai, Ax, x, b, r)
  integer :: n, Ap(:), Ai(:), j, i, p
  double precision :: Ax(:), x(:), b(:), r(:), rmax, aij

  r(1:n) = -b(1:n)

  do j = 1,n
    do p = Ap(j) + 1, Ap(j+1)
        i = Ai(p) + 1
        aij = Ax(p)
        r(i) = r(i) + aij * x(j)
    enddo
  enddo
  rmax = 0
  do i = 1, n
    rmax = max(rmax, r(i))
  enddo

  write(*,*) 'norm (A*x-b): ', rmax
end subroutine


end module mix_fem
