module geometry
!*Description:* This module handles all geometry read/write/generation.
!
  use other_consts
  implicit none

  !Module parameters

  !Module types

  !Module variables

  !Interfaces
  private
  public add_matching_mesh
  public append_units
  public calc_capillary_unit_length
  public define_1d_elements
  public define_node_geometry
  public define_rad_from_geom
  public element_connectivity_1d
  public evaluate_ordering
  public read_icem_msh
  public read_k_file
  public get_final_real

contains
!
!###################################################################################
!
  subroutine add_matching_mesh()
  !*Description:* adds a matching venous mesh to an arterial mesh
    use arrays,only: dp,elems,elem_cnct,elem_direction,elem_field,&
         elem_nodes,elem_ordrs,elem_symmetry,elems_at_node,&
         nodes,node_xyz,num_elems,&
         num_nodes,num_units,units,num_arterial_elems
    use indices
    use other_consts,only: PI
    use diagnostics, only: enter_exit,get_diagnostics_level
    implicit none
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_ADD_MATCHING_MESH" :: ADD_MATCHING_MESH  
    !Parameters to become inputs
    real(dp) :: offset(3)
    logical :: REVERSE=.TRUE.
    !local variables
    integer :: num_nodes_new,num_elems_new,ne,ne_global,np,np_global,np0,nonode,np_m
    integer :: nj,ne_m,noelem,ne0,n,nindex,ne1,noelem0,nu,cap_conns,cap_term,np1,np2,counter,i,j
    integer, allocatable :: np_map(:)
    character(len=60) :: sub_name
    integer :: diagnostics_level

    sub_name = 'add_matching_mesh'
    call enter_exit(sub_name,1)
    call get_diagnostics_level(diagnostics_level)
    !Ultimately offset should be an input argument
    offset(1)=0.0_dp
    offset(2)=1e-6_dp
    offset(3)=0.0_dp


    allocate(np_map(num_nodes))
!!! increase the size of node and element arrays to accommodate the additional elements
    ! the number of nodes after adding mesh will be:
    num_nodes_new = 2*num_nodes
    ! the number of elems after adding mesh will be:
    num_elems_new = 2*num_elems + num_units
  
    call reallocate_node_elem_arrays(num_elems_new,num_nodes_new)
    noelem0=0
    ne0 = num_elems ! the starting local element number
    ne_global = elems(ne0) ! assumes this is the highest element number (!!!)
    np0 = num_nodes ! the starting local node number
    np_global = nodes(np0) ! assumes this is the highest node number (!!!)

    do nonode=1,num_nodes
       np=np_global+nonode
       np_m=nodes(nonode)
       np_map(np_m)=np !maps new to old node numbering
       nodes(np0+nonode)=np
       do nj=1,3
         node_xyz(nj,np)=node_xyz(nj,np_m)+offset(nj)
       enddo
       elems_at_node(np,0)=0 !initialise
     !Doesnt map versions, would be added here
    enddo
    
    do noelem=1,num_elems
        ne=ne_global+noelem
        elem_field(ne_group,ne)=2.0_dp!VEIN
        ne_m=elems(noelem)
        elem_field(ne_group,ne_m)=0.0_dp!ARTERY
        elems(ne0+noelem)=ne
        if(.NOT.REVERSE)then
          elem_nodes(1,ne)=np_map(elem_nodes(1,ne_m))
          elem_nodes(2,ne)=np_map(elem_nodes(2,ne_m))
          elem_cnct(1,0,ne)=elem_cnct(1,0,ne_m)
          elem_cnct(-1,0,ne)=elem_cnct(-1,0,ne_m)
          do n=1,elem_cnct(1,0,ne)
            elem_cnct(1,n,ne)=elem_cnct(1,n,ne_m)+ne0
          enddo
          do n=1,elem_cnct(-1,0,ne)
            elem_cnct(-1,n,ne)=elem_cnct(-1,n,ne_m)+ne0
          enddo
        else
          elem_nodes(1,ne)=np_map(elem_nodes(2,ne_m))
          elem_nodes(2,ne)=np_map(elem_nodes(1,ne_m))
          elem_cnct(-1,0,ne)=elem_cnct(1,0,ne_m)
          elem_cnct(1,0,ne)=elem_cnct(-1,0,ne_m)
          do n=1,elem_cnct(-1,0,ne)        
            elem_cnct(-1,n,ne)=elem_cnct(1,n,ne_m)+ne0
          enddo
          do n=1,elem_cnct(1,0,ne)
            elem_cnct(1,n,ne)=elem_cnct(-1,n,ne_m)+ne0
          enddo
        endif
        !if worrying about regions and versions do it here
        elems_at_node(elem_nodes(1,ne),0)=elems_at_node(elem_nodes(1,ne),0)+1
        elems_at_node(elem_nodes(1,ne),elems_at_node(elem_nodes(1,ne),0))=ne
        elems_at_node(elem_nodes(2,ne),0)=elems_at_node(elem_nodes(2,ne),0)+1
        elems_at_node(elem_nodes(2,ne),elems_at_node(elem_nodes(2,ne),0))=ne
        nindex=no_gen
        elem_ordrs(nindex,ne)=elem_ordrs(nindex,ne_m)
        nindex=no_sord
        elem_ordrs(nindex,ne)=elem_ordrs(nindex,ne_m)
        nindex=no_hord
        elem_ordrs(nindex,ne)=elem_ordrs(nindex,ne_m)
      enddo

     !update current no of nodes and elements to determine connectivity
     np0=np !current highest node
     ne1=ne !current highest element
     noelem0=num_elems+noelem0

     cap_conns=0
     cap_term=0
     do nu=1,num_units
         ne=units(nu)
         cap_term=cap_term+1
         np1=elem_nodes(2,ne)
         np2=np_map(np1)
         noelem0=noelem0+1
         ne1=ne1+1
         elems(noelem0)=ne1
         elem_nodes(1,ne1)=np1
         elem_nodes(2,ne1)=np2
         elems_at_node(np1,0)=elems_at_node(np1,0)+1
         elems_at_node(np1,elems_at_node(np1,0))=ne1
         elems_at_node(np2,0)=elems_at_node(np2,0)+1
         elems_at_node(np2,elems_at_node(np2,0))=ne1
         elem_cnct(1,elem_cnct(1,0,ne)+1,ne)=ne1
         elem_cnct(1,0,ne)=elem_cnct(1,0,ne)+1
         elem_cnct(-1,elem_cnct(-1,0,ne+ne_global)+1,ne+ne_global)=ne1
         elem_cnct(-1,0,ne+ne_global)=elem_cnct(-1,0,ne+ne_global)+1
         elem_cnct(-1,0,ne1)=1
         elem_cnct(1,0,ne1)=1
         elem_cnct(-1,1,ne1)=ne
         elem_cnct(1,1,ne1)=ne+ne0
         nindex=no_gen
         elem_ordrs(nindex,ne1)=elem_ordrs(nindex,ne_m)
         nindex=no_sord
         elem_ordrs(nindex,ne1)=elem_ordrs(nindex,ne_m)
         nindex=no_hord
         elem_ordrs(nindex,ne1)=elem_ordrs(nindex,ne_m)
         elem_field(ne_group,ne1)=1.0_dp!connection between meshes
    enddo
    if(diagnostics_level.GT.1)then
      print *, 'Number of connections', cap_term
    endif
 
    num_nodes=num_nodes_new
    num_arterial_elems = num_elems
    num_elems=num_elems_new
  
    !calculate the element lengths and directions for the venous elements
    do ne=num_arterial_elems+1,num_elems
       np1=elem_nodes(1,ne)
       np2=elem_nodes(2,ne)
       elem_field(ne_length,ne) = DSQRT((node_xyz(1,np2) - &
            node_xyz(1,np1))**2 + (node_xyz(2,np2) - &
            node_xyz(2,np1))**2 + (node_xyz(3,np2) - &
            node_xyz(3,np1))**2)   
       do j=1,3
          elem_direction(j,ne) = (node_xyz(j,np2) - &
               node_xyz(j,np1))/elem_field(ne_length,ne)           
       enddo
    enddo

  
    if(diagnostics_level.GT.1)then 
        print *, "num_nodes=",num_nodes
        print *, "num_arterial_elems=",num_arterial_elems
        print *, "num_elems=",num_elems
   	 	!print out new node geometry, number of connected elements for each element and element connectivity array
   		print *,"element nodes:"
   		DO ne=1,num_elems
   	   		DO nj=1,2 !each element has 2 nodes in 1D tree
   	     		print *, "elem_nodes(nj,ne)",nj,ne,"=",elem_nodes(nj,ne)
   	   		ENDDO
   		ENDDO
   
    		DO n=1,num_nodes
       		DO nj=1,3
         		print *,"node_xyz(nj,n)",nj,n,"=",node_xyz(nj,n)
       		ENDDO
    		ENDDO
   
    		DO n=1,num_nodes
       		print *," "
       		print *,"node",n
       		print *, "total number of elements connected",elems_at_node(n,0)
       		DO ne=1,elems_at_node(n,0)
          		print *, "element",elems_at_node(n,ne)
       		ENDDO
    		ENDDO
 
		! total count of upstream elements connected to element ne elem_cnct(-1,0,ne)
		! upstream elements elem_cnct(-1,counter,ne)
		! total count of downstream elements connected to element ne elem_cnct(1,0,ne)
		! downstream elements elem_cnct(1,counter,ne)
   		DO ne=1,num_elems
   	    		print *,""
   	    		print *,"element",ne 
       		IF(elem_cnct(-1,0,ne).gt.0)THEN
       	    		print *, "total number of upstream elements:",elem_cnct(-1,0,ne)
       			DO counter=1,elem_cnct(-1,0,ne)
          			print *, "upstream element",elem_cnct(-1,counter,ne)
       	    		ENDDO
       		ENDIF
       		IF(elem_cnct(1,0,ne).gt.0)THEN
       	    		print *, "total number of downstream elements:",elem_cnct(1,0,ne)
       			DO counter=1,elem_cnct(1,0,ne)
          			print *, "downstream element",elem_cnct(1,counter,ne)
       	    		ENDDO
       		ENDIF    
    		ENDDO   
    		
    		DO ne=1,num_elems
    			do counter=1,3		
    		  		print *, "elem_ordrs(",counter,",",ne,")=",elem_ordrs(counter,ne)
    			ENDDO
    		enddo
    		print *, elem_ordrs
    		 
    endif !diagnostics_level
       
    deallocate(np_map)
    call enter_exit(sub_name,2)

  end subroutine add_matching_mesh
!
!###################################################################################
!
  subroutine append_units()
  !*Description:* Appends terminal units at the end of a tree structure
    use arrays,only: dp, elem_cnct,elem_symmetry,elem_units_below,&
         num_elems,num_units,units,unit_field
    use indices,only: num_nu
    use diagnostics, only: enter_exit,get_diagnostics_level
    implicit none
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_APPEND_UNITS" :: APPEND_UNITS

    integer :: ne,ne0,nu
    character(len=60) :: sub_name
    integer:: diagnostics_level

    sub_name = 'append_units'
    call enter_exit(sub_name,1)
    call get_diagnostics_level(diagnostics_level)

    num_units = 0
    DO ne=1,num_elems
       IF(elem_cnct(1,0,ne).eq.0)THEN! terminal element
          num_units=num_units+1
       ENDIF
    ENDDO
    
    if(diagnostics_level.GT.1)then
		print *,"num_units=",num_units
	endif

    if(allocated(units))then !increasing the array size; just overwrite
       deallocate(units)
       deallocate(unit_field)
    endif
    allocate(units(num_units))
    allocate(unit_field(num_nu,num_units))

    unit_field=0.0_dp
    units=0
    elem_units_below(1:num_elems) = 0 !initialise the number of terminal units below a branch

    nu=0
    DO ne=1,num_elems
       IF(elem_cnct(1,0,ne).eq.0)THEN
          nu=nu+1
          units(nu)=ne     !Set up units array containing terminals
          elem_units_below(ne)=1
       ENDIF
    ENDDO

    ! count the effective number of elements below each branch
    do ne=num_elems,2,-1
       ne0=elem_cnct(-1,1,ne)
       elem_units_below(ne0) = elem_units_below(ne0) &
            + elem_units_below(ne)*elem_symmetry(ne)
    enddo !ne
	
	if(diagnostics_level.GT.1)then
		do ne=1,num_elems
			print *,"elem_units_below(",ne,")",elem_units_below(ne)
		enddo
	endif

    call enter_exit(sub_name,2)

  end subroutine append_units

!
!###################################################################################
!
  subroutine calc_capillary_unit_length(num_convolutes,num_generations)
  !*Description:* Calculates the effective length of a capillary unit based on its total resistance
  ! and assumed radius, given the number of terminal convolute connections and the number of
  ! generations of symmetric intermediate villous trees 
    use arrays,only: dp,num_units,units,elem_field,elem_direction, &
                     node_xyz,elem_nodes,elem_cnct
    use diagnostics, only: enter_exit,get_diagnostics_level
    use other_consts, only: PI
    use indices, only: ne_length,ne_radius
    implicit none
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_CALC_CAPILLARY_UNIT_LENGTH" :: CALC_CAPILLARY_UNIT_LENGTH

    integer, intent(inout) :: num_convolutes,num_generations

    real(dp) :: int_length,int_radius,cap_length,cap_radius,seg_length,viscosity, &
                seg_resistance,cap_resistance,terminal_resistance,total_resistance,cap_unit_radius,total_length
    real(dp) :: int_radius_in, int_radius_out
    real(dp),allocatable :: resistance(:)
    integer :: ne,nu,i,j,np1,np2,nc,nv
    integer :: AllocateStatus
    character(len=60) :: sub_name
    integer:: diagnostics_level

    sub_name = 'calc_terminal_unit_length'
    call enter_exit(sub_name,1)
    call get_diagnostics_level(diagnostics_level)

    !check number of capillary convolutes and number of intermediate villous tree generations
    if (num_convolutes.LE.0)then
      num_convolutes = 6
    endif
    if (num_generations.LE.0)then
      num_generations = 3
    endif
    if (diagnostics_level.GE.1)then
      print *, "num_convolutes=",num_convolutes
      print *, "num_generations=",num_generations
    endif

    allocate (resistance(num_convolutes+1), STAT = AllocateStatus)
    if (AllocateStatus /= 0)then
       STOP "*** Not enough memory for resistance array ***"
    endif

    ne =units(1) !Get a terminal unit
    nc = elem_cnct(1,1,ne) !capillary unit is downstream of a terminal unit
    nv =  elem_cnct(1,1,nc) !vein is downstream of the capillary
    int_radius_in = (elem_field(ne_radius,ne)+elem_field(ne_radius,nv))/2.0_dp ! mm radius of inlet_faces intermediate villous (average of artery and vein)
    int_radius_out=(0.03_dp + 0.03_dp/2.0_dp)/2.0_dp ! mm radius of mature intermediate villous (average of artery and vein)
    int_length=1.5_dp !mm Length of each intermediate villous
    cap_length=3_dp/num_convolutes !mm length of capillary convolutes
    cap_radius=0.0144_dp/2.0_dp !radius of capillary convolutes
    seg_length=int_length/num_convolutes !lengh of each intermediate villous segment
    viscosity=0.33600e-02_dp !Pa.s !viscosity: fluid viscosity
    cap_unit_radius = 0.03_dp
    cap_resistance=(8.d0*viscosity*cap_length)/(PI*cap_radius**4) !resistance of each capillary convolute segment
	total_resistance = 0
	

    do j=1,num_generations
      int_radius = int_radius_in - (int_radius_in-int_radius_out)/num_generations*j

      seg_resistance=(8.d0*viscosity*seg_length)/(PI*int_radius**4) !resistance of each intermediate villous segment
      !calculate total resistance of terminal capillary conduits
      i=1
      resistance(i)= cap_resistance + 2.d0*seg_resistance
      do i=2,num_convolutes
        resistance(i)=2.d0*seg_resistance + 1/(1/cap_resistance + 1/resistance(i-1))
      enddo
      terminal_resistance = resistance(num_convolutes) !Pa . s per mm^3 total resistance of terminal capillary conduits
	
      !We have symmetric generations of intermediate villous trees so we can calculate the total resistance
      !of the system by summing the resistance of each generation

      total_resistance = total_resistance + terminal_resistance/2**j
    enddo

    total_length = total_resistance*(PI*cap_unit_radius**4)/(8.d0*viscosity)

    if(diagnostics_level.GE.1)then
      print *, "terminal_resistance=",terminal_resistance
      print *, "total_resistance=",total_resistance
    endif

    !set the effective length of each capillary unit based the total resistance of capillary convolutes   
    do nu=1,num_units
      ne =units(nu) !Get a terminal unit
      nc = elem_cnct(1,1,ne) !capillary unit is downstream of a terminal unit
      !update element radius
      elem_field(ne_radius,nc) = cap_unit_radius
      !update element length   
      elem_field(ne_length,nc) = total_length
      !update element direction
      np1=elem_nodes(1,nc)
      np2=elem_nodes(2,nc)
      do j=1,3
        elem_direction(j,nc) = (node_xyz(j,np2) - &
               node_xyz(j,np1))/elem_field(ne_length,nc)           
      enddo
    enddo

    call enter_exit(sub_name,2)

  end subroutine calc_capillary_unit_length
!
!###################################################################################
!
  subroutine define_1d_elements(ELEMFILE)
  !*Description:* Reads in an element ipelem file to define a geometry
    use arrays,only: dp, elem_direction,elem_field,elems,elem_cnct,elem_nodes,&
         elem_ordrs,elem_symmetry,elems_at_node,elem_units_below,&
         node_xyz,num_elems,num_nodes
    use indices
    use diagnostics, only: enter_exit,get_diagnostics_level
    implicit none
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_DEFINE_1D_ELEMENTS" :: DEFINE_1D_ELEMENTS

    character(len=MAX_FILENAME_LEN), intent(in) :: ELEMFILE
    !     Local Variables
    integer :: ibeg,iend,ierror,i_ss_end,j,ne,ne_global,&
         nn,np,np1,np2,np_global
    character(LEN=132) :: ctemp1
    character(LEN=40) :: sub_string
    character(len=60) :: sub_name
    integer :: diagnostics_level

    sub_name = 'define_1d_elements'
    call enter_exit(sub_name,1)
    call get_diagnostics_level(diagnostics_level)

    open(10, file=ELEMFILE, status='old')

    read_number_of_elements : do
       read(unit=10, fmt="(a)", iostat=ierror) ctemp1
       if(index(ctemp1, "elements")> 0) then
          call get_final_integer(ctemp1,num_elems)
          if(diagnostics_level.GT.1)then
          	print *, "num_elems", num_elems
          endif
          exit read_number_of_elements
       endif
    enddo read_number_of_elements

!!! allocate memory for element arrays
    if(allocated(elems)) deallocate(elems)
    allocate(elems(num_elems))
    if(allocated(elem_cnct)) deallocate(elem_cnct)
    allocate(elem_cnct(-1:1,0:2,0:num_elems))
    if(allocated(elem_nodes)) deallocate(elem_nodes)
    allocate(elem_nodes(2,num_elems))
    if(allocated(elem_ordrs)) deallocate(elem_ordrs)
    allocate(elem_ordrs(num_ord,num_elems))
    if(allocated(elem_symmetry)) deallocate(elem_symmetry)
    allocate(elem_symmetry(num_elems))
    if(allocated(elem_units_below)) deallocate(elem_units_below)
    allocate(elem_units_below(num_elems))
    if(allocated(elems_at_node)) deallocate(elems_at_node)
    allocate(elems_at_node(num_nodes,0:3))
    if(allocated(elem_field)) deallocate(elem_field)
    allocate(elem_field(num_ne,num_elems))
    if(allocated(elem_direction)) deallocate(elem_direction)
    allocate(elem_direction(3,num_elems))

!!! initialise element arrays
    elems=0
    elem_nodes=0
    elem_symmetry = 1
    elem_field = 0.0_dp

    ne=0
    !each element has 2 nodes in 1D tree - this is defined in array elem_nodes: elem_nodes (1,element_number ne) = node np"	
    read_an_element : do
       !.......read element number
       read(unit=10, fmt="(a)", iostat=ierror) ctemp1
       if(index(ctemp1, "Element")> 0) then
          call get_final_integer(ctemp1,ne_global) !get element number
          ne=ne+1
          elems(ne)=ne_global
             read_element_nodes : do
             read(unit=10, fmt="(a)", iostat=ierror) ctemp1
             if(index(ctemp1, "global")> 0) then !found the correct line
                iend=len(ctemp1)
                ibeg=index(ctemp1,":")+1 !get location of first integer in string
                sub_string = adjustl(ctemp1(ibeg:iend)) ! get the characters beyond : remove leading blanks
                i_ss_end=len(sub_string) !get the end location of the sub-string
                ibeg=1
                do nn=1,2
                   iend=index(sub_string," ") !get location of first blank in sub-string
                   read (sub_string(ibeg:iend-1), '(i7)' ) np_global
                   call get_local_node(np_global,np) ! get local node np for global node
                   elem_nodes(nn,ne)=np ! the local node number, not global             
                   if(diagnostics_level.GT.1)then
                   		print *,"elem_nodes(nn,ne)", nn, ne, "= np", np
                   endif
                   sub_string = adjustl(sub_string(iend:i_ss_end)) ! get chars beyond blank, remove leading blanks
                enddo
                exit read_element_nodes
             endif !index
          enddo read_element_nodes
          if(ne.ge.num_elems) exit read_an_element
       endif

    enddo read_an_element

    close(10)

    ! calculate the element lengths and directions
    do ne=1,num_elems
       np1=elem_nodes(1,ne)
       np2=elem_nodes(2,ne)
       elem_field(ne_length,ne) = DSQRT((node_xyz(1,np2) - &
            node_xyz(1,np1))**2 + (node_xyz(2,np2) - &
            node_xyz(2,np1))**2 + (node_xyz(3,np2) - &
            node_xyz(3,np1))**2)   
       do j=1,3
          elem_direction(j,ne) = (node_xyz(j,np2) - &
               node_xyz(j,np1))/elem_field(ne_length,ne)           
       enddo !j
    enddo

    call element_connectivity_1d
    call evaluate_ordering

    call enter_exit(sub_name,2)

  END subroutine define_1d_elements
!
!###################################################################################
!
  subroutine define_node_geometry(NODEFILE)
  !*Description:* Reads in an ipnode file to define a tree geometry
    use arrays,only: dp,nodes,node_field,node_xyz,num_nodes
    use diagnostics, only: enter_exit,get_diagnostics_level
    use indices
    use other_consts, only: MAX_FILENAME_LEN
    implicit none
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_DEFINE_NODE_GEOMETRY" :: DEFINE_NODE_GEOMETRY

    character(len=MAX_FILENAME_LEN), intent(in) :: NODEFILE !Input nodefile
    !     Local Variables
    integer :: i,ierror,np,np_global,&
         num_versions,nv,NJT
    character(LEN=132) :: ctemp1
    LOGICAL :: versions
    real(dp) :: point
    character(len=60) :: sub_name
    integer :: diagnostics_level
  
    sub_name = 'define_node_geometry'
    call enter_exit(sub_name,1)
    call get_diagnostics_level(diagnostics_level)
   

    versions = .TRUE.
    NJT = 0
    open(10, file=NODEFILE, status='old')

    !.....read in the total number of nodes. read each line until one is found
    !.....that has the correct keyword (nodes). then return the integer that is
    !.....at the end of the line
    read_number_of_nodes : do !define a do loop name
       read(unit=10, fmt="(a)", iostat=ierror) ctemp1 !read a line into ctemp1
       if(index(ctemp1, "nodes")> 0) then !keyword "nodes" is found in ctemp1
          call get_final_integer(ctemp1,num_nodes) !return the final integer
          if(diagnostics_level.GT.1)then
       	  	print *, "num_nodes",num_nodes
       	  endif
          exit read_number_of_nodes !exit the named do loop
       endif

    enddo read_number_of_nodes

    if(allocated(nodes)) deallocate (nodes)
    allocate (nodes(num_nodes))
    if(allocated(node_xyz)) deallocate (node_xyz)
    allocate (node_xyz(3,num_nodes))
    if(allocated(node_field)) deallocate (node_field)
    allocate (node_field(num_nj,num_nodes))
    nodes = 0 !initialise node index values
    node_xyz = 0.0_dp !initialise
    node_field = 0.0_dp !initialise

    !.....read in the number of coordinates
    read_number_of_coords : do !define a do loop name
       read(unit=10, fmt="(a)", iostat=ierror) ctemp1 !read a line into ctemp1
       if(index(ctemp1, "coordinates")> 0) then !keyword "coordinates" is found
          call get_final_integer(ctemp1,NJT) !return the final integer
          exit read_number_of_coords !exit the named do loop
       endif
    enddo read_number_of_coords

    !.....check whether versions are prompted (>1)
    read_versions : do !define a do loop name
       read(unit=10, fmt="(a)", iostat=ierror) ctemp1 !read a line into ctemp1
       if(index(ctemp1, "different")> 0) then !keyword "different" is found
          if(index(ctemp1, " N")> 0) then !keyword " N" is found
             versions=.false.
          endif
          exit read_versions !exit the named do loop
       endif
    enddo read_versions

!!! WARNING :: following should be in general code
    ! note that only the first version of coordinate is currently read in

    !.....read the coordinate, derivative, and version information for each node.
    np=0
    read_a_node : do !define a do loop name
       !.......read node number
       read(unit=10, fmt="(a)", iostat=ierror) ctemp1
       if(index(ctemp1, "Node")> 0) then
          call get_final_integer(ctemp1,np_global) !get node number
          np=np+1
          nodes(np)=np_global
          !.......read coordinates and derivatives
          do i=1,NJT ! for the NJT coordinates
             !...........coordinate
             num_versions=1
             if(versions)then
                read(unit=10, fmt="(a)", iostat=ierror) ctemp1
                call get_final_integer(ctemp1,num_versions)
             endif
             if(num_versions > 1)then
                read(unit=10, fmt="(a)", iostat=ierror) ctemp1 !temporary line
                read(unit=10, fmt="(a)", iostat=ierror) ctemp1
                call get_final_real(ctemp1,point)
                do nv=2,num_versions
                   read(unit=10, fmt="(a)", iostat=ierror) ctemp1 !temporary line
                   read(unit=10, fmt="(a)", iostat=ierror) ctemp1
                enddo
             else
                read(unit=10, fmt="(a)", iostat=ierror) ctemp1
                call get_final_real(ctemp1,point)
             endif
             node_xyz(i,np)=point
             if(diagnostics_level.GT.1)then
             	print *, "node_xyz(i,np)",i," ",np,"=", point
             endif
          enddo !i

       endif !index
                  
       if(np.ge.num_nodes) exit read_a_node
    enddo read_a_node

    close(10)
    
    call enter_exit(sub_name,2)

  END subroutine define_node_geometry
!
!###################################################################################
!
  subroutine define_rad_from_geom(ORDER_SYSTEM, CONTROL_PARAM, START_FROM, START_RAD, group_type_in, group_option_in)
  !*Description:* Defines vessel radius based on their geometric structure
    use arrays,only: dp,num_elems,elem_field,elem_ordrs,maxgen,elem_cnct,num_arterial_elems
    use indices
    use diagnostics, only: enter_exit,get_diagnostics_level
    implicit none
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_DEFINE_RAD_FROM_GEOM" :: DEFINE_RAD_FROM_GEOM
    
   character(LEN=100), intent(in) :: ORDER_SYSTEM,START_FROM
   character(LEN=100), optional :: group_type_in, group_option_in
   real(dp), intent(in) :: CONTROL_PARAM,START_RAD
   !Input options ORDER_SYSTEM=STRAHLER (CONTROL_PARAM=RDS), HORSFIELD (CONTROL_PARAM=RDH)

   !Local variables
   character(LEN=100) :: group_type, group_options
   integer :: ne_min,ne_max,nindex,ne,n_max_ord,n,ne_start,&
      inlet_faces_count
   real(dp) :: radius
   character(len=60) :: sub_name
   integer :: diagnostics_level

   sub_name = 'define_rad_from_geom'
    call enter_exit(sub_name,1)
    call get_diagnostics_level(diagnostics_level)
    !define list of elements you are going to operate on
    if(present(group_type_in))then
      group_type = group_type_in
    else!default to all
      group_type='all'
    endif
    if(group_type.eq.'all')then
       ne_min=1
       ne_max=num_elems
    elseif(group_type.eq.'efield')then
    elseif(group_type.eq.'arterial')then
       ne_min=1
       ne_max=num_arterial_elems
    elseif(group_type.eq.'venous')then
       ne_min=num_arterial_elems + 1
       ne_start=ne_min
       ne_max=num_elems          
    elseif(group_type.eq.'list')then
      read (START_FROM,'(I10)') ne_min
      read (group_option_in,'(I10)') ne_max
    endif
    if(diagnostics_level.GT.1)then
    	  print *,"ne_min=", ne_min, "ne_max=",ne_max
    	endif
    !Define start element
    if(group_type.ne.'venous')then
    		if(START_FROM.eq.'inlet_faces')then
      	inlet_faces_count=0
      		do ne=ne_min,ne_max
         		if(elem_cnct(-1,0,ne).eq.0)then
           			inlet_faces_count=inlet_faces_count+1
           			ne_start=ne
         		endif
         		if(inlet_faces_count.gt.1)then
            			WRITE(*,*) ' More than one inlet_faces in this group, using last found, ne = ',ne
         		endif
      		enddo
    		else!element number defined
       		read (START_FROM,'(I10)') ne_start
    		endif
	endif
	if(diagnostics_level.GT.1)then
		print *, "ne_start=",ne_start
	endif

    !Strahler and Horsfield ordering system
    if(ORDER_SYSTEM(1:5).EQ.'strah')THEN
      nindex=no_sord !for Strahler ordering
    else if(ORDER_SYSTEM(1:5).eq.'horsf')then
      nindex = no_hord !for Horsfield ordering
    endif

	if(diagnostics_level.GT.1)then
		print *, "ordering system= ",ORDER_SYSTEM(1:5)
	endif

    ne=ne_start
    n_max_ord=elem_ordrs(nindex,ne)
    elem_field(ne_radius,ne)=START_RAD

	if(diagnostics_level.GT.1)then
		print *, "start radius START_RAD = ",START_RAD
	endif

    do ne=ne_min,ne_max
     radius=10.0_dp**(log10(CONTROL_PARAM)*dble(elem_ordrs(nindex,ne)-n_max_ord)&
        +log10(START_RAD))
     elem_field(ne_radius,ne)=radius  
     elem_field(ne_radius_in,ne)=radius
     elem_field(ne_radius_out,ne)=radius
     if(diagnostics_level.GT.1)then
     	print *,"radius for element",ne,"=",elem_field(ne_radius,ne)
     	print *,"radius in for element",ne,"=",elem_field(ne_radius_in,ne)
     	print *,"radius out for element",ne,"=",elem_field(ne_radius_out,ne)
     endif
    enddo

    call enter_exit(sub_name,2)

  END subroutine define_rad_from_geom
!
!###########################################################################
!
  subroutine element_connectivity_1d()
  !*Description:* Calculates element connectivity in 1D and stores in elelem_cnct
    use arrays,only: elem_cnct,elem_nodes,elems_at_node,num_elems,num_nodes
    use diagnostics, only: enter_exit,get_diagnostics_level
    implicit none
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_ELEMENT_CONNECTIVITY_1D" :: ELEMENT_CONNECTIVITY_1D
  
    !     Local Variables
    integer :: ne,ne2,nn,noelem,np,np2,np1,counter,orphan_counter
    integer,parameter :: NNT=2
    character(len=60) :: sub_name
    integer :: orphan_nodes(num_nodes)
    integer :: diagnostics_level

    sub_name = 'element_connectivity_1d'
    call enter_exit(sub_name,1)
    call get_diagnostics_level(diagnostics_level)

    elem_cnct = 0 !initialise

    ! calculate elems_at_node array: stores the elements that nodes are in
    ! elems_at_node(node np,0)= total number of elements connected to this node
    ! elems_at_node(node np, index of each connected element starting at 1) = connected element
    elems_at_node = 0 !initialise

    DO ne=1,num_elems
       DO nn=1,2
          np=elem_nodes(nn,ne)
          elems_at_node(np,0)=elems_at_node(np,0)+1
          elems_at_node(np,elems_at_node(np,0))=ne ! local element that np is in
        ENDDO !nn
    ENDDO !noelem

    if(diagnostics_level.GT.1)then
    		DO nn=1,num_nodes
       		print *," "
       		print *,"node",nn
       		print *,"total number of elements connected",elems_at_node(nn,0)
       		DO ne=1,elems_at_node(nn,0)
          		print *,"element",elems_at_node(nn,ne)
       		ENDDO
    		ENDDO
    endif
    
    !check for nodes with 0 elements - exit if any are found
    orphan_counter = 0
    DO nn=1,num_nodes
		if(elems_at_node(nn,0).EQ.0)then
			orphan_counter = orphan_counter + 1
			orphan_nodes(orphan_counter) = nn
		endif
	ENDDO
	if(orphan_counter.GT.0)then
		print *, "found",orphan_counter,"node(s) not connected to any elements"
		do counter=1,orphan_counter
			print *,"node",orphan_nodes(counter),"is not connected to any elements"
		enddo
		call exit(0)
	endif

    ! calculate elem_cnct array: stores the connectivity of all elements
    
    elem_cnct=0 !initialise all elem_cnct

    DO ne=1,num_elems
       !     ne_global=elems(noelem)
       IF(NNT == 2) THEN !1d
          np1=elem_nodes(1,ne) !first local node
          np2=elem_nodes(2,ne) !second local node

          DO noelem=1,elems_at_node(np2,0) !for each element connected to node np2
             ne2=elems_at_node(np2,noelem) !get the element number connected to node np2
             IF(ne2 /= ne)THEN !if element connected to node np2 is not the current element ne
                elem_cnct(-1,0,ne2)=elem_cnct(-1,0,ne2)+1
                elem_cnct(-1,elem_cnct(-1,0,ne2),ne2)=ne !previous element              
                elem_cnct(1,0,ne)=elem_cnct(1,0,ne)+1
                elem_cnct(1,elem_cnct(1,0,ne),ne)=ne2
             ENDIF !ne2
          ENDDO !noelem2


       ENDIF
    ENDDO

	! total count of upstream elements connected to element ne elem_cnct(-1,0,ne)
	! upstream elements elem_cnct(-1,counter,ne)
	! total count of downstream elements connected to element ne elem_cnct(1,0,ne)
	! downstream elements elem_cnct(1,counter,ne)
	if(diagnostics_level.GT.1)then
   		DO ne=1,num_elems
   	    		print *,""
   	    		print *,"element",ne 
       		IF(elem_cnct(-1,0,ne).gt.0)THEN
       	    		print *,"total number of upstream elements:",elem_cnct(-1,0,ne)
       			DO counter=1,elem_cnct(-1,0,ne)
          			print *,"upstream element",elem_cnct(-1,counter,ne)
       	    		ENDDO
       		ENDIF
       		IF(elem_cnct(1,0,ne).gt.0)THEN
       	    		print *,"total number of downstream elements:",elem_cnct(1,0,ne)
       			DO counter=1,elem_cnct(1,0,ne)
          			print *,"downstream element",elem_cnct(1,counter,ne)
       	    		ENDDO
       		ENDIF
    		ENDDO
    endif

    call enter_exit(sub_name,2)

  END subroutine element_connectivity_1d

!
!###################################################################################
!
  subroutine evaluate_ordering()
  !*Description:* calculates generations, Horsfield orders, Strahler orders for a given tree
   
    use arrays,only: elem_cnct,elem_nodes,elem_ordrs,elem_symmetry,&
         elems_at_node,num_elems,num_nodes,maxgen
    use diagnostics, only: enter_exit,get_diagnostics_level
    implicit none
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_EVALUATE_ORDERING" :: EVALUATE_ORDERING

    integer :: inlet_faces,ne,ne0,ne2,noelem2,np,np2, &
         num_attach,n_children,n_generation, &
         n_horsfield,outlet_faces,STRAHLER,STRAHLER_ADD,temp1
    LOGICAL :: DISCONNECT,DUPLICATE
    character(len=60) :: sub_name
    integer :: diagnostics_level
    
    sub_name = 'evaluate_ordering'    
    call enter_exit(sub_name,1)
    call get_diagnostics_level(diagnostics_level)

    !Calculate generations, Horsfield orders, Strahler orders
    !.....Calculate branch generations

	elem_ordrs = 0
    maxgen=1
    DO ne=1,num_elems
       ne0=elem_cnct(-1,1,ne) !parent
       IF(ne0.NE.0)THEN
          n_generation=elem_ordrs(1,ne0) !parent generation
          IF(elem_cnct(1,0,ne0).EQ.1)THEN !single daughter
             elem_ordrs(1,ne)=n_generation + (elem_symmetry(ne)-1)
          ELSE IF(elem_cnct(1,0,ne0).GE.2)THEN
             elem_ordrs(1,ne)=n_generation+1
          ENDIF
       ELSE
          elem_ordrs(1,ne)=1 !generation 1
       ENDIF
       maxgen=max(maxgen,elem_ordrs(1,ne))
    ENDDO !noelem
	if(diagnostics_level.GT.1)then
		print *,"branch generations - maxgen",maxgen
	endif

    !.....Calculate the branch orders
    DO ne=num_elems,1,-1
       n_horsfield=MAX(elem_ordrs(2,ne),1)
       n_children=elem_cnct(1,0,ne) !number of child branches
       IF(n_children.EQ.1)THEN
          IF(elem_ordrs(1,elem_cnct(1,1,ne)).EQ.0)  n_children=0
       ENDIF
       STRAHLER=0
       STRAHLER_ADD=1
       IF(n_children.GE.2)THEN !branch has two or more daughters
          STRAHLER=elem_ordrs(3,elem_cnct(1,1,ne)) !first daughter
          DO noelem2=1,n_children !for all daughters
             ne2=elem_cnct(1,noelem2,ne) !global element # of daughter
             temp1=elem_ordrs(2,ne2) !Horsfield order of daughter
             IF(temp1.GT.n_horsfield) n_horsfield=temp1
             IF(elem_ordrs(3,ne2).LT.STRAHLER)THEN
                STRAHLER_ADD=0
             ELSE IF(elem_ordrs(3,ne2).GT.STRAHLER)THEN
                STRAHLER_ADD=0
                STRAHLER=elem_ordrs(3,ne2) !highest daughter
             ENDIF
          ENDDO !noelem2 (ne2)
          n_horsfield=n_horsfield+1 !Horsfield ordering
       ELSE IF(n_children.EQ.1)THEN
          ne2=elem_cnct(1,1,ne) !local element # of daughter
          n_horsfield=elem_ordrs(2,ne2)+(elem_symmetry(ne)-1)
          STRAHLER_ADD=elem_ordrs(3,ne2)+(elem_symmetry(ne)-1)
       ENDIF !elem_cnct
       elem_ordrs(2,ne)=n_horsfield !store the Horsfield order
       elem_ordrs(3,ne)=STRAHLER+STRAHLER_ADD !Strahler order
    ENDDO !noelem

    !       Check for disconnected nodes and number of inlet_faces and outlet_faces
    DUPLICATE=.FALSE.
    DO ne=1,num_elems
       np=elem_nodes(1,ne)
       np2=elem_nodes(2,ne)
       IF(np.EQ.np2)THEN
          DUPLICATE=.TRUE.
       ENDIF
    ENDDO

    DISCONNECT=.FALSE.
    inlet_faces=0
    outlet_faces=0
    DO np=1,num_nodes
       num_attach=elems_at_node(np,0)
       IF(num_attach.EQ.0)THEN
          DISCONNECT=.TRUE.
       ELSEIF(num_attach.EQ.1)THEN
          ne=elems_at_node(np,1)
          IF(elem_cnct(1,0,ne).EQ.0) outlet_faces=outlet_faces+1
         IF(elem_cnct(-1,0,ne).EQ.0) inlet_faces=inlet_faces+1
       ELSEIF(num_attach.GT.3)THEN
          WRITE(*,*) ' Node ',np,' attached to',num_attach,' elements'
       ENDIF
    ENDDO
  
    call enter_exit(sub_name,2)

  end subroutine evaluate_ordering
!
!###################################################################################
!
  subroutine reallocate_node_elem_arrays(num_elems_new,num_nodes_new)
  !*Description:* Reallocates the size of arrays when modifying geometries
  
    use arrays,only: dp,elems,elem_cnct,elem_direction,elem_field,&
         elem_ordrs,elem_nodes,&
         elem_symmetry,elem_units_below,elems_at_node,&
         nodes,node_field,node_xyz,num_elems,num_nodes
    use indices
    use diagnostics, only: enter_exit
    implicit none

!!! Parameters
    integer,intent(in) :: num_elems_new,num_nodes_new

!!! Local variables
    integer,allocatable :: nodelem_temp(:),enodes_temp(:,:),enodes_temp2(:,:,:)
    real(dp),allocatable :: xyz_temp(:,:),rnodes_temp(:,:)
    character(len=60) :: sub_name

    sub_name = 'reallocate_node_elem_arrays'
    call enter_exit(sub_name,1)

    allocate(nodelem_temp(num_nodes))
    nodelem_temp = nodes ! copy to temporary array
    deallocate(nodes) !deallocate initially allocated memory
    allocate(nodes(num_nodes_new))
    nodes(1:num_nodes)=nodelem_temp(1:num_nodes)
    deallocate(nodelem_temp) !deallocate the temporary array

    allocate(xyz_temp(3,num_nodes))
    xyz_temp=node_xyz
    deallocate(node_xyz)
    allocate(node_xyz(3,num_nodes_new))
    node_xyz(1:3,1:num_nodes)=xyz_temp(1:3,1:num_nodes)

    allocate(nodelem_temp(num_elems))
    nodelem_temp = elems ! copy to temporary array
    deallocate(elems) !deallocate initially allocated memory
    allocate(elems(num_elems_new))
    elems(1:num_elems)=nodelem_temp(1:num_elems)
    deallocate(nodelem_temp) !deallocate the temporary array

    allocate(enodes_temp(2,num_elems))
    enodes_temp=elem_nodes
    deallocate(elem_nodes)
    allocate(elem_nodes(2,num_elems_new))
    elem_nodes(1:2,1:num_elems)=enodes_temp(1:2,1:num_elems)
    deallocate(enodes_temp)

    allocate(rnodes_temp(num_ne,num_elems))
    rnodes_temp=elem_field
    deallocate(elem_field)
    allocate(elem_field(num_ne,num_elems_new))
    elem_field(1:num_ne,1:num_elems)=rnodes_temp(1:num_ne,1:num_elems)
    deallocate(rnodes_temp)
    elem_field(1:num_ne,num_elems+1:num_elems_new) = 0.0_dp

    allocate(rnodes_temp(3,num_elems))
    rnodes_temp=elem_direction
    deallocate(elem_direction)
    allocate(elem_direction(3,num_elems_new))
    elem_direction(1:3,1:num_elems)=rnodes_temp(1:3,1:num_elems)
    deallocate(rnodes_temp)
    elem_direction(1:3,num_elems+1:num_elems_new) = 0.0_dp

    allocate(rnodes_temp(num_nj,num_nodes))
    rnodes_temp=node_field
    deallocate(node_field)
    allocate(node_field(num_nj,num_nodes_new))
    node_field(1:num_nj,1:num_nodes)=rnodes_temp(1:num_nj,1:num_nodes)
    deallocate(rnodes_temp)
    node_field(1:num_nj,num_nodes+1:num_nodes_new)=0.0_dp

    allocate(nodelem_temp(num_elems))
    nodelem_temp = elem_symmetry ! copy to temporary array
    deallocate(elem_symmetry) !deallocate initially allocated memory
    allocate(elem_symmetry(num_elems_new))
    elem_symmetry(1:num_elems)=nodelem_temp(1:num_elems)
    deallocate(nodelem_temp) !deallocate the temporary array
    elem_symmetry(num_elems+1:num_elems_new)=1

    allocate(enodes_temp2(-1:1,0:2,0:num_elems))
    enodes_temp2=elem_cnct
    deallocate(elem_cnct)
    allocate(elem_cnct(-1:1,0:2,0:num_elems_new))
    elem_cnct(-1:1,0:2,0:num_elems)=enodes_temp2(-1:1,0:2,0:num_elems)
    deallocate(enodes_temp2)
    elem_cnct(-1:1,0:2,num_elems+1:num_elems_new) = 0

    allocate(enodes_temp(num_ord,num_elems))
    enodes_temp=elem_ordrs
    deallocate(elem_ordrs)
    allocate(elem_ordrs(num_ord,num_elems_new))
    elem_ordrs(1:num_ord,1:num_elems)=enodes_temp(1:num_ord,1:num_elems)
    deallocate(enodes_temp)
    elem_ordrs(1:num_ord,num_elems+1:num_elems_new) = 0

    allocate(nodelem_temp(num_elems))
    nodelem_temp=elem_units_below
    deallocate(elem_units_below)
    allocate(elem_units_below(num_elems_new))
    elem_units_below(1:num_elems)=nodelem_temp(1:num_elems)
    deallocate(nodelem_temp)
    elem_units_below(num_elems+1:num_elems_new)=0

    allocate(enodes_temp(num_nodes,0:3))
    enodes_temp=elems_at_node
    deallocate(elems_at_node)
    allocate(elems_at_node(num_nodes_new,0:3))
    elems_at_node(1:num_nodes,0:3)=enodes_temp(1:num_nodes,0:3)
    deallocate(enodes_temp)
    elems_at_node(num_nodes+1:num_nodes_new,0:3)=0

    call enter_exit(sub_name,2)

  end subroutine reallocate_node_elem_arrays

!
!###################################################################################
!
  subroutine get_final_real(string,rtemp)
    use arrays,only: dp
    implicit none
    character, intent(in) :: string*(132)
    integer :: ibeg,iend
    real(dp), intent(out) :: rtemp
    real(dp) :: rsign
    character :: sub_string*(40)

    iend=len(string) !get the length of the string
    ibeg=index(string,":")+1 !get location of real in string
    sub_string = adjustl(string(ibeg:iend)) ! get the characters beyond :
    iend=len(sub_string) !get the length of the sub-string
    if(sub_string(1:1).eq.'-')then !check whether negative
       rsign=-1.0_dp
       ibeg=2
    else
       rsign=1.0_dp
       ibeg=1
    endif
    read (sub_string(ibeg:iend), * ) rtemp !get real value
    rtemp=rtemp*rsign !apply sign to number

  end subroutine get_final_real

!
!###################################################################################
!
  subroutine get_final_string(string,rtemp)
    use arrays,only: dp
    implicit none
    character, intent(in) :: string*(132)
    integer :: ibeg,iend
    real(dp), intent(out) :: rtemp
    real(dp) :: rsign
    character :: sub_string*(40)

    iend=len(string) !get the length of the string
    ibeg=index(string,":")+1 !get location of real in string
    sub_string = adjustl(string(ibeg:iend)) ! get the characters beyond :
    iend=len(sub_string) !get the length of the sub-string
    if(sub_string(1:1).eq.'-')then !check whether negative
       rsign=-1.0_dp
       ibeg=2
    else
       rsign=1.0_dp
       ibeg=1
    endif
    read (sub_string(ibeg:iend), '(D25.17)' ) rtemp !get real value
    rtemp=rtemp*rsign !apply sign to number

  end subroutine get_final_string

!
!###################################################################################
!
  subroutine get_local_node(np_global,np_local)
    use arrays,only: nodes,num_nodes
    implicit none

    integer,intent(in) :: np_global
    integer,intent(out) :: np_local

    integer :: np
    logical :: found

    np=1   
    found=.false.
    do while ((.not.found).AND.(np.le.num_nodes))
       if(nodes(np).eq.np_global)then
          found=.true.
       else
          np=np+1
       endif
    enddo

    if(.not.found)then
       np = 0          
       write(*,'('' Global node '',I6,'' not in node list'')') np_global
    endif

    np_local = np
    return

  end subroutine get_local_node

!
!###################################################################################
!
  subroutine get_final_integer(string,num)
    implicit none
    character,intent(in) :: string*(132)
    integer,intent(out) :: num
    integer :: ibeg,iend,nsign,ntemp
    character :: sub_string*(40)

    iend=len(string) !get the length of the string
    ibeg=index(string,":")+1 !get location of integer in string, follows ":"
    sub_string = adjustl(string(ibeg:iend)) ! get the characters beyond ":"
    iend=len(sub_string) !length of the sub-string
    if(sub_string(1:1).eq.'-')then !check for negative sign
       nsign=-1
       ibeg=2
    else
       nsign=1
       ibeg=1
    endif
    read (sub_string(ibeg:iend), '(i10)' ) ntemp !get integer values
    ntemp=ntemp*nsign !apply sign to number

    num=ntemp !return the integer value

  end subroutine get_final_integer

!
! ##########################################################################      
!

  function inlist(item,ilist)
!!! dummy arguments
    integer :: item,ilist(:)
! local variables
    integer :: n
    logical :: inlist

    inlist = .false.
    do n=1,size(ilist)
       if(item == ilist(n)) inlist = .true.
    enddo

  end function inlist
!
!###########################################################################################
!

!------------------------------------------------
subroutine read_icem_msh(filename)
    use arrays,only: dp,nodes,node_xyz,num_nodes,internal_faces,num_faces,&
      inlet_faces,num_inlet_faces,num_outlet_faces,num_wall_faces,&
      outlet_faces,wall_faces,all_faces,num_all_faces
    use diagnostics, only: enter_exit,get_diagnostics_level
    use indices
    use other_consts, only: MAX_FILENAME_LEN
    implicit none
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_READ_ICEM_MSH" :: READ_ICEM_MSH

    character(len=MAX_FILENAME_LEN), intent(in) :: filename !Input nodefile
    !     Local Variables
    character (100) :: cur_line, no_data
    character*4:: hexaval1,hexaval2,hexaval3,hexaval4,hexaval5,hexaval6
    integer, parameter :: line_buf_len= 1024*4
    character(LEN=line_buf_len) :: InS
    integer :: input, size1, status1, startunit
    integer :: narray(32), decimaln,i,j,k,m
    integer :: input_stat, ii, results, indx,Startstate, lineNum
    real(dp) :: x, y, z
    character(len=60) :: sub_name
    integer :: diagnostics_level

    sub_name = 'read_icem_msh'
    call enter_exit(sub_name,1)
    call get_diagnostics_level(diagnostics_level)

    lineNum=0
    input=1
    if(diagnostics_level.gt.1)write(*,*) 'Opening file : ', filename

    open ( unit = input, file = filename, status = 'old', iostat = input_stat )
    if ( input_stat .ne. 0 ) then
        write ( *, '(a)' ) ''
        write ( *, '(a)' ) 'GMSH FILE READ  - Fatal error!'
        write ( *, '(a)' ) '  Could not open input file:' &
         // trim (filename)
        stop 1
    end if

    read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
    lineNum=lineNum+1

    do while ( size1<6 .OR. cur_line(size1-1:size1).NE.')(' )
        read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
        lineNum=lineNum+1
    enddo

    !*******Start reading the vertices(x,y,z)*******************************
    read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
    num_nodes=0 !was I
    do while ( size1>6 .OR. cur_line(size1-1:size1).NE.'))' )
        num_nodes=num_nodes+1
        read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
    enddo

    rewind(input)
    if(diagnostics_level.gt.1)write(*,*) 'Allocating vertices (node_xyz) : ', num_nodes
    allocate(node_xyz(num_nodes,3)) !was vertices in original code

    do ii=1, lineNum
        read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) no_data
    enddo
    do ii=1,num_nodes
        read(input, *) node_xyz(ii,:)
        lineNum=lineNum+1
        if(diagnostics_level.gt.1)write(*,*) 'new_node : ',ii, node_xyz(ii,:)
    enddo
    !******* END reading the vertices(x,y,z)***********************************


    do while ( size1<6 .OR. cur_line(size1-1:size1).NE.')(' )
        read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
        lineNum=lineNum+1
    enddo
    !*******Read internal faces(nod1,node2,node3,node4,RightElement,Leftelement)********
    read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
    num_faces=0
    do while ( size1>6 .OR. cur_line(size1:size1).NE.')' )
        num_faces=num_faces+1
        read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
    enddo
    rewind(input)
    if(diagnostics_level.gt.1)write(*,*) 'Allocating internal_faces : ', num_faces
    allocate(internal_faces(num_faces,6))
    do ii=1, lineNum
        read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) no_data
    enddo
    do ii=1,num_faces
        read(input, *) hexaval1,hexaval2,hexaval3,hexaval4,hexaval5,hexaval6
        read(hexaval1, '(z4)') internal_faces(ii,1)
        read(hexaval2, '(z4)') internal_faces(ii,2)
        read(hexaval3, '(z4)') internal_faces(ii,3)
        read(hexaval4, '(z4)') internal_faces(ii,4)
        read(hexaval5, '(z4)') internal_faces(ii,5)
        read(hexaval6, '(z4)') internal_faces(ii,6)
        lineNum=lineNum+1
        if(diagnostics_level.gt.1)write(*,*) 'new face : ',ii, internal_faces(ii,:)
    enddo
    !******* End reading internal faces *************************************************
    read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
    lineNum=lineNum+1
    read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
    lineNum=lineNum+1
    read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
    lineNum=lineNum+1

    do while ( size1<6 .OR. cur_line(size1-1:size1).NE.')(' )
        read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
        lineNum=lineNum+1
    enddo

    !*******Read inlet_faces faces(nod1,node2,node3,node4,RightElement,Leftelement)********
    read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
    num_inlet_faces=0
    do while ( size1>6 .OR. cur_line(size1:size1).NE.')' )
        num_inlet_faces=num_inlet_faces+1
        read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
    enddo
    rewind(input)
    if(diagnostics_level.gt.1)write(*,*) 'Allocating inlet_faces : ', num_inlet_faces
    allocate(inlet_faces(num_inlet_faces,6))
    do ii=1, lineNum
        read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) no_data
    enddo
    do ii=1,num_inlet_faces
        read(input, *) hexaval1,hexaval2,hexaval3,hexaval4,hexaval5,hexaval6
        read(hexaval1, '(z4)') inlet_faces(ii,1)
        read(hexaval2, '(z4)') inlet_faces(ii,2)
        read(hexaval3, '(z4)') inlet_faces(ii,3)
        read(hexaval4, '(z4)') inlet_faces(ii,4)
        read(hexaval5, '(z4)') inlet_faces(ii,5)
        read(hexaval6, '(z4)') inlet_faces(ii,6)
        lineNum=lineNum+1
        if(diagnostics_level.gt.1)write(*,*) 'inlet face : ',ii, inlet_faces(ii,:)

    enddo


    !******* End reading inlet_faces faces *************************************************
    read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
    lineNum=lineNum+1
    read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
    lineNum=lineNum+1
    read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
    lineNum=lineNum+1
    do while ( size1<6 .OR. cur_line(size1-1:size1).NE.')(' )
        read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
        lineNum=lineNum+1
    enddo
    !*******Read outlet_faces faces(nod1,node2,node3,node4,RightElement,Leftelement)********
    read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
    num_outlet_faces=0
    do while ( size1>6 .OR. cur_line(size1:size1).NE.')' )
        num_outlet_faces=num_outlet_faces+1
        read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
    enddo
    rewind(input)
    if(diagnostics_level.gt.1)write(*,*) 'Allocating outlet_faces : ', num_outlet_faces

    allocate(outlet_faces(num_outlet_faces,6))
    do ii=1, lineNum
        read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) no_data
    enddo
    do ii=1,num_outlet_faces
        read(input, *) hexaval1,hexaval2,hexaval3,hexaval4,hexaval5,hexaval6
        read(hexaval1, '(z4)') outlet_faces(ii,1)
        read(hexaval2, '(z4)') outlet_faces(ii,2)
        read(hexaval3, '(z4)') outlet_faces(ii,3)
        read(hexaval4, '(z4)') outlet_faces(ii,4)
        read(hexaval5, '(z4)') outlet_faces(ii,5)
        read(hexaval6, '(z4)') outlet_faces(ii,6)
        if(diagnostics_level.gt.1)write(*,*) 'outlet face : ',ii, outlet_faces(ii,:)
        lineNum=lineNum+1
    enddo
    !******* End reading outlet_faces faces *************************************************
    read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
    lineNum=lineNum+1
    read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
    lineNum=lineNum+1
    read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
    lineNum=lineNum+1

    do while ( size1<6 .OR. cur_line(size1-1:size1).NE.')(' )
        read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
        lineNum=lineNum+1
    enddo
    !*******Read wall_faces faces(nod1,node2,node3,node4,RightElement,Leftelement)********
    read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
    num_wall_faces=0
    do while ( size1>6 .OR. cur_line(size1:size1).NE.')' )
        num_wall_faces=num_wall_faces+1
        read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
    enddo
    rewind(input)
    if(diagnostics_level.gt.1)write(*,*) 'Allocating wall_faces : ', num_wall_faces

    allocate(wall_faces(num_wall_faces,6))
    do ii=1, lineNum
        read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) no_data

    enddo
    do ii=1,num_wall_faces
        read(input, *) hexaval1,hexaval2,hexaval3,hexaval4,hexaval5,hexaval6
        read(hexaval1, '(z4)') wall_faces(ii,1)
        read(hexaval2, '(z4)') wall_faces(ii,2)
        read(hexaval3, '(z4)') wall_faces(ii,3)
        read(hexaval4, '(z4)') wall_faces(ii,4)
        read(hexaval5, '(z4)') wall_faces(ii,5)
        read(hexaval6, '(z4)') wall_faces(ii,6)
        if(diagnostics_level.gt.1)write(*,*) 'wall face : ',ii, wall_faces(ii,:)
        lineNum=lineNum+1
    enddo
    !******* End reading wall_faces faces *************************************************

    close (input)

    num_all_faces = num_wall_faces + num_faces + num_inlet_faces + num_outlet_faces
    allocate(all_faces(num_all_faces,6))

    do i=1,num_faces
       all_faces(i,:)= internal_faces(i,:)
    enddo
    i = num_faces
    do j=1,num_inlet_faces
       all_faces(i+j,:)= inlet_faces(j,:)
    enddo
    j= num_inlet_faces
    do k=1,num_outlet_faces
       all_faces(i+j+k,:)= outlet_faces(k,:)
    enddo
    k = num_outlet_faces
    do m=1,num_wall_faces
       all_faces(i+j+k+m,:)= wall_faces(m,:)
    enddo


    call enter_exit(sub_name,2)
end subroutine read_icem_msh

subroutine read_k_file(filename)
    use arrays,only: dp,node_3d,elem_3d,num_nodes,num_elems
    use diagnostics, only: enter_exit,get_diagnostics_level
    use indices
    use other_consts, only: MAX_FILENAME_LEN
    implicit none
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_READ_K_FILE" :: READ_K_FILE

    character(len=MAX_FILENAME_LEN), intent(in) :: filename !Input nodefile

    character (100) :: cur_line, no_data
    integer :: input,status1,size1,lineNum,e
    integer :: input_stat,I, ii, jj

    character(len=60) :: sub_name
    integer :: diagnostics_level
    integer, allocatable :: ELEMENT(:,:)

    sub_name = 'read_k_file'
    call enter_exit(sub_name,1)
    call get_diagnostics_level(diagnostics_level)

    lineNum=0
    input=1
    if(diagnostics_level.gt.1)write(*,*) 'Opening file : ', filename
    open ( unit = input, file = filename, status = 'old', &
     iostat = input_stat )
    if ( input_stat .ne. 0 ) then
        write ( *, '(a)' ) ''
        write ( *, '(a)' ) 'KMSH_DATA_READ - Fatal error!'
        write ( *, '(a)' ) '  Could not open input file:' &
         // trim ( filename )
        stop 1
    end if

    read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
    lineNum=lineNum+1
    do while ( cur_line(1:size1).NE.'*NODE' )
        read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
        lineNum=lineNum+1
    end do

    read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
    lineNum=lineNum+1

    !*******Start reading the NODES(#,x,y,z)*******************************
    read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
    I=0
    do while ( cur_line(size1:size1).NE.'$' )
        I=I+1
        read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
        !print *, "!!", cur_line
    end do
    rewind(input)

    if(diagnostics_level.gt.1)write(*,*) 'Allocating nodes : ', I !double allocation?
    allocate(node_3d(I,4))
    do ii=1, lineNum
        read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) no_data
    enddo
    do ii=1,I
        read(input, *) node_3d(ii,:)
        if(diagnostics_level.gt.1)write(*,*) 'new node (double up?) : ',ii, node_3d(ii,:)
        lineNum=lineNum+1
    enddo
    !******* END reading the NODEs(#,x,y,z)***********************************

    read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
    lineNum=lineNum+1
    do while ( cur_line(1:size1).NE.'*ELEMENT_SOLID' )
        read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
        lineNum=lineNum+1
    end do

    !*******Start reading the ELEMENT(node1,node2...node8)*******************************
    read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
    I=0
    do while ( cur_line(1:size1).NE.'*ELEMENT_SHELL' )
        I=I+1
        read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
    end do

    rewind(input)
    num_elems = I/2
    if(diagnostics_level.gt.1)write(*,*) 'Allocating elems : ', I/2
    allocate(elem_3d(I/2,8)) !nodes associated with elements
    do ii=1, lineNum
        read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) no_data
    enddo
    jj=0
    do ii=1,I
    if (mod(ii,2)==0)then
        jj=jj+1
        read(input, *)  elem_3d(jj,:)
        if(diagnostics_level.gt.1)write(*,*) 'elem : ',jj, elem_3d(jj,:)
        lineNum=lineNum+1
    else
        read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) no_data
        lineNum=lineNum+1
    endif
    enddo
    !******* END reading the ELEMENTs(x,y,z)***********************************

    close (input)

    ! Some elements have nodes defined in wrong order - note this is not robust

    allocate(ELEMENT(size(elem_3d,1),8))
    do e=1,size(ELEMENT,1)
        if (e<=1800 .OR. e>2700) then
            ELEMENT(e,1)=elem_3d(e,1)
            ELEMENT(e,2)=elem_3d(e,5)
            ELEMENT(e,3)=elem_3d(e,6)
            ELEMENT(e,4)=elem_3d(e,2)
            ELEMENT(e,5)=elem_3d(e,4)
            ELEMENT(e,6)=elem_3d(e,8)
            ELEMENT(e,7)=elem_3d(e,7)
            ELEMENT(e,8)=elem_3d(e,3)
        else
            ELEMENT(e,1)=elem_3d(e,1)
            ELEMENT(e,2)=elem_3d(e,2)
            ELEMENT(e,3)=elem_3d(e,3)
            ELEMENT(e,4)=elem_3d(e,4)
            ELEMENT(e,5)=elem_3d(e,5)
            ELEMENT(e,6)=elem_3d(e,6)
            ELEMENT(e,7)=elem_3d(e,7)
            ELEMENT(e,8)=elem_3d(e,8)
        endif
    enddo
    elem_3d=ELEMENT
    deallocate(ELEMENT)

    write(*,*) filename
    call enter_exit(sub_name,2)

end subroutine read_k_file
end module geometry

