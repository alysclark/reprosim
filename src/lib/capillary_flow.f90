module capillary_flow

  implicit none
  !Module parameters

  !Module types

  !Module depvar

  !Interfaces
  private
  public export_capillary_pressures
contains


subroutine export_capillary_pressures(filename)

    use indices
    use arrays,only: dp,elem_field,node_field,num_units,units,unit_field,elem_nodes,elem_cnct
    use other_consts, only: PI
    use diagnostics, only: enter_exit,get_diagnostics_level

    character(len=60), intent(in) :: filename
    

    integer :: nu,ne,np1,nc,np2
    integer :: nv,i,j
    real(dp) :: int_radius_in,int_radius_out, int_radius,int_length,&
      seg_length,viscosity,cap_resistance,terminal_resistance,&
      cap_resistance2,resistance(6),seg_resistance,pressure_in,pressure_out,q
    character(len=60) :: sub_name
    integer :: diagnostics_level,num_convolutes,num_generations
    
    sub_name = 'export_capillary_pressures'
    call enter_exit(sub_name,1)
    call get_diagnostics_level(diagnostics_level)


      num_convolutes = 6
      num_generations = 3
	!num_units - number of terminal elements
    do nu=1,num_units    
      ne=units(nu)!get terminal element
      nc = elem_cnct(1,1,ne) !capillary unit is downstream of a terminal unit
      np1=elem_nodes(1,nc)    
      np2=elem_nodes(2,nc)
      nv =  elem_cnct(1,1,nc) !vein is downstream of the capillary
      write(*,*) ne,nc,nv
      int_radius_in = (elem_field(ne_radius,ne)+elem_field(ne_radius,nv))/2.0_dp ! mm radius of inlet intermediate villous (average of artery and vein)
      int_radius_out=(0.03_dp + 0.03_dp/2.0_dp)/2.0_dp ! mm radius of mature intermediate villous (average of artery and vein)
      int_length=1.5_dp !mm Length of each intermediate villous
      seg_length=int_length/num_convolutes !lengh of each intermediate villous segment
      viscosity=0.33600e-02_dp !Pa.s !viscosity: fluid viscosity
      cap_resistance = 50000
      terminal_resistance = 0.0_dp
      
      q = elem_field(ne_Qdot,ne)
      
      pressure_in = node_field(nj_bv_press,np1)
      pressure_out = pressure_in
      
      write(*,*) nu, node_field(nj_bv_press,np1)-node_field(nj_bv_press,np2),elem_field(ne_Qdot,ne) 
      
      do j=1,num_generations
      !arterial connection
      q=q/2.0_dp
      int_radius = int_radius_in - (int_radius_in-int_radius_out)/num_generations*j
      seg_resistance=(8.0_dp*viscosity*seg_length)/(PI*int_radius**4.0_dp) !resistance of each intermediate villous segment
      pressure_out = pressure_in - q*seg_resistance
      !calculate total resistance of terminal capillary conduits
      write(*,*) pressure_in, pressure_out
      i=1
      resistance(i)= cap_resistance + 2.0_dp*seg_resistance
      do i=2,num_convolutes
        resistance(i)=2.0_dp*seg_resistance + 1.0_dp/(1.0_dp/cap_resistance + 1.0_dp/resistance(i-1))
      enddo
      cap_resistance2 = resistance(num_convolutes) !Pa . s per mm^3 total resistance of terminal capillary conduits
	   write(*,*) resistance
      !We have symmetric generations of intermediate villous trees so we can calculate the total resistance
      !of the system by summing the resistance of each generation

      terminal_resistance = terminal_resistance + cap_resistance2/2**j
    enddo
          
    enddo
    call enter_exit(sub_name,2)
end subroutine export_capillary_pressures


end module capillary_flow