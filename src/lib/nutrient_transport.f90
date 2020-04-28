module nutrient_transport
!*Description:* This module contains tools that are used to solve systems of equations nutrient transport .
!
! Descriptions for subroutines that are not included in the subroutine:
  use indices
  use arrays
  implicit none
  !Module parameters

  !Module types

  !Module depvar

  !Interfaces
  private
  public evaluate_transport

contains

subroutine evaluate_transport()
    use diagnostics, only: enter_exit
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_EVALUATE_TRANSPORT" :: EVALUATE_TRANSPORT

    integer :: nu,noelem,nnod1,nnod2
    character(len=60) :: sub_name

    sub_name = 'evaluate_transport'

    call enter_exit(sub_name,1)

    do nu = 1,num_units
      noelem = units(nu)
      nnod1 = elem_nodes(1,noelem)
      nnod2 = elem_nodes(2,noelem)
      !write(*,*) units(nu),nnod1,nnod2,node_field(nnod1,nj_bv_press) !Get a terminal unit
      elem_field(ne_dp,noelem) = node_field(nnod1,nj_bv_press)  - node_field(nnod2,nj_bv_press)
      write(*,*) units(nu), elem_field(ne_dp,noelem),elem_field(ne_Qdot,noelem)
      !Total delta P across element !was 40Pa in paper
      call calc_capillary_pdrop(noelem,6,3,node_field(nnod1,nj_bv_press),node_field(nnod2,nj_bv_press),elem_field(ne_Qdot,noelem))

    enddo

    !Calculate Damkohler
    !Dt, L, R, B, dltaP
    !L is a length scale specific to the villus (basically the relating to diffusion distance, soln to Laplaces for a villous of theat type )
    !R is network resistance
    !B is a constant (models facilitated transport)
    !deltaP is pressure drop



    !Calculate mu
    !Dt, L, Dp, Lc



    call enter_exit(sub_name,2)
end subroutine evaluate_transport


!
!###################################################################################
!
  subroutine calc_capillary_pdrop(ne,num_convolutes,num_generations,press_in,press_out,flow)
  !*Description:*
    use diagnostics, only: enter_exit,get_diagnostics_level
    use other_consts, only: PI
    implicit none

    integer, intent(in) :: ne,num_convolutes,num_generations
    real(dp), intent(in) :: press_in,press_out,flow

    real(dp) :: int_length,int_radius,seg_length,viscosity, &
                seg_resistance,cap_unit_radius, cap_length,&
                capillary_unit_vol, capillary_unit_area
    real(dp) :: cap_resistance2
    real(dp) :: int_radius_in,int_radius_out
    real(dp),allocatable :: resistance(:)
    integer :: nu,i,j,np1,np2,nc,nv
    integer :: AllocateStatus
    real(dp) :: press_in_cap,press_out_cap,press_in_gen,press_out_gen,Qin_gen
    character(len=60) :: sub_name
    integer:: diagnostics_level

    sub_name = 'calc_capillary_pdrop'
    call enter_exit(sub_name,1)
    call get_diagnostics_level(diagnostics_level)

    num_conv = num_convolutes
    num_conv_gen = num_generations

    if (diagnostics_level.GE.2)then
      print *, "num_convolutes=",num_convolutes
      print *, "num_generations=",num_generations
    endif

    allocate (resistance(num_convolutes+1), STAT = AllocateStatus)
    if (AllocateStatus /= 0)then
       STOP "*** Not enough memory for resistance array ***"
    endif

    resistance = 0.0_dp


    nc = elem_cnct(1,1,ne) !capillary unit is downstream of a terminal unit
    nv =  elem_cnct(1,1,nc) !vein is downstream of the capillary
    int_radius_in = (elem_field(ne_radius,ne)+elem_field(ne_radius,nv))/2.0_dp ! mm radius of inlet intermediate villous (average of artery and vein)
    int_radius_out=(0.03_dp + 0.03_dp/2.0_dp)/2.0_dp ! mm radius of mature intermediate villous (average of artery and vein)
    int_length=1.5_dp !mm Length of each intermediate villous
    cap_length=3.0_dp/num_convolutes !mm length of capillary convolutes
    cap_radius=0.0144_dp/2.0_dp !radius of capillary convolutes
    seg_length=int_length/num_convolutes !lengh of each intermediate villous segment
    viscosity=0.33600e-02_dp !Pa.s !viscosity: fluid viscosity
    cap_unit_radius = 0.03_dp
    cap_resistance=(8.d0*viscosity*cap_length)/(PI*cap_radius**4) !resistance of each capillary convolute segment
    cap_resistance = 50000.0_dp
    terminal_resistance = 0.0_dp

    press_in_gen = press_in
    press_out_gen = press_out
    Qin_gen = flow

    do j=1,num_generations
      Qin_gen = Qin_gen/2.0_dp
      int_radius = int_radius_in - (int_radius_in-int_radius_out)/num_generations*j
      seg_resistance=(8.0_dp*viscosity*seg_length)/(PI*int_radius**4.0_dp) !resistance of each intermediate villous segment
      !capillary unit volume and surface area calculation - adding intermediate villous volume
      i=1
      press_in_cap =press_in -seg_resistance*Qin_gen
      press_out_cap = press_out + seg_resistance*Qin_gen
      write(*,*)ne,j,i,press_in_cap,press_out_cap,press_in_cap - press_out_cap,Qin_gen
      resistance(i)= cap_resistance + 2.0_dp*seg_resistance
      do i=2,num_convolutes
        Qin_gen = Qin_gen*seg_resistance/cap_resistance
        resistance(i)=2.0_dp*seg_resistance + 1.0_dp/(1.0_dp/cap_resistance + 1.0_dp/resistance(i-1))
        press_in_cap =press_in_cap -seg_resistance*Qin_gen
        press_out_cap = press_out_cap + seg_resistance*Qin_gen

        write(*,*)ne,j,i,press_in_cap,press_out_cap,press_in_cap-press_out_cap,Qin_gen
      enddo
      cap_resistance2 = resistance(num_convolutes) !Pa . s per mm^3 total resistance of terminal capillary conduits
      !We have symmetric generations of intermediate villous trees so we can calculate the total resistance
      !of the system by summing the resistance of each generation

      terminal_resistance = terminal_resistance + cap_resistance2/2**j
    enddo

    terminal_length = terminal_resistance*(PI*cap_unit_radius**4.0_dp)/(8.0_dp*viscosity)

    total_cap_volume = capillary_unit_vol * num_units/1000 !in cm**3
    total_cap_surface_area = capillary_unit_area * num_units/100 !in cm**2

    if(diagnostics_level.GE.2)then
      print *, "Resistance of capillary conduits=",cap_resistance
      print *, "Resistance of all generations of capillaries per terminal unit=",terminal_resistance
      print *, "Effective length of each capillary unit (mm)",terminal_length
      print *, "Total capillary length for the vasculature (cm)",(terminal_length*num_units)/10
      print *, "Total capillary volume (cm**3) = ",total_cap_volume
      print *, "Total capillary surface area (cm**2) = ", total_cap_surface_area
      write(*,*) 'overall cap resistance', terminal_resistance
      write(*,*) 'Effective length',terminal_length
    endif


    deallocate(resistance)
    call enter_exit(sub_name,2)

  end subroutine calc_capillary_pdrop

end module nutrient_transport
