module abm_forces

  implicit none
  !Module parameters

  !Module types

  !Module depvar

  !Interfaces
  private
  public calc_random_forces
  public calc_saghian_chemo_forces
  public calc_saghian_cell_cell
contains
!
!###################################################################################
!
subroutine calc_random_forces(cell_population, force_field, force_magnitude)
    use arrays, only: dp,num_cells,cell_list,cell_stat,cell_field
    use diagnostics, only: enter_exit,get_diagnostics_level
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_CALC_RANDOM_FORCES" :: CALC_RANDOM_FORCES

    
    integer, intent(in) :: cell_population
    integer, intent(in) :: force_field
    real(dp), intent(in) :: force_magnitude
    !local variables
    integer :: kcell,dir
    integer :: diagnostics_level
    real(dp) :: direction(3)
    character(len=60) :: sub_name
    
    sub_name = 'calc_random_forces'
    call enter_exit(sub_name,1)
    call get_diagnostics_level(diagnostics_level)
    
    do kcell = 1,num_cells
      if(cell_list(kcell)%ctype.eq.cell_population.and.cell_list(kcell)%state.eq.cell_stat%ALIVE)then
        if(diagnostics_level.gt.1)then
          write(*,*) kcell
          write(*,*) cell_list(kcell)%ctype
        endif
        call get_random_dir(direction)
         if(diagnostics_level.gt.1)write(*,*) direction
        do dir = 1,3
          cell_field(kcell, force_field, dir) = force_magnitude*direction(dir)
        enddo
        if(diagnostics_level.gt.1)then
          write(*,*) cell_field(kcell, force_field, 1) ,cell_field(kcell, force_field, 2),cell_field(kcell, force_field, 3)
        endif
      endif
    enddo

    call enter_exit(sub_name,2)
end subroutine calc_random_forces
!
!###################################################################################
!
subroutine calc_saghian_chemo_forces(cell_population,force_field,fradial,cradial,faxial,caxial)
    use arrays, only: dp,num_cells,cell_list,cell_stat,cell_field,plug_params
    use other_consts, only: PI
    use diagnostics, only: enter_exit,get_diagnostics_level
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_CALC_SAGHIAN_CHEMO_FORCES" :: CALC_SAGHIAN_CHEMO_FORCES
    integer, intent(in) :: cell_population
    integer, intent(in) :: force_field
    real(dp), intent(in) :: fradial,cradial,faxial,caxial
    
    integer :: kcell
    real(dp) :: ccentre(3),d,r,theta
    integer :: diagnostics_level
    character(len=60) :: sub_name
    
    sub_name = 'calc_saghian_chemo_forces'
    call enter_exit(sub_name,1)
    call get_diagnostics_level(diagnostics_level)
    
    
    do kcell = 1,num_cells
      if(cell_list(kcell)%ctype.eq.cell_population.and.cell_list(kcell)%state.eq.cell_stat%ALIVE)then
        if(diagnostics_level.gt.1)then
          write(*,*) kcell
          write(*,*) cell_list(kcell)%ctype
        endif
        ccentre = cell_list(kcell)%centre(:,1)
        d = ccentre(3)/(plug_params%tube_length)
        r = sqrt(dot_product(ccentre,ccentre))/(plug_params%tube_length)
        theta=ATAN2(ccentre(2),ccentre(1))
        if (theta<0) then
            theta=theta+2*PI
        endif
    
        cell_field(kcell, force_field, :) = faxial *[0, 0, 1]*d*EXP(caxial*r) + &
            fradial *[1, 0, 0]*r*cos(theta)*exp(cradial*r) + &
            fradial *[0, 1, 0]*r*sin(theta)*exp(cradial*r)
        if(diagnostics_level.gt.1)then
          write(*,*) cell_field(kcell, force_field, :)
        endif
      endif 
    enddo

    call enter_exit(sub_name,2)
end subroutine calc_saghian_chemo_forces

!
!###################################################################################
!
subroutine calc_saghian_cell_cell(cell_population,force_field,r0,r1,a,b)
    use arrays, only: dp,num_cells,cell_list,cell_stat,cell_field,plug_params
    use other_consts, only: PI
    use math_utilities, only: unit_vector
    use diagnostics, only: enter_exit,get_diagnostics_level
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_CALC_SAGHIAN_CELL_CELL" :: CALC_SAGHIAN_CELL_CELL
    integer, intent(in) :: cell_population
    integer, intent(in) :: force_field
    real(dp), intent(in) :: r0,r1,a,b

    integer :: kcell,knbr,nbridx
    real(dp) :: dx, c, delta, xcross,F,x,Fdir(3),unitFdir(3),rad1,rad2

    integer :: diagnostics_level
    character(len=60) :: sub_name

    sub_name = 'calc_saghian_cell_cell'
    call enter_exit(sub_name,1)
    call get_diagnostics_level(diagnostics_level)

    dx = r1 - r0
    c = -b*8.0_dp - 4.0_dp*a/dx**2
    delta = dx**2.0_dp + 4.0_dp*a/c
    xcross = (r0+r1)/2.0_dp + 0.5_dp*sqrt(delta)

    do kcell = 1,num_cells
     !write(*,*) cell_list(kcell)%nbrs, cell_list(kcell)%nearest_dist,xcross
     do knbr = 1,cell_list(kcell)%nbrs
       nbridx = cell_list(kcell)%nbrlist(knbr)%indx
       Fdir =cell_list(nbridx)%centre(:,1)-cell_list(kcell)%centre(:,1)
       unitFdir = unit_vector(Fdir)
       rad1 = cell_list(kcell)%radius(1)
       rad2 = cell_list(nbridx)%radius(1)
       x= cell_list(kcell)%nbrlist(knbr)%distance/(rad1+rad2)
       if (x > xcross) then
         F = 0.0_dp
       elseif (x < r0) then ! should not happen with adaptive time stepping
         write(*,'(a,3e12.3,2f8.3,a,L2)') 'Error: get_force: x < x0: '!,R1,R2,d,x,x0_force,' incontact: ',incontact
         F = 0.0_dp
       else
         F = a/((x-r0)*(r1-x)) + c
       endif

       cell_field(kcell, force_field, :) = F*unitFdir
      ! write(*,*) cell_field(kcell, force_field, :),rad1,rad2,x

     enddo

    enddo

    call enter_exit(sub_name,2)
end subroutine calc_saghian_cell_cell
!
!###################################################################################
!
subroutine get_random_dir(dr)
  use arrays, only: dp
  use math_utilities, only: unit_vector
  real(dp) :: dr(3)
  real(dp) :: r
  integer :: i
  
  do i = 1,3
    call RANDOM_NUMBER(r)
    dr(i) = 2*(r - 0.5)
  enddo
  dr = unit_vector(dr)
 
end subroutine get_random_dir

end module abm_forces
