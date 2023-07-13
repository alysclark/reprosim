module fetal
!*Description:* This module contains fetal models
! Descriptions for subroutines that are not included in the subroutine:

  use arrays
  use diagnostics
  use indices
  use other_consts

  implicit none
  !Module parameters

  !Module types

  !Module depvar

  !Interfaces
  private
  public fetal_model
  public assign_fetal_arrays

  real(dp),parameter,private :: T_beat = 0.43_dp         ! heart beat period (s)
  real(dp),parameter,private :: T_vs  = 0.215_dp ! Time period of ventricular contraction (s)
  real(dp), parameter,private :: T_as = 0.1075_dp !Time period of atrial contraction (s)
  real(dp), parameter,private :: T_v_delay = 0.1075_dp !delay in ventrial contraction (compare to atria) (s)
  real(dp), parameter, private :: U0RV = 5332.89_dp !Pa
  real(dp), parameter, private :: EsysRV = 0.399967_dp !Pa/mm3
  real(dp), parameter, private :: EdiaRV = 0.0399967_dp !Pa/mm3
  real(dp), parameter, private :: RvRV = 0.010665 !Pa.s/mm3
  real(dp), parameter, private :: U0LV = 5332.89_dp !Pa
  real(dp), parameter, private :: EsysLV = 0.399967_dp !Pa/mm3
  real(dp), parameter, private :: EdiaLV = 0.0399967_dp !Pa/mm3
  real(dp), parameter, private :: RvLV = 0.010665_dp !Pa.s/mm3
  real(dp), parameter, private :: U0A = 399.967_dp !Pa


contains
    subroutine fetal_model
        use diagnostics, only: enter_exit,get_diagnostics_level

    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_FETAL_MODEL" :: FETAL_MODEL

        real(dp) :: time                  !current time (s)
        real(dp) :: dt                     !Time step (s)
        integer :: num_heart_beats         !num heart beats
        real(dp) :: T_interval            ! the total length of the heat beat (s)

        integer :: n    !current heart beat
        real(dp) :: ttime !time within current heart beat
        integer :: np,ne,np_in,np_out

        real(dp) :: Avent !Ventricular activation (no units)
        real(dp) :: Aatria !Atrial activation (no units)
        real(dp) :: dpress,Pgrad,Qnod,Qnew,dQ,Vnod
        logical :: continue
        character(len=60) :: sub_name
        integer :: diagnostics_level

        !------
        sub_name = 'assign_fetal_arrays'
        call enter_exit(sub_name,1)
        call get_diagnostics_level(diagnostics_level)

        dt = 0.0001_dp
        num_heart_beats = 1
        T_interval = num_heart_beats * T_beat
        write(*,*) "simulating for" , T_interval, " s"


        !!! HARD CODED - TO BE READ IN AS FIELDS
        !Compartment 1 - Right ventricle
        node_field_fetal(njf_type,1) = 1.0_dp !The right ventricle
        node_field_fetal(njf_press,1) = 869.82_dp
        node_field_fetal(njf_comp,1) = 1.0_dp/EdiaRV


        !Compartment 2 - Left ventricle
        node_field_fetal(njf_type,2) = 2.0_dp !The right ventricle
        node_field_fetal(njf_press,2) = 869.82_dp
        node_field_fetal(njf_comp,2) = 1.0_dp/EdiaRV
        !Compartment 3 - Right atrium
        node_field_fetal(njf_type,3) = 3.0_dp !Any atrium
        node_field_fetal(njf_press,3) = 2.1_dp*133.0_dp
        node_field_fetal(njf_comp,3) = 2.0_dp*1000.0_dp/133.0_dp !ml/mmHg to mm3/Pa.
        !Compartment 4 - Left atrium
        node_field_fetal(njf_type,4) = 3.0_dp !Any atrium
        node_field_fetal(njf_press,4) = 3.0_dp*133.0_dp !m
        node_field_fetal(njf_comp,4) = 2.0_dp*1000.0_dp/133.0_dp !ml/mmHg to mm3/Pa.
        !Compartment 5 PA1
        node_field_fetal(njf_press,5) = 44.0_dp*133.0_dp !m
        node_field_fetal(njf_comp,5) = 0.08*1000.0_dp/133.0_dp !ml/mmHg to mm3/Pa
        do np = 5,num_nodes_fetal
            node_field_fetal(njf_type,np) = 4.0_dp !Generic compartment
        end do
        !AA
        node_field_fetal(njf_press,6) = 44.4_dp*133.0_dp !Pa
        node_field_fetal(njf_comp,6) = 0.37593984962406_dp !mm3/Pa
        !AO1
        node_field_fetal(njf_press,7) = 43.1_dp*133.0_dp !Pa
        node_field_fetal(njf_comp,7) = 0.601503759398496_dp
        !AO2
        node_field_fetal(njf_press,8) = 5506.2_dp!Pa
        node_field_fetal(njf_comp,8) = 0.526315789473684_dp
        !Ao3
        node_field_fetal(njf_press,9) = 5426.4_dp!Pa
        node_field_fetal(njf_comp,9) = 0.300751879699248_dp
        !AO4
        node_field_fetal(njf_press,10) = 5359.9_dp!Pa
        node_field_fetal(njf_comp,10) = 0.37593984962406_dp
        !PA2
        node_field_fetal(njf_press,11) = 5732.3_dp!Pa
        node_field_fetal(njf_comp,11) = 0.601503759398496_dp
        !Lung
        node_field_fetal(njf_press,12) = 1463_dp!Pa
        node_field_fetal(njf_comp,12) = 3.00751879699248_dp
        !CA
        node_field_fetal(njf_press,13) = 5599.3_dp!Pa
        node_field_fetal(njf_comp,13) = 0.075187969924812_dp
        !BR
        node_field_fetal(njf_press,14) = 4309.2_dp!Pa
        node_field_fetal(njf_comp,14) = 2.25563909774436_dp
        !SVC
        node_field_fetal(njf_press,15) = 625.1_dp!Pa
        node_field_fetal(njf_comp,15) = 7.5187969924812_dp
        !UB
        node_field_fetal(njf_press,16) = 2566.9_dp!Pa
        node_field_fetal(njf_comp,16) = 6.39097744360902_dp
        !HE
        node_field_fetal(njf_press,17) = 678.3_dp!Pa
        node_field_fetal(njf_comp,17) = 22.5563909774436_dp
        !INTE
        node_field_fetal(njf_press,18) = 1476.3_dp!Pa
        node_field_fetal(njf_comp,18) = 1.8796992481203_dp
        !KID
        node_field_fetal(njf_press,19) = 4269.3_dp!Pa
        node_field_fetal(njf_comp,19) =  0.150375939849624_dp
        !IVC
        node_field_fetal(njf_press,20) = 571.9!Pa
        node_field_fetal(njf_comp,20) = 4.51127819548872_dp
        !PLAC
        node_field_fetal(njf_press,21) = 2979.2_dp!Pa
        node_field_fetal(njf_comp,21) = 11.2781954887218_dp
        !UV
        node_field_fetal(njf_press,22) = 893.2598656_dp!Pa
        node_field_fetal(njf_comp,22) = 2.25563909774436_dp
        !LE
        node_field_fetal(njf_press,23) = 1356.6_dp!Pa
        node_field_fetal(njf_comp,23) = 30.0751879699248_dp !mm3/Pa


        do np = 1,num_nodes_fetal
            node_field_fetal(njf_vol,np) = node_field_fetal(njf_press,np)/node_field_fetal(njf_comp,np)
            write(*,*) np, node_field_fetal(njf_vol,np)
        end do




        !ELement 1, 1-5 RV-PA1, one way flow should occur when RV pressure > PA pressure
        elem_field_fetal(ne_group,1) = 1.0_dp !One way valve
        elem_field_fetal(ne_resist,1) = 0.0_dp
        elem_field_fetal(nef_K,1) = 0.001_dp *133.0_dp/(1000.0_dp*1000.0_dp)!mmHg s2/ml2 -> Pa.s2/mm6
        elem_field_fetal(nef_L,1) = 0.0_dp
        ! 2-6 LV-AA, one way flow, should occur when LV pressure > AA pressure
        elem_field_fetal(ne_group,2) = 1.0_dp !One way valve
        elem_field_fetal(ne_resist,2) = 0.0_dp
        elem_field_fetal(nef_K,2) = 0.001_dp *133.0_dp/(1000.0_dp*1000.0_dp)!mmHg s2/ml2 -> Pa.s2/mm6
        elem_field_fetal(nef_L,2) = 0.0_dp
        ! 3-15 RA-SVC, standard 2 2way flow
        elem_field_fetal(ne_group,3) = 2.0_dp !Simple R-Q unit
        elem_field_fetal(ne_resist,3) = 0.0266644736_dp ! Pa s /mm3
        elem_field_fetal(nef_K,3) = 0.0_dp
        elem_field_fetal(nef_L,3) = 0.0_dp
        !3-1 RA-RV !One way flow from atrium to ventricle
        elem_field_fetal(ne_group,4) = 1.0_dp !One way valve
        elem_field_fetal(ne_resist,4) = 0.0_dp
        elem_field_fetal(nef_K,4) = 0.002_dp *133.0_dp/(1000.0_dp*1000.0_dp)!mmHg s2/ml2 -> Pa.s2/mm6
        elem_field_fetal(nef_L,4) = 0.0016_dp*133.0_dp/1000.0_dp !mmHg s2/ml - Pa . s2/mm3

        !4-12 LA to lung, standard 2 way flow
        elem_field_fetal(ne_group,5) = 2.0_dp !Simple R-Q unit
        elem_field_fetal(ne_resist,5) = 0.266644736_dp ! Pa s /mm3
        elem_field_fetal(nef_K,5) = 0.0_dp
        elem_field_fetal(nef_L,5) = 0.0_dp

        !4-2 LA-LV One way flow from atrium to ventricle
        elem_field_fetal(ne_group,6) = 1.0_dp !One way valve
        elem_field_fetal(ne_resist,6) = 0.0_dp
        elem_field_fetal(nef_K,6) = 0.002_dp *133.0_dp/(1000.0_dp*1000.0_dp)!mmHg s2/ml2 -> Pa.s2/mm6
        elem_field_fetal(nef_L,6) = 0.0016_dp*133.0_dp/1000.0_dp !mmHg s2/ml - Pa . s2/mm3

        !5-11 PA1-PA2, standard 2 way flow
         elem_field_fetal(ne_group,7) = 2.0_dp !R_Q-L-unit
        elem_field_fetal(ne_resist,7) = 0.00933256576_dp ! Pa s /mm3 RPA
        elem_field_fetal(nef_K,7) = 0.0_dp
        elem_field_fetal(nef_L,7) = 0.002_dp*133.0_dp/1000.0_dp

        !6-7 AA-AO1, standard 2 way flow, not simple R-Q
        elem_field_fetal(ne_group,8) = 3.0_dp !R-Q-L unit
        elem_field_fetal(ne_resist,8) = 0.01599868416_dp ! Pa s /mm3 !AA Resistance
        elem_field_fetal(nef_K,8) = 0.0_dp
        elem_field_fetal(nef_L,8) = 0.002_dp*133.0_dp/1000.0_dp !mmHg s2/ml - Pa . s2/mm3

        !7-8 AO1-AO2 standard 2 way flow
        elem_field_fetal(ne_group,9) = 2.0_dp !Simple R-Q unit
        elem_field_fetal(ne_resist,9) =0.0533289472_dp ! Pa s /mm3 RIsthm
        elem_field_fetal(nef_K,9) = 0.0_dp
        elem_field_fetal(nef_L,9) = 0.0_dp

        !7-13 AO1-CA, CARO
        elem_field_fetal(ne_group,10) = 3.0_dp !R-Q-L unit
        elem_field_fetal(ne_resist,10) = 0.0399967104_dp ! Pa s /mm3
        elem_field_fetal(nef_K,10) = 0.0_dp
        elem_field_fetal(nef_L,10) = 0.08_dp*133.0_dp/1000.0_dp

        !7-16 AO1-UB - RUBA
        elem_field_fetal(ne_group,11) = 2.0_dp !Simple R-Q unit
        elem_field_fetal(ne_resist,11) = 1.066578944_dp ! Pa s /mm3
        elem_field_fetal(nef_K,11) = 0.0_dp
        elem_field_fetal(nef_L,11) = 0.0_dp

        !8-9 AO2-AO3,
        elem_field_fetal(ne_group,12) = 2.0_dp !Simple R-Q unit
        elem_field_fetal(ne_resist,12) = 0.00533289472_dp ! Pa s /mm3 R_DTAO
        elem_field_fetal(nef_K,12) = 0.0_dp
        elem_field_fetal(nef_L,12) = 0.0_dp

        !8-11 AO2-PA2, PA2-AO2, Ductus arteriosus
        elem_field_fetal(ne_group,13) = 6.0_dp !R-Q-K-L unit, beta=2
        elem_field_fetal(ne_resist,13) = 0.00133322368_dp ! Pa s /mm3
        elem_field_fetal(nef_K,13) = 0.0009_dp*133.0_dp/(1000.0_dp*1000.0_dp)!
        elem_field_fetal(nef_L,13) = 0.006_dp*133.0_dp/1000.0_dp

        !9-10 AO3-AO4
        elem_field_fetal(ne_group,14) = 2.0_dp !Simple R-Q unit
        elem_field_fetal(ne_resist,14) = 0.00799934208_dp ! Pa s /mm3 !R_DAAO
        elem_field_fetal(nef_K,14) = 0.0_dp
        elem_field_fetal(nef_L,14) = 0.0_dp

        !9-17 AO3-He FIG
        elem_field_fetal(ne_group,15) = 2.0_dp !Simple R-Q unit
        elem_field_fetal(ne_resist,15) = 10.799111808_dp ! Pa s /mm3
        elem_field_fetal(nef_K,15) = 0.0_dp
        elem_field_fetal(nef_L,15) = 0.0_dp

        !9-18 AO3-Inte MEA
        elem_field_fetal(ne_group,16) = 2.0_dp !Simple R-Q unit
        elem_field_fetal(ne_resist,16) = 4.532960512_dp ! Pa s /mm3
        elem_field_fetal(nef_K,16) = 0.0_dp
        elem_field_fetal(nef_L,16) = 0.0_dp

        !9-19 AO3-Kid REA
        elem_field_fetal(ne_group,17) = 2.0_dp !Simple R-Q unit
        elem_field_fetal(ne_resist,17) = 0.466628288_dp ! Pa s /mm3
        elem_field_fetal(nef_K,17) = 0.0_dp
        elem_field_fetal(nef_L,17) = 0.0_dp

        !10-21 AO4-Plac UV
        elem_field_fetal(ne_group,18) = 2.0_dp !Simple R-Q unit
        elem_field_fetal(ne_resist,18) = 0.5199572352_dp ! Pa s /mm3
        elem_field_fetal(nef_K,18) = 0.0_dp
        elem_field_fetal(nef_L,18) = 0.0_dp

        !10-23 AO4-Leg Fa
        elem_field_fetal(ne_group,19) = 2.0_dp !Simple R-Q unit
        elem_field_fetal(ne_resist,19) = 0.466628288_dp ! Pa s /mm3
        elem_field_fetal(nef_K,19) = 0.0_dp
        elem_field_fetal(nef_L,19) = 0.0_dp

        !11-12 PA2-Lung,
        elem_field_fetal(ne_group,20) = 2.0_dp !Simple R-Q unit
        elem_field_fetal(ne_resist,20) = 1.799851968_dp ! Pa s /mm3 !Rlung
        elem_field_fetal(nef_K,20) = 0.0_dp
        elem_field_fetal(nef_L,20) = 0.0_dp

        !12-4 One way lung-LA (LA)
        elem_field_fetal(ne_group,21) = 2.0_dp !Simple R-Q unit
        elem_field_fetal(ne_resist,21) = 0.266644736_dp ! Pa s /mm3 !
        elem_field_fetal(nef_K,21) = 0.0_dp
        elem_field_fetal(nef_L,21) = 0.0_dp


        !13-14 CA-BR
        elem_field_fetal(ne_group,22) = 2.0_dp !Simple R-Q unit
        elem_field_fetal(ne_resist,22) =0.399967104_dp ! Pa s /mm3 !RMCA
        elem_field_fetal(nef_K,22) = 0.0_dp
        elem_field_fetal(nef_L,22) = 0.0_dp

        !14-15 BR-SVC
        elem_field_fetal(ne_group,23) = 2.0_dp !Simple R-Q unit
        elem_field_fetal(ne_resist,23) =1.133240128_dp ! Pa s /mm3 !RBR
        elem_field_fetal(nef_K,23) = 0.0_dp
        elem_field_fetal(nef_L,23) = 0.0_dp

        !15-16 SVC-UB - UBV
        elem_field_fetal(ne_group,24) = 2.0_dp !Simple R-Q unit
        elem_field_fetal(ne_resist,24) =0.6532796032_dp ! Pa s /mm3
        elem_field_fetal(nef_K,24) = 0.0_dp
        elem_field_fetal(nef_L,24) = 0.0_dp


        !15-3 one way SVC to RA
        elem_field_fetal(ne_group,25) = 2.0_dp !Simple R-Q unit
        elem_field_fetal(ne_resist,25) =0.0266644736_dp ! Pa s /mm3
        elem_field_fetal(nef_K,25) = 0.0_dp
        elem_field_fetal(nef_L,25) = 0.0_dp

        !17-18 He-Inte PORV
        elem_field_fetal(ne_group,26) = 2.0_dp !Simple R-Q unit
        elem_field_fetal(ne_resist,26) =0.933256576_dp ! Pa s /mm3
        elem_field_fetal(nef_K,26) = 0.0_dp
        elem_field_fetal(nef_L,26) = 0.0_dp

        !17-20 He-IVC HV
        elem_field_fetal(ne_group,27) = 2.0_dp !Simple R-Q unit
        elem_field_fetal(ne_resist,27) =0.02133157888_dp ! Pa s /mm3
        elem_field_fetal(nef_K,27) = 0.0_dp
        elem_field_fetal(nef_L,27) = 0.0_dp

        !19-20 Kid-IVC REV
        elem_field_fetal(ne_group,28) = 2.0_dp !Simple R-Q unit
        elem_field_fetal(ne_resist,28) =1.866513152_dp ! Pa s /mm3
        elem_field_fetal(nef_K,28) = 0.0_dp
        elem_field_fetal(nef_L,28) = 0.0_dp

        !19-3 Kid to RA, one way

        !20-23 IVC-Leg FV
        elem_field_fetal(ne_group,30) = 2.0_dp !Simple R-Q unit
        elem_field_fetal(ne_resist,30) =0.0799934208_dp ! Pa s /mm3
        elem_field_fetal(nef_K,30) = 0.0_dp
        elem_field_fetal(nef_L,30) = 0.0_dp

        !20-22 IVC-UV DV
        elem_field_fetal(ne_group,31) = 5.0_dp !R-K-Q unit
        elem_field_fetal(ne_resist,31) =0.1733190784_dp ! Pa s /mm3
        elem_field_fetal(nef_K,31) = 0.0_DP!0.26_dp*133.0_dp/(1000.0_dp*1000.0_dp)!
        elem_field_fetal(nef_L,31) = 0.0_dp

        !20-4 IVC-LA,LA-IVC R-K-Q unit, FO beta special
        elem_field_fetal(ne_group,32) = 4.0_dp !R-K-Q unit
        elem_field_fetal(ne_resist,32) =0.1733190784_dp ! Pa s /mm3
        elem_field_fetal(nef_K,32) = 0.4_dp*133.0_dp/(1000.0_dp*1000.0_dp)!
        elem_field_fetal(nef_L,32) = 0.0_dp


        !21-22 Plac-UV RPLAC
        elem_field_fetal(ne_group,33) = 2.0_dp !Simple R-Q unit
        elem_field_fetal(ne_resist,33) =0.4532960512_dp ! Pa s /mm3
        elem_field_fetal(nef_K,33) = 0.0_dp
        elem_field_fetal(nef_L,33) = 0.0_dp

        !22-17 UV-He - HA
        elem_field_fetal(ne_group,34) = 2.0_dp !Simple R-Q unit
        elem_field_fetal(ne_resist,34) =0.066661184_dp ! Pa s /mm3
        elem_field_fetal(nef_K,34) = 0.0_dp
        elem_field_fetal(nef_L,34) = 0.0_dp

        !3-20 RA-IVC
        elem_field_fetal(ne_group,35) = 2.0_dp !Simple R-Q unit
        elem_field_fetal(ne_resist,35) =0.01599868416_dp ! Pa s /mm3
        elem_field_fetal(nef_K,35) = 0.0_dp
        elem_field_fetal(nef_L,35) = 0.0_dp
        !



        !print *,  node_field_fetal(njf_press,:)

        time = 0.0_dp !initialise the simulation time.
        open(10, file='results_volume.out', status='replace')
        open(20, file='results_pressure.out', status='replace')
        open(30, file='results_flow.out', status='replace')
        continue = .true.
        n = 0
        do while (continue)
            n = n + 1 ! increment the heart beat number
            ttime = 0.0_dp ! each breath starts with ttime=0
            !endtime = T_interval * n - 0.5_dp * dt ! the end time of this breath
            do while (ttime.lt.T_beat)
                ttime = ttime + dt ! increment the heartbeat time
                time = time + dt ! increment the whole simulation time
                if ((ttime.ge.T_v_delay).and.(ttime.le.T_vs+T_v_delay)) then
                    Avent = sin(pi/T_vs*(ttime-T_v_delay))
                else
                    Avent = 0.0_dp
                end if
                if ((ttime.le.T_as))then
                    Aatria = sin(pi/T_as*(ttime))
                else
                    Aatria = 0.0_dp
                end if

                do ne = 1, num_elems_fetal
                     np_in = elem_nodes_fetal(1,ne)
                     np_out = elem_nodes_fetal(2,ne)
                     Pgrad = node_field_fetal(njf_press,np_in)-node_field_fetal(njf_press,np_out)
                    if(elem_field_fetal(ne_group,ne).eq.1.0_dp)then!One way valve
                        !write(*,*) np_in,np_out,Pgrad
                        call one_way_valve(dt,dQ,elem_field_fetal(ne_Qdot,ne),Pgrad,elem_field_fetal(ne_resist,ne),&
                        elem_field_fetal(nef_K,ne),elem_field_fetal(nef_L,ne))

                        !print *, ne,dt,dQ,elem_field_fetal(ne_Qdot,ne)
                    else if (elem_field_fetal(ne_group,ne).eq.2.0_dp)then!R-Q unit
                        call rq_unit(dt,elem_field_fetal(ne_Qdot,ne),Pgrad,elem_field_fetal(ne_resist,ne))
                    else if (elem_field_fetal(ne_group,ne).eq.3.0_dp)then!R-Q-L unit
                        call rql_unit(dt,dQ,elem_field_fetal(ne_Qdot,ne),Pgrad,elem_field_fetal(ne_resist,ne),&
                                elem_field_fetal(nef_L,ne))
                    else if (elem_field_fetal(ne_group,ne).eq.4.0_dp)then!R-Q-K unit, FO
                           call rqk_unit(dt,elem_field_fetal(ne_Qdot,ne),Pgrad,elem_field_fetal(ne_resist,ne),&
                           elem_field_fetal(nef_K,ne), 0.625_dp)
                        if (elem_field_fetal(ne_Qdot,ne).le.0.0_dp) elem_field_fetal(ne_Qdot,ne)=0.0_dp
                    else if (elem_field_fetal(ne_group,ne).eq.5.0_dp)then!R-Q-K unit, DA and DV
                           call rqk_unit(dt,elem_field_fetal(ne_Qdot,ne),Pgrad,elem_field_fetal(ne_resist,ne),&
                           elem_field_fetal(nef_K,ne), 2.0_dp)
                    else if (elem_field_fetal(ne_group,ne).eq.6.0_dp)then!R-Q-K unit, DA and DV
                           call rqkl_unit(dt,dQ,elem_field_fetal(ne_Qdot,ne),Pgrad,elem_field_fetal(ne_resist,ne),&
                           elem_field_fetal(nef_K,ne), 2.0_dp,elem_field_fetal(nef_L,ne))
                    end if
                    !write(*,*) ne, np_in,np_out,elem_field_fetal(ne_Qdot,ne),elem_field_fetal(ne_group,ne),Pgrad
                end do

                do np =1,num_nodes_fetal
                    Qnod = 0.0_dp !Net flow passing through the node
                    do ne = 1, elems_at_node_fetal(np,0)
                        if (elem_nodes_fetal(1,elems_at_node_fetal(np,ne)).eq.np)then
                            !print *, 'first node of element, +ve flow leaves the node'
                            !#print *, elem_field_fetal(ne_Qdot,elems_at_node_fetal(np,ne))
                            Qnod = Qnod-elem_field_fetal(ne_Qdot,elems_at_node_fetal(np,ne))
                        else
                            !print *, 'second node of element, +ve flow enters the node'
                            !print *, elem_field_fetal(ne_Qdot,elems_at_node_fetal(np,ne))
                            Qnod = Qnod+elem_field_fetal(ne_Qdot,elems_at_node_fetal(np,ne))
                        end if
                   end do
                    dQ = Qnod - node_field_fetal(njf_netQ,np)
                    node_field_fetal(njf_netQ,np) = Qnod
                    !write(*,*) np,Qnod, dQ

                    Vnod = node_field_fetal(njf_vol,np)
                !
                    if(node_field_fetal(njf_type,np).eq.1.0_dp)then!right ventricle
                        call ventricle_pressure_step(dpress,dt,Avent,U0RV,EsysRV,EdiaRV,RvRV,Qnod,dQ,Vnod)
                        node_field_fetal(njf_press,np) = node_field_fetal(njf_press,np) + dpress
                    elseif(node_field_fetal(njf_type,np).eq.2.0_dp)then!left ventricle)
                        call ventricle_pressure_step(dpress,dt,Avent,U0LV,EsysLV,EdiaLV,RvLV,Qnod,dQ,Vnod)
                        node_field_fetal(njf_press,np) = node_field_fetal(njf_press,np) + dpress
                    elseif(node_field_fetal(njf_type,np).eq.3.0_dp)then!Its an atrium
                        call atrium_pressure_step(dpress,dt,Aatria,U0A,node_field_fetal(njf_comp,np),Qnod,dQ,Vnod)
                        node_field_fetal(njf_press,np) = node_field_fetal(njf_press,np) + dpress
                    else ! This is a standard node
                        call compartment_pressure_step(dpress,dt,node_field_fetal(njf_comp,np), Qnod,Vnod)
                        node_field_fetal(njf_press,np) = node_field_fetal(njf_press,np) + dpress
                    end if
                    node_field_fetal(njf_vol,np) = node_field_fetal(njf_vol,np)+dQ*dt
                enddo


                WRITE(10,'(27(F15.4,X))')&
                time, ttime,Avent,Aatria,node_field_fetal(njf_vol,1),&
                        node_field_fetal(njf_vol,2), node_field_fetal(njf_vol,3),node_field_fetal(njf_vol,4),&
                node_field_fetal(njf_vol,5), node_field_fetal(njf_vol,6),node_field_fetal(njf_vol,7),&
                node_field_fetal(njf_vol,8), node_field_fetal(njf_vol,9),node_field_fetal(njf_vol,10),&
                node_field_fetal(njf_vol,11), node_field_fetal(njf_vol,12),node_field_fetal(njf_vol,13),&
                node_field_fetal(njf_vol,14), node_field_fetal(njf_vol,15),node_field_fetal(njf_vol,16),&
                node_field_fetal(njf_vol,17), node_field_fetal(njf_vol,18),node_field_fetal(njf_vol,19),&
                node_field_fetal(njf_vol,20), node_field_fetal(njf_vol,21),node_field_fetal(njf_vol,22),&
                node_field_fetal(njf_vol,23)

                WRITE(20,'(27(F15.4,X))')&
                time, ttime,Avent,Aatria,node_field_fetal(njf_press,1),&
                        node_field_fetal(njf_press,2), node_field_fetal(njf_press,3),node_field_fetal(njf_press,4),&
                node_field_fetal(njf_press,5), node_field_fetal(njf_press,6),node_field_fetal(njf_press,7),&
                node_field_fetal(njf_press,8), node_field_fetal(njf_press,9),node_field_fetal(njf_press,10),&
                node_field_fetal(njf_press,11), node_field_fetal(njf_press,12),node_field_fetal(njf_press,13),&
                node_field_fetal(njf_press,14), node_field_fetal(njf_press,15),node_field_fetal(njf_press,16),&
                node_field_fetal(njf_press,17), node_field_fetal(njf_press,18),node_field_fetal(njf_press,19),&
                node_field_fetal(njf_press,20), node_field_fetal(njf_press,21),node_field_fetal(njf_press,22),&
                node_field_fetal(njf_press,23)

                WRITE(30,'(27(F15.4,X))')&
                time, ttime,Avent,Aatria,node_field_fetal(njf_netQ,1),&
                        node_field_fetal(njf_netQ,2), node_field_fetal(njf_netQ,3),node_field_fetal(njf_netQ,4),&
                node_field_fetal(njf_netQ,5), node_field_fetal(njf_netQ,6),node_field_fetal(njf_netQ,7),&
                node_field_fetal(njf_netQ,8), node_field_fetal(njf_netQ,9),node_field_fetal(njf_netQ,10),&
                node_field_fetal(njf_netQ,11), node_field_fetal(njf_netQ,12),node_field_fetal(njf_netQ,13),&
                node_field_fetal(njf_netQ,14), node_field_fetal(njf_netQ,15),node_field_fetal(njf_netQ,16),&
                node_field_fetal(njf_netQ,17), node_field_fetal(njf_netQ,18),node_field_fetal(njf_netQ,19),&
                node_field_fetal(njf_netQ,20), node_field_fetal(njf_netQ,21),node_field_fetal(njf_netQ,22),&
                node_field_fetal(njf_netQ,23)

            end do
            !'(8(F10.4),8(F10.2))')
            !if (n.eq.num_heart_beats) then


            continue = .false.
            !endif
        end do
        close(10)
        close(20)
        close(30)


        call enter_exit(sub_name,2)
    end subroutine fetal_model


    subroutine assign_fetal_arrays

        use arrays,only: dp,elem_field_fetal, num_elems_fetal, elem_nodes_fetal, nodes_fetal, elems_fetal, num_nodes_fetal,&
          node_field_fetal,node_xyz_fetal, elem_cnct_fetal,elem_direction_fetal,elems_at_node_fetal,elem_field, num_elems,&
          elem_nodes, nodes, elems, num_nodes, node_field, elem_cnct,elem_direction,elems_at_node,node_xyz
        use diagnostics, only: enter_exit,get_diagnostics_level

    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_ASSIGN_FETAL_ARRAYS" :: ASSIGN FETAL_ARRAYS

        character(len=60) :: sub_name
        integer :: diagnostics_level

        !------
        sub_name = 'assign_fetal_arrays'
        call enter_exit(sub_name,1)
        call get_diagnostics_level(diagnostics_level)

        num_nodes_fetal = num_nodes
        if(allocated(nodes_fetal)) deallocate (nodes_fetal)
        allocate (nodes_fetal(num_nodes_fetal))
        nodes_fetal = nodes
        if(allocated(node_xyz_fetal)) deallocate (node_xyz_fetal)
        allocate (node_xyz_fetal(3,num_nodes_fetal))
        node_xyz_fetal = node_xyz
        if(allocated(node_field_fetal)) deallocate (node_field_fetal)
        allocate (node_field_fetal(num_nj_fetal,num_nodes_fetal))
        node_field_fetal = node_field_fetal
        num_elems_fetal = num_elems
        if(allocated(elems_fetal)) deallocate(elems_fetal) !Array that defines nodal connections between elements
        allocate(elems_fetal(num_elems_fetal))
        elems_fetal = elems
        if(allocated(elem_cnct_fetal)) deallocate(elem_cnct_fetal) !Array that defines connections between elements
        allocate(elem_cnct_fetal(-1:1,0:10,0:num_elems_fetal))!Allows up to 10 elements per node
        elem_cnct_fetal = elem_cnct
        if(allocated(elem_nodes_fetal)) deallocate(elem_nodes_fetal)
        allocate(elem_nodes_fetal(2,num_elems_fetal)) !defines in and out nodes at each element
        elem_nodes_fetal = elem_nodes
        if(allocated(elems_at_node_fetal)) deallocate(elems_at_node_fetal)
        allocate(elems_at_node_fetal(num_nodes_fetal,0:10)) !Allows up to 10 elements per node
        elems_at_node_fetal = elems_at_node
        if(allocated(elem_field_fetal)) deallocate(elem_field_fetal)
        allocate(elem_field_fetal(num_ne,num_elems_fetal))
        elem_field_fetal = elem_field
        if(allocated(elem_direction_fetal)) deallocate(elem_direction_fetal)
        allocate(elem_direction_fetal(3,num_elems_fetal))
        elem_direction_fetal = elem_direction



        if(allocated(nodes)) deallocate (nodes)
        if(allocated(node_xyz)) deallocate (node_xyz)
        if(allocated(node_field)) deallocate (node_field)
        if(allocated(elems)) deallocate(elems)
        if(allocated(elem_cnct)) deallocate(elem_cnct)
        if(allocated(elem_nodes)) deallocate(elem_nodes)
        if(allocated(elems_at_node)) deallocate(elems_at_node)
        if(allocated(elem_field)) deallocate(elem_field)
        if(allocated(elem_direction)) deallocate(elem_direction)

        call enter_exit(sub_name,2)

    end subroutine assign_fetal_arrays

    subroutine ventricle_pressure_step(dpress,dt,Avent,U0,Edia,Esys,Rv,Q,dQ,V)
        use diagnostics, only: enter_exit,get_diagnostics_level

    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_VENTRICLE_PRESSURE_STEP" :: VENTRICLE_PRESSURE_STEP
        real(dp), intent(out) :: dpress
        real(dp), intent(in) :: dt
        real(dp), intent(in) :: Avent
        real(dp), intent(in) :: U0
        real(dp), intent(in) :: Edia
        real(dp), intent(in) :: Esys
        real(dp), intent(in) :: Rv
        real(dp), intent(in) :: Q
        real(dp), intent(in) :: dQ
        real(dp), intent(in) :: V

        character(len=60) :: sub_name
        integer :: diagnostics_level

        !------
        sub_name = 'ventricle pressure step'
        call enter_exit(sub_name,1)
        call get_diagnostics_level(diagnostics_level)

        dpress = dt*(U0*Avent + (Edia + Esys*Avent)*V + Rv*Q)
        call enter_exit(sub_name,2)

    end subroutine ventricle_pressure_step

subroutine atrium_pressure_step(dpress,dt,Aatria,U0,comp,Q,dQ,V)
        use diagnostics, only: enter_exit,get_diagnostics_level

    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_ATRIUM_PRESSURE_STEP" :: ATRIUM_PRESSURE_STEP
        real(dp), intent(out) :: dpress
        real(dp), intent(in) :: dt
        real(dp), intent(in) :: Aatria
        real(dp), intent(in) :: U0
        real(dp), intent(in) :: comp
        real(dp), intent(in) :: Q
        real(dp), intent(in) :: dQ
        real(dp), intent(in) :: V
        character(len=60) :: sub_name
        integer :: diagnostics_level

        !------
        sub_name = 'atrium pressure step'
        call enter_exit(sub_name,1)
        call get_diagnostics_level(diagnostics_level)

        dpress = dt*U0!(U0*Aatria + V/(comp))
        call enter_exit(sub_name,2)

    end subroutine atrium_pressure_step


subroutine compartment_pressure_step(dpress,dt,comp,Q,V)
        use diagnostics, only: enter_exit,get_diagnostics_level

    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_COMPARTMENT_PRESSURE_STEP" :: COMPARTMENT_PRESSURE_STEP
        real(dp), intent(out) :: dpress
        real(dp), intent(in) :: dt
        real(dp), intent(in) :: comp
        real(dp), intent(in) :: Q
        real(dp), intent(in) :: V

        character(len=60) :: sub_name
        integer :: diagnostics_level

        !------
        sub_name = 'compartment pressure step'
        call enter_exit(sub_name,1)
        call get_diagnostics_level(diagnostics_level)

        dpress = dt*Q/comp
        call enter_exit(sub_name,2)

    end subroutine compartment_pressure_step

subroutine one_way_valve(dt,dQ,Q,Pgrad, R, K, L)
        use diagnostics, only: enter_exit,get_diagnostics_level

    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_ONE_WAY_VALVE :: ONE_WAY_VALVE
        real(dp), intent(in) :: dt
        real(dp), intent(inout) :: dQ
        real(dp), intent(inout) :: Q
        real(dp), intent(in) :: Pgrad
        real(dp), intent(in) :: R
        real(dp), intent(in) :: K
        real(dp), intent(in) :: L

        real(dp) :: check_sign
        character(len=60) :: sub_name
        integer :: diagnostics_level

        !------
        sub_name = 'one_way_valve'
        call enter_exit(sub_name,1)
        call get_diagnostics_level(diagnostics_level)
        !write(*,*) Q
        if(L.gt.0)then
            dQ = dt*(Pgrad-K*Q**2.0_dp)/L
            Q = Q+dQ
        else
            if(Pgrad.gt.0) then
                Q = sqrt(Pgrad/K)
            else
                !write(*,*) 'doing this one'
                Q = -1.0_dp* sqrt(abs(Pgrad)/K)
            end if
        end if
        !write(*,*) Q,Pgrad,abs(Pgrad),sqrt(abs(Pgrad)/K),K
        if(Q.lt.0.0_dp)then
            Q=0.0_dp
        end if


        call enter_exit(sub_name,2)

    end subroutine one_way_valve

    subroutine rq_unit(dt,Q,Pgrad, R)
        use diagnostics, only: enter_exit,get_diagnostics_level

    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_RQ_UNIT :: RQ_UNIT
        real(dp), intent(out) :: dt
        real(dp), intent(inout) :: Q
        real(dp), intent(in) :: Pgrad
        real(dp), intent(in) :: R


        real(dp) :: check_sign
        character(len=60) :: sub_name
        integer :: diagnostics_level

        !------
        sub_name = 'rq_unit'
        call enter_exit(sub_name,1)
        call get_diagnostics_level(diagnostics_level)

        Q = Pgrad/R

        call enter_exit(sub_name,2)

    end subroutine rq_unit

        subroutine rql_unit(dt,dQ,Q,Pgrad, R,L)
        use diagnostics, only: enter_exit,get_diagnostics_level

    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_RQL_unit :: RQL_UNIT
        real(dp), intent(out) :: dt
        real(dp), intent(inout) :: dQ
        real(dp), intent(inout) :: Q
        real(dp), intent(in) :: Pgrad
        real(dp), intent(in) :: R
        real(dp), intent(in) :: L


        real(dp) :: check_sign
        character(len=60) :: sub_name
        integer :: diagnostics_level

        !------
        sub_name = 'rql_unit'
        call enter_exit(sub_name,1)
        call get_diagnostics_level(diagnostics_level)

        dQ = dt*(Pgrad - R*Q)/L

        Q = Q+dQ


        call enter_exit(sub_name,2)

    end subroutine rql_unit

       subroutine rqk_unit(dt,Q,Pgrad, R,K,beta)
        use diagnostics, only: enter_exit,get_diagnostics_level

    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_RQK_UNIT :: RQK_UNIT
        real(dp), intent(out) :: dt
        real(dp), intent(inout) :: Q
        real(dp), intent(in) :: Pgrad
        real(dp), intent(in) :: R
        real(dp), intent(in) :: K
        real(dp), intent(in) :: beta


        real(dp) :: check_sign
        character(len=60) :: sub_name
        integer :: diagnostics_level

        !------
        sub_name = 'rqk_unit'
        call enter_exit(sub_name,1)
        call get_diagnostics_level(diagnostics_level)


        write(*,*) Pgrad, R, K, Q, beta, Q**beta, Q**2.0_dp

        Q = Pgrad/R - K*Q**beta/R


        write(*,*) Pgrad, R, K, Q, beta, Q**beta, 1e-7**beta

        call enter_exit(sub_name,2)

        end subroutine rqk_unit

       subroutine rqkl_unit(dt,dQ,Q,Pgrad, R,K,beta,L)
        use diagnostics, only: enter_exit,get_diagnostics_level

    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_RQKL_UNIT :: RQKL_UNIT
        real(dp), intent(out) :: dt
        real(dp), intent(inout) :: dQ
        real(dp), intent(inout) :: Q
        real(dp), intent(in) :: Pgrad
        real(dp), intent(in) :: R
        real(dp), intent(in) :: K
        real(dp), intent(in) :: beta
        real(dp), intent(in) :: L


        real(dp) :: check_sign
        character(len=60) :: sub_name
        integer :: diagnostics_level

        !------
        sub_name = 'rqk_unit'
        call enter_exit(sub_name,1)
        call get_diagnostics_level(diagnostics_level)

        dQ = dt*(Pgrad - k*Q**2.0_dp - R*Q)/L

        Q = Q+dQ

        !write(*,*) Pgrad, R, K, Q, beta, Q**beta, 1e-7**beta

        call enter_exit(sub_name,2)

    end subroutine rqkl_unit




end module fetal