!
!  isf_driving_stress.h - Compute SSA driving stress at PSI points.
!
!  Computes taudx_psi = -rho_i * g * H_psi * dS/dx
!          taudy_psi = -rho_i * g * H_psi * dS/dy
!
!  For each active PSI node (i,j), the driving stress is computed from
!  the surface elevation gradient of the 4 surrounding RHO cells.
!  Ice-front boundary pressure is added where one side has no ice.
!
!  Reference: STREAMICE streamice_driving_stress.F
!
!  Initialize driving stress to zero.
!
      DO j=JstrP,JendP
        DO i=IstrP,IendP
          taudx_psi(i,j)=0.0_r8
          taudy_psi(i,j)=0.0_r8
        END DO
      END DO
!
      rhoi_g=rho_i*gravity
      rhow_g=rho_w_local*gravity
!
!  Compute driving stress at active PSI nodes.
!  PSI node (i,j) is surrounded by RHO cells (i-1,j-1), (i,j-1),
!  (i-1,j), (i,j).  Surface gradient is computed by averaging
!  the differences across the two pairs of cells in each direction.
!
      DO j=JstrP,JendP
        DO i=IstrP,IendP
          IF (pmask_ssa(i,j).eq.1.0_r8) THEN
!
!  Check ice presence in 4 surrounding RHO cells.
!  Cell layout:  (i-1,j)  | (i,j)
!                (i-1,j-1) | (i,j-1)
!
            has_ice_w=0
            has_ice_e=0
            has_ice_s=0
            has_ice_n=0
            IF (Si(i-1,j,ksa_old,isIhis).gt.0.0_r8) has_ice_w=1
            IF (Si(i,j,ksa_old,isIhis).gt.0.0_r8)   has_ice_e=1
            IF (Si(i,j-1,ksa_old,isIhis).gt.0.0_r8) has_ice_s=1
            IF (Si(i-1,j-1,ksa_old,isIhis).gt.0.0_r8) has_ice_s=1
            IF (Si(i-1,j,ksa_old,isIhis).gt.0.0_r8) has_ice_n=1
            IF (Si(i,j,ksa_old,isIhis).gt.0.0_r8)   has_ice_n=1
!
!  X-direction driving stress.
!  Use cells (i-1,j) and (i,j) for the upper row,
!  and cells (i-1,j-1) and (i,j-1) for the lower row.
!
            IF ((has_ice_w.eq.1).and.(has_ice_e.eq.1)) THEN
!
!  Both cells have ice: standard surface gradient.
!
              S_grad_x=0.5_r8*                                            &
     &          ((Fi(i,j,icIzsf)-Fi(i-1,j,icIzsf))+                      &
     &           (Fi(i,j-1,icIzsf)-Fi(i-1,j-1,icIzsf)))*pm(i,j)
              H_sum=0.25_r8*                                              &
     &          (Si(i,j,ksa_old,isIhis)+Si(i-1,j,ksa_old,isIhis)+         &
     &           Si(i,j-1,ksa_old,isIhis)+Si(i-1,j-1,ksa_old,isIhis))
              taudx_psi(i,j)=-rhoi_g*H_sum*S_grad_x
!
            ELSE IF (has_ice_w.eq.1) THEN
!
!  Ice only on west side: ice-front boundary pressure.
!
              H_rho=Si(i-1,j,ksa_old,isIhis)
              D_rho=-h(i-1,j)
              IF (rmask_ist_scr(i-1,j).eq.1.0_r8) THEN
                front_press=rhoi_g*                                       &
     &            (1.0_r8-rho_i/rho_w_local)*H_rho*H_rho
              ELSE
                front_press=rhoi_g*H_rho*H_rho -                         &
     &            rhow_g*MAX(D_rho,0.0_r8)*MAX(D_rho,0.0_r8)
              END IF
              taudx_psi(i,j)=0.25_r8*front_press
!
            ELSE IF (has_ice_e.eq.1) THEN
!
!  Ice only on east side: ice-front boundary pressure (reversed sign).
!
              H_rho=Si(i,j,ksa_old,isIhis)
              D_rho=-h(i,j)
              IF (rmask_ist_scr(i,j).eq.1.0_r8) THEN
                front_press=rhoi_g*                                       &
     &            (1.0_r8-rho_i/rho_w_local)*H_rho*H_rho
              ELSE
                front_press=rhoi_g*H_rho*H_rho -                         &
     &            rhow_g*MAX(D_rho,0.0_r8)*MAX(D_rho,0.0_r8)
              END IF
              taudx_psi(i,j)=-0.25_r8*front_press
            END IF
!
!  Y-direction driving stress.
!  Use cells (i,j-1) and (i,j) for the right column,
!  and cells (i-1,j-1) and (i-1,j) for the left column.
!
            IF ((has_ice_s.eq.1).and.(has_ice_n.eq.1)) THEN
!
!  Both cells have ice: standard surface gradient.
!
              S_grad_y=0.5_r8*                                            &
     &          ((Fi(i,j,icIzsf)-Fi(i,j-1,icIzsf))+                      &
     &           (Fi(i-1,j,icIzsf)-Fi(i-1,j-1,icIzsf)))*pn(i,j)
              H_sum=0.25_r8*                                              &
     &          (Si(i,j,ksa_old,isIhis)+Si(i-1,j,ksa_old,isIhis)+         &
     &           Si(i,j-1,ksa_old,isIhis)+Si(i-1,j-1,ksa_old,isIhis))
              taudy_psi(i,j)=-rhoi_g*H_sum*S_grad_y
!
            ELSE IF (has_ice_s.eq.1) THEN
!
!  Ice only on south side: ice-front boundary pressure.
!
              H_rho=Si(i,j-1,ksa_old,isIhis)
              D_rho=-h(i,j-1)
              IF (rmask_ist_scr(i,j-1).eq.1.0_r8) THEN
                front_press=rhoi_g*                                       &
     &            (1.0_r8-rho_i/rho_w_local)*H_rho*H_rho
              ELSE
                front_press=rhoi_g*H_rho*H_rho -                         &
     &            rhow_g*MAX(D_rho,0.0_r8)*MAX(D_rho,0.0_r8)
              END IF
              taudy_psi(i,j)=0.25_r8*front_press
!
            ELSE IF (has_ice_n.eq.1) THEN
!
!  Ice only on north side: ice-front boundary pressure (reversed sign).
!
              H_rho=Si(i,j,ksa_old,isIhis)
              D_rho=-h(i,j)
              IF (rmask_ist_scr(i,j).eq.1.0_r8) THEN
                front_press=rhoi_g*                                       &
     &            (1.0_r8-rho_i/rho_w_local)*H_rho*H_rho
              ELSE
                front_press=rhoi_g*H_rho*H_rho -                         &
     &            rhow_g*MAX(D_rho,0.0_r8)*MAX(D_rho,0.0_r8)
              END IF
              taudy_psi(i,j)=-0.25_r8*front_press
            END IF
!
          END IF  ! pmask_ssa
        END DO
      END DO
!
!  Store driving stress to diagnostic output fields.
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          Fi(i,j,icItdx)=taudx_psi(i,j)
          Fi(i,j,icItdy)=taudy_psi(i,j)
        END DO
      END DO
