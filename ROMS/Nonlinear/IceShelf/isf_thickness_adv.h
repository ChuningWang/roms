!  isf_thickness_adv.h - Ice thickness advection (mass conservation).
!
!  Solves the ice thickness evolution equation in flux form:
!
!    dH/dt = -div(u*H) - melt_rate
!
!  where u is the depth-averaged ice velocity (from SSA solver),
!  H is ice thickness, and melt_rate is the basal melt rate from
!  the 3-equation thermodynamics (positive = melting).
!
!  Numerical scheme:
!    - PSI velocity interpolated to RHO cell faces (2-point average)
!    - First-order upwind flux for monotonicity
!    - CFL-based sub-stepping for stability
!    - Positivity preservation (H >= 0)
!
!  Grid layout (B-grid / PSI velocity):
!    PSI(i,j) is the SW corner of RHO cell (i,j).
!    x-face between RHO(i-1,j) and RHO(i,j) has velocity:
!      u_face(i,j) = 0.5*(u_psi(i,j) + u_psi(i,j+1))
!    y-face between RHO(i,j-1) and RHO(i,j) has velocity:
!      v_face(i,j) = 0.5*(v_psi(i,j) + v_psi(i+1,j))
!
!  Reference: STREAMICE streamice_advect_thickness.F,
!             PISM mass transport (upwind FD)
!
!  Initialize work arrays.
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          H_new(i,j)=0.0_r8
        END DO
      END DO
!
!  Compute face velocities from PSI-point velocity.
!  u_xface(i,j): x-direction velocity on the x-face between
!                RHO(i-1,j) and RHO(i,j).
!  v_yface(i,j): y-direction velocity on the y-face between
!                RHO(i,j-1) and RHO(i,j).
!
      DO j=Jstr-1,Jend+1
        DO i=Istr-1,Iend+1
          u_xface(i,j)=0.5_r8*                                          &
     &        (Si(i,j,ksa_new,isIuis)+Si(i,j+1,ksa_new,isIuis))
          v_yface(i,j)=0.5_r8*                                          &
     &        (Si(i,j,ksa_new,isIvis)+Si(i+1,j,ksa_new,isIvis))
        END DO
      END DO
!
!  Compute maximum CFL number for sub-stepping.
!  icfl = max(|u_face|/dx + |v_face|/dy) over all active cells.
!
      max_icfl=1.0E-20_r8
      DO j=Jstr,Jend
        DO i=Istr,Iend
          IF (Si(i,j,ksa_old,isIhis).gt.0.0_r8) THEN
            cfl_u=MAX(ABS(u_xface(i,j)),                              &
     &                ABS(u_xface(i+1,j)))*pm(i,j)
            cfl_v=MAX(ABS(v_yface(i,j)),                              &
     &                ABS(v_yface(i,j+1)))*pn(i,j)
            loc_icfl=cfl_u+cfl_v
            IF (loc_icfl.gt.max_icfl) max_icfl=loc_icfl
          END IF
        END DO
      END DO
!
!  Global maximum for distributed memory.
!
# ifdef DISTRIBUTE
      isf_dot_prod=max_icfl
      op_handle='MAX'
      CALL mp_reduce(ng,iNLM,1,isf_dot_prod,op_handle)
      max_icfl=isf_dot_prod
# endif
!
!  CFL sub-stepping loop.
!
      dt_full=dt(ng)*REAL(N_ssa(ng)*N_thick(ng),r8)
      dt_rem=dt_full
!
      DO WHILE (dt_rem.gt.1.0E-10_r8)
!
        dt_sub=MIN(cfl_factor(ng)/max_icfl, dt_rem)
!
!  Compute upwind fluxes and update thickness.
!  Flux through x-face: Fx = u_face * H_upwind
!  Flux through y-face: Fy = v_face * H_upwind
!
!  H_new = H_old - dt_sub/area * (Fx_out - Fx_in + Fy_out - Fy_in)
!          - dt_sub * melt_rate
!
        DO j=Jstr,Jend
          DO i=Istr,Iend
            H_old=Si(i,j,ksa_old,isIhis)
!
            IF (H_old.gt.0.0_r8) THEN
!
!  X-direction fluxes.
!  Left x-face (i): between RHO(i-1,j) and RHO(i,j).
!  Right x-face (i+1): between RHO(i,j) and RHO(i+1,j).
!
              H_left=Si(i-1,j,ksa_old,isIhis)
              H_right=Si(i+1,j,ksa_old,isIhis)
!
!  Upwind: if u > 0, H comes from left (i-1); if u < 0, from right (i).
!
              IF (u_xface(i,j).ge.0.0_r8) THEN
                Fx_left=u_xface(i,j)*MAX(H_left,0.0_r8)
              ELSE
                Fx_left=u_xface(i,j)*H_old
              END IF
!
              IF (u_xface(i+1,j).ge.0.0_r8) THEN
                Fx_right=u_xface(i+1,j)*H_old
              ELSE
                Fx_right=u_xface(i+1,j)*MAX(H_right,0.0_r8)
              END IF
!
!  Y-direction fluxes.
!  Bottom y-face (j): between RHO(i,j-1) and RHO(i,j).
!  Top y-face (j+1): between RHO(i,j) and RHO(i,j+1).
!
              H_below=Si(i,j-1,ksa_old,isIhis)
              H_above=Si(i,j+1,ksa_old,isIhis)
!
              IF (v_yface(i,j).ge.0.0_r8) THEN
                Fy_bottom=v_yface(i,j)*MAX(H_below,0.0_r8)
              ELSE
                Fy_bottom=v_yface(i,j)*H_old
              END IF
!
              IF (v_yface(i,j+1).ge.0.0_r8) THEN
                Fy_top=v_yface(i,j+1)*H_old
              ELSE
                Fy_top=v_yface(i,j+1)*MAX(H_above,0.0_r8)
              END IF
!
!  Thickness update: dH/dt = -div(uH) - melt.
!  On a uniform grid: div(uH) ≈ (Fx_right - Fx_left)*pm
!                              + (Fy_top - Fy_bottom)*pn
!
              H_new(i,j)=H_old                                         &
     &          -dt_sub*((Fx_right-Fx_left)*pm(i,j)                   &
     &                  +(Fy_top-Fy_bottom)*pn(i,j))                  &
     &          -dt_sub*Fi(i,j,icImel)
!
!  Positivity preservation and safety cap.
!
              IF (H_new(i,j).lt.0.0_r8) H_new(i,j)=0.0_r8
              IF (H_new(i,j).gt.2000.0_r8) H_new(i,j)=2000.0_r8
!
            ELSE
              H_new(i,j)=0.0_r8
            END IF
          END DO
        END DO
!
!  Update H for next sub-step (copy H_new back to ksa_old for
!  multi-substep advection).
!
        DO j=Jstr,Jend
          DO i=Istr,Iend
            Si(i,j,ksa_old,isIhis)=H_new(i,j)
          END DO
        END DO
!
!  Exchange halo for multi-substep advection (only if more
!  sub-steps remain).
!
        IF (dt_rem-dt_sub.gt.1.0E-10_r8) THEN
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_p2d_tile(ng,tile,LBi,UBi,LBj,UBj,          &
     &                             Si(:,:,ksa_old,isIhis))
          END IF
# ifdef DISTRIBUTE
          CALL mp_exchange2d(ng,tile,iNLM,1,LBi,UBi,LBj,UBj,         &
     &                       NghostPoints,EWperiodic(ng),              &
     &                       NSperiodic(ng),                           &
     &                       Si(:,:,ksa_old,isIhis))
# endif
        END IF
!
        dt_rem=dt_rem-dt_sub
!
      END DO  ! sub-stepping
!
!  Store final thickness at ksa_new.
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          Si(i,j,ksa_new,isIhis)=H_new(i,j)
        END DO
      END DO
!
!  Exchange halo for updated thickness.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_p2d_tile(ng,tile,LBi,UBi,LBj,UBj,              &
     &                         Si(:,:,ksa_new,isIhis))
      END IF
# ifdef DISTRIBUTE
      CALL mp_exchange2d(ng,tile,iNLM,1,LBi,UBi,LBj,UBj,             &
     &                   NghostPoints,EWperiodic(ng),NSperiodic(ng),    &
     &                   Si(:,:,ksa_new,isIhis))
# endif
