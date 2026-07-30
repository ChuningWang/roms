!
!  isf_cg_solve.h - SSA Picard iteration with matrix-free PCG solver.
!
!  Solves the SSA stress balance equation A(u) * u = taud using:
!    - Picard outer loop (freeze viscosity, solve linear, update viscosity)
!    - Matrix-free preconditioned conjugate gradient (PCG) inner solve
!    - Jacobi (diagonal) preconditioner
!
!  Reference: STREAMICE streamice_cg_solve_matfree.F, streamice_cg_functions.F
!
!  Copy old velocity as initial guess for new time level.
!
      DO j=JstrP,JendP
        DO i=IstrP,IendP
          Si(i,j,ksa_new,isIuis)=Si(i,j,ksa_old,isIuis)
          Si(i,j,ksa_new,isIvis)=Si(i,j,ksa_old,isIvis)
          IF (pmask_ssa(i,j).eq.0.0_r8) THEN
            Si(i,j,ksa_new,isIuis)=0.0_r8
            Si(i,j,ksa_new,isIvis)=0.0_r8
          END IF
        END DO
      END DO
!
!  Enforce lateral boundary conditions (free-slip / ice divide):
!    - v = 0 at southern/northern domain boundaries (no cross-flow)
!    - u = 0 at western domain boundary (ice divide / symmetry)
!    - Eastern boundary (calving front) is left free.
!
      IF (Jstr.eq.1) THEN
        DO i=IstrP,IendP
          Si(i,JstrP,ksa_new,isIvis)=0.0_r8
        END DO
      END IF
      IF (Jend.eq.Mm(ng)) THEN
        DO i=IstrP,IendP
          Si(i,JendP,ksa_new,isIvis)=0.0_r8
        END DO
      END IF
      IF (Istr.eq.1) THEN
        DO j=JstrP,JendP
          Si(IstrP,j,ksa_new,isIuis)=0.0_r8
        END DO
      END IF
!
!  Warm-start: if initial velocity is near zero, use a simple
!  basal drag balance (tau_d = C^2 * u) to provide a physically
!  reasonable initial guess.  This avoids the Picard "high-viscosity
!  trap" where zero velocity -> eps_e = eps_glen_min -> huge nuH
!  -> tiny velocity -> convergence to the wrong fixed point.
!
      C_basal_sq=C_basal(ng)*C_basal(ng)
      DO j=JstrP,JendP
        DO i=IstrP,IendP
          IF (pmask_ssa(i,j).eq.1.0_r8) THEN
            u_slide=SQRT(Si(i,j,ksa_new,isIuis)**2+                    &
     &                   Si(i,j,ksa_new,isIvis)**2)
            IF (u_slide.lt.1.0E-5_r8) THEN
              avg_ff=0.25_r8*(float_frac(i,j)+                        &
     &                        float_frac(i-1,j)+                       &
     &                        float_frac(i,j-1)+                       &
     &                        float_frac(i-1,j-1))
              jac_avg=0.25_r8*(                                        &
     &          1.0_r8/(pm(i-1,j-1)*pn(i-1,j-1))+                     &
     &          1.0_r8/(pm(i,  j-1)*pn(i,  j-1))+                     &
     &          1.0_r8/(pm(i-1,j  )*pn(i-1,j  ))+                     &
     &          1.0_r8/(pm(i,  j  )*pn(i,  j  )))
              cff=(1.0_r8-avg_ff)/(C_basal_sq*jac_avg*0.25_r8)
              Si(i,j,ksa_new,isIuis)=taudx_psi(i,j)*cff
              Si(i,j,ksa_new,isIvis)=taudy_psi(i,j)*cff
            END IF
          END IF
        END DO
      END DO
!
!  Re-enforce lateral BCs after warm-start.
!
      IF (Jstr.eq.1) THEN
        DO i=IstrP,IendP
          Si(i,JstrP,ksa_new,isIvis)=0.0_r8
        END DO
      END IF
      IF (Jend.eq.Mm(ng)) THEN
        DO i=IstrP,IendP
          Si(i,JendP,ksa_new,isIvis)=0.0_r8
        END DO
      END IF
      IF (Istr.eq.1) THEN
        DO j=JstrP,JendP
          Si(IstrP,j,ksa_new,isIuis)=0.0_r8
        END DO
      END IF
!
!  Initialize CG work arrays to zero over full scratch range.
!
      DO j=JminS,JmaxS
        DO i=IminS,ImaxS
          Ru(i,j)=0.0_r8
          Rv(i,j)=0.0_r8
          Zu(i,j)=0.0_r8
          Zv(i,j)=0.0_r8
          Du(i,j)=0.0_r8
          Dv(i,j)=0.0_r8
          Au(i,j)=0.0_r8
          Av(i,j)=0.0_r8
          diag_u(i,j)=0.0_r8
          diag_v(i,j)=0.0_r8
          Ru_old(i,j)=0.0_r8
          Rv_old(i,j)=0.0_r8
          Zu_old(i,j)=0.0_r8
          Zv_old(i,j)=0.0_r8
          u_scr(i,j)=0.0_r8
          v_scr(i,j)=0.0_r8
        END DO
      END DO
!
!  NOTE: Q1 shape functions (DPhi, phi_val, Xquad, jacq) are precomputed
!  in isf_mom.F before this include file.
!
!=======================================================================
!  Picard iteration loop.
!=======================================================================
!
      cgtol=cg_tol(ng)
      converged=.FALSE.
!
      DO nl_iter=1,nl_max_iter(ng)
!
!  Compute effective viscosity and basal drag from current velocity
!  estimate (ksa_new).  Moved inside the Picard loop so that the
!  viscosity is always consistent with the current velocity, avoiding
!  the use of uninitialized or stale values from a pre-loop call.
!
#   include <isf_viscosity.h>
!
!  Compute Jacobi diagonal (diag_u, diag_v) from nuH and beta.
!  Element loop: for each RHO cell, accumulate diagonal contribution
!  to its 4 PSI corner nodes.
!
        DO j=JstrP,JendP
          DO i=IstrP,IendP
            diag_u(i,j)=0.0_r8
            diag_v(i,j)=0.0_r8
          END DO
        END DO
!
        DO j=Jstr,Jend
          DO i=Istr,Iend
            IF (Si(i,j,ksa_old,isIhis).gt.0.0_r8) THEN
              dx_rho=1.0_r8/pm(i,j)
              dy_rho=1.0_r8/pn(i,j)
              jac=1.0_r8/(pm(i,j)*pn(i,j))
!
              DO nq=1,4
!
!  Diagonal: set u_node=1 for node k, compute its contribution.
!
                DO k_node=1,4
                  IF (k_node.eq.1) THEN
                    i1=i;   j1=j
                  ELSE IF (k_node.eq.2) THEN
                    i1=i+1; j1=j
                  ELSE IF (k_node.eq.3) THEN
                    i1=i;   j1=j+1
                  ELSE
                    i1=i+1; j1=j+1
                  END IF
!
                  IF (pmask_ssa(i1,j1).eq.1.0_r8) THEN
                    dphi_x=DPhi(k_node,nq,1)*2.0_r8/dx_rho
                    dphi_y=DPhi(k_node,nq,2)*2.0_r8/dy_rho
                    phiv=phi_val(k_node,nq)
!
!  Diagonal: set u_node=1, all others=0, all v=0.
!  Then exx=dphi_x, eyy=0, exy=0.5*dphi_y.
!  u-diag = DPhi_x*(4*exx+2*eyy) + DPhi_y*(2*exy)
!         = DPhi_x*(4*dphi_x) + DPhi_y*(dphi_y)
!  v-diag: set v_node=1, all others=0, all u=0.
!  Then exx=0, eyy=dphi_y, exy=0.5*dphi_x.
!  v-diag = DPhi_y*(4*eyy+2*exx) + DPhi_x*(2*exy)
!         = DPhi_y*(4*dphi_y) + DPhi_x*(dphi_x)
!
                    diag_u(i1,j1)=diag_u(i1,j1)+                            &
     &                jacq(nq)*jac*nuH_rho(i,j)*                               &
     &                (4.0_r8*dphi_x*dphi_x+dphi_y*dphi_y)+                &
     &                jacq(nq)*jac*beta_rho(i,j)*phiv*phiv
!
                    diag_v(i1,j1)=diag_v(i1,j1)+                            &
     &                jacq(nq)*jac*nuH_rho(i,j)*                               &
     &                (4.0_r8*dphi_y*dphi_y+dphi_x*dphi_x)+                &
     &                jacq(nq)*jac*beta_rho(i,j)*phiv*phiv
                  END IF
                END DO  ! k_node
              END DO  ! nq
            END IF
          END DO  ! i
        END DO  ! j
!
!  Exchange halo for diagonal.
!
        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
          CALL exchange_p2d_tile(ng,tile,LBi,UBi,LBj,UBj,diag_u)
          CALL exchange_p2d_tile(ng,tile,LBi,UBi,LBj,UBj,diag_v)
        END IF
# ifdef DISTRIBUTE
        CALL mp_exchange2d(ng,tile,iNLM,2,LBi,UBi,LBj,UBj,               &
     &                     NghostPoints,EWperiodic(ng),NSperiodic(ng),     &
     &                     diag_u,diag_v)
# endif
!
!  Compute initial residual: r = b - A*u
!  First compute A*u, then subtract from taud.
!
        DO j=JstrP-1,JendP+1
          DO i=IstrP-1,IendP+1
            Au(i,j)=0.0_r8
            Av(i,j)=0.0_r8
          END DO
        END DO
!
        DO j=Jstr,Jend
          DO i=Istr,Iend
            IF (Si(i,j,ksa_old,isIhis).gt.0.0_r8) THEN
              dx_rho=1.0_r8/pm(i,j)
              dy_rho=1.0_r8/pn(i,j)
!
              DO nq=1,4
!
!  Interpolate velocity and compute strain rates at quad point.
!
                uq=0.0_r8
                vq=0.0_r8
                DO k_node=1,4
                  IF (k_node.eq.1) THEN
                    i1=i;   j1=j
                  ELSE IF (k_node.eq.2) THEN
                    i1=i+1; j1=j
                  ELSE IF (k_node.eq.3) THEN
                    i1=i;   j1=j+1
                  ELSE
                    i1=i+1; j1=j+1
                  END IF
                  uq=uq+Si(i1,j1,ksa_new,isIuis)*phi_val(k_node,nq)
                  vq=vq+Si(i1,j1,ksa_new,isIvis)*phi_val(k_node,nq)
                END DO
!
                exx=0.0_r8
                eyy=0.0_r8
                exy=0.0_r8
                DO k_node=1,4
                  IF (k_node.eq.1) THEN
                    i1=i;   j1=j
                  ELSE IF (k_node.eq.2) THEN
                    i1=i+1; j1=j
                  ELSE IF (k_node.eq.3) THEN
                    i1=i;   j1=j+1
                  ELSE
                    i1=i+1; j1=j+1
                  END IF
                  dphi_x=DPhi(k_node,nq,1)*2.0_r8/dx_rho
                  dphi_y=DPhi(k_node,nq,2)*2.0_r8/dy_rho
                  exx=exx+Si(i1,j1,ksa_new,isIuis)*dphi_x
                  eyy=eyy+Si(i1,j1,ksa_new,isIvis)*dphi_y
                  exy=exy+0.5_r8*(Si(i1,j1,ksa_new,isIuis)*dphi_y+        &
     &                            Si(i1,j1,ksa_new,isIvis)*dphi_x)
                END DO
!
!  Scatter contributions to 4 corner nodes.
!
                DO k_node=1,4
                  IF (k_node.eq.1) THEN
                    i1=i;   j1=j
                  ELSE IF (k_node.eq.2) THEN
                    i1=i+1; j1=j
                  ELSE IF (k_node.eq.3) THEN
                    i1=i;   j1=j+1
                  ELSE
                    i1=i+1; j1=j+1
                  END IF
!
                  IF (pmask_ssa(i1,j1).eq.1.0_r8) THEN
                    dphi_x=DPhi(k_node,nq,1)*2.0_r8/dx_rho
                    dphi_y=DPhi(k_node,nq,2)*2.0_r8/dy_rho
                    phiv=phi_val(k_node,nq)
!
                    Au(i1,j1)=Au(i1,j1)+                                    &
     &                jacq(nq)*jac*nuH_rho(i,j)*                                 &
     &                (dphi_x*(4.0_r8*exx+2.0_r8*eyy)+                       &
     &                 dphi_y*(2.0_r8*exy))+                                  &
     &                jacq(nq)*jac*beta_rho(i,j)*phiv*uq
!
                    Av(i1,j1)=Av(i1,j1)+                                    &
     &                jacq(nq)*jac*nuH_rho(i,j)*                                 &
     &                (dphi_y*(4.0_r8*eyy+2.0_r8*exx)+                       &
     &                 dphi_x*(2.0_r8*exy))+                                  &
     &                jacq(nq)*jac*beta_rho(i,j)*phiv*vq
                  END IF
                END DO  ! k_node
              END DO  ! nq
            END IF
          END DO  ! i
        END DO  ! j
!
!  Exchange halo for A*u.
!
        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
          CALL exchange_p2d_tile(ng,tile,LBi,UBi,LBj,UBj,Au)
          CALL exchange_p2d_tile(ng,tile,LBi,UBi,LBj,UBj,Av)
        END IF
# ifdef DISTRIBUTE
        CALL mp_exchange2d(ng,tile,iNLM,2,LBi,UBi,LBj,UBj,               &
     &                     NghostPoints,EWperiodic(ng),NSperiodic(ng),     &
     &                     Au,Av)
# endif
!
!  Initial residual: r = b - A*u
!  b = taudx_psi/taudy_psi is the FE-assembled driving stress nodal
!  force vector (computed in isf_driving_stress.h).
!
        dot_zr=0.0_r8
        DO j=JstrP,JendP
          DO i=IstrP,IendP
            IF (pmask_ssa(i,j).eq.1.0_r8) THEN
              Ru(i,j)=taudx_psi(i,j)-Au(i,j)
              Rv(i,j)=taudy_psi(i,j)-Av(i,j)
!  Enforce BCs: zero residual for constrained velocity components.
              IF ((Jstr.eq.1.and.j.eq.JstrP).or.                        &
     &            (Jend.eq.Mm(ng).and.j.eq.JendP)) THEN
                Rv(i,j)=0.0_r8
              END IF
              IF (Istr.eq.1.and.i.eq.IstrP) THEN
                Ru(i,j)=0.0_r8
              END IF
              IF (diag_u(i,j).gt.eps) THEN
                Zu(i,j)=Ru(i,j)/diag_u(i,j)
              ELSE
                Zu(i,j)=0.0_r8
              END IF
              IF (diag_v(i,j).gt.eps) THEN
                Zv(i,j)=Rv(i,j)/diag_v(i,j)
              ELSE
                Zv(i,j)=0.0_r8
              END IF
              Du(i,j)=Zu(i,j)
              Dv(i,j)=Zv(i,j)
              dot_zr=dot_zr+Zu(i,j)*Ru(i,j)+Zv(i,j)*Rv(i,j)
            ELSE
              Ru(i,j)=0.0_r8
              Rv(i,j)=0.0_r8
              Zu(i,j)=0.0_r8
              Zv(i,j)=0.0_r8
              Du(i,j)=0.0_r8
              Dv(i,j)=0.0_r8
            END IF
          END DO
        END DO
!
!  Global reduction for initial residual norm.
!
# ifdef DISTRIBUTE
        isf_tile_count=0
        isf_dot_prod=dot_zr
        op_handle='SUM'
        CALL mp_reduce(ng,iNLM,1,isf_dot_prod,op_handle)
        dot_zr=isf_dot_prod
# endif
        resid_0=SQRT(ABS(dot_zr))
        resid=resid_0
!
!  CG iteration loop.
!
        DO cg_iter=1,cg_max_iter(ng)
          IF (resid.le.cgtol*resid_0) EXIT
!
!  Compute A*d.
!  First exchange Si halo so that A*D uses current solution at
!  tile-boundary nodes (critical for MPI correctness).
!
          DO j=JstrP,JendP
            DO i=IstrP,IendP
              u_scr(i,j)=Si(i,j,ksa_new,isIuis)
              v_scr(i,j)=Si(i,j,ksa_new,isIvis)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_p2d_tile(ng,tile,LBi,UBi,LBj,UBj,u_scr)
            CALL exchange_p2d_tile(ng,tile,LBi,UBi,LBj,UBj,v_scr)
          END IF
# ifdef DISTRIBUTE
          CALL mp_exchange2d(ng,tile,iNLM,2,LBi,UBi,LBj,UBj,             &
     &                       NghostPoints,EWperiodic(ng),NSperiodic(ng),   &
     &                       u_scr,v_scr)
# endif
          DO j=JstrP-1,JendP+1
            DO i=IstrP-1,IendP+1
              IF (i.lt.IstrP.or.i.gt.IendP.or.                        &
     &            j.lt.JstrP.or.j.gt.JendP) THEN
                Si(i,j,ksa_new,isIuis)=u_scr(i,j)
                Si(i,j,ksa_new,isIvis)=v_scr(i,j)
              END IF
            END DO
          END DO
!
          DO j=JstrP-1,JendP+1
            DO i=IstrP-1,IendP+1
              Au(i,j)=0.0_r8
              Av(i,j)=0.0_r8
            END DO
          END DO
!
          DO j=Jstr,Jend
            DO i=Istr,Iend
              IF (Si(i,j,ksa_old,isIhis).gt.0.0_r8) THEN
                dx_rho=1.0_r8/pm(i,j)
                dy_rho=1.0_r8/pn(i,j)
                DO nq=1,4
                  uq=0.0_r8
                  vq=0.0_r8
                  exx=0.0_r8
                  eyy=0.0_r8
                  exy=0.0_r8
                  DO k_node=1,4
                    IF (k_node.eq.1) THEN
                      i1=i;   j1=j
                    ELSE IF (k_node.eq.2) THEN
                      i1=i+1; j1=j
                    ELSE IF (k_node.eq.3) THEN
                      i1=i;   j1=j+1
                    ELSE
                      i1=i+1; j1=j+1
                    END IF
                    phiv=phi_val(k_node,nq)
                    dphi_x=DPhi(k_node,nq,1)*2.0_r8/dx_rho
                    dphi_y=DPhi(k_node,nq,2)*2.0_r8/dy_rho
                    uq=uq+Du(i1,j1)*phiv
                    vq=vq+Dv(i1,j1)*phiv
                    exx=exx+Du(i1,j1)*dphi_x
                    eyy=eyy+Dv(i1,j1)*dphi_y
                    exy=exy+0.5_r8*(Du(i1,j1)*dphi_y+Dv(i1,j1)*dphi_x)
                  END DO
!
                  DO k_node=1,4
                    IF (k_node.eq.1) THEN
                      i1=i;   j1=j
                    ELSE IF (k_node.eq.2) THEN
                      i1=i+1; j1=j
                    ELSE IF (k_node.eq.3) THEN
                      i1=i;   j1=j+1
                    ELSE
                      i1=i+1; j1=j+1
                    END IF
                    IF (pmask_ssa(i1,j1).eq.1.0_r8) THEN
                      dphi_x=DPhi(k_node,nq,1)*2.0_r8/dx_rho
                      dphi_y=DPhi(k_node,nq,2)*2.0_r8/dy_rho
                      phiv=phi_val(k_node,nq)
                      Au(i1,j1)=Au(i1,j1)+                                  &
     &                  jacq(nq)*jac*nuH_rho(i,j)*                               &
     &                  (dphi_x*(4.0_r8*exx+2.0_r8*eyy)+                     &
     &                   dphi_y*(2.0_r8*exy))+                                &
     &                  jacq(nq)*jac*beta_rho(i,j)*phiv*uq
                      Av(i1,j1)=Av(i1,j1)+                                  &
     &                  jacq(nq)*jac*nuH_rho(i,j)*                               &
     &                  (dphi_y*(4.0_r8*eyy+2.0_r8*exx)+                     &
     &                   dphi_x*(2.0_r8*exy))+                                &
     &                  jacq(nq)*jac*beta_rho(i,j)*phiv*vq
                    END IF
                  END DO
                END DO
              END IF
            END DO
          END DO
!
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_p2d_tile(ng,tile,LBi,UBi,LBj,UBj,Au)
            CALL exchange_p2d_tile(ng,tile,LBi,UBi,LBj,UBj,Av)
          END IF
# ifdef DISTRIBUTE
          CALL mp_exchange2d(ng,tile,iNLM,2,LBi,UBi,LBj,UBj,             &
     &                       NghostPoints,EWperiodic(ng),NSperiodic(ng),   &
     &                       Au,Av)
# endif
!
!  alpha = (z,r) / (d,A*d)
!
          dot_dad=0.0_r8
          dot_zr_old=0.0_r8
          DO j=JstrP,JendP
            DO i=IstrP,IendP
              IF (pmask_ssa(i,j).eq.1.0_r8) THEN
                dot_dad=dot_dad+Du(i,j)*Au(i,j)+Dv(i,j)*Av(i,j)
                dot_zr_old=dot_zr_old+Zu(i,j)*Ru(i,j)+Zv(i,j)*Rv(i,j)
              END IF
            END DO
          END DO
!
# ifdef DISTRIBUTE
          isf_dot_prod=dot_dad
          op_handle='SUM'
          CALL mp_reduce(ng,iNLM,1,isf_dot_prod,op_handle)
          dot_dad=isf_dot_prod
          isf_dot_prod=dot_zr_old
          CALL mp_reduce(ng,iNLM,1,isf_dot_prod,op_handle)
          dot_zr_old=isf_dot_prod
# endif
!
          IF (ABS(dot_dad).lt.eps) EXIT
          alpha=dot_zr_old/dot_dad
!
!  Update solution, residual, and preconditioned residual.
!
          dot_zr=0.0_r8
          DO j=JstrP,JendP
            DO i=IstrP,IendP
              IF (pmask_ssa(i,j).eq.1.0_r8) THEN
                Si(i,j,ksa_new,isIuis)=Si(i,j,ksa_new,isIuis)+alpha*Du(i,j)
                Si(i,j,ksa_new,isIvis)=Si(i,j,ksa_new,isIvis)+alpha*Dv(i,j)
                Ru_old(i,j)=Ru(i,j)
                Rv_old(i,j)=Rv(i,j)
                Zu_old(i,j)=Zu(i,j)
                Zv_old(i,j)=Zv(i,j)
                Ru(i,j)=Ru(i,j)-alpha*Au(i,j)
                Rv(i,j)=Rv(i,j)-alpha*Av(i,j)
!  Enforce Dirichlet BCs on residual (prevent leakage).
                IF ((Jstr.eq.1.and.j.eq.JstrP).or.                    &
     &              (Jend.eq.Mm(ng).and.j.eq.JendP)) THEN
                  Rv(i,j)=0.0_r8
                END IF
                IF (Istr.eq.1.and.i.eq.IstrP) THEN
                  Ru(i,j)=0.0_r8
                END IF
                IF (diag_u(i,j).gt.eps) THEN
                  Zu(i,j)=Ru(i,j)/diag_u(i,j)
                ELSE
                  Zu(i,j)=0.0_r8
                END IF
                IF (diag_v(i,j).gt.eps) THEN
                  Zv(i,j)=Rv(i,j)/diag_v(i,j)
                ELSE
                  Zv(i,j)=0.0_r8
                END IF
                dot_zr=dot_zr+Zu(i,j)*Ru(i,j)+Zv(i,j)*Rv(i,j)
              END IF
            END DO
          END DO
!
# ifdef DISTRIBUTE
          isf_dot_prod=dot_zr
          op_handle='SUM'
          CALL mp_reduce(ng,iNLM,1,isf_dot_prod,op_handle)
          dot_zr=isf_dot_prod
# endif
!
!  beta = (z_new,r_new) / (z_old,r_old)
!
          IF (ABS(dot_zr_old).lt.eps) EXIT
          beta_cg=dot_zr/dot_zr_old
!
!  New search direction.
!
          DO j=JstrP,JendP
            DO i=IstrP,IendP
              IF (pmask_ssa(i,j).eq.1.0_r8) THEN
                Du(i,j)=Zu(i,j)+beta_cg*Du(i,j)
                Dv(i,j)=Zv(i,j)+beta_cg*Dv(i,j)
!  Enforce Dirichlet BCs on search direction.
                IF ((Jstr.eq.1.and.j.eq.JstrP).or.                    &
     &              (Jend.eq.Mm(ng).and.j.eq.JendP)) THEN
                  Dv(i,j)=0.0_r8
                END IF
                IF (Istr.eq.1.and.i.eq.IstrP) THEN
                  Du(i,j)=0.0_r8
                END IF
              END IF
            END DO
          END DO
!
          resid=SQRT(ABS(dot_zr))
!
        END DO  ! cg_iter
!
!  Exchange halo for solution.
!
        DO j=JstrP,JendP
          DO i=IstrP,IendP
            u_scr(i,j)=Si(i,j,ksa_new,isIuis)
            v_scr(i,j)=Si(i,j,ksa_new,isIvis)
          END DO
        END DO
        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
          CALL exchange_p2d_tile(ng,tile,LBi,UBi,LBj,UBj,u_scr)
          CALL exchange_p2d_tile(ng,tile,LBi,UBi,LBj,UBj,v_scr)
        END IF
# ifdef DISTRIBUTE
        CALL mp_exchange2d(ng,tile,iNLM,2,LBi,UBi,LBj,UBj,               &
     &                     NghostPoints,EWperiodic(ng),NSperiodic(ng),     &
     &                     u_scr,v_scr)
# endif
        DO j=JstrP,JendP
          DO i=IstrP,IendP
            Si(i,j,ksa_new,isIuis)=u_scr(i,j)
            Si(i,j,ksa_new,isIvis)=v_scr(i,j)
          END DO
        END DO
!  Re-enforce lateral BCs after halo exchange.
        IF (Jstr.eq.1) THEN
          DO i=IstrP,IendP
            Si(i,JstrP,ksa_new,isIvis)=0.0_r8
          END DO
        END IF
        IF (Jend.eq.Mm(ng)) THEN
          DO i=IstrP,IendP
            Si(i,JendP,ksa_new,isIvis)=0.0_r8
          END DO
        END IF
        IF (Istr.eq.1) THEN
          DO j=JstrP,JendP
            Si(IstrP,j,ksa_new,isIuis)=0.0_r8
          END DO
        END IF
!
!  PISM velocity cap.
!
        DO j=JstrP,JendP
          DO i=IstrP,IendP
            IF (pmask_ssa(i,j).eq.1.0_r8) THEN
              u_mid=Si(i,j,ksa_new,isIuis)
              v_mid=Si(i,j,ksa_new,isIvis)
              u_slide=SQRT(u_mid*u_mid+v_mid*v_mid)
              IF (u_slide.gt.ssa_max_speed(ng)) THEN
                cff=ssa_max_speed(ng)/u_slide
                Si(i,j,ksa_new,isIuis)=u_mid*cff
                Si(i,j,ksa_new,isIvis)=v_mid*cff
              END IF
            END IF
          END DO
        END DO
!
!  Recompute viscosity with new velocity for next Picard iteration.
!
#   include <isf_viscosity.h>
!
!  Picard relaxation: u = alpha*u_new + (1-alpha)*u_old.
!  Prevents oscillation between high/low viscosity regimes.
!  Applied BEFORE convergence check so the check reflects the
!  actual iterate that will be used in the next iteration.
!
        IF (nl_iter.gt.1) THEN
          DO j=JstrP,JendP
            DO i=IstrP,IendP
              IF (pmask_ssa(i,j).eq.1.0_r8) THEN
                Si(i,j,ksa_new,isIuis)=picard_alpha(ng)*               &
     &            Si(i,j,ksa_new,isIuis)+(1.0_r8-picard_alpha(ng))*    &
     &            Si(i,j,ksa_old,isIuis)
                Si(i,j,ksa_new,isIvis)=picard_alpha(ng)*               &
     &            Si(i,j,ksa_new,isIvis)+(1.0_r8-picard_alpha(ng))*    &
     &            Si(i,j,ksa_old,isIvis)
              END IF
            END DO
          END DO
        END IF
!
!  Picard convergence check: relative velocity change.
!  ||u_relaxed - u_old|| / ||u_relaxed|| < nl_tol
!
        vel_change=0.0_r8
        vel_norm=0.0_r8
        DO j=JstrP,JendP
          DO i=IstrP,IendP
            IF (pmask_ssa(i,j).eq.1.0_r8) THEN
              vel_change=vel_change+                                    &
     &          (Si(i,j,ksa_new,isIuis)-Si(i,j,ksa_old,isIuis))**2+    &
     &          (Si(i,j,ksa_new,isIvis)-Si(i,j,ksa_old,isIvis))**2
              vel_norm=vel_norm+                                        &
     &          Si(i,j,ksa_new,isIuis)**2+                             &
     &          Si(i,j,ksa_new,isIvis)**2
            END IF
          END DO
        END DO
# ifdef DISTRIBUTE
        isf_dot_prod=vel_change
        op_handle='SUM'
        CALL mp_reduce(ng,iNLM,1,isf_dot_prod,op_handle)
        vel_change=isf_dot_prod
        isf_dot_prod=vel_norm
        op_handle='SUM'
        CALL mp_reduce(ng,iNLM,1,isf_dot_prod,op_handle)
        vel_norm=isf_dot_prod
# endif
        IF (vel_norm.gt.eps) THEN
          nl_resid=SQRT(vel_change/vel_norm)
        ELSE
          nl_resid=0.0_r8
        END IF
!
        IF (nl_resid.le.nl_tol(ng)) THEN
          converged=.TRUE.
        END IF
!
!  DEBUG: Picard iteration diagnostics (tile 0 only, sparse output).
!
        IF (tile.eq.0) THEN
          IF (converged.or.nl_iter.eq.1.or.                            &
     &        MOD(nl_iter,10).eq.0) THEN
            u_slide=0.0_r8
            DO j=JstrP,JendP
              DO i=IstrP,IendP
                IF (pmask_ssa(i,j).eq.1.0_r8) THEN
                  cff=SQRT(Si(i,j,ksa_new,isIuis)**2+                  &
     &                     Si(i,j,ksa_new,isIvis)**2)
                  IF (cff.gt.u_slide) u_slide=cff
                END IF
              END DO
            END DO
  #ifdef ISF_DEBUG
          WRITE(*,'(A,I3,A,ES12.4,A,ES12.4,A,L1)')                   &
     &        '  Picard iter ',nl_iter,' nl_resid=',nl_resid,          &
     &        ' max|u|=',u_slide,' conv=',converged
#endif
          END IF
        END IF
!
        IF (converged) EXIT
!
!  Update ksa_old for next Picard iteration (use new velocity as guess).
!
        DO j=JstrP,JendP
          DO i=IstrP,IendP
            Si(i,j,ksa_old,isIuis)=Si(i,j,ksa_new,isIuis)
            Si(i,j,ksa_old,isIvis)=Si(i,j,ksa_new,isIvis)
          END DO
        END DO
!
      END DO  ! nl_iter (Picard)
!
!  Post-Picard cleanup: zero velocity at ice-margin PSI points
!  (min surrounding H < 10 m).  Applied once after convergence
!  to remove garbage without affecting the nonlinear iteration.
!
      DO j=JstrP,JendP
        DO i=IstrP,IendP
          IF (pmask_ssa(i,j).eq.1.0_r8) THEN
            IF (MIN(Si(i,j,ksa_old,isIhis),                           &
     &            Si(i-1,j,ksa_old,isIhis),                           &
     &            Si(i,j-1,ksa_old,isIhis),                           &
     &            Si(i-1,j-1,ksa_old,isIhis)).lt.10.0_r8) THEN
              Si(i,j,ksa_new,isIuis)=0.0_r8
              Si(i,j,ksa_new,isIvis)=0.0_r8
            END IF
          END IF
        END DO
      END DO
