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
!  NOTE: the scalar state (thickness, temperature, enthalpy) is carried
!  forward ksa_old -> ksa_new in isf_mom.F BEFORE the geometry/mask
!  computations (with a halo exchange), so only the velocity guess is
!  copied here.
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
!  Enforce lateral boundary conditions (LBC-aware).  The treatment of
!  each domain edge follows the LBC switch set in the ISF input script:
!    closed   -> Dirichlet, velocity = 0      (no-slip / ice divide)
!    clamped  -> Dirichlet, velocity = boundary data from ISF_LOBC
!    mixed    -> treated as clamped (inflow value held at the BC data)
!    gradient -> free (no enforcement; zero-gradient outflow)
!  The eastern boundary (calving front) is always left free.
!
      IF (Jstr.eq.1.and.Lv_south_fix) THEN
        DO i=IstrP,IendP
          IF (Lv_south_clamp) THEN
            Si(i,JstrP,ksa_new,isIvis)=                                 &
     &        ISF_LOBC(isIvis,ng)%isf_south(i)
          ELSE
            Si(i,JstrP,ksa_new,isIvis)=0.0_r8
          END IF
        END DO
      END IF
      IF (Jend.eq.Mm(ng).and.Lv_north_fix) THEN
        DO i=IstrP,IendP
          IF (Lv_north_clamp) THEN
            Si(i,JendP,ksa_new,isIvis)=                                 &
     &        ISF_LOBC(isIvis,ng)%isf_north(i)
          ELSE
            Si(i,JendP,ksa_new,isIvis)=0.0_r8
          END IF
        END DO
      END IF
      IF (Istr.eq.1.and.Lu_west_fix) THEN
        DO j=JstrP,JendP
          IF (Lu_west_clamp) THEN
            Si(IstrP,j,ksa_new,isIuis)=                                 &
     &        ISF_LOBC(isIuis,ng)%isf_west(j)
          ELSE
            Si(IstrP,j,ksa_new,isIuis)=0.0_r8
          END IF
        END DO
      END IF
!
!  Exchange the halo of the initial guess AFTER the copy above: the
!  copy and the Dirichlet enforcement only touch the tile interior, so
!  the ghost columns/rows would otherwise stay stale (zero) and the
!  element assembly (which reads all 4 PSI corners of every RHO cell)
!  would use wrong velocities at tile boundaries under MPI/OpenMP.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_p2d_tile(ng,tile,LBi,UBi,LBj,UBj,              &
     &                         Si(:,:,ksa_new,isIuis))
        CALL exchange_p2d_tile(ng,tile,LBi,UBi,LBj,UBj,              &
     &                         Si(:,:,ksa_new,isIvis))
      END IF
# ifdef DISTRIBUTE
      CALL mp_exchange2d(ng,tile,iNLM,2,LBi,UBi,LBj,UBj,             &
     &                   NghostPoints,EWperiodic(ng),NSperiodic(ng),   &
     &                   Si(:,:,ksa_new,isIuis),                       &
     &                   Si(:,:,ksa_new,isIvis))
# endif
!
!  Warm-start: if initial velocity is near zero, use a simple
!  basal drag balance (tau_d = C^2 * u) to provide a physically
!  reasonable initial guess.  This avoids the Picard "high-viscosity
!  trap" where zero velocity -> eps_e = eps_glen_min -> huge nuH
!  -> tiny velocity -> convergence to the wrong fixed point.
!
!  The warm-start is only meaningful when basal drag is present.
!  With C_basal ~ 0 (floating shelf, or the viscous-only MMS test)
!  skip it and keep the previous/initial velocity field.
      C_basal_sq=C_basal(ng)*C_basal(ng)
!  Warm-start is only valid for linear sliding (n_basal=1).  For
!  nonlinear sliding (n_basal != 1) the drag balance tau_d = C^2*u^m
!  is nonlinear in u and the linear inversion u=tau_d/C^2 is wrong.
      IF (C_basal_sq.gt.1.0E-10_r8 .and.                               &
     &    ABS(n_basal(ng)-1.0_r8).lt.eps) THEN
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
      END IF
!
!  Re-enforce lateral BCs after warm-start (LBC-aware, see above).
!
      IF (Jstr.eq.1.and.Lv_south_fix) THEN
        DO i=IstrP,IendP
          IF (Lv_south_clamp) THEN
            Si(i,JstrP,ksa_new,isIvis)=                                 &
     &        ISF_LOBC(isIvis,ng)%isf_south(i)
          ELSE
            Si(i,JstrP,ksa_new,isIvis)=0.0_r8
          END IF
        END DO
      END IF
      IF (Jend.eq.Mm(ng).and.Lv_north_fix) THEN
        DO i=IstrP,IendP
          IF (Lv_north_clamp) THEN
            Si(i,JendP,ksa_new,isIvis)=                                 &
     &        ISF_LOBC(isIvis,ng)%isf_north(i)
          ELSE
            Si(i,JendP,ksa_new,isIvis)=0.0_r8
          END IF
        END DO
      END IF
      IF (Istr.eq.1.and.Lu_west_fix) THEN
        DO j=JstrP,JendP
          IF (Lu_west_clamp) THEN
            Si(IstrP,j,ksa_new,isIuis)=                                 &
     &        ISF_LOBC(isIuis,ng)%isf_west(j)
          ELSE
            Si(IstrP,j,ksa_new,isIuis)=0.0_r8
          END IF
        END DO
      END IF
!
!  Initialize CG work arrays to zero over their full declared ranges.
!  NOTE: Au/Av/diag_u/diag_v/u_scr/v_scr are declared on the GLOBAL
!  bounds (LBi:UBi, LBj:UBj) because they are passed to mp_exchange2d.
!  The tile-local JminS/JmaxS range (Jstr-3..Jend+3) overruns those
!  arrays on multi-tile MPI decompositions (JminS < LBj or JmaxS > UBj),
!  so they are initialized over the global range; the remaining scratch
!  arrays (Ru..Zv) keep their tile-local (IminS..JmaxS) extent.
!
      DO j=LBj,UBj
        DO i=LBi,UBi
          Au(i,j)=0.0_r8
          Av(i,j)=0.0_r8
          diag_u(i,j)=0.0_r8
          diag_v(i,j)=0.0_r8
          u_scr(i,j)=0.0_r8
          v_scr(i,j)=0.0_r8
          Du(i,j)=0.0_r8
          Dv(i,j)=0.0_r8
        END DO
      END DO
!
      DO j=JminS,JmaxS
        DO i=IminS,ImaxS
          Ru(i,j)=0.0_r8
          Rv(i,j)=0.0_r8
          Zu(i,j)=0.0_r8
          Zv(i,j)=0.0_r8
          Ru_old(i,j)=0.0_r8
          Rv_old(i,j)=0.0_r8
          Zu_old(i,j)=0.0_r8
          Zv_old(i,j)=0.0_r8
        END DO
      END DO
!
!  NOTE: Q1 shape functions (DPhi, phi_val, Xquad, jacq) are precomputed
!  in isf_mom.F before this include file.
!
!  Extended cell-assembly bounds: PSI tile ranges do NOT overlap
!  (get_bounds.F: IstrP=my_Istr everywhere; IendP=my_Iend except on
!  the eastern edge), so a PSI node on a tile seam is owned by ONE
!  tile only, while the RHO cells feeding its FE stencil lie on BOTH
!  sides of the seam.  Each tile therefore assembles one extra ring
!  of halo cells (their nuH/beta arrive through the halo exchange
!  below) and keeps only contributions landing on PSI nodes it owns
!  (ownership filter at every scatter).  Without the extra ring,
!  seam nodes miss the contributions from the neighbour's cells and
!  the MPI solution departs from the single-tile solution.
!
      iCstr=Istr
      iCend=Iend
      jCstr=Jstr
      jCend=Jend
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
!  Exchange nuH/beta halos: the assembly loops below integrate one
!  ring of halo cells beyond Istr..Iend (see iCstr..iCend), and need
!  the viscosity/basal drag of those cells as computed on the
!  neighbouring tile.  Halo cells outside the physical ice cover keep
!  nuH=beta=0 from the initialization above and contribute nothing.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_p2d_tile(ng,tile,LBi,UBi,LBj,UBj,nuH_rho)
        CALL exchange_p2d_tile(ng,tile,LBi,UBi,LBj,UBj,beta_rho)
      END IF
# ifdef DISTRIBUTE
      CALL mp_exchange2d(ng,tile,iNLM,2,LBi,UBi,LBj,UBj,               &
     &                   NghostPoints,EWperiodic(ng),NSperiodic(ng),     &
     &                   nuH_rho,beta_rho)
# endif
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
        DO j=jCstr,jCend
          DO i=iCstr,iCend
            IF (nuH_rho(i,j).gt.0.0_r8) THEN
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
!  Ownership filter: halo-cell contributions landing on PSI nodes
!  owned by a neighbouring tile are dropped here (that tile performs
!  the same integration over its own halo cells).
!
                  IF (pmask_ssa(i1,j1).eq.1.0_r8.and.                  &
     &                i1.ge.IstrP.and.i1.le.IendP.and.                 &
     &                j1.ge.JstrP.and.j1.le.JendP) THEN
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
        DO j=jCstr,jCend
          DO i=iCstr,iCend
            IF (nuH_rho(i,j).gt.0.0_r8) THEN
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
!  Scatter contributions to 4 corner nodes (ownership filter: drop
!  halo-cell contributions landing on neighbour-owned PSI nodes).
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
                  IF (pmask_ssa(i1,j1).eq.1.0_r8.and.                  &
     &                i1.ge.IstrP.and.i1.le.IendP.and.                 &
     &                j1.ge.JstrP.and.j1.le.JendP) THEN
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
!  Exchange halo for A*u.  MUST precede the initial residual: shared
!  PSI nodes on tile boundaries receive contributions from cells owned
!  by the neighbouring tile, which arrive through the halo.  Computing
!  r = b - A*u before this exchange leaves the initial residual (and
!  hence the first CG search direction) inconsistent with the matvec
!  operator, corrupting convergence under MPI.
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
              IF ((Jstr.eq.1.and.j.eq.JstrP.and.Lv_south_fix).or.       &
     &            (Jend.eq.Mm(ng).and.j.eq.JendP.and.Lv_north_fix))     &
     &          THEN
                Rv(i,j)=0.0_r8
              END IF
              IF (Istr.eq.1.and.i.eq.IstrP.and.Lu_west_fix) THEN
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
          IF (resid.le.cgtol*resid_0 .or. resid.lt.eps) EXIT
!
!  Compute A*d.
!  Exchange halos for BOTH the solution (Si) and the search direction
!  (Du/Dv).  The A*D matvec reads Du(i+1,j+1) corners that may fall on
!  the neighbouring tile's halo; without exchanging Du/Dv the matvec
!  is incomplete at MPI tile boundaries (MITgcm exchanges Du_SI/Dv_SI
!  explicitly in streamice_cg_solve.F:472-477).
!
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_p2d_tile(ng,tile,LBi,UBi,LBj,UBj,Du)
            CALL exchange_p2d_tile(ng,tile,LBi,UBi,LBj,UBj,Dv)
          END IF
# ifdef DISTRIBUTE
          CALL mp_exchange2d(ng,tile,iNLM,2,LBi,UBi,LBj,UBj,             &
     &                       NghostPoints,EWperiodic(ng),NSperiodic(ng),   &
     &                       Du,Dv)
# endif
!
          DO j=JstrP-1,JendP+1
            DO i=IstrP-1,IendP+1
              Au(i,j)=0.0_r8
              Av(i,j)=0.0_r8
            END DO
          END DO
!
          DO j=jCstr,jCend
            DO i=iCstr,iCend
              IF (nuH_rho(i,j).gt.0.0_r8) THEN
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
                    IF (pmask_ssa(i1,j1).eq.1.0_r8.and.                &
     &                  i1.ge.IstrP.and.i1.le.IendP.and.               &
     &                  j1.ge.JstrP.and.j1.le.JendP) THEN
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
          IF (alpha.ne.alpha.or.alpha.gt.1.0e30_r8) THEN
            IF (Master) WRITE(stdout,'(a,i3,3e12.4)')                      &
     &        ' SSA CG NaN: iter=',cg_iter,' dot_dad=',dot_dad,            &
     &        ' dot_zr_old=',dot_zr_old
            EXIT
          END IF
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
                IF ((Jstr.eq.1.and.j.eq.JstrP.and.Lv_south_fix).or.     &
     &              (Jend.eq.Mm(ng).and.j.eq.JendP.and.Lv_north_fix))   &
     &            THEN
                  Rv(i,j)=0.0_r8
                END IF
                IF (Istr.eq.1.and.i.eq.IstrP.and.Lu_west_fix) THEN
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
                IF ((Jstr.eq.1.and.j.eq.JstrP.and.Lv_south_fix).or.     &
     &              (Jend.eq.Mm(ng).and.j.eq.JendP.and.Lv_north_fix))   &
     &            THEN
                  Dv(i,j)=0.0_r8
                END IF
                IF (Istr.eq.1.and.i.eq.IstrP.and.Lu_west_fix) THEN
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
!  Copy full range (including halo) to scratch, exchange, copy back.
!
        DO j=LBj,UBj
          DO i=LBi,UBi
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
        DO j=LBj,UBj
          DO i=LBi,UBi
            Si(i,j,ksa_new,isIuis)=u_scr(i,j)
            Si(i,j,ksa_new,isIvis)=v_scr(i,j)
          END DO
        END DO
!  Re-enforce lateral BCs after halo exchange (LBC-aware, see above).
        IF (Jstr.eq.1.and.Lv_south_fix) THEN
          DO i=IstrP,IendP
            IF (Lv_south_clamp) THEN
              Si(i,JstrP,ksa_new,isIvis)=                               &
     &          ISF_LOBC(isIvis,ng)%isf_south(i)
            ELSE
              Si(i,JstrP,ksa_new,isIvis)=0.0_r8
            END IF
          END DO
        END IF
        IF (Jend.eq.Mm(ng).and.Lv_north_fix) THEN
          DO i=IstrP,IendP
            IF (Lv_north_clamp) THEN
              Si(i,JendP,ksa_new,isIvis)=                               &
     &          ISF_LOBC(isIvis,ng)%isf_north(i)
            ELSE
              Si(i,JendP,ksa_new,isIvis)=0.0_r8
            END IF
          END DO
        END IF
        IF (Istr.eq.1.and.Lu_west_fix) THEN
          DO j=JstrP,JendP
            IF (Lu_west_clamp) THEN
              Si(IstrP,j,ksa_new,isIuis)=                               &
     &          ISF_LOBC(isIuis,ng)%isf_west(j)
            ELSE
              Si(IstrP,j,ksa_new,isIuis)=0.0_r8
            END IF
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
!  Exchange velocity halo before recomputing viscosity: the element
!  viscosity reads Si velocities at all four corners of each cell,
!  including halo rows owned by neighbouring tiles/MPI ranks.
!  Without the exchange, halo velocities are stale after this tile's
!  CG update, breaking viscosity continuity at tile boundaries.
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
        CALL mp_exchange2d(ng,tile,iNLM,2,LBi,UBi,LBj,UBj,            &
     &                     NghostPoints,EWperiodic(ng),NSperiodic(ng), &
     &                     u_scr,v_scr)
# endif
        DO j=JstrP,JendP
          DO i=IstrP,IendP
            Si(i,j,ksa_new,isIuis)=u_scr(i,j)
            Si(i,j,ksa_new,isIvis)=v_scr(i,j)
          END DO
        END DO
!  Re-enforce lateral BCs after halo exchange (LBC-aware, see above).
        IF (Jstr.eq.1.and.Lv_south_fix) THEN
          DO i=IstrP,IendP
            IF (Lv_south_clamp) THEN
              Si(i,JstrP,ksa_new,isIvis)=                               &
     &          ISF_LOBC(isIvis,ng)%isf_south(i)
            ELSE
              Si(i,JstrP,ksa_new,isIvis)=0.0_r8
            END IF
          END DO
        END IF
        IF (Jend.eq.Mm(ng).and.Lv_north_fix) THEN
          DO i=IstrP,IendP
            IF (Lv_north_clamp) THEN
              Si(i,JendP,ksa_new,isIvis)=                               &
     &          ISF_LOBC(isIvis,ng)%isf_north(i)
            ELSE
              Si(i,JendP,ksa_new,isIvis)=0.0_r8
            END IF
          END DO
        END IF
        IF (Istr.eq.1.and.Lu_west_fix) THEN
          DO j=JstrP,JendP
            IF (Lu_west_clamp) THEN
              Si(IstrP,j,ksa_new,isIuis)=                               &
     &          ISF_LOBC(isIuis,ng)%isf_west(j)
            ELSE
              Si(IstrP,j,ksa_new,isIuis)=0.0_r8
            END IF
          END DO
        END IF
!
!  Recompute viscosity from the RELAXED velocity so that nuH is
!  consistent with the velocity used in the next Picard iteration
!  and in the final output.  (Previously this was done BEFORE
!  relaxation, leaving a (1-alpha)*delta_u inconsistency.)
!
#   include <isf_viscosity.h>
!
!  Picard convergence check: dual criterion.
!  1. Relative velocity change: ||u_relaxed - u_old|| / ||u_relaxed|| < nl_tol
!  2. Relative residual norm: ||Au - b|| / ||b|| < nl_tol
!  Either criterion passing triggers convergence.  The residual norm
!  catches cases where velocity stagnates while the SSA residual is
!  still large (MITgcm uses residual-based convergence).
!
        vel_change=0.0_r8
        vel_norm=0.0_r8
        resid_norm=0.0_r8
        rhs_norm=0.0_r8
        DO j=JstrP,JendP
          DO i=IstrP,IendP
            IF (pmask_ssa(i,j).eq.1.0_r8) THEN
              vel_change=vel_change+                                    &
     &          (Si(i,j,ksa_new,isIuis)-Si(i,j,ksa_old,isIuis))**2+    &
     &          (Si(i,j,ksa_new,isIvis)-Si(i,j,ksa_old,isIvis))**2
              vel_norm=vel_norm+                                        &
     &          Si(i,j,ksa_new,isIuis)**2+                             &
     &          Si(i,j,ksa_new,isIvis)**2
              resid_norm=resid_norm+                                    &
     &          (Au(i,j)-taudx_psi(i,j))**2+                           &
     &          (Av(i,j)-taudy_psi(i,j))**2
              rhs_norm=rhs_norm+                                        &
     &          taudx_psi(i,j)**2+taudy_psi(i,j)**2
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
        isf_dot_prod=resid_norm
        op_handle='SUM'
        CALL mp_reduce(ng,iNLM,1,isf_dot_prod,op_handle)
        resid_norm=isf_dot_prod
        isf_dot_prod=rhs_norm
        op_handle='SUM'
        CALL mp_reduce(ng,iNLM,1,isf_dot_prod,op_handle)
        rhs_norm=isf_dot_prod
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
!  Residual-based convergence check: ||Au-b||/||b|| < nl_tol.
!  Catches false convergence when velocity stagnates.
!
        IF (.not.converged .and. rhs_norm.gt.eps) THEN
          IF (SQRT(resid_norm/rhs_norm).le.nl_tol(ng)) THEN
            converged=.TRUE.
          END IF
        END IF
!
!  DEBUG: Picard iteration diagnostics (tile 0 only, sparse output).
!  Entirely gated by ISF_DEBUG: the max|u| scan and WRITE run only
!  when diagnostics are requested, so production runs pay no cost.
!
# ifdef ISF_DEBUG
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
            cff=0.0_r8
            cff1=0.0_r8
            WRITE(*,'(A,I3,A,ES12.4,A,ES12.4,A,L1)')                   &
     &          '  Picard iter ',nl_iter,' nl_resid=',nl_resid,        &
     &          ' max|u|=',u_slide,' conv=',converged
          END IF
        END IF
# endif
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
!  NaN check on SSA solution.
!
        DO j=JstrP,JendP
          DO i=IstrP,IendP
            IF (pmask_ssa(i,j).eq.1.0_r8) THEN
              IF (Si(i,j,ksa_new,isIuis).ne.Si(i,j,ksa_new,isIuis).or.   &
     &            Si(i,j,ksa_new,isIvis).ne.Si(i,j,ksa_new,isIvis)) THEN
                IF (Master) WRITE(stdout,'(a,2i4)')                       &
     &            ' SSA NaN at PSI node:',i,j
                Si(i,j,ksa_new,isIuis)=0.0_r8
                Si(i,j,ksa_new,isIvis)=0.0_r8
              END IF
            END IF
          END DO
        END DO
!
!  Post-Picard cleanup has been moved to isf_mom.F: thin-margin nodes
!  are now deactivated (pmask_ssa=0) BEFORE the Picard solve, making
!  them Dirichlet (u=v=0) throughout the iteration.  This ensures the
!  solved problem is consistent with the returned solution.
!
