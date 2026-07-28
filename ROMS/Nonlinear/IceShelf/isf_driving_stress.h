!  isf_driving_stress.h - SSA driving stress RHS via FE assembly.
!
!  Computes the consistent (FE-integrated) nodal force vector for the
!  SSA driving stress using 2x2 Gauss quadrature on each RHO cell:
!
!    b_k = ∫∫ τ_d · φ_k dA
!
!  where τ_d = -ρ_i g H ∇S is the driving stress vector, H and S are
!  interpolated from the 4 RHO corner nodes using Q1 shape functions,
!  and φ_k are the bilinear basis functions at PSI velocity nodes.
!
!  The surface elevation S is smoothed before gradient computation to
!  regularize the grounding line "kink" (discontinuity in dS/dn between
!  grounded S=-h+H and floating S=(1-ρ_i/ρ_w)H surfaces). Without
!  smoothing, the FE gradient produces spurious driving stress at the
!  grounding line that is O(100x) too large.
!
!  Smoothing: N iterations of 5-point Laplacian filter on zs_smooth.
!  Reference: PISM compute_nuH_everywhere (slope regularization),
!             STREAMICE streamice_driving_stress.F
!
!  Initialize FE-integrated nodal forces to zero.
!
      DO j=JstrP,JendP
        DO i=IstrP,IendP
          taudx_psi(i,j)=0.0_r8
          taudy_psi(i,j)=0.0_r8
        END DO
      END DO
!
      rhoi_g=rho_i*gravity
!
!-----------------------------------------------------------------------
!  Smooth surface elevation for driving stress gradient computation.
!  This regularizes the grounding line kink where zs transitions from
!  grounded (S = -h + H) to floating (S = (1-rho_i/rho_w)*H).
!-----------------------------------------------------------------------
!
!  Initialize smoothed field.
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          IF (Si(i,j,ksa_old,isIhis).gt.0.0_r8) THEN
            zs_smooth(i,j)=Fi(i,j,icIzsf)
          ELSE
            zs_smooth(i,j)=0.0_r8
          END IF
        END DO
      END DO
!
!  Apply 5-point Laplacian smoothing (N iterations).
!  Only smooth ice-covered nodes; ice-free neighbors use center value.
!
      n_smooth=5
      DO is_iter=1,n_smooth
        DO j=Jstr,Jend
          DO i=Istr,Iend
            IF (Si(i,j,ksa_old,isIhis).gt.0.0_r8) THEN
              zs_sum=zs_smooth(i,j)
              n_ice=1
              IF (Si(i+1,j,ksa_old,isIhis).gt.0.0_r8) THEN
                zs_sum=zs_sum+zs_smooth(i+1,j)
                n_ice=n_ice+1
              ELSE
                zs_sum=zs_sum+zs_smooth(i,j)
                n_ice=n_ice+1
              END IF
              IF (Si(i-1,j,ksa_old,isIhis).gt.0.0_r8) THEN
                zs_sum=zs_sum+zs_smooth(i-1,j)
                n_ice=n_ice+1
              ELSE
                zs_sum=zs_sum+zs_smooth(i,j)
                n_ice=n_ice+1
              END IF
              IF (Si(i,j+1,ksa_old,isIhis).gt.0.0_r8) THEN
                zs_sum=zs_sum+zs_smooth(i,j+1)
                n_ice=n_ice+1
              ELSE
                zs_sum=zs_sum+zs_smooth(i,j)
                n_ice=n_ice+1
              END IF
              IF (Si(i,j-1,ksa_old,isIhis).gt.0.0_r8) THEN
                zs_sum=zs_sum+zs_smooth(i,j-1)
                n_ice=n_ice+1
              ELSE
                zs_sum=zs_sum+zs_smooth(i,j)
                n_ice=n_ice+1
              END IF
              zs_smooth(i,j)=zs_sum/REAL(n_ice,r8)
            END IF
          END DO
        END DO
      END DO
!
!  Exchange halo for smoothed field (needed for tile-boundary elements).
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_p2d_tile(ng,tile,LBi,UBi,LBj,UBj,zs_smooth)
      END IF
# ifdef DISTRIBUTE
      CALL mp_exchange2d(ng,tile,iNLM,1,LBi,UBi,LBj,UBj,              &
     &                   NghostPoints,EWperiodic(ng),NSperiodic(ng),    &
     &                   zs_smooth)
# endif
!
!-----------------------------------------------------------------------
!  Element loop: assemble driving stress RHS for each active RHO cell.
!  RHO cell (i,j) has 4 PSI corner nodes where velocity is defined.
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          IF (Si(i,j,ksa_old,isIhis).gt.0.0_r8) THEN
            dx_rho=1.0_r8/pm(i,j)
            dy_rho=1.0_r8/pn(i,j)
            jac=1.0_r8/(pm(i,j)*pn(i,j))
!
!  Get smoothed surface elevation for the 4 corner nodes.
!  Ice-free nodes are assigned the mean of ice-covered corners.
!
            zs_sum=0.0_r8
            n_ice=0
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
              H_corner(k_node)=Si(i1,j1,ksa_old,isIhis)
              IF (H_corner(k_node).gt.0.0_r8) THEN
                zs_corner(k_node)=zs_smooth(i1,j1)
                zs_sum=zs_sum+zs_corner(k_node)
                n_ice=n_ice+1
              ELSE
                zs_corner(k_node)=0.0_r8
              END IF
            END DO
            IF (n_ice.gt.0) THEN
              zs_mean=zs_sum/REAL(n_ice,r8)
            ELSE
              zs_mean=0.0_r8
            END IF
            DO k_node=1,4
              IF (H_corner(k_node).le.0.0_r8) THEN
                zs_corner(k_node)=zs_mean
              END IF
            END DO
!
            DO nq=1,4
!
!  Interpolate H and compute ∇S at quadrature point from 4 RHO
!  corner values using Q1 shape function derivatives.
!
              H_q=0.0_r8
              dSdx_q=0.0_r8
              dSdy_q=0.0_r8
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
                H_q=H_q+Si(i1,j1,ksa_old,isIhis)*phi_val(k_node,nq)
                dSdx_q=dSdx_q+zs_corner(k_node)*                      &
     &              DPhi(k_node,nq,1)*2.0_r8/dx_rho
                dSdy_q=dSdy_q+zs_corner(k_node)*                      &
     &              DPhi(k_node,nq,2)*2.0_r8/dy_rho
              END DO
!
!  Scatter driving force to 4 corner PSI nodes.
!  b_k += w_q * |J| * (-ρ_i g H ∇S) · φ_k
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
                  phiv=phi_val(k_node,nq)
                  taudx_psi(i1,j1)=taudx_psi(i1,j1)-                  &
     &              jacq(nq)*jac*rhoi_g*H_q*dSdx_q*phiv
                  taudy_psi(i1,j1)=taudy_psi(i1,j1)-                  &
     &              jacq(nq)*jac*rhoi_g*H_q*dSdy_q*phiv
                END IF
              END DO  ! k_node
            END DO  ! nq
          END IF  ! H > 0
        END DO  ! i
      END DO  ! j
!
!  Exchange halo for assembled RHS (needed for shared PSI nodes
!  at tile boundaries).
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_p2d_tile(ng,tile,LBi,UBi,LBj,UBj,taudx_psi)
        CALL exchange_p2d_tile(ng,tile,LBi,UBi,LBj,UBj,taudy_psi)
      END IF
# ifdef DISTRIBUTE
      CALL mp_exchange2d(ng,tile,iNLM,2,LBi,UBi,LBj,UBj,              &
     &                   NghostPoints,EWperiodic(ng),NSperiodic(ng),    &
     &                   taudx_psi,taudy_psi)
# endif
!
!  Store FE-integrated driving force to diagnostic output fields.
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          Fi(i,j,icItdx)=taudx_psi(i,j)
          Fi(i,j,icItdy)=taudy_psi(i,j)
        END DO
      END DO
