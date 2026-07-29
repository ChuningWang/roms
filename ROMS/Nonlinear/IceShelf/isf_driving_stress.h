!  isf_driving_stress.h - SSA driving stress RHS computation.
!
!  Computes the nodal force vector for the SSA driving stress using
!  finite differences on a smoothed surface elevation field:
!
!    b_k = -rho_i * g * H * dS/dx * (nodal area)
!
!  where S is the smoothed surface elevation, H is ice thickness
!  interpolated to PSI nodes, and nodal area = jac/4 (mass-lumped Q1).
!
!  The surface elevation S is smoothed (20 iterations of 5-point
!  Laplacian) before gradient computation to regularize the grounding
!  line "kink" (discontinuity in dS/dn between grounded S=-h+H and
!  floating S=(1-rho_i/rho_w)*H surfaces).
!
!  Ice margin taper: driving stress is linearly tapered to zero within
!  n_taper cells of the ice margin (where any neighbor has H=0).
!  The missing resistance should be provided by an ice front pressure
!  boundary condition (future work).
!
!  Known limitation: the smoothing creates small artifacts at MPI tile
!  boundaries due to the local nature of the Laplacian filter. This
!  affects the diagnostic output but has minimal impact on the solution.
!
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
!  Initialize smoothed field (including halo for tile-boundary elements).
!
      DO j=LBj,UBj
        DO i=LBi,UBi
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
      n_smooth=20
      DO is_iter=1,n_smooth
# ifdef DISTRIBUTE
        CALL mp_exchange2d(ng,tile,iNLM,1,LBi,UBi,LBj,UBj,            &
     &                     NghostPoints,EWperiodic(ng),NSperiodic(ng),  &
     &                     zs_smooth)
# endif
        DO j=Jstr-1,Jend+1
          DO i=Istr-1,Iend+1
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
!  Compute driving stress RHS at PSI points using finite differences
!  on the smoothed surface elevation. This avoids FE assembly artifacts
!  at tile boundaries and ice margins.
!
!  At PSI point (i,j), the 4 surrounding RHO points are:
!    (i-1,j-1), (i,j-1), (i-1,j), (i,j)
!
!  dS/dx = [(zs(i,j-1)+zs(i,j)) - (zs(i-1,j-1)+zs(i-1,j))] / (2*dx)
!  dS/dy = [(zs(i-1,j)+zs(i,j)) - (zs(i-1,j-1)+zs(i,j-1))] / (2*dy)
!
!  Nodal area = jac/4 (mass-lumped Q1, 4 elements share each node).
!-----------------------------------------------------------------------
!
      n_taper=3
      DO j=JstrP,JendP
        DO i=IstrP,IendP
          IF (pmask_ssa(i,j).eq.1.0_r8) THEN
!
!  Get H and zs from 4 surrounding RHO points.
!
            H11=Si(i-1,j-1,ksa_old,isIhis)
            H21=Si(i,  j-1,ksa_old,isIhis)
            H12=Si(i-1,j,  ksa_old,isIhis)
            H22=Si(i,  j,  ksa_old,isIhis)
!
!  Skip if any surrounding RHO point is ice-free (margin handling).
!
            IF (H11.gt.0.0_r8.and.H21.gt.0.0_r8.and.                  &
     &          H12.gt.0.0_r8.and.H22.gt.0.0_r8) THEN
!
!  Ice margin taper: check 8 neighbors for ice-free cells.
!
              dist_margin=n_taper
              DO dj=-1,1
                DO di=-1,1
                  IF (di.ne.0.or.dj.ne.0) THEN
                    IF (Si(i+di,j+dj,ksa_old,isIhis).le.0.0_r8) THEN
                      dist_margin=MIN(dist_margin,                    &
     &                    MAX(ABS(di),ABS(dj)))
                    END IF
                  END IF
                END DO
              END DO
              taper_w=MIN(REAL(dist_margin,r8)/                       &
     &                    REAL(n_taper,r8),1.0_r8)
!
!  FD gradient of smoothed surface at PSI point.
!
              dx_rho=0.5_r8*(1.0_r8/pm(i-1,j)+1.0_r8/pm(i,j))
              dy_rho=0.5_r8*(1.0_r8/pn(i,j-1)+1.0_r8/pn(i,j))
              jac_avg=0.25_r8*(1.0_r8/(pm(i-1,j-1)*pn(i-1,j-1))+      &
     &                         1.0_r8/(pm(i,  j-1)*pn(i,  j-1))+      &
     &                         1.0_r8/(pm(i-1,j  )*pn(i-1,j  ))+      &
     &                         1.0_r8/(pm(i,  j  )*pn(i,  j  )))
!
              dSdx=((zs_smooth(i,j-1)+zs_smooth(i,j))-                 &
     &              (zs_smooth(i-1,j-1)+zs_smooth(i-1,j)))/            &
     &             (2.0_r8*dx_rho)
              dSdy=((zs_smooth(i-1,j)+zs_smooth(i,j))-                 &
     &              (zs_smooth(i-1,j-1)+zs_smooth(i,j-1)))/            &
     &             (2.0_r8*dy_rho)
!
!  Average H to PSI point.
!
              H_psi=0.25_r8*(H11+H21+H12+H22)
!
!  Driving stress RHS (FE-consistent nodal force).
!  b_k = -rho_i * g * H * dS/dx * (nodal area)
!  nodal area = jac/4 for mass-lumped Q1.
!
              taudx_psi(i,j)=-taper_w*rhoi_g*H_psi*dSdx*               &
     &                       jac_avg*0.25_r8
              taudy_psi(i,j)=-taper_w*rhoi_g*H_psi*dSdy*               &
     &                       jac_avg*0.25_r8
            END IF
          END IF
        END DO
      END DO
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
