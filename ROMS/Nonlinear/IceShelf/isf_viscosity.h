!
!  isf_viscosity.h - Compute SSA effective viscosity and basal drag.
!
!  Computes at RHO points:
!    nuH_rho = depth-integrated effective viscosity (Glen law)
!    beta_rho = basal drag coefficient (Weertman sliding law)
!
!  Includes PISM robustness measures:
!    - Additive nuH regularization (nuH += nuH_reg)
!    - Strength extension for thin ice (H < min_thickness)
!
!  Reference: STREAMICE streamice_visc_beta.F, PISM compute_nuH_everywhere
!
!  Initialize.
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          nuH_rho(i,j)=0.0_r8
          beta_rho(i,j)=0.0_r8
        END DO
      END DO
!
      B_glen_sq=B_glen(ng)*B_glen(ng)
      C_basal_sq=C_basal(ng)*C_basal(ng)
!
!  Compute viscosity and basal drag at RHO points where ice is present.
!  RHO cell (i,j) has 4 PSI corner nodes: (i,j), (i+1,j), (i,j+1),
!  (i+1,j+1).
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          H_rho=Si(i,j,ksa_old,isIhis)
          IF (H_rho.gt.0.0_r8) THEN
!
!  Get 4 corner velocities (PSI points).
!
            u1=Si(i,j,ksa_new,isIuis)
            u2=Si(i+1,j,ksa_new,isIuis)
            u3=Si(i,j+1,ksa_new,isIuis)
            u4=Si(i+1,j+1,ksa_new,isIuis)
            v1=Si(i,j,ksa_new,isIvis)
            v2=Si(i+1,j,ksa_new,isIvis)
            v3=Si(i,j+1,ksa_new,isIvis)
            v4=Si(i+1,j+1,ksa_new,isIvis)
!
!  Velocity at cell center (4-corner average).
!
            u_mid=0.25_r8*(u1+u2+u3+u4)
            v_mid=0.25_r8*(v1+v2+v3+v4)
!
!  Strain rates at RHO point (centered differences from 4 corners).
!
            eps_xx=0.5_r8*((u2-u1)+(u4-u3))*pm(i,j)
            eps_yy=0.5_r8*((v3-v1)+(v4-v2))*pn(i,j)
            eps_xy=0.5_r8*(0.5_r8*((v2-v1)+(v4-v3))*pm(i,j)+             &
     &                     0.5_r8*((u3-u1)+(u4-u2))*pn(i,j))
!
!  Glen law effective viscosity.
!    eps_e = sqrt(eps_xx^2 + eps_yy^2 + eps_xx*eps_yy + eps_xy^2
!                + eps_glen_min^2)
!    nuH = 0.5 * B^2 * eps_e^((1-n)/n) * H
!
            eps_e=SQRT(eps_xx*eps_xx+eps_yy*eps_yy+                      &
     &                 eps_xx*eps_yy+eps_xy*eps_xy+                      &
     &                 eps_glen_min(ng)*eps_glen_min(ng))
!
            IF (H_rho.lt.min_thickness(ng)) THEN
!
!  PISM strength extension: constant nuH for thin ice to maintain
!  ellipticity at the margin.
!
              nuH_rho(i,j)=0.5_r8*B_glen_sq*min_thickness(ng)*          &
     &            eps_e**((1.0_r8-n_glen(ng))/n_glen(ng))
            ELSE
              nuH_rho(i,j)=0.5_r8*B_glen_sq*                            &
     &            eps_e**((1.0_r8-n_glen(ng))/n_glen(ng))*H_rho
            END IF
!
!  PISM additive regularization.
!
            nuH_rho(i,j)=nuH_rho(i,j)+nuH_reg(ng)
!
!  Weertman basal drag: beta = C^2 * |u|^(m-1) * (1 - float_frac)
!  Only grounded ice has basal drag (float_frac = 0 -> full drag,
!  float_frac = 1 -> no drag).
!
            u_slide=SQRT(u_mid*u_mid+v_mid*v_mid+                       &
     &                   eps_u_min(ng)*eps_u_min(ng))
            beta_rho(i,j)=C_basal_sq*                                   &
     &        u_slide**(n_basal(ng)-1.0_r8)*                            &
     &        (1.0_r8-float_frac(i,j))
!
          END IF  ! H > 0
        END DO
      END DO
!
!  Store to diagnostic output fields.
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          Fi(i,j,icInuh)=nuH_rho(i,j)
          Fi(i,j,icIbta)=beta_rho(i,j)
        END DO
      END DO
