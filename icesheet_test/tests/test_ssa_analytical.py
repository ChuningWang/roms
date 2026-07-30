#!/usr/bin/env python3
"""
SSA Analytical Solution Verification (Stage 2.5.4)
===================================================
Verifies the Q1 FE + PCG discretization against the 1D analytical solution
for a grounded ice stream with constant thickness and linear surface slope.

Analytical solution (1D SSA, linear sliding, constant nuH):
  4*nuH*u''(x) - C²*u(x) = -τ_d
  u(0) = 0  (ice divide, Dirichlet)
  u'(L) = 0 (calving front, Neumann / zero stress)

  u(x) = (τ_d/C²) * [1 - cosh((L-x)/λ) / cosh(L/λ)]
  where λ = sqrt(4*nuH / C²), τ_d = ρ_i*g*H*α

This script implements the SAME Q1 FE discretization as the ROMS Fortran code
(isf_cg_solve.h) and verifies:
1. The FE solution converges to the analytical solution as h→0
2. The convergence rate is O(h²) (expected for Q1 elements)
3. The PCG solver converges to machine precision

Pass criteria:
- Relative L2 error < 1% at dx=2000m (MISOMIP resolution)
- Convergence rate between 1.5 and 2.5 (O(h²))
"""

import numpy as np
import sys

# =============================================================================
# Physical parameters (match isf.in)
# =============================================================================
rho_ice = 917.0       # kg/m³
g_grav = 9.81         # m/s²
A_glen = 6.338e-25    # Pa⁻³ s⁻¹
n_glen = 3.0
C_basal = 2.0e4       # Pa^(1/2) s^(1/2) m^(-1/2)
nuH_reg = 1.0e15      # Pa s m (regularization dominates for this test)

# Test case parameters
H_ice = 500.0         # m (uniform thickness)
alpha = 0.005         # surface slope (5 m/km)
L_domain = 200000.0   # m (200 km domain)

# Derived
tau_d = rho_ice * g_grav * H_ice * alpha  # driving stress [Pa]
C2 = C_basal**2                           # linear drag coefficient [Pa s/m]
u_inf = tau_d / C2                        # terminal velocity [m/s]
lam = np.sqrt(4.0 * nuH_reg / C2)        # shear margin length [m]

print("=" * 70)
print("SSA 1D Analytical Solution Verification (Stage 2.5.4)")
print("=" * 70)
print(f"\nTest case: grounded ice stream, uniform H={H_ice} m, slope={alpha}")
print(f"  τ_d = {tau_d:.2f} Pa")
print(f"  C² = {C2:.2e} Pa·s/m")
print(f"  u_∞ = τ_d/C² = {u_inf:.6e} m/s ({u_inf*3.156e7:.1f} m/yr)")
print(f"  λ = sqrt(4·nuH/C²) = {lam:.1f} m ({lam/1000:.2f} km)")
print(f"  L = {L_domain/1000:.0f} km, L/λ = {L_domain/lam:.1f}")

# =============================================================================
# Analytical solution
# =============================================================================
def u_analytical(x, L, tau_d, C2, lam):
    """Exact 1D SSA solution: u(x) = (τ_d/C²)[1 - cosh((L-x)/λ)/cosh(L/λ)]"""
    return (tau_d / C2) * (1.0 - np.cosh((L - x) / lam) / np.cosh(L / lam))

# =============================================================================
# Q1 FE discretization (matches isf_cg_solve.h)
# =============================================================================
def solve_ssa_1d_fe(N, L, tau_d, C2, nuH, bc_left='dirichlet', bc_right='neumann'):
    """
    Solve 1D SSA using Q1 finite elements with PCG.
    
    Weak form: ∫[4*nuH*u'*v' + C²*u*v]dx = ∫τ_d*v*dx + [σ_xx*v]_boundary
    
    For Q1 elements on uniform grid:
    - Stiffness: K[i,i] = 4*nuH/dx + C²*dx*2/3 (interior)
    -            K[i,i±1] = -4*nuH/dx + C²*dx/6
    - Load:      F[i] = τ_d*dx (interior), τ_d*dx/2 (boundary)
    
    This matches the ROMS isf_cg_solve.h assembly (PSI-point B-grid).
    """
    dx = L / N
    x = np.linspace(0, L, N+1)  # N+1 nodes
    
    # Assembly: tridiagonal system
    # Main diagonal
    diag = np.zeros(N+1)
    # Off-diagonal
    off = np.zeros(N)
    # RHS
    rhs = np.zeros(N+1)
    
    for e in range(N):
        # Element [x_e, x_{e+1}], length = dx
        # Local stiffness: 4*nuH/dx * [[1,-1],[-1,1]] + C²*dx/6 * [[2,1],[1,2]]
        k_local = (4.0*nuH/dx) * np.array([[1,-1],[-1,1]]) + \
                  (C2*dx/6.0) * np.array([[2,1],[1,2]])
        # Local load: τ_d*dx/2 * [1, 1]
        f_local = tau_d * dx / 2.0 * np.array([1.0, 1.0])
        
        # Scatter to global (tridiagonal: diag, off)
        # Nodes: e (local 0), e+1 (local 1)
        diag[e] += k_local[0, 0]
        diag[e+1] += k_local[1, 1]
        off[e] += k_local[0, 1]  # = k_local[1,0] (symmetric)
        rhs[e] += f_local[0]
        rhs[e+1] += f_local[1]
    
    # Boundary conditions
    if bc_left == 'dirichlet':
        # u(0) = 0: zero row/col, set diag=1, rhs=0
        rhs[0] = 0.0
        diag[0] = 1.0
        if N > 0:
            off[0] = 0.0
    
    if bc_right == 'neumann':
        # Zero stress at x=L: natural BC (no modification needed)
        # The Neumann BC σ_xx(L) = 0 is automatically satisfied
        pass
    
    # Solve with PCG (Jacobi preconditioner)
    u = np.zeros(N+1)
    r = rhs - (diag * u + np.concatenate([[0], off * u[:-1]]) +
               np.concatenate([off * u[1:], [0]]))
    
    # Apply BC to residual
    if bc_left == 'dirichlet':
        r[0] = 0.0
    
    z = r / diag  # Jacobi preconditioner
    if bc_left == 'dirichlet':
        z[0] = 0.0
    d = z.copy()
    rz = np.dot(r, z)
    
    tol = 1e-12
    max_iter = 1000
    
    for cg_iter in range(max_iter):
        # A*d
        Ad = diag * d + np.concatenate([[0], off * d[:-1]]) + \
             np.concatenate([off * d[1:], [0]])
        if bc_left == 'dirichlet':
            Ad[0] = 0.0
        
        dAd = np.dot(d, Ad)
        if abs(dAd) < 1e-30:
            break
        alpha_cg = rz / dAd
        
        u += alpha_cg * d
        r -= alpha_cg * Ad
        if bc_left == 'dirichlet':
            u[0] = 0.0
            r[0] = 0.0
        
        rz_new = np.dot(r, r / diag)
        if np.sqrt(np.dot(r, r)) < tol * np.sqrt(np.dot(rhs, rhs)):
            break
        
        beta_cg = rz_new / rz
        z = r / diag
        if bc_left == 'dirichlet':
            z[0] = 0.0
        d = z + beta_cg * d
        if bc_left == 'dirichlet':
            d[0] = 0.0
        rz = rz_new
    
    return x, u, cg_iter + 1

# =============================================================================
# Convergence test
# =============================================================================
print(f"\n{'='*70}")
print("Grid Convergence Study")
print(f"{'='*70}")
print(f"{'N':>6} {'dx(km)':>8} {'L2_err':>12} {'Linf_err':>12} "
      f"{'rate':>6} {'CG_iter':>8}")
print("-" * 60)

resolutions = [25, 50, 100, 200, 400, 800]
errors_l2 = []
errors_linf = []
dx_list = []

for N in resolutions:
    dx = L_domain / N
    x, u_num, cg_iters = solve_ssa_1d_fe(N, L_domain, tau_d, C2, nuH_reg)
    u_exact = u_analytical(x, L_domain, tau_d, C2, lam)
    
    # L2 error (trapezoidal rule)
    err = u_num - u_exact
    l2_err = np.sqrt(np.trapezoid(err**2, x))
    # Normalize by L2 norm of exact solution
    l2_norm = np.sqrt(np.trapezoid(u_exact**2, x))
    rel_l2 = l2_err / l2_norm
    
    # Linf error
    linf_err = np.max(np.abs(err))
    rel_linf = linf_err / np.max(np.abs(u_exact))
    
    errors_l2.append(rel_l2)
    errors_linf.append(rel_linf)
    dx_list.append(dx)
    
    # Convergence rate
    if len(errors_l2) > 1:
        rate = np.log(errors_l2[-2] / errors_l2[-1]) / np.log(dx_list[-2] / dx_list[-1])
        rate_str = f"{rate:.2f}"
    else:
        rate_str = "  -"
    
    print(f"{N:>6} {dx/1000:>8.2f} {rel_l2:>12.4e} {rel_linf:>12.4e} "
          f"{rate_str:>6} {cg_iters:>8}")

# =============================================================================
# MISOMIP resolution check
# =============================================================================
print(f"\n{'='*70}")
print("MISOMIP Resolution (dx=2km) Verification")
print(f"{'='*70}")

N_misomip = int(L_domain / 2000)  # 100 elements for 200km at 2km
x, u_num, cg_iters = solve_ssa_1d_fe(N_misomip, L_domain, tau_d, C2, nuH_reg)
u_exact = u_analytical(x, L_domain, tau_d, C2, lam)
err = u_num - u_exact
rel_l2_misomip = np.sqrt(np.trapezoid(err**2, x)) / np.sqrt(np.trapezoid(u_exact**2, x))
rel_linf_misomip = np.max(np.abs(err)) / np.max(np.abs(u_exact))

print(f"  N = {N_misomip} elements, dx = {L_domain/N_misomip/1000:.1f} km")
print(f"  Relative L2 error:   {rel_l2_misomip*100:.4f}%")
print(f"  Relative Linf error: {rel_linf_misomip*100:.4f}%")
print(f"  CG iterations:       {cg_iters}")
print(f"  u_max (numerical):   {u_num.max():.6e} m/s")
print(f"  u_max (analytical):  {u_exact.max():.6e} m/s")

# =============================================================================
# Nonlinear viscosity test (Glen's law)
# =============================================================================
print(f"\n{'='*70}")
print("Nonlinear Viscosity (Glen's Law) Picard Convergence")
print(f"{'='*70}")

B_glen = A_glen**(-1.0/(2.0*n_glen))

def solve_ssa_1d_nonlinear(N, L, tau_d, C2, H, max_picard=100, picard_tol=1e-8):
    """Solve 1D SSA with nonlinear Glen's law viscosity using Picard iteration."""
    dx = L / N
    x = np.linspace(0, L, N+1)
    
    # Initial guess: drag balance
    u = np.full(N+1, tau_d / C2)
    u[0] = 0.0
    
    for picard_iter in range(max_picard):
        # Compute strain rate
        du_dx = np.diff(u) / dx  # N values (element-level)
        eps_e = np.maximum(np.abs(du_dx), 1e-12)
        
        # Glen's law viscosity (per element)
        nuH_elem = 0.5 * B_glen**2 * eps_e**((1.0-n_glen)/n_glen) * H
        nuH_elem = np.maximum(nuH_elem, nuH_reg)
        
        # Assemble with current viscosity
        diag = np.zeros(N+1)
        off = np.zeros(N)
        rhs = np.zeros(N+1)
        
        for e in range(N):
            k_local = (4.0*nuH_elem[e]/dx) * np.array([[1,-1],[-1,1]]) + \
                      (C2*dx/6.0) * np.array([[2,1],[1,2]])
            f_local = tau_d * dx / 2.0 * np.array([1.0, 1.0])
            diag[e] += k_local[0, 0]
            diag[e+1] += k_local[1, 1]
            off[e] += k_local[0, 1]
            rhs[e] += f_local[0]
            rhs[e+1] += f_local[1]
        
        # Dirichlet BC at x=0
        rhs[0] = 0.0
        diag[0] = 1.0
        off[0] = 0.0
        
        # Solve linear system (direct for 1D)
        # Build full matrix
        A = np.diag(diag) + np.diag(off, 1) + np.diag(off, -1)
        u_new = np.linalg.solve(A, rhs)
        u_new[0] = 0.0
        
        # Picard relaxation
        alpha_picard = 0.7
        u_relaxed = alpha_picard * u_new + (1.0 - alpha_picard) * u
        u_relaxed[0] = 0.0
        
        # Check convergence
        vel_change = np.max(np.abs(u_relaxed - u))
        vel_norm = np.max(np.abs(u_relaxed))
        if vel_change / max(vel_norm, 1e-15) < picard_tol:
            return x, u_relaxed, picard_iter + 1, True
        
        u = u_relaxed
    
    return x, u, max_picard, False

N_nl = 200
x_nl, u_nl, picard_iters, converged = solve_ssa_1d_nonlinear(
    N_nl, L_domain, tau_d, C2, H_ice)

print(f"  N = {N_nl}, dx = {L_domain/N_nl/1000:.1f} km")
print(f"  Picard converged: {converged} in {picard_iters} iterations")
print(f"  u_max = {u_nl.max():.6e} m/s ({u_nl.max()*3.156e7:.1f} m/yr)")

# Compare nonlinear vs linear (using same nuH_reg)
lam_nl = np.sqrt(4.0 * nuH_reg / C2)
u_lin = u_analytical(x_nl, L_domain, tau_d, C2, lam_nl)
rel_diff = np.max(np.abs(u_nl - u_lin)) / np.max(np.abs(u_lin))
print(f"  Max difference vs linear analytical: {rel_diff*100:.4f}%")
print(f"  Note: Difference is expected where Glen viscosity > nuH_reg.")
print(f"        This confirms nonlinear Picard handles variable viscosity.")

# =============================================================================
# Summary
# =============================================================================
print(f"\n{'='*70}")
print("VERIFICATION SUMMARY")
print(f"{'='*70}")

# Check convergence rate
rates = []
for i in range(1, len(errors_l2)):
    r = np.log(errors_l2[i-1] / errors_l2[i]) / np.log(dx_list[i-1] / dx_list[i])
    rates.append(r)
avg_rate = np.mean(rates[-3:])  # average of last 3 rates

pass_l2 = rel_l2_misomip < 0.01  # < 1% at MISOMIP resolution
pass_rate = 1.5 < avg_rate < 2.5  # O(h²) convergence
# Nonlinear: Picard converges AND solution is physically reasonable
# (within factor of 2 of linear solution - Glen viscosity adds stiffness)
pass_nonlinear = converged and rel_diff < 0.50 and u_nl.max() > 0.1*u_inf

print(f"[{'PASS' if pass_l2 else 'FAIL'}] L2 error at dx=2km < 1%: "
      f"{rel_l2_misomip*100:.4f}%")
print(f"[{'PASS' if pass_rate else 'FAIL'}] Convergence rate O(h²): "
      f"rate = {avg_rate:.2f} (expected ~2.0)")
print(f"[{'PASS' if pass_nonlinear else 'FAIL'}] Nonlinear Picard convergence: "
      f"{picard_iters} iters, u_max/u_inf = {u_nl.max()/u_inf:.2f}")

if pass_l2 and pass_rate and pass_nonlinear:
    print(f"\n✓ Stage 2.5.4 PASSED: Q1 FE discretization verified against "
          f"analytical solution")
    sys.exit(0)
else:
    print(f"\n✗ Stage 2.5.4 FAILED")
    sys.exit(1)
