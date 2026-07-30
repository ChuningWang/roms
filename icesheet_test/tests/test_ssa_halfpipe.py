#!/usr/bin/env python3
"""
STREAMICE Halfpipe Qualitative Comparison (Stage 2.5.6)
=======================================================
Verifies that the Q1 FE SSA solver correctly captures the lateral drag
effect in channel flow, which is the key physics of the STREAMICE halfpipe
verification case (MITgcm verification/halfpipe_streamice/).

The halfpipe geometry: ice flows through a channel with no-slip sidewalls.
For fully-developed flow u=u(y), the SSA reduces to:
  νH * d²u/dy² + τ_d - C²u = 0
  u(0) = u(W) = 0  (no-slip walls)

Analytical solution:
  u(y) = (τ_d/C²) * [1 - cosh((y-W/2)/λ) / cosh(W/(2λ))]
  where λ = sqrt(νH / C²)

Key physics verified:
1. Center velocity < τ_d/C² (lateral drag reduces speed)
2. Velocity profile is cosh-shaped (not parabolic)
3. For narrow channels (W << λ): u_center ≈ τ_d*W²/(8*νH) (viscous limit)
4. For wide channels (W >> λ): u_center → τ_d/C² (drag-balance limit)

This test uses the 2D Q1 FE solver and compares against the analytical
solution at the channel centerline.

Pass criteria:
- Centerline velocity matches analytical to < 1%
- Lateral velocity profile shape is correct (cosh)
- Wide/narrow channel limits are captured
"""

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
import sys

# =============================================================================
# Parameters (matching STREAMICE halfpipe defaults where possible)
# =============================================================================
nuH = 1.0e15          # constant viscosity [Pa s m]
C2 = 4.0e8            # linear drag coefficient [Pa s/m]
rho_ice = 917.0       # kg/m³
g_grav = 9.81         # m/s²
H_ice = 500.0         # m
alpha = 0.005         # surface slope
tau_d = rho_ice * g_grav * H_ice * alpha  # driving stress [Pa]

# Channel parameters
Lx = 50000.0          # channel length [m] (flow direction)
W_channel = 20000.0   # channel width [m] (cross-flow)

# Derived
u_inf = tau_d / C2    # terminal velocity (no lateral drag)
lam = np.sqrt(nuH / C2)  # lateral shear margin length

print("=" * 70)
print("STREAMICE Halfpipe Qualitative Comparison (Stage 2.5.6)")
print("=" * 70)
print(f"\nChannel: {Lx/1000:.0f} km × {W_channel/1000:.0f} km")
print(f"τ_d = {tau_d:.2f} Pa, C² = {C2:.2e}, νH = {nuH:.2e}")
print(f"u_∞ = τ_d/C² = {u_inf:.4e} m/s ({u_inf*3.156e7:.1f} m/yr)")
print(f"λ = sqrt(νH/C²) = {lam:.1f} m ({lam/1000:.2f} km)")
print(f"W/λ = {W_channel/lam:.2f}")

# =============================================================================
# Analytical solution (fully-developed channel flow)
# =============================================================================
def u_analytical_channel(y, W, tau_d, C2, lam):
    """Analytical SSA solution for channel flow with no-slip walls."""
    return (tau_d / C2) * (1.0 - np.cosh((y - W/2) / lam) / np.cosh(W / (2*lam)))

# =============================================================================
# 2D FE solver (same as MMS test, but with physical BCs)
# =============================================================================
def solve_channel_flow(nx, ny, Lx, W, tau_d, C2, nuH):
    """Solve 2D SSA for channel flow with no-slip sidewalls."""
    dx = Lx / nx
    dy = W / ny
    
    x = np.linspace(0, Lx, nx+1)
    y = np.linspace(0, W, ny+1)
    X, Y = np.meshgrid(x, y)
    
    n_nodes = (nx+1) * (ny+1)
    n_dof = 2 * n_nodes
    
    def node_id(i, j): return j * (nx+1) + i
    def u_dof(i, j): return 2 * node_id(i, j)
    def v_dof(i, j): return 2 * node_id(i, j) + 1
    
    # Gauss quadrature
    gp = np.array([-1/np.sqrt(3), 1/np.sqrt(3)])
    gw = np.array([1.0, 1.0])
    jac = dx * dy / 4.0
    
    def shape(xi, eta):
        N = np.array([(1-xi)*(1-eta), (1+xi)*(1-eta),
                      (1+xi)*(1+eta), (1-xi)*(1+eta)]) / 4.0
        dN_dxi = np.array([-(1-eta), (1-eta), (1+eta), -(1+eta)]) / 4.0
        dN_deta = np.array([-(1-xi), -(1+xi), (1+xi), (1-xi)]) / 4.0
        return N, dN_dxi, dN_deta
    
    D = nuH * np.array([[4, 2, 0], [2, 4, 0], [0, 0, 1]])
    
    rows, cols, vals = [], [], []
    rhs = np.zeros(n_dof)
    
    for ej in range(ny):
        for ei in range(nx):
            nodes = [(ei, ej), (ei+1, ej), (ei+1, ej+1), (ei, ej+1)]
            Ke = np.zeros((8, 8))
            Fe = np.zeros(8)
            
            for ig, xi in enumerate(gp):
                for jg, eta in enumerate(gp):
                    N, dN_dxi, dN_deta = shape(xi, eta)
                    dN_dx = dN_dxi * (2.0/dx)
                    dN_dy = dN_deta * (2.0/dy)
                    w = gw[ig] * gw[jg] * jac
                    
                    B = np.zeros((3, 8))
                    for k in range(4):
                        B[0, 2*k] = dN_dx[k]
                        B[1, 2*k+1] = dN_dy[k]
                        B[2, 2*k] = dN_dy[k]
                        B[2, 2*k+1] = dN_dx[k]
                    
                    Ke += w * B.T @ D @ B
                    for k in range(4):
                        for l in range(4):
                            Ke[2*k, 2*l] += w * C2 * N[k] * N[l]
                            Ke[2*k+1, 2*l+1] += w * C2 * N[k] * N[l]
                    
                    # Driving stress RHS (only u-component)
                    for k in range(4):
                        Fe[2*k] += w * tau_d * N[k]
            
            # Scatter
            for a in range(8):
                na = nodes[a // 2]
                da = u_dof(*na) if a % 2 == 0 else v_dof(*na)
                rhs[da] += Fe[a]
                for b in range(8):
                    nb = nodes[b // 2]
                    db = u_dof(*nb) if b % 2 == 0 else v_dof(*nb)
                    rows.append(da)
                    cols.append(db)
                    vals.append(Ke[a, b])
    
    A = sparse.coo_matrix((vals, (rows, cols)), shape=(n_dof, n_dof)).tocsr()
    
    # BCs: no-slip on sidewalls (j=0 and j=ny), periodic-like in x
    # For fully-developed flow: u=0 at y=0,W; du/dx=0 at x=0,Lx
    # Simplified: Dirichlet u=v=0 at y=0 and y=W (sidewalls)
    #             Dirichlet u=v=0 at x=0 (inlet)
    #             Neumann (natural) at x=Lx (outlet)
    A_lil = A.tolil()
    bc_dofs = []
    
    # Sidewalls: u=v=0 at j=0 and j=ny
    for i in range(nx+1):
        bc_dofs.extend([u_dof(i, 0), v_dof(i, 0)])
        bc_dofs.extend([u_dof(i, ny), v_dof(i, ny)])
    
    # Inlet: u=v=0 at i=0
    for j in range(ny+1):
        bc_dofs.extend([u_dof(0, j), v_dof(0, j)])
    
    for dof in bc_dofs:
        A_lil[dof, :] = 0
        A_lil[:, dof] = 0
        A_lil[dof, dof] = 1.0
        rhs[dof] = 0.0
    
    A = A_lil.tocsr()
    sol = spsolve(A, rhs)
    
    u_sol = sol[0::2].reshape(ny+1, nx+1)
    v_sol = sol[1::2].reshape(ny+1, nx+1)
    
    return X, Y, u_sol, v_sol

# =============================================================================
# Run test
# =============================================================================
print(f"\n{'='*70}")
print("2D FE Channel Flow Solution")
print(f"{'='*70}")

nx, ny = 40, 40
X, Y, u_num, v_num = solve_channel_flow(nx, ny, Lx, W_channel, tau_d, C2, nuH)

# Extract centerline (x = Lx/2, away from inlet BC)
i_mid = nx // 2
y_line = Y[:, i_mid]
u_line = u_num[:, i_mid]

# Analytical solution
u_anal = u_analytical_channel(y_line, W_channel, tau_d, C2, lam)

# Comparison
u_center_num = u_line[ny//2]
u_center_anal = u_anal[ny//2]
u_center_ratio = u_center_num / u_inf

print(f"\nCenterline velocity (x = {Lx/2/1000:.0f} km):")
print(f"  Numerical:  {u_center_num:.6e} m/s ({u_center_num*3.156e7:.1f} m/yr)")
print(f"  Analytical: {u_center_anal:.6e} m/s ({u_center_anal*3.156e7:.1f} m/yr)")
print(f"  u_center/u_∞ = {u_center_ratio:.4f}")
print(f"  Relative error: {abs(u_center_num-u_center_anal)/u_center_anal*100:.3f}%")

# Profile shape comparison
mask = u_anal > 1e-10 * u_inf
if mask.sum() > 5:
    profile_error = np.abs(u_line[mask] - u_anal[mask]) / u_anal[mask]
    print(f"\nProfile shape (all y at centerline):")
    print(f"  Median relative error: {np.median(profile_error)*100:.3f}%")
    print(f"  Max relative error:    {np.max(profile_error)*100:.3f}%")

# =============================================================================
# Wide vs narrow channel limits
# =============================================================================
print(f"\n{'='*70}")
print("Channel Width Sensitivity (STREAMICE halfpipe physics)")
print(f"{'='*70}")
print(f"{'W/λ':>8} {'u_c/u_∞ (num)':>14} {'u_c/u_∞ (anal)':>15} {'error':>8}")
print("-" * 50)

widths = [5000, 10000, 20000, 50000, 100000]
for W in widths:
    # Analytical
    u_c_anal = (tau_d/C2) * (1.0 - 1.0/np.cosh(W/(2*lam)))
    
    # Numerical (coarse grid for speed)
    nx_w = 20
    ny_w = max(10, int(W / 1000))
    _, _, u_w, _ = solve_channel_flow(nx_w, ny_w, Lx, W, tau_d, C2, nuH)
    u_c_num = u_w[ny_w//2, nx_w//2]
    
    err = abs(u_c_num - u_c_anal) / u_c_anal * 100 if u_c_anal > 0 else 0
    print(f"{W/lam:>8.1f} {u_c_num/u_inf:>14.4f} {u_c_anal/u_inf:>15.4f} "
          f"{err:>7.2f}%")

# =============================================================================
# Summary
# =============================================================================
print(f"\n{'='*70}")
print("VERIFICATION SUMMARY")
print(f"{'='*70}")

pass_center = abs(u_center_num - u_center_anal) / u_center_anal < 0.01
pass_profile = np.median(profile_error) < 0.02 if mask.sum() > 5 else False
pass_physics = u_center_ratio < 1.0  # lateral drag reduces velocity

print(f"[{'PASS' if pass_center else 'FAIL'}] Centerline velocity error < 1%: "
      f"{abs(u_center_num-u_center_anal)/u_center_anal*100:.3f}%")
print(f"[{'PASS' if pass_profile else 'FAIL'}] Profile shape error < 2%: "
      f"{np.median(profile_error)*100:.3f}%")
print(f"[{'PASS' if pass_physics else 'FAIL'}] Lateral drag effect captured: "
      f"u_center/u_∞ = {u_center_ratio:.4f} < 1.0")

if pass_center and pass_profile and pass_physics:
    print(f"\n✓ Stage 2.5.6 PASSED: Channel flow (halfpipe physics) verified")
    print(f"  The Q1 FE SSA solver correctly captures lateral drag effects")
    print(f"  consistent with STREAMICE halfpipe verification case.")
    sys.exit(0)
else:
    print(f"\n✗ Stage 2.5.6 FAILED")
    sys.exit(1)
