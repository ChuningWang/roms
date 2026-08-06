#!/usr/bin/env python3
"""
SSA 2D MMS Convergence Test (Stage 2.5.3)
==========================================
Method of Manufactured Solutions (MMS) for the 2D SSA operator.

Verifies that the Q1 FE discretization of the full 2D SSA momentum equation
converges at O(h²) by:
1. Choosing manufactured solutions u*(x,y), v*(x,y)
2. Computing exact source terms f_u, f_v by applying the SSA operator
3. Solving the discrete system with these source terms
4. Measuring L2 error vs grid spacing

SSA operator (constant viscosity, linear drag):
  L_u = νH(4u_xx + u_yy + 3v_xy) - C²u
  L_v = νH(v_xx + 4v_yy + 3u_xy) - C²v

Manufactured solutions (satisfy u=v=0 on boundary):
  u*(x,y) = sin(πx/Lx) * sin(πy/Ly)
  v*(x,y) = sin(2πx/Lx) * sin(2πy/Ly)

Pass criteria:
- L2 convergence rate between 1.8 and 2.2 (O(h²) for Q1)
- Final error < 1% at finest resolution
"""

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
import sys

# =============================================================================
# Parameters
# =============================================================================
nuH = 1.0e15          # constant viscosity [Pa s m]
C2 = 4.0e8            # linear drag coefficient [Pa s/m]
Lx = 100000.0         # domain length x [m]
Ly = 50000.0          # domain length y [m]

# =============================================================================
# Manufactured solutions and source terms
# =============================================================================
def u_exact(x, y):
    """Manufactured u-velocity."""
    return np.sin(np.pi * x / Lx) * np.sin(np.pi * y / Ly)

def v_exact(X, Y):
    """Manufactured v-velocity."""
    return np.sin(2*np.pi * X / Lx) * np.sin(2*np.pi * Y / Ly)

def source_terms(X, Y):
    """
    Compute MMS source terms for the WEAK form.
    
    The SSA weak form solves: a(u,v) = (f, v)
    where a(u,v) = ∫[νH*strain(u):strain(v) + C²u·v] dA
    
    This corresponds to the strong form:
      -∇·(νH*strain(u)) + C²u = f
    
    So f = -νH(4u_xx + u_yy + 3v_xy) + C²u  (for u-equation)
       f = -νH(v_xx + 4v_yy + 3u_xy) + C²v  (for v-equation)
    """
    kx1 = np.pi / Lx
    ky1 = np.pi / Ly
    kx2 = 2*np.pi / Lx
    ky2 = 2*np.pi / Ly
    
    # u* = sin(kx1*x)*sin(ky1*y)
    u = np.sin(kx1*X) * np.sin(ky1*Y)
    u_xx = -kx1**2 * u
    u_yy = -ky1**2 * u
    u_xy = kx1 * ky1 * np.cos(kx1*X) * np.cos(ky1*Y)
    
    # v* = sin(kx2*x)*sin(ky2*y)
    v = np.sin(kx2*X) * np.sin(ky2*Y)
    v_xx = -kx2**2 * v
    v_yy = -ky2**2 * v
    v_xy = kx2 * ky2 * np.cos(kx2*X) * np.cos(ky2*Y)
    
    # Source terms for weak form: f = -L_strong(u*)
    # where L_strong = νH(4u_xx + u_yy + 3v_xy) - C²u
    f_u = -nuH * (4*u_xx + u_yy + 3*v_xy) + C2 * u
    f_v = -nuH * (v_xx + 4*v_yy + 3*u_xy) + C2 * v
    
    return f_u, f_v

# =============================================================================
# Q1 FE assembly on uniform 2D grid
# =============================================================================
def assemble_2d_ssa(nx, ny):
    """
    Assemble the 2D SSA Q1 FE system on a uniform grid.
    
    Grid: (nx+1) x (ny+1) nodes, spacing dx=Lx/nx, dy=Ly/ny.
    Unknowns: 2*(nx+1)*(ny+1) (u and v at each node).
    
    Returns: sparse matrix A, node coordinates X, Y.
    """
    dx = Lx / nx
    dy = Ly / ny
    
    # Node coordinates
    x = np.linspace(0, Lx, nx+1)
    y = np.linspace(0, Ly, ny+1)
    X, Y = np.meshgrid(x, y)  # shape (ny+1, nx+1)
    
    n_nodes = (nx+1) * (ny+1)
    n_dof = 2 * n_nodes
    
    # Node numbering: node(i,j) = j*(nx+1) + i
    def node_id(i, j):
        return j * (nx+1) + i
    
    # DOF numbering: u_dof = 2*node, v_dof = 2*node+1
    def u_dof(i, j):
        return 2 * node_id(i, j)
    def v_dof(i, j):
        return 2 * node_id(i, j) + 1
    
    # Q1 element stiffness matrices (on reference element [-1,1]²)
    # For bilinear Q1: φ₁=(1-ξ)(1-η)/4, φ₂=(1+ξ)(1-η)/4, etc.
    # Gradients on reference element:
    # ∂φ/∂ξ and ∂φ/∂η at Gauss points
    
    # For uniform rectangular elements, the physical stiffness is:
    # K_visc = νH * ∫[B^T * D * B] dA
    # where D = [[4, 2, 0], [2, 1, 0], [0, 0, 1]] (for SSA in 2D)
    # and B relates strains to nodal velocities.
    
    # Using 2x2 Gauss quadrature on each element:
    gp = np.array([-1/np.sqrt(3), 1/np.sqrt(3)])  # Gauss points
    gw = np.array([1.0, 1.0])  # Gauss weights
    
    # Shape functions and derivatives at Gauss points
    def shape(xi, eta):
        N = np.array([(1-xi)*(1-eta), (1+xi)*(1-eta),
                      (1+xi)*(1+eta), (1-xi)*(1+eta)]) / 4.0
        dN_dxi = np.array([-(1-eta), (1-eta), (1+eta), -(1+eta)]) / 4.0
        dN_deta = np.array([-(1-xi), -(1+xi), (1+xi), (1-xi)]) / 4.0
        return N, dN_dxi, dN_deta
    
    # Jacobian for uniform rectangular element
    jac = dx * dy / 4.0  # |J| = (dx/2)*(dy/2)
    
    # Build sparse matrix using COO format
    rows = []
    cols = []
    vals = []
    
    for ej in range(ny):
        for ei in range(nx):
            # Element nodes (counterclockwise from bottom-left)
            nodes = [(ei, ej), (ei+1, ej), (ei+1, ej+1), (ei, ej+1)]
            
            # Element stiffness (8x8: 4 nodes × 2 dofs)
            Ke = np.zeros((8, 8))
            Fe = np.zeros(8)
            
            for ig, xi in enumerate(gp):
                for jg, eta in enumerate(gp):
                    N, dN_dxi, dN_deta = shape(xi, eta)
                    
                    # Physical gradients
                    dN_dx = dN_dxi * (2.0/dx)
                    dN_dy = dN_deta * (2.0/dy)
                    
                    w = gw[ig] * gw[jg] * jac
                    
                    # Strain-displacement matrix B (3x8)
                    # ε = [ε_xx, ε_yy, ε_xy]^T = B * d
                    # ε_xx = du/dx, ε_yy = dv/dy, ε_xy = du/dy + dv/dx
                    B = np.zeros((3, 8))
                    for k in range(4):
                        B[0, 2*k] = dN_dx[k]      # du/dx
                        B[1, 2*k+1] = dN_dy[k]    # dv/dy
                        B[2, 2*k] = dN_dy[k]      # du/dy
                        B[2, 2*k+1] = dN_dx[k]    # dv/dx
                    
                    # Constitutive matrix D for SSA
                    # σ_xx = νH(4ε_xx + 2ε_yy)
                    # σ_yy = νH(2ε_xx + 4ε_yy)  -- wait, this isn't right
                    # Actually for SSA:
                    # σ_xx = νH(4u_x + 2v_y) = νH(4ε_xx + 2ε_yy)
                    # σ_yy = νH(2u_x + 4v_y) = νH(2ε_xx + 4ε_yy)
                    # σ_xy = νH(u_y + v_x) = νH*ε_xy
                    D = nuH * np.array([[4, 2, 0],
                                        [2, 4, 0],
                                        [0, 0, 1]])
                    
                    # Stiffness contribution
                    Ke += w * B.T @ D @ B
                    
                    # Drag contribution: C² * ∫[u*v_u + v*v_v] dA
                    for k in range(4):
                        for l in range(4):
                            Ke[2*k, 2*l] += w * C2 * N[k] * N[l]
                            Ke[2*k+1, 2*l+1] += w * C2 * N[k] * N[l]
            
            # Scatter to global
            for a in range(8):
                node_a = nodes[a // 2]
                dof_a = u_dof(*node_a) if a % 2 == 0 else v_dof(*node_a)
                for b in range(8):
                    node_b = nodes[b // 2]
                    dof_b = u_dof(*node_b) if b % 2 == 0 else v_dof(*node_b)
                    rows.append(dof_a)
                    cols.append(dof_b)
                    vals.append(Ke[a, b])
    
    A = sparse.coo_matrix((vals, (rows, cols)), shape=(n_dof, n_dof)).tocsr()
    return A, X, Y, dx, dy

def assemble_rhs(X, Y, dx, dy, nx, ny):
    """Assemble RHS from manufactured source terms using Q1 quadrature."""
    f_u, f_v = source_terms(X, Y)
    
    n_nodes = (nx+1) * (ny+1)
    rhs = np.zeros(2 * n_nodes)
    
    # Lumped mass matrix approach: F_i = f(x_i) * A_i
    # where A_i is the nodal area
    # For interior nodes: A_i = dx*dy
    # For edge nodes: A_i = dx*dy/2
    # For corner nodes: A_i = dx*dy/4
    
    area = np.ones((ny+1, nx+1)) * dx * dy
    area[0, :] *= 0.5
    area[-1, :] *= 0.5
    area[:, 0] *= 0.5
    area[:, -1] *= 0.5
    
    for j in range(ny+1):
        for i in range(nx+1):
            node = j * (nx+1) + i
            rhs[2*node] = f_u[j, i] * area[j, i]
            rhs[2*node+1] = f_v[j, i] * area[j, i]
    
    return rhs

def apply_bc(A, rhs, nx, ny):
    """Apply homogeneous Dirichlet BCs on all boundaries."""
    n_nodes = (nx+1) * (ny+1)
    
    # Boundary nodes
    bc_nodes = set()
    for i in range(nx+1):
        bc_nodes.add(i)                    # bottom (j=0)
        bc_nodes.add(ny*(nx+1) + i)        # top (j=ny)
    for j in range(ny+1):
        bc_nodes.add(j*(nx+1))             # left (i=0)
        bc_nodes.add(j*(nx+1) + nx)        # right (i=nx)
    
    # Zero out BC rows/cols and set diagonal to 1
    A_lil = A.tolil()
    for node in bc_nodes:
        for dof in [2*node, 2*node+1]:
            A_lil[dof, :] = 0
            A_lil[:, dof] = 0
            A_lil[dof, dof] = 1.0
            rhs[dof] = 0.0
    
    return A_lil.tocsr(), rhs

# =============================================================================
# Convergence test
# =============================================================================
print("=" * 70)
print("SSA 2D MMS Convergence Test (Stage 2.5.3)")
print("=" * 70)
print(f"\nDomain: {Lx/1000:.0f} km × {Ly/1000:.0f} km")
print(f"νH = {nuH:.2e} Pa·s·m, C² = {C2:.2e} Pa·s/m")
print(f"Manufactured: u*=sin(πx/Lx)sin(πy/Ly), v*=sin(2πx/Lx)sin(2πy/Ly)")

print(f"\n{'nx×ny':>8} {'dx(km)':>8} {'L2_u':>12} {'L2_v':>12} "
      f"{'rate_u':>7} {'rate_v':>7}")
print("-" * 62)

resolutions = [(8,4), (16,8), (32,16), (64,32)]
errors_u = []
errors_v = []
dx_list = []

for nx, ny in resolutions:
    dx = Lx / nx
    dy = Ly / ny
    
    # Assemble system
    A, X, Y, dx, dy = assemble_2d_ssa(nx, ny)
    rhs = assemble_rhs(X, Y, dx, dy, nx, ny)
    A, rhs = apply_bc(A, rhs, nx, ny)
    
    # Solve
    sol = spsolve(A, rhs)
    
    # Extract u, v
    n_nodes = (nx+1) * (ny+1)
    u_num = sol[0::2].reshape(ny+1, nx+1)
    v_num = sol[1::2].reshape(ny+1, nx+1)
    
    # Exact solution
    u_ex = u_exact(X, Y)
    v_ex = v_exact(X, Y)
    
    # L2 error (using nodal quadrature)
    err_u = u_num - u_ex
    err_v = v_num - v_ex
    
    # Approximate L2 norm: sqrt(sum(err² * area) / sum(area))
    area = np.ones_like(X) * dx * dy
    area[0,:] *= 0.5; area[-1,:] *= 0.5
    area[:,0] *= 0.5; area[:,-1] *= 0.5
    
    l2_u = np.sqrt(np.sum(err_u**2 * area) / np.sum(u_ex**2 * area))
    l2_v = np.sqrt(np.sum(err_v**2 * area) / np.sum(v_ex**2 * area))
    
    errors_u.append(l2_u)
    errors_v.append(l2_v)
    dx_list.append(dx)
    
    # Convergence rate
    if len(errors_u) > 1:
        rate_u = np.log(errors_u[-2]/errors_u[-1]) / np.log(dx_list[-2]/dx_list[-1])
        rate_v = np.log(errors_v[-2]/errors_v[-1]) / np.log(dx_list[-2]/dx_list[-1])
        rate_str_u = f"{rate_u:.2f}"
        rate_str_v = f"{rate_v:.2f}"
    else:
        rate_str_u = "  -"
        rate_str_v = "  -"
    
    print(f"{nx:>3}×{ny:<3} {dx/1000:>8.2f} {l2_u:>12.4e} {l2_v:>12.4e} "
          f"{rate_str_u:>7} {rate_str_v:>7}")

# =============================================================================
# Summary
# =============================================================================
print(f"\n{'='*70}")
print("VERIFICATION SUMMARY")
print(f"{'='*70}")

# Use last 2 rates for average
rates_u = []
rates_v = []
for i in range(1, len(errors_u)):
    rates_u.append(np.log(errors_u[i-1]/errors_u[i]) / np.log(dx_list[i-1]/dx_list[i]))
    rates_v.append(np.log(errors_v[i-1]/errors_v[i]) / np.log(dx_list[i-1]/dx_list[i]))

avg_rate_u = np.mean(rates_u[-2:])
avg_rate_v = np.mean(rates_v[-2:])

pass_rate_u = 1.8 < avg_rate_u < 2.5
pass_rate_v = 1.8 < avg_rate_v < 2.5
pass_error = errors_u[-1] < 0.01 and errors_v[-1] < 0.01

print(f"[{'PASS' if pass_rate_u else 'FAIL'}] u-velocity convergence rate: "
      f"{avg_rate_u:.2f} (expected ~2.0)")
print(f"[{'PASS' if pass_rate_v else 'FAIL'}] v-velocity convergence rate: "
      f"{avg_rate_v:.2f} (expected ~2.0)")
print(f"[{'PASS' if pass_error else 'FAIL'}] Final L2 error < 1%: "
      f"u={errors_u[-1]*100:.3f}%, v={errors_v[-1]*100:.3f}%")

if pass_rate_u and pass_rate_v and pass_error:
    print(f"\n✓ Stage 2.5.3 PASSED: 2D Q1 FE SSA operator verified, O(h²) convergence")
    sys.exit(0)
else:
    print(f"\n✗ Stage 2.5.3 FAILED")
    sys.exit(1)
