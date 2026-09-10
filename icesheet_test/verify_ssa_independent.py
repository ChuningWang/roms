#!/usr/bin/env python3
"""Independent verification of the ROMS isf_ssa first-timestep velocity field.

Everything here is derived from the continuum SSA weak form; no ROMS code
is transcribed.  Inputs are the raw grid/initial/boundary files plus the
isf.in parameters.  The solve is a fresh Q1-FE assembly with a sparse
direct linear solver and Picard iteration on the Glen viscosity.

Weak form (u-row test phi):
  a(u,w) = int 2*nuH*[(2 exx^u + eyy^u) exx^w
                     + (2 eyy^u + exx^u) eyy^w
                     + 2 exy^u exy^w] dA + int beta*(u.w) dA
  rhs_k  = int tau_d . phi_k dA   (mass-lumped, tau_d = -rhoi*g*H_psi*grad S)
"""
import netCDF4 as nc
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl

# ---------------- parameters (from isf.in / defaults) ----------------
n_glen   = 3.0
A_glen   = 6.338e-25
B_hard   = A_glen**(-1.0/n_glen)
eps2_gm  = (1.0e-12)**2          # eps_glen_min^2 added to the invariant
nuH_reg  = 1.0e15
H_min    = 10.0
C_basal  = 2.0e4
beta_g   = C_basal**2            # n_basal = 1  ->  constant beta when grounded
rho_i    = 916.7
rho_w    = 1025.0
grav     = 9.81
Dcrit    = 0.10
n_smooth = 20
n_taper  = 3
eps_f    = 1e-10

# ---------------- input fields ----------------
grd = nc.Dataset('calving_test/misomip_grd_calvmask.nc')
ini = nc.Dataset('data/misomip_ini.nc')
bry = nc.Dataset('data/roms_bry.nc')
his = nc.Dataset('outputs_smoke/roms_his.nc')

# RHO grid is [eta=40, xi=400]; PSI nodes are [39, 399] (interior corners).
h   = np.array(grd.variables['h'][:])            # bed depth, positive down
pm  = np.array(grd.variables['pm'][:])
pn  = np.array(grd.variables['pn'][:])
H   = np.array(ini.variables['is_h'][0])
# global index layout: RHO (j=0..39, i=0..399); cell (j,i) j=1..38,i=1..398
# PSI (jp=0..38, ip=0..398) is the SW corner of RHO cell (j=i+1 indexes).

# ---------------- masks (AND convention) ----------------
ice   = H > 0.0                                   # rmask_ist
grounded = (rho_i*H > rho_w*(h - Dcrit)) & ice    # beta flotation test, zeta=0

def pmask_of(mask):
    P = np.zeros((39, 399))
    P[:] = (mask[0:39,0:399] & mask[0:39,1:400] &
            mask[1:40,0:399] & mask[1:40,1:400])
    return P
Pmask = pmask_of(ice)

# ---------------- surface elevation and smoothing ----------------
zsf = np.where(grounded, -h + Dcrit + eps_f + H, (1.0-rho_i/rho_w)*H)
zsf = np.where(ice, zsf, 0.0)
zs = zsf.copy()
for _ in range(n_smooth):
    zsum = zs.copy()
    for (dj,di) in ((0,1),(0,-1),(1,0),(-1,0)):
        nb = np.roll(np.roll(ice, dj, axis=0), di, axis=1)
        zv = np.roll(np.roll(zs , dj, axis=0), di, axis=1)
        add = np.where(nb, zv, zs)
        zsum = zsum + add
    zs = np.where(ice, zsum/5.0, 0.0)

# ---------------- driving stress RHS at PSI nodes ----------------
tau_x = np.zeros((39, 399)); tau_y = np.zeros((39, 399))
for jp in range(39):
    for ip in range(399):
        if Pmask[jp,ip]:
            dist = n_taper
            for dj in (-1,0,1):
                for di in (-1,0,1):
                    if (dj or di):
                        jj, ii = jp+dj, ip+di
                        if 0 <= jj < 39 and 0 <= ii < 399 and Pmask[jj,ii] == 0:
                            dist = min(dist, max(abs(dj),abs(di)))
            taper = min(dist/n_taper, 1.0)
            dxr = 0.5*(1.0/pm[jp,ip] + 1.0/pm[jp+1,ip])
            dyr = 0.5*(1.0/pn[jp,ip] + 1.0/pn[jp,ip+1])
            jac = 0.25*sum(1.0/(pm[a,b]*pn[a,b]) for (a,b) in
                           ((jp,ip),(jp,ip+1),(jp+1,ip),(jp+1,ip+1)))
            dSdx = ((zs[jp,ip+1]+zs[jp+1,ip+1]) -
                    (zs[jp,ip]  +zs[jp+1,ip]  ))/(2*dxr)
            dSdy = ((zs[jp+1,ip]+zs[jp+1,ip+1]) -
                    (zs[jp,ip]  +zs[jp,ip+1]  ))/(2*dyr)
            Hpsi = 0.25*sum(H[a,b] for (a,b) in
                            ((jp,ip),(jp,ip+1),(jp+1,ip),(jp+1,ip+1)))
            cff = -taper*rho_i*grav*Hpsi*jac
            tau_x[jp,ip] = cff*dSdx
            tau_y[jp,ip] = cff*dSdy

# ---------------- beta field ----------------
beta = np.where(grounded, beta_g, 0.0)

# ---------------- quadrature data (derived here) ----------------
# Q1 on the reference square [-1,1]^2, 2x2 Gauss.
gq = 1.0/np.sqrt(3.0)
qp = [(-gq,-gq),(gq,-gq),(-gq,gq),(gq,gq)]
def shape(xi, eta):
    # returns (phi[4], dphi_dxi[4], dphi_deta[4]) node order SW,SE,NW,NE
    s = [(-1,-1),(1,-1),(-1,1),(1,1)]
    phi = []; dx_ = []; de_ = []
    for (sx, se) in s:
        phi.append(0.25*(1+sx*xi)*(1+se*eta))
        dx_.append(0.25*sx*(1+se*eta))
        de_.append(0.25*se*(1+sx*xi))
    return np.array(phi), np.array(dx_), np.array(de_)

# ---------------- dof numbering ----------------
NP = 39*399
def dof(jp, ip, c):
    return (jp*399 + ip)*2 + c

rows = []; cols = []; vals = []
bvec = np.zeros(2*NP)

# Dirichlet: every boundary and inactive node carries u = v = 0,
# exactly what the model does (Stage 1 zeroing; the clamped boundary# data is not imposed because those nodes are inactive).
fixed = np.zeros(2*NP, dtype=bool); fval = np.zeros(2*NP)
for jp in range(39):
    for ip in range(399):
        if Pmask[jp,ip] == 0:
            for c in (0,1):
                fixed[dof(jp,ip,c)] = True

active = ~fixed

def nuH_of(u, v):
    nuh = np.zeros((40,400))
    for j in range(1,39):
        for i in range(1,399):
            if H[j,i] > 0:
                # PSI corners of RHO cell (j,i): (j-1,i-1),(j-1,i),(j,i-1),(j,i)
                exx = 0.5*((u[j,i-1]-u[j-1,i-1]) +
                           (u[j,i]-u[j-1,i]))*pm[j,i]
                eyy = 0.5*((v[j-1,i]-v[j-1,i-1]) +
                           (v[j,i]-v[j,i-1]))*pn[j,i]
                exy = 0.5*(0.5*((v[j,i-1]-v[j-1,i-1]) +
                                (v[j,i]-v[j-1,i]))*pm[j,i] +
                          0.5*((u[j-1,i]-u[j-1,i-1]) +
                               (u[j,i]-u[j,i-1]))*pn[j,i])
                gam = exx*exx + eyy*eyy + exx*eyy + exy*exy + eps2_gm
                nuh[j,i] = (0.5*B_hard*gam**((1.0-n_glen)/(2.0*n_glen))
                            * max(H[j,i], H_min) + nuH_reg)
    return nuh

def assemble(nuh):
    rows = []; cols = []; vals = []
    b = np.zeros(2*NP)
    b.reshape(39,399,2)[...,0] = tau_x
    b.reshape(39,399,2)[...,1] = tau_y
    for j in range(1,39):
        for i in range(1,399):
            if nuh[j,i] <= 0:
                continue
            dx = 1.0/pm[j,i]; dy = 1.0/pn[j,i]; area = dx*dy
            nds = [(j-1,i-1),(j-1,i),(j,i-1),(j,i)]   # psi (jp,ip) SW,SE,NW,NE
            for (xi,eta) in qp:
                phi, dxi, deta = shape(xi,eta)
                dpx = dxi*2.0/dx                       # dphi/dx
                dpy = deta*2.0/dy
                # strain contribution of basis a in components c:
                # c=0 (u): exx=dpx, eyy=0, exy=0.5*dpy
                # c=1 (v): exx=0, eyy=dpy, exy=0.5*dpx
                w = 0.25*area
                nu = nuh[j,i]
                for a in range(4):
                    for ca in (0,1):
                        if ca == 0: Exa, Eya, Ega = dpx[a], 0.0, 0.5*dpy[a]
                        else:       Exa, Eya, Ega = 0.0, dpy[a], 0.5*dpx[a]
                        da = dof(*nds[a], ca)
                        for bt in range(4):
                            for cb in (0,1):
                                if cb == 0: Exb, Eyb, Egb = dpx[bt], 0.0, 0.5*dpy[bt]
                                else:       Exb, Eyb, Egb = 0.0, dpy[bt], 0.5*dpx[bt]
                                k = w*2.0*nu*((2*Exa+Eya)*Exb +
                                              (2*Eya+Exa)*Eyb + 2*Ega*Egb)
                                if ca == cb:
                                    k += w*beta[j,i]*phi[a]*phi[bt]
                                rows.append(da); cols.append(dof(*nds[bt],cb))
                                vals.append(k)
    K = sp.coo_matrix((vals,(rows,cols)), shape=(2*NP,2*NP)).tocsr()
    return K, b

# ---------------- Picard with direct solves ----------------
uprev = np.zeros(2*NP)
for it in range(400):
    u3 = uprev.reshape(39,399,2)
    nuh = nuH_of(u3[...,0].copy(), u3[...,1].copy())
    K, b = assemble(nuh)
    # Dirichlet elimination
    rhs = b - K @ (fval*fixed)
    A = K.tolil()
    idx = np.where(fixed)[0]
    for d in idx:
        A[d,:] = 0.0; A[d,d] = 1.0
    A = A.tocsr()
    rhs[idx] = fval[idx]
    usol = spl.spsolve(A, rhs)
    alpha = 0.7                      # same under-relaxation as the model
    usol = alpha*usol + (1.0-alpha)*uprev
    rel = np.linalg.norm(usol-uprev)/max(np.linalg.norm(usol),1e-30)
    uprev = usol
    if rel < 1e-9:
        break
print('picard iterations:', it+1, 'rel change:', rel)

# ---------------- compare with model record 1 ----------------
mu0 = his.variables['is_u'][1]; mv0 = his.variables['is_v'][1]
print('model is_u shape:', np.shape(mu0))
mu = np.ma.filled(mu0, np.nan); mv = np.ma.filled(mv0, np.nan)
uind = uprev.reshape(39,399,2)

# compare on model-valid, active nodes
ok = (~np.isnan(mu)) & (Pmask > 0)
du = uind[...,0][ok] - mu[ok]
dv = uind[...,1][ok] - mv[ok]
un = np.linalg.norm(mu[ok]); vn = np.linalg.norm(mv[ok])
print('compare on %d nodes' % ok.sum())
print('u: rel L2 diff = %.3e   max|diff| = %.3e   |u|max model = %.3e'
      % (np.linalg.norm(du)/un, np.abs(du).max(), np.nanmax(np.abs(mu))))
print('v: rel L2 diff = %.3e   max|diff| = %.3e   |v|max model = %.3e'
      % (np.linalg.norm(dv)/vn, np.abs(dv).max(), np.nanmax(np.abs(mv))))
