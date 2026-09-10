#!/usr/bin/env python3
"""Independent verification of the first SSA solve (corrected geometry).

ROMS convention: the PSI value stored at python (jp,ip) is the SW
corner of RHO cell (jp,ip); its four surrounding cells are rows
jp-1..jp and cols ip-1..ip.  Element (cell) (j,i) has corners at
psi slots (j,i),(j,i+1),(j+1,i),(j+1,i+1)  [SW,SE,NW,NE].
"""
import netCDF4 as nc
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl

n_glen, A_glen = 3.0, 6.338e-25
B_hard = A_glen**(-1.0/n_glen)
eps2_gm, nuH_reg, H_min = 1e-24, 1.0e15, 10.0
C_basal = 2.0e4; beta_g = C_basal**2
rho_i, rho_w, grav, Dcrit = 916.7, 1025.0, 9.81, 0.10
n_smooth, n_taper = 20, 3

grd = nc.Dataset('calving_test/misomip_grd_calvmask.nc')
ini = nc.Dataset('data/misomip_ini.nc')
his = nc.Dataset('outputs_smoke/roms_his.nc')
h  = np.array(grd.variables['h'][:]); pm = np.array(grd.variables['pm'][:]); pn = np.array(grd.variables['pn'][:])
H  = np.array(ini.variables['is_h'][0])
um = np.ma.filled(his.variables['is_u'][1], 0.0)
vm = np.ma.filled(his.variables['is_v'][1], 0.0)
T  = np.ma.filled(his.variables['is_taudx'][1], 0.0)[:39,:399]
TY = np.ma.filled(his.variables['is_taudy'][1], 0.0)[:39,:399]

ice = H > 0
grounded = (rho_i*H > rho_w*(h-Dcrit)) & ice
# PSI node (jp,ip): cells rows jp-1..jp, cols ip-1..ip
P = (ice[0:39,0:399] & ice[0:39,1:400] & ice[1:40,0:399] & ice[1:40,1:400])
# NOTE the mask expression equals the shifted-window boolean only after
# relabeling: cells (jp-1,ip-1),(jp-1,ip),(jp,ip-1),(jp,ip) with node
# label shifted by (-1,-1); define it explicitly instead:
P = np.zeros((39,399), dtype=bool)
P[:] = (ice[0:39,0:399] & ice[0:39,1:400] & ice[1:40,0:399] & ice[1:40,1:400])
Pshift = np.zeros((39,399), dtype=bool)
for jp in range(39):
    for ip in range(399):
        c = ice[max(jp-1,0):jp+1, max(ip-1,0):ip+1]
        Pshift[jp,ip] = c.all() if c.size==4 else False

# ---- surface elevation, smoothing ----
zsf = np.where(grounded, -h+Dcrit+H, (1-rho_i/rho_w)*H)
zs = np.where(ice, zsf, 0.0)
for _ in range(n_smooth):
    zsum = zs.copy()
    for (dj,di) in ((0,1),(0,-1),(1,0),(-1,0)):
        nb = np.roll(np.roll(ice,dj,0),di,1); zv = np.roll(np.roll(zs,dj,0),di,1)
        zsum = zsum + np.where(nb, zv, zs)
    zs = np.where(ice, zsum/5.0, 0.0)

dxr_ = 1.0/pm; dyr_ = 1.0/pn; A_ = (1.0/pm)*(1.0/pn)

# ---- RHS at node (jp,ip) from cells rows jp-1..jp, cols ip-1..ip ----
tau_x = np.zeros((39,399)); tau_y = np.zeros((39,399))
for jp in range(39):
    for ip in range(399):
        if Pshift[jp,ip]:
            dist = n_taper
            for dj in (-1,0,1):
                for di in (-1,0,1):
                    if (dj or di):
                        jj,ii = jp+dj, ip+di
                        if 0<=jj<39 and 0<=ii<399 and not Pshift[jj,ii]:
                            dist = min(dist, max(abs(dj),abs(di)))
            taper = min(dist/n_taper, 1.0)
            # node (jp,ip): cells (jp-1,ip-1),(jp-1,ip),(jp,ip-1),(jp,ip)
            dxr = 0.5*(dxr_[jp,ip-1]+dxr_[jp,ip])
            dyr = 0.5*(dyr_[jp-1,ip]+dyr_[jp,ip])
            jac = 0.25*(A_[jp-1,ip-1]+A_[jp-1,ip]+A_[jp,ip-1]+A_[jp,ip])
            dSx = ((zs[jp,ip]+zs[jp-1,ip])-(zs[jp,ip-1]+zs[jp-1,ip-1]))/(2*dxr)
            dSy = ((zs[jp,ip]+zs[jp,ip-1])-(zs[jp-1,ip]+zs[jp-1,ip-1]))/(2*dyr)
            Hpsi = 0.25*(H[jp-1,ip-1]+H[jp-1,ip]+H[jp,ip-1]+H[jp,ip])
            cff = -taper*rho_i*grav*Hpsi*jac
            tau_x[jp,ip] = cff*dSx
            tau_y[jp,ip] = cff*dSy

print('=== RHS check (all active nodes, corrected incidence):')
print('tau_x vs T : corr=%+.5f relL2=%.3e' %
      (np.corrcoef(tau_x[Pshift],T[Pshift])[0,1],
       np.linalg.norm(tau_x[Pshift]-T[Pshift])/np.linalg.norm(T[Pshift])))
print('tau_y vs TY: corr=%+.5f relL2=%.3e' %
      (np.corrcoef(tau_y[Pshift],TY[Pshift])[0,1],
       np.linalg.norm(tau_y[Pshift]-TY[Pshift])/np.linalg.norm(TY[Pshift])))

# ---- assembly: cell (j,i) corners psi (j,i),(j,i+1),(j+1,i),(j+1,i+1) ----
NP = 39*399
def dof(jp,ip,c): return (jp*399+ip)*2+c
gq = 1/np.sqrt(3); qp=[(-gq,-gq),(gq,-gq),(-gq,gq),(gq,gq)]
sc=[(-1,-1),(1,-1),(-1,1),(1,1)]
def shp(xi,eta):
    return (np.array([0.25*(1+sx*xi)*(1+se*eta) for sx,se in sc]),
            np.array([0.25*sx*(1+se*eta) for sx,se in sc]),
            np.array([0.25*se*(1+sx*xi) for sx,se in sc]))

beta = np.where(grounded, beta_g, 0.0)

def nuH_of(u, v):
    nuh = np.zeros((40,400))
    for j in range(1,39):
        for i in range(1,399):
            if H[j,i] > 0:
                # cell (j,i) corners: psi (j-1..j, i-1..i)
                exx = 0.5*((u[j,i]-u[j,i-1]) + (u[j-1,i]-u[j-1,i-1]))*pm[j,i]
                eyy = 0.5*((v[j-1,i]-v[j-1,i-1]) + (v[j,i]-v[j,i-1]))*pn[j,i]
                exy = 0.5*(0.5*((v[j,i]-v[j,i-1]) + (v[j-1,i]-v[j-1,i-1]))*pm[j,i] +
                          0.5*((u[j-1,i]-u[j-1,i-1]) + (u[j,i]-u[j,i-1]))*pn[j,i])
                gam = exx*exx+eyy*eyy+exx*eyy+exy*exy+eps2_gm
                nuh[j,i] = 0.5*B_hard*gam**((1.0-n_glen)/(2.0*n_glen)) \
                           *max(H[j,i],H_min)+nuH_reg
    return nuh

def assemble(nuh):
    rows=[];cols=[];vals=[]
    b = np.zeros(2*NP)
    b.reshape(39,399,2)[...,0] = tau_x
    b.reshape(39,399,2)[...,1] = tau_y
    for j in range(1,39):
        for i in range(1,399):
            if nuh[j,i] <= 0: continue
            dx=1/pm[j,i]; dy=1/pn[j,i]; area=dx*dy
            nds=[(j-1,i-1),(j-1,i),(j,i-1),(j,i)]   # SW,SE,NW,NE
            for (xi,eta) in qp:
                phi,dxi,deta = shp(xi,eta)
                dpx=dxi*2/dx; dpy=deta*2/dy
                Ex=list(dpx); Ey=list(dpy); Gx=[0.5*x for x in dpy]; Gy=[0.5*x for x in dpx]
                w=0.25*area; nu=nuh[j,i]
                for a in range(4):
                    for ca in (0,1):
                        Exa=Ex[a] if ca==0 else 0.0; Eya=Ey[a] if ca==1 else 0.0
                        Ga=Gx[a] if ca==0 else Gy[a]
                        for bt in range(4):
                            for cb in (0,1):
                                Exb=Ex[bt] if cb==0 else 0.0; Eyb=Ey[bt] if cb==1 else 0.0
                                Gb=Gx[bt] if cb==0 else Gy[bt]
                                k=w*2*nu*((2*Exa+Eya)*Exb+(2*Eya+Exa)*Eyb+2*Ga*Gb)
                                if ca==cb: k+=w*beta[j,i]*phi[a]*phi[bt]
                                rows.append(dof(*nds[a],ca)); cols.append(dof(*nds[bt],cb)); vals.append(k)
    return sp.coo_matrix((vals,(rows,cols)),shape=(2*NP,2*NP)).tocsr(), b

fixed = np.zeros(2*NP, dtype=bool)
for jp in range(39):
    for ip in range(399):
        if not Pshift[jp,ip]:
            fixed[dof(jp,ip,0)] = True; fixed[dof(jp,ip,1)] = True
idx = np.where(fixed)[0]

uprev = np.zeros(2*NP)
for it in range(400):
    u3 = uprev.reshape(39,399,2)
    nuh = nuH_of(u3[...,0].copy(), u3[...,1].copy())
    K, b = assemble(nuh)
    rhs = b - K @ (uprev*fixed)
    A = K.tolil()
    for d0 in idx: A[d0,:] = 0.0; A[d0,d0] = 1.0
    A = A.tocsr(); rhs[idx] = 0.0
    usol = spl.spsolve(A, rhs)
    usol = 0.7*usol + 0.3*uprev
    rel = np.linalg.norm(usol-uprev)/max(np.linalg.norm(usol),1e-30)
    uprev = usol
    if rel < 1e-8: break
print('picard: %d iterations, rel=%.2e' % (it+1, rel))

uind = uprev.reshape(39,399,2)
du = uind[...,0][Pshift]-um[Pshift]; dv = uind[...,1][Pshift]-vm[Pshift]
print('=== independent solve vs model (all %d active nodes):' % Pshift.sum())
print('u: relL2=%.4e max|d|=%.3e (|u|max=%.3e)' %
      (np.linalg.norm(du)/np.linalg.norm(um[Pshift]), np.abs(du).max(), np.abs(um[Pshift]).max()))
print('v: relL2=%.4e max|d|=%.3e' %
      (np.linalg.norm(dv)/np.linalg.norm(vm[Pshift]), np.abs(dv).max()))
