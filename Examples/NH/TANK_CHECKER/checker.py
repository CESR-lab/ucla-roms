# Growth of the 2dx checkerboard in zeta (and w) for a tank run directory with joined tank_his*.nc
import sys, glob, numpy as np, netCDF4 as nc
d=sys.argv[1]
import re; files=sorted(f for f in glob.glob(d+'/tank_his*.nc') if not re.search(r'\.[0-9]\.nc$',f))
rows=[]
for f in files:
    ds=nc.Dataset(f); t=ds['ocean_time'][:]
    for r in range(len(t)):
        z=np.ma.filled(ds['zeta'][r][1:-1,1:-1],0.); w=np.ma.filled(ds['w'][r][:,1:-1,1:-1],0.)
        ny,nx=z.shape; J,I=np.meshgrid(np.arange(ny),np.arange(nx),indexing='ij'); s=1-2*((I+J)%2)
        A=(z*s).mean(); rows.append((t[r],A,np.sqrt((z**2).mean()),abs(z).max(),abs(w).max(),np.sqrt((w**2).mean())))
print(f"{'t[s]':>6} {'checker amp [m]':>15} {'rms zeta':>10} {'max|zeta|':>10} {'max|w|':>10} {'rms w':>10}")
for t,A,rz,mz,mw,rw in rows: print(f"{t:6.1f} {A:15.3e} {rz:10.3e} {mz:10.3e} {mw:10.3e} {rw:10.3e}")
