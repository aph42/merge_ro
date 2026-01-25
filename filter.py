from . import read 
import pygeode as pyg
import numpy as np, pylab as pyl

from ECMWF import ERAI as ei

import os

def rotlon(v, orig):
# {{{
   import pygeode as pyg
   
   if not v.hasaxis('lon'): return v

   l = v.lon[:]
   o = orig % 360
   off = o - orig
   lp = pyg.Lon(values=(l - o) % 360 + o - off)
   return v.replace_axes(lon=lp).sorted('lon')
# }}}

def var_filter(v, maxwn = 30, maxfreq = 1., nfft = 128):
# {{{
   import wnf

   v = v.transpose('pres', 'lon', 'time')

   dx = np.pi/180.
   dt = 0.25

   ix = 1
   it = 2

   times = v(i_time=(0, -nfft))(hour=0).time
   itimes = (times[:]/dt).astype('i') 
   itimes -= itimes[0]
   
   tcen = v(i_time=(nfft/2, -nfft/2))(hour=0).time
   pres = v.pres

   nt = len(times)
   npres = len(v.pres)
   nks = maxwn
   nws = nfft + 1

   V = v[:]
   nans = np.sum(np.isnan(V))
   if nans > 0: raise ValueError('NaNs present: %d, aborting.' % nans)

   VV = np.zeros((nt, npres, nks, nws), 'd')

   for i, t in enumerate(itimes):
      tsl = slice(t, t+nfft)
      ks, ws, vv = wnf.wnf(V[:, :, tsl], V[:, :, tsl], dx, dt, ix, it, nfft)

      VV[i, :, :, :] = vv[:, :maxwn, :]

   ak = pyg.NamedAxis(1+np.arange(maxwn), 'wn')
   aw = pyg.NamedAxis(-ws/(2*np.pi), 'freq')

   vv = pyg.Var((tcen, pres, ak, aw), values=VV, name='TT')

   return vv(freq=(-maxfreq, maxfreq))
# }}} 

def var_filter2d(v, maxwn = 30, maxfreq = 1., nfft = 128):
# {{{
   import wnf

   v = v.transpose('lon', 'time')

   dx = 2.5 * np.pi/180.
   dt = 1.

   ix = 0
   it = 1

   times = v(i_time=(0, -nfft))(hour=0).time
   itimes = (times[:] / 24).astype('i') 
   itimes -= itimes[0]
   
   tcen = v(i_time=(nfft/2, -nfft/2))(hour=0).time

   nt = len(times)
   nks = maxwn
   nws = nfft + 1

   V = v[:]
   nans = np.sum(np.isnan(V))
   if nans > 0: raise ValueError('NaNs present: %d, aborting.' % nans)

   VV = np.zeros((nt, nks, nws), 'd')

   for i, t in enumerate(itimes):
      tsl = slice(t, t+nfft)
      ks, ws, vv = wnf.wnf(V[:, tsl], V[:, tsl], dx, dt, ix, it, nfft)

      VV[i, :, :] = vv[:maxwn, :]

   ak = pyg.NamedAxis(1+np.arange(maxwn), 'wn')
   aw = pyg.NamedAxis(-ws/(2*np.pi), 'freq')

   vv = pyg.Var((tcen, ak, aw), values=VV, name='TT')

   return vv(freq=(-maxfreq, maxfreq))
# }}} 

def var_filter4d(v, maxwn = 30, maxfreq = 1., nfft = 128):
# {{{
   import wnf

   if v.hasaxis('Ds'):
      v  = v.transpose('Ds', 'pres', 'lon', 'time')
      ds = v.Ds
   elif v.hasaxis('Ts'):
      v  = v.transpose('Ts', 'pres', 'lon', 'time')
      ds = v.Ts
   else:
      raise 'Unknown dimension present'

   dx = np.pi/180.
   dt = 0.25

   ix = 2
   it = 3

   times = v(i_time=(0, -nfft))(hour=0).time
   itimes = (times[:]/dt).astype('i') 
   itimes -= itimes[0]
   
   tcen = v(i_time=(nfft/2, -nfft/2))(hour=0).time
   pres = v.pres

   nt = len(times)
   nds = len(ds)
   npres = len(pres)
   nks = maxwn
   nws = nfft + 1

   V = v[:]
   nans = np.sum(np.isnan(V))
   if nans > 0: raise ValueError('NaNs present: %d, aborting.' % nans)

   VV = np.zeros((nt, nds, npres, nks, nws), 'd')

   for i, t in enumerate(itimes):
      tsl = slice(t, t+nfft)
      ks, ws, vv = wnf.wnf(V[:, :, :, tsl], V[:, :, :, tsl], dx, dt, ix, it, nfft)

      VV[i, :, :, :, :] = vv[:, :, :maxwn, :]

   ak = pyg.NamedAxis(1+np.arange(maxwn), 'wn')
   aw = pyg.NamedAxis(-ws/(2*np.pi), 'freq')

   vv = pyg.Var((tcen, ds, pres, ak, aw), values=VV, name='TT')

   return vv(freq=(-maxfreq, maxfreq))
# }}} 

def var_filterlat(v, maxwn = 30, maxfreq = 1., nfft = 128, mlat = 10):
# {{{
   import wnf

   v  = v.transpose('pres', 'lat', 'lon', 'time')

   vN = v(lat=(0, 90))
   vS = v(lat=(-90, 0)).sorted(lat=-1).replace_axes(lat=vN.lat)
   vs = 0.5 * (vN + vS)(m_lat=(0, mlat))
   va = 0.5 * (vN - vS)(m_lat=(0, mlat))

   dx = (v.lon[1] - v.lon[0]) * np.pi/180.
   dt = 0.25

   ix = 1
   it = 2

   times = v(i_time=(0, -nfft))(hour=0).time
   itimes = (times[:]/dt).astype('i')
   itimes -= itimes[0]

   tcen = v(i_time=(nfft/2, -nfft/2))(hour=0).time
   pres = v.pres

   nt = len(times)
   npres = len(pres)
   nks = maxwn
   nws = nfft + 1

   Vs = vs[:]
   Va = va[:]
   #return Vs, Va
   nans = np.sum(np.isnan(Vs)) + np.sum(np.isnan(Va))
   if nans > 0: raise ValueError('NaNs present: %d, aborting.' % nans)

   VVs = np.zeros((nt, npres, nks, nws), 'd')
   VVa = np.zeros((nt, npres, nks, nws), 'd')

   for i, t in enumerate(itimes):
      tsl = slice(t, t+nfft)
      ks, ws, vv = wnf.wnf(Vs[:, :, tsl], Vs[:, :, tsl], dx, dt, ix, it, nfft)
      VVs[i, :, :, :] = vv[:, :maxwn, :]

      ks, ws, vv = wnf.wnf(Va[:, :, tsl], Va[:, :, tsl], dx, dt, ix, it, nfft)
      VVa[i, :, :, :] = vv[:, :maxwn, :]

   ak = pyg.NamedAxis(1+np.arange(maxwn), 'wn')
   aw = pyg.NamedAxis(-ws/(2*np.pi), 'freq')

   vvs = pyg.Var((tcen, pres, ak, aw), values=VVs, name='TTsym')
   vva = pyg.Var((tcen, pres, ak, aw), values=VVa, name='TTasym')

   return pyg.asdataset([vvs, vva])(freq=(-maxfreq, maxfreq))
# }}}

def covar_filter(u, v, maxwn = 30, maxfreq = 1., nfft = 128):
# {{{
   import wnf

   u = u.transpose('pres', 'lon', 'time')
   v = v.transpose('pres', 'lon', 'time')

   dx = np.pi/180.
   dt = 0.25

   ix = 1
   it = 2

   times = v(i_time=(0, -nfft))(hour=0).time
   itimes = (times[:]/dt).astype('i') 
   itimes -= itimes[0]
   
   tcen = v(i_time=(nfft/2, -nfft/2))(hour=0).time
   pres = v.pres

   nt = len(times)
   npres = len(v.pres)
   nks = maxwn
   nws = nfft + 1

   V = v[:]
   U = u[:]
   nans = np.sum(np.isnan(V)) + np.sum(np.isnan(U))
   if nans > 0: raise ValueError('NaNs present: %d, aborting.' % nans)

   UV = np.zeros((nt, npres, nks, nws), 'd')

   for i, t in enumerate(itimes):
      tsl = slice(t, t+nfft)
      ks, ws, uv = wnf.wnf(U[:, :, tsl], V[:, :, tsl], dx, dt, ix, it, nfft)

      UV[i, :, :, :] = uv[:, :maxwn, :]

   ak = pyg.NamedAxis(1+np.arange(maxwn), 'wn')
   aw = pyg.NamedAxis(-ws/(2*np.pi), 'freq')

   uv = pyg.Var((tcen, pres, ak, aw), values=UV, name='TQ')

   return uv
# }}} 

def do_ei_tt_filter(yr):
# {{{
   maxwn = 30
   nfft  = 128

   times = ('16 Dec %04d' % (yr - 1), '17 Jan %04d' % (yr + 1))
   pres = (900, 0.1)

   # Interpolate onto same pressure grid
   aonp = ei.on_pres(ei.an, ei.pres, ext=True)
   tei = aonp.t(time=times, m_lat=(-10, 10)).load()
   ttei = var_filter(tei, maxwn = maxwn, maxfreq=2., nfft = nfft)(year=yr)

   fn = 'tt_ei_full_%04d-v4.nc' % yr
   path = read.path + 'gridded/kw/'
   print('Writing %s (%s).' % (fn, ttei.name))
   pyg.save(path + fn, ttei, version=4, compress=True)
# }}}

def do_ei_ttsym_filter(yr):
# {{{
   maxwn = 16
   maxfreq = 1.
   nfft  = 128

   times = ('16 Dec %04d' % (yr - 1), '17 Jan %04d' % (yr + 1))
   pres = (400, 1.)

   # Interpolate onto same pressure grid
   aonp = ei.on_pres(ei.an, ei.pres, ext=True)
   tei = aonp.t(time=times, lat=(-12, 12)).load()
   ttei = var_filterlat(tei, maxwn = maxwn, maxfreq = maxfreq, nfft = nfft, mlat=10)(year=yr)

   fn = 'tts_ei_full_%04d.nc' % yr
   path = read.path + 'gridded/kw2/'
   print('Writing %s.' % (fn, ))
   pyg.save(path + fn, ttei, version=4, compress=True)
# }}}

def do_ei_uw_filter(yr):
# {{{
   maxwn = 30
   nfft  = 128

   times = ('16 Dec %04d' % (yr - 1), '17 Jan %04d' % (yr + 1))
   pres = (900, 0.1)

   # Interpolate onto same pressure grid
   aonp = ei.on_pres(ei.an, ei.pres, ext=True)
   uei = aonp.u(time=times, m_lat=(-10, 10)).load()
   wei = aonp.w(time=times, m_lat=(-10, 10)).load()
   uwei = covar_filter(uei, wei, maxwn = maxwn, maxfreq=2., nfft = nfft)(year=yr)
   uwei = uwei.rename('UW')

   fn = 'uw_ei_full_%04d-v4.nc' % yr
   path = read.path + 'gridded/kw/'
   print('Writing %s (%s).' % (fn, uwei.name))
   pyg.save(path + fn, uwei, version=4, compress=True)
# }}}

def do_mls_tt_filter(yr):
# {{{
   import MLS as mls

   maxwn = 30
   nfft  = 128

   times = ('16 Dec %04d' % (yr - 1), '17 Jan %04d' % (yr + 1))
   pres = (900, 0.1)

   # Interpolate onto same pressure grid
   tm = mls.d.t(time=times, m_lat=(-10, 10)).load()
   ttm = var_filter(tm, maxwn = maxwn, nfft = nfft)(year=yr)

   fn = 'tt_mls_full_%04d.nc' % yr
   path = read.path + 'gridded/kw/'
   print('Writing %s (%s).' % (fn, ttm.name))
   pyg.save(path + fn, ttm, version=4, compress=True)
# }}}

def do_dsens_filter(yr):
# {{{
   maxwn = 15
   nfft  = 128

   times = ('16 Dec %04d' % (yr - 1), '17 Jan %04d' % (yr + 1))
   pres = (900, 0.3)

   # Read in GPS gridded data
   ds = read.open_gridded('dsens1d')(time=times, pres=pres).load()

   # Interpolate onto same pressure grid
   aonp = ei.on_pres(ei.an, ei.pres, ext=True)
   tei = aonp.t(time=times, m_lat=(-10, 10), pres=pres).load()

   def ei_filter(v, D, T):
   # {{{
      dt = np.arange(-3., 3.1, 0.25)
      tfilt = np.exp( -(dt / T)**2)

      dl = np.arange(-30, 30, 1.)
      lfilt = np.exp( -(dl / D)**2)

      vr = pyg.concatenate([rotlon(v, -540), v, rotlon(v, 180)])

      vs = vr.smooth('lon', lfilt)(lon=(-180, 179)).smooth('time', tfilt)
      return vs.extend(0, pyg.NamedAxis(name='Ds', values=[D]))
   # }}}

   tei_filt = pyg.concatenate([ei_filter(tei, d, 0.5) for d in ds.Ds[:]])

   ttei    = var_filter  (tei,               maxwn = maxwn, maxfreq=2., nfft = nfft)(year=yr).rename('TTei')
   tteiflt = var_filter4d(tei_filt,          maxwn = maxwn, maxfreq=2., nfft = nfft)(year=yr).rename('TTeiflt')
   tteismp = var_filter4d(ds.Tei,            maxwn = maxwn, maxfreq=2., nfft = nfft)(year=yr).rename('TTeismp')
   ttgps   = var_filter4d(ds.Tbl,            maxwn = maxwn, maxfreq=2., nfft = nfft)(year=yr).rename('TTgps')

   tteg_ef = var_filter4d(tei      - ds.Tei, maxwn = maxwn, maxfreq=2., nfft = nfft)(year=yr).rename('TTeges')
   ttef_es = var_filter4d(tei_filt - ds.Tei, maxwn = maxwn, maxfreq=2., nfft = nfft)(year=yr).rename('TTefes')
   ttes_gp = var_filter4d(ds.Tei   - ds.Tbl, maxwn = maxwn, maxfreq=2., nfft = nfft)(year=yr).rename('TTesgp')

   vs = [ttei, tteiflt, tteismp, ttgps, tteg_ef, ttef_es, ttes_gp]

   fn = 'tt_dsamp1d_%04d-v5.nc' % yr
   path = read.path + 'gridded/kw2/'
   print('Writing %s.' % (fn,))
   pyg.save(path + fn, vs, version=4, compress=True)
# }}}

def do_tsens_filter(yr):
# {{{
   maxwn = 45
   nfft  = 128

   times = ('16 Dec %04d' % (yr - 1), '17 Jan %04d' % (yr + 1))
   pres = (900, 0.3)

   # Read in GPS gridded data
   ds = read.open_gridded('tsens')(time=times, pres=pres).load()

   Ts = pyg.NamedAxis(name = 'Ts', values = [1., 0.5, 0.25])
   ds = ds.replace_axes(Ds = Ts)

   # Interpolate onto same pressure grid
   aonp = ei.on_pres(ei.an, ei.pres, ext=True)
   tei = aonp.t(time=times, m_lat=(-10, 10), pres=pres).load()

   def ei_filter(v, D, T):
   # {{{
      dt = np.arange(-3., 3.1, 0.25)
      tfilt = np.exp( -(dt / T)**2)

      dl = np.arange(-30, 30, 1.)
      lfilt = np.exp( -(dl / D)**2)

      vr = pyg.concatenate([rotlon(v, -540), v, rotlon(v, 180)])

      vs = vr.smooth('lon', lfilt)(lon=(-180, 179)).smooth('time', tfilt)
      return vs.extend(0, pyg.NamedAxis(name='Ts', values=[T]))
   # }}}

   tei_filt = pyg.concatenate([ei_filter(tei, 1., d) for d in ds.Ts[:]])

   ttei    = var_filter  (tei,               maxwn = maxwn, maxfreq=2., nfft = nfft)(year=yr).rename('TTei')
   tteiflt = var_filter4d(tei_filt,          maxwn = maxwn, maxfreq=2., nfft = nfft)(year=yr).rename('TTeiflt')
   tteismp = var_filter4d(ds.Tei,            maxwn = maxwn, maxfreq=2., nfft = nfft)(year=yr).rename('TTeismp')
   ttgps   = var_filter4d(ds.Tbl,            maxwn = maxwn, maxfreq=2., nfft = nfft)(year=yr).rename('TTgps')

   tteg_ef = var_filter4d(tei      - ds.Tei, maxwn = maxwn, maxfreq=2., nfft = nfft)(year=yr).rename('TTeges')
   ttef_es = var_filter4d(tei_filt - ds.Tei, maxwn = maxwn, maxfreq=2., nfft = nfft)(year=yr).rename('TTefes')
   ttes_gp = var_filter4d(ds.Tei   - ds.Tbl, maxwn = maxwn, maxfreq=2., nfft = nfft)(year=yr).rename('TTesgp')

   vs = [ttei, tteiflt, tteismp, ttgps, tteg_ef, ttef_es, ttes_gp]

   fn = 'tt_tsamp_%04d-v5.nc' % yr
   path = read.path + 'gridded/kw2/'
   print('Writing %s.' % (fn,))
   pyg.save(path + fn, vs, version=4, compress=True)
# }}}

def do_olr_filter(yr):
# {{{
   import NOAA_OLR as olr

   maxwn  = 30
   maxfreq = 1.
   nfft  = 32

   times = ('16 Dec %04d' % (yr - 1), '17 Jan %04d' % (yr + 1))

   olr = olr.olr.olr.fill(0)(time=times, m_lat=(-10, 10)).load()
   oo = var_filter2d(olr, maxwn = maxwn, maxfreq=maxfreq, nfft = nfft)(year=yr)
   oo = oo.rename('OLR2')

   fn = 'olr_%04d.nc' % yr
   path = read.path + 'gridded/olr/'
   print('Writing %s (%s).' % (fn, oo.name))
   pyg.save(path + fn, oo, version=4, compress=True)
# }}}

def do_tt_filter(yr):
# {{{
   maxwn = 30
   maxfreq = 1.
   nfft  = 128

   times = ('16 Dec %04d' % (yr - 1), '17 Jan %04d' % (yr + 1))
   pres = (900, 0.1)

   #d = read.open_gridded('eq_rad')(time=times, pres=pres)
   d = read.open_gridded('kw2')(time=times, pres=pres)
   tei = d.Tei
   ttei = var_filter(tei, maxwn = maxwn, maxfreq = maxfreq, nfft = nfft)(year=yr)

   #fn = 'tt_%d_ei_%04d.nc' % (nfft, yr)
   fn = 'tt_ei_%04d.nc' % (yr)
   path = read.path + 'gridded/kw2/'
   print('Writing %s (%s).' % (fn, ttei.name))
   pyg.save(path + fn, ttei)

   tbl = d.Tbl
   ttbl = var_filter(tbl, maxwn = maxwn, maxfreq = maxfreq, nfft = nfft)(year=yr)

   #fn = 'tt_%d_gps_%04d.nc' % (nfft, yr)
   fn = 'tt_gps_%04d.nc' % (yr)
   path = read.path + 'gridded/kw2/'
   print('Writing %s (%s).' % (fn, ttei.name))
   pyg.save(path + fn, ttei, version=4, compress=True)
# }}}

def do_trop_tt_filter(yr):
# {{{
   maxwn = 16
   maxfreq = 1.
   nfft  = 128

   times = ('16 Dec %04d' % (yr - 1), '17 Jan %04d' % (yr + 1))
   pres = (400, 1.)

   #d = read.open_gridded('eq_rad')(time=times, pres=pres)
   d = read.open_gridded('trop')(time=times, pres=pres)

   tbl = d.Tbl
   ttbl = var_filterlat(tbl, maxwn = maxwn, maxfreq = maxfreq, nfft = nfft, mlat=10)(year=yr)

   #fn = 'tt_%d_gps_%04d.nc' % (nfft, yr)
   fn = 'tt_bl_%04d.nc' % (yr)
   path = read.path + 'gridded/kw2/'
   print('Writing %s.' % (fn,))
   pyg.save(path + fn, ttbl, version=4, compress=True)
# }}}

def do_tq_filter(yr, q = 'ei'):
# {{{
   maxwn = 30
   nfft  = 128

   times = ('16 Dec %04d' % (yr - 1), '17 Jan %04d' % (yr + 1))
   pres = (900, 0.1)

   #d = read.open_gridded('eq_rad')(time=times, pres=pres)
   d = read.open_gridded('kw2')(time=times, pres=pres)

   if q == 'ei':
      t = d.Tei
      qlw = d.lwhrei
      qsw = d.swhrei

      tql = covar_filter(t, qlw, maxwn = maxwn, nfft = nfft)(year=yr)
      tqs = covar_filter(t, qsw, maxwn = maxwn, nfft = nfft)(year=yr)

      vs = [tql, tqs]
      fn = 'tq_ei_%04d.nc' % yr
   elif q == 'gps':
      t = d.Tbl
      qlw = d.lwhr
      qsw = d.swhr
      tql = covar_filter(t, qlw, maxwn = maxwn, nfft = nfft)(year=yr)
      tqs = covar_filter(t, qsw, maxwn = maxwn, nfft = nfft)(year=yr)

      vs = [tql, tqs]
      fn = 'tq_gps_%04d.nc' % yr
   elif q == 'spc':
      t = d.Tbl
      qclw = d.lwhrco2
      qolw = d.lwhro3
      qosw = d.swhro3
      qhlw = d.lwhrh2o
      tqlc = covar_filter(t, qclw, maxwn = maxwn, nfft = nfft)(year=yr).rename('TQlwco2')
      tqlo = covar_filter(t, qolw, maxwn = maxwn, nfft = nfft)(year=yr).rename('TQlwo3')
      tqso = covar_filter(t, qosw, maxwn = maxwn, nfft = nfft)(year=yr).rename('TQswo3')
      tqlh = covar_filter(t, qhlw, maxwn = maxwn, nfft = nfft)(year=yr).rename('TQlwh2o')

      vs = [tqlc, tqlo, tqso, tqlh]
      fn = 'tq_spc_%04d.nc' % yr
   else:
      raise ValueError('Unrecognized q type %s' % q)

   path = read.path + 'gridded/kw2/'
   print('Writing %s.' % fn)
   pyg.save(path + fn, vs, version=4, compress=True)
# }}}

def do_to3_filter(yr):
# {{{
   maxwn = 30
   nfft  = 128

   times = ('16 Dec %04d' % (yr - 1), '17 Jan %04d' % (yr + 1))
   pres = (900, 0.1)

   d = read.open_gridded('eq_rad')(time=times, pres=pres)

   t   = d.Tbl
   o3  = d.O3bl
   h2o = d.H2Obl
   to = covar_filter(t, o3 , maxwn = maxwn, nfft = nfft)(year=yr)
   th = covar_filter(t, h2o, maxwn = maxwn, nfft = nfft)(year=yr)

   vs = [to, th]
   fn = 'to3_gps_%04d.nc' % yr

   path = read.path + 'gridded/kw/'
   print('Writing %s.' % fn)
   pyg.save(path + fn, vs, version=4, compress=True)
# }}}

def do_chi_filter(yr):
# {{{
   maxwn = 30
   nfft  = 128

   times = ('16 Dec %04d' % (yr - 1), '17 Jan %04d' % (yr + 1))
   pres = (900, 0.1)

   #d = read.open_gridded('kw1')(time=times, pres=pres)
   d = read.open_gridded('kw2')(time=times, pres=pres)

   o3  = d.O3bl
   h2o = d.H2Obl
   oo = var_filter(o3 , maxwn = maxwn, maxfreq=1., nfft = nfft)(year=yr)
   hh = var_filter(h2o, maxwn = maxwn, maxfreq=1., nfft = nfft)(year=yr)

   vs = [oo, hh]
   fn = 'chi_mls_%04d.nc' % yr

   path = read.path + 'gridded/kw2/'
   print('Writing %s.' % fn)
   pyg.save(path + fn, vs, version=4, compress=True)
# }}}

def profile_stats(yr, mn, dy):
# {{{
   d = read.open_rad(yr, mn, dy)
   if d is None: return
   
   d = d(lat=(-10, 10)).load()

   def ave(v): return v.mean('time')
   def frac(v):  return (v.sum('time') / float(len(v.time))).rename(v.name + 'frac')
   def count(v): return (v.sum('time')).rename(v.name + 'cnt')
   def std(v): return v.stdev('time').rename(v.name + 'std')

   bns = pyg.cldict(cdelt=5, cidelt=0.5, range=70)['clines']
   bax = pyg.NamedAxis(name='Bin', values=(bns[:-1] + bns[1:]) / 2.)

   def hist(v): 
      nbins = len(bax)
      npres = len(v.pres)
      out = np.zeros((nbins, npres), 'i')
      for p in range(npres):
         out[:, p] = np.histogram(v[:, p], bns, range=(-70, 70))[0]

      return pyg.Var((bax, v.pres), name = v.name + 'hist', values=out)

   Tbl      = ave(d.Tbl)
   Tei      = ave(d.Tei)
   Tbl_std  = std(d.Tbl)
   Tei_std  = std(d.Tei)
   Td_std   = std(d.Tbl - d.Tei).rename('dTstd')
   Td_hist  = hist(d.Tbl - d.Tei).rename('dThist')
   vs = [Tbl, Tei, Tbl_std, Tei_std, Td_std, Td_hist]

   H = 7e3
   p00 = 1000.
   kappa = 2./7.
   g0 = 9.81

   zs = -H * pyg.log(d.pres / p00)
   th = d.Tbl  * (d.pres  / p00) ** (-kappa)
   thz= th.deriv('pres', dx=zs)

   N2 = (g0 * thz / th).rename('N2')

   N2g = N2[:]
   N2_0, N2_2, N2_10 = np.percentile(N2g, [0., 2., 10.], axis = 0)

   N2m = ave(N2)
   N2_min = pyg.Var(N2m.axes, values=N2_0, name = 'N2min')
   N2_2 = pyg.Var(N2m.axes, values=N2_2, name = 'N22')
   N2_10 = pyg.Var(N2m.axes, values=N2_10, name = 'N210')
   Thz_cnt = count((1.*(N2 < 0.)).rename('DryUnstbl')) 
   vs += [N2m, N2_min, N2_2, N2_10, Thz_cnt]

   O3bl     = ave(d.O3bl)
   O3w      = ave(d.O3w)
   O3bl_std = std(d.O3bl)
   O3w_std  = std(d.O3w)
   O3d_std  = std(d.O3bl - d.O3w).rename('dO3std')
   vs += [O3bl, O3w, O3bl_std, O3w_std, O3d_std]

   H2Obl     = ave(d.H2Obl)
   H2Ow      = ave(d.H2Ow)
   H2Obl_std = std(d.H2Obl)
   H2Ow_std  = std(d.H2Ow)
   H2Od_std  = std(d.H2Obl - d.H2Ow).rename('dH2Ostd')
   vs += [H2Obl, H2Ow, H2Obl_std, H2Ow_std, H2Od_std]

   lwhr       = ave(d.lwhr)
   lwhrei     = ave(d.lwhrei*86400.).rename('lwhrei')
   lwhr_std   = std(d.lwhr)
   lwhrei_std = std(d.lwhrei*86400.).rename('lwhreistd')
   lwhrd_std  = std(d.lwhr - d.lwhrei*86400.).rename('dlwhrstd')
   vs += [lwhr, lwhrei, lwhr_std, lwhrei_std, lwhrd_std]

   swhr       = ave(d.swhr)
   swhrei     = ave(d.swhrei*86400.).rename('swhrei')
   swhr_std   = std(d.swhr)
   swhrei_std = std(d.swhrei*86400.).rename('swhreistd')
   swhrd_std  = std(d.swhr - d.swhrei*86400.).rename('dswhrstd')
   vs += [swhr, swhrei, swhr_std, swhrei_std, swhrd_std]

   lwhrtei       = ave(d.lwhrtei).rename('lwhrtei')
   lwhrteid_std  = std(d.lwhr - d.lwhrtei).rename('dlwhrtstd')
   lwhro3        = ave(d.lwhro3w).rename('lwhro3w')
   lwhro3d_std   = std(d.lwhr - d.lwhro3w).rename('dlwhro3std')
   swhro3        = ave(d.swhro3w).rename('swhro3w')
   swhro3d_std   = std(d.swhr - d.swhro3w).rename('dswhro3std')

   vs += [lwhrtei, lwhrteid_std, lwhro3, lwhro3d_std, swhro3, swhro3d_std]

   ds = pyg.Dataset(vs)
   ds.atts['profiles'] = len(d.time)

   pth = read.path + 'rad/all/%04d/' % yr
   if not os.path.exists(pth):
      print('Making %s' % pth)
      os.makedirs(pth)

   pyg.save(pth + 'mls_stats_%04d-%02d-%02d.nc' % (yr, mn, dy), ds)
# }}}

def profile_stats_by_mission(yr, mn, dy, m):
# {{{
   d = read.open_merged(yr, mn, dy, miss=[m], pres=None)
   if d is None: return
   
   d = d(lat=(-10, 10)).load()

   def ave(v):  return v.mean('time')
   def frac(v):  return (v.sum('time') / float(len(v.time))).rename(v.name + 'frac')
   def count(v): return (v.sum('time')).rename(v.name + 'cnt')
   def std(v):  return v.stdev('time').rename(v.name + 'std')

   bns = pyg.cldict(cdelt=5, cidelt=0.5, range=70)['clines']
   bax = pyg.NamedAxis(name='Bin', values=(bns[:-1] + bns[1:]) / 2.)

   def hist(v): 
      nbins = len(bax)
      npres = len(v.pres)
      out = np.zeros((nbins, npres), 'i')
      for p in range(npres):
         out[:, p] = np.histogram(v[:, p], bns, range=(-70, 70))[0]

      return pyg.Var((bax, v.pres), name = v.name + 'hist', values=out)

   Tbl      = ave(d.Tbl)
   Tei      = ave(d.Tei)
   Tbl_std  = std(d.Tbl)
   Tei_std  = std(d.Tei)
   Td_std   = std(d.Tbl - d.Tei).rename('dTstd')
   vs = [Tbl, Tei, Tbl_std, Tei_std, Td_std]

   H = 7e3
   p00 = 1000.
   kappa = 2./7.
   g0 = 9.81

   zs = -H * pyg.log(d.pres / p00)
   th = d.Tbl  * (d.pres  / p00) ** (-kappa)
   thz= th.deriv('pres', dx=zs)

   N2 = (g0 * thz / th).rename('N2')

   N2g = N2[:]
   N2_0, N2_2, N2_10 = np.percentile(N2g, [0., 2., 10.], axis = 0)

   N2m = ave(N2)
   N2_min = pyg.Var(N2m.axes, values=N2_0, name = 'N2min')
   N2_2 = pyg.Var(N2m.axes, values=N2_2, name = 'N22')
   N2_10 = pyg.Var(N2m.axes, values=N2_10, name = 'N210')
   Thz_cnt = count((1.*(N2 < 0.)).rename('DryUnstbl')) 
   vs += [N2m, N2_min, N2_2, N2_10, Thz_cnt]

   ds = pyg.Dataset(vs)
   ds.atts['profiles'] = len(d.time)

   return ds

   pth = read.path + 'rad/%s/' % m
   if not os.path.exists(pth):
      print('Making %s' % pth)
      os.makedirs(pth)

   pyg.save(pth + 'mls_%s_stats_%04d-%02d-%02d.nc' % (m, yr, mn, dy), ds)
# }}}

def profile_variance(yr, mn, dy, grid):
# {{{
   d = read.open_rad(yr, mn, dy)
   if d is None: return
   
   d = d(lat=(-30, 30), pres=(300, 5)).load()

   dg = read.open_gridded(grid)
   tm0 = dg(year = yr, month = mn, day = dy).time[0]
   dg = dg(time = (tm0, tm0+1))

   dtbs = []
   dtes = []
   for t in range(len(d.time)):
      dgi = dg.Tbl(s_lat=d.OLat[t], s_lon=d.OLon[t], s_time = dg.time[0] + d.time[t]).interpolate('pres', d.pres, inx=pyg.log(dg.pres), outx = pyg.log(d.pres))
      dTb = (d.Tbl(i_time = t) - dgi).load()
      dTe = (d.Tei(i_time = t) - dgi).load()
      dtbs.append(dTb)
      dtes.append(dTe)

   dTb = pyg.concatenate(dtbs).rename('dTbl')
   dTe = pyg.concatenate(dtes).rename('dTei')

   def ave(v): return v.mean('time')
   def std(v): return v.stdev('time').rename(v.name + 'std')

   dTbl_std  = std(dTb)
   dTei_std  = std(dTe)

   vs = [dTbl_std, dTei_std]

   ds = pyg.Dataset(vs)
   ds.atts['profiles'] = len(d.time)

   pth = read.path + 'merged/variance/%04d/' % yr
   if not os.path.exists(pth):
      print('Making %s' % pth)
      os.makedirs(pth)

   pyg.save(pth + 'profile_%s_variance_%04d-%02d-%02d.nc' % (grid, yr, mn, dy), ds)
# }}}

def eisubgrid_variance(yr, mn, dy, grid):
# {{{
   pres = (300, 5)

   dg = read.open_gridded(grid)
   tbl = dg.Tbl(year = yr, month = mn, day = dy, pres = pres)
   print(tbl.pres)

   # Interpolate onto same pressure grid
   aonp = ei.on_pres(ei.an, tbl.pres, ext=True)
   tei = aonp.t(year = yr, month = mn, day = dy, lat = (-28, 28)).replace_axes(time=tbl.time)

   tbl = tbl.interpolate('lat', tei.lat)
   dt = (tbl - tei).stdev('time', 'lat', 'lon')

   ds = pyg.Dataset([dt])

   pth = read.path + 'merged/variance/%04d/' % yr
   if not os.path.exists(pth):
      print('Making %s' % pth)
      os.makedirs(pth)

   pyg.save(pth + 'EI_%s_variance_%04d-%02d-%02d.nc' % (grid, yr, mn, dy), ds)
# }}}

if __name__ == '__main__':
   import sys
   if len(sys.argv) > 1:
      cmd = sys.argv[1]

      if cmd == 'stats':
         year   = int(sys.argv[2])
         month  = int(sys.argv[3])
         day    = int(sys.argv[4])

         for m in read.missions:
            profile_stats_by_mission(year, month, day, m)
         profile_stats(year, month, day)

      if cmd == 'variance':
         year   = int(sys.argv[2])
         month  = int(sys.argv[3])
         day    = int(sys.argv[4])

         #profile_variance(year, month, day, 'trop')
         eisubgrid_variance(year, month, day, 'trop')

      elif cmd == 'tt':
         year   = int(sys.argv[2])

         do_tt_filter(year)

      elif cmd == 'chi2':
         year   = int(sys.argv[2])

         do_chi_filter(year)

      elif cmd == 'trop':
         year   = int(sys.argv[2])

         #do_trop_tt_filter(year)
         do_ei_ttsym_filter(year)

      elif cmd == 'ttei':
         year   = int(sys.argv[2])

         do_ei_tt_filter(year)

      elif cmd == 'uwei':
         year   = int(sys.argv[2])

         do_ei_uw_filter(year)

      elif cmd == 'dsens':
         year   = int(sys.argv[2])

         do_dsens_filter(year)

      elif cmd == 'tsens':
         year   = int(sys.argv[2])

         do_tsens_filter(year)

      elif cmd == 'olr':
         year   = int(sys.argv[2])

         do_olr_filter(year)

      elif cmd == 'tqei':
         year   = int(sys.argv[2])

         do_tq_filter(year, q = 'ei')

      elif cmd == 'tqgps':
         year   = int(sys.argv[2])

         do_tq_filter(year, q = 'gps')

      elif cmd == 'tqspc':
         year   = int(sys.argv[2])

         do_tq_filter(year, q = 'spc')
      
      elif cmd == 'to3':
         year   = int(sys.argv[2])

         do_to3_filter(year)

      else:
         raise ValueError('Unknown command %s' % cmd)
