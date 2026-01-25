#from . import read
import read
import numpy as np, pylab as pyl
import pygeode as pyg

from ro import ro
from ECMWF import ERAI as ei

def eq_grid(yr, mn):
# {{{
   #pres = ei.an_on_p.pres
   pres = pyg.Pres(np.exp(np.linspace(np.log(450), np.log(2.5), 200)))

   if mn == 1: 
      yrm = yr - 1
      mnm = 12
   else:
      yrm = yr
      mnm = mn - 1

   if mn == 12:
      yrp = yr + 1
      mnp = 1
   else:
      yrp = yr
      mnp = mn + 1

   dmm = np.arange(-2, 1) + ei.an_on_p.time.days_in_month(yrm, mnm)
   dmp = np.arange(1, 4)

   dmm = read.open_nc_on_pres(yrm, mnm, pres, days = dmm)
   dm  = read.open_nc_on_pres(yr,  mn,  pres)
   dmp = read.open_nc_on_pres(yrp, mnp, pres, days = dmp)

   dm = pyg.concatenate([d for d in [dmm, dm, dmp] if d is not None])

   d0 = '1 %s %04d' % (pyg.timeaxis.months[mn], yr)
   d1 = '1 %s %04d' % (pyg.timeaxis.months[mnp], yrp)

   #tm = pyg.StandardTime(startdate=dm.time.startdate, units='days', values=dm.time[:], lat=dm.OLat[:], lon=dm.OLon[:])
   #pres = pyg.Pres(10**np.linspace(2.6, 0.5, 200))

   #dmt = dm.replace_axes(time=tm)
   dm = dm(lat = (-10, 10)).load()

   dgt, dgtc = ro.grid_kw_gauss(dm.Temp, d0, d1, 1, 0.25, 10., 1.)
   dgz, dgzc = ro.grid_kw_gauss(dm.Alt,  d0, d1, 1, 0.25, 10., 1.)
   
   ds = pyg.Dataset([dgt, dgz, dgtc.rename('count')])

   fn = '/data2/RO/gridded/eq/ro_1deg_D10deg_T1d_kw_gauss_gridded_hvr_%04d-%02d.nc' % (yr, mn)
   print('Writing %s.' % fn)
   pyg.save(fn, ds)
# }}}

def eq_merged_grid(yr, mn):
# {{{
   tax = pyg.standardtimen('1 Jan 2002', 1)

   if mn == 1: 
      yrm = yr - 1
      mnm = 12
   else:
      yrm = yr
      mnm = mn - 1

   if mn == 12:
      yrp = yr + 1
      mnp = 1
   else:
      yrp = yr
      mnp = mn + 1

   dmm = np.arange(-9, 1) + tax.days_in_month(yrm, mnm)
   dmp = np.arange(1, 10)

   dmm = read.open_rad(yrm, mnm, days = dmm)
   dm  = read.open_rad(yr,  mn)
   dmp = read.open_rad(yrp, mnp, days = dmp)

   dm = pyg.concatenate([d for d in [dmm, dm, dmp] if d is not None])

   d0 = '1 %s %04d' % (pyg.timeaxis.months[mn], yr)
   d1 = '1 %s %04d' % (pyg.timeaxis.months[mnp], yrp)

   dm = dm(lat = (-10, 10)).load()

   def scl(v): return (86400.*v).rename(v.name)

   vars = [dm.Tbl, dm.Tei, \
         scl(dm.lwhrei), scl(dm.swhrei), dm.lwhr, dm.swhr, \
         dm.lwhrco2, dm.lwhro3, dm.lwhrh2o, \
         dm.swhrco2, dm.swhro3, dm.swhrh2o, \
         #dm.lwhrtei, dm.lwhro3w, dm.swhro3w, \
         dm.CO2w, dm.O3w, dm.H2Ow, dm.O3bl, dm.H2Obl]

   vs, cs = ro.grid_kw_gauss(vars,  d0, d1, 1, 0.25, 5., 0.25)
   
   ds = pyg.Dataset(vs + cs)

   fn = read.path + 'gridded/eq_rad_D5_T6h2/ro_1deg_D5deg_T6h_kw_gauss_gridded_hvr_%04d-%02d.nc' % (yr, mn)
   print('Writing %s.' % fn)
   pyg.save(fn, ds, version=4, compress=True)
   return ds
# }}}

def eq_merged_grid_e5(yr, mn):
# {{{
   tax = pyg.standardtimen('1 Jan 2002', 1)

   if mn == 1: 
      yrm = yr - 1
      mnm = 12
   else:
      yrm = yr
      mnm = mn - 1

   if mn == 12:
      yrp = yr + 1
      mnp = 1
   else:
      yrp = yr
      mnp = mn + 1

   dmm = np.arange(-9, 1) + tax.days_in_month(yrm, mnm)
   dmp = np.arange(1, 10)

   rset = 'mls-e5'
   dset = 'mls-e5'

   dmm = read.open_rad(yrm, mnm, days = dmm, rset=rset, dset=dset)
   dm  = read.open_rad(yr,  mn,              rset=rset, dset=dset)
   dmp = read.open_rad(yrp, mnp, days = dmp, rset=rset, dset=dset)

   dm = pyg.concatenate([d for d in [dmm, dm, dmp] if d is not None])

   d0 = '1 %s %04d' % (pyg.timeaxis.months[mn], yr)
   d1 = '1 %s %04d' % (pyg.timeaxis.months[mnp], yrp)

   dm = dm(lat = (-10, 10)).load()

   def scl(v): return (86400.*v).rename(v.name)

   vars = [dm.Tbl, dm.Te5, \
         dm.lwhr, dm.swhr, \
         dm.lwhrco2, dm.lwhro3, dm.lwhrh2o, \
         dm.swhrco2, dm.swhro3, dm.swhrh2o, \
         dm.lwhre5h2o, dm.swhre5h2o, \
         #scl(dm.lwhrei), scl(dm.swhrei), 
         dm.CO2w, dm.O3w, dm.H2Ow, dm.O3bl, dm.H2Obl]

   vs, cs = ro.grid_kw_gauss(vars,  d0, d1, 1, 0.25, 5., 0.25)
   
   ds = pyg.Dataset(vs + cs)

   fn = read.path + 'gridded/eq_rad_e5_D5_T6h/ro_1deg_D5deg_T6h_kw_gauss_gridded_%04d-%02d.nc' % (yr, mn)
   print('Writing %s.' % fn)
   pyg.save(fn, ds, version=4, compress=True)
   return ds
# }}}

def trop_grid(yr, mn):
# {{{
   if mn == 1: 
      yrm = yr - 1
      mnm = 12
   else:
      yrm = yr
      mnm = mn - 1

   if mn == 12:
      yrp = yr + 1
      mnp = 1
   else:
      yrp = yr
      mnp = mn + 1

   dmm = np.arange(-9, 1) + ei.an_on_p.time.days_in_month(yrm, mnm)
   dmp = np.arange(1, 10)

   dmm = read.open_rad(yrm, mnm, days = dmm, pres=ei.an_on_p.pres)
   dm  = read.open_rad(yr,  mn,              pres=ei.an_on_p.pres)
   dmp = read.open_rad(yrp, mnp, days = dmp, pres=ei.an_on_p.pres)

   dm = pyg.concatenate([d for d in [dmm, dm, dmp] if d is not None])

   d0 = '1 %s %04d' % (pyg.timeaxis.months[mn], yr)
   d1 = '1 %s %04d' % (pyg.timeaxis.months[mnp], yrp)

   dm = dm(lat = (-30, 30))

   def scl(v): return (86400.*v).rename(v.name)

   vars = [v.load() for v in [dm.Tbl, dm.Tei]]

   vs, ds, ts, ns = ro.grid_trop_gauss(vars, d0, d1, 1, 2, 0.25, 10., 1.0)
   #vs, ds, ts, ns = ro.grid_trop_gauss(vars, d0, d1, 1, 4, 0.25, 5., 0.25)
   
   ds = pyg.Dataset(vs + ds + ts + ns)

   #fn = read.path + 'gridded/trop_D5_T6h/ro_1deg_D5deg_T6h_gridded_%04d-%02d.nc' % (yr, mn)
   fn = read.path + 'gridded/trop_D10_T24h/ro_1deg_D10deg_T24h_gridded_%04d-%02d.nc' % (yr, mn)
   print('Writing %s.' % fn)
   pyg.save(fn, ds, version=4, compress=True)
   return ds
# }}}

def eq_grid_sens(yr, mn):
# {{{

   if mn == 1: 
      yrm = yr - 1
      mnm = 12
   else:
      yrm = yr
      mnm = mn - 1

   if mn == 12:
      yrp = yr + 1
      mnp = 1
   else:
      yrp = yr
      mnp = mn + 1

   dmm = np.arange(-5, 1) + ei.an_on_p.time.days_in_month(yrm, mnm)
   dmp = np.arange(1, 6)

   dmm = read.open_rad(yrm, mnm, days = dmm, pres=ei.an_on_p.pres)
   dm  = read.open_rad(yr,  mn,              pres=ei.an_on_p.pres)
   dmp = read.open_rad(yrp, mnp, days = dmp, pres=ei.an_on_p.pres)

   dm = pyg.concatenate([d for d in [dmm, dm, dmp] if d is not None])

   d0 = '1 %s %04d' % (pyg.timeaxis.months[mn], yr)
   d1 = '1 %s %04d' % (pyg.timeaxis.months[mnp], yrp)

   dm = dm(lat = (-10, 10)).load()

   def scl(v): return (86400.*v).rename(v.name)

   vars = [dm.Tbl, dm.Tei]

   #Ds = np.array([30., 20., 10.,  5.,  2.,  1.])
   #Ts = np.array([0.5, 0.5, 0.5, 0.5, 0.5, 0.5])

   #Ds = np.array([1., 1.,  1.])
   #Ts = np.array([1, 0.5, 0.25])

   #Ds = np.array([30., 10.,  5.])
   #Ts = np.array([0.25, 0.25, 0.25])

   Ds = np.array([30., 10.,  5.])
   Ts = np.array([1., 1., 1.])

   vs, cs = ro.grid_kw_gauss_sens(vars,  d0, d1, 1, 0.25, Ds, Ts)
   
   ds = pyg.Dataset(vs + cs)

   fn = read.path + 'gridded/eq_grid/ro_kw_gauss_gridded_Dsens_T1d_%04d-%02d.nc' % (yr, mn)
   print('Writing %s.' % fn)
   pyg.save(fn, ds, version=4, compress=True)
   return ds
# }}}

def zm_merged_grid(yr, mn):
# {{{
   if mn == 1: 
      yrm = yr - 1
      mnm = 12
   else:
      yrm = yr
      mnm = mn - 1

   if mn == 12:
      yrp = yr + 1
      mnp = 1
   else:
      yrp = yr
      mnp = mn + 1

   dmm = np.arange(-2, 1) + ei.an_on_p.time.days_in_month(yrm, mnm)
   dmp = np.arange(1, 4)

   dmm = read.open_rad(yrm, mnm, days = dmm)
   dm  = read.open_rad(yr,  mn)
   dmp = read.open_rad(yrp, mnp, days = dmp)

   dm = pyg.concatenate([d for d in [dmm, dm, dmp] if d is not None])

   d0 = '1 %s %04d' % (pyg.timeaxis.months[mn], yr)
   d1 = '1 %s %04d' % (pyg.timeaxis.months[mnp], yrp)

   dm = dm.load()

   def scl(v): return (86400.*v).rename(v.name)

   vars = [dm.Tbl, dm.Tei, \
         scl(dm.lwhrei), scl(dm.swhrei), dm.lwhr, dm.swhr, \
         dm.lwhrco2, dm.lwhro3, dm.lwhrh2o, \
         dm.swhrco2, dm.swhro3, dm.swhrh2o, \
         dm.lwhrtei, dm.lwhro3w, dm.swhro3w, \
         dm.CO2w, dm.O3w, dm.H2Ow, dm.O3bl, dm.H2Obl]

   vs, cs = ro.grid_zm_gauss(vars,  d0, d1, 2, 0.25, 5., 1.0)
   
   ds = pyg.Dataset(vs + cs)

   fn = read.path + 'gridded/zm/ro_1deg_D5deg_T0p5d_zm_gauss_gridded_%04d-%02d.nc' % (yr, mn)
   print('Writing %s.' % fn)
   pyg.save(fn, ds, version=4, compress=True)
   return ds
# }}}

def grid_profiles_theta(ds, theta):
# {{{
   p00   = 1000.
   kappa = 287.04 / 1009.
   g0    = 9.81

   Tblh = ds.Tbl.interpolate('eta', ds.etah, interp_type='linear')
   Teih = ds.Tei.interpolate('eta', ds.etah, interp_type='linear')

   thf_gps = ds.Tbl  * (ds.p  / p00) ** (-kappa)
   thh_gps =    Tblh * (ds.ph / p00) ** (-kappa)
   thf_ei  = ds.Tei  * (ds.p  / p00) ** (-kappa)
   thh_ei  =    Teih * (ds.ph / p00) ** (-kappa)
   
   #thck_gps = -ds.ph.diff(axis='etah') / thh_gps.diff(axis='etah') / g0
   #thck_ei  = -ds.ph.diff(axis='etah') / thh_ei.diff(axis='etah') / g0

   thck_gps = -ds.p.deriv('eta', dx=thf_gps) / g0
   thck_ei  = -ds.p.deriv('eta', dx=thf_ei)  / g0

   thck_gps = thck_gps * ((thck_gps > 0.).unfill(False))
   thck_ei  = thck_ei  * ((thck_ei  > 0.).unfill(False))

   thck_gps = thck_gps.replace_axes(eta = ds.eta).rename('Thkbl')
   thck_ei  = thck_ei. replace_axes(eta = ds.eta).rename('Thkei')

   dst = pyg.Dataset([thck_gps, thck_ei])
   
   torm = ['p', 'ph']

   eivars = ['Tei', 'lwhrei', 'swhrei']

   def interp(v, theta):
      if v.name in eivars:
         thf = thf_ei
         thh = thh_ei
      else:
         thf = thf_gps
         thh = thh_gps

      if v.hasaxis('etah'):
         return v.interpolate('etah', theta, inx = thh, omit_nonmonotonic=True, interp_type='linear')# d_above=0., d_below=0.)
      elif v.hasaxis('eta'):
         return v.interpolate('eta',  theta, inx = thf, omit_nonmonotonic=True, interp_type='linear')# d_above=0., d_below=0.)
      else:
         return v

      #if v.hasaxis('etah'):
         #return v.interpolate('etah', theta, inx = thh, omit_nonmonotonic=True)# d_above=0., d_below=0.)
      #elif v.hasaxis('eta'):
         #return v.interpolate('eta',  theta, inx = thf, omit_nonmonotonic=True)# d_above=0., d_below=0.)

   return ds.remove(*torm).map(interp, theta) + dst.map(interp, theta)
# }}}

def zm_merged_theta_grid(yr, mn):
# {{{
   if mn == 1: 
      yrm = yr - 1
      mnm = 12
   else:
      yrm = yr
      mnm = mn - 1

   if mn == 12:
      yrp = yr + 1
      mnp = 1
   else:
      yrp = yr
      mnp = mn + 1

   dmm = np.arange(-2, 1) + ei.an_on_p.time.days_in_month(yrm, mnm)
   dmp = np.arange(1, 4)

   dmm = read.open_rad(yrm, mnm, days = dmm, pres=False)
   dm  = read.open_rad(yr,  mn,              pres=False)
   dmp = read.open_rad(yrp, mnp, days = dmp, pres=False)

   dm = pyg.concatenate([d for d in [dmm, dm, dmp] if d is not None])

   d0 = '1 %s %04d' % (pyg.timeaxis.months[mn], yr)
   d1 = '1 %s %04d' % (pyg.timeaxis.months[mnp], yrp)

   dm = grid_profiles_theta(dm, read.thdef).load()

   def scl(v): return (86400.*v).rename(v.name)

   vars = [dm.Tbl, dm.Tei, dm.Thkbl, dm.Thkei, \
         scl(dm.lwhrei), scl(dm.swhrei), dm.lwhr, dm.swhr]#, \
        #dm.lwhrco2, dm.lwhro3, dm.lwhrh2o, \
        #dm.swhrco2, dm.swhro3, dm.swhrh2o, \
        #dm.lwhrtei, dm.lwhro3w, dm.swhro3w, \
        #dm.CO2w, dm.O3w, dm.H2Ow, dm.O3bl, dm.H2Obl]

   # Add thickness weighted variables

   thTbl  = (dm.Tbl * dm.Thkbl).rename('thTbl')
   thTei  = (dm.Tei * dm.Thkei).rename('thTei')
   thlwei = (scl(dm.lwhrei) * dm.Thkei).rename('thlwei')
   thswei = (scl(dm.swhrei) * dm.Thkei).rename('thswei')
   thlwhr = (dm.lwhr * dm.Thkbl).rename('thlwhr')
   thswhr = (dm.swhr * dm.Thkbl).rename('thswhr')

   vars += [thTbl, thTei, thlwei, thswei, thlwhr, thswhr]

   vs, cs = ro.grid_zm_gauss(vars,  d0, d1, 2, 1, 5., 1.0)
   
   ds = pyg.Dataset(vs + cs)

   fn = read.path + 'gridded/zm_th/ro_D5deg_T1p0d_zm_gauss_gridded_%04d-%02d.nc' % (yr, mn)
   print('Writing %s.' % fn)
   pyg.save(fn, ds, version=4, compress=True)
   return ds
# }}}

if __name__ == '__main__':
   import sys
   if len(sys.argv) > 1:
      cmd = sys.argv[1]

      if cmd == 'eq':
         year   = int(sys.argv[2])
         month  = int(sys.argv[3])

         eq_grid(year, month)

      elif cmd == 'eqrad':
         year   = int(sys.argv[2])
         month  = int(sys.argv[3])

         eq_merged_grid_e5(year, month)

      elif cmd == 'trop':
         year   = int(sys.argv[2])
         month  = int(sys.argv[3])

         trop_grid(year, month)

      elif cmd == 'sens':
         year   = int(sys.argv[2])
         month  = int(sys.argv[3])

         eq_grid_sens(year, month)

      elif cmd == 'zm':
         year   = int(sys.argv[2])
         month  = int(sys.argv[3])

         zm_merged_grid(year, month)

      elif cmd == 'zmth':
         year   = int(sys.argv[2])
         month  = int(sys.argv[3])

         zm_merged_theta_grid(year, month)

      else:
         raise ValueError('Unknown command %s' % cmd)
