import pygeode as pyg
import numpy as np
import os, glob

#path = '/UTLS/phitch/RO/'
#if not os.path.exists(path):
   #path = '/net/modeling2' + path
path = '/local/storage/RO/'

class Theta(pyg.axis.ZAxis):
# {{{
   name='Theta'
   plotatts=pyg.axis.ZAxis.plotatts.copy()
   units='K'
   plotatts['plotname'] = 'Theta'
   formatstr = '%g'
# }}}

thdef = Theta([  302.,   304.,   306.,   308.,   310.,   313.,   315.,   317.,
                 319.,   321.,   323.,   324.,   326.,   327.,   328.,   329.,
                 331.,   332.,   333.,   335.,   336.,   337.,   338.,   340.,
                 341.,   343.,   345.,   346.,   348.,   350.,   352.,   354.,
                 355.,   357.,   359.,   361.,   363.,   365.,   367.,   368.,
                 370.,   372.,   374.,   377.,   379.,   381.,   383.,   386.,
                 388.,   391.,   393.,   396.,   399.,   402.,   405.,   409.,
                 412.,   416.,   420.,   424.,   428.,   432.,   437.,   441.,
                 445.,   450.,   454.,   459.,   464.,   468.,   473.,   478.,
                 483.,   487.,   492.,   497.,   502.,   507.,   512.,   517.,
                 522.,   527.,   532.,   537.,   542.,   548.,   553.,   558.,
                 564.,   569.,   575.,   580.,   586.,   592.,   597.,   603.,
                 609.,   615.,   621.,   627.,   634.,   640.,   646.,   653.,
                 659.,   666.,   673.,   679.,   686.,   693.,   700.,   708.,
                 715.,   723.,   730.,   738.,   746.,   754.,   763.,   771.,
                 780.,   789.,   799.,   809.,   819.,   829.,   840.,   851.,
                 863.,   876.,   890.,   904.,   920.,   936.,   955.,   975.,
                 998.,  1024.,  1055.,  1092.,  1137.,  1194.,  1268.,  1363.,
                1480.,  1629.,  1806.,  2022.,  2278.,  2572.,  3209.])

missions = ['cnofs', 'cosmic2021', 'metopa2016', 'metopb2016', 'metopa', 'metopb', 'metopc', \
            'champ2016', 'grace', 'sacc', 'tsx', 'paz', 'kompsat5', 'tdx', 'spire', 'cosmic2', \
            'geoopt', 'planetiq']
missions.sort()

def grid_profiles_pres(ds, pres):
# {{{
   if hasattr(ds, 'Pres'):
      pr = ds.Pres

      def interp(v, pres):
         if v.hasaxis('level'):
            return v.interpolate('level', pres, inx = pr)#, d_above=0., d_below=0.)
         else:
            return v

      torm = ['Pres']
   else:
      pf = ds.p
      ph = ds.ph

      def interp(v, pres):
         if v.hasaxis('layer'):
            return v.interpolate('layer', pres, inx = pf)#, d_above=0., d_below=0.)
         elif v.hasaxis('level'):
            return v.interpolate('level', pres, inx = ph)#, d_above=0., d_below=0.)
         else:
            return v

      torm = ['p', 'ph']

   return ds.remove(*torm).map(interp, pres)
# }}}
 
def grid_profiles_theta(ds, theta):
# {{{
   p00 = 1000.
   kappa = 2./7.

   Tblh = ds.Tbl.interpolate('eta', ds.etah)#, d_above=0., d_below=0.)
   Teih = ds.Tei.interpolate('eta', ds.etah)#, d_above=0., d_below=0.)

   thf_gps = ds.Tbl  * (ds.p  / p00) ** (-kappa)
   thh_gps =    Tblh * (ds.ph / p00) ** (-kappa)
   thf_ei  = ds.Tei  * (ds.p  / p00) ** (-kappa)
   thh_ei  =    Teih * (ds.ph / p00) ** (-kappa)
   
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

   return ds.remove(*torm).map(interp, theta)
# }}}

def find_lr_tp(ds):
# {{{
   R = 287.04
   cp = 1009.
   g = 9.81
   kappa = R/cp

   dz = 0.
# }}}

def open_nc(mission, year, month, days=None):
# {{{
   import glob

   ds = []

   m = mission
   if days is not None and not hasattr(days, '__len__'):
      days = [days]

   if days is None:
      fns = glob.glob(path + 'raw/%s/%s_atmPrf_%d-%02d-??.nc' % (m, m, year, month))
   else:
      fns = []
      for d in days:
         fns.extend(glob.glob(path + 'raw/%s/%s_atmPrf_%d-%02d-%02d.nc' % (m, m, year, month, d)))

   def opener(fn):
      p = pyg.open(fn)
      t = p.time.withnewvalues(np.round(p.time[:]))
      t.rtol = 1e-9
      return p.replace_axes(time=t)

   d = pyg.openall(fns, opener=opener)
   tm = pyg.StandardTime(startdate=d.time.startdate, units='days', values=d.time[:], lat=d.OLat[:], lon=d.OLon[:])
   dt = d.replace_axes(time=tm)

   return dt
# }}}

def open_nc_on_pres(year, month, pres, days=None):
# {{{
   import glob

   ds = []

   for m in missions[:]:
      if days is None:
         fns = glob.glob(path + 'raw/%s/%s_atmPrf_%d-%02d-??.nc' % (m, m, year, month))
      else:
         fns = []
         for d in days:
            fns.extend(glob.glob(path + 'raw/%s/%s_atmPrf_%d-%02d-%02d.nc' % (m, m, year, month, d)))

      if len(fns) == 0:
         continue

      d = pyg.openall(fns)
      tm = pyg.StandardTime(startdate=d.time.startdate, units='days', values=d.time[:], lat=d.OLat[:], lon=d.OLon[:])
      dt = d.replace_axes(time=tm)
      dp = grid_profiles_pres(dt, pres) 
         
      ds.append(dp)

   if len(ds) == 0: return None

   dm = pyg.concatenate(ds).sorted('time')

   return dm
# }}}

def open_merged(year, month, days=None, dset='mls', miss=None, pres=False):
# {{{
   import glob

   ds = []

   if days is not None and not hasattr(days, '__len__'):
      days = [days]

   fns = []
   if miss is None:
      miss = missions
   for m in miss:
      pth = path + 'merged/%s/%04d-%02d/' % (dset, year, month)
      if days is None:
         fns.extend(glob.glob(pth + '%s_*merged_%04d-%02d-??.nc' % (m, year, month)))
      else:
         for d in days:
            fns.extend(glob.glob(pth + '%s_*merged_%04d-%02d-%02d.nc' % (m, year, month, d)))

   if len(fns) == 0: return None

   d = pyg.openall(fns, sorted=False).sorted('time')
   tm = pyg.StandardTime(startdate=d.time.startdate, units='days', values=d.time[:], lat=d.OLat[:], lon=d.OLon[:])
   dt = d.replace_axes(time=tm)
   
   #dt = dt.rename_axes(eta='layer', etah='level')

   if pres is False: 
      return dt

   if pres is None:
      pres = pyg.Pres(dt.layer[:] * 1013.25)

   return grid_profiles_pres(dt, pres) 
# }}}

def open_rad(year, month, days=None, rset='pyr', dset='mls', pres=None, theta=None):
#U {{{
   import glob

   ds = []

   if days is not None and not hasattr(days, '__len__'):
      days = [days]

   fns = []
   pth = path + 'merged/%s/%04d-%02d/' % (dset, year, month)
   for m in missions:
      if days is None:
         fns.extend(glob.glob(pth + '%s_*merged_%04d-%02d-??.nc' % (m, year, month)))
      else:
         for d in days:
            fns.extend(glob.glob(pth + '%s_*merged_%04d-%02d-%02d.nc' % (m, year, month, d)))

   rfns = []
   rpth = path + 'rad/%s/%04d/' % (rset, year)
   if days is None:
      rfns.extend(glob.glob(rpth + 'rrtm_%04d-%02d-??.nc' % (year, month)))
   else:
      for d in days:
         rfns.extend(glob.glob(rpth + 'rrtm_%04d-%02d-%02d.nc' % (year, month, d)))

   if len(fns) == 0 or len(rfns) == 0: return None

   #rvlist = ['cosz', 'lwhr', 'swhr', 'dflxlw', 'uflxlw', 'dflxsw', 'uflxsw']
   #dr = pyg.openall(fns, varlist=rvlist, sorted=False).sorted('time')

   dm = pyg.openall(fns, sorted=False).sorted('time')
   dr = pyg.openall(rfns, sorted=False).sorted('time')

   tm = pyg.StandardTime(startdate=dm.time.startdate, units='days', values=dm.time[:], lat=dm.OLat[:], lon=dm.OLon[:])

   if len(tm) != len(dr.time): 
      print('Warning; %d %d %s has mismatch between merged and radiative dataset.' % (year, month, str(days)))
      return None

   dt = dm.replace_axes(time=tm) + dr.replace_axes(time=tm)

   if pres == False:
      return dt

   if theta is not None:
      return grid_profiles_theta(dt, theta)

   if pres is None:
      pres = pyg.Pres(dt.eta[:] * 1013.25)

   dt = dt.rename_axes(eta='layer', etah='level')

   return grid_profiles_pres(dt, pres) 
# }}}

def open_gridded(dset, vers=0):
# {{{
   pths = dict(kw1 = 'eq_rad',
               kw2 = 'eq_rad_5_p5',
               trop = 'trop_D10_T12h',
               trop2 = 'eq_rad_D5_T6h2',
               dsens = 'eq_grid',
               dsens6h = 'eq_grid',
               dsens1d = 'eq_grid',
               tsens = 'eq_grid',
               zm  = 'zm')

   fns = dict(kw1 = 'ro_1deg_D10deg_T1d_kw_gauss_gridded_hvr_????-??',
              kw2 = 'ro_1deg_D5deg_T0p5d_kw_gauss_gridded_hvr_????-??',
              trop = 'ro_1deg_D10deg_T12h_gridded_????-??',
              trop2 = 'ro_1deg_D5deg_T6h_kw_gauss_gridded_hvr_????-??',
              dsens = 'ro_kw_gauss_gridded_Dsens_????-??',
              dsens6h = 'ro_kw_gauss_gridded_Dsens_T6h_????-??',
              dsens1d = 'ro_kw_gauss_gridded_Dsens_T1d_????-??',
              tsens = 'ro_kw_gauss_gridded_Tsens_????-??',
              zm  = 'ro_1deg_D5deg_T0p5d_zm_gauss_gridded_????-??')

   ver = dict(kw1 = ['-v3', '-v2'],
              kw2 = [''],
              trop = [''],
              trop2 = [''],
              dsens = [''],
              dsens6h = [''],
              dsens1d = [''],
              tsens = [''],
              zm = [''])

   root = path + 'gridded/%s/' % pths[dset]
   fns = [root + fns[dset] + ver[dset][vers] + '.nc']
   return pyg.open_multi(fns, pattern='$Y-$m')
# }}}

def open_filtered(dset, var, vers='-v2'):
# {{{
   root = path + 'gridded/%s/' % dset
   fns = [root + '%s_????%s.nc' % (var, vers)]
   #fns = [root + '%s_200[789]%s.nc' % (var, vers)]
   #fns  = [root + '%s_201[01345]%s.nc' % (var, vers)]
   return pyg.openall(fns)
# }}}

def open_stats(dset, miss=None):
# {{{
   def opener(fn):
      import re
      patt = '(?P<year>[0-9]{4})-(?P<month>[0-9]{2})-(?P<day>[0-9]{2})'
      rx = re.compile(patt)
      d = rx.search(fn)
      assert d is not None, "can't use the pattern on the filenames?"
      d = d.groupdict()
      d = dict([k,[int(v)]] for k,v in d.items() if v is not None)

      tm = pyg.StandardTime(units = 'days', startdate = dict(year = 2000, month = 1, day = 1), **d)
      ds = pyg.open(fn).extend(0, tm)
      ds += pyg.Var((tm,), values = [ds.atts['profiles']], name = 'profiles')
      return ds

   if miss is None:
      root = path + 'rad/all/'
      fns = [root + dset + '_stats_????-??-??.nc']
   else:
      root = path + 'rad/%s/' % miss
      fns = [root + '_'.join([dset, miss, 'stats_????-??-??.nc'])]

   return pyg.open_multi(fns, pattern='$Y-$m-$d', opener=opener)
# }}}

def open_var():
# {{{
   def opener(fn):
      import re
      patt = '(?P<year>[0-9]{4})-(?P<month>[0-9]{2})-(?P<day>[0-9]{2})'
      rx = re.compile(patt)
      d = rx.search(fn)
      assert d is not None, "can't use the pattern on the filenames?"
      d = d.groupdict()
      d = dict([k,[int(v)]] for k,v in d.items() if v is not None)

      tm = pyg.StandardTime(units = 'days', startdate = dict(year = 2000, month = 1, day = 1), **d)
      ds = pyg.open(fn).extend(0, tm)
      ds += pyg.Var((tm,), values = [ds.atts['profiles']], name = 'profiles')
      return ds

   root = path + 'merged/variance/'
   fns = [root + '%04d/*_variance_????-??-??.nc' % y for y in range(2007, 2017)]

   return pyg.open_multi(fns, pattern='$Y-$m-$d', opener=opener)
# }}}

def open_summary(mission):
# {{{
    spath = path + 'raw/%s/%s_summary*.nc' % (mission, mission)
    return pyg.open_multi(spath, pattern='$Y-$m')
# }}}

def write_summary(raw=False, merged=False, rrtm=False, reprocess=False):
# {{{
   t0 = pyg.standardtimen('2000-01-01', 1)

   if raw:
      for m in missions:
         print('============== %s ==============' % m)
         for yr in np.arange(2006, 2019):
            for mn in range(1, 13):
               summary_file = path + 'raw/%s/%s_summary_%d-%02d.nc' % (m, m, yr, mn)
               if not reprocess and os.path.exists(summary_file):
                   print('Summary file exists. Skipping.')
                   continue

               days = int(t0.days_in_month(yr, mn))
               tax = pyg.standardtimen('%04d-%02d-01' % (yr, mn), days)
               ntot = np.zeros(days, 'i')
               ntrop = np.zeros(days, 'i')
               for dy in range(days):
                  try:
                     ds = open_nc(m, yr, mn, dy + 1)
                     ntot[dy] = ds.time.shape[0]
                     ntrop[dy] = ds.time(lat=(-10, 10)).shape[0]
                  except AssertionError as e:
                     pass
      
               pr = pyg.Var((tax,), values = ntot, name = f'profiles')
               prt = pyg.Var((tax,), values = ntrop, name = f'tropical_profiles')
               pyg.save(summary_file, [pr, prt])
   if merged:
      for m in missions:
         print('============== %s ==============' % m)
         for yr in np.arange(2001, 2017):
            for mn in range(1, 13):
               d = open_merged(yr, mn, dset='mls', miss=[m])
               if d is not None:
                  print('%04d-%02d: %d profiles.' % (yr, mn, len(d.time)))

      print('============== all ==============')
      for yr in np.arange(2001, 2017):
         for mn in range(1, 13):
            d = open_merged(yr, mn, dset='mls')
            if d is not None:
               print('%04d-%02d: %d profiles.' % (yr, mn, len(d.time)))
   if rrtm:
      print('============== all ==============')
      for yr in np.arange(2012, 2017):
         for mn in range(1, 13):
            d = open_rad(yr, mn, rset='mls', dset='mls')
            if d is not None:
               print('%04d-%02d: %d profiles.' % (yr, mn, len(d.time)))
# }}}

def summary(raw=False, merged=False, rrtm=False, reprocess=False):
# {{{
   t0 = pyg.standardtimen('2000-01-01', 1)

   if raw:
      for m in missions:
         patt = path + f'raw/{m}/{m}_atmPrf_*.nc'
         fns = glob.glob(patt)

         if fns == None or len(fns) == 0: 
             print(f'{m:<13}: No data.')
             continue

         dts = [f[-13:-3] for f in fns]
         dts.sort()

         ndates = len(dts)
         yr1 = dts[0][:4]
         v1 = t0.str_as_val('', dts[0])
         doy1 = int(v1 - t0.str_as_val('', f'{yr1}-01-01') + 1)

         yr2 = dts[-1][:4]
         v2 = t0.str_as_val('', dts[-1])
         doy2 = int(v2 - t0.str_as_val('', f'{yr2}-01-01') + 1)
         ndays = int(v2 - v1)

         print(f'{m:<13}: {yr1}.{doy1:03d} - {yr2}.{doy2:03d}  '
               f'{dts[0]} - {dts[-1]}; {ndays} of which {ndays - ndates} are missing.')

# }}}

if __name__ == '__main__':
   import sys
   if len(sys.argv) > 1:
      cmd = sys.argv[1]

      if cmd == 'summarizeraw':
          write_summary(raw = True)
