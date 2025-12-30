import pygeode as pyg
import numpy as np

all_dates = set()
hashed_dates = {}

rpath = '/local/storage/RO/'

def open_profile(fn, N=3100):
# {{{
   dt = dict(MSL_alt = pyg.Height)
   nm = dict(MSL_alt = 'alt')
   vl = ['Lat', 'Lon', 'Pres', 'Temp']
   attl = ['occulting_sat_id', 'irs', 'lat', 'lon', 'rgeoid', 'stop_time', \
         'trhcp', 'trhwmo', 'trhwmo2', 'trtcp', 'trtwmo', 'trtwmo2']

   ds = pyg.open(fn, dimtypes = dt, namemap = nm, varlist = vl, format = 'netcdf')

   if int(ds.atts['bad']) == 1:
      #print 'Rejecting %s; %s' % (fn, ds.atts['errstr'])
      return None

   ncom = len(ds.alt) - len(np.unique(ds.alt))
   if ncom > 0:
      print('Warning: %s contains %d non-unique altitudes. Omitting duplicate gridpoints.' % (fn, ncom))
      icom = np.where(np.diff(ds.alt[:]) == 0.)[0]
      alts = ds.alt[:].copy()
      alts[icom] = -100. - (icom - icom[0])
      dt['MSL_alt'] = pyg.Height(name='alt', values=alts, atts=ds.alt.atts)
      rng = (ds.alt[0], ds.alt[-1])
      ds = pyg.open(fn, dimtypes = dt, namemap = nm, varlist = vl, format = 'netcdf')
      ds = ds(alt = rng)

   ds = ds.sorted('alt')
   ds.alt.units = ds.alt.atts['units']
   Alt = (ds.alt + 0.).as_type('f').rename('Alt')

   nlevs = len(ds.alt)

   if nlevs >= N:
      print('Warning: truncating %s from %d to %d levels.' % (fn, nlevs, N))
      level = pyg.NamedAxis(name='level', values=np.arange(N))
      def reshape(v): 
         V = v(i_alt = (0, N)).replace_axes(alt = level)
         V = V.extend(0, t)
         V.atts.update(v.atts)
         return V.rename(v.name)

   else:
      lvl0 = pyg.NamedAxis(name='level', values=np.arange(nlevs))
      lvl1 = pyg.NamedAxis(name='level', values=np.arange(nlevs, N))

      def reshape(v): 
         V = pyg.concat.concat(v.replace_axes(alt=lvl0), (np.nan * lvl1).as_type(v.dtype))
         V = V.extend(0, t)
         V.atts.update(v.atts)
         return V.rename(v.name)

   def encode_date(year, month, day, hour, minute, second):
      return '%04d%02d%02d%02d%02d%02d' % (year[0], month[0], day[0], hour[0], minute[0], second[0])

   date = {k: np.array([int(ds.atts.get(k, 0))]) for k in ['year', 'month', 'day', 'hour', 'minute', 'second']}
   t = pyg.StandardTime(units='days', starttime='1 Jan 2000', **date)

   dstr = encode_date(**date)
   hsh = 0
   while dstr in all_dates:
      date['second'] += 1
      date = t.val_as_date(t.date_as_val(date))
      dstr = encode_date(**date)
      hsh += 1

   if hsh > 0: 
      print('Offsetting %s by %d second(s) to %s' % (fn, hsh, dstr))
      hashed_dates[fn] = hsh

   all_dates.add(dstr)

   t = pyg.StandardTime(units='days', starttime='1 Jan 2000', **date)

   def make_var(a):
      name = a
      vals = ds.atts[a]
      if a == 'occulting_sat_id':
         name = 'GPSsatID'
      elif a == 'lat':
         name = 'OLat'
      elif a == 'lon':
         name = 'OLon'
      elif a == 'stop_time':
         name = 'ODuration'
         vals = vals - ds.atts['start_time']

      return pyg.Var((t,), name = name, values=[vals])

   vs = pyg.Dataset([make_var(a) for a in attl])

   ds.atts = {}

   return ds.map(reshape) + pyg.Dataset([reshape(Alt)]) + vs 
# }}} 

def write_date_hash(mission, year, doy, dates, hashes):
# {{{
   import os
   import datetime

   dhpath = rpath + 'raw/date_hash/'
   dhfn = '%04d-%03d_%s.txt' % (year, doy, mission)

   today = datetime.date.today().isoformat()

   if not os.path.exists(dhpath): os.makedirs(dhpath)

   dates = list(dates)
   dates.sort()

   with open(dhpath + dhfn, 'a') as f:
      f.write('%s processed %s\n' % (mission, today))

      for d in dates:
         f.write(d + '\n')

      for k, v in hashes.items():
         f.write('%s delayed %d second(s).\n' % (k, v))

      f.write('---------------------\n')
# }}}

def read_date_hash(year, doy, mission=None):
# {{{
   import os
   import datetime
   import glob

   dhpath = rpath + 'raw/date_hash/'
   fngl = '%04d-%03d_*.txt' % (year, doy)

   fns = glob.glob(dhpath + fngl)

   if len(fns) == 0:
      print('No previous missions found.')
      return set()

   missions = [os.path.basename(f)[9:-4] for f in fns]
   if mission is not None and mission in missions:
      raise ValueError("%s already processed for %04d-%03d" % (mission, year, doy))

   print('%d previously processed mission(s): ' % len(missions) + ', '.join(missions))

   today = datetime.date.today().isoformat()

   dates = set()

   for m in missions:
      dhfn = '%04d-%03d_%s.txt' % (year, doy, m)

      with open(dhpath + dhfn, 'r') as f:
         for line in f:
            if line[:4] == str(year): 
               line = line.strip()
               if line in dates:
                  print('Hash fail:', line)
               dates.add(line)

   return dates
# }}}

def process_date(mission, year, doy, Nz = 3100, files = rpath + 'dl/{mission}/atmPrf/{year:04d}.{doy:03d}/*_nc'):
# {{{
   import sys

   global all_dates
   global hashed_dates

   all_dates = read_date_hash(year, doy, mission)
   old_dates = set()
   old_dates.update(all_dates)
   hashed_dates = {}

   print('# of old dates: %d' % len(old_dates))

   fns = files.format(mission=mission, year=year, doy=doy)

   log_fn = rpath + 'raw/{mission}/logs/{mission}_atmPrf_{year:04d}-{doy:03d}_log.txt'.format(mission=mission, year=year, doy=doy)

   with open(log_fn, 'w') as f:
      sys.stdout = f
      try:
         ds = pyg.openall(fns, format='netcdf', opener=open_profile, N = Nz)
      except AssertionError as e:
         print(fns)
         sys.stdout = sys.__stdout__
         print(e, 'Aborting.')
         return
      sys.stdout = sys.__stdout__

   new_dates = all_dates - old_dates
   print('# of new dates: %d' % len(new_dates))

   month = ds.time.month[0]
   day   = ds.time.day[0]

   out_fn = rpath + 'raw/{mission}/{mission}_atmPrf_{year:04d}-{month:02d}-{day:02d}.nc'.format(mission=mission, year=year, month=month, day=day)
   print('Writing %s.' % out_fn)
   pyg.save(out_fn, ds, version = 4, compress = True)

   write_date_hash(mission, year, doy, new_dates, hashed_dates)
# }}}

def open_nc(year, month):
# {{{
   ds = []
   missions = ['cnofs', 'cosmic2013', 'metopa2016', 'grace', 'sacc']#, 'tsx']

   for m in missions[1:]:
      try:
         #fn = '%s/%s_atmPrf_%d-%02d-%02d.nc' % (m, m, year, month, 1)
         fn = '%s/%s_atmPrf_%d-%02d-??.nc' % (m, m, year, month)
         d = pyg.openall(fn)
      except AssertionError as e:
         print(fn, 'not found.')
         continue
      ds.append(d)

   dm = pyg.dataset.concat(ds).sorted('time')

   return dm
# }}}

def open_nc_on_pres(year, month, pres, days=None):
# {{{
   from RO.merge_profiles_rrtm import grid_profiles_pres
   ds = []
   missions = ['cnofs', 'cosmic2013', 'metopa2016', 'grace', 'sacc', 'tsx']

   for m in missions[1:]:
      try:
         #fn = '%s/%s_atmPrf_%d-%02d-%02d.nc' % (m, m, year, month, 1)
         if days is None:
            fns = '%s/%s_atmPrf_%d-%02d-??.nc' % (m, m, year, month)
         else:
            fns = ['%s/%s_atmPrf_%d-%02d-%02d.nc' % (m, m, year, month, d) for d in days]
         d = pyg.openall(fns)
         tm = pyg.StandardTime(startdate=d.time.startdate, units='days', values=d.time[:], lat=d.OLat[:], lon=d.OLon[:])
         dt = d.replace_axes(time=tm)
         dp = grid_profiles_pres(dt, pres) 
         
      except AssertionError as e:
         print(fns, 'not found.')
         continue
      ds.append(dp)

   if len(ds) == 0: return None

   dm = pyg.concatenate(ds).sorted('time')

   return dm
# }}}

def test_grid():
# {{{
   from ro import ro
   from RO.merge_profiles_rrtm import grid_profiles_pres

   yr = 2010
   for mn in range(12, 13):
      dm = open_nc(yr, mn)

      tm = pyg.StandardTime(startdate=dm.time.startdate, units='days', values=dm.time[:], lat=dm.OLat[:], lon=dm.OLon[:])
      pres = pyg.Pres(10**np.linspace(2.6, 0.5, 200))

      dmt = dm.replace_axes(time=tm)
      dmpt = grid_profiles_pres(dmt, pres) 
      dmpt.load()

      dgt, dgtc = ro.grid2(dmpt.Temp, 10, 20, 1)
      dgz, dgzc = ro.grid2(dmpt.Alt, 10, 20, 1)

      ds = pyg.Dataset([dgt, dgz, dgtc.rename('count')])
      fn = rpath + 'raw/grid_test/ro_10by20_gridded_%04d-%02d.nc' % (yr, mn)
      print('Writing %s.' % fn)
      pyg.save(fn, ds)
# }}}

def test_kw_grid():
# {{{
   from ro import ro
   from RO.merge_profiles_rrtm import grid_profiles_pres
   import ERAInter as ei

   #pres = ei.an_on_p.pres
   pres = pyg.Pres(np.exp(np.linspace(np.log(450), np.log(4.5), 200)))

   yr = 2009
   for mn in range(1, 13):
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
      dmp = list(range(1, 4))

      dmm = open_nc_on_pres(yrm, mnm, pres, days = dmm)
      dm  = open_nc_on_pres(yr, mn, pres)
      dmp = open_nc_on_pres(yrp, mnp, pres, days = dmp)

      dm = pyg.concatenate([d for d in [dmm, dm, dmp] if d is not None])

      d0 = '1 %s %04d' % (pyg.timeaxis.months[mn], yr)
      d1 = '1 %s %04d' % (pyg.timeaxis.months[mnp], yrp)

      #tm = pyg.StandardTime(startdate=dm.time.startdate, units='days', values=dm.time[:], lat=dm.OLat[:], lon=dm.OLon[:])
      #pres = pyg.Pres(10**np.linspace(2.6, 0.5, 200))

      #dmt = dm.replace_axes(time=tm)
      #dmpt = grid_profiles_pres(dmt, pres) 
      dm = dm(lat = (-10, 10)).load()

      dgt, dgtc = ro.grid_kw_gauss(dm.Temp, d0, d1, 1, 0.25, 10., 1.)
      dgz, dgzc = ro.grid_kw_gauss(dm.Alt,  d0, d1, 1, 0.25, 10., 1.)
      
      ds = pyg.Dataset([dgt, dgz, dgtc.rename('count')])

      fn = rpath + 'raw/grid_test/ro_1deg_D10deg_T1d_kw_gauss_gridded_hvr_%04d-%02d.nc' % (yr, mn)
      print('Writing %s.' % fn)
      pyg.save(fn, ds)
# }}}

