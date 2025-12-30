import numpy as np, pylab as pyl
import pygeode as pyg
from rrtm import rrtmg
from ECMWF import ERAI as ei, ozone as o3
import MLS as mls
import WACCM as wc
import os

import consts as c
import read

def merge(year, month, day, nprof=-1):
# {{{
   out_pth = read.path + 'merged/mls/%04d-%02d/'
   out_fn = '%s_merged_%04d-%02d-%02d.nc'

   p0 = 1013.25

   # Construct hybrid level grid on which to construct merged data
   ei_etaf = ei.an.eta
   ei_ef = ei_etaf[:]
   ei_eh = (p0*ei.hyb_B_on_half + ei.hyb_A_on_half/100.) / p0
   ei_etah = pyg.Hybrid(ei_eh, A=ei.hyb_A_on_half, B=ei.hyb_B_on_half)

   ei_Af = ei_etaf.auxasvar('A') / 100.
   ei_Bf = ei_etaf.auxasvar('B')
   ei_Bh = ei_etah.auxasvar('B')

   etaf, etah = make_eta_grid(ei_ef, ei_eh)
   ef = pyg.Hybrid(etaf, A=etaf, B=0.)
   eh = pyg.Hybrid(etah, A=etah, B=0.).rename('etah')

   Bf = ei_Bf.interpolate('eta', ef)[:]
   Bh = ei_Bh.interpolate('eta', eh)[:]

   ef = pyg.Hybrid(etaf, A=etaf-Bf, B=Bf)
   eh = pyg.Hybrid(etah, A=etah-Bh, B=Bh).rename('etah')

   def tstr(tax, y, mn, d, h=0, mi=0, s=0, fmt='$H:$M:$S $d $b $a'):
      dt = dict(year=y, month=mn, day=d, hour=h, minute=mi, second=s)
      return tax.formatvalue(tax.date_as_val(dt), fmt)
   ltimes = (tstr(ei.an.time, year, month, day), tstr(ei.an.time, year, month, day+1))

   # Load ERA Interim data
   def load_ei(v): 
      v = v(time=ltimes)
      if v.hasaxis('eta'): v = v.sorted(eta = 1)
      return v.load()

   ei_ps  = pyg.exp(ei.lnsp.lnsp)/100.
   ei_ps  = load_ei(ei_ps)

   ei_t   = load_ei(ei.an.t)
   ei_qlw = load_ei(ei.fc_tend.qtendCSLW)
   ei_qsw = load_ei(ei.fc_tend.qtendCSSW)
   ei_skt = load_ei(ei.fc.skt)
   ei_fal = load_ei(ei.fc.fal)
 
   # Load WACCM data
   def load_wc(v): 
      tm = pyg.modeltime365range(ltimes[0], ltimes[1], 0.25, inc=True)
      v = (v + 0*tm).rename(v.name)
      if v.hasaxis('eta'): v = v.sorted(eta = 1)
      return v.load()

   wc_co2 = load_wc(wc.cld.CO2)
   wc_o3  = load_wc(wc.cld.O3 )
   wc_h2o = load_wc(wc.cld.H2O)
 
   # Load MLS data
   def load_mls(v): 
      v = v(year=year, month=month, day=day)
      if v.hasaxis('eta'): v = v.sorted(eta = 1)
      return v.load()

   mls_h2o = load_mls(1e-6*mls.h2o.H2O.unfill(-999.9))
   mls_o3  = load_mls(1e-6*mls.o3.O3.unfill(-999.9))
   mls_valid = True

   if len(mls_h2o.time) == 0 or len(mls_o3.time) == 0:
      mls_valid = False

   def sel_ei(v): 
      return v(s_time=tstr, s_lat=lat, s_lon=lon)

   def sel_mls(v): 
      ln = lon % 360.
      return v(s_time=tstr, s_lat=lat, s_lon=ln)

   def sel_wc(v): 
      return v(s_time=tstr, s_lat=lat)

   # Blending profiles (all in units of log-pressure km w/ p0 = 1000 hPa, H = 7 km)
   zt_T  = 40.   # Temperature
   zb_T  = 10.
   dzt_T = 3.
   dzb_T = 0.6

   zt_X  = 40.   # Constituents
   zb_X  = 11.
   dzt_X = 3.
   dzb_X = 0.6

   #missions = ['cnofs', 'cosmic2013', 'metopa2016', 'metopb2016', 'champ2016', 'grace', 'sacc', 'tsx']

   for m in read.missions:
      # Open GPS data
      try:
         d = read.open_nc(m, year, month, day)
      except AssertionError as e:
         print('No Matches: %s' % m)
         continue

      nlevs = len(ef)
      if nprof == -1: 
         npr = len(d.time)
      else:
         npr = nprof

      # Allocate output arrays
      P_f     = np.zeros((npr, nlevs),     'd')
      P_h     = np.zeros((npr, nlevs + 1), 'd')
      T_gps   = np.zeros((npr, nlevs),     'd')
      T_ei    = np.zeros((npr, nlevs),     'd')
      T_bl    = np.zeros((npr, nlevs),     'd')
      CO2_wc  = np.zeros((npr, nlevs),     'd')
      O3_wc   = np.zeros((npr, nlevs),     'd')
      H2O_wc  = np.zeros((npr, nlevs),     'd')
      O3_mls  = np.zeros((npr, nlevs),     'd')
      H2O_mls = np.zeros((npr, nlevs),     'd')
      O3_bl   = np.zeros((npr, nlevs),     'd')
      H2O_bl  = np.zeros((npr, nlevs),     'd')
      qlw_ei  = np.zeros((npr, nlevs),     'd')
      qsw_ei  = np.zeros((npr, nlevs),     'd')
      tsf_ei  = np.zeros((npr,),           'd')
      alb_ei  = np.zeros((npr,),           'd')

      for i in range(npr):
         p = d.Pres(si_time = i)
         # Discard data which are too closely spaced in pressure
         imsk = np.where((pyg.log(p).diff('level') < -1e-4)[:])[0]

         p = d.Pres(si_time = i, li_level=imsk)
         t = d.Temp(si_time = i, li_level=imsk) + c.zeroC
         lat = d.OLat[i]
         lon = d.OLon[i]
         tstr   = d.time.formatvalue(d.time[i], '$H:$M:$S $d $b $a')

         # Construct pressure grid
         ps = sel_ei(ei_ps)[()]
         ei_pr = ps * ei_Bf + ei_Af

         pf = pyg.Pres( p0 * etaf + (ps - p0) * Bf )
         ph = pyg.Pres( p0 * etah + (ps - p0) * Bh )

         P_f[i, :] = pf[:]
         P_h[i, :] = ph[:]

         zs = -7*pyg.log(pf/1000.)

         # Interpolate temperature data
         T_ei [i, :] = sel_ei(ei_t).interpolate('eta'  , pf, inx=pyg.log(ei_pr), outx=pyg.log(pf), d_below=0., d_above=0.)[:]
         T_gps[i, :] = t           .interpolate('level', pf, inx=pyg.log(p),     outx=pyg.log(pf), omit_nonmonotonic=True)[:]

         # Merge temperatures
         bl = (0.5 * (pyg.tanh((zs - zb_T)/dzb_T) + pyg.tanh((zt_T - zs)/dzt_T)))[:]

         msk = ~np.isnan(T_gps[i, :])
         T_bl[i, :  ] = T_ei[i, :]
         T_bl[i, msk] = T_bl[i, msk] * (1 - bl[msk]) + T_gps[i, msk] * bl[msk]

         # Interpolate constituent
         CO2_wc [i, :] = sel_wc(wc_co2 ).interpolate('eta',  pf, inx=pyg.log(wc_co2.eta), outx=pyg.log(pf), d_below=0., d_above=0.)[:]
         O3_wc  [i, :] = sel_wc(wc_o3  ).interpolate('eta',  pf, inx=pyg.log(wc_o3 .eta), outx=pyg.log(pf), d_below=0., d_above=0.)[:]
         H2O_wc [i, :] = sel_wc(wc_h2o ).interpolate('eta',  pf, inx=pyg.log(wc_h2o.eta), outx=pyg.log(pf), d_below=0., d_above=0.)[:]
    
         if mls_valid:
            O3_mls [i, :] = sel_mls(mls_o3 ).interpolate('pres', pf, inx=pyg.log(mls_o3 .pres), outx=pyg.log(pf), d_below=0., d_above=0.)[:]
            H2O_mls[i, :] = sel_mls(mls_h2o).interpolate('pres', pf, inx=pyg.log(mls_h2o.pres), outx=pyg.log(pf), d_below=0., d_above=0.)[:]
         else:
            O3_mls [i, :] = O3_wc[i, :]
            H2O_mls[i, :] = H2O_wc[i, :]

         # Merge constituents
         bl = (0.5 * (pyg.tanh((zs - zb_X)/dzb_X) + pyg.tanh((zt_X - zs)/dzt_X)))[:]

         msk = ~(np.isnan(O3_mls[i, :]) | (O3_mls[i, :] < 0.))
         O3_bl[i, :  ] = O3_wc[i, :]
         O3_bl[i, msk] = O3_bl[i, msk] * (1 - bl[msk]) + O3_mls[i, msk] * bl[msk]

         msk = ~(np.isnan(H2O_mls[i, :]) | (H2O_mls[i, :] < 0.))
         H2O_bl[i, :  ] = H2O_wc[i, :]
         H2O_bl[i, msk] = H2O_bl[i, msk] * (1 - bl[msk]) + H2O_mls[i, msk] * bl[msk]

         # Interpolate heating rates
         qlw_ei[i, :] = sel_ei(ei_qlw).interpolate('eta', pf, inx=pyg.log(ei_pr), outx=pyg.log(pf), d_below=0., d_above=0.)[:]
         qsw_ei[i, :] = sel_ei(ei_qsw).interpolate('eta', pf, inx=pyg.log(ei_pr), outx=pyg.log(pf), d_below=0., d_above=0.)[:]

         # Store skin temperature and albedo for appropriate gridpoint/time
         tsf_ei[i] = sel_ei(ei_skt)[()]
         alb_ei[i] = sel_ei(ei_fal)[()]

      time = d.time(i_time = (0, npr))

      pf      = pyg.Var((time, ef), values = P_f,     name = 'p')
      ph      = pyg.Var((time, eh), values = P_h,     name = 'ph')
      tei     = pyg.Var((time, ef), values = T_ei,    name = 'Tei')
      tgps    = pyg.Var((time, ef), values = T_gps,   name = 'Tgps')
      tbl     = pyg.Var((time, ef), values = T_bl,    name = 'Tbl')
      CO2_wc  = pyg.Var((time, ef), values = CO2_wc,  name = 'CO2w')
      O3_wc   = pyg.Var((time, ef), values = O3_wc,   name = 'O3w')
      O3_mls  = pyg.Var((time, ef), values = O3_mls,  name = 'O3mls')
      O3_bl   = pyg.Var((time, ef), values = O3_bl,   name = 'O3bl')
      H2O_wc  = pyg.Var((time, ef), values = H2O_wc,  name = 'H2Ow')
      H2O_mls = pyg.Var((time, ef), values = H2O_mls, name = 'H2Omls')
      H2O_bl  = pyg.Var((time, ef), values = H2O_bl,  name = 'H2Obl')
      qLW_ei  = pyg.Var((time, ef), values = qlw_ei,  name = 'lwhrei')
      qSW_ei  = pyg.Var((time, ef), values = qsw_ei,  name = 'swhrei')
      tsfc    = pyg.Var((time,),    values = tsf_ei,  name = 'Tsfc')
      alb     = pyg.Var((time,),    values = alb_ei,  name = 'alb')

      vs = [pf, ph, tei, tgps, tbl, CO2_wc, O3_wc, O3_mls, O3_bl, H2O_wc, H2O_mls, H2O_bl, \
            qLW_ei, qSW_ei, tsfc, alb]
      
      vs += [v(i_time = (0, npr)) for v in [d.OLat, d.OLon, d.GPSsatID]]

      ds = pyg.Dataset(vs)

      pth = out_pth % (year, month)
      if not os.path.exists(pth): os.makedirs(pth)
      fn = pth + out_fn % (m, year, month, day)
      print('Writing %s.' % fn)
      pyg.save(fn, ds, version = 4, compress = True)
# }}}

def merge_EI_comp(year, month, day, nprof=-1):
# {{{
   out_pth = read.path + 'merged/ei/%04d-%02d/'
   out_fn = '%s_ei_comp_merged_%04d-%02d-%02d.nc'

   p0 = 1013.25

   # Construct hybrid level grid on which to construct merged data
   ei_etaf = ei.an.eta
   ei_ef = ei_etaf[:]
   ei_eh = (p0*ei.hyb_B_on_half + ei.hyb_A_on_half/100.) / p0
   ei_etah = pyg.Hybrid(ei_eh, A=ei.hyb_A_on_half, B=ei.hyb_B_on_half)

   ei_Af = ei_etaf.auxasvar('A') / 100.
   ei_Bf = ei_etaf.auxasvar('B')
   ei_Bh = ei_etah.auxasvar('B')

   etaf, etah = make_eta_grid(ei_ef, ei_eh)
   ef = pyg.Hybrid(etaf, A=etaf, B=0.)
   eh = pyg.Hybrid(etah, A=etah, B=0.).rename('etah')

   Bf = ei_Bf.interpolate('eta', ef)[:]
   Bh = ei_Bh.interpolate('eta', eh)[:]

   ef = pyg.Hybrid(etaf, A=etaf-Bf, B=Bf)
   eh = pyg.Hybrid(etah, A=etah-Bh, B=Bh).rename('etah')

   def tstr(tax, y, mn, d, h=0, mi=0, s=0, fmt='$H:$M:$S $d $b $a'):
      dt = dict(year=y, month=mn, day=d, hour=h, minute=mi, second=s)
      return tax.formatvalue(tax.date_as_val(dt), fmt)
   ltimes = (tstr(ei.an.time, year, month, day), tstr(ei.an.time, year, month, day+1))

   # Load ERA Interim data
   def load_ei(v): 
      v = v(time=ltimes)
      if v.hasaxis('eta'): v = v.sorted(eta = 1)
      return v.load()

   ei_ps  = pyg.exp(ei.lnsp.lnsp)/100.
   ei_ps  = load_ei(ei_ps)

   ei_h2o = load_ei(c.mdw*ei.an.q) # Convert specific humidity to volume mixing ratio)
 
   # Load ERAInterim ozone data
   def load_ei_o3(v): 
      tm = pyg.modeltime365range(ltimes[0], ltimes[1], 0.25, inc=True)
      v = (v + 0*tm).rename(v.name)
      if v.hasaxis('pres'): v = v.sorted(pres = 1)
      return v.load()

   ei_o3 = load_ei_o3(o3.o3cd.OZO*1e-6) # Convert to volume mixing ratio (from ppmv)
 
   def sel_ei(v): 
      return v(s_time=tstr, s_lat=lat, s_lon=lon)

   def sel_ei_o3(v): 
      return v(s_time=tstr, s_lat=lat)

   # Blending profiles (all in units of log-pressure km w/ p0 = 1000 hPa, H = 7 km)
   zt_T  = 40.   # Temperature
   zb_T  = 10.
   dzt_T = 3.
   dzb_T = 0.6

   zt_X  = 40.   # Constituents
   zb_X  = 11.
   dzt_X = 3.
   dzb_X = 0.6

   #missions = ['cnofs', 'cosmic2013', 'metopa2016', 'metopb2016', 'champ2016', 'grace', 'sacc', 'tsx']

   for m in read.missions:
      # Open GPS data
      try:
         d = read.open_nc(m, year, month, day)
      except AssertionError as e:
         print('No Matches: %s' % m)
         continue

      nlevs = len(ef)
      if nprof == -1: 
         npr = len(d.time)
      else:
         npr = nprof

      # Allocate output arrays
      P_f     = np.zeros((npr, nlevs),     'd')
      P_h     = np.zeros((npr, nlevs + 1), 'd')
      O3_ei   = np.zeros((npr, nlevs),     'd')
      H2O_ei  = np.zeros((npr, nlevs),     'd')

      for i in range(npr):
         p = d.Pres(si_time = i)
         # Discard data which are too closely spaced in pressure
         imsk = np.where((pyg.log(p).diff('level') < -1e-4)[:])[0]

         p = d.Pres(si_time = i, li_level=imsk)
         t = d.Temp(si_time = i, li_level=imsk) + c.zeroC
         lat = d.OLat[i]
         lon = d.OLon[i]
         tstr   = d.time.formatvalue(d.time[i], '$H:$M:$S $d $b $a')

         # Construct pressure grid
         ps = sel_ei(ei_ps)[()]
         ei_pr = ps * ei_Bf + ei_Af

         pf = pyg.Pres( p0 * etaf + (ps - p0) * Bf )
         ph = pyg.Pres( p0 * etah + (ps - p0) * Bh )

         P_f[i, :] = pf[:]
         P_h[i, :] = ph[:]

         # Interpolate constituent
         O3_ei  [i, :] = sel_ei_o3(ei_o3  ).interpolate('pres', pf, inx=pyg.log(ei_o3.pres), outx=pyg.log(pf), d_below=0., d_above=0.)[:]
         H2O_ei [i, :] = sel_ei   (ei_h2o ).interpolate('eta',  pf, inx=pyg.log(ei_pr)     , outx=pyg.log(pf), d_below=0., d_above=0.)[:]

      time = d.time(i_time = (0, npr))

      pf      = pyg.Var((time, ef), values = P_f,     name = 'p')
      ph      = pyg.Var((time, eh), values = P_h,     name = 'ph')
      O3_ei   = pyg.Var((time, ef), values = O3_ei,   name = 'O3ei')
      H2O_ei  = pyg.Var((time, ef), values = H2O_ei,  name = 'H2Oei')

      vs = [pf, ph, O3_ei, H2O_ei]
      
      vs += [v(i_time = (0, npr)) for v in [d.OLat, d.OLon, d.GPSsatID]]

      ds = pyg.Dataset(vs)

      pth = out_pth % (year, month)
      if not os.path.exists(pth): os.makedirs(pth)
      fn = pth + out_fn % (m, year, month, day)
      print('Writing %s.' % fn)
      pyg.save(fn, ds, version = 4, compress = True)
# }}}

def make_eta_grid(eta_r, eta_rh):
# {{{
   N = 166
   Ne = len(eta_r)
   eta = np.zeros(N, 'd')
   etah = np.zeros(N + 1, 'd')

   # Attach ERA grid below 400 and above 3
   e0_s = 0.003
   e0_t = 0.55
   i_eit = Ne - np.where(eta_r > e0_t)[0][0]
   i_eim = np.where(eta_r < e0_s)[0][-1]

   #print i_eit, eta_r[-i_eit]
   #print i_eim, eta_r[i_eim]

   eta[-i_eit:]  = eta_r[-i_eit:]
   eta[:i_eim]   = eta_r[:i_eim]
   etah[-i_eit:] = eta_rh[-i_eit:]
   etah[:i_eim]  = eta_rh[:i_eim]

   Ni = N - (i_eit + i_eim)

   Ng = Ni
   #iter = 0

   kap = 1.
   nkap = 1.

   # Target resolution: 600 m at 400 hPa, 300 m from 150 hPa to 10 hPa, 1600 m at 3 hPa 
   z_tp = 8.
   z_sp = 40. 
   z_wt = 2.1
   z_ws = 5.
   dzt = 0.5
   dzg = 0.2
   dzm = 2.0

   z0 = -7 * np.log(eta_r[-i_eit])
   z1 = -7 * np.log(eta_r[i_eim - 1])

   Z = lambda z : dzt + 0.5 * (dzg - dzt) * (1 + np.tanh((z - z_tp)/ (z_wt))) + 0.5 * (dzm - dzg) * (1 + np.tanh((z - z_sp) / z_ws))

   #while (Ng != Ni or iter == 0) and iter < 10:
   #   kap = nkap

   #   Z = lambda z : dzt + 0.5 * (dzg - dzt) * (1 + np.tanh((z - z_tp)/ (nkap*z_wt))) + 0.5 * (dzm - dzg) * (1 + np.tanh((z - z_sp) / z_ws))

   #   z = z1
   #   Ng = 0
   #   while z > z0:
   #      Ng += 1
   #      z -= Z(z)

   #   nkap = ((z1 - z) *  Ni) / (Ng * (z1 - z0))
   #   iter += 1
   #   print iter, Ng, Ni, kap, nkap, z

   z = z1
   i = i_eim - 1
   #print np.exp(-z/7.), eta[i], e0_s
   while z > z0 + 0.7*dzt:
      i += 1
      z -= Z(z)
      e = np.exp(-z/7.)
      eta[i] = e
   #print Ni, i, eta[i], eta[i+1]

   etah[0] = 1e-5
   etah[i_eim:i+1] = np.sqrt(eta[i_eim-1:i] * eta[i_eim:i+1])

   return eta, etah

   zs = -7*np.log(eta)
   zh = -7*np.log(etah)
   f = pyl.figure(4)
   f.clf()
   ax = f.add_subplot(121)
   ax.plot(Z(zs), zs, 'r')
   ax.plot(-np.diff(zs), zs[:-1], 'ko')
   ax = f.add_subplot(122)
   ax.plot(0*zs, zs, 'ro')
   ax.plot(0*zh, zh, 'k+')

   return eta, etah, Z
# }}}

def run_rrtm(year, month, day, nprof=-1, iprof=0):
# {{{
   from rrtm import rrtmg, astr

   #path = '/data2/RO/merged/mls/'
   #d = pyg.openall(path + '*_%04d-%02d-%02d.nc' % (year, month, day), sorted=False).sorted('time')

   d = read.open_merged(year, month, day, dset='mls')

   fnout = read.path + 'rad/mls/rrtm_%04d-%02d-%02d.nc' % (year, month, day)

   def sanitize(a): return np.asfortranarray(a, 'd')

   if nprof == -1:  nprof = len(d.time)
   print(nprof)

   d2 = d(i_time=(iprof, iprof + nprof))
   time = d2.time

   ph   = sanitize(d2.ph [:, :])
   pf   = sanitize(d2.p  [:, :])
   tf   = sanitize(d2.Tbl[:, :])
   tei  = sanitize(d2.Tei[:, :])

   tsfc = sanitize(d2.Tsfc[:])
   alb  = sanitize(d2.alb [:])

   co2  = sanitize(d2.CO2w [:, :])
   h2o  = sanitize(d2.H2Obl[:, :])
   o3   = sanitize(d2.O3bl [:, :])
   o3w  = sanitize(d2.O3w  [:, :])
   zero = sanitize(d2.O3bl [:, :]*0.)

   lats = sanitize(d2.OLat [:])
   lons = sanitize(d2.OLon [:])

   t    = pyg.timeutils.reltime(time, dict(year = 2000, month = 1, day = 1))
   hr   = 24 * np.fmod(t, 1.)
   dy   = t.astype('i')

   cosz = sanitize(np.clip(astr.zenith(lats, lons, hr, dy), 0., 1.))
   emi  = sanitize(c.emi * np.ones(nprof, 'd'))

   rrtmg.init(c.cpair)

   ret_lw     = rrtmg.rrtmg_lw(pf, ph, tf, tsfc, emi, co2,  h2o,  o3)
   ret_lw_co2 = rrtmg.rrtmg_lw(pf, ph, tf, tsfc, emi, zero, h2o,  o3)
   ret_lw_h2o = rrtmg.rrtmg_lw(pf, ph, tf, tsfc, emi, co2,  zero, o3)
   ret_lw_o3  = rrtmg.rrtmg_lw(pf, ph, tf, tsfc, emi, co2,  h2o,  zero)

   ret_lw_ei  = rrtmg.rrtmg_lw(pf, ph, tei, tsfc, emi, co2, h2o,  o3)
   ret_lw_o3w = rrtmg.rrtmg_lw(pf, ph, tf , tsfc, emi, co2, h2o,  o3w)

   ret_sw     = rrtmg.rrtmg_sw(pf, ph, tf, tsfc, c.scon, cosz, alb, co2,  h2o,  o3)
   ret_sw_co2 = rrtmg.rrtmg_sw(pf, ph, tf, tsfc, c.scon, cosz, alb, zero, h2o,  o3)
   ret_sw_h2o = rrtmg.rrtmg_sw(pf, ph, tf, tsfc, c.scon, cosz, alb, co2,  zero, o3)
   ret_sw_o3  = rrtmg.rrtmg_sw(pf, ph, tf, tsfc, c.scon, cosz, alb, co2,  h2o,  zero)
   ret_sw_o3w = rrtmg.rrtmg_sw(pf, ph, tf, tsfc, c.scon, cosz, alb, co2,  h2o,  o3w)

   cosz = pyg.Var((time,),      values=cosz,          name='cosz')

   Pf       = pyg.Var((time, d.eta) , name='pf'   , values=pf)
   Ph       = pyg.Var((time, d.etah), name='ph'   , values=ph)
   Tf       = pyg.Var((time, d.eta) , name='tf'   , values=tf)
   Tsfc     = pyg.Var((time,)       , name='tsfc' , values=tsfc)
   Alb      = pyg.Var((time,)       , name='alb'  , values=alb)
   CosZ     = pyg.Var((time,)       , name='cosz' , values=cosz)
   Emis     = pyg.Var((time,)       , name='emis' , values=emi)
   Day      = pyg.Var((time,)       , name='day'  , values=dy)
   Hr       = pyg.Var((time,)       , name='hour' , values=hr)
   Lat      = pyg.Var((time,)       , name='lat'  , values=lats)
   Lon      = pyg.Var((time,)       , name='lon'  , values=lons)
   CO2      = pyg.Var((time, d.eta) , name='co2'  , values=co2)
   O3       = pyg.Var((time, d.eta) , name='o3'   , values=o3)
   H2O      = pyg.Var((time, d.eta) , name='h2o'  , values=h2o)

   lwhr    = pyg.Var((time, d.eta) , name='lwhr'   , values=ret_lw['lwhr']  )
   dflw    = pyg.Var((time, d.etah), name='dflxlw' , values=ret_lw['dflxlw'])
   uflw    = pyg.Var((time, d.etah), name='uflxlw' , values=ret_lw['uflxlw'])
   swhr    = pyg.Var((time, d.eta) , name='swhr'   , values=ret_sw['swhr']  )
   dfsw    = pyg.Var((time, d.etah), name='dflxsw' , values=ret_sw['dflxsw'])
   ufsw    = pyg.Var((time, d.etah), name='uflxsw' , values=ret_sw['uflxsw'])

   lwhrco2 = pyg.Var((time, d.eta) , name='lwhrco2', values=ret_lw['lwhr'] - ret_lw_co2['lwhr'])
   lwhro3  = pyg.Var((time, d.eta) , name='lwhro3' , values=ret_lw['lwhr'] - ret_lw_o3 ['lwhr'])
   lwhrh2o = pyg.Var((time, d.eta) , name='lwhrh2o', values=ret_lw['lwhr'] - ret_lw_h2o['lwhr'])
   lwhrtei = pyg.Var((time, d.eta) , name='lwhrtei', values=ret_lw_ei ['lwhr'])
   lwhro3w = pyg.Var((time, d.eta) , name='lwhro3w', values=ret_lw_o3w['lwhr'])

   swhrco2 = pyg.Var((time, d.eta) , name='swhrco2', values=ret_sw['swhr'] - ret_sw_co2['swhr'])
   swhro3  = pyg.Var((time, d.eta) , name='swhro3' , values=ret_sw['swhr'] - ret_sw_o3 ['swhr'])
   swhrh2o = pyg.Var((time, d.eta) , name='swhrh2o', values=ret_sw['swhr'] - ret_sw_h2o['swhr'])
   swhro3w = pyg.Var((time, d.eta) , name='swhro3w', values=ret_sw_o3w['swhr'])

   #ds = pyg.Dataset([cosz, lwhr, dflw, uflw])
   #ds = pyg.Dataset([cosz, swhr, dfsw, ufsw])
   ds = pyg.Dataset([cosz, lwhr, swhr, dflw, uflw, dfsw, ufsw,\
                     lwhrco2, lwhro3, lwhrh2o, swhrco2, swhro3, swhrh2o, \
                     lwhrtei, lwhro3w, swhro3w])

   ds = pyg.Dataset([Pf, Ph, Tf, Tsfc, Emis, CO2, H2O, O3, CosZ, Alb, Day, Hr, Lat, Lon,\
                     lwhr, swhr, dflw, uflw, dfsw, ufsw,\
                     lwhrco2, lwhro3, lwhrh2o, swhrco2, swhro3, swhrh2o, \
                     lwhrtei, lwhro3w, swhro3w])

   return ds
   print('Writing %s' % fnout)

   pyg.save(fnout, ds)
# }}}

def run_pyracc(year, month, day, nprof=-1, iprof=0):
# {{{
   from rrtm import rrtmg, astr
   import pyracc as pyr

   d   = read.open_merged(year, month, day, dset='mls')
   dei = read.open_merged(year, month, day, dset='ei')

   fnout = read.path + 'rad/pyr/%04d/rrtm_%04d-%02d-%02d.nc' % (year, year, month, day)

   nlev = len(d.eta)
   if nprof == -1:  nprof = len(d.time)

   d2 = (d + dei.O3ei + dei.H2Oei)(i_time=(iprof, iprof + nprof))
   time = d2.time

   #alb  = sanitize(d2.alb [:])

   #o3w  = sanitize(d2.O3w  [:, :])

   #lats = sanitize(d2.OLat [:])
   #lons = sanitize(d2.OLon [:])

   #t    = pyg.timeutils.reltime(time, dict(year = 2000, month = 1, day = 1))
   #hr   = 24 * np.fmod(t, 1.)
   #dy   = t.astype('i')

   #cosz = sanitize(np.clip(astr.zenith(lats, lons, hr, dy), 0., 1.))

   # Set up pyraccoons LW calculation
   lw   = pyr.lw('RRTMG', nlev, nprof, cpair = c.cpair)
   lw.pres  = d2.p [()]
   lw.phalf = d2.ph[()]

   lw.T     = d2.Tbl [()]
   lw.Tsfc  = d2.Tsfc[()]

   lw.emis  = c.emi

   # Base case: blended temps, EI constituents
   lw.CO2 = d2.CO2w[()];  lw.H2O = d2.H2Obl[()];  lw.O3 = d2.O3bl[()]
   lw.H2O[lw.H2O < 0.] = 1e-6
   ret_lw = lw.run()

   # Species dependent rates: CO2
   lw.CO2 = 0.;           lw.H2O = d2.H2Obl[()];  lw.O3 = d2.O3bl[()]
   lw.H2O[lw.H2O < 0.] = 1e-6
   ret_lw_co2 = lw.run()

   # Species dependent rates: H2O
   lw.CO2 = d2.CO2w[()];  lw.H2O = 0.;            lw.O3 = d2.O3bl[()]
   ret_lw_h2o = lw.run()

   # Species dependent rates: O3
   lw.CO2 = d2.CO2w[()];  lw.H2O = d2.H2Obl[()];  lw.O3 = 0.
   lw.H2O[lw.H2O < 0.] = 1e-6
   ret_lw_o3  = lw.run()

   # Wrap into pygeode variables
   lwhr    = pyg.Var((time, d.eta) , name='lwhr'   , values=ret_lw.lwhr  )
   dflw    = pyg.Var((time, d.etah), name='dflxlw' , values=ret_lw.dflxlw)
   uflw    = pyg.Var((time, d.etah), name='uflxlw' , values=ret_lw.uflxlw)

   lwhrco2 = pyg.Var((time, d.eta) , name='lwhrco2', values=ret_lw.lwhr - ret_lw_co2.lwhr)
   lwhro3  = pyg.Var((time, d.eta) , name='lwhro3' , values=ret_lw.lwhr - ret_lw_o3 .lwhr)
   lwhrh2o = pyg.Var((time, d.eta) , name='lwhrh2o', values=ret_lw.lwhr - ret_lw_h2o.lwhr)

   #lwhrtei = pyg.Var((time, d.eta) , name='lwhrtei', values=ret_lw_ei .lwhr)
   #lwhro3w = pyg.Var((time, d.eta) , name='lwhro3w', values=ret_lw_o3w.lwhr)


   # Set up pyraccoons SW calculation
   sw   = pyr.sw('RRTMG', nlev, nprof, cpair = c.cpair)
   sw.pres  = d2.p [()]
   sw.phalf = d2.ph[()]

   sw.T     = d2.Tbl [()]
   sw.Tsfc  = d2.Tsfc[()]

   sw.scon  = c.scon
   sw.alb   = d2.alb[()]

   lats = d2.OLat[()]
   lons = d2.OLon[()]

   t    = pyg.timeutils.reltime(time, dict(year = 2000, month = 1, day = 1))
   hr   = 24 * np.fmod(t, 1.)
   dy   = t.astype('i')

   sw.cosz  = np.clip(astr.zenith(lats, lons, hr, dy), 1e-6, 1.)


   # Base case: blended temps, constituents
   sw.CO2 = d2.CO2w[()];  sw.H2O = d2.H2Obl[()];  sw.O3 = d2.O3bl[()]
   sw.H2O[sw.H2O < 0.] = 1e-6
   ret_sw = sw.run()

   # Species dependent rates: CO2
   sw.CO2 = 0.;           sw.H2O = d2.H2Obl[()];  sw.O3 = d2.O3bl[()]
   sw.H2O[sw.H2O < 0.] = 1e-6
   ret_sw_co2 = sw.run()

   # Species dependent rates: H2O
   sw.CO2 = d2.CO2w[()];  sw.H2O = 0.;            sw.O3 = d2.O3bl[()]
   ret_sw_h2o = sw.run()

   # Species dependent rates: O3
   sw.CO2 = d2.CO2w [()]; sw.H2O = d2.H2Obl[()];  sw.O3 = 0.
   sw.H2O[sw.H2O < 0.] = 1e-6
   ret_sw_o3  = sw.run()

   # Wrap into pygeode dataset
   cosz     = pyg.Var((time,),        name='cosz'   , values=sw.cosz)

   swhr    = pyg.Var((time, d.eta) , name='swhr'   , values=ret_sw.swhr  )
   dfsw    = pyg.Var((time, d.etah), name='dflxsw' , values=ret_sw.dflxsw)
   ufsw    = pyg.Var((time, d.etah), name='uflxsw' , values=ret_sw.uflxsw)

   swhrco2 = pyg.Var((time, d.eta) , name='swhrco2', values=ret_sw.swhr - ret_sw_co2.swhr)
   swhro3  = pyg.Var((time, d.eta) , name='swhro3' , values=ret_sw.swhr - ret_sw_o3 .swhr)
   swhrh2o = pyg.Var((time, d.eta) , name='swhrh2o', values=ret_sw.swhr - ret_sw_h2o.swhr)
   #swhro3w = pyg.Var((time, d.eta) , name='swhro3w', values=ret_sw_o3w.swhr)

   ds = pyg.Dataset([cosz, lwhr, swhr, dflw, uflw, dfsw, ufsw,\
                     lwhrco2, lwhro3, lwhrh2o, swhrco2, swhro3, swhrh2o])#, \
                     #lwhrtei, lwhro3w, swhro3w])

   print('Writing %s' % fnout)

   pyg.save(fnout, ds)
# }}}

def run_profile(phalf, pfull, temp, co2, h2o, o3, tsfc, alb, cosz):
# {{{
   from rrtm import rrtmg, astr

   def order(v):
      assert v.naxes <= 2
      assert v.hasaxis('time')
      if v.naxes == 2:
         assert v.hasaxis('pres')
         return v.sorted(pres = 1).transpose('time', 'pres')
      else:
         return v

   def sanitize(v): 
      return np.asfortranarray(v[:], 'd')

   pfull = order(pfull)
   phalf = order(phalf)
   temp  = order(temp)
   co2   = order(co2)
   h2o   = order(h2o)
   o3    = order(o3)
   tsfc  = order(tsfc)
   alb   = order(alb)

   time = temp.time

   nprof = len(time)

   emi   = c.emi * np.ones(nprof, 'd')
   cosz  = cosz  * np.ones(nprof, 'd')

   ph = sanitize(phalf)
   pf = sanitize(pfull)
   tf = sanitize(temp)
   co2 = sanitize(co2)
   h2o = sanitize(h2o)
   o3  = sanitize(o3)
   tsfc = sanitize(tsfc)
   emi  = sanitize(emi)
   alb  = sanitize(alb)
   cosz = sanitize(cosz)

   #print pf.shape, ph.shape, tf.shape, tsfc.shape, o3.shape
   #print cosz.shape, alb.shape, emi.shape
   #print pf.shape, ph.shape, tf.shape, tsfc.shape, o3.shape, tsfc.shape
   #return pf, ph, tf, tsfc, co2, h2o, o3

   rrtmg.init(c.cpair)

   ret_lw = rrtmg.rrtmg_lw(pf, ph, tf, tsfc, emi, co2, h2o, o3)

   ret_sw = rrtmg.rrtmg_sw(pf, ph, tf, tsfc, c.scon, cosz, alb, co2, h2o, o3)

   pfl = pfull.pres.rename('pres')
   phf = phalf.pres.rename('phalf')

   pf   = pyg.Var((time, pfl), values=pf,            name='pf')

   ph   = pyg.Var((time, phf), values=ph,            name='ph')
   tf   = pyg.Var((time, pfl), values=tf,            name='t')
   tsfc = pyg.Var((time,),      values=tsfc,          name='skt')
   alb  = pyg.Var((time,),      values=alb,           name='alb')
   cosz = pyg.Var((time,),      values=cosz,          name='cosz')
   co2  = pyg.Var((time, pfl), values=co2,           name='CO2')
   h2o  = pyg.Var((time, pfl), values=h2o,           name='H2O')
   o3   = pyg.Var((time, pfl), values=o3,            name='O3')

   lwhr = pyg.Var((time, pfl), values=ret_lw['lwhr'],   name='lwhr')
   swhr = pyg.Var((time, pfl), values=ret_sw['swhr'],   name='swhr')
   dflw = pyg.Var((time, phf), values=ret_lw['dflxlw'], name='dflxlw')
   uflw = pyg.Var((time, phf), values=ret_lw['uflxlw'], name='uflxlw')
   dfsw = pyg.Var((time, phf), values=ret_sw['dflxsw'], name='dflxsw')
   ufsw = pyg.Var((time, phf), values=ret_sw['uflxsw'], name='uflxsw')

   ds = pyg.Dataset([pf, ph, tf, tsfc, alb, cosz, co2, h2o, o3, lwhr, swhr, dflw, uflw, dfsw, ufsw])

   return ds
# }}}

def load_ei_hr_profile(pf, year, month, day, hour, lat, lon):
# {{{
   # Load ERA Interim data
   def load_ei(v): 
      v = v(year=year, month=month, day=day)
      if v.hasaxis('eta'): v = v.sorted('eta')
      return v

   ei_ps = pyg.exp(ei.lnsp.lnsp)/100.
   ei_p = ei_ps * ei.an.eta.auxasvar('B') + ei.an.eta.auxasvar('A')/100.
   ei_p = ei_p.transpose('time', 'eta', 'lat', 'lon').rename('pres')

   ei_p   = load_ei(ei_p)
   ei_t   = load_ei(ei.an.t)
   ei_qlw = load_ei(ei.fc_tend.qtendCSLW)
   ei_qsw = load_ei(ei.fc_tend.qtendCSSW)
   ei_skt = load_ei(ei.fc.skt)
   ei_fal = load_ei(ei.fc.fal)

   timestr = '%d:00 %d %s %d' % (hour, day, pyg.timeaxis.months[month], year)

   def sel_ei(v): 
       return v(s_time=timestr, s_lat=lat, s_lon=lon).load()

   def interp_ei(v):
       vp = sel_ei(v)
       return vp.interpolate('eta', pf.eta, inx=pyg.log(ei_pr), outx=pyg.log(pf), d_below=0.).load()

   ei_pr = sel_ei(ei_p)
   
   t_ei   = interp_ei(ei_t)
   qlw_ei = interp_ei(ei_qlw)
   qsw_ei = interp_ei(ei_qsw)

   return pyg.Dataset([t_ei, qlw_ei, qsw_ei])
# }}}

def test_profile():
# {{{
   zh = np.arange(0, 57, 0.3)
   zf = (zh[:-1] + zh[1:]) / 2.

   nlevs = len(zf)

   phalf = pyg.Pres(1000.*np.exp(-zh/7.))
   pfull = pyg.Pres(1000.*np.exp(-zf/7.))

   year = 2010
   month = 6
   day = 15
   iprf = 9

   # Open GPS data
   #d = read.open_nc(year, month, [day])
   d = read.open_nc_on_pres(year, month, pfull, [day])

   # Open WACCM data
   wc = pyg.openall('/local/storage/WACCM/waccm.ccmi79to14.clim.*.nc')

   co2 = wc.CO2(month=month).squeeze()
   o3  = wc.O3 (month=month).squeeze()
   h2o = wc.H2O(month=month).squeeze()

   # Load ERA Interim data
   def load_ei(v): 
      v = v(year=year, month=month, day=day)
      if v.hasaxis('eta'): v = v.sorted('eta')
      return v

   ei_ps = pyg.exp(ei.lnsp.lnsp)/100.
   ei_p = ei_ps * ei.an.eta.auxasvar('B') + ei.an.eta.auxasvar('A')/100.
   ei_p = ei_p.transpose('time', 'eta', 'lat', 'lon').rename('pres')

   ei_p   = load_ei(ei_p)
   ei_t   = load_ei(ei.an.t)
   ei_qlw = load_ei(ei.fc_tend.qtendCSLW)
   ei_qsw = load_ei(ei.fc_tend.qtendCSSW)
   ei_skt = load_ei(ei.fc.skt)
   ei_fal = load_ei(ei.fc.fal)

   timestr = d.time.formatvalue(d.time[iprf], '$H:$M:$S $d $b $Y')
   lat = d.OLat[iprf]
   lon = d.OLon[iprf]

   def sel_ei(v): return v(s_time=timestr, s_lat=lat, s_lon=lon).load()
   def sel_wc(v): return v(s_lat=lat).load()

   ei_pr = sel_ei(ei_p)
   
   t_gps = d.Temp(si_time=iprf)[:] + c.zeroC
   t_ei = sel_ei(ei_t).interpolate('eta', pfull, inx=pyg.log(ei_pr), outx=pyg.log(pfull), d_below=0.).load()


   zb = 7.*np.log(250./1000.)
   zt = 7.*np.log( 10./1000.)
   dz = 2.

   imsk = ~np.isnan(t_gps)

   bl = 1. + 0.5*(np.tanh((zf - zb) / dz) - np.tanh((zf - zb) / dz))

   temp = t_ei[:].copy()
   temp[imsk] = temp[imsk] * (1. - bl[imsk]) + t_gps[imsk] * bl[imsk]

   tax = d.time(i_time=iprf)

   co2_wc = sel_wc(co2).interpolate('eta', pfull, inx=pyg.log(wc.eta/1000.), outx=pyg.log(pfull/1000.), d_below=0., d_above=0.)
   o3_wc  = sel_wc(o3 ).interpolate('eta', pfull, inx=pyg.log(wc.eta/1000.), outx=pyg.log(pfull/1000.), d_below=0., d_above=0.)
   h2o_wc = sel_wc(h2o).interpolate('eta', pfull, inx=pyg.log(wc.eta/1000.), outx=pyg.log(pfull/1000.), d_below=0., d_above=0.)

   temp = pyg.Var((tax, pfull), values=temp   .reshape(1, -1), name = 't')
   phalf = 0*tax + phalf
   pfull = 0*tax + pfull
   co2 = co2_wc.extend(0, tax)
   o3  =  o3_wc.extend(0, tax)
   h2o = h2o_wc.extend(0, tax)
   tsfc = pyg.Var((tax,), values=[sel_ei(ei_skt)[()]], name = 'tsfc')
   alb  = pyg.Var((tax,), values=[sel_ei(ei_fal)[()]], name = 'fal')
   
   ds = run_profile(phalf, pfull, temp, co2, h2o, o3, tsfc, alb, 0.6)

   tgps = pyg.Var((tax, ds.pres), values=t_gps  .reshape(1, -1), name = 'tgps')
   tei  = pyg.Var((tax, ds.pres), values=t_ei[:].reshape(1, -1), name = 'tei')
   qlw_ei = sel_ei(ei_qlw).interpolate('eta', ds.pres, inx=pyg.log(ei_pr), outx=pyg.log(ds.pres), d_below=0.).load()
   qsw_ei = sel_ei(ei_qsw).interpolate('eta', ds.pres, inx=pyg.log(ei_pr), outx=pyg.log(ds.pres), d_below=0.).load()
   qlw_ei = qlw_ei.extend(0, tax)
   qsw_ei = qsw_ei.extend(0, tax)

   return ds + pyg.Dataset([tgps, tei, qlw_ei, qsw_ei])
# }}}

if __name__ == '__main__':
   import sys
   if len(sys.argv) > 1:
      cmd = sys.argv[1]

      if cmd == 'merge':
         year = int(sys.argv[2])
         day  = int(sys.argv[3])

         tm = pyg.standardtimen('1 Jan 2000', 1)
         dt = tm.val_as_date(tm.date_as_val(dict(year=year, day=day, month=1)))

         merge(dt['year'], dt['month'], dt['day'])

      elif cmd == 'mergeei':
         year = int(sys.argv[2])
         day  = int(sys.argv[3])

         tm = pyg.standardtimen('1 Jan 2000', 1)
         dt = tm.val_as_date(tm.date_as_val(dict(year=year, day=day, month=1)))

         merge_EI_comp(dt['year'], dt['month'], dt['day'])

      elif cmd == 'rrtm':
         year = int(sys.argv[2])
         day  = int(sys.argv[3])

         tm = pyg.standardtimen('1 Jan 2000', 1)
         dt = tm.val_as_date(tm.date_as_val(dict(year=year, day=day, month=1)))

         #run_rrtm(dt['year'], dt['month'], dt['day'])
         run_pyracc(dt['year'], dt['month'], dt['day'])

      else:
         raise ValueException('Unknown command %s' % cmd)
