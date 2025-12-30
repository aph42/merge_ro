from . import process as pr
import gc

import sys

Nzs = dict(champ2016  = 3500,
           cnofs      = 6500,
           cnofsrt    = 6500,
           cosmic     = 3500,
           cosmic2013 = 3500,
           grace      = 7800,
           kompsat5   = 6000,
           paz        = 6000,
           metopa     = 3500,
           metopb     = 3500,
           metopa2016 = 3500,
           metopb2016 = 3500,
           sacc       = 7000,
           saccrt     = 7000,
           tsx        = 7800)

missions = [
            'champ2016'
            #'cnofs' 
            #'cosmic2013' 
            #'metopa2016', 
            #'grace',
            #'sacc' 
            #'tsx'
            ]



m = sys.argv[1]
year = int(sys.argv[2])
d = int(sys.argv[3])

assert m in list(Nzs.keys())

print('Processing %s %d-%03d' % (m, year, d))

pr.process_date(m, year, d, Nz=Nzs[m])
