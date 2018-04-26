import datetime


qtdbprefix = 'sel'

source_details ={
    'sel100': {
        'srcdb': 'mitdb',
        'srcstart': datetime.timedelta(minutes=7),
        'pattern': '(p)(N)t)'
    },
    'sel102': {
        'srcdb': 'mitdb',
        'srcstart': datetime.timedelta(minutes=6),
        'pattern': '(N)t)u)'
    },
    'sel103': {
        'srcdb': 'mitdb',
        'srcstart': datetime.timedelta(minutes=0),
        'pattern': '(p)(N)t)u)'
    },
    'sel104': {
        'srcdb': 'mitdb',
        'srcstart': datetime.timedelta(minutes=15),
        'pattern': '(N)t)'
    },
    'sel114': {
        'srcdb': 'mitdb',
        'srcstart': datetime.timedelta(minutes=0),
        'pattern': '(p)(N)t)'
    },
    'sel116': {
        'srcdb': 'mitdb',
        'srcstart': datetime.timedelta(minutes=0),
        'pattern': '(p)(N)t)u)'
    },
    'sel117': {
        'srcdb': 'mitdb',
        'srcstart': datetime.timedelta(minutes=15),
        'pattern': '(p)(N)t)u)'
    },
    'sel123': {
        'srcdb': 'mitdb',
        'srcstart': datetime.timedelta(minutes=3),
        'pattern': '(p)(N)t)u)'
    },
    'sel213': {
        'srcdb': 'mitdb',
        'srcstart': datetime.timedelta(minutes=0),
        'pattern': '(p)(N)t)'
    },
    'sel221': {
        'srcdb': 'mitdb',
        'srcstart': datetime.timedelta(minutes=0),
        'pattern': '(N)t)'
    },
    'sel223': {
        'srcdb': 'mitdb',
        'srcstart': datetime.timedelta(minutes=13),         #note: actual offset is different from laguna paper
        'pattern': '(p)(N)(t)'
    },
    'sel230': {
        'srcdb': 'mitdb',
        'srcstart': datetime.timedelta(minutes=7),
        'pattern': '(p)(N)t)'
    },
    'sel231': {
        'srcdb': 'mitdb',
        'srcstart': datetime.timedelta(minutes=3),
        'pattern': '(p)(N)t)'
    },
    'sel232': {
        'srcdb': 'mitdb',
        'srcstart': datetime.timedelta(minutes=15),
        'pattern': '(A)t)'
    },
    'sel233': {
        'srcdb': 'mitdb',
        'srcstart': datetime.timedelta(minutes=15),
        'pattern': '(p)(N)t)'
    },
    'sel301': {
        'srcdb': 'stdb',
        'srcstart': datetime.timedelta(minutes=14),
        'pattern': '(p)(N)(t)'
    },
    'sel302': {
        'srcdb': 'stdb',
        'srcstart': datetime.timedelta(minutes=8),
        'pattern': '(p)(N)t)'
    },
    'sel306': {
        'srcdb': 'stdb',
        'srcstart': datetime.timedelta(minutes=0),
        'pattern': '(p)(N)t)'
    },
    'sel307': {
        'srcdb': 'stdb',
        'srcstart': datetime.timedelta(minutes=0),
        'pattern': '(p)(N)t)'
    },
    'sel308': {
        'srcdb': 'stdb',
        'srcstart': datetime.timedelta(minutes=13),
        'pattern': '(p)(N)t)'
    },
    'sel310': {
        'srcdb': 'stdb',
        'srcstart': datetime.timedelta(minutes=4),
        'pattern': '(N)t)u)'
    },
   'sel803': {
        'srcdb': 'svdb',
        'srcstart': datetime.timedelta(minutes=0),
        'pattern': '(p)(N)t)'
    },
   'sel808': {
        'srcdb': 'svdb',
        'srcstart': datetime.timedelta(minutes=0),
        'pattern': '(p)(N)t)u)'
    },
   'sel811': {
        'srcdb': 'svdb',
        'srcstart': datetime.timedelta(minutes=7,seconds=30),
        'pattern': '(p)(N)t)'
    },
   'sel820': {
        'srcdb': 'svdb',
        'srcstart': datetime.timedelta(minutes=0),
        'pattern': '(p)(N)t)'
    },
   'sel821': {
        'srcdb': 'svdb',
        'srcstart': datetime.timedelta(minutes=0),
        'pattern': '(p)(N)(t)'
    },
   'sel840': {
        'srcdb': 'svdb',
        'srcstart': datetime.timedelta(minutes=14,seconds=50),
        'pattern': '(p)(N)t)'
    },
   'sel847': {
        'srcdb': 'svdb',
        'srcstart': datetime.timedelta(minutes=14,seconds=22),
        'pattern': '(p)(N)t)'
    },
   'sel853': {
        'srcdb': 'svdb',
        'srcstart': datetime.timedelta(minutes=5),
        'pattern': '(p)(N)t)'
    },
   'sel871': {
        'srcdb': 'svdb',
        'srcstart': datetime.timedelta(minutes=4),
        'pattern': '(p)(N)t)'
    },
   'sel872': {
        'srcdb': 'svdb',
        'srcstart': datetime.timedelta(minutes=0),
        'pattern': '(p)(N)t)'
    },
   'sel873': {
        'srcdb': 'svdb',
        'srcstart': datetime.timedelta(minutes=0),
        'pattern': '(p)(N)t)'
    },
   'sel883': {
        'srcdb': 'svdb',
        'srcstart': datetime.timedelta(minutes=13,seconds=30),
        'pattern': '(p)(N)t)'
    },
   'sel891': {
        'srcdb': 'svdb',
        'srcstart': datetime.timedelta(minutes=14,seconds=59),
        'pattern': '(p)(N)t)'
    },
   'sel16265': {
        'srcdb': 'nsrdb',
        'srcstart': datetime.timedelta(hours=10,minutes=37,seconds=30),
        'pattern': '(p)(N)t)'
    },
   'sel16272': {
        'srcdb': 'nsrdb',
        'srcstart': datetime.timedelta(hours=11,minutes=13,seconds=10),
        'pattern': '(p)(N)t)'
    },
   'sel16273': {
        'srcdb': 'nsrdb',
        'srcstart': datetime.timedelta(hours=11,minutes=22,seconds=36.5),
        'pattern': '(p)(N)t)'
    },
   'sel16420': {
        'srcdb': 'nsrdb',
        'srcstart': datetime.timedelta(hours=14,minutes=24,seconds=0),
        'pattern': '(p)(N)t)'
    },
   'sel16483': {
        'srcdb': 'nsrdb',
        'srcstart': datetime.timedelta(hours=13,minutes=32,seconds=0),
        'pattern': '(p)(N)t)'
    },
   'sel16539': {
        'srcdb': 'nsrdb',
        'srcstart': datetime.timedelta(hours=20,minutes=56,seconds=0),
        'pattern': '(p)(N)(t)'
    },
   'sel16773': {
        'srcdb': 'nsrdb',
        'srcstart': datetime.timedelta(hours=10,minutes=23,seconds=40),
        'pattern': '(p)(N)t)u)'
    },
   'sel16786': {
        'srcdb': 'nsrdb',
        'srcstart': datetime.timedelta(hours=18,minutes=22,seconds=0),
        'pattern': '(p)(N)(t)'
    },
   'sel16795': {
        'srcdb': 'nsrdb',
        'srcstart': datetime.timedelta(hours=15,minutes=0,seconds=0),
        'pattern': '(p)(N)(t)'
    },
   'sel17453': {
        'srcdb': 'nsrdb',
        'srcstart': datetime.timedelta(hours=14,minutes=15,seconds=0),
        'pattern': '(p)(N)(t)u)'
    },
   'sele0104': {
        'srcdb': 'edb',
        'srcstart': datetime.timedelta(hours=1,minutes=35,seconds=0),
        'pattern': '(p)(N)t)'
    },
   'sele0106': {
        'srcdb': 'edb',
        'srcstart': datetime.timedelta(hours=1,minutes=31,seconds=0),
        'pattern': '(p)(N)(t)u)'
    },
   'sele0107': {
        'srcdb': 'edb',
        'srcstart': datetime.timedelta(hours=1,minutes=19,seconds=0),
        'pattern': '(p)(N)t)'
    },
   'sele0110': {
        'srcdb': 'edb',
        'srcstart': datetime.timedelta(hours=1,minutes=32,seconds=0),
        'pattern': '(p)(N)t)'
    },
   'sele0111': {
        'srcdb': 'edb',
        'srcstart': datetime.timedelta(hours=1,minutes=19,seconds=0),
        'pattern': '(p)(N)t)u)'
    },
   'sele0112': {
        'srcdb': 'edb',
        'srcstart': datetime.timedelta(hours=1,minutes=32,seconds=10),
        'pattern': '(p)(N)t)u)'
    },
   'sele0114': {
        'srcdb': 'edb',
        'srcstart': datetime.timedelta(hours=1,minutes=38,seconds=0),
        'pattern': '(p)(N)(t)(u)'
    },
   'sele0116': {
        'srcdb': 'edb',
        'srcstart': datetime.timedelta(hours=0,minutes=37,seconds=30),
        'pattern': '(p)(N)(t)'
    },
   'sele0121': {
        'srcdb': 'edb',
        'srcstart': datetime.timedelta(hours=1,minutes=7,seconds=30),
        'pattern': '(p)(N)t)'
    },
   'sele0122': {
        'srcdb': 'edb',
        'srcstart': datetime.timedelta(hours=0,minutes=40,seconds=0),
        'pattern': '(p)(N)t)u)'
    },
   'sele0124': {
        'srcdb': 'edb',
        'srcstart': datetime.timedelta(hours=1,minutes=41,seconds=0),
        'pattern': '(p)(N)(t)'
    },
   'sele0126': {
        'srcdb': 'edb',
        'srcstart': datetime.timedelta(hours=1,minutes=36,seconds=40),
        'pattern': '(p)(N)t)u)'
    },
   'sele0129': {
        'srcdb': 'edb',
        'srcstart': datetime.timedelta(hours=0,minutes=26,seconds=30),
        'pattern': '(p)(N)t)'
    },
   'sele0133': {
        'srcdb': 'edb',
        'srcstart': datetime.timedelta(hours=0,minutes=37,seconds=30),
        'pattern': '(p)(N)t)'
    },
   'sele0136': {
        'srcdb': 'edb',
        'srcstart': datetime.timedelta(hours=0,minutes=56,seconds=0),
        'pattern': '(p)(N)(t)'
    },
   'sele0166': {
        'srcdb': 'edb',
        'srcstart': datetime.timedelta(hours=0,minutes=52,seconds=30),
        'pattern': '(p)(N)(t)'
    },
   'sele0170': {
        'srcdb': 'edb',
        'srcstart': datetime.timedelta(hours=0,minutes=22,seconds=30),
        'pattern': '(p)(N)t)'
    },
   'sele0203': {
        'srcdb': 'edb',
        'srcstart': datetime.timedelta(hours=0,minutes=52,seconds=30),
        'pattern': '(p)(N)t)'
    },
   'sele0210': {
        'srcdb': 'edb',
        'srcstart': datetime.timedelta(hours=1,minutes=45,seconds=0),
        'pattern': '(p)(N)(t)'
    },
   'sele0211': {
        'srcdb': 'edb',
        'srcstart': datetime.timedelta(hours=1,minutes=45,seconds=0),
        'pattern': '(p)(N)(t)'
    },
   'sele0303': {
        'srcdb': 'edb',
        'srcstart': datetime.timedelta(hours=0,minutes=34,seconds=30),
        'pattern': '(p)(N)(t)'
    },
   'sele0405': {
        'srcdb': 'edb',
        'srcstart': datetime.timedelta(hours=0,minutes=22,seconds=30),
        'pattern': '(p)(N)t)'
    },
   'sele0406': {
        'srcdb': 'edb',
        'srcstart': datetime.timedelta(hours=1,minutes=0,seconds=0),
        'pattern': '(p)(N)t)'
    },
   'sele0409': {
        'srcdb': 'edb',
        'srcstart': datetime.timedelta(hours=1,minutes=20,seconds=0),
        'pattern': '(p)(N)(t)'
    },
   'sele0411': {
        'srcdb': 'edb',
        'srcstart': datetime.timedelta(hours=1,minutes=30,seconds=0),
        'pattern': '(p)(N)t)'
    },
   'sele0509': {
        'srcdb': 'edb',
        'srcstart': datetime.timedelta(hours=0,minutes=15,seconds=0),
        'pattern': '(p)(N)t)'
    },
   'sele0603': {
        'srcdb': 'edb',
        'srcstart': datetime.timedelta(hours=0,minutes=17,seconds=0),
        'pattern': '(p)(N)t)u)'
    },
   'sele0604': {
        'srcdb': 'edb',
        'srcstart': datetime.timedelta(hours=0,minutes=0,seconds=30),
        'pattern': '(p)(N)t)u)'
    },
   'sele0606': {
        'srcdb': 'edb',
        'srcstart': datetime.timedelta(hours=1,minutes=31,seconds=0),
        'pattern': '(p)(N)t)'
    },
   'sele0607': {
        'srcdb': 'edb',
        'srcstart': datetime.timedelta(hours=1,minutes=0,seconds=10),
        'pattern': '(p)(N)t)'
    },
   'sele0609': {
        'srcdb': 'edb',
        'srcstart': datetime.timedelta(hours=0,minutes=35,seconds=30),
        'pattern': '(p)(N)(t)'
    },
   'sele0612': {
        'srcdb': 'edb',
        'srcstart': datetime.timedelta(hours=1,minutes=35,seconds=0),
        'pattern': '(p)(N)(t)'
    },
   'sele0704': {
        'srcdb': 'edb',
        'srcstart': datetime.timedelta(hours=1,minutes=16,seconds=30),
        'pattern': '(p)(N)t)u)'
    },
   'sel30': {
        'srcdb': 'sddb',
        'srcstart': datetime.timedelta(hours=7,minutes=39,seconds=30),
        'pattern': '(p)(N)(t)'
    },
   'sel31': {
        'srcdb': 'sddb',
        'srcstart': datetime.timedelta(hours=12,minutes=57,seconds=30),
        'pattern': '(p)(N)(t)'
    },
   'sel32': {
        'srcdb': 'sddb',
        'srcstart': datetime.timedelta(hours=20,minutes=52,seconds=17.68),
        'pattern': '(p)(N)(t)'
    },
   'sel33': {
        'srcdb': 'sddb',
        'srcstart': datetime.timedelta(hours=14,minutes=57,seconds=40),
        'pattern': '(p)(N)(t)'
    },
   'sel34': {
        'srcdb': 'sddb',
        'srcstart': datetime.timedelta(hours=5,minutes=53,seconds=0),
        'pattern': '(p)(N)(t)'
    },
   'sel35': {
        'srcdb': 'sddb',
        'srcstart': datetime.timedelta(hours=4,minutes=52,seconds=0),
        'pattern': '(N)'
    },
   'sel36': {
        'srcdb': 'sddb',
        'srcstart': datetime.timedelta(hours=18,minutes=43,seconds=0),
        'pattern': '(B)(t)u)'
    },
   'sel37': {
        'srcdb': 'sddb',
        'srcstart': datetime.timedelta(hours=1,minutes=14,seconds=10),
        'pattern': '(Q)'
    },
   'sel38': {
        'srcdb': 'sddb',
        'srcstart': datetime.timedelta(hours=5,minutes=0,seconds=0),
        'pattern': '(p)(N)(t)'
    },
   'sel39': {
        'srcdb': 'sddb',
        'srcstart': datetime.timedelta(hours=0,minutes=4,seconds=30),
        'pattern': '(p)(N)(t)'
    },
   'sel40': {
        'srcdb': 'sddb',
        'srcstart': datetime.timedelta(hours=7,minutes=7,seconds=0),
        'pattern': '(p)(N)(t)'
    },
   'sel41': {
        'srcdb': 'sddb',
        'srcstart': datetime.timedelta(hours=0,minutes=5,seconds=0),
        'pattern': '(p)(N)(t)'
    },
   'sel42': {
        'srcdb': 'sddb',
        'srcstart': datetime.timedelta(hours=10,minutes=3,seconds=0),
        'pattern': '(p)(N)(t)'
    },
   'sel43': {
        'srcdb': 'sddb',
        'srcstart': datetime.timedelta(hours=8,minutes=0,seconds=0),
        'pattern': '(p)(N)(t)'
    },
   'sel44': {
        'srcdb': 'sddb',
        'srcstart': datetime.timedelta(hours=21,minutes=30,seconds=0),
        'pattern': '(p)(N)(t)'
    },
   'sel45': {
        'srcdb': 'sddb',
        'srcstart': datetime.timedelta(hours=17,minutes=30,seconds=30),
        'pattern': '(p)(N)(t)u)'
    },
   'sel46': {
        'srcdb': 'sddb',
        'srcstart': datetime.timedelta(hours=3,minutes=10,seconds=0),
        'pattern': '(p)(N)(t)'
    },
   'sel47': {
        'srcdb': 'sddb',
        'srcstart': datetime.timedelta(hours=8,minutes=43,seconds=40),
        'pattern': '(p)(N)(t)(u)'
    },
   'sel48': {
        'srcdb': 'sddb',
        'srcstart': datetime.timedelta(hours=7,minutes=11,seconds=0),
        'pattern': '(p)(N)(t)'
    },
   'sel49': {
        'srcdb': 'sddb',
        'srcstart': datetime.timedelta(hours=23,minutes=10,seconds=40),
        'pattern': '(p)(N)(t)'
    },
   'sel50': {
        'srcdb': 'sddb',
        'srcstart': datetime.timedelta(hours=10,minutes=53,seconds=00),
        'pattern': '(N)(t)u)'
    },
   'sel50': {
        'srcdb': 'sddb',
        'srcstart': datetime.timedelta(hours=10,minutes=53,seconds=00),
        'pattern': '(N)(t)u)'
    },
   'sel51': {
        'srcdb': 'sddb',
        'srcstart': datetime.timedelta(hours=0,minutes=55,seconds=00),
        'pattern': '(p)(N)(t)'
    },
   'sel52': {
        'srcdb': 'sddb',
        'srcstart': datetime.timedelta(hours=0,minutes=8,seconds=00),
        'pattern': '(p)(N)(t)'
    },
   'sel14046': {
        'srcdb': 'ltdb',
        'srcstart': datetime.timedelta(hours=9,minutes=14,seconds=13),
        'pattern': '(p)(N)t)'
    },
   'sel14157': {
        'srcdb': 'ltdb',
        'srcstart': datetime.timedelta(hours=7,minutes=58,seconds=12),
        'pattern': '(p)(N)t)u)'
    },
   'sel14172': {
        'srcdb': 'ltdb',
        'srcstart': datetime.timedelta(hours=9,minutes=14,seconds=13),
        'pattern': '(p)(N)(t)(u)'
    },
   'sel15814': {
        'srcdb': 'ltdb',
        'srcstart': datetime.timedelta(hours=21,minutes=57,seconds=0),
        'pattern': '(p)(N)t)'
    },
}





