#!/usr/bin/env python3

import sys

help = """
  sequence_extract.py busy|idle

  Given a sequence of start/end times on stdin like:
    3.9408 3.9430
    3.9447 3.9466
    3.9487 3.9506
  Output either the busy times (end-start)
  or idle times (this.start - prev.end)

"""

"""
Source of input:
grep 'X_POSIX[ 0]*write' hdf5_large_file_details.txt | head -954 | cutfast -w -f 6,7
grep 'X_POSIX[ 0]*write' meshio_large_file_details.txt | head -954 | cutfast -w -f 6,7
"""

if len(sys.argv) != 2:
  print(help)
  sys.exit(1)

mode = sys.argv[1]
if mode != 'busy' and mode != 'idle':
  print(help)
  sys.exit(1)

prev_end = None

for line in sys.stdin:
  fields = line.split()[:2]
  s = float(fields[0])
  e = float(fields[1])
  
  if mode == 'busy':
    print('%.6e' % (e-s))
  else:
    if prev_end != None:
      print('%.6e' % (s - prev_end))

  prev_end = e
