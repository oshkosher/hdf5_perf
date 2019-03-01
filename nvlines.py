#!/usr/bin/env python3

"""
parse lines like this:

hdf5_open_write np=128 stripe_count=1 stripe_len_mb=1 open_sec=0.093
hdf5_write np=128 size_mb=1024.000 stripe_count=1 stripe_len_mb=1 write_sec=4.692 write_mbps=218.260
hdf5_write np=128 size_mb=1024.000 stripe_count=1 stripe_len_mb=1 write_sec=3.914 write_mbps=261.631

Use the first word on the line as a key, to specify which lines to ignore.
name=value pairs are parsed from the rest of the words.

  nvlines [opts]
  opts:
    ^foo : only select lines with this prefix. Multiple prefixes can be
           specified, and any will match.
    foo=bar : select only lines where field foo is bar
    foo : select field foo
    -n : no header row

  Each selected field will be a column in the output.
  
"""

import sys

prefix_vals = set()
output_list = []
# map column name to 0-based index
output_pos = {}
# map column index to required value
required = {}

output_header = True


def main(args):
  parseArgs(args)
  """
  print('prefix_vals=' + repr(prefix_vals))
  print('output_list=' + repr(output_list))
  print('output_pos=' + repr(output_pos))
  print('required=' + repr(required))
  """

  if len(output_list) == 0:
    sys.stderr.write('No columns selected\n')
    return 0

  if output_header:
    print('\t'.join(output_list))

  empty_line = [''] * len(output_list)

  line_no = 0
  for line in sys.stdin:
    line_no += 1
    # print('line %d' % line_no)
    fields = line.rstrip().split(' ')
    if len(fields) == 0: continue

    # screen out non-matching prefixes
    if len(prefix_vals) > 0 and fields[0] not in prefix_vals:
      continue

    # print(repr(fields))
    
    out = empty_line.copy()

    # extract fields, assign them to out[]
    out_empty = True
    for field in fields[1:]:
      eq = field.find('=')
      if eq == -1: continue
      name = field[:eq]
      value = field[eq+1:]
      # print("  %s=%s, %s" % (name, value, repr(output_pos.get(name, -1))))
      pos = output_pos.get(name, -1)
      
      if pos != -1:
        out_empty = False
        out[pos] = value

    if out_empty: continue

    # check that all required values are set
    if len(required) > 0:
      required_missing = False
      for (pos,value) in required.items():
        if out[pos] != value:
          required_missing = True
          break

      if required_missing: continue

    # output the row
    print('\t'.join(out))
    


def parseArgs(args):
  global prefix_vals, output_list, output_pos, required, output_header

  for arg in args:
    
    if arg == '-n':
      output_header = False

    elif len(arg) > 0 and arg[0] == '^':
      prefix_vals.add(arg[1:])
    
    else:
      eq = arg.find('=')
      if eq != -1:
        name = arg[:eq]
        val = arg[eq+1:]
        pos = len(output_list)
        output_pos[name] = len(output_list)
        output_list.append(name)
        required[pos] = val

      else:
        output_pos[arg] = len(output_list)
        output_list.append(arg)


if __name__ == '__main__':
  sys.exit(main(sys.argv[1:]))
