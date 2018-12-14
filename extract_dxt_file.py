#!/usr/bin/env python3

import sys, re


def main(args):
  if len(args) > 2 or '-h' in args:
    sys.stderr.write("""
    extract_dxt_file.py [filename [rank]]
    This takes the output of darshan-dxt-parser as input.
    With no arguments, this outputs all the file names listed in
    the Darshan data.
    Given a filename as an argument, it extracts the X_POSIX events for
    that file.
    Given a rank, it only output data from that rank.

  """)
    return 1

  if len(args) == 0:
    extractFilenames()
    return 0

  filename = args[0]
  rank = None

  if len(args) > 1:
    rank = args[1]

  extractOneFile(filename, rank)

file_id_regex = re.compile('^# DXT, file_id: [0-9]+, file_name: (.*)$')
rank_regex = re.compile('^# DXT, rank: ([0-9]+), ')

def extractFilenames():

  filename_set = set()
  
  for line in sys.stdin:
    match = file_id_regex.search(line)
    if match:
      filename = match.group(1)
      if filename not in filename_set:
        filename_set.add(filename)
        print(filename)


def extractOneFile(filename, rank):
  current_file = None
  current_rank = None

  for line in sys.stdin:
    # print('line ' + line)
    match = file_id_regex.search(line)
    if match:
      current_file = match.group(1)
      current_rank = None
      # print('  current_file = ' + current_file)
      continue

    match = rank_regex.search(line)
    if match:
      current_rank = match.group(1)
      # print('  current_rank = ' + current_rank)
      continue

    # print('sw %s, %s, %s, %s' % 
    #       (repr(line.startswith(' X_POSIX')), repr(current_file == filename),
    #        rank, current_rank))
                                 
    if (line.startswith(' X_POSIX') and current_file == filename
        and (rank == None or rank == current_rank)):
      sys.stdout.write(line)
    

if __name__ == "__main__":
  sys.exit(main(sys.argv[1:]))
