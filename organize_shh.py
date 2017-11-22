#!/usr/bin/env python

from __future__ import print_function
import sys
import os
import re
import glob
import shutil
import shlex, subprocess

if __name__=='__main__':
  helpmess="""Usage:
organize_shh dir

Converts the output of the Leica SB configuration into single tiffs.
"""
  # Inputs
  if len(sys.argv)<2:
      print(helpmess)
      sys.exit(0)
  else:
      indir=sys.argv[1]

  indir=os.path.realpath(indir)

  count = 1;
  dpa = '-100';
  isgfp = False;

  if os.path.exists(indir) and os.path.isdir(indir):
    files = os.listdir(indir)
    for fname in files:
      res = re.match('SB25_(\d+)_(-?\d+)dpa_(Z\d+)_(\d+)(.*)\.tif', fname)
      if res != None:
        ids = res.group(2)
        index = ''

        if ids != dpa:
          if (res.group(5) == '_GFP'):
            dpa = ids
            count = 1
        else:
          index = chr(96+count)
          if (res.group(5) == '_GFP'):
            count += 1

        print('{name} -> SB25_{z}{i}{g}_{d}dpa_{D}'.format(name=fname, z=res.group(3), i=index, g=res.group(5), d=ids, D=res.group(1)))
        os.rename(os.path.join(indir, fname), os.path.join(indir, 'SB25_{z}{i}{g}_{d}dpa_{D}.tif'.format(z=res.group(3), i=index, g=res.group(5), d=ids, D=res.group(1))))
