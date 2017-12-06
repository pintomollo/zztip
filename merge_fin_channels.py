#!/usr/bin/env python

from __future__ import print_function
import sys
import os
import re
import glob
import shlex, subprocess

if __name__=='__main__':
  helpmess="""Usage:
merge_fin_channels dir

Converts three independent channels into one RGB image.
"""
  # Inputs
  if len(sys.argv)<2:
      print(helpmess)
      sys.exit(0)
  else:
      indir=sys.argv[1]

  indir=os.path.realpath(indir)

  if os.path.exists(indir) and os.path.isdir(indir):
    for root, dirs, files in os.walk(indir):

      files=glob.glob('{mydir}'.format(mydir=os.path.join(root, '*_ch00.tif')))

      if len(files)>0:
        outdir=os.path.join(root, 'export')
        if not os.path.exists(outdir):
          os.mkdir(outdir)

      for infile in files:
        basedir,basename=os.path.split(infile)
        basename,ext=os.path.splitext(basename)
        parts=basename.split('_')
        parts[-1]='*'

        pattern=os.path.join(basedir, '{mypattern}{ext}'.format(mypattern='_'.join(parts),ext=ext))
        group=glob.glob(pattern)

        outfile=os.path.join(outdir, '{myfile}{ext}'.format(myfile='_'.join(parts[0:-1]),ext=ext))
        cmd='convert {file1} {file2} {file3} -normalize -quiet -combine {outfile}'.format(file1=group[1],file2=group[3],file3=group[2], outfile=outfile)
        print(cmd)

        args = shlex.split(cmd)
        subprocess.call(args)
