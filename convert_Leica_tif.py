#!/usr/bin/env python

from __future__ import print_function
import sys
import os
import re
import glob
import shlex, subprocess

if __name__=='__main__':
  helpmess="""Usage:
convert_Leica_tif dir

Converts the output of the Leica SB configuration into single tiffs.
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
        cmd='convert {files} -evaluate-sequence Mean -type Grayscale -auto-level {outfile}'.format(files=' '.join(group), outfile=outfile)
        print(cmd)

        args = shlex.split(cmd)
        subprocess.call(args)

  #cmd='Trinity --seqType {ftype} {data}--max_memory 8G --CPU 4 --trimmomatic ' \
 #     '--verbose --FORCE_INCHWORM_KMER_METHOD --quality_trimming_params "ILLUMINACLIP:{trimmpath}adapters/TruSeq3-PE.fa:2:30:12:1:true MAXINFO:40:0.4 LEADING:3 TRAILING:3 MINLEN:40" ' \
 #     '--output {outfile}_trinity --full_cleanup'.format(ftype=ftype, data=data, trimmpath=trimmpath, outfile=outfile)
#
#  if ftype!=None and data!=None:
#    print('\n*******************************************')
#    print(cmd)
#    print('*******************************************\n')
#    args = shlex.split(cmd)
#    subprocess.call(args)
