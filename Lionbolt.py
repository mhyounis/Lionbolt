#!/usr/bin/env python3
import sys
import os
import glob
import subprocess
import multiprocessing
import argparse
import textwrap
from pathlib import Path

ScriptDir = os.path.dirname(os.path.abspath(__file__))

# MHY LATER - PUT A CLEANUP FUNCTION THAT CLEANS TEMPORARY FILES AND ALL THAT.

def Parser ():
    
    parser = argparse.ArgumentParser(description=textwrap.dedent('''
        ====================================================================================
                                 This is the Lionbolt submit script.                        
          It handles building, cleaning, and running Lionbolt in either the standard input  
          file mode or the driver program mode.                                             
                                                                                            
          It is by no means necessary to use Lionbolt, however, it makes life a lot easier  
          when running simple jobs and developing new code.                                 
        ====================================================================================
        
        If only -c is specified, then Lionbolt is cleaned, and nothing is built.
        If only -b is specified, then Lionbolt is built.
        If an input file is provided, then a job is run after the other commands are executed. If -c is specified here, then after cleaning, Lionbolt is built, and then run.
        '''), formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('Input file(s)', help='Input file(s) (.in) or one driver program (.f90)', nargs='*')
    parser.add_argument('-n', type=int, default=multiprocessing.cpu_count(), help='Number of threads you wish to make available to Lionbolt. Default is all.')
    parser.add_argument('--scratch', type=str, default=None, help='Custom scratch directory to use. You can also set the scratch directory permanently using environment variable $OPENRPS_SCRATCH')
    builder = parser.add_mutually_exclusive_group(required=False)
    builder.add_argument('-c', action='store_true', help="Clean Lionbolt's build and temporary directories")
    builder.add_argument('-b', action='store_true', help='Build Lionbolt.')
    builder.add_argument('-m', action='store_true', help='Build Lionbolt with developer compile options.')
    
    args = parser.parse_args()
    
    inputs = getattr(args, 'Input file(s)')
    
    # If no inputs given, manually assign None, as nargs='*' does not allow default value
    if len(inputs) == 0:
        inputs = None
    
    return inputs, args.n, args.scratch, args.c, args.b, args.m

def Verify (inputs):
    
    if inputs is None:
        return 0
    
    # Check if all files exist
    for i in inputs:
        path = Path (i)
        if not path.is_file():
            print ('  Lionbolt.py ---')
            print ('  ERROR : A provided file does not exist.')
            print ('          The file is : ' + i)
            sys.exit(1)
    
    if len(inputs) > 1:
        # If there are multiple inputs check their extensions.
        # Can only have one driver file, but arbitrarily many
        # input files
        for i in inputs:
            parts = i.split('.')
            
            # Extension is the last part
            ext = parts[-1]
            
            if ext == 'f90':
                print ('  Lionbolt.py ---')
                print ('  ERROR : Multiple files were provided, at least one of which is a driver.')
                print ('          The file is : ' + i)
                print ('          You can provide multiple input files but only one driver at a time.')
                sys.exit(1)
            elif ext != 'in':
                print ('  Lionbolt.py ---')
                print ('  ERROR : A provided file is neither an input (.in) nor a driver (.f90).')
                print ('          The file is : ' + i)
                sys.exit(1)
        
        return 1
    else:
        parts = inputs[0].split('.')
        ext = parts[-1]
        if ext == 'in':
            return 1
        elif ext =='f90':
            return 2

def Clean ():
    
    subprocess.call(['rm', '-rf', 'build'], cwd=ScriptDir)
    
    print ('    ')
    print ('    ===========================================')
    print ('      Lionbolt has been cleaned successfully.  ')
    print ('    ===========================================')
    print ('    ')

def Build (dev=False, driverfile=None):
    
    if driverfile is None:
        if dev:
            rc = subprocess.call(['cmake', '-S', '.', '-B', 'build', '-DENABLE_DEV_COMPILE=ON'], cwd=ScriptDir)
        else:
            rc = subprocess.call(['cmake', '-S', '.', '-B', 'build'], cwd=ScriptDir)
    else:
        # Build using full driver path
        abspath = os.path.abspath(driverfile)
        
        if dev:
            rc = subprocess.call(['cmake', '-S', '.', '-B', 'build', '-DLBDRIVER=' + abspath, '-DENABLE_DEV_COMPILE=ON'], cwd=ScriptDir)
        else:
            rc = subprocess.call(['cmake', '-S', '.', '-B', 'build', '-DLBDRIVER=' + abspath], cwd=ScriptDir)
    
    if rc != 0:
        sys.exit(1)
    
    rc = subprocess.call(['cmake', '--build', 'build', '-j'], cwd=ScriptDir)
    if rc != 0:
        sys.exit(1)
    
    print ('    ')
    print ('    =========================================')
    print ('      Lionbolt has been built successfully.  ')
    print ('    =========================================')
    print ('    ')

def Run (inp, nthreads, scratch, isinputfile):
    
    abspath = os.path.abspath(inp)
    inpdir  = os.path.dirname(abspath)
    inpfile = os.path.basename(abspath)
    
    env = os.environ.copy()
    env['OMP_NUM_THREADS'] = str(nthreads)
    env['GFORTRAN_UNBUFFERED_ALL'] = 'y'
    
    if scratch is None:
        # This tries to use $OPENRPS_SCRATCH, but if that isn't defined, it tries TMPDIR, then if that isn't defined, it uses /tmp
        scratch = (os.environ.get('OPENRPS_SCRATCH') or os.environ.get('TMPDIR') or '/tmp')
    
    env['OPENRPS_SCRATCH'] = scratch
    
    LionboltExec = os.path.join(ScriptDir, 'build', 'Lionbolt')
    
    if isinputfile:
        subprocess.run([LionboltExec, inpfile], cwd=inpdir, env=env)
    else:
        subprocess.run([LionboltExec], cwd=inpdir, env=env)

##################################################################################################

try:
    inputs, nthreads, scratch, clean, make, dev = Parser ()
    
    # First check if the inputs are input files or a driver file.
    info = Verify (inputs)
    
    # info = 0, no input was given, this is a pure build/clean
    # info = 1, input file(s)
    # info = 2, driver file
    
    if clean:
        Clean ()
    
    if info == 0:
        # Make if requested then exit
        if make or dev:
            Build (dev=dev, driverfile=None)
        
        # If nothing had been requested then state that this script did nothing
        if not (make or dev or clean):
            print ('  Lionbolt.py --- Did nothing.')
    
    elif info == 1:
        # Make if requested (or required) then loop through input files and run
        if not (make or dev):
            # Check if we NEED to make regardless
            # NOTE - this will not necessarily mean that the built Lionbolt ISN'T a driver.
            # This means users will have to be cognizant to clean a driver build before feeding an input file.
            # I could get around this by making the CMakeLists use a different executable name for the driver case
            # vs. the input file case.
            executable = Path (ScriptDir + '/build/Lionbolt')
            make = not executable.is_file()
        
        if make or dev:
            Build (dev=dev, driverfile=None)
        
        if len(inputs) > 1:
            print('  Lionbolt.py --- Submitting multiple input files.')
            for i in inputs:
                print(f'    {os.path.abspath(i)}')
        for i in inputs:
            Run (i, nthreads, scratch, isinputfile=True)
    
    elif info == 2:
        # Always must make
        Build (dev=dev, driverfile=inputs[0])
        
        Run (inputs[0], nthreads, scratch, isinputfile=False)

except KeyboardInterrupt:
    # Might later add an error message for keyboard interrupts
    pass
