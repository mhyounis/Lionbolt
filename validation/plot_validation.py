
# A quick Terpdose-based script for plotting validation cases.
# Just to make it easy for the user (and myself) to create/re-create plots.
# Comparison is Lockwood experimental data (if present) vs. Lionbolt vs. DOSXYZnrc (mesh) / DOSRZnrc (slab)
# To run, just provide the folder which contains all the validation data you are interested in plotting.

import os
import numpy as np
import matplotlib
matplotlib.use ('Agg')
import matplotlib.pyplot as plt
from Terpdose import *
from argparse import ArgumentParser

#  =====================
#    Read in arguments  
#  =====================

parser = ArgumentParser ()
parser.add_argument ('folder', help="The validation case's folder")
args = parser.parse_args ()

input_dir = args.folder

colors = ['k', '#009CDE', '#E76F51' ]

#  ========================
#    Is this a slab case?  
#  ========================
# I'll determine this by seeing whether or not I can find DOSXYZnrc.3ddose in the folder

if os.path.isfile(f'{input_dir}/DOSXYZnrc.3ddose'):
    slab = False
else:
    slab = True

if slab:
    
    #  =================================
    #    Determine Lionbolt dose curve  
    #  =================================
    
    D = Lionbolt (f'{input_dir}/lb.h5')
    
    dmap = D.dose_deposition ()
    
    # fl = D.electrons.fluence()
    # # dmap = sum(fl)
    # for g in range(len(fl)):
    #     print(0.5 * (D.electrons.energy_grid[g] + D.electrons.energy_grid[g + 1]),
    #           fl[g][len(fl[0][:]) // 2] / (D.electrons.energy_grid[g] - D.electrons.energy_grid[g + 1]))
    # raise ValueError ('temp')
    
    fig, ax = plot_slab (
        D.mesh, 
        dmap, 
        # title  = r'521 keV electrons normally incident on Aluminum slab', 
        xlabel = r'Depth (cm)', 
        ylabel = r'Dose per fluence (MeV-cm$^{2}$/g)',
        color  = colors[1],
        pretty = False
    )
    
    #  ===========================
    #    Lockwood (if available)  
    #  ===========================
    
    L = f'{input_dir}/Lockwood.txt'
    try:
        with open(L, 'r') as f:
            # Nm = np.fromstring(f.readline(), sep=' ')    # Number of materials (for mean range reading)
            MR = np.fromstring(f.readline(), sep=' ')[0] # First line is just the mean range
            zL = np.fromstring(f.readline(), sep=' ')    # This is read in as fraction of a mean range, scaled later
            dL = np.fromstring(f.readline(), sep=' ')    # Dose data
        Lockwoodavail = True
    except:
        Lockwoodavail = False
    
    if Lockwoodavail:
        zL = zL * MR
        ax.plot (zL, dL, '-o', markersize=2.0, linewidth=0, markerfacecolor='none', color=colors[0])
    
    #  ============================
    #    Determine EGS dose curve  
    #  ============================
    
    def ReadDOSRZnrc(fname):
        with open(fname, 'r') as f:
            # for _ in range(632):
            for _ in range(16):
                next(f)
            
            z   = []
            EGS = []
            
            done = False
            while not done:
                line = f.readline()
                done = line[0] == '&'
                
                if not done:
                    zi, EGSi, _ = map(float, line.split())
                    
                    z.append(zi)
                    EGS.append(EGSi)
        
        return np.array(z), np.array(EGS)
    
    # Load EGS file
    z, EGS = ReadDOSRZnrc(f'{input_dir}/DOSRZnrc.plotdat')
    EGS = EGS * 6.242e+9
    
    ax.plot (z, EGS, '-', markersize=0.0, linewidth=0.75, color=colors[2])
    
    # Add legend
    if Lockwoodavail:
        # (Reorder so Lockwood is first)
        lines = ax.get_lines()
        ax.legend ([lines[1], lines[0], lines[2]],
                   [r'Lockwood (Exp.)', r'Lionbolt', r'EGSnrc'])
    else:
        ax.legend ([r'Lionbolt', r'EGSnrc'])
    
    # Make pretty using Terpdose's pretty_fig routine
    pretty_fig (fig, ax)
    
    # Save
    fig.savefig(f'{input_dir}/comparison.png', transparent=False, format='png', bbox_inches='tight', dpi=600)

else:
    
    #  ==============================================================================
    #    The below is highly specific to the water tank phantom test cases.
    #    Specifically, field size (via incident fluence factor),
    #    geo objects (Line, Plane), and plot shapes.
    #    At the moment the thinking is that validation of physics was largely
    #    accomplished in the slab case, so the general case doesn't need to be
    #    as extensive, just needs to show a correct X-ray->electron calculation to 
    #    validate source, PI-EI, sweeping (which was also independently validated),
    #    etc. The other general case solves were not compared to DOSXYZnrc.
    #  ==============================================================================
    
    import matplotlib.tri as tri
    
    #  ==================================
    #    Determine Lionbolt dose curves  
    #  ==================================
    
    D = Lionbolt (f'{input_dir}/lb.h5')
    
    dmap  = D.dose_deposition ()
    nodes = D.mesh.nodes ()
    
    # First want to plot YZ slice
    origin = [ HALF * (np.max(nodes[:,0]) - np.min(nodes[:,0])) + np.min(nodes[:,0]),
               HALF * (np.max(nodes[:,1]) - np.min(nodes[:,1])) + np.min(nodes[:,1]),
               HALF * (np.max(nodes[:,2]) - np.min(nodes[:,2])) + np.min(nodes[:,2]) ]
    
    sides  = [ np.max(nodes[:,1]) - np.min(nodes[:,1]),
                           np.max(nodes[:,2]) - np.min(nodes[:,2]) ]
    ax1 = np.array([0.0, 1.0, 0.0]) * sides[0] / TWO
    ax2 = np.array([0.0, 0.0, 1.0]) * sides[1] / TWO
    
    geo = Plane ( n=[500, 500], origin=origin, ax1=ax1, ax2=ax2 )
    fig2D1, ax2D1 = plot_2D (D.mesh, geo, dmap,
                             title   = r'Lionbolt',
                             xlabel  = r'$y$ (cm)',
                             ylabel  = r'$z$ (cm)',
                             cblabel = r'Dose per fluence (MeV-cm$^{2}$/g)',
                             prune   = True)
    
    # Then want to plot a depth profile
    # By convention we're taking this to be high -> low on the central axis
    x0 = np.array([ ZERO, ZERO, np.max(nodes[:,2]) ])
    x1 = np.array([ ZERO, ZERO, np.min(nodes[:,2]) ])
    
    geo = Line ( n=100, x0=x0, x1=x1 )
    fig1D, ax1D = plot_1D (
        D.mesh, 
        geo, 
        dmap,
        color  = colors[1],
        xlabel = r'Depth (cm)', 
        ylabel = r'Dose per fluence (MeV-cm$^{2}$/g)',
        pretty = False
    )
    
    #  =============================
    #    Determine EGS dose curves  
    #  =============================
    
    def ReadDOSXYZnrc(fname):
        with open(fname, 'r') as f:
            # Get the voxel numbers
            nx, ny, nz = map(int, f.readline().split())
            nvox = nx * ny * nz
            
            # Get the voxel boundaries
            xb = np.fromstring(f.readline(), sep=' ')
            yb = np.fromstring(f.readline(), sep=' ')
            zb = np.fromstring(f.readline(), sep=' ')
            
            # Plot points are voxel centers
            x = 0.5 * (xb[:-1] + xb[1:])
            y = 0.5 * (yb[:-1] + yb[1:])
            z = 0.5 * (zb[:-1] + zb[1:])
            
            # Read dose
            data = np.fromstring(f.read(), sep=' ')
            
            dose  = data[0:nvox]
            # error = data[nvox : 2 * nvox] # Dont really need error honestly
            
            dose = dose.reshape((nx, ny, nz), order='F')
            # error = error.reshape((nx, ny, nz), order='F')
        
        return x, y, z, dose
    
    # Load EGS file
    x, y, z, dose = ReadDOSXYZnrc(f'{input_dir}/DOSXYZnrc.3ddose')
    dose = dose * 6.242e+9
    
    nx = len(x)
    ny = len(y)
    nz = len(z)
    
    # For 2D plot, just plot it side by side.
    mid_x_idx = nx // 2
    yz_slice = dose[mid_x_idx, :, :]
    extent = (y[0], y[-1], z[-1], z[0])
    
    # Plot the EGS 2D Heatmap
    fig2D2, ax2D2 = plt.subplots(figsize=(7.1, 4))
    
    im2 = ax2D2.imshow(yz_slice.T, extent=extent, cmap='jet', aspect='equal', interpolation='nearest')
    fig2D2.colorbar (im2, ax=ax2D2, label=r'Dose per fluence (MeV-cm$^{2}$/g)')
    ax2D2.set_title (r'EGSnrc', fontsize=14)
    ax2D2.set_xlabel (r'$y$ (cm)', fontsize=14)
    ax2D2.set_ylabel (r'$z$ (cm)', fontsize=14)
    
    # For 1D plot, get the central axis values.
    # If nx is even, then average between the two voxels that hug the central axis.
    # Note, this is not an integration. This is just a numerical average
    if nx % 2 == 0:
        tmp = HALF * (dose[nx // 2 - 1,:,:] + dose[nx // 2,:,:])
    else:
        # Just take the middle value
        tmp = dose[nx // 2,:,:]
    # Repeat for ny
    if ny % 2 == 0:
        EGS1D = HALF * (tmp[ny // 2 - 1,:] + tmp[ny // 2,:])
    else:
        EGS1D = tmp[ny // 2,:]
    
    # Add EGS data to plot in the form of bins
    edges = np.empty(len(z) + 1)
    edges[1:-1] = (z[:-1] + z[1:]) / 2
    edges[0] = z[0] - (z[1] - z[0]) / 2
    edges[-1] = z[-1] + (z[-1] - z[-2]) / 2
    
    ax1D.stairs (EGS1D, edges, linewidth=0.75, color=colors[2])
    
    # Add legend
    ax1D.legend ([r'Lionbolt', r'EGSnrc'])
    
    # Make pretty using Terpdose's pretty_fig routine
    pretty_fig (fig1D, ax1D)
    
    # Save
    fig1D.savefig(f'{input_dir}/depth_comparison.png',        transparent=False, format='png', bbox_inches='tight', dpi=600)
    fig2D1.savefig(f'{input_dir}/Lionbolt_YZ.png', transparent=False, format='png', bbox_inches='tight', dpi=600)
    fig2D2.savefig(f'{input_dir}/EGSnrc_YZ.png',   transparent=False, format='png', bbox_inches='tight', dpi=600)
