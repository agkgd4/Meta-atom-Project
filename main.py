from cProfile import label
from lib2to3.pgen2.literals import simple_escapes
from pickle import NONE
from re import L
import meep as mp
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

def build_geometry(params, radius): 
    geometry = []
    for element in params["geometry"]:
        if params["geometry"][element]["block"] == True:
            element = params["geometry"][element]

            xdim = element["xdim"]
            ydim = element["ydim"]
            zdim = element["zdim"]

            xcen = element["xcen"]
            ycen = element["ycen"]
            zcen = element["zcen"]

            size = mp.Vector3(xdim,ydim,zdim)
            center = mp.Vector3(xcen,ycen,zcen)

            n = element["n"]
            material = mp.Medium(index=n)
            
            geometry.append(mp.Block(size=size,
                                center=center,
                                material=material))

        # for params["geometry"][element]["block"] == False (This is a cylinder)
        else:
            element = params["geometry"][element]
            height = element["height"]
            radius = radius

            xcen = element["xcen"]
            ycen = element["ycen"]
            zcen = element["zcen"]

            axis = mp.Vector3(0,1,0)

            center = mp.Vector3(xcen,ycen,zcen)

            n = element["n"]
            material = mp.Medium(index=n)

            geometry.append(mp.Cylinder(radius=radius,
                                        axis=axis,
                                        height=height,
                                        center=center,
                                        material=material))
            
    return geometry

def build_source(params):
    source = []
    for element in params["source"]:
        element = params["source"][element]
        freq = element["freq"]

        xcen = element["xcen"]
        ycen = element["ycen"]
        zcen = element["zcen"]
        
        xdim = element["xdim"]
        ydim = element["ydim"]
        zdim = element["zdim"]

        n = element["n"]

        source.append(mp.EigenModeSource(mp.ContinuousSource(frequency=freq,is_integrated=True),            # is_integrated = True for bloch-periodic boundary conditions for planar wavefront
                                    center=mp.Vector3(xcen,ycen,zcen),
                                    size=mp.Vector3(xdim,ydim,zdim),
                                    direction=mp.AUTOMATIC,                                                 # this will need to be changed if not normally incident wave
                                    eig_kpoint=mp.Vector3(freq*n).rotate(mp.Vector3(z=1),np.radians(0)),
                                    eig_parity=mp.EVEN_Y+mp.ODD_Z,                                          # this will need to be updated if not normally incident wave
                                    eig_match_freq=True))
                          
    return source

def build_sim(params,source,geometry=None):          
    params = params["cell_config"]
    cell_x = params["cell_x"]
    cell_y = params["cell_y"]
    cell_z = params["cell_z"]
    cell_size = mp.Vector3(cell_x,cell_y,cell_z)

    pml_thickness = params["pml_thickness"]
    pml_layers = [mp.PML(thickness=pml_thickness,direction=mp.Y)]

    resolution = params["resolution"]

    freq = params["freq"]
    n = params["n"]
    rot_angle = params["rot_angle"]

    #k_point=mp.Vector3(freq*n).rotate(mp.Vector3(z=1),rot_angle)    # planewaves in homogenous media: https://meep.readthedocs.io/en/latest/Python_Tutorials/Eigenmode_Source/ Also has PML only in x-direction and specifies k_point just like this.
    k_point=mp.Vector3(0,0,0)
    if geometry is None:
        sim = mp.Simulation(cell_size=cell_size,
                    boundary_layers=pml_layers,
                    sources=source,
                    k_point=k_point,
                    resolution=resolution)
    else:
        sim = mp.Simulation(cell_size=cell_size,
                    boundary_layers=pml_layers,
                    geometry=geometry,
                    sources=source,
                    k_point=k_point,
                    resolution=resolution)
    
    return sim

def build_flux_region(params,sim):
    params = params["flux_region"]
    xdim = params["xdim"]
    ydim = params["ydim"]
    zdim = params["zdim"]

    xcen = params["xcen"]
    ycen = params["ycen"]
    zcen = params["zcen"]

    fcen = params["fcen"]
    df = params["df"]
    nfreq = params["nfreq"]

    fr = mp.FluxRegion(center=mp.Vector3(xcen,ycen,zcen),
                                size=mp.Vector3(xdim,ydim,zdim))

    flux_region = sim.add_flux(fcen,df,nfreq,fr)

    return flux_region 

def get_phases(sim,flux_obj):
    result = sim.get_eigenmode_coefficients(flux_obj,[1],eig_parity=mp.NO_PARITY)
    coeffs = result.alpha
          
    phase = np.angle(coeffs[0,0,0]) 
    if phase > 0:
      phase -= 2*np.pi

    return phase

def get_fluxes(input_flux,sim,flux_object):
    res = sim.get_eigenmode_coefficients(flux_object,[1],eig_parity=mp.NO_PARITY)
    coeffs = res.alpha
    mode_tran = abs(coeffs[0,0,0])**2

    return mode_tran

def interpret_data(fluxes,phases):
    delta_phases = []
    percent_transmitted = []
    initial_flux = fluxes[0]
    #initial_phase = phases[0]

    for flux in fluxes[1:]:
        percent_transmitted.append(flux/initial_flux)
     
    for phase in phases[1:]:
        delta_phases.append(phase)

    return percent_transmitted,delta_phases

def build_plots(radii,fluxes,phases):
    plt.style.use('science')
    fig,ax = plt.subplots(1,2,figsize=(5,5))
    tickfontsize=14
    labelfontsize=16
    titlefontsize=22
    ax[0].set_title("Transmission Magnitude",fontsize=titlefontsize)
    ax[0].plot(radii,fluxes,'b',label='Transmission')
    ax[0].set_xlabel("Radius (nm)",fontsize=labelfontsize)
    ax[0].set_xticks([0.075,0.100,0.125,0.150,0.175,0.200,0.225,0.250])
    ax[0].set_xticklabels([75,100,125,150,175,200,225,250],fontsize=tickfontsize)
    ax[0].set_ylim([0,1])
    ax[0].set_ylabel(r'Percent Transmitted ($\%$)',fontsize=16)
    ax[0].set_yticks([0,0.2,0.4,0.6,0.8,1])
    ax[0].set_yticklabels([0,20,40,60,80,100],fontsize=tickfontsize)
    ax[0].legend(loc="upper right")
    

    ax[1].set_title("Transmission Phase",fontsize=titlefontsize)
    ax[1].plot(radii,phases,'r',label='Phase')
    ax[1].set_xlabel("Radius (nm)",fontsize=labelfontsize)
    ax[1].set_xticks([0.075,0.100,0.125,0.150,0.175,0.200,0.225,0.250])
    ax[1].set_xticklabels([75,100,125,150,175,200,225,250],fontsize=tickfontsize)
    ax[1].set_ylim([-2*np.pi,0])
    ax[1].set_ylabel(r'Phase Delay (rad)',fontsize=labelfontsize)
    ax[1].set_yticks([-2*np.pi, -(3*np.pi) / 2, -np.pi, -np.pi / 2, 0])
    ax[1].set_yticklabels([r'-2$\pi$',r'-$\frac{3\pi}{2}$',r'-$\pi$',r'-$\frac{\pi}{2}$',r'0'],fontsize=tickfontsize+4)
    ax[1].legend(loc="upper right")

    plt.tight_layout()
    plt.show()

def plot_geometry(sim):
    plt.figure(dpi=100)
    sim.plot2D() 
    plt.show()

def run(params):
    num=5
    fluxes = []
    phases = []
    radii = []
    mp.verbosity(0)

    #wave propagating without meta-atom
    sources = build_source(params)
    sim = build_sim(params,sources)
    flux_obj = build_flux_region(params,sim)
    sim.run(until=39)

    fluxes.append(mp.get_fluxes(flux_obj)[0])
    phases.append(get_phases(sim,flux_obj))
    
    # sim.plot2D(fields=mp.Ez)
    # plt.show()
    sim.reset_meep()
    # input()

    #build geometry, sources, and simulation
    for radius in tqdm(np.linspace(0.075,0.25,num=num),leave=False):
        geometry = build_geometry(params,radius)
        sources = build_source(params)
        sim = build_sim(params,sources,geometry)
    
        #get flux from simulation
        flux_obj = build_flux_region(params,sim)

        sim.run(until=39)

        fluxes.append(get_fluxes(fluxes[0],sim,flux_obj))
        phases.append(get_phases(sim,flux_obj))
        radii.append(radius)

        sim.reset_meep()

    percent_transmitted,delta_phases = interpret_data(fluxes,phases)
    build_plots(radii,percent_transmitted,delta_phases)

    import csv

    fields = ['radii','percent transmitted','delta phases']
    rows=[radii,percent_transmitted,delta_phases]

    with open('data.csv','w') as csvfile:
        writer = csv.writer(csvfile,delimiter=',')
        writer.writerow(fields)
        #np.transpose(radii)
        writer.writerows(rows)
        # writer.writerows(percent_transmitted)
        # writer.writerows(delta_phases)