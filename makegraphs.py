from cProfile import label
from lib2to3.pgen2.literals import simple_escapes
from pickle import NONE
from re import L
import meep as mp
import matplotlib.pyplot as plt 
import numpy as np


def plot_location_data(data):
    plt.style.use('science')
    tickfontsize=14
    labelfontsize=16
    titlefontsize=22
    fig,ax = plt.subplots(1, 2, figsize=(11,8)) # may want to change. Also see DPI keyword

    fig.suptitle("Meta-atom Transmission",fontsize=titlefontsize)

    ax[0].set_xlabel("radii",fontsize=labelfontsize)
    ax[0].set_ylabel("phase",fontsize=labelfontsize)

    ax[1].set_xlabel("radii",fontsize=labelfontsize)
    ax[1].set_ylabel("trasmission",fontsize=labelfontsize)

    temp = np.linspace(1,2.5,7)
    for i, location in enumerate(data):
        radii,phases,transmitted = location
        ax[0].plot(radii,phases,label=f"location {temp[i]:.2f}")
        ax[1].plot(radii,transmitted,label=f"location {round(temp[i],1)}")
    ax[0].legend(loc="upper right")
    ax[1].legend(loc="lower right")
    plt.show()

if __name__=="__main__":
    data = np.load('data.npy',allow_pickle=True)

    radii = []
    transmitted = []
    phases = []

    for n in range(7):
        flux_location = data[n]
        r=np.asarray(data[n][0])
        radii.append(r)
        t=np.asarray(data[n][1])
        transmitted.append(t)
        p=np.asarray(data[n][2])
        phases.append(p)

    print(r)
    print(t)
    print
    # radii=np.delete(radii[0],0)
    # phases[0]=np.delete(phases[0],0)
    # transmitted=transmitted

    # temp = np.zeros((7,14))
    # temp[:] = radii
    # radii = temp
    # radii=np.asarray(radii)
    # phases=np.asarray(phases)
    # transmitted=np.asarray(transmitted)
    # print(radii.shape)
    # print(phases.shape)
    # print(transmitted.shape)
    # data = np.zeros((3,7,14)) 
    # data[0] = np.asarray(radii)
    # data[1] = np.asarray(phases)
    # data[2] = np.asarray(transmitted)

    # data = data.swapaxes(0,1)
    # from IPython import embed
    # embed()
    # exit()

    #plot_location_data(data)

    #plot_data(radii,phases,transmitted)
    
