
import csv
import matplotlib.pyplot as plt
import numpy as np
plt.style.use(['science'])


def read_csv(path):
    data = []
    with open(path) as csvFile:
        reader = csv.reader(csvFile, delimiter=',')
        for row in reader:
            data.append(row)
    return data



def plot_data(radii, fluxes, phases):
    plt.style.use('science')
    tickfontsize=14
    labelfontsize=16
    titlefontsize=22
    fig,ax = plt.subplots(1, figsize=(5,5)) # you may want to change. Also see DPI keyword

    plt.title("Meta-atom Transmission",fontsize=titlefontsize)
    ax.plot(radii, transmitted, c='b', label='Transmission')  #legend thing
    ax.set_xlabel("Radius (nm)",fontsize=labelfontsize)
    ax.set_xticks([0.075,0.100,0.125,0.150,0.175,0.200,0.225,0.250])
    ax.set_xticklabels([75,100,125,150,175,200,225,250],fontsize=tickfontsize)   
    ax.set_ylabel(r'Percent Transmitted ($\%$)',fontsize=16)
    ax.set_ylim([0,1])
    ax.set_yticks([0,0.2,0.4,0.6,0.8,1])
    ax.set_yticklabels([0,20,40,60,80,100],fontsize=tickfontsize)

    ax2 = ax.twinx()
    ax2.plot(radii, phases, c='r', label='Phase Delay')
    ax2.set_ylabel("Radius (nm)",fontsize=labelfontsize)
    ax2.set_xticks([0.075,0.100,0.125,0.150,0.175,0.200,0.225,0.250])
    ax2.set_xticklabels([75,100,125,150,175,200,225,250],fontsize=tickfontsize)
    ax2.set_ylim([-2*np.pi,0])
    ax2.set_ylabel(r'Phase Delay (rad)',fontsize=labelfontsize)
    ax2.set_yticks([-2*np.pi, -(3*np.pi) / 2, -np.pi, -np.pi / 2, 0])
    ax2.set_yticklabels([r'-2$\pi$',r'-$\frac{3\pi}{2}$',r'-$\pi$',r'-$\frac{\pi}{2}$',r'0'],fontsize=tickfontsize+4)
    
    plt.show()

def build_separate_plots(radii,fluxes,phases):
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

if __name__=="__main__":
    
    path = "datafor200to250radwithYat1.csv"

    #Read the CSV
    data = read_csv(path)

    #Split the data into lists
    radii = np.asarray(data[1])
    transmitted = np.asarray(data[2])
    phases = np.asarray(data[3])

    #Need to convert the strings from CSV to floats
    phases = [float(phase) for phase in phases]
    transmitted = [float(tran) for tran in transmitted]
    radii = [float(rad) for rad in radii]

    #print(radii)
    print(transmitted)
    #print(phases)
    #Plot it
    plot_data(radii, transmitted, phases)
    

