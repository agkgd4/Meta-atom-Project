import matplotlib.pyplot as plt



def plot_data(eps, ez):

    fig, ax = plt.subplots(1,2)

    ax[0].imshow(eps.transpose(), interpolation='spline36', cmap='binary')
    ax[1].imshow(ez.transpose(), interpolation='spline36', cmap='RdBu', alpha=0.9)
    fig.show()
    