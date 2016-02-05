import matplotlib.pyplot as plt
import matplotlib.lines as mlines

def plot_principal_components(eigenvectors, output_filename):
    f, ax = plt.subplots(1, 1, figsize=(8, 8))
    f.suptitle('Population Structure: PCA')
    ax.set_xlabel('First Eigenvector')
    ax.set_ylabel('Second Eigenvector')
    #generational_eigenvectors = np.split(eigenvectors, number_of_generations, axis=1)
    #for i in range(number_of_generations):
    ax.plot(eigenvectors[:100, 0], eigenvectors[:100, 1], 'r*', linewidth=0.0, markersize=5, alpha=0.5)
    ax.plot(eigenvectors[100:200, 0], eigenvectors[100:200, 1], 'c^', linewidth=0.0, markersize=5, alpha=0.5)
    ax.plot(eigenvectors[200:300, 0], eigenvectors[200:300, 1], 'go', linewidth=0.0, markersize=5, alpha=0.5)
    ax.plot(eigenvectors[300:400, 0], eigenvectors[300:400, 1], 'y8', linewidth=0.0, markersize=5, alpha=0.5)
    ax.plot(eigenvectors[400:, 0], eigenvectors[400:, 1], 'mv', linewidth=0.0, markersize=5, alpha=0.5)
    plt.savefig(output_filename, dpi=300)

def plot_multiple_replicate_population_structure(eigenvectors, output_filename):
    """
    Generates a plot of population stucture of multiple replicates of the meta-population. Parameter
    eigenvectors is a list of (pop_size)*(2) numpy.arrays because only the first two eigenvectors
    are considered. Generations are given different markers and colors for purposes of comparison.
    """
    import matplotlib.pyplot as plt
    import matplotlib.lines as mlines
    f, axes = plt.subplots(2, 2, figsize=(16, 16))
    f.suptitle("Population Structure Under Five Different Scenarios")

    itervecs = iter(eigenvectors)
    for i in range(2):
        for j in range(2):
            e = next(itervecs)
            axes[i, j].plot(e[:100, 0], e[:100, 1], 'r*',
                            label='Generation One', linewidth=0.0, markersize=5, alpha=0.5)
            axes[i, j].plot(e[100:200, 0], e[100:200, 1], 'c^',
                            label='Generation Two', linewidth=0.0, markersize=5, alpha=0.5)
            axes[i, j].plot(e[200:300, 0], e[200:300, 1], 'go',
                            label='Generation Three', linewidth=0.0, markersize=5, alpha=0.5)
            axes[i, j].plot(e[300:400, 0], e[300:400, 1], 'y8',
                            label='Generation Four', linewidth=0.0, markersize=5, alpha=0.5)
            axes[i, j].plot(e[400:, 0], e[400:, 1], 'mv',
                            label='Generation Five', linewidth=0.0, markersize=5, alpha=0.5)
    axes[0, 0].set_xlabel('First Eigenvector')
    axes[0, 0].set_ylabel('Second Eigenvector')

    red_stars = mlines.Line2D([], [], color='red', marker='*', markersize=8, label='Generation One')
    cyan_triangles = mlines.Line2D([], [], color='cyan', marker='^', markersize=8, label='Generation Two')
    green_circles = mlines.Line2D([], [], color='green', marker='o', markersize=8, label='Generation Three')
    yellow_octagons = mlines.Line2D([], [], color='yellow', marker='8', markersize=8, label='Generation Four')
    inverted_triangles = mlines.Line2D([], [], color='m', marker='v', markersize=8, label='Generation Five')

    labels = ['Generation One', 'Generation Two', 'Generation Three', 'Generation Four', 'Generation Five']
    plt.figlegend(handles=[red_stars, cyan_triangles,
                           green_circles, yellow_octagons,
                           inverted_triangles], labels=labels, loc='upper right')

    axes[0, 0].set_title('Selection with Random Mating')
    axes[0, 1].set_title('Drift with Random Mating, Largest')
    axes[1, 0].set_title('Selection without Random Mating')
    axes[1, 1].set_title('Drift without Random Mating')
    plt.savefig(output_filename, dpi=300)
