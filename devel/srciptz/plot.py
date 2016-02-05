import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def plot_single_generation(pop, output_filename):
    """
    Function is intended mainly as a memory aid in using mapltotlib/pylab/seaborn etc.
    I will probably end up having to modify the source code for every plot I want to make until I have
    a better grasp of matplotlib.

    Genotypic effect scores are stored in pop.dvars().parental_ge by generation. This code generates
    a pdf file of
    :param gens_to_plot:
    :return:
    """
    sns.set(context='paper', style='darkgrid', palette='deep')
    f, axes = plt.subplots(1, 1, figsize=(10, 10))
    f.suptitle("Single Generations of Selection")
    G_zero = np.array(pop.dvars().parental_ge[7])
    sns.distplot(G_zero, ax=axes)
    axes.set_title(r"$G_0$")
    axes.set_xlabel(r"$G$")
    axes.set_ylabel(r"$frequency$")
    zero_x_bar = G_zero.mean()
    zero_sigma = G_zero.var()
    text_str = '$\mu=%.2f$\n$\sigma=%.2f$' % (zero_x_bar, zero_sigma)
    axes.text(120, 0.07, text_str, fontsize=14)
    plt.savefig(output_filename, dpi=300)

plot_single_generation(pop, "single_gen_selection.pdf", dpi=300)
