import numpy
import pareto
import pandas
import matplotlib.figure
import matplotlib.backends.backend_agg as agg

data = pandas.read_table("data.txt", header=None, sep=" ")
sets = {}

fig = matplotlib.figure.Figure(figsize=(15,15))
agg.FigureCanvasAgg(fig)

counter = 0

resolutions = [1e-9, 0.03, 0.06, 0.1, 0.15, 0.2, 0.25, 0.4, 1.0]
for resolution in resolutions:
    sets[resolution] = pandas.DataFrame(data=pareto.eps_sort([data.itertuples(False)], [0,1], [resolution]*2).archive)
    
for resolution in resolutions:
    counter += 1
    ax = fig.add_subplot(3,3,counter)
    ax.scatter(data[0], data[1], lw=0, facecolor=(0.7, 0.7, 0.7), zorder=-1)
    ax.scatter(sets[resolution][0], sets[resolution][1], facecolor=(1.0, 1.0, 0.4), zorder=1, s=50)
    if resolution < 0.1:
        ax.set_xticks(numpy.arange(0, 1.5, 0.2))
        ax.set_yticks(numpy.arange(0, 1.5, 0.2))
    else:
        ax.set_xticks(numpy.arange(0, 1.2, resolution))
        ax.set_yticks(numpy.arange(0, 1.2, resolution))
    
    if resolution > 0.001:
        ax.hlines(numpy.arange(0, 1.4, resolution), 0, 1.4, colors=(0.9, 0.9, 0.9), zorder=-2)
        ax.vlines(numpy.arange(0, 1.4, resolution), 0, 1.4, colors=(0.9, 0.9, 0.9), zorder=-2)
    ax.set_xlim(0, 1.2)
    ax.set_ylim(0, 1.2)
    ax.set_title("Epsilon resolution: {0:.2g}".format(resolution))
                      
fig.savefig("picture")
"""
needs shading for dominated boxes
"""
