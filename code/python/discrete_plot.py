import ipyvolume as ipv
from scipy.ndimage import zoom
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

def imshow_with_contours(x_list, y_list, f, imshow_kwargs={}, contour_kwargs={}, cbar_kwargs={}):
    # imshow
    extent = (min(x_list), max(x_list), min(y_list), max(y_list))
    plt.imshow(f.T, origin='lower', **imshow_kwargs)
    cbar = plt.colorbar()

    # contours
    x_grid, y_grid = np.meshgrid(x_list, y_list, indexing='ij')
    contour_kwargs = {'cmap':'magma_r', **contour_kwargs}
    cont = plt.contour(f.T, **contour_kwargs)
    cbar.add_lines(cont, **cbar_kwargs)

    # Ticks
    xpos = np.arange(len(x_list))
    convert_to_str = lambda x: '{:.2f}'.format(x) if isinstance(x, float) else str(x)
    xlabels = map(convert_to_str, x_list)
    ypos = np.arange(len(y_list))
    ylabels = map(convert_to_str, y_list)
    plt.xticks(xpos, xlabels)
    plt.yticks(ypos, ylabels)

def double_n(arr, n):
    for i in range(n):
        arr = zoom(arr, 2, order=0)
    return arr

def imshow_with_contours_and_zoom(x_list, y_list, f, zoom_factor, imshow_kwargs={}, contour_kwargs={}, cbar_kwargs={}, log_data=False):
    # zoom_factor is the number of times the grid is doubled
    # (actually a zoom by 2 ** zoom_factor)

    # Assume constant spacing in each dimension
    dx = np.mean(np.diff(x_list))
    dy = np.mean(np.diff(y_list))

    # zoom
    x_grid, y_grid = np.meshgrid(x_list, y_list, indexing='ij')
    new_f = double_n(f, zoom_factor)

    # extent
    # extent = (
    #     min(x_list)-dx/2,
    #     max(x_list)+dx/2,
    #     min(y_list)-dy/2,
    #     max(y_list)+dy/2
    # )

    extent = (
        0, len(x_list),
        0, len(y_list)
    )

    # imshow
    ax = plt.gca()
    im = plt.imshow(new_f.T, origin='lower', extent=extent, **imshow_kwargs)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)

    cbar = plt.colorbar(cax=cax, **cbar_kwargs)

    # contours
    contour_kwargs = {'cmap':'magma_r', **contour_kwargs}
    cont = ax.contour(new_f.T, extent=extent, **contour_kwargs)
    cbar.add_lines(cont)

    # Ticks
    convert_to_str = lambda x: '{:.2f}'.format(x) if isinstance(x, float) else str(x)
    x_labels = map(convert_to_str, x_list)
    y_labels = map(convert_to_str, y_list)
    ax.set_xticks(np.arange(len(x_list))+0.5)
    ax.set_xticklabels(x_labels)
    ax.set_yticks(np.arange(len(y_list))+0.5)
    ax.set_yticklabels(y_labels)

    if log_data:
        plt.yscale('log')

# Due to a bug in IPyVolume (#117, fixed but not released),
# there are two versions: One with the correct scale
# and one with isosurfaces.

def volshow_with_isoplanes_and_zoom(x_list, y_list, z_list, f, zoom_factor):
    # zoom_factor is the number of times the grid is doubled
    # (actually a zoom by 2 ** zoom_factor)

    # Assume constant spacing in each dimension
    dx = np.mean(np.diff(x_list))
    dy = np.mean(np.diff(y_list))
    dz = np.mean(np.diff(z_list))

    # zoom
    x_grid, y_grid, z_grid = np.meshgrid(x_list, y_list, z_list, indexing='ij')
    new_f = double_n(f, zoom_factor)

    # data extent
    extent = (
        (min(x_list)-dx/2, max(x_list)+dx/2),
        (min(y_list)-dy/2, max(y_list)+dy/2),
        (min(z_list)-dz/2, max(z_list)+dz/2)
    )

    ipv.figure()

    # isosurfaces
    ipv.plot_isosurface(new_f)

    # volshow
    ipv.volshow(new_f.T)

    ipv.show()

def volshow_zoom_correct_scale(x_list, y_list, z_list, f, zoom_factor=0):
    # zoom_factor is the number of times the grid is doubled
    # (actually a zoom by 2 ** zoom_factor)

    # Assume constant spacing in each dimension
    dx = np.mean(np.diff(x_list))
    dy = np.mean(np.diff(y_list))
    dz = np.mean(np.diff(z_list))

    # zoom
    x_grid, y_grid, z_grid = np.meshgrid(x_list, y_list, z_list, indexing='ij')
    new_f = double_n(f, zoom_factor)

    # data extent
    extent = (
        (min(x_list)-dx/2, max(x_list)+dx/2),
        (min(y_list)-dy/2, max(y_list)+dy/2),
        (min(z_list)-dz/2, max(z_list)+dz/2)
    )

    ipv.figure()

    # volshow
    ipv.volshow(new_f.T)

    # figure extent
    ipv.xlim(min(x_list)-dx/2, max(x_list)+dx/2)
    ipv.ylim(min(y_list)-dy/2, max(y_list)+dy/2)
    ipv.zlim(min(z_list)-dz/2, max(z_list)+dz/2)

    ipv.show()
