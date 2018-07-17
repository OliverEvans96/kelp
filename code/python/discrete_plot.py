import ipyvolume as ipv
from scipy.ndimage import zoom
import numpy as np
import matplotlib.pyplot as plt

def imshow_with_contours(x_list, y_list, f):
    # imshow
    extent = (min(x_list), max(x_list), min(y_list), max(y_list))
    plt.imshow(f.T, origin='lower')
    cbar = plt.colorbar()
    
    # contours
    x_grid, y_grid = np.meshgrid(x_list, y_list, indexing='ij')
    cont = plt.contour(f.T, cmap='magma_r')
    cbar.add_lines(cont)
    
    # Ticks
    xpos = np.arange(len(x_list))
    xlabels = map(str, x_list)
    ypos = np.arange(len(y_list))
    ylabels = map(str, y_list)
    plt.xticks(xpos, xlabels)
    plt.yticks(ypos, ylabels)

def double_n(arr, n):
    for i in range(n):
        arr = zoom(arr, 2, order=0)
    return arr

def imshow_with_contours_and_zoom(x_list, y_list, f, zoom_factor):
    # zoom_factor is the number of times the grid is doubled
    # (actually a zoom by 2 ** zoom_factor)
    
    # Assume constant spacing in each dimension
    dx = np.mean(np.diff(x_list))
    dy = np.mean(np.diff(y_list))
    
    # zoom
    x_grid, y_grid = np.meshgrid(x_list, y_list, indexing='ij')
    new_f = double_n(f, zoom_factor)
    
    # extent
    extent = (
        min(x_list)-dx/2, 
        max(x_list)+dx/2, 
        min(y_list)-dy/2, 
        max(y_list)+dy/2
    )
    
    # imshow
    plt.imshow(new_f.T, origin='lower', extent=extent)
    cbar = plt.colorbar()
    
    # contours
    cont = plt.contour(new_f.T, cmap='magma_r', extent=extent)
    cbar.add_lines(cont)
    
    # Ticks
    plt.xticks(x_list)
    plt.yticks(y_list)

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

def volshow_zoom_correct_scale(x_list, y_list, z_list, f, zoom_factor):
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