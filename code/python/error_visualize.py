import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

import kelp_analyze
import discrete_plot

def plot_two_edges(rel_err_arr, abs_err_arr, ns_list, na_list, ylim=None):
    plt.figure(figsize=(6,6))
    
    ax1 = plt.subplot(1,1,1)
    ax1.semilogy(na_list, rel_err_arr[-1,:], 'C1o-', label='na')
    ax1.set_xlabel('na')
    ax1.legend(loc='right')

    ax2 = ax1.twiny()
    ax2.semilogy(ns_list, rel_err_arr[:,-1], 'C0o-', label='ns')
    ax2.set_xlabel('ns')
    ax2.legend(loc='upper right')
    
    ax1.set_ylabel('Average relative error (Perceived Irradiance)')
    if ylim:
        ax1.set_ylim(*ylim)
    
    # ax3 = plt.subplot(1,2,2)
    # ax3.semilogy(na_list, abs_err_arr[-1,:], 'C1o-', label='na')
    # ax3.set_xlabel('na')
    # ax3.legend(loc='lower right')

    # ax4 = ax3.twiny()
    # ax4.semilogy(ns_list, abs_err_arr[:,-1], 'C0o-', label='ns')
    # ax4.set_xlabel('ns')
    # ax4.legend(loc='upper right')
    # ax3.set_ylabel('Average absolute error (Perceived Irradiance)')

def plot_2d_resolution_grid(rel_err_arr, abs_err_arr, ns_list, na_list, vlim=None, xlabel='ns', ylabel='na'):
    # ns max
    plt.figure(figsize=(6,6))
    ax1 = plt.subplot(1,1,1)
    
    imshow_kwargs = {'norm': LogNorm()}
    
    if vlim:
        imshow_kwargs = {
            'vmin': vlim[0],
            'vmax': vlim[1],
            **imshow_kwargs
        }
    
    plt.title('log10 rel err (Perceived Irradiance)')
    discrete_plot.imshow_with_contours_and_zoom(
        ns_list, 
        na_list, 
        rel_err_arr, 
        zoom_factor=4,
        imshow_kwargs=imshow_kwargs,
        contour_kwargs={
            #'levels': np.logspace(-3,0,7),
            'norm': LogNorm()
        },
        cbar_kwargs={
            'label': 'average relative error',
        }

    )
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    
    # ax2 = plt.subplot(1,2,2)

    # plt.title('log10 abs err (Perceived Irradiance)')
    # discrete_plot.imshow_with_contours_and_zoom(
    #     ns_list, 
    #     na_list, 
    #     abs_err_arr, 
    #     zoom_factor=4,
    #     imshow_kwargs={
    #         'norm': LogNorm(),
    #         #'vmin': 1e-3,
    #         #'vmax': 1e0
    #     },
    #     contour_kwargs={
    #         'levels': np.logspace(-3,0,7),
    #         'norm': LogNorm()
    #     },
    #     cbar_kwargs={
    #         'label': 'average relative error',
    #     }
    # )
    # ax2.set_xlabel(xlabel)
    # ax2.set_ylabel(ylabel)
    
def plot_gs_edges_and_grid(rel_err_arr, abs_err_arr, ns_list, na_list):
    plot_two_edges(rel_err_arr, abs_err_arr, ns_list, na_list)
    plot_2d_resolution_grid(rel_err_arr, abs_err_arr, ns_list, na_list)