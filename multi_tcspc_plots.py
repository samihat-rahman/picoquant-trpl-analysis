import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from phu_to_df import read_phu

# Time axis
def make_time_axis(bin_number, resolution, units, half_index):
    bin_list = np.arange(0, bin_number, 1) # list from 0 of the total number of bins
    time_data = np.multiply(bin_list, resolution) # multiply the list by the resolution to get the time
    # change the units
    if units == 'us':
        time_axis = np.multiply(time_data, 1E-6)
    if units == 'ns':
        time_axis = np.multiply(time_data, 1E-3)
    if units == 'ps':
        time_axis = time_data
    time_corrected = time_axis - time_axis[max_index] # adjust the time scale so time 0 is the max count
    return time_corrected

# Background subtraction
def background_subtraction(data):
    uncorrected_counts = data
    max_counts = np.max(data)
    background_counts = []
    # calculates the background as 0.5% of the peak (can change if needed)
    for count in data[:np.argmax(data)]:
        if count <= 0.005 * max_counts:
            background_counts.append(count)
    background = np.mean(background_counts) # the avg background
    bgcorrected_counts = uncorrected_counts - background # subtraction
    return bgcorrected_counts

# Normalize to peak
def normalization(bg_corrected_data):
    max_count = np.max(bg_corrected_data)
    corrected_data = np.divide(bg_corrected_data, max_count)
    return corrected_data

# Find average lifetime
def riemann_sum(time, counts):
    max_index = np.argmax(counts)
    lifetime = np.sum(np.multiply(time[max_index:], counts[max_index:]))/np.sum(counts[max_index:])
    return lifetime
  
# Finds max index of half peak (for a group)
def find_half_peak(data_list):
    half_max_index = 0
    for data in data_list:
        peak_height = np.max(data) - np.mean(data[:50])
        half_peak_index = (np.abs(data[:np.argmax(data)] - peak_height/2)).argmin()
        if half_max_index < half_peak_index:
            half_max_index = half_peak_index
    return half_max_index

# At half max
def roll_counts_half(data, half_index):
    peak_height = np.max(data) - np.mean(data[:50])
    half_peak_index = (np.abs(data[:np.argmax(data)] - peak_height/2)).argmin()
    rolled = np.roll(data, (half_index - half_peak_index))
    return rolled

# Plots a data set on the same canvas: this is rolled at half max
def trpl_multi_plot(file_list, labels, resolution, units, xlims, ylims, color_palette='default', lifetimes=False, output_file='No Output'):
    # I like using this font
    #plt.style.use('dark_background')
    plt.rcParams['font.family'] = 'Open Sans'
    # Auto change the labels of units
    unit_label = 'µs'
    if units == 'ns':
        unit_label = 'ns'
    if units == 'ps':
        unit_label = 'ps'
    # Make the canvas
    plt.figure(figsize=(8, 8))
    # bins
    bin_number = 65536
    # Store the counts
    count_list = []
    # read data
    for index, file in enumerate(file_list):
        # counts
        uncorrected_counts = read_phu(file)[0]
        bg_subtracted_counts = background_subtraction(uncorrected_counts)
        norm_counts = normalization(bg_subtracted_counts)
        count_list.append(norm_counts)
    # from the counts, find the latest time we get the half peak
    half_max_index = find_half_peak(count_list)
    
    # spacing between text
    spacing = 0
    # colors for different counts
    colors = list(mcolors.TABLEAU_COLORS.keys())
    if color_palette != 'default':
      colors = sns.color_palette(color_palette, n_colors=len(file_list))

    text_positions = []
    for index in range(len(file_list)):
        # text positions
        text_x = 0.55 - 0.1
        text_y = 0.97
        text_p = (text_x, text_y - spacing)
        text_positions.append(text_p)
        spacing += 0.05
      
    # time
    time_axis = make_time_axis(bin_number, resolution, units, max_index)
    for index, counts in enumerate(count_list):
        # plot
        counts = roll_counts_half(counts, half_ max_index)
        if lifetimes:
            lifetime = riemann_sum(time_axis, counts)
            plt.text(text_positions[index][0], text_positions[index][1], f'$\\langle$\N{GREEK SMALL LETTER TAU}$\\rangle$ = {lifetime:0.3g} {unit_label}', transform=plt.gca().transAxes, fontsize=22, fontweight='regular', ha='left', va='top', color=colors[index])
        plt.semilogy(time_axis, counts, 'o', label = labels[index], mfc='none', mec=colors[index])
        
    # input x and y limits as a list, feel free to comment this out
    plt.ylim(ylims[0], ylims[1])
    plt.xlim(xlims[0], xlims[1])
    plt.xlabel(f'Time ({unit_label})', size=26)
    plt.ylabel('PL Intensity (Norm.)',  size=26)
    
    # Only change this to change the fontsize of everything
    plt.tick_params(axis='both', labelsize=22)

    # Legend
    plt.legend(markerscale=2, fontsize=16)
    plt.tight_layout()

    # If you want to save
    if output_file != 'No Output':
        plt.savefig(f'{output_file}.png', dpi=300)
      
    # open a preview
    plt.show()
    plt.close()
    return time_axis, count_list
