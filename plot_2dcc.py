#!/usr/local/bin/python3
# coding: utf-8

# In[27]:

import nmrglue as ng
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm
import matplotlib


def plot_2d_cc(ddir, sigma=4):
    """
    Ploter for typical Protein CC data.
    """
    dic, data = ng.sparky.read_2D(ddir)

    uc_13c_a = ng.sparky.make_uc(dic, data, dim=1)
    ppm_13c_a = uc_13c_a.ppm_scale()
    uc_13c_b = ng.sparky.make_uc(dic, data, dim=0)
    ppm_13c_b = uc_13c_b.ppm_scale()

    noise = np.std(data[uc_13c_a(200, 'ppm'):uc_13c_a(190, 'ppm'),
                        uc_13c_b(210, 'ppm'):uc_13c_b(200, 'ppm')]) * sigma

    # plot parameters
    cmap = matplotlib.cm.Blues_r    # contour map (colors to use for contours)
    contour_start = noise           # contour level start value
    contour_num = 16                # number of contour levels
    contour_factor = 1.40          # scaling factor between contour levels

    # calculate contour levels
    cl = contour_start * contour_factor ** np.arange(contour_num)

    fig=plt.figure()
    #plt.subplot(1, 2, 1)
    # plt.contour(ppm_13c_a, ppm_13c_b, data, cl, cmap=cmap)
    # plt.xlim([165, 185])
    # plt.ylim([10, 80])
    # plt.gca().invert_xaxis()
    # plt.gca().invert_yaxis()

    # plt.subplot(1, 2, 2)
    # plt.contour(ppm_13c_a, ppm_13c_b, data, cl, cmap=cmap)
    # plt.xlim([10, 80])
    # plt.ylim([10, 80])
    # plt.gca().invert_xaxis()
    # plt.gca().invert_yaxis()

    ax1=plt.subplot2grid((14,19), (0,0), colspan=4,rowspan=14)

    plt.contour(ppm_13c_a, ppm_13c_b, data, cl, cmap=cmap)
    plt.xlim([165, 185])
    plt.ylim([10, 80])
    plt.gca().set_xticks([180, 170])
    plt.gca().invert_xaxis()
    plt.gca().invert_yaxis()
    plt.ylabel("13C/ppm")

    ax2=plt.subplot2grid((14,19), (0,5), colspan=14,rowspan=14)
    plt.contour(ppm_13c_a, ppm_13c_b, data, cl, cmap=cmap,linewidth=1.0)
    plt.xlim([10, 80])
    plt.ylim([10, 80])
    plt.gca().invert_xaxis()
    plt.gca().invert_yaxis()
    plt.xlabel("13C/ppm")
    plt.title(os.path.basename(ddir))

    fig.savefig(os.path.basename(ddir) + '.pdf')

# def plot_2d_cc_region(ddir, x_range=[60,70], y_range=[60,70],sigma=5):
#     """
#     Ploter for  Protein CC data for certain ranges
#     """
#     dic, data = ng.sparky.read_2D(ddir)
#
#     uc_13c_a = ng.sparky.make_uc(dic, data, dim=1)
#     ppm_13c_a = uc_13c_a.ppm_scale()
#     uc_13c_b = ng.sparky.make_uc(dic, data, dim=0)
#     ppm_13c_b = uc_13c_b.ppm_scale()
#
#     noise = np.std(data[uc_13c_a(200, 'ppm'):uc_13c_a(190, 'ppm'),
#                         uc_13c_b(210, 'ppm'):uc_13c_b(200, 'ppm')]) * sigma
#
#     # plot parameters
#     cmap = matplotlib.cm.Blues_r    # contour map (colors to use for contours)
#     contour_start = noise           # contour level start value
#     contour_num = 16                # number of contour levels
#     contour_factor = 1.40          # scaling factor between contour levels
#
#     # calculate contour levels
#     cl = contour_start * contour_factor ** np.arange(contour_num)
#
#     fig=plt.figure()
#     plt.subplot(1, 1, 1)
#     plt.contour(ppm_13c_a, ppm_13c_b, data, cl, cmap=cmap)
#     plt.xlim(x_range)
#     plt.ylim(y_range)
#     plt.gca().invert_xaxis()
#     plt.gca().invert_yaxis()
#
#     # plt.subplot(1, 2, 2)
#     # plt.contour(ppm_13c_a, ppm_13c_b, data, cl, cmap=cmap)
#     # plt.xlim([10, 80])
#     # plt.ylim([10, 80])
#     # plt.gca().invert_xaxis()
#     # plt.gca().invert_yaxis()
#
#     plt.ylabel("13C/ppm")
#     plt.xlabel("13C/ppm")
#     plt.title(os.path.basename(ddir))
#
#     fig.savefig(os.path.basename(ddir) + '.pdf')


if __name__ == '__main__':
    import os
    import os.path
    import subprocess
    import argparse

    #matplotlib.rcParams['pdf.fonttype']='none'

    # def spec_range(s):
    #     """
    #     Type conversion function for argument parser
    #     """
    #     x_min, x_max, y_min, y_max = map(float, s.split(', '))
    #     return x_min, x_max, y_min, y_max
    #     # try:
    #     #     x_min, x_max, y_min, y_max = map(float, s.split(', '))
    #     #     return x_min, x_max, y_min, y_max
    #     # except:
    #     #     msg = 'Range must be x_min, x_max, y_min, y_max'
    #     #     raise argparse.ArgumentTypeError(msg)

    parser = argparse.ArgumentParser(description='Plot 2D 13C-13C spectrum.')
    parser.add_argument('file_name',type=str, help='file path to data')
    parser.add_argument('-w', '--whole_region',
                        action='store',
                        default= True,
                        type=bool, help='True to plot the whole region')
    # parser.add_argument('-r', '--range',
    #                     type=spec_range,
    #                     help='x_min, x_max, y_min, y_max')

    parser_dict = vars(parser.parse_args())
    ddir = parser_dict['file_name']

    # if a bruker file transfer into a sparky format
    for dirpath, dirnames, filenames in os.walk(ddir):
        # if 'Sep0922' in dirpath:
            for filename in [f for f in filenames if f.endswith("2rr")]:
                #for i in dirpath.split('/'):
                    # if 'UKcsA'in i:
                        name = filename + '.ucsf'
                        dirna = os.path.join(dirpath, filename)
                        dirna2 = os.path.join(dirpath, name)
                        # print(dirna)
                        subprocess.Popen("/Applications/Sparky.app/Contents/Resources/bin/bruk2ucsf %s %s" % (dirna, dirna2),
                                         shell=True, executable='/bin/sh')
                        # then read into nmrglue

                        #print(parser_dict['range'])
                        if parser_dict['whole_region']:
                             plot_2d_cc(dirna2)
                             print(dirna2)
                        # else:
                        #     plot_2d_cc_region(dirna2,
                        #     x_range = parser_dict['range'][0],
                        #     y_range= x_range=parser_dict['range'][1])
                        plt.show()