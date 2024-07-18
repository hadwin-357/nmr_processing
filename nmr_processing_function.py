# various functions for processing NMR spectra, examining NMR simulation results
# generate artificial NCOCX CANCO  NCACX spectra and compare them to experimental data
import nmrglue as ng
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import os
import seaborn as sns
import matplotlib.ticker as ticker
from matplotlib.ticker import AutoMinorLocator
from matplotlib import gridspec


# function to convert SIMPSON simulation output to spectrum-like
def convert_simp2spec(path, mag, ref_ppm):
    with open(path,'r') as file:
        firstNlines = file.readlines()[0:6]
        for line in firstNlines:
            # print(line)
            if 'NP' in line:
                # print(line.rstrip("\n"))
                NP = int(line.rstrip("\n").split('=')[1])
                # [print(s) for s in re.split('= ',line) if s.isdigit()]
            elif 'SW' in line:
                SW = int(line.rstrip("\n").split('=')[1])
    carrier_freq = mag
    span_ppm = SW / carrier_freq
    start_x = ref_ppm - span_ppm / 2
    end_x = ref_ppm + span_ppm / 2
    xxais = np.linspace(start_x, end_x, num=NP, endpoint=True)
    df = pd.read_csv(path, skiprows=5, skipfooter=1, index_col=False, sep=' ', header=None, names=['real', 'image'], engine='python')
    return (xxais, df['real'])

def read_1Ddata(path):
    new_path='/'+path.strip('pdata/1')
    print(new_path)
    dic, data=ng.bruker.read(new_path)
    # generate the ppm/ Hz scale for the data
    udic=ng.bruker.guess_udic(dic,data)
    uc = ng.fileiobase.uc_from_udic(udic)

#generate region for fitting
    ppm=uc.ppm_scale()
    return [udic, data, ppm, uc]

def read_1Dpdata(path):
    '''read 1D processed NMR spectrum
    input： path to bruker pdata folder
    '''
    dic,data = ng.bruker.read_pdata(path)
    # generate the ppm/ Hz scale for the data
    udic=ng.bruker.guess_udic(dic,data)
    uc = ng.fileiobase.uc_from_udic(udic)
    ppm=uc.ppm_scale()
    return [udic, data, ppm, uc]

def convert_spc2simp(Spec,file_name):
    '''
    :param Spec: Spec contains three inputs from read_pdata function and file_name
    :param file_name: output simposn file name
    :return: simpson file to be used for simpson simulation
    '''
    np=Spec[0][0]['size']
    sw=Spec[0][0]['sw']
    ref=Spec[0][0]['car']/Spec[0][0]['obs']
# write spe file for SIMPSON simulation
    with open(file_name,'w+') as file:
        file.write('SIMP\n')
        file.write('NP={}\n'.format(np))
        file.write('SW={}\n'.format(sw))
        file.write('#REF={}off_set_for_record\n'.format(ref))
        file.write('TYPE=SPE\n')
        file.write('DATA\n')
# write the data part
    df=pd.DataFrame(len(Spec[1])*[0],Spec[1][::-1])
    df.to_csv(file_name, sep=' ',mode='a', header=False)
    with open(file_name,'a') as file:
        file.write('END\n')

from scipy.stats import chisquare

def reduced_chi_square(observed, expected, degrees_of_freedom):
    """
    Calculate the reduced chi-square value.

    Parameters:
    observed (array-like): Observed data points.
    expected (array-like): Expected data points.
    degrees_of_freedom (int): Degrees of freedom (number of observations - number of fitted parameters - 1).

    Returns:
    float: Reduced chi-square value.
    """
    chi2_stat, p_val = chisquare(observed, expected)
    reduced_chi2 = chi2_stat / degrees_of_freedom
    return reduced_chi2



def plot_exp_simp(Simp_path, Spec):
    '''
    compare experimental and fitted spectra
    inputs: simp_path-path for simpson output
    Spec: experiment spectra output from read_1Dpdate function
    '''

    REF = Spec[0][0]['car'] / Spec[0][0]['obs']

    obs = Spec[0][0]['obs']
    span_ppm = SW / obs
    start_x = REF - span_ppm / 2
    # print(start_x)
    end_x = REF + span_ppm / 2
    # print(end_x)
    xxais = np.linspace(start_x, end_x, num=NP, endpoint=True)
    fig = plt.figure(figsize=(14, 6))
    ax = plt.subplot2grid((4, 4), (0, 0), colspan=4, rowspan=3)
    ppm_scale = Spec[2]
    # noise_std=np.std(Spec[1][Spec[3](-200,'ppm'):Spec[3](-300,'ppm')]/Spec[1].max())
    noise_std = np.std(Spec[1][Spec[3](-200, 'ppm'):Spec[3](-213, 'ppm')])

    df = pd.read_csv(Simp_path, skiprows=5, skipfooter=1, index_col=False, sep=' ', header=None,
                     names=['real', 'image'], engine='python')
    # y_data= Spec[1][1024:1024+2048][::-1]/Spec[1][1024:1024+2048][::-1].max()
    # model_data=df['real']/df['real'].max()
    # exp_data=Spec[1]/Spec[1].max()
    model_data,xxais=convert_simp2spec(Simp_path,obs,REF)
    exp_data = Spec[1]

    ax.plot(xxais, model_data, label='Simulated')
    ax.plot(ppm_scale, exp_data, label='Experimental')

    ax2 = plt.subplot2grid((4, 4), (3, 0), colspan=4, sharex=ax)
    # ax2.plot(xxais,df['real']-Spec[1][1024:1024+2048][::-1])
    ax2.plot(xxais, exp_data[::-1] - model_data, color='r')
    #ax2.set_ylim(-0.2 * 10 ** 10, 0.2 * 10 ** 10)

    #ax.set_xlim(300, -200)
    #ax.set_xlim(100, 40)
    #ax2.set_xlabel('31P chemical shift/ppm')
    #ax.set_title("31P spectra of Co6Se8(PEt3)5CO at 20KHz")
    ax.legend()
    ax.xaxis.set_major_locator(ticker.AutoLocator())
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    # fig.savefig("31P spectra of Co6Se8(PEt3)5CO at 20KHz_SimpvsExp.pdf",format='pdf', dpi=1000)

    reduced_chi = reduced_chi_square(exp_data[::-1], model_data, degrees_of_freedom=len(model_data))
    return reduced_chi

def plot_2dsparky_cc(ddir, sigma=4):
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


def shiftx_to_df(file_path):
    from io import StringIO
    """
    Read heavy atom SHIFTX2 tabular data
    """
    # The data is in several differnt tables
    frames = []
    with open(file_path, 'r') as fid:
        frame=''
        for line in fid:
            frame+=line
            if not line.strip():
                if frame.strip():
                    frames.append(frame)
                frame=''
        if frame.strip():
            frames.append(frame)
    
    # Parse the deperate tables and merge the tables
    df = pd.read_csv(StringIO(frames[0]), header=1, delim_whitespace=True,
            na_values='****')
    for frame in frames[1:]:
        dfi = pd.read_csv(StringIO(frame), header=1, delim_whitespace=True,
            na_values='****')
        df = pd.merge(df, dfi, on=['Num', 'RES'])
    
    # Change the  Dataframe header
    header = list(df.columns)
    header[4] = 'C'
    df.columns = header
    return df

#function for 3D plot
def ncocx_plot(df, cs_nitrogen, off_set=0.4, c_shift=0, n_shift=0, **kwargs):
    '''
    generate computed 3D plane using shfitx value for NCOCX
    paras: df- df output from shiftx_to_df
    cs_nitrogen- nitrogen plane, off_set: nitrogen chemical shift center at cs_nitrogen within offset wil be inluded 
    '''
    test={}
    res_limx=df['Num'][0]
    res_limy=df['Num'].iloc[-1]
    #read NCOCX into a dictionary
    for i in range(res_limx,res_limy):
        str='{}{}{}'.format(df[df.Num==i]['RES'].values[0],df[df.Num==i+1]['Num'].values[0],df[df.Num==i+1]['RES'].values[0])
        test[str]=[df[df.Num==i]['N'].values[0],df[df.Num==i+1]['C'].values[0],df[df.Num==i+1]['CA'].values[0]]
    # convert the dictionary to a dataframe
    new_pd=pd.DataFrame(test,index=['N','CO','CA'])
    #transpose the dataframe
    new_pd_Tran=new_pd.T
    #select slice based on nitrogen chemical shift
    sele=new_pd_Tran[abs(new_pd_Tran['N']-cs_nitrogen)<=off_set]
    #plot the data using seaborn
    sns.scatterplot(data=sele,x='CA', y = 'CO')
    for i in range(sele.shape[0]):
        plt.text(x=sele.CA[i]+0.1,y=sele.CO[i]+0.1,s=sele.index[i], 
          fontdict=dict(color= 'red',size=10))


def canco_plot(df, cs_nitrogen, off_set=0.4, c_shift=0, n_shift=0, **kwargs):
    test1={}
    res_limx=df['Num'][0]
    res_limy=df['Num'].iloc[-1]
    #read CANCO(i-1) into a dictionary
    for i in range(res_limx,res_limy):
            str1='{}{}{}'.format(df[df.Num==i-1]['RES'].values[0],df[df.Num==i-1]['Num'].values[0],df[df.Num==i]['RES'].values[0])
            test1[str1]=[df[df.Num==i-1]['N'].values[0],df[df.Num==i]['C'].values[0],df[df.Num==i-1]['CA'].values[0]]
    # convert the dictionary to a dataframe
    new_pd=pd.DataFrame(test1,index=['N','CO','CA'])
    #transpose the dataframe
    new_pd_Tran=new_pd.T
    #select slice based on nitrogen chemical shift
    sele=new_pd_Tran[abs(new_pd_Tran['N']-(cs_nitrogen+nshift))<=off_set]
    #plot the data using seaborn
    sns.scatterplot(data=sele,x='CA', y = 'CO', color='cyan')
    for i in range(sele.shape[0]):
        plt.text(x=sele.CA[i]+0.1,y=sele.CO[i]+0.1,s=sele.index[i], 
          fontdict=dict(color= 'blue',size=7))
#correct n(i)cocx(i-1)

def ncocx_plot(df, cs_nitrogen, off_set=0.4, c_shift=0, n_shift=0, **kwargs):
    test={}
    res_limx=df['Num'][0]+1
    res_limy=df['Num'].iloc[-1]
    #read NCOCX into a dictionary
    for i in range(res_limx,res_limy):
        str='{}{}{}'.format(df[df.Num==i]['RES'].values[0],df[df.Num==i-1]['Num'].values[0],df[df.Num==i-1]['RES'].values[0])
        test[str]=[df[df.Num==i]['N'].values[0],df[df.Num==i-1]['C'].values[0],df[df.Num==i-1]['CA'].values[0]]
    # convert the dictionary to a dataframe
    new_pd=pd.DataFrame(test,index=['N','CO','CA'])
    #transpose the dataframe
    new_pd_Tran=new_pd.T
    #select slice based on nitrogen chemical shift
    sele=new_pd_Tran[abs(new_pd_Tran['N']-(cs_nitrogen+n_shift))<=off_set]
    #plot the data using seaborn
    sns.scatterplot(data=sele,x='CA', y = 'CO', marker='x')
    for i in range(sele.shape[0]):
        plt.text(x=sele.CA[i]+0.1,y=sele.CO[i]+0.1,s=sele.index[i], 
          fontdict=dict(color= 'green',size=7))
        
def read_2Ddata(dir, data_type):
    if 'sparky' in data_type:
        dic, data = ng.sparky.read(dir)
    elif 'bruker' in data_type:
        dic, data = ng.bruker.read_pdata(dir)
    
    udic= ng.bruker.guess_udic(dic, data)
    uc_1=ng.fileiobase.uc_from_udic(udic, dim=1) 
    uc_0=ng.fileiobase.uc_from_udic(udic, dim=0) 
    #uc_13co_2=ng.fileiobase.uc_from_udic(udic_2, dim=1)
    #ppm_1 = uc_1.ppm_scale()
    #ppm_0 = uc_0.ppm_scale()

    return data, uc_0, uc_1

def read_3Ddata(dir, data_type):
    if 'sparky' in data_type:
        dic, data = ng.sparky.read(dir)
    elif 'bruker' in data_type:
        dic, data = ng.bruker.read_pdata(dir)
    
    udic= ng.bruker.guess_udic(dic, data)
    uc_2=ng.fileiobase.uc_from_udic(udic, dim=2) 
    uc_1=ng.fileiobase.uc_from_udic(udic, dim=1) 
    uc_0=ng.fileiobase.uc_from_udic(udic, dim=0) 
    #uc_13co_2=ng.fileiobase.uc_from_udic(udic_2, dim=1)
    #ppm_2 = uc_2.ppm_scale()
    #ppm_1 = uc_1.ppm_scale()
    #ppm_0 = uc_0.ppm_scale()

    return data, uc_0, uc_1, uc_2



        

def ncocx_connect(dir1, dir2, cs_nitrogen,df1, N_off1=-2, N_off2=-2,span=0.8,**kwargs):
    '''
    plot function to show generate peaks based on chemical shift list to experimental spectra
    params: dir1- NCOCX spectrum, dir2- 2D NCO spectrum
    cs_nitrogen: n chemical shift plane
    df1: df output from shiftx_to_df
    
    '''
    #read 2D, #_3 to label 2d data
    data_3, uc_15n_3, uc_13co_3 = read_2Ddata(dir2,'bruker')
    ppm_15n_3, ppm_13co_3=uc_15n_3.ppm_scale(), uc_13co_3.ppm_scale()
    #read 3D
    data_1, uc_13co, uc_15n, uc_13cx = read_3Ddata(dir1,'sparky')
    ppm_13co = uc_13co.ppm_scale()
    ppm_15n= uc_15n.ppm_scale()
    ppm_13cx= uc_13cx.ppm_scale()

    
    #calcualte noise in 3D
    noise = np.std(data_1[uc_13co(200, 'ppm'):uc_13co(190, 'ppm'), uc_15n(144,'ppm'),
                        uc_13cx(210, 'ppm'):uc_13cx(200, 'ppm')])*3 

    noise_3 = np.std(data_3[uc_15n_3(144,'ppm'):uc_15n_3(140,'ppm'),uc_13co_3(190, 'ppm'):uc_13co_3(185, 'ppm')])*3

    
    #countour information
    contour_start_s1 = noise
    contour_step_s1 = 1.3

    contour_start_s3 = noise_3
    contour_step_s3 = 1.3

    # contour_start_s3 = 3.0e5
    # contour_step_s3 = 1.20

    #colors_s1 = 'blue'
    #colors_s2 = 'green'
    #colors_s3 = 'red'

    cl_s1 = contour_start_s1 * contour_step_s1 ** np.arange(20)
    #cl_s2 = contour_start_s2 * contour_step_s2 ** np.arange(20)
    cl_s3 = contour_start_s3 * contour_step_s3 ** np.arange(20)
    #genereate 2D plane at cs_nitrogen 
    data_slice=data_1[:,uc_15n(cs_nitrogen,'ppm'),:]
    data_slice_Tran=np.array([*zip(*data_slice)])
    #data_slice_2=data_2[uc_15n_2(cs_nitrogen,'ppm'),:,:]
    data_3_Tran= np.array([*zip(*data_3)])
    #plot
    fig=plt.figure(figsize=(10,4))
    gs = gridspec.GridSpec(2, 2, width_ratios=[1, 3]) 
    ax1=plt.subplot(gs[0])
    ax1.contour(ppm_13co_3, ppm_15n_3, data_3, cl_s3)
    ax1.set_xlim([184,168])
    ax1.set_ylim([140,90])
    # mark nitrogen cs in NCO plane
    xslice = data_3[uc_15n_3(cs_nitrogen, 'ppm'), :]
    ax1.plot(ppm_13co_3, -xslice / 2.e5 + cs_nitrogen)
    line=np.array([0]*len(ppm_13co_3))
    ax1.plot(ppm_13co_3, line + cs_nitrogen, 'r--')
    ax1.xaxis.set_minor_locator(AutoMinorLocator())
    ax1.yaxis.set_minor_locator(AutoMinorLocator())


    ax2=plt.subplot(gs[3])
    ncocx_plot(df1,cs_nitrogen,n_shift=N_off1,off_set=span)
    #ncocx_e_plot(df2,cs_nitrogen,n_shift=N_off2,off_set=span)
    
    ax2.contour(ppm_13cx, ppm_13co, data_slice, cl_s1, colors='red',linewidths=0.5)
    ax2.set_ylim([184,168])
    ax2.set_xlim([72,14])
    ax2.xaxis.set_minor_locator(AutoMinorLocator())
    ax2.yaxis.set_minor_locator(AutoMinorLocator())
    #ax2.get_yaxis().set_visible(False)
    ax3=plt.subplot(gs[2])
    ax3.contour(ppm_13cx, ppm_13co, data_slice, cl_s1, colors='red',linewidths=0.5)
    ax3.set_xlim([184,168])
    ax3.set_ylim([184,168])
    ax3.xaxis.set_minor_locator(AutoMinorLocator())
    ax3.yaxis.set_minor_locator(AutoMinorLocator())
    
    fig.savefig('/Users/yunyao_1/Dropbox/KcsA_DNP/assign/assign_slice2/nitrogen_{}.pdf'.format(cs_nitrogen),dpi=300)












