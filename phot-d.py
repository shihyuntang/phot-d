#-----------------------------------------------------------------------
# phot-d v0.2
# 2018/10/16: initiated by shih-yun tang, at NCU, Taiwan
# 2025/02/19: update code for Github release by shih-yun tang, at Rice, Houston, USA
#-----------------------------------------------------------------------
from module.importmodule import *
import module.mag_transfer as mag_transfer
import module.estimation   as estimation
import module.query        as query

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

class Indata():
    '''
    Input data
    '''
    def __init__(self, in_name):
        self.format = None
        self.in_name = in_name
        self.result = None
    def open(self):
        self.format = os.path.splitext(self.in_name)[-1][1:]
        self.result = Table.read('./estimate_data/' + self.in_name, format=self.format)
        # Add tracking sequence
        self.result['_q'] = np.arange(len(self.result), dtype='int') + 1
        return self.result

class Outdata():
    '''
    Output data
    '''
    def __init__(self, in_name, out_name, inputdata, resultdata):
        self.format     = None
        self.in_name    = in_name
        self.inputdata  = inputdata
        self.resultdata = resultdata
        self.yy = strftime("%Y", localtime())
        self.dd = strftime(self.yy + "%m%d-%H%M", localtime())
        if out_name == '':
            self.out_name = '_phot-d_' + in_name
        else:
            self.out_name = out_name
    def write(self):
        if os.path.splitext(self.out_name)[-1][1:] == '':
            self.format = os.path.splitext(self.in_name)[-1][1:]
        else:
            self.format = os.path.splitext(self.out_name)[-1][1:]
        outputcolumn = ['ps1_g', 'ps1_r', 'ps1_i', 'ps1_z', 'ps1_y', 'ma_j', 'ma_h',
                        'ma_k', 'w1', 'w2',
                        'ps1_eg', 'ps1_er', 'ps1_ei', 'ps1_ez', 'ps1_ey', 'ma_ej', 'ma_eh',
                        'ma_ek', 'ew1', 'ew2',
                        'phot_sp1', 'phot_chi1', 'phot_sp2', 'phot_chi2',
                        'phot_sp', 'phot_nsp', 'phot_chi',
                        'phot_dist', 'phot_dist_err']
        outputdata = hstack((self.inputdata, self.resultdata[outputcolumn]))
        outputdata.write('result/' + self.dd + self.out_name, format=self.format)

class Incolor_tab():
    '''
    Read in the color table
    '''
    def __init__(self, in_name):
        self.in_name = in_name
        self.format = None
        self.result = None
    def open(self):
        self.format = os.path.splitext(self.in_name)[-1][1:]
        if self.format == 'txt':
            self.result = Table.read('./color_table/' + self.in_name, format='ascii')
        else:
            self.result = Table.read('./color_table/' + self.in_name, format=self.format)
        return self.result


def SpTy_estimation1(color_correction, mag_cols, mer_cols, magmer):
    ref_mag            = estimation.m_bt1(magmer[mag_cols[1:]], magmer[mer_cols[1:]], color_correction)
    idx, mychi, near_S = estimation.sp_typing1(magmer[mag_cols[1:]], magmer[mer_cols[1:]],ref_mag, color_correction)
    return np.array([magmer['_q'],idx, mychi, near_S])

def SpTy_estimation2(color_correction, mag_cols, mer_cols, magmer):
    ref_mag            = estimation.m_bt2(magmer[mag_cols[1:]], magmer[mer_cols[1:]], color_correction)
    idx, mychi, near_S = estimation.sp_typing2(magmer[mag_cols[1:]], magmer[mer_cols[1:]],ref_mag, color_correction)
    return np.array([magmer['_q'],idx, mychi, near_S])

def Distance_estimation(absM_sque, absmags, d):
    seque = d['phot_nsp']
    doST  = np.where(seque == np.array(absM_sque))[0]
    absg  = absmags['g'][doST].data[0]
    absr  = absmags['r'][doST].data[0]
    absi  = absmags['i'][doST].data[0]
    absz  = absmags['z'][doST].data[0]
    absy  = absmags['y'][doST].data[0]
    absj  = absmags['J'][doST].data[0]
    absk  = absmags['H'][doST].data[0]
    absh  = absmags['K'][doST].data[0]
    absw1 = absmags['W1'][doST].data[0]
    absw2 = absmags['W2'][doST].data[0]
    ab    = np.array([absg, absr, absi, absz, absy, absj, absh, absk, absw1, absw2])
    obs   = np.array(list([d['ps1_g'], d['ps1_r'], d['ps1_i'], d['ps1_z'],
                         d['ps1_y'], d['ma_j'], d['ma_h'], d['ma_k'], d['w1'], d['w2']]))
    ddtemp  = np.power(10, 1. + (obs - ab) / 5.)
    dd      = np.nanmedian(ddtemp)
    dder    = np.nanstd(ddtemp)
    max_dif = np.nanmax(ddtemp) - np.nanmin(ddtemp)
    return np.array([d['_q'], dd, dder, max_dif])

def abs_tab_999(column):
    for i in range(len(column)):
        if column[i] > 900:
            column[i] += np.nan
    return column

def tcd_drew(tcdtab, data):
    plt.style.use('seaborn-v0_8-deep')
    tcd = plt.figure(figsize=(10,3), facecolor='white', dpi=200)
    gs1 = gridspec.GridSpec(1,1)
    gs2 = gridspec.GridSpec(1,1)
    gs3 = gridspec.GridSpec(1,1)
    gs1.update(left=0.02, right=0.26,  bottom=0.05, top=0.9)
    gs2.update(left=0.32, right=0.56, bottom=0.05, top=0.9)
    gs3.update(left=0.62, right=0.90, bottom=0.05, top=0.9)
    ax1 = plt.subplot(gs1[0,0]) ; ax2 = plt.subplot(gs2[0,0]) ; ax3 = plt.subplot(gs3[0,0])


    a = ax1.scatter(
        tcdtab['PS1zmag']-tcdtab['PS1ymag'], tcdtab['PS1ymag']-tcdtab['Jmag'], 
        c=tcdtab['Type'], cmap='rainbow'
        )
    b = ax2.scatter(
        tcdtab['PS1imag']-tcdtab['PS1zmag'], tcdtab['PS1zmag']-tcdtab['PS1ymag'], 
        c=tcdtab['Type'], cmap='rainbow'
        )
    c = ax3.scatter(
        tcdtab['PS1zmag']-tcdtab['PS1ymag'], tcdtab['W1mag']-tcdtab['W2mag'], 
        c=tcdtab['Type'], cmap='rainbow'
        )

    ax1.plot(data['ps1_z']-data['ps1_y'], data['ps1_y']-data['ma_j'], '+k', mew=2, ms=15)
    ax2.plot(data['ps1_i']-data['ps1_z'], data['ps1_z']-data['ps1_y'], '+k', mew=2, ms=15)
    ax3.plot(data['ps1_z']-data['ps1_y'], data['w1']-data['w2'],      '+k', mew=2, ms=15)

    d = plt.colorbar(c, ax=ax3)
    d.set_label(label='SpTy', size='small', )
    d.ax.tick_params(labelsize='small')
    d.ax.set_yticklabels(['M0','M5','L0','L5','T0','T5'])

    ax1.set_xlabel('z$_{ps1}-$y$_{ps1}$', size='small', style='normal', family='serif')
    ax1.set_ylabel('y$_{ps1}-$J$_{2M}$',  size='small', style='normal', family='serif')

    ax2.set_xlabel('i$_{ps1}-$z$_{ps1}$', size='small', style='normal', family='serif')
    ax2.set_ylabel('z$_{ps1}-$y$_{2M}$',  size='small', style='normal', family='serif')

    ax3.set_xlabel('z$_{ps1}-$y$_{ps1}$', size='small', style='normal', family='serif')
    ax3.set_ylabel('W1$-$W2',             size='small', style='normal', family='serif')

    ax1.tick_params(axis='both', colors='k', direction='in', right=True, top=True, labelsize='small')
    ax2.tick_params(axis='both', colors='k', direction='in', right=True, top=True, labelsize='small')
    ax3.tick_params(axis='both', colors='k', direction='in', right=True, top=True, labelsize='small')
    
    tcd.savefig('./result/tcd/tcd_'+str(data['_q'])+'.png', format='png', bbox_inches='tight')


#-----------------------------------------------------------------------
timer_start = timeit.default_timer()
#-----------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
                                     prog        = 'phot-d',
                                     description = '''
                                     An spectral type & distance estimator base on the muti-band photometry. \n
                                     Please place the target list into the *estimate_data* folder, \n
                                     and check if the header is named as requested. \n
                                     The requested header example can be found in the *estimate_data* folder. \n
                                     We only take target list with fits/csv/ascii format in current version. \n
                                     ''',
                                     epilog = "Contact the author through sytang@g.ncu.edu.tw")
    parser.add_argument(    "TableName",                    action="store",
                        help="target list name that you want to estimate. !no need to add the dir!",        type=str)
    parser.add_argument('-o',     dest="Fout",              action="store",
                        help="output filename. the input file name will be used if not specified",
                        type=str, default='')
    parser.add_argument('-c',     dest="CPUuse",            action="store",
                        help="# of cpu to use, default is 1/2 # avalible (you have %i cpu avaliable)"%(mp.cpu_count()),
                        type=int, default=int(mp.cpu_count()//2) )
    parser.add_argument('-q',     dest="auto_query",        action="store_true",
                        help="If set, will do astro-query for PanSTARRs, ALLWISE, and 2MASS (Warning!! not good for large dataset!)")
    parser.add_argument('-qr',    dest="query_r",            action="store",
                        help="astro-query searching radius, default = 5 arcsec",
                        type=int, default = 5 )
    parser.add_argument('-tcd',     dest="TCD_drew",         action="store_true",
                        help="If set, will plot the two-color diagram")
    parser.add_argument('--version',            action='version',  version='%(prog)s 0.2')
    args = parser.parse_args()

    print('\n-------------------------------------------')
    print('Program started with using %i cpus\n'%args.CPUuse)
    #-----------------------------------------------------------------------
    # Readin data
    inputdata = Indata(args.TableName).open()
    data = inputdata[0:]

    #----------------------------------------------------------------
    # Astroquery
    #----------------------------------------------------------------
    if args.auto_query:
        data = query.astro_query(data, args.query_r).query_all()

    # ************************************************************************ #
    #----------------------------------------------------------------
    # Drew Two-color diagram
    #----------------------------------------------------------------
    if args.TCD_drew:

        print('Start to plot the two-color diagram...')
        tcdtab = Table.read('./color_table/best18_ps1_MLTcolor.fits', format='fits')
        tcdtab = tcdtab[ 
            (tcdtab['Binary'].filled('999')=='999') &  (tcdtab['Young']=='    ')
            ] # clean out binary & YSO

        pool = mp.Pool(processes=args.CPUuse)
        func = partial(tcd_drew, tcdtab)
        pool.map(func, data)
        pool.close()
        pool.join()

    # ************************************************************************ #
    #----------------------------------------------------------------
    # Estimaiton B1 -- M4
    #----------------------------------------------------------------
    # Read in color table1
    color_cols=['g-r','r-i','i-z','z-J','J-H','H-K','K-W1']
    colors = Incolor_tab('B1M4_colors.csv').open()
    sp = colors['Sp']
    colors = colors[color_cols]
    color_correction = estimation.c_bt1(colors)

    # take out photometry needed
    mag_cols=['_q','ps1_g','ps1_r','ps1_i','ps1_z','ma_j','ma_h','ma_k','w1']
    mer_cols=['_q','ps1_eg','ps1_er','ps1_ei','ps1_ez','ma_ej','ma_eh','ma_ek','ew1']
    mags=data[mag_cols]
    mers=data[mer_cols]

    # photometry transformation
    pool = mp.Pool(processes=args.CPUuse)
    sdssmags = pool.map( mag_transfer.ps1tosdss, data['_q','ps1_g','ps1_r','ps1_i','ps1_z'])
    pool.close()
    pool.join()

    sdssmags = np.array(sdssmags)
    sdssmags = sdssmags[np.argsort(sdssmags[:,0])]

    mags[mag_cols[1]] = np.around(sdssmags[:,1], decimals=3)
    mags[mag_cols[2]] = np.around(sdssmags[:,2], decimals=3)
    mags[mag_cols[3]] = np.around(sdssmags[:,3], decimals=3)
    mags[mag_cols[4]] = np.around(sdssmags[:,4], decimals=3)
    magmer = join(mags, mers,  keys='_q')
    #-----------------------------------------------
    #-----------------------------------------------
    print('SpTy estimation Started (1/2)')
    pool = mp.Pool(processes=args.CPUuse)
    func = partial(SpTy_estimation1, color_correction, mag_cols, mer_cols)
    est_rest1 = pool.map(func, magmer)
    pool.close()
    pool.join()
    #print('\r', '%1.2f%% completed (1/3)' % (100. * (i / len(mags))), end='')

    est_rest1 = np.array(est_rest1)
    est_rest1 = est_rest1[np.argsort(est_rest1[:,0])]
    temp_sp = [sp[int(i)] for i in est_rest1[:,1]]
    sp_col    = Column(data = np.array(temp_sp),        name='phot_sp1',   dtype ='5U')
    nsp_col   = Column(data = np.array(est_rest1[:,1]), name='phot_nsp1',  )
    chi_col   = Column(data = np.array(est_rest1[:,2]), name='phot_chi1',  format='1.2f')
    Nears_col = Column(data = np.array(est_rest1[:,3]), name='phot_slop1', format='1.3f')
    data.add_columns([sp_col, nsp_col, chi_col, Nears_col])

    #----------------------------------------------------------------
    # Estimaiton M0 -- T2
    #----------------------------------------------------------------

    # Read in color table1
    color_cols=['i-z','z-y','y-J','J-H','H-Ks','Ks-W1','W1-W2']
    colors = Incolor_tab('M_colors.txt').open()
    sp     = colors[0:23]['SpT']
    colors = colors[0:23][color_cols]
    color_correction = estimation.c_bt2(colors)


    # take out photometry needed
    mag_cols=['_q','ps1_i','ps1_z','ps1_y','ma_j','ma_h','ma_k','w1','w2']
    mer_cols=['_q','ps1_ei','ps1_ez','ps1_ey','ma_ej','ma_eh','ma_ek','ew1','ew2']
    magmer = data[ mag_cols+mer_cols[1:] ]

    #-----------------------------------------------
    #-----------------------------------------------

    print('SpTy estimation Started (2/2)')
    pool = mp.Pool(processes=args.CPUuse)
    func = partial(SpTy_estimation2, color_correction, mag_cols, mer_cols)
    est_rest2 = pool.map(func, magmer)
    pool.close()
    pool.join()
    #print('\r', '%1.2f%% completed (1/3)' % (100. * (i / len(mags))), end='')

    est_rest2 = np.array(est_rest2)
    est_rest2 = est_rest2[np.argsort(est_rest2[:,0])]
    temp_sp = [sp[int(i)] for i in est_rest2[:,1]]
    sp_col    = Column(data = np.array(temp_sp),        name='phot_sp2',   dtype ='5U')
    nsp_col   = Column(data = np.array(est_rest2[:,1]), name='phot_nsp2',  )
    chi_col   = Column(data = np.array(est_rest2[:,2]), name='phot_chi2',  format='1.2f')
    Nears_col = Column(data = np.array(est_rest2[:,3]), name='phot_slop2', format='1.3f')
    data.add_columns([sp_col, nsp_col, chi_col, Nears_col])

    print('-------------------------------------------')
    print('  SpTy estimation DONE')
    print('  Merging the results..')

    #----------------------------------------------------------------
    # Merge the result phot_sp1, phot_sp2 >> phot_sp
    #----------------------------------------------------------------

    # merge phot_sp result
    data['phot_sp']  = data['phot_sp1']              # Final result "phot_sp" >> 1st,
    data['phot_chi'] = data['phot_chi1']             #       put all targets in
    data['phot_slop'] = data['phot_slop1']           #       phot_sp1's that have results

    tar_nan    = np.isnan(np.array(data['phot_chi'])) # find targets in phot_sp that
    tar_ok_put = np.nonzero(tar_nan)                  #      have no estimated result

    for i in tar_ok_put[0]:                           # Final result >> 2st, put phot_sp2's results
        data['phot_sp'][i]   = data['phot_sp2'][i]    #       in for those have no estimated results
        data['phot_chi'][i]  = data['phot_chi2'][i]   #       in phot_sp (phot_sp1 have no results)
        data['phot_slop'][i] = data['phot_slop2'][i]


# Deal with targets that have estimates results in both sp1 & sp2
#      by finding which have the most significant low chi2
    tar_not_nan  = np.logical_not(np.isnan(np.array(data['phot_chi'])))  # Find targets have results in
    tar_care_put = np.nonzero(tar_not_nan)                               #      phot_sp now
    for i in tar_care_put[0]:
        if np.isnan(data['phot_chi2'][i]) == 0:                          # If the same target have phot_sp2 result:
            if ( (data['phot_chi2'][i]< data['phot_chi1'][i]) or
                 ( (data['phot_slop2'][i]/data['phot_chi2'][i]) > (data['phot_slop1'][i]/data['phot_chi1'][i])) ):
                data['phot_sp'][i]   = data['phot_sp2'][i]
                data['phot_chi'][i]  = data['phot_chi2'][i]
                data['phot_slop'][i] = data['phot_slop2'][i]

    # making phot_nsp
    spnum = [ mag_transfer.SpTy(data['phot_sp'][i], 0) for i in range(len(data['phot_sp'])) ]
    data['phot_nsp'] = spnum

    # Let output to be nan if no estimated result
    tar_nan    = np.isnan(np.array(data['phot_chi']))
    tar_ok_put = np.nonzero(tar_nan)
    for i in tar_ok_put[0]:
        data['phot_sp'][i]   = np.nan
        data['phot_chi'][i]  = np.nan
        data['phot_slop'][i] = np.nan


    canot = sum(np.isnan(np.array(data['phot_chi']))) / len(data['phot_chi'])
    persent = canot * 100.
    print('-------------------------------------------')
    print('There are %1.2f %% can not be classicfied' % persent)
    print('-------------------------------------------')

    # print(data['phot_sp1', 'phot_chi1', 'phot_sp2', 'phot_chi2','phot_sp', 'phot_chi' ])
    #----------------------------------------------------------------
    # Distance estimation (1)
    #----------------------------------------------------------------

    wherenan = np.isnan(np.array(data['phot_chi']))
    okdata   = np.logical_not(wherenan)
    Tar_distance_E = data[np.nonzero(okdata)]

    # reading abs mag table
    absmags = Incolor_tab('absmag.csv').open()

    absmags['g'] = abs_tab_999(absmags['g'])
    absmags['r'] = abs_tab_999(absmags['r'])
    absmags['i'] = abs_tab_999(absmags['i'])
    absmags['z'] = abs_tab_999(absmags['z'])
    absmags['y'] = abs_tab_999(absmags['y'])
    absmags['J'] = abs_tab_999(absmags['J'])
    absmags['H'] = abs_tab_999(absmags['H'])
    absmags['K'] = abs_tab_999(absmags['K'])
    absmags['W1'] = abs_tab_999(absmags['W1'])
    absmags['W2'] = abs_tab_999(absmags['W2'])
    Tar_distance_E['ma_j'].format='1.3f'
    Tar_distance_E['ma_h'].format='1.3f'
    Tar_distance_E['ma_k'].format='1.3f'

    Tar_distance_E['sd_g'] = np.zeros(len(Tar_distance_E['ma_j']))
    Tar_distance_E['sd_r'] = np.zeros(len(Tar_distance_E['ma_j']))
    Tar_distance_E['sd_i'] = np.zeros(len(Tar_distance_E['ma_j']))
    Tar_distance_E['sd_z'] = np.zeros(len(Tar_distance_E['ma_j']))

    mag_cols=['_q','ps1_g','ps1_r','ps1_i','ps1_z']
    temp_mags=Tar_distance_E[mag_cols]

    #-----------------------------------------------
    #-----------------------------------------------

    # photometry transfomation
    pool = mp.Pool(processes=args.CPUuse)
    sdssmags = pool.map( mag_transfer.ps1tosdss, Tar_distance_E['_q','ps1_g','ps1_r','ps1_i','ps1_z'])
    pool.close()
    pool.join()

    sdssmags = np.array(sdssmags)
    sdssmags = sdssmags[np.argsort(sdssmags[:,0])]

    temp_mags['ps1_g'] = np.around(sdssmags[:,1], decimals=3)
    temp_mags['ps1_r'] = np.around(sdssmags[:,2], decimals=3)
    temp_mags['ps1_i'] = np.around(sdssmags[:,3], decimals=3)
    temp_mags['ps1_z'] = np.around(sdssmags[:,4], decimals=3)

    for x,y in zip(Tar_distance_E,temp_mags):
        if x['sd_g']>900:
            x['sd_g'] = y['ps1_g']
            x['sd_r'] = y['ps1_r']
            x['sd_i'] = y['ps1_i']
            x['sd_z'] = y['ps1_z']

    absM_sque = [mag_transfer.SpTy(i, 0) for i in absmags['SpTy']]

    #----------------------------------------------------------------
    # Distance estimation (2)
    #----------------------------------------------------------------

    print('Start the distance estimation...')

    pool = mp.Pool(processes=args.CPUuse)
    func = partial(Distance_estimation, absM_sque, absmags)
    est_rest3 = pool.map(func, Tar_distance_E)
    pool.close()
    pool.join()

    est_rest3 = np.array(est_rest3)
    est_rest3 = est_rest3[np.argsort(est_rest3[:,0])]

    phot_dist     = Column(data=np.array(est_rest3[:,1]), name='phot_dist',     format='1.2f')
    phot_dist_err = Column(data=np.array(est_rest3[:,2]), name='phot_dist_err', format='1.2f')
    phot_dist_max = Column(data=np.array(est_rest3[:,3]), name='phot_dist_max', format='1.2f')
    Tar_distance_E.add_columns([phot_dist, phot_dist_err, phot_dist_max])

    print('  Distance estimation DONE')


    #----------------------------------------------------------------
    # merging phot_d result
    #----------------------------------------------------------------
    # Combind Target, B_M7_data, A_M7_data
    data['phot_dist']     = np.zeros(len(data['_q']))
    data['phot_dist_err'] = np.zeros(len(data['_q']))
    data['phot_dist_max'] = np.zeros(len(data['_q']))

    for i in range(len(Tar_distance_E['_q'])):
        targetno = np.where(Tar_distance_E['_q'][i] == data['_q'])
        data['phot_dist'][targetno]     = Tar_distance_E['phot_dist'][i]
        data['phot_dist_err'][targetno] = Tar_distance_E['phot_dist_err'][i]
        data['phot_dist_max'][targetno] = Tar_distance_E['phot_dist_max'][i]

    for i in range(len(data['phot_dist'])):
        if data['phot_dist'][i] == 0:
            data['phot_dist'][i]     = np.nan
            data['phot_dist_err'][i] = np.nan
            data['phot_dist_max'][i] = np.nan
            data['phot_nsp'][i]      = 99

    data['phot_dist'].format     = '1.2f'
    data['phot_dist_err'].format = '1.2f'
    data['phot_dist_max'].format = '1.2f'

    print('Outputting the results...')

    #----------------------------------------------------------------
    # Outputting the result
    #----------------------------------------------------------------
    Outdata(args.TableName, args.Fout, inputdata, data).write()

    timer_end = timeit.default_timer()
    spendtime = timer_end - timer_start
    print('-------------------------------------------')
    print('Finished\nTotal riun time:', strftime("%H:%M:%S", gmtime(spendtime)))

if __name__ == '__main__':
    main()
