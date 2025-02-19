from module.importmodule import *

def uktotm(ps1_row):
    uk_mags=dict(zip(ps1_row.colnames,list(ps1_row)))
    tm_j=uk_mags['uk_j'] - 0.01 + 0.073*(uk_mags['uk_j'] - uk_mags['uk_h'])
    tm_h=uk_mags['uk_h'] - 0.01 + 0.069*(uk_mags['uk_h'] - uk_mags['uk_k1'])
    tm_k=uk_mags['uk_k1'] - 0.073*(uk_mags['uk_h'] - uk_mags['uk_k1'])
    return np.array([tm_j,tm_h,tm_k])

def ps1tosdss(ps1_row):
    #transforming ps1 mags to SDSS mags (only for griz bands)
    #reference: Finkbeiner et al. 2016
    #http://iopscience.iop.org/article/10.3847/0004-637X/822/2/66/pdf
    ps1mag=dict(zip(ps1_row.colnames,list(ps1_row)))
    gc=np.array([-0.01808, -0.13595,  0.01941, -0.00183])
    rc=np.array([-0.01836, -0.03577,  0.02612, -0.00558])
    ic=np.array([ 0.01170, -0.00400,  0.00066, -0.00058])
    zc=np.array([-0.01062,  0.07529, -0.03592,  0.00890])
    x=ps1mag['ps1_g']-ps1mag['ps1_i']
    gsd=ps1mag['ps1_g']-(gc[0]+gc[1]*x+gc[2]*np.power(x,2)+gc[3]*np.power(x,3))
    rsd=ps1mag['ps1_r']-(rc[0]+rc[1]*x+rc[2]*np.power(x,2)+rc[3]*np.power(x,3))
    isd=ps1mag['ps1_i']-(ic[0]+ic[1]*x+ic[2]*np.power(x,2)+ic[3]*np.power(x,3))
    zsd=ps1mag['ps1_z']-(zc[0]+zc[1]*x+zc[2]*np.power(x,2)+zc[3]*np.power(x,3))
    return np.array([ps1_row['_q'],gsd,rsd,isd,zsd])

sequence = 0
SpTy_lib = {}
for alpha in ['O','B','A','F','G','K','M','L','T','C','W']:
    if alpha == 'O':
        for num in np.arange(0, 10):
            SpTy_lib_1 = {sequence:alpha+'%i'%num}
            SpTy_lib.update(SpTy_lib_1)
            sequence+=1
    elif  alpha == 'B':
        for num in np.arange(0, 10):
            SpTy_lib_1 = {sequence:alpha+'%i'%num}
            SpTy_lib.update(SpTy_lib_1)
            sequence+=1
    elif  alpha == 'A':
        for num in np.arange(0, 10):
            SpTy_lib_1 = {sequence:alpha+'%i'%num}
            SpTy_lib.update(SpTy_lib_1)
            sequence+=1
    elif  alpha == 'F':
        for num in np.arange(0, 10):
            SpTy_lib_1 = {sequence:alpha+'%i'%num}
            SpTy_lib.update(SpTy_lib_1)
            sequence+=1
    elif  alpha == 'G':
        for num in np.arange(0, 10):
            SpTy_lib_1 = {sequence:alpha+'%i'%num}
            SpTy_lib.update(SpTy_lib_1)
            sequence+=1
    elif  alpha == 'K':
        for num in np.arange(0, 10):
            SpTy_lib_1 = {sequence:alpha+'%i'%num}
            SpTy_lib.update(SpTy_lib_1)
            sequence+=1
    elif  alpha == 'M':
        for num in np.arange(0, 10):
            SpTy_lib_1 = {sequence:alpha+'%i'%num}
            SpTy_lib.update(SpTy_lib_1)
            sequence+=1
    elif  alpha == 'L':
        for num in np.arange(0, 10):
            SpTy_lib_1 = {sequence:alpha+'%i'%num}
            SpTy_lib.update(SpTy_lib_1)
            sequence+=1
    elif  alpha == 'T':
        for num in np.arange(0, 10):
            SpTy_lib_1 = {sequence:alpha+'%i'%num}
            SpTy_lib.update(SpTy_lib_1)
            sequence+=1
    elif  alpha == 'C':
        for num in np.arange(0, 10):
            SpTy_lib_1 = {sequence:alpha+'V'}
            SpTy_lib.update(SpTy_lib_1)
            sequence+=1
    elif  alpha == 'W':
        for num in np.arange(0, 10):
            SpTy_lib_1 = {sequence:alpha+'D'}
            SpTy_lib.update(SpTy_lib_1)
            sequence+=1
inver_SpTy_lib = dict([[v,k] for k,v in SpTy_lib.items()])

def SpTy(target, num=0):
    if num==0:
        # SpTy to num
        result = inver_SpTy_lib[target]
    else:
        # number to SpTy
        result = SpTy_lib[target]
    return result
