from module.importmodule import *

def c_bt1(colors):
    '''
    colors: table object with color information (i-z, z-Y, Y-J.... W1-W2)
    result: numpy array of color difference anchor on J band (i-J, z-J, Y-J, J-J.... J-W2)
    '''
    arr=colors.as_array()
    gj=arr['g-r']+arr['r-i']+arr['i-z']+arr['z-J']
    rj=arr['r-i']+arr['i-z']+arr['z-J']
    ij=arr['i-z']+arr['z-J']
    zj=arr['z-J']
    jj=arr['z-J']-arr['z-J']
    hj=0.-arr['J-H']
    kj=0.-(arr['J-H']+arr['H-K'])
    w1j=0.-(arr['J-H']+arr['H-K']+arr['K-W1'])
    result=np.array([gj,rj,ij,zj,jj,hj,kj,w1j])
    return result

def m_bt1(mag,mer,color_correction):
    '''
    mag: mags of ONE star from i to W2 (i,z,y,J,H,K,W1,W2)
    mer: errors of mags above
    color_correction: color difference anchor on J band, 24 templates
    result: return numpy array, each star has 24 m_bt, crossponding to each template
    '''
    result=[]
    y=np.array(list(mag))
    yer=np.array(list(mer))
    yyer=np.power(yer,2)
    for i in range(color_correction.shape[1]):
        top = (y-color_correction[:,i])/yyer
        low = 1./yyer
        result.append(np.nansum(top)/np.nansum(low))
    return np.array(result)

def sp_typing1(mag,mer,ref_mag,color_correction):
    chi=[]
    y=np.array(list(mag))
    ok = np.where(y<900)
    for j in range(color_correction.shape[1]):
        chi.append(
            np.nansum(
                np.power(
                    (np.array(list(mag))[ok]-ref_mag[j]-color_correction[ok][:,j])/
                    np.array(list(mer))[ok],2)))

    chi=np.array(chi)
    type_index=np.nanargmin(chi)
    if (type_index != 0) & (type_index < 53):
        near_slop = ((chi[type_index-1]-chi[type_index]) +
                     (chi[type_index+1]-chi[type_index]) )
    else: near_slop = np.nan
    my_chi = chi[type_index]/len(ok[0])
    total_chi_B1M4 = chi/len(ok[0])

    #if len(ok[0]) < (2*(len(mag_cols)//3)): my_chi = np.nan
    if (type_index == 0) or (type_index == 53):              # If == B1,M4-> show nan
        my_chi = np.nan
    #if ((my_chi) > 1000): my_chi = np.nan
    return type_index, my_chi, near_slop


def c_bt2(colors):
    '''
    colors: table object with color information (g-r,i-z, z-Y, Y-J.... W1-W2)
    result: numpy array of color difference anchor on J band (i-J, z-J, Y-J, J-J.... J-W2)
    '''
    arr=colors.as_array()
    #PS1 grizy
    #gj=arr['g-r']+arr['r-i']+arr['i-z']+arr['z-y']+arr['y-J']
    #rj=arr['r-i']+arr['i-z']+arr['z-y']+arr['y-J']
    ij=arr['i-z']+arr['z-y']+arr['y-J']
    zj=arr['z-y']+arr['y-J']
    yj=arr['y-J']
    jj=arr['y-J']-arr['y-J']
    hj=0.-arr['J-H']
    kj=0.-(arr['J-H']+arr['H-Ks'])
    w1j=0.-(arr['J-H']+arr['H-Ks']+arr['Ks-W1'])
    w2j=0.-(arr['J-H']+arr['H-Ks']+arr['Ks-W1']+arr['W1-W2'])
    result=np.array([ij,zj,yj,jj,hj,kj,w1j,w2j]) # result=np.array([ij,zj,yj,jj,hj,kj,w1j,w2j])
    #result=np.array([gj,rj,ij,zj,yj,jj,hj,kj])
    return result

def m_bt2(mag,mer,color_correction2):
    '''
    mag: mags of ONE star from i to W2 (i,z,y,J,H,K,W1,W2)
    mer: errors of mags above
    color_correction: color difference anchor on J band, 24 templates
    result: return numpy array, each star has 24 m_bt, crossponding to each template
    '''
    result=[]
    y=np.array(list(mag))
    yer=np.array(list(mer))
    yyer=np.power(yer,2)
    for i in range(color_correction2.shape[1]):
        top = (y-color_correction2[:,i])/yyer
        low = 1./yyer
        result.append(np.nansum(top)/np.nansum(low))
    return np.array(result)

def sp_typing2(mag,mer,ref_mag,color_correction):
    chi=[]
    y=np.array(list(mag))
    ok = np.where(y<900)
    for j in range(color_correction.shape[1]):
        chi.append(
            np.nansum(
                np.power(
                    (np.array(list(mag))-ref_mag[j]-color_correction[:,j])/
                    np.array(list(mer)),2)))
    chi=np.array(chi)
    type_index=np.nanargmin(chi)
    if (type_index != 0) & (type_index < 22):
        near_slop = ((chi[type_index-1]-chi[type_index]) +
                     (chi[type_index+1]-chi[type_index]) )
    else: near_slop = np.nan
    my_chi = chi[type_index]/len(ok[0])
    total_chi = chi/len(ok[0])

    #if len(ok[0]) < (2*(len(mag_cols)//3)): my_chi = np.nan
    if (type_index == 0) or (type_index == 22):              #If == M0,T2 -> show nan
        my_chi = np.nan
    #if ((my_chi) > 400): my_chi = np.nan
    return type_index, my_chi, near_slop
