from module.importmodule import *

class astro_query():
    def __init__(self, data, match_r=5):
        self.data    = data
        self.targetc = None
        self.match_r = match_r
        self.ps1_result = None
        self.wise_result = None
        self.mass_result = None
    def query_all(self):
        self.targetc = coord.SkyCoord(ra  = self.data['_RAJ2000'] * u.deg,
                                      dec = self.data['_DEJ2000'] * u.deg,
                                      frame = 'icrs')
        print('-------------------------------------------')
        print('Query ps1 data...')
        Vizier.ROW_LIMIT = 1
        Vizier.TIMEOUT = 300
        v_ps1 = Vizier(columns=["*", "+_r"], column_filters={"Nd": ">12"})
        result_ps = v_ps1.query_region(self.targetc, radius=self.match_r * u.arcsec, catalog="II/349/ps1")
        self.ps1_result = result_ps[0]['_q', '_r', 'Nd', 'gmag', 'e_gmag', 'rmag',
                                       'e_rmag', 'imag', 'e_imag', 'zmag', 'e_zmag', 'ymag', 'e_ymag']
        self.ps1_result = unique(self.ps1_result, keys='_q', keep='first')
        #if ps1_result['_q'].size > inputdata['_RAJ2000'].size:
        #    print('ps1 have duplications')
        self.ps1_result.rename_column('gmag', 'ps1_g')
        self.ps1_result.rename_column('e_gmag', 'ps1_eg')
        self.ps1_result.rename_column('rmag', 'ps1_r')
        self.ps1_result.rename_column('e_rmag', 'ps1_er')
        self.ps1_result.rename_column('imag', 'ps1_i')
        self.ps1_result.rename_column('e_imag', 'ps1_ei')
        self.ps1_result.rename_column('zmag', 'ps1_z')
        self.ps1_result.rename_column('e_zmag', 'ps1_ez')
        self.ps1_result.rename_column('ymag', 'ps1_y')
        self.ps1_result.rename_column('e_ymag', 'ps1_ey')
        #----------------------------------------------------------------

        #----------------------------------------------------------------
        print('Query allwise data...')
        # ALLWISE Query
        v_allwise = Vizier(columns=["*", "+_r", 'e_Jmag', 'e_Hmag', 'e_Kmag'])
        result_al = v_allwise.query_region(self.targetc, radius=self.match_r * u.arcsec, catalog="II/328/allwise")
        # wise_result = result_al[0]['_q','_r','W1mag','e_W1mag','W2mag','e_W2mag','W3mag','e_W3mag',
        #                        'W4mag','e_W4mag', 'Jmag','e_Jmag','Hmag','e_Hmag','Kmag','e_Kmag','pmRA','e_pmRA','pmDE','e_pmDE','qph']
        self.wise_result = result_al[0]['_q', '_r', 'W1mag', 'e_W1mag', 'W2mag', 'e_W2mag', 'W3mag', 'e_W3mag',
                                        'W4mag', 'e_W4mag']

        self.wise_result = unique(self.wise_result, keys='_q', keep='first')
        #if wise_result['_q'].size > inputdata['_RAJ2000'].size:
        #    print('WISE have duplications')
        self.wise_result.rename_column('W1mag', 'w1')
        self.wise_result.rename_column('e_W1mag', 'ew1')
        self.wise_result.rename_column('W2mag', 'w2')
        self.wise_result.rename_column('e_W2mag', 'ew2')
        self.wise_result.rename_column('W3mag', 'w3')
        self.wise_result.rename_column('e_W3mag', 'ew3')
        self.wise_result.rename_column('W4mag', 'w4')
        self.wise_result.rename_column('e_W4mag', 'ew4')

        #----------------------------------------------------------------

        #----------------------------------------------------------------
        print('Query 2MASS data...')
        # 2MASS Query
        v_2mass = Vizier(columns=["*", "+_r"])
        result_ma = v_2mass.query_region(self.targetc, radius=self.match_r * u.arcsec, catalog="II/246/out")
        self.mass_result = result_ma[0]['_q', '_r', 'Jmag', 'e_Jmag', 'Hmag', 'e_Hmag', 'Kmag', 'e_Kmag', 'Qflg']

        self.mass_result = unique(self.mass_result, keys='_q', keep='first')
        #if mass_result['_q'].size > inputdata['_RAJ2000'].size:
        #    print('2MASS have duplications')
        Vizier.ROW_LIMIT = 50
        self.mass_result.rename_column('Jmag', 'ma_j')
        self.mass_result.rename_column('e_Jmag', 'ma_ej')
        self.mass_result.rename_column('Hmag', 'ma_h')
        self.mass_result.rename_column('e_Hmag', 'ma_eh')
        self.mass_result.rename_column('Kmag', 'ma_k')
        self.mass_result.rename_column('e_Kmag', 'ma_ek')


        self.data = join(self.data, self.ps1_result, keys='_q', join_type='left')
        self.data = join(self.data, self.wise_result, keys='_q', join_type='left')
        self.data = join(self.data, self.mass_result, keys='_q', join_type='left')
        print('-------------------------------------------')
        return self.data
