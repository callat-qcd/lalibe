import os, sys
import make_secsource_h5_xml_ini as xml

flavs = ['UU','DD']
spins = ['up','dn']
corrs = ['proton','proton_np']
t_src = [3,6]
t_snks = {
    'proton'   :['<t_sep>4</t_sep>',''],
    'proton_np':['<t_sep>-4</t_sep>',''],
}
seqprop_tsnk = {
    '<t_sep>4</t_sep>':'tsep_4',
    '<t_sep>-4</t_sep>':'tsep_m4',
    '':'t_all'
}
t_sep   = {'proton':4  ,'proton_np':-4}
t_sep_s = {'proton':'4','proton_np':'m4'}

params = {'QUARK_SPIN':'FULL'}
for t in range(8):
    t0 = str(t)
    params['T0'] = t0
    for z in range(4):
        z0 = str(z)
        params['Z0'] = z0
        for y in range(4):
            y0 = str(y)
            params['Y0'] = y0
            for x in range(4):
                x0 = str(x)
                params['X0'] = x0
                prop = 'pt_prop_x'+x0+'_y'+y0+'_z'+z0+'_t'+t0
                src  = 'pt_src_x'+x0+'_y'+y0+'_z'+z0+'_t'+t0
                params['PROP'] = prop
                params['SOURCE'] = src
                if not os.path.exists('xml'):
                    os.makedirs('xml')
                f = open('xml/'+prop+'.ini.xml','w')
                f.write(xml.head)
                f.write(xml.random_pt_src % params)
                f.write(xml.prop % params)
                params['OBJECT_ID'] = prop
                params['OBJECT_TYPE'] = 'LatticePropagator'
                params['H5_FILE'] = 'all_pt_props.h5'
                params['H5_PATH'] = 'props'
                params['H5_NAME'] = prop
                f.write(xml.h5_write % params)
                f.write(xml.end)
                f.close()
