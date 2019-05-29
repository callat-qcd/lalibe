import os, sys
import make_seqsource_h5_xml_ini as xml

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


f = open('seqsource_h5.ini.xml','w')
f.write(xml.head)
params = {}
coherent = {}
multiadd = {}
for flav in flavs:
    coherent[flav] = dict()
    multiadd[flav] = dict()
    for spin in spins:
        coherent[flav][spin] = dict()
        multiadd[flav][spin] = dict()
        for corr in corrs:
            coherent[flav][spin][corr] = []
            multiadd[flav][spin][corr] = []

for t0 in t_src:
    params['T0'] = t0
    params['SOURCE'] = 'pt_source_'+str(t0)
    f.write(xml.pt_src % params)
    ''' write source '''
    params['OBJECT_ID'] = 'pt_source_'+str(t0)
    params['OBJECT_TYPE'] = 'LatticePropagator'
    params['H5_FILE'] = './test_propagator.h5'
    params['H5_PATH'] = 'pt'
    params['H5_NAME'] = 'pt_source_'+str(t0)
    f.write(xml.h5_write % params)
    params['PROP'] = 'PP_prop_'+str(t0)
    params['QUARK_SPIN'] = 'FULL'
    f.write(xml.prop % params)
    params['OBJECT_ID'] = 'PP_prop_'+str(t0)
    params['H5_NAME'] = 'PP_prop_'+str(t0)
    f.write(xml.h5_write % params)
    #params['POINT_SINK_PROP'] = 'PS_prop_'+str(t0)
    #params['SMEARED_SINK_PROP'] = 'SS_prop_'+str(t0)
    #f.write(xml.snk % params)
    #params['OBJECT_ID'] = 'SS_prop_'+str(t0)
    #params['H5_NAME'] = 'SS_prop_'+str(t0)
    #f.write(xml.h5_write % params)
    params['PROP'] = 'PP_prop_'+str(t0)
    params['H5_PATH'] = 'PP'
    f.write(xml.lalibe_2pt % params)
    for flav in flavs:
        params['FLAV'] = flav
        for spin in spins:
            params['SRC_SPIN'] = spin
            params['SNK_SPIN'] = spin
            for corr in corrs:
                params['PARTICLE'] = corr
                for ti,ts in enumerate(t_snks[corr]):
                    ''' SEQSOURCE '''
                    params['PROP']   = 'PP_prop_'+str(t0)
                    params['T_SINK'] = ts
                    seqsource  = 'PP_sink_'+corr+'_'+flav+'_'+spin+'_'
                    seqsource += spin+'_tsrc_'+str(t0)+'_'
                    seqsource += seqprop_tsnk[ts]
                    if seqprop_tsnk[ts] == 't_all':
                        coherent[flav][spin][corr].append(seqsource)
                    else:
                        multiadd[flav][spin][corr].append(seqsource)
                    params['SEQSOURCE'] = seqsource
                    f.write(xml.seqsource % params)
                    ''' WRITE '''
                    params['OBJECT_ID'] = seqsource
                    params['H5_FILE']   = './test_seqsource.h5'
                    params['H5_PATH']   = ''
                    params['H5_NAME']   = seqsource
                    f.write(xml.h5_write % params)
                    ''' SEQPROP '''
                    if seqprop_tsnk[ts] != 't_all':
                        seqprop = seqsource.replace('_sink_','_seqprop_')
                        params['PROP']   = seqprop
                        if corr == 'proton':
                            params['QUARK_SPIN'] = 'UPPER'
                        elif corr == 'proton_np':
                            params['QUARK_SPIN'] = 'LOWER'
                        params['SOURCE'] = seqsource
                        f.write(xml.prop % params)
                        params['OBJECT_ID'] = seqprop
                        params['H5_FILE']   = './test_seqprop.h5'
                        params['H5_PATH']   = ''
                        params['H5_NAME']   = seqprop
                        f.write(xml.h5_write % params)
                        ''' 3pt '''
                        params['PROP']    = 'PP_prop_'+str(t0)
                        params['SEQPROP'] = seqprop
                        params['H5_PATH'] = 'half'
                        f.write(xml.bar3pt % params)
                        f.write(xml.bar3pt_4d % params)

                        ''' DO FULL SPIN SOLVES ALSO '''
                        seqprop_full = seqprop+'_full'
                        params['QUARK_SPIN'] = 'FULL'
                        params['PROP']   = seqprop_full
                        f.write(xml.prop % params)
                        params['OBJECT_ID'] = seqprop_full
                        params['H5_FILE']   = './test_seqprop.h5'
                        params['H5_PATH']   = ''
                        params['H5_NAME']   = seqprop_full
                        f.write(xml.h5_write % params)
                        params['PROP']    = 'PP_prop_'+str(t0)
                        params['SEQPROP'] = seqprop_full
                        params['H5_PATH'] = 'full'
                        f.write(xml.bar3pt % params)
                        f.write(xml.bar3pt_4d % params)

for flav in flavs:
    for spin in spins:
        for corr in corrs:
            ''' COHERENT_SEQSOURCE '''
            params['T_SEP'] = t_sep[corr]
            params['J_DECAY'] = 3
            params['SINKS'] = ''
            for si,sink in enumerate(coherent[flav][spin][corr]):
                params['SINKS'] += '            <elem>'+sink+'</elem>'
                if si < len(coherent[flav][spin][corr])-1:
                    params['SINKS'] += '\n'
            params['RESULT_SINK'] =  'PP_coherent_sink_'+corr+'_'+flav+'_'
            params['RESULT_SINK'] += spin+'_tsep_'+t_sep_s[corr]
            f.write(xml.coherent_seqsource % params)
            params['OBJECT_ID']   = params['RESULT_SINK']
            params['OBJECT_TYPE'] = 'LatticePropagator'
            params['H5_FILE']     = './test_coherent_seqsource.h5'
            params['H5_PATH']     = ''
            params['H5_NAME']     = params['RESULT_SINK']
            f.write(xml.h5_write % params)

            ''' MULTI_PROP_ADD SEQSOURCE '''
            params['PROPS'] = ''
            params['PROP_WEIGHTS'] = ''
            for si,sink in enumerate(multiadd[flav][spin][corr]):
                params['PROPS'] += '            <elem>'+sink+'</elem>'
                params['PROP_WEIGHTS'] +='1'
                if si < len(multiadd[flav][spin][corr])-1:
                    params['PROPS'] += '\n'
                    params['PROP_WEIGHTS'] += ' '
            params['RESULT_PROP'] =  'PP_multiadd_sink_'+corr+'_'+flav+'_'
            params['RESULT_PROP'] += spin+'_tsep_'+t_sep_s[corr]
            f.write(xml.multi_prop_add % params)
            params['OBJECT_ID']   = params['RESULT_PROP']
            params['OBJECT_TYPE'] = 'LatticePropagator'
            params['H5_FILE']     = './test_coherent_seqsource.h5'
            params['H5_PATH']     = ''
            params['H5_NAME']     = params['RESULT_PROP']
            f.write(xml.h5_write % params)


f.write(xml.end)
f.close()
