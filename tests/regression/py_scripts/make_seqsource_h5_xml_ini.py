head = '''<?xml version="1.0"?>
<lalibe>
<annotation>
;
; Test input file for lalibe
;
</annotation>
<Param>
<InlineMeasurements>

'''

src = '''<elem>
<Name>MAKE_SOURCE</Name>
<Frequency>1</Frequency>
<Param>
    <version>6</version>
    <Source>
        <version>2</version>
        <SourceType>SHELL_SOURCE</SourceType>
        <j_decay>3</j_decay>
        <t_srce>0 0 0 %(T0)s</t_srce>
        <SmearingParam>
            <wvf_kind>GAUGE_INV_GAUSSIAN</wvf_kind>
            <wvf_param>2.0</wvf_param>
            <wvfIntPar>5</wvfIntPar>
            <no_smear_dir>3</no_smear_dir>
        </SmearingParam>
    </Source>
</Param>
<NamedObject>
    <gauge_id>default_gauge_field</gauge_id>
    <source_id>%(SOURCE)s</source_id>
</NamedObject>
</elem>

'''

random_pt_src = '''<elem>
<Name>MAKE_SOURCE</Name>
<Frequency>1</Frequency>
<Param>
    <version>6</version>
    <Source>
        <version>2</version>
        <SourceType>POINT_SOURCE</SourceType>
        <j_decay>3</j_decay>
        <t_srce>%(X0)s %(Y0)s %(Z0)s %(T0)s</t_srce>
    </Source>
</Param>
<NamedObject>
    <gauge_id>default_gauge_field</gauge_id>
    <source_id>%(SOURCE)s</source_id>
</NamedObject>
</elem>

'''


pt_src = '''<elem>
<Name>MAKE_SOURCE</Name>
<Frequency>1</Frequency>
<Param>
    <version>6</version>
    <Source>
        <version>2</version>
        <SourceType>POINT_SOURCE</SourceType>
        <j_decay>3</j_decay>
        <t_srce>0 0 0 %(T0)s</t_srce>
    </Source>
</Param>
<NamedObject>
    <gauge_id>default_gauge_field</gauge_id>
    <source_id>%(SOURCE)s</source_id>
</NamedObject>
</elem>

'''

prop='''<elem>
<Name>PROPAGATOR</Name>
<Frequency>1</Frequency>
<Param>
    <version>10</version>
    <quarkSpinType>%(QUARK_SPIN)s</quarkSpinType>
    <obsvP>false</obsvP>
    <numRetries>1</numRetries>
    <FermionAction>
        <FermAct>UNPRECONDITIONED_CLOVER</FermAct>
        <Mass>0.5</Mass>
        <clovCoeff>1.17</clovCoeff>
        <FermionBC>
            <FermBC>SIMPLE_FERMBC</FermBC>
            <boundary>1 1 1 -1</boundary>
        </FermionBC>
    </FermionAction>
    <InvertParam>
        <invType>CG_INVERTER</invType>
        <RsdCG>1.0e-9</RsdCG>
        <MaxCG>1000</MaxCG>
    </InvertParam>
</Param>
<NamedObject>
    <gauge_id>default_gauge_field</gauge_id>
    <source_id>%(SOURCE)s</source_id>
    <prop_id>%(PROP)s</prop_id>
</NamedObject>
</elem>

'''

snk = '''<elem>
<Name>SINK_SMEAR</Name>
<Frequency>1</Frequency>
<Param>
    <version>5</version>
    <Sink>
        <version>2</version>
        <SinkType>SHELL_SINK</SinkType>
        <j_decay>3</j_decay>
        <SmearingParam>
            <wvf_kind>GAUGE_INV_GAUSSIAN</wvf_kind>
            <wvf_param>2.0</wvf_param>
            <wvfIntPar>5</wvfIntPar>
            <no_smear_dir>3</no_smear_dir>
        </SmearingParam>
    </Sink>
</Param>
<NamedObject>
    <gauge_id>default_gauge_field</gauge_id>
    <prop_id>%(POINT_SINK_PROP)s</prop_id>
    <smeared_prop_id>%(SMEARED_SINK_PROP)s</smeared_prop_id>
</NamedObject>
</elem>

'''

lalibe_2pt = '''<elem>
<Name>BARYON_CONTRACTIONS</Name>
    <Frequency>1</Frequency>
    <BaryonParams>
        <ng_parity>true</ng_parity>
        <h5_file_name>./lalibe_2pt_spectrum.h5</h5_file_name>
        <path>/%(H5_PATH)s</path>
        <mom_list>
            <elem>0 0 0</elem>
            <elem>0 0 1</elem>
            <elem>0 0 2</elem>
        </mom_list>
        <particle_list>
            <elem>proton</elem>
        </particle_list>
    </BaryonParams>
<NamedObject>
    <up_quark>%(PROP)s</up_quark>
    <down_quark>%(PROP)s</down_quark>
</NamedObject>
</elem>


'''

h5_write = '''<elem>
<Name>HDF5_WRITE_NAMED_OBJECT</Name>
<Frequency>1</Frequency>
<NamedObject>
    <object_id>%(OBJECT_ID)s</object_id>
    <object_type>%(OBJECT_TYPE)s</object_type>
</NamedObject>
<File>
    <file_name>%(H5_FILE)s</file_name>
    <path>/%(H5_PATH)s</path>
    <obj_name>%(H5_NAME)s</obj_name>
</File>
</elem>

'''

seqsource = '''<elem>
<Name>LALIBE_SEQSOURCE</Name>
<Frequency>1</Frequency>
<SeqSourceParams>
    <particle>%(PARTICLE)s</particle>
    <flavor>%(FLAV)s</flavor>
    <source_spin>%(SRC_SPIN)s</source_spin>
    <sink_spin>%(SNK_SPIN)s</sink_spin>
    <sink_mom>0 0 0</sink_mom>
    %(T_SINK)s
</SeqSourceParams>
<NamedObject>
    <gauge_id>default_gauge_field</gauge_id>
    <up_quark>%(PROP)s</up_quark>
    <down_quark>%(PROP)s</down_quark>
    <seqsource_id>%(SEQSOURCE)s</seqsource_id>
</NamedObject>
</elem>

'''

multi_prop_add='''<elem>
<Name>MULTI_PROP_ADD</Name>
    <Frequency>1</Frequency>
    <PropWeights>
        <delete_props>false</delete_props>
        <weights>%(PROP_WEIGHTS)s</weights>
    </PropWeights>
    <NamedObject>
        <prop_ids>
%(PROPS)s
        </prop_ids>
        <result_prop>%(RESULT_PROP)s</result_prop>
    </NamedObject>
</elem>

'''

coherent_seqsource='''<elem>
<Name>COHERENT_SEQSOURCE</Name>
    <Frequency>1</Frequency>
    <SinkParams>
        <t_sep>%(T_SEP)s</t_sep>
        <j_decay>%(J_DECAY)s</j_decay>
    </SinkParams>
    <NamedObject>
        <sink_ids>
%(SINKS)s
        </sink_ids>
        <result_sink>%(RESULT_SINK)s</result_sink>
    </NamedObject>
</elem>

'''

bar3pt = '''<elem>
<Name>LALIBE_BAR3PTFN</Name>
    <Frequency>1</Frequency>
    <Param>
        <version>7</version>
        <j_decay>3</j_decay>
        <currents>
            <elem>V4</elem>
            <elem>A3</elem>
            <elem>S</elem>
            <elem>P</elem>
            <elem>T34</elem>
        </currents>
        <mom_list>
            <elem>0 0 0</elem>
            <elem>0 0 1</elem>
            <elem>0 0 2</elem>
        </mom_list>
        <h5_file_name>./lalibe_3ptfn.h5</h5_file_name>
        <path>/%(H5_PATH)s</path>
    </Param>
    <NamedObject>
        <gauge_id>default_gauge_field</gauge_id>
        <prop_id>%(PROP)s</prop_id>
        <seqprops>
            <elem>
                <seqprop_id>%(SEQPROP)s</seqprop_id>
                <gamma_insertion>0</gamma_insertion>
            </elem>
        </seqprops>
    </NamedObject>
</elem>

'''

bar3pt_4d = '''<elem>
<Name>LALIBE_BAR3PTFN</Name>
    <Frequency>1</Frequency>
    <Param>
        <version>7</version>
        <j_decay>3</j_decay>
        <currents>
            <elem>V4</elem>
            <elem>A3</elem>
            <elem>S</elem>
            <elem>P</elem>
            <elem>T34</elem>
        </currents>
        <h5_file_name>./lalibe_3ptfn_4d.h5</h5_file_name>
        <path>/%(H5_PATH)s</path>
    </Param>
    <NamedObject>
        <gauge_id>default_gauge_field</gauge_id>
        <prop_id>%(PROP)s</prop_id>
        <seqprops>
            <elem>
                <seqprop_id>%(SEQPROP)s</seqprop_id>
                <gamma_insertion>0</gamma_insertion>
            </elem>
        </seqprops>
    </NamedObject>
</elem>

'''

end = '''</InlineMeasurements>
<nrow>4 4 4 8</nrow>
</Param>

<RNG>
  <Seed>
    <elem>11</elem>
    <elem>11</elem>
    <elem>11</elem>
    <elem>0</elem>
  </Seed>
</RNG>

<Cfg>
  <cfg_type>WEAK_FIELD</cfg_type>
  <cfg_file>dummy</cfg_file>
</Cfg>
</lalibe>
'''
