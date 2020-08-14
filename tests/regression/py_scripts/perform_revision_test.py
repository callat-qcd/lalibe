import argparse
import compare_h5_files as diff_h5

''' When you add a new routine, which generates new data, you have to verify
    that data with some other means, then generate the appropriate h5 data file
    to place in the "known_results" folder, and generate the xml that will
    create this file, such that it will be compared in future revision tests
'''

''' NOTE - prop files fail the test due to the .PropagatorDouble3 data type
           the final spectrum files used newly generated propagator files so
           this indirectly tests those as well.
    #'test_propagator.h5',
    #'test_fh_propagator.h5',
    #'test_coherent_sink.h5',
    #'test_seqprop.h5',

'''
list_of_h5_files = [
    'lalibe_2pt_spectrum.h5',
    'lalibe_fh_proton.h5',
    'lalibe_3ptfn.h5',
    'lalibe_3ptfn_coherent_sink.h5',
    'lalibe_3ptfn_2src_coherent_sink.h5'
]

diff_1_2 = dict()
diff_1_2['lalibe_2pt_spectrum.h5'] = 'lalibe_2pt_spectrum_lime.h5'
diff_1_2['lalibe_3ptfn.h5'] = 'lalibe_3ptfn_coherent_sink.h5'

PARSER = argparse.ArgumentParser(description='Perform revision test after files are generated')
PARSER.add_argument('known_file_path', type=str, help='known_result_folder')
PARSER.add_argument('new_file_path',   type=str, help='new_result_folder')
PARSER.add_argument('-a','--atol',   type=float, default=0., help='absolute tolerance for comparison')
PARSER.add_argument('-r','--rtol',   type=float, default=5.e-9, help='relative tolerance for comparison')
PARSER.add_argument('-v','--verbose',default=False,action='store_const',const=True,help='verbose? [%(default)s]')
args = PARSER.parse_args()

''' perform test '''
revision_tests = dict()
for f5 in list_of_h5_files:
    f_known = args.known_file_path+'/'+f5
    f_new   = args.new_file_path+'/'+f5

    revision_tests[f5] = diff_h5.assert_h5files_equal(
        f_known,
        f_new,
        atol=args.atol,
        rtol=args.rtol,
        verbose=args.verbose
        )
    if revision_tests[f5]:
        print('PASS:    ',f5)
    else:
        print('FAIL:    ',f5)

for f5 in diff_1_2:
    f_known = args.known_file_path+'/'+f5
    f_new   = args.new_file_path+'/'+diff_1_2[f5]

    revision_tests[f5] = diff_h5.assert_h5files_equal(
        f_known,
        f_new,
        atol=args.atol,
        rtol=args.rtol,
        verbose=args.verbose
        )
    if revision_tests[f5]:
        print('PASS:    ',f5, ' = ',diff_1_2[f5])
    else:
        print('FAIL:    ',f5, ' != ',diff_1_2[f5])
