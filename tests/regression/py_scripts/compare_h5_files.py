import os, sys
import argparse
import h5py
from typing import Union
from typing import Optional
import numpy as np

def get_dsets(container: Union[h5py.File, h5py.Group], parent_name: Optional[str] = None, verbose=True):
    if isinstance(container, h5py.File) and verbose:
        print('Locating all dsets of h5 file %s' %container.filename)

    dsets = {}
    for key in container:
        obj = container[key]
        address = os.path.join(parent_name, key) if parent_name else key

        if isinstance(obj, h5py.Dataset):
            dsets[address] = obj
        elif isinstance(obj, h5py.Group):
            dsets.update(get_dsets(obj, parent_name=address))

    return dsets

def assert_h5files_equal(known_result, new_result, atol=0., rtol=1.e-10, verbose=False):
    with h5py.File(known_result,'r') as h5_known:
        dsets_known = get_dsets(h5_known, verbose=verbose)

        with h5py.File(new_result,'r') as h5_new:
            dsets_new = get_dsets(h5_new, verbose=verbose)

            known_keys = set(dsets_known.keys())
            new_keys   = set(dsets_new.keys())

            if new_keys != known_keys:
                data_equal = False
                if verbose:
                    raise AssertionError(
                        (
                            "Files have different datasets:"
                            "\n---Dsets in known but not in new---\n\t%s"
                            "\n---Dsets in new but not in known---\n\t%s"
                        )
                        %(
                            "\n\t".join(known_keys.difference(new_keys)),
                            "\n\t".join(new_keys.difference(known_keys)),
                        )
                    )
            else:
                try:
                    for key in known_keys:
                        if verbose:
                            error_message = "Dataset %s has unequal values" %key
                        else:
                            error_message = None
                        np.testing.assert_allclose(
                            dsets_known[key],
                            dsets_new[key],
                            atol=atol,
                            rtol=rtol,
                            err_msg=error_message,
                        )
                    data_equal = True
                except Exception as e:
                    data_equal = False
                    if verbose:
                        print(e)
    return data_equal

def main():
    PARSER = argparse.ArgumentParser()
    PARSER.add_argument('known_file',    type=str, help='known_result_file')
    PARSER.add_argument('new_file',      type=str, help='new_result_file')
    PARSER.add_argument('-a','--atol',   type=float, default=0., help='absolute tolerance for comparison')
    PARSER.add_argument('-r','--rtol',   type=float, default=1.e-10, help='relative tolerance for comparison')
    PARSER.add_argument('-v','--verbose',default=False,action='store_const',const=True,help='verbose? [%(default)s]')
    args = PARSER.parse_args()

    h5_files_equal = assert_h5files_equal(args.known_file, args.new_file, atol=args.atol, rtol=args.rtol, verbose=True)
    if h5_files_equal:
        print("[+] Files are equal!")
    else:
        print("[-] Files are not equal :(")

if __name__ == "__main__":
    main()
