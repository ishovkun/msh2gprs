#!/usr/bin/env python

from Postprocessor import Postprocessor
import sys                      # argument handling
import os

def main(argv):
    usage_msg = """
    Postprocessor takes either zero or a single argument.
    If zero arguments is passed, the script will use the invocation directory
    to run the program.
    If a single argument is passed, it will be used as the working directory.
    The argument must be a valid directory and contain a
    preprocessor_config.yaml file made by msh2gprs.

    Examples:
    - While being in directory ~/simulation_dir/case1, invoke
        $ python ~/preprocessor_path/main.py
    to use the current directory as a working directory.

    - While being in any directory, invoke
    $ python ~/preprocessor_path/main.py /some/work/path
    to use /some/work/path as a working directory.
    """
    # print(os.getcwd())
    # print(sys.argv[0])
    # print(sys.argv[1])
    case_path = ""
    if len(argv) == 1: # use current path
        case_path = os.getcwd() + "/"
    elif len(argv) == 2:        # use argument
        case_path = argv[1]
        if case_path[-1] != "/":
            case_path += "/"
        print(case_path)
    else:
        print ("Wrong number of arguments.\n")
        print (usage_msg)
        exit(1)

    print("Running in '%s'" % case_path)
    postprocessor = Postprocessor(case_path)
    postprocessor.run()
    # except Exception as e:
    #     print(e)

main(sys.argv)
