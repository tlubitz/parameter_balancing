#!/usr/bin/env python
import argparse

try: from . import parameter_balancing_core
except: import parameter_balancing_core

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('sbml', help='Path to an SBML file.')
    parser.add_argument('--sbtab_data', help='Path to an SBtab data file.')
    parser.add_argument('--sbtab_prior', help='Path to an SBtab prior file.')
    parser.add_argument('--sbtab_options', help='Path to an SBtab options file.')
    parser.add_argument('--output_name', help='Choose a name for the output files.')
    parser.add_argument('-l', '--pb_log', help='Flag to print a log file.', action='store_true')
    parser.add_argument('-p', '--no_pseudo_values', help='Flag for disabling the usage of pseudo values.', action='store_true')
    parser.add_argument('-v', '--verbose', help='Flag to display script messages.', action='store_true')

    args = parser.parse_args()

    parameter_balancing_core.parameter_balancing_wrapper(args.sbml,
                                                         args.sbtab_data,
                                                         args.sbtab_prior,
                                                         args.sbtab_options,
                                                         args.verbose,
                                                         args.no_pseudo_values,
                                                         args.output_name,
                                                         args.pb_log)

