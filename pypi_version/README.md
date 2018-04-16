# Parameter Balancing: PyPi version

Parameter balancing is a tool for metabolic modelling in systems biology. It is implemented in Python3 and its code underlies the PEP8 guidelines. This subdirectory holds the code for the PyPi python installer. This code is not meant to be downloaded but to be installed via pip. To do this, first of all you need Python3. Next, you will need the pip3 installer. You can find information
on how to do this here: https://pip.pypa.io/en/stable/installing/. Then you can type on your commandline

> sudo pip3 install pbalancing

This will also install libsbml and tablib on your computer if these libraries are missing. You can now employ
Parameter Balancing as a Python3 package by, e.g., writing a script such as

> from pbalancing import parameter_balancing
> 
> parameter_balancing.parameter_balancing_wrapper('model.xml')

In this example case, 'model.xml' is the file name of an SBML model. Further optional arguments are an SBtab parameter
file, an SBtab prior distribution file, and an SBtab configuration file. You will find examples for all these
file types in parameter_balancing/standalone_version/files/.

After installation, you can also employ parameter balancing from your commandline by typing

> python3 -m pbalancing.parameter_balancing model.xml

where model.xml corresponds to the path of your SBML model. It is also possible to provide further input files, such as
an SBtab parameter files (.tsv), an SBtab prior information file (.tsv), and an SBtab options file (.tsv) for the
configuration of parameter balancing. Providing complete file information would look like this:

> python3 -m pbalancing.parameter_balancing model.xml --sbtab_data data_file.tsv --sbtab_prior prior_file.tsv --sbtab_options options_file.tsv

You can create a log file by setting the flag -l, you can use pseudo values to account for a lack of data by setting the flag -p, you can
watch program outputs on your commandline by setting the flag -v.

> usage: parameter_balancing.py [-h] [--sbtab_data SBTAB_DATA]
>                               [--sbtab_prior SBTAB_PRIOR]
>                               [--sbtab_options SBTAB_OPTIONS]
>                               [--output_name OUTPUT_NAME] [-l] [-p] [-v]
>                               sbml
> 
> positional arguments:
>   sbml                  Path to an SBML file.
> 
> optional arguments:
>   -h, --help            show this help message and exit
>   --sbtab_data SBTAB_DATA
>                         Path to an SBtab data file.
>   --sbtab_prior SBTAB_PRIOR
>                         Path to an SBtab prior file.
>   --sbtab_options SBTAB_OPTIONS
>                         Path to an SBtab options file.
>   --output_name OUTPUT_NAME
>                         Choose a name for the output files.
>   -l, --pb_log          Flag to print a log file.
>   -p, --pb_pseudos      Flag for usage of pseudo values.
>   -v, --verbose         Flag to display script messages.

Information on the SBtab format can be found on www.sbtab.net, more information
on the mentioned file types can be found in the parameter balancing manual in parameter_balancing/standalone_version/files/manual.pdf,
and example files can be found in parameter_balancing/standalone_version/files/example_files/.

If you use parameter balancing, please cite http://pubs.acs.org/doi/abs/10.1021/jp108764b for details.

For questions and feedback,please consult the repository admin.
