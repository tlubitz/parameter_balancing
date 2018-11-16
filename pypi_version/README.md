# Parameter Balancing: PyPi version

Parameter Balancing is a tool for metabolic modelling in systems biology. It is implemented in Python3 and its code underlies the PEP8 guidelines. This subdirectory holds the code for the PyPi python installer. This code is not meant to be downloaded but to be installed via pip. To do this, first of all you need Python3. Next, you will need the pip3 installer. You can find information
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

After installation, you can also employ Parameter Balancing from your commandline by typing

> python3 -m pbalancing.parameter_balancing model.xml

where model.xml corresponds to the path of your SBML model. It is also possible to provide further input files, such as
an SBtab parameter files (.tsv), an SBtab prior information file (.tsv), and an SBtab options file (.tsv) for the
configuration of Parameter Balancing. Providing complete file information would look like this:

> python3 -m pbalancing.parameter_balancing model.xml --sbtab_data data_file.tsv --sbtab_prior prior_file.tsv --sbtab_options options_file.tsv

You can create a log file by setting the flag -l, you can use pseudo values to account for a lack of data by setting the flag -p, you can
watch program outputs on your commandline by setting the flag -v.

> usage: parameter_balancing.py [-h] [--sbtab_data SBTAB_DATA] <br>
>                               [--sbtab_prior SBTAB_PRIOR] <br>
>                               [--sbtab_options SBTAB_OPTIONS] <br>
>                               [--output_name OUTPUT_NAME] [-l] [-p] [-v] <br>
>                               sbml <br>
>  <br>
> positional arguments: <br>
>   sbml                  Path to an SBML file. <br>
>  <br>
> optional arguments: <br>
>   -h, --help            show this help message and exit <br>
>   --sbtab_data SBTAB_DATA <br>
>                         Path to an SBtab data file. <br>
>   --sbtab_prior SBTAB_PRIOR <br>
>                         Path to an SBtab prior file. <br>
>   --sbtab_options SBTAB_OPTIONS <br>
>                         Path to an SBtab options file. <br>
>   --output_name OUTPUT_NAME <br>
>                         Choose a name for the output files. <br>
>   -l, --pb_log          Flag to print a log file. <br>
>   -p, --pb_pseudos      Flag for usage of pseudo values. <br>
>   -v, --verbose         Flag to display script messages. <br>

Information on the SBtab format can be found on www.sbtab.net, more information
on the mentioned file types can be found in the Parameter Balancing manual in parameter_balancing/standalone_version/files/manual.pdf,
and example files can be found in parameter_balancing/standalone_version/files/example_files/.

If you use Parameter Balancing, please cite http://pubs.acs.org/doi/abs/10.1021/jp108764b for details.

For questions and feedback,please consult the repository admin.
