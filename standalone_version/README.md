# Parameter Balancing: Standalone version

Parameter balancing is a tool for metabolic modelling in systems biology. It is implemented in Python3 and its code underlies the PEP8 guidelines. This subdirectory of the project holds the files for the standalone Python3 version. It can be downloaded and used from commandline as is, the call being

> python3 parameter_balancing.py model.xml

where model.xml corresponds to the path of your SBML model. It is also possible to provide further input files, such as
an SBtab parameter files (.tsv), an SBtab prior information file (.tsv), and an SBtab options file (.tsv) for the
configuration of parameter balancing. Providing complete file information would look like this:

> python3 parameter_balancing model.xml --sbtab_data data_file.tsv --sbtab_prior prior_file.tsv --sbtab_options options_file.tsv

You can create a log file by setting the flag -l, you can use pseudo values to account for a lack of data by setting the flag -p, you can
watch program outputs on your commandline by setting the flag -v. Information on the SBtab format can be found on www.sbtab.net, more information
on the mentioned file types can be found in the parameter balancing manual in parameter_balancing/standalone_version/files/manual.pdf,
and example files can be found in parameter_balancing/standalone_version/files/example_files/.

You will be required to install several Python packages, though. A list of these packages can be found in parameter_balancing/requirements.txt.

If you use parameter balancing, please cite http://pubs.acs.org/doi/abs/10.1021/jp108764b for details.

For questions and feedback,please consult the repository admin.