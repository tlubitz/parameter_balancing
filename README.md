# Parameter Balancing

Parameter balancing is a tool for metabolic modelling in systems biology. It is implemented in Python3 and its code underlies the PEP8 guidelines. There are 4 major ways of employing parameter balancing for your project:

1. Online version

The tool can be employed via www.parameterbalancing.net. All required knowledge can be found on the webpage.

2. Python PyPi installer

To install Parameter Balancing as a Python3 package, first of all you need Python3. Next, you will need the pip3 installer. You can find information
on how to do this here: https://pip.pypa.io/en/stable/installing/
Afterwards, install Parameter Balancing by typing in your command line:

> sudo pip3 install pbalancing

This will also install libsbml and tablib on your computer if these libraries are missing. You can now employ
Parameter Balancing as a Python3 package by, e.g., writing a script such as

> from pbalancing import parameter_balancing
> 
> parameter_balancing.parameter_balancing_wrapper('model.xml')

In this example case, 'model.xml' is the file name of an SBML model. Further optional arguments are an SBtab parameter
file, an SBtab prior distribution file, and an SBtab configuration file. You will find examples for all these
file types in parameter_balancing/standalone_version/files/.

3. Standalone commandline version

To run parameter balancing as a commandline tool, the package needs to be installed as explained in (2). Then,
it can be executed in the commandline as follows:

> python3 -m pbalancing.parameter_balancing model.xml

where model.xml corresponds to the path of your SBML model. It is also possible to provide further input files, such as
an SBtab parameter files (.tsv), an SBtab prior information file (.tsv), and an SBtab options file (.tsv) for the
configuration of parameter balancing. Providing complete file information would look like this:

> python3 -m pbalancing.parameter_balancing model.xml --sbtab_data data_file.tsv --sbtab_prior prior_file.tsv --sbtab_options options_file.tsv

You can create a log file by setting the flag -l, you can use pseudo values to account for a lack of data by setting the flag -p, you can
watch program outputs on your commandline by setting the flag -v. Information on the SBtab format can be found on www.sbtab.net, more information
on the mentioned file types can be found in the parameter balancing manual in parameter_balancing/standalone_version/files/manual.pdf,
and example files can be found in parameter_balancing/standalone_version/files/example_files/.

If you do not want to install the pip package, you can still use the commandline modules in the subdirectory parameter_balancing/standalone_version.
The usage works as a standard call of a Python module:

> python3 parameter_balancing.py model.xml

Here, as well, you can use the optional file provision like explained above. You will be required to install several Python packages, though. A list of
these packages can be found in parameter_balancing/requirements.txt.

4. Standalone server version

The online tool, which is hosted on www.parameterbalancing.net, is open source and can also be hosted as an offline server. Thus, you can tailor it to your specific needs and use it via your browser of choice. You will need to download a version of the web framework [web2py](http://www.web2py.com/). Then, you can directly download the parameter balancing application into the applications folder of your web2py server (web2py/applications/). You will find the parameter balancing application in this repository (parameter_balancing/web_version). For manipulating the code and user interface, you will require basic knowledge about the web2py web framework.

If you use parameter balancing, please cite http://pubs.acs.org/doi/abs/10.1021/jp108764b for details.

For questions and feedback,please consult the repository admin.