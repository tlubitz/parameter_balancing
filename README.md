# Parameter Balancing

Parameter balancing is a tool for metabolic modelling in systems biology. It is implemented in Python3 and its code underlies the PEP8 guidelines. There are 4 major ways of employing parameter balancing for your project:

1. Online version

The tool can be employed via www.parameterbalancing.net. All required knowledge can be found on the webpage.

2. Standalone server version

The online tool, which is hosted on www.parameterbalancing.net, is open source and can also be hosted as an offline server. Thus, you can tailor it to your specific needs and use it via your browser of choice. You will to download a version of the web framework [web2py](http://www.web2py.com/). Then, you can directly download the parameter balancing application into the applications folder of your web2py server (web2py/applications/). You will find the parameter balancing application in this repository (parameter_balancing/web_version). For manipulating the code and user interface, you will require basic knowledge about the web2py web framework.

3. Python PyPi installer

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
file, an SBtab prior distribution file, and an SBtab configuration file. You will find information on all these
file types in this repository's examples section.

4. Standalone commandline version

To run parameter balancing as a commandline tool, the package needs to be installed as explained in (2). Then,
it can be executed in the commandline as follows:

> python3 -m pbalancing.parameter_balancing model.xml

where model.xml corresponds to the path to your SBML model. It is also possible to add further input files, such as
an SBtab parameter files (.tsv), an SBtab prior information file (.tsv), and an SBtab options file (.tsv) for the
configuration of parameter balancing. Information on the SBtab format can be found on www.sbtab.net, more information
on the mentioned file types can be found in the parameter balancing manual on https://rumo.biologie.hu-berlin.de/pb/static/documentation.html,
and example files can be found in the parameter balancing download section on https://rumo.biologie.hu-berlin.de/pb/static/download.html.
The optional files can be provided as further arguments in the commandline, such as

> python3 -m pbalancing.parameter_balancing model.xml parameter_file.tsv prior.tsv options.tsv

where the first argument MUST be the SBML file and the order and amount of SBtab files is optional, as long as they
are correct in their file paths and SBtab file syntax.

If you use parameter balancing, please cite http://pubs.acs.org/doi/abs/10.1021/jp108764b for details.

For questions and feedback,please consult the repository admin.