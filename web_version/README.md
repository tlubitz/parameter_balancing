# Parameter Balancing

Parameter balancing is a tool for estimating parameters in kinetic modelling. It can be employed

__ 1. On our website www.parameterbalancing.net __  
Parameter balancing can be used on our web interface by uploading the necessary files.

__ 2. Installed as a python3 package __  
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
file types in this repository's examples section: /pb/static/files/

__ 3. Installed as a standalone web2py application __  
This bitbucket repository fully corresponds to a web2py application directory. It needs to be embedded in 
a web2py server. Web2py can be downloaded on http://www.web2py.com/init/default/download and this repository 
needs to be embedded in /web2py/applications/pb to make it function properly. 

__ 4. As a commandline tool __
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