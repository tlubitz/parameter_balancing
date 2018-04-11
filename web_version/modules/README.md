### Python modules for parameter balancing ###

This directory contains the python modules required for parameter balancing.

To perform parameter balancing in the command line, call balance_reaction_parameters.py by typing

> python parameter_balancing.py argument

where argument = path to SBML model (obligatory).

If the provided SBML model is valid, the script is starting and you can follow the instructions on the command line.
You will be asked to provide the following, optional files:

- SBtab parameter file
- SBtab file with information on prior distributions
- SBtab file with individual configurations for the parameter balancing.

There are examples for these types of SBtab files in the repository path /pb/static/files.

The SBML model can be validated beforehand on http://sbml.org/Facilities/Validator/
The SBtab files can be validated using the validation tool on www.sbtab.net.

More information on parameter balancing and the input files can be found on www.parameterbalancing.net.