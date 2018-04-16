### Python modules for parameter balancing ###

This is a standalone library to use parameter balancing in kinetic models.

To perform parameter balancing in the command line, call parameter_balancing.py by typing

> python parameter_balancing.py argument1 argument2 argument3 argument4

where argument1 = path to SBML model (obligatory). The other three arguments are optional: argument2 = SBtab data file;
argument3 = path to a file with prior distributions; argument4 = path to a file with configuration options.

Examples for these file types are found in the static directory.

If the provided SBML model is valid, the script is started automatically. After successful balancing of the model it
will ask the user for a name for the output files. This is optional; if no name is provided, the file will be named
after the input SBML file.

The SBML model can be validated beforehand on http://sbml.org/Facilities/Validator/
The SBtab files can be validated using the validation tool on www.sbtab.net.

More information on parameter balancing and the input files can be found on www.parameterbalancing.net.
