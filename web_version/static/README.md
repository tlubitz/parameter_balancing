### Static files of the Parameter Balancing Repository ###

The static files comprise static HTML files for the website (*.html), images for the website visualisation (/static/images/*), and a variety of default and exemplary data files (/static/files/*).

The exemplary data files include SBtab kinetic parameter files, SBML model files (*.xml), and the following SBtab default files, which are crucial to the parameter balancing workflow:

- pb_options.tsv (customary configuration options of the parameter balancing; numerous settings of the tool can be changed via this file)
- definitions.tsv (definitions of the SBtab table types and allowed columns; this is required for the validation of SBtab files and can be extended according to individual needs. Please see www.sbtab.org for more details)
- pb_prior.tsv (file that holds information on the prior distributions of the different kinetic parameter types; these are required for the Bayesian estimation within the parameter balancing)