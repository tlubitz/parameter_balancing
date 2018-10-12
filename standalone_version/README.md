# Parameter Balancing: Standalone version

Parameter balancing is a tool for metabolic modelling in systems biology. It is implemented in Python3 and its code underlies the PEP8 guidelines. This subdirectory of the project holds the files for the standalone Python3 version. It can either be downloaded and used from commandline as is or embedded in your own Python3 projects. In both cases, you will be required to install Python packages to make the scripts work. A list of these packages can be found in parameter_balancing/requirements.txt.

<h3>Parameter Balancing on the commandline</h3>

Parameter Balancing can be employed on the commandline by

> python3 cl_balancing.py model.xml

where model.xml corresponds to the path of your SBML model. It is also possible to provide further input files, such as
an SBtab parameter files (.tsv), an SBtab prior information file (.tsv), and an SBtab options file (.tsv) for the
configuration of parameter balancing. Providing complete file information would look like this:

> python3 cl_balancing.py model.xml --sbtab_data data_file.tsv --sbtab_prior prior_file.tsv --sbtab_options options_file.tsv

You can create a log file by setting the flag -l, you can use pseudo values to account for a lack of data by setting the flag -p, you can watch program outputs on your commandline by setting the flag -v. Information on the SBtab format can be found on www.sbtab.net, more information on the mentioned file types can be found in the parameter balancing manual in this repository's parameter_balancing/standalone_version/files/manual.pdf, and example files can be found in parameter_balancing/standalone_version/files/example_files/.

<h3>Embedding Parameter Balancing in your Python3 package</h3>

You can embed the modules of Parameter Balancing in your own Python3 workflow.

```python
  import parameter_balancing
  
  (balanced_sbml, balanced_sbtab) = parameter_balancing.parameter_balancing_wrapper(sbml, sbtab_data, sbtab_prior, sbtab_options, verbose, no_pseudo_values, output_name, pb_log)
```
<strong>Input arguments</strong>
<ul>
  <li>sbml (path to SBML file, REQUIRED)</li>
  <li>sbtab_data (path to SBtab data file, OPTIONAL)</li>
  <li>sbtab_prior (path to SBtab prior file, OPTIONAL)</li>
  <li>sbtab_options (path to SBtab options file, OPTIONAL)</li>
  <li>verbose (Boolean, enable messages on commandline, OPTIONAL)</li>
  <li>no_pseudo_values (Boolean, disable usage of pseudo values, OPTIONAL)</li>
  <li>output_name (name for the output files, OPTIONAL)</li>
  <li>pb_log (Boolean, enable writing of a log file, OPTIONAL)</li>
</ul>

<strong>Output parameters</strong>
<ul>
  <li>balanced_sbml (SBML file object with balanced content)</li>
  <li>balanced_sbtab (SBtab object with balanced content)</li>
</ul>


<h3>Citation and Contact</h3>

If you use parameter balancing, please cite http://pubs.acs.org/doi/abs/10.1021/jp108764b for details.

If you are encountering trouble with any of the above, please file a bug report in Github. You can also feel free to file feature requests in the same manner.
