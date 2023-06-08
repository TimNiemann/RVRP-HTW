## Requirements

To run, this codes needs an installation of SCIP version 8.0.1 under the variable $(SCIP_PATH).
See the website for an installation guide:
https://scipopt.org/#scipoptsuite

## Generating Testdata

To generate the test instances, run the following script:

```markdown
$ ./data_generation_script/test-modeldata-generation.sh
```

## Usage
 
$ make
$ ./bin/columnGeneration <path-to-modeldata-file.dat> 
    [-w <appointment window length in sec (mandatory)>] 
    [-c <objective function weigth in [0,1] (mandatory)>]
    [-o <output json file (optional; default used if no name is provided)>]
    [-s <input solution json file (mandatory)>]
    [-a <objective function parameters (delay, travel time, encounter probability) (default: 1, 0, 0) (mandatory) >]
    [-g <value for Gamma (mandatory)>]

Example usage with modeldata "bayern0_r20_d5_w0.2.dat" in "../test-data/model-data-paper/paper-evaluation-small-tw-0/bayern" 
and an appointment window length of 60m (=3600s) with the creation of an output file:

```markdown
$ make
$ /bin/columnGeneration ../test-data/model-data-paper/paper-evaluation-small-tw-0/bayern/bayern0_r20_d5_w0.2.dat -w 3600 -o
```

If no values for the objective parameters are given, the default values are (1, 0, 0) is used.
If no appointment window length is given, the default value of 0 is used.
If no outputfile name is given, "outputfile.json" is used.
Configuration for debugging and parameters for heuristics can be set in "src/tools_vrp.h".


## License

Shield: [![CC BY 4.0][cc-by-shield]][cc-by]

This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg
