# Readme

## Building

 - Install docker on your machine
 - Run `docker build -t ssp .` in the directory where the code is extracted.
   If you have Gurobi 9.5.1 and Java version >= 17 installed, you can also run the code directly on your machine by replicating the steps from the Dockerfile.
 - Obtain a WLS licence from Gurobi (see https://license.gurobi.com/, free for academic use) and save the `gurobi.lic` somewhere

## Running

 - To run a `command`, execute `docker run -it -v <path to gurobi.lic>:/opt/gurobi/gurobi.lic ssp <command>`
 - Optional: When running experiments, also add `-v <some path>:/app/output` to the docker call to obtain the output files on your machine.
   So, the complete command is
   ```
   docker run -it -v <path to gurobi.lic>:/opt/gurobi/gurobi.lic -v <some path>:/app/output ssp <command>
   ```
 - The models can be found in the folder `/app/models`
 - Running our implementation always has the form
   ```
   bin/ssp-risk <model transitions> <model labels> <goal label> reach <methods: vi,lp> <threshold> output/<output name>.json
   ```
 - All relevant data can be found in `<some path>/<name>.json`
 - To specify memory limits, ensure that docker has enough memory available and add `-e JAVA_OPTS="-Xmx10G"` directly after `docker run`.

## Evaluation

Run the following commands for the respective models:

 - Grid4: `bin/ssp-cvar -m models/gridworld/gridworld.04 --goal station --method LP,VI -t 0.10 -o output/gridworld.04.json`
 - Grid8: `bin/ssp-cvar -m models/gridworld/gridworld.08.d --goal station --method LP,VI -t 0.10 -o output/gridworld.08.json`
 - Grid16: `bin/ssp-cvar -m models/gridworld/gridworld.16.d --goal station --method LP,VI -t 0.10 -o output/gridworld.16.json`
 - Grid32: `bin/ssp-cvar -m models/gridworld/gridworld.32.d --goal station --method VI -t 0.10 -o output/gridworld.32.json`
 - FireWire: `bin/ssp-cvar -m models/firewire/firewire --goal absorbing --method VI -t 0.10 -o output/firewire.json`
 - WLAN: `bin/ssp-cvar -m models/wlan/wlan2 --goal absorbing --method LP,VI -t 0.10 -o output/wlan.json`
 - Walk: `bin/ssp-cvar -m models/random_walk/random_walk --goal goal --method VI -t "range:0.001,0.1,100|logrange:-7,-3,100" -o output/walk-multi.json`

Variants which time out or run out of memory are excluded, to include them simply change `VI` / `LP` to `VI,LP` in the respective commands.
To view the results, simply open the created JSON files.

Run `python3 eval.py walk-multi.json` to output a CSV for the different threshold performances.

By increasing the logging level of Java's default logging facility, further details of the computation can be viewed (note that the statistics computation may increase runtime by a few percent).

## Modifying

 - The `sources.tar.gz` contain all sources of our implementation (both Java and a python prototype).
   In fact, the binary executed in the evaluation is compiled from the contained sources.
   If you want to modify the sources, simply re-package the `sources.tar.gz` and re-build the docker image.
 - The models are written in the PRISM modelling language.
   They can be modified and then exported using PRISM (a release can be obtained at https://www.prismmodelchecker.org/), using the call
   ```
   prism <model>.nm -exportmodel .all -const <constants>
   ```
   This produces (among others) the required `.tra` and `.lab` files.
   The constants we used to create the models can be found in the `.const` files in each model folder.