# Readme

## Building

 - Install docker on your machine
 - Run `docker build -t ssp .` in the directory where the code is extracted.
   If you have Gurobi 9.1 and Java version >= 17 installed, you can also run the code directly on your machine by replicating the steps from the Dockerfile.
 - Obtain a WLS licence from Gurobi (see https://license.gurobi.com/, free for academic use) and save the `gurobi.lic` somewhere

## Running

 - To run a `command`, execute `docker run -it -v <path to gurobi.lic>:/opt/gurobi/gurobi.lic bin/ssp-cvar <command>`
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
 - If the execution stops because of OOM errors, ensure that docker has enough memory available and potentially add `-e JAVA_OPTS="-Xmx16G"` directly after `docker run` to increase the memory available to the JVM.

## Evaluation

Run the following commands for the respective models:

 - Grid5: `bin/ssp-cvar models/gridworld/gridworld.5.p.tra models/gridworld/gridworld.5.p.lab station reach lp 0.10 output/gridworld.5.p.json`
 - Grid20: `bin/ssp-cvar models/gridworld/gridworld.20.d.tra models/gridworld/gridworld.20.d.lab station reach vi 0.10 output/gridworld.20.d.json`
 - FireWire: `bin/ssp-cvar models/firewire/firewire.tra models/firewire/firewire.lab absorbing reach vi 0.10 output/firewire.json`
 - WLAN: `bin/ssp-cvar models/wlan/wlan2.tra models/wlan/wlan2.lab absorbing reach vi,lp 0.10 output/wlan.json`
 - Walk: `bin/ssp-cvar models/random_walk/random_walk.tra models/random_walk/random_walk.lab goal reach vi,lp range:0.01,0.1,19 output/walk-multi.json`

We have excluded the variants which time out, to include them simply change `vi` / `lp` to `vi,lp` in the respective commands.
To view the results, simply open the created JSON files.

Run `python3 eval.py figure walk-multi.json` to output a CSV for the different threshold performances.

By increasing the logging level of Java's default logging facility, further details of the computation can be viewed.

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