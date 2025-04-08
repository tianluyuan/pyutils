Originally by robj
Last edited by tianlu, 15/12/05


Hi! This is a README for how to use these scripts to generate:
  GENIE vector files (numc / gnmc)
  Cherry-Picked bundled vector files (nucp)
  ND280 production files (nd280)

Executables are in `bin/` or `flux/bin` within this production package.

Here's a [page](http://www.t2k.org/nd280/datacomp/meetings/2015/soft_comp_nov23/genie6pissues) with some useful caveats learned from P6 GENIE generation.

# GENIE

Tom Dealtry has written a nice python script that installs all prereqs for GENIE and GENIE itself.  It is located in the `nd280Computing` directory in the T2K Repository under `GENIE_install_scripts`.

The command to run `installGENIE.py` is:
```
./installGENIE.py -g 2.8.0 -l 5.9.1 -p 6_428 -r 5.30.06 -x 2.9.0 -m 'Make.include' -u 'UserPhysicsOptions.xml' -d /nfs/hepcode/t2k/ -o
```
For an explanation of the options do `./installGENIE.py -h`.

Compiling ROOT can be problematic due to `/usr/local/lib/libfftw3.a` on being a statically linked library instead of a dynamically linked lib.  Options seem to be either recompile fftw3 library with `-fPIC` option, or disable fftw3 in ROOT when configuring via `./configure --disable-fftw3`.

Once GENIE has been installed, there is a setup script created, and the scripts need to know the location of that (currently genieSetupPath).

# GENIE Flux Probs Production
Download the flux release from the `nd280/beam/MC` dir on the lcg.  Untar.  You will find many files that contain the flux from the beam group.  It's best to hadd these files into a single large flux file that's used for the vector production.

For an example, the `/nfs/data44/t2k/GNProd/production006/flux` directory contains the hadded flux files that are used to produce genie flux probs.  The flux probs are saved in the `genie_evtrate/` directory.  In general, when there is a new flux-release, download all the flux files from the grid for ND5/ND6 and hadd them each into a large root file using haddFlux.py.  Then run `makeGENIEsetup.py` to create the nd280.cfg files necessary to generate the flux probs.  Finally, run submit_scripts.py (in `~tianlu/software/Python/tianlu/bin/`) to submit the flux prob generation to the cluster.  For flux files on the order of a few GB this will take a very long time (we're tlaking days here).

edited T.Y. 5/12/14
python executables needed are:
`haddFlux.py`, `makeGENIEsetup.py`, `submit_scripts.py`

All the executables can be found in `flux/bin/`.  Symlink or put in your $PATH to use.

## Note on flux-prob generation for TPC-Gas sample
Since interactions on Ar are extremely rare, a large PoT value is needed if we want to cherry-pick a statistically significant sample of Ar-interactions.  A better way is to set `master_volume = Half` in the cfg file and regenerate the GENIE-MC chain (i.e. nd280Control modules `genieSetup`, `genieMC`, and `genieConvert`).  Several issues arose while doing this for prod 6B.  In the genieSetup and genieMC modules the `master_volume` is read and passed to `$GENIE/bin/gT2Kevgen` with the `-t` option.  If there is no '+' prefix on the master_volume, the resulting Ar-interactions are restricted to half of TPC and the interaction rate is also incorrect.  We want a uniform distribution across all three TPCs, and this can be done by appending a '+' in front of the `master_volume` (e.g. `gT2Kevgen -t +Half`).  I have hardcoded this into both `genieSetup.py` and `genieMC.py` and the tpcgas production can proceed.  However, the processing of flux-probs is finicky with '+' prefixed as genieSetup jobs with certain flux files will get stuck and seemingly never end.  This means if you hadd a large number of flux files, most likely the job to process the flux-probs will get stuck and never finish.  I waited for more than a day for a job to finish that should have finished in a few hours.  This is most likely a bug in GENIE, as the '+' prefix is explained in detail in `gT2Kevgen.cxx` and should be included.  Something somewhere is getting caught in an infinite loop.  I am leaving `genieSetup.py` and `genieMC.py` as is because it _should_ work and to enforce consistency with the NEUT modules, but for future productions, this may be an issue that needs to be checked deeper.

# Directory Creation, Directory Layout, and Cards
The directory structure more or less follows the style laid out [here](http://www.t2k.org/nd280/datacomp/howtoaccessdata/directorystructure)
This means that the structure more or less is the same as it will be in the LFC on the GRID, e.g. <base>/production005/C/mcp/genie/2010-11-air/magnet/beama/anal reco ...

There is an additional subdirectory, 'staging', which is the 'work area' to  facilitate production and organization. This staging subdirectory contains subfolders 'numc', 'nd280', and 'nucp', for use in creation of numc, nucp, or nd280 files. Each of these subdirectories has subfolders 'cfg', 'scripts', and 'logs' which should be self explanatory.

Creation of the directory structure is handled via the 'CreateDirs.py' script.  It requires a 'card' file (the sample card is named 'dirs.card') that lets the
script know which respin, the mc type, the sotware version, and whether it's a verify run or not. It is run from the production direcotry (e.g. <base>/production005/CreateDirs.py), and it creates all of the relevant subdirectories for production for a given respin. N.B.: not all subdirectories are necessarily used. The 'nucp' directories are only used for cherry picking, for example.  It will also copy the runInfo.card to all the staging directories.

A separate type of card is required within each staging directory, giving the desired output file type (numc, nucp, nd280) and the desired input file type (None, numc, nucp), as well as the software version. This could be extended to more comlicated scenarios for smaller respins, such using cali files to respin repo and analysis. The name of this card is runInfo.card

# GENIE Vector Production
First ensure that the `flux/` directory containing the GENIE flux probs exists in the correct location.  By default it's based off `self.prodDir` in MCJobClass.py but you can also set the "fluxdir" parameter in the runInfo.card to override the default value.

The production is entirely based on teh `runND280` executable that is built when building the ND280 software.  This takes as input a cfg file with numerous different control parameters.  For an explanation see this [link](http://www.hep.lancs.ac.uk/nd280Doc/devel/invariant/nd280Control/nd280ControlOptions.html).  Everything regarding the production MUST be run within the `staging/` directory!

What the scripts `CreateCfgFile.py` does is finds the `runInfo.card` file in your cwd and creates cfg files based on the values set there.  For instance, you can modify the card file to have `output: numc`, which will tell `CreateCfgFile.py` to create a cfg file suitable for making the vector files.  You can run it by calling `CreateCfgFile.py 0 0` which creates the cfg file for run 0, subrun 0.  The cfg file is stored in `blah/staging/numc/cfg`.  Inspect this to ensure it's correct.

To then create the script that is submitted to the SGE batch, use `SubmitND280Analysis.py 0 0` and inspect the output in `blah/staging/numc/scripts`.

# ND280 MC production
By now things should be clearer as to how these scripts are meant to be used.  Once the vector files have been produced, the full nd280mc will be based off those vector files.  Modify the runInfo.card to reflect the fact that `output: nd280` and ensure the `numcthroughdir` is correctly set.  By default it will look in the `production00n/A` directory.  The scripts will complain if it cannot file the correct numc file for the run-subrun.  The Create and Submit scripts are run with the same arguments.

Note that creation of the .cfg file is based off of both the runInfo.card, which sets some basic requirements, _and_ the current directory full-path.  I.e. the directory structure is important for creating the cfg files (e.g. 2010-11-air/magnet/run3/ will tell the configuration class to set the correct settings for run3).  Since these settings are hard-coded, and they sometimes change, it's important to ensure that things like the interactions-per-spill is correctly set.  The int/spill _must_ be updated whenever the flux probs change.  The output logs from the GENIE flux creation contains info the interaction rate.  For magnet files, the interactions per spill is calculated as Int/Sp = (Evts/POT) * (POT/Sp).  The POT/Sp is calculated and listed on the [MCProdSummary page](http://www.t2k.org/nd280/datacomp/Production006/mcp/mcProdSummary), and the Evts/POT is written out to the log files of the `genie_setup` module from generating the flux-probs.

One thing to be aware of is that sometimes the jobs will fail, or there will be errors with the passthrough.  These are bugs within the code that seem to be triggered by the random seed or may be due to job failures on the cluster itself.  We must not probe too deeply.  The workaround is to use `CountDirMissing.py` on the output `anal/` directory and see if things are amiss and use `Purge.py` to remove traces of missing/incomplete sets of files.  These scripts have helper functions that explain their usage `Purge.py -h`.

# Cherry Picking
Cherry picked files start with basket beam-numc vector files and select specific interactions from the beam gRooTracker.  It utilizes [oaCherryPicker](http://www.hep.lancs.ac.uk/nd280Doc/devel/invariant/oaCherryPicker/) by setting the `[cherry_picker]` module in cfg and passing that to runND280.  To create nucp files, first generate the RooTracker numc files with adequate PoT.  With that as input, set output in runInfo.card to "nucp" and provide the correct "passthrudir" that contains the necessary numc file.  Outputs are stored in `nucp/` and are used when running the full production chain.
