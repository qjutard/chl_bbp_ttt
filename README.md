# DMMC : FCHLA and BBP delayed mode
Software for processing the Fchla and bbp variables from the BGC Argo netcdf files and create Delayed Mode (DM) files following the work done by M. Cornec in Bellacicco et al. 2019 (http://dx.doi.org/10.1029/2019GL084078)

## Recommended "installation" process :
* Clone the repository where you want it ( `$ git clone https://github.com/qjutard/chl_bbp_ttt` ) or download the release you want to use
* Get the latest profile index and greylist if necessary
* Adapt pathway definitions in *pathways.R* and *start_DMMC.R*
* Create an alias in your *.bashrc* or *.bash_aliases* ( `alias DMMC="~/path/to/repository/DMMC.sh"` )
* Go to the working directory (where you want log outputs to be written)
* READ THE HELP ( `$ DMMC -h` )
* Updates can then be obtained with `$ git pull` if you cloned the repository
* Feel free to contact me directly for help and to open issues for bugs or desired features


## Quick file descriptions
* **./argo_merge-profile_index.txt** : index of argo profiles in coriolis
* **./ar_greylist.txt** : greylist
* **./DMMC.sh** : shell script that manages input options
* **./error_message.R** : gives the error message associated to an error index
* **./increment_N_CALIB** : increments the N_CALIB dimension in a BGC Argo necdf file
* **./pathways.R** : local pathways to define
* **./process_files.R** : computes delayed mode for one profile (range test, offset, NPQ, F2 for the chla ; range test, drift for the bbp)
* **./start_DMMC.R** : is called by *DMMC.sh* and translates the options to call the delayed mode operations
* **./write_DM.R** : writes DM information for a given parameter on a file
* **./write_DM_with_mcornec.R** : takes one profile as input, calls *process_files.R* to obtain DM data, then creates surrounding DM informations to write with *write_DM.R*
* **./Functions** : functions called by *process_files.R* while computing the delayed mode operations
  * **/Darkoz.R** : calculate a dark offset in the regions of Oxygen Minimum Zone
  * **/DarkXing.R** : calculate a dark offset in the regions where the chla show an increase at depth
  * **/Dark_Fchla_Corr.R** : correct the dark offset of the chla
  * **/Dark_MLD_table_coriolis.R** : correct the dark offset in cases of deep vertical mixing
  * **/MLD_calc.R** : calculate the MLD
  * **/NPQ_corX12_XB18.R** : correct the profile from the NPQ (in presence of PAR measured in situ)
  * **/NPQ_cor_P18.R** : correct the profile from the NPQ (in absence of PAR measured in situ)
  * **/Outliars_med.R** : test outliars in a time series
  * **/RunningFilter.R** : smooth profiles (running median/mean on a given window)
  * **/Zone.R** : regional criterion on the profile
* **./Tests** : tools for automated testing of the software, intended for development purposes
  * **/attach_commit.sh** : writes the latest commit sha to an error report as produced by *run_tests.sh*
  * **/post_report.sh** : sends an error report to github
  * **/run_tests.sh** : calls the commands from *test_runs.t* and writes an error report, can be immediately linked to the current commit
  * **/test_runs.t** : defines a series of DMMC calls to try








