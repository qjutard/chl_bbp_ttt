# DMMC : FCHLA and BBP delayed mode
<br> Code for processing the Fchla and bbp variables from the BGC Argo netcdf files and create Delayed Mode (DM) files following the work done by M. Cornec in Bellacicco et al. 2019 (http://dx.doi.org/10.1029/2019GL084078)
<br> Usage : Use with DMMC.sh (type './DMMC.sh -h' for help), usage will require adapting some pathways to the local file system (see pathways.R)
<br> pathways.R : local pathways to define
<br> start_DMMC.R : is called by DMMC.sh to itself call the delayed mode operations with the correct options
<br> write_DM_with_mcornec.R : takes one profile as input and writes the DM
<br> process_files.R : computes delayed mode for one profile (range test, offset, NPQ, F2 for the chla ; range test, drift for the bbp) 
<br> increment_N_CALIB : increments the N_CALIB dimension in a BGC Argo necdf file
<br> error_message.R : gives the error message associated to an error index
<br> ar_greylist.txt : greylist
<br> argo_merge-profile_index.txt : index of argo profiles in coriolis
<br> /Functions :
<ul type = "disc">
    <br> <li> DarXing.R : calculate a dark offset in the regions where the chla show an increase at depth </li>
    <br> <li> Dark_Fchla_Corr.R : correct the dark offset of the chla </li>
    <br> <li> Dark_MLD_table_coriolis.R : correct the dark offset in cases of deep vertical mixing  </li>
    <br> <li> Darkoz.R : calculate a dark offset in the regions of Oxygen Minimum Zone </li>
    <br> <li> MLD_calc.R : calculate the MLD </li>
    <br> <li> NPQ_cor_P18.R : correct the profile from the NPQ (in absence of PAR measured in situ) </li>
    <br> <li> NPQ_corX12_XB18.R : correct the profile from the NPQ (in presence of PAR measured in situ) </li>
    <br> <li> Outliars_med.R : test outliars in a time series </li>
    <br> <li> RunningFilter.R : smooth profiles (running median/mean on a given window)  </li>
    <br> <li> Zone.R : regional criterion on the profile </li>
</ul>








