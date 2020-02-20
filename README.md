# DMMC : FCHLA and BBP delayed mode
<br> Software for processing the Fchla and bbp variables from the BGC Argo netcdf files and create Delayed Mode (DM) files following the work done by M. Cornec in Bellacicco et al. 2019 (http://dx.doi.org/10.1029/2019GL084078)

<br> Recommended "installation" process :
<ul>
    <li> Clone the repository where you want it (git clone https://github.com/qjutard/chl_bbp_ttt) or dowload the release you want to use </li>
    <li> Get the latest profile index and greylist if necessary </li>
    <li> Adapt pathway definitions in 'pathways.R' and 'start_DMMC.R' </li>
    <li> Create an alias in your '.bashrc' or '.bash_aliases' (alias DMMC="~/path/to/repository/DMMC.sh") </li>
    <li> Go to the working directory (where you want log outputs to be written) </li>
    <li> READ THE HELP (DMMC -h) </li>
    <li> Updates can then be obtained with git pull if you cloned the repository </li>
    <li> Fell free to contact me directly for help and to open issues for bugs or wanted features </li>
</ul>

<br> File descriptions :
<ul type = "none">
<li> ./pathways.R : local pathways to define </li>
<li> ./start_DMMC.R : is called by DMMC.sh and translates the options to call the delayed mode operations </li>
<li> ./write_DM_with_mcornec.R : takes one profile as input, calls process_files.R to obtain DM data, then creates surrounding DM informations to write with write_DM.R </li>
<li> ./write_DM.R : writes DM information for a given parameter on a file </li>
<li> ./process_files.R : computes delayed mode for one profile (range test, offset, NPQ, F2 for the chla ; range test, drift for the bbp) </li>
<li> ./increment_N_CALIB : increments the N_CALIB dimension in a BGC Argo necdf file </li>
<li> ./error_message.R : gives the error message associated to an error index </li>
<li> ./ar_greylist.txt : greylist </li>
<li> ./argo_merge-profile_index.txt : index of argo profiles in coriolis </li>
<li> ./Functions
    <ul type = "none">
        <li> /DarXing.R : calculate a dark offset in the regions where the chla show an increase at depth </li>
        <li> /Dark_Fchla_Corr.R : correct the dark offset of the chla </li>
        <li> /Dark_MLD_table_coriolis.R : correct the dark offset in cases of deep vertical mixing  </li>
        <li> /Darkoz.R : calculate a dark offset in the regions of Oxygen Minimum Zone </li>
        <li> /MLD_calc.R : calculate the MLD </li>
        <li> /NPQ_cor_P18.R : correct the profile from the NPQ (in absence of PAR measured in situ) </li>
        <li> /NPQ_corX12_XB18.R : correct the profile from the NPQ (in presence of PAR measured in situ) </li>
        <li> /Outliars_med.R : test outliars in a time series </li>
        <li> /RunningFilter.R : smooth profiles (running median/mean on a given window)  </li>
        <li> /Zone.R : regional criterion on the profile </li>
    </ul> </li>
</ul>








