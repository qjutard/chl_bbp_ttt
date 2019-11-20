# FCHLA and BBP TREATMENT
<br> Code for processing the Fchla and bbp variables from the BGC Argo netcdf files
<br> Main code: processing.files.R (range test, offset, NPQ, F2 for the chla ; range test, drift for the bbp) 
<br> Functions:
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








