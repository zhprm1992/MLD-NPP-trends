% Used for: MLD of one T/S profile, 0.03 kg/m3 sigma0 threshold from 10 m
% Original: Code\mld_ECCO\Compute_MLD\findMLD_sigma003.m (data paths updated to I:\Prj; computation unchanged)
function mld_dens = findMLD_sigma003( temp, sal, pres )
% This subroutine calculates the mixed layer depth (MLD) of the input profile.  
% Output:
% - mldepthdensp: Density threshold 0.03 kg/m3
% - mldepthptmpp: Temperature threshold 2 C
% - gdmldp: Density gradient method
% - mixeddp: Holte's method
%
% Rich Signell
% 6/12/2015
%
% Haipeng Zhao
% 2/12/2025


% Calculate the MLD using a threshold method with de Boyer Montegut et al's
% criteria; a density difference of .03 kg/m^3 or a temperature difference
% of .2 degrees C.  The measurement closest to 10 dbar is used as the 
% reference value.  The threshold MLDs are interpolated to exactly match 
% the threshold criteria.
%

% Calculate the index of the reference value
m = length( sal ); 
starti = find( (pres - 10) .^ 2 == min( (pres - 10) .^ 2 ) , 1 );  
pres = pres( starti : m );
sal = sal( starti : m );
temp = temp( starti : m );

starti = 1;
m = length( sal ); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the potential density anomaly, with a reference pressure of 0 
% pden = sw_pden(sal,temp,pres,0)-1000;
pden = gsw_sigma0( sal, temp );  % potential density anomaly referenced to surface (rrom TEOS10)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Search for the first level that exceeds the potential density threshold
mld_dens = m;  % 
for j = starti : m
    if abs( pden( starti ) - pden( j )) > 0.03 
        mld_dens = j;
        break;
    end
end

% Interpolate to exactly match the potential density threshold
% clear pdenseg presseg presinterp pdenthreshold
presseg = [pres(mld_dens - 1) pres(mld_dens)];
pdenseg = [pden(starti) - pden(mld_dens - 1) pden(starti) - pden(mld_dens)];
P = polyfit( presseg, pdenseg, 1 );
presinterp = presseg(1) : .5 : presseg(2);
pdenthreshold = polyval( P, presinterp );

% The potential density threshold MLD value:
mld_dens = presinterp(find( abs( pdenthreshold ) < 0.03 , 1, 'last' ) );

% Search for the first level that exceeds the temperature threshold
% mld_tmp = m;
% for j = starti : m
%     if abs( temp(starti)-temp(j) )>.2 
%         mld_tmp = j;
%         break;
%     end
% end
% 
% % Interpolate to exactly match the temperature threshold
% clear tempseg presseg presinterp tempthreshold
% presseg = [pres(mld_tmp-1) pres(mld_tmp)];
% tempseg = [temp(starti)-temp(mld_tmp-1) temp(starti) - temp(mld_tmp)];
% P = polyfit(presseg,tempseg,1);
% presinterp = presseg(1) : .5 : presseg(2);
% tempthreshold = polyval(P,presinterp);
% 
% % The temperature threshold MLD value:
% mld_tmp = presinterp(find(abs(tempthreshold)<.2, 1, 'last' ));

end