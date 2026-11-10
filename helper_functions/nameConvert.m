% Used for: biome name conversion (used by the biomes_*_ot.m scripts)
% Original: Code\N2_NEMO\nameConvert.m (data paths updated to I:\Prj; computation unchanged)
function name = nameConvert( area )
% convert area name to biome name
% Haipeng Zhao
% 10/10

if strcmp( area, 'Atlantic EQU' )
    name = 'A EQU';
elseif strcmp( area, 'Atlantic STPS' )
    name = 'NA STPS';
elseif strcmp( area, 'Atlantic STSS' )
    name = 'NA STSS';
elseif strcmp( area, 'Eest Pacific EQU' )
    name = 'EP EQU';
elseif strcmp( area, 'Indian Ocean' )
    name = 'IND STPS';
elseif strcmp( area, 'North Atlantic' )
    name = 'NA SPSS';
elseif strcmp( area, 'North Pacific')
    name = 'NP SPSS';
elseif strcmp( area, 'Pacific Arctic' )
    name = 'NP ICE';
elseif strcmp( area, 'Pacific STPS' )
    name = 'NP STPS';
elseif strcmp( area, 'Pacific STSS' )
    name = 'NP STSS';
elseif strcmp( area, 'S Pacific  STPS' )
     name = 'SP STPS';
elseif strcmp( area, 'SH SPTS')
    name = 'SA STPS';
elseif strcmp( area, 'SH ICE' )
    name = 'SO ICE';
elseif strcmp( area, 'SH SSPS' )
    name = 'SO SPSS';
elseif strcmp( area, 'SH STSS' )
    name = 'SO STSS';
elseif strcmp( area, 'West Pacfic EQU')
    name = 'WP EQU';
else
    name = area;
end