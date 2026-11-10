% Used for: monthly mean of the daily LLC270 MLD
% Original: Code\mld_ECCO\Compute_MLD\monthly_ot.m (data paths updated to I:\Prj; computation unchanged)
% Output monthly MLD
% Haipeng Zhao
% 8/26/2025

clear
clc


for yyyy = 2003 : 2022
    dataFOLDER = ['I:\ECCO-Darwin_extension\daily\mld_003\Daily', '\', num2str( yyyy )];
    saveFOLDER = ['I:\ECCO-Darwin_extension\daily\mld_003\Monthly', '\', num2str( yyyy )];
    
    numFaces = 13;
    nx = 270;
    ny = 270 .* numFaces;
    nz = 50;
    
    fileList = dir( fullfile( dataFOLDER, 'mld*.data' ) );
    nFile = length( fileList );
    
    daysInMonth = arrayfun( @(m) eomday( yyyy, m ), 1 : 12 );
    monthEdges = [0 cumsum( daysInMonth )]; %
    
    for i = 1 : 12
        startDay = monthEdges(i) + 1;
        endDay   = monthEdges(i+1);
    
        mld = nan( nx, ny, endDay - startDay + 1 );
        k = 1;
        for j = startDay + 1 : endDay
            tmp = readbin( fullfile( fileList(j).folder, fileList(j).name ), [nx ny], 1, 'real*4' );
            mld(:, :, k) = tmp;
            k = k + 1;
        end
        mld_avg = nanmean( mld, 3 );
    
        writebin( fullfile( saveFOLDER,  sprintf('mld_%d_%02d.data', yyyy, i) ), mld_avg );
    end
end

