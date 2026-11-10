% Used for: daily MLD on the LLC270 grid (0.03 kg/m3 threshold, uses gcmfaces)
% Original: Code\mld_ECCO\Compute_MLD\compute_MXLDEPTH_NAS.m (data paths updated to I:\Prj; computation unchanged)
function [] = compute_MXLDEPTH_NAS(dataDir,saveDir)

cd(saveDir)
delete *

% gridDir = '<cluster_dir>/grid/LLC_270/';
gridDir = 'I:\Prj\MLD_ECCO\LLC_270\grid\';

%%
%load llc270 grid

numFaces = 13;
nx = 270;
ny = 270 .* numFaces;
nz = 50;

gcmfaces_global;

global mygrid;

mygrid = [];

grid_load( gridDir, 5, 'compact' ); % Load llc270 grid


llc270.YC = readbin([gridDir 'YC.data'],[nx ny],1,'real*4');
llc270.hFacC = readbin([gridDir 'hFacC.data'],[nx ny nz],1,'real*4');

llc270.hFacC(llc270.hFacC == 0) = nan;

llc270.RAC = readbin([gridDir 'RAC.data'],[nx ny],1,'real*4');
llc270.drC = repmat(readbin([gridDir 'DRC.data'],[nz],1,'real*4'),1,nx,ny);
drC = repmat(readbin([gridDir 'DRC.data'],[nz],1,'real*4'),1,nx,ny);
llc270.drC = permute(drC,[2 3 1]);

thetaFiles = dir([dataDir,'THETA.0*.data']);
saltFiles = dir([dataDir,'SALTanom*.data']);

c = 1;

for i = 1:length(thetaFiles)
    
    fileName = thetaFiles(i).name;

    THETA = convert2gcmfaces(readbin([dataDir thetaFiles(i).name],[nx ny nz]));
    SALT = convert2gcmfaces(readbin([dataDir saltFiles(i).name],[nx ny nz]) + 35);

    fldT = THETA .* mygrid.mskC; 
    fldS = SALT .* mygrid.mskC;
       
    %prepare to compute potential density
    fldP = 0 .* mygrid.mskC; 
        
    for kk = 1:length(mygrid.RC)
        
    	fldP(:,:,kk) = -mygrid.RC(kk); 
        
    end
        
    T = convert2vector(fldT);
    S = convert2vector(fldS);
    msk = convert2vector(mygrid.mskC);
    P = convert2vector(fldP);
        
    %compute potential density    
    RHO = 0 .* msk; 
    alpha = 0 .* msk;
        
    tmp1 = find(~isnan(msk));
    RHO(tmp1) = density(T(tmp1),S(tmp1),P(tmp1));
        
    fldRhoPot = convert2vector(RHO);
    alpha(tmp1) = density(T(tmp1)+1e-4,S(tmp1),P(tmp1));
    fldAlpha = (convert2vector(alpha)-fldRhoPot)/1e-4;

    clear T S P msk RHO RHOis tmp1;

    %compute mld
    tmp1 = NaN .* mygrid.mskC(:,:,1);
        
    for kk = 1:50
        
    	tmp2 = fldRhoPot(:,:,kk) - fldRhoPot(:,:,1);
          
    	%if we pass RHO(1)+0.03 for the first time (or we reach the bottom)
    	%then mld is the velocity point above RC(kk), which is RF(kk)
          
    	jj = find((tmp2>0.03 | isnan(tmp2)) & isnan(tmp1));
    	tmp1(jj) = -mygrid.RF(kk);
        
    end
        
    fldMldBoyer = convert2gcmfaces(tmp1);
    clear tmp1 tmp2

    %compute mld
    tmp1 = NaN .* mygrid.mskC(:,:,1);
    
    for kk = 1:50;
    	
        tmp2 = fldRhoPot(:,:,kk) - fldRhoPot(:,:,1);
        
        %if we pass RHO(1)+0.125 for the first time (or we reach the bottom)
        %then mld is the velocity point above RC(kk), which is RF(kk)
        
        jj = find((tmp2 > 0.125 | isnan(tmp2)) & isnan(tmp1));
      
        tmp1(jj) = -mygrid.RF(kk);
   
    end

    fldMldSuga=convert2gcmfaces(tmp1);
    clear tmp1 tmp2

    %compute mld
    tmp1 = NaN .* mygrid.mskC(:,:,1);
   
    fldRhoPotMax = fldRhoPot(:,:,1) - 0.8 .* fldAlpha(:,:,1);
    
    for kk = 1:50;
    
  	tmp2 = fldRhoPot(:,:,kk) - fldRhoPotMax;
     
        %if we pass RHO(1)+0.8*alpha(1) for the first time (or we reach the bottom)
        %then mld is the velocity point above RC(kk), which is RF(kk)
    
        jj = find((tmp2 > 0 | isnan(tmp2)) & isnan(tmp1));
        tmp1(jj) = -mygrid.RF(kk);
    
    end

    fldMldKara = convert2gcmfaces(tmp1);
    clear tmp1 temp2

    writebin([saveDir 'MXLDEPTH_boyer.' fileName(7:end-5) '.data'],fldMldBoyer);
    writebin([saveDir 'MXLDEPTH_suga.' fileName(7:end-5) '.data'],fldMldSuga);
    writebin([saveDir 'MXLDEPTH_kara.' fileName(7:end-5) '.data'],fldMldKara);

    c = c + 1;
 
    disp(num2str(c));
    
    clear THETA SALT fldT fldS fldMldBoyer fldMldSuga fldMldKara tmp1 tmp2

end

end

