function [Lx,Lz,mpt] = readmpt(filename)
  %Loads mean-pol-tor formulation from file
  alpha = ncreadatt(filename,'/','alpha');
  gamma = ncreadatt(filename,'/','gamma');
  mpt = ncread(filename,'mpt');
  Lx = 2*pi/alpha;
  Lz = 2*pi/gamma;
end
