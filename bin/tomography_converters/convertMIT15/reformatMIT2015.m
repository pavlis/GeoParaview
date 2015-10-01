% Converter script for MIT 2015 model to output regular grid.  That output
% is converted to format similar to earlier MIT models for input into
% Pavlis's GCLgrid format with a C++ conversion program written 
% for earlier MIT models
path('./mbin',path);

% Define lat, lon, and depth grids 
lonsraw = [-150, -50];
lons=lonsraw+360;
lats = [10, 70];
depths = [0,1000.0];


% Set up lon and lat grids in radians (from Scott original script)
ax = [lons lats]*pi/180; 
longr=linspace(ax(1),ax(2),round(4*(ax(2)-ax(1))/pi*180+1))';
latgr=linspace(ax(3),ax(4),round(4*(ax(4)-ax(3))/pi*180+1))';
d=linspace(0,1000,21);
nz=max(size(d));
% load model
dV=loadmitmodel15;
% Scott's getlayers script returns a vector of matrices given
% an input vector radii.   L contains that result after the next 
% two lines
r=6370.0-d;
L=getlayers(r,dV,longr,latgr);
% All layers are same size matrix so we use the first one to get 
% dimensions
[nrow,ncol]=size(L{1}.A);
% We write results to this fixed file name
fd=fopen('MIT2015.txt','w');

% nrow is same as size latgr so rows are latitude grid
% ncol is then longitude grid variable
% Old MIT models had depth range as fastest index followed
% by latitude from low to high and last by lon low to high
% We replicate that here
for j=1:ncol
    for i=1:nrow
        fprintf(fd,'%f %f ',longr(j)*180.0/pi,latgr(i)*180.0/pi);
        for k=1:nz
            val=L{k}.A(i,j);
            fprintf(fd,' %f',val);
            fprintf(fd,'\n');
        end
    end
end
fclose(fd);

