alpha = 5000000;
beta = 0.00000000;
junk=0.0; samprate=10; dt=1/samprate;

cd synth; load radial.asc; load vert.asc; cd ..;
vert = vert(:,2); rad = radial(:,2);

% What is Q model?
load q.dat; qs = q(:,1); qdz = q(:,2);
nlayer = length(qs);

% Setup time axis
t = 1/samprate:1/samprate:(2*length(vert)-1)/samprate;

% normalize the input data by dividing by the rms
dum1=0.0;
for i=1:length(vert)
  dum1=dum1+vert(i).^2;
end
rms_vert=sqrt(dum1)./length(vert); vert=vert./rms_vert;
dum2=0.0;
for i=1:length(rad)
  dum2=dum2+rad(i).^2;
end
rms_rad = sqrt(dum2)./length(rad); rad =rad./rms_rad;

% build matrices for inversion. R is from the vertical and built every time.
disp('building inversion matrix')

n=2*(length(vert))-1;

% calculate frequencies at which fft will be evaluated
cnt=1;
for ii=n/2-1:-1:-n/2
  freq(cnt) = ii/(n.*dt);
  cnt=cnt+1;
end

% vertical is padded with zeros to build inversion matrix
newvert=zeros(1,n);
for i=1:(length(vert))
  newvert(i)=vert(i);
end

% memory allocations	
toe = zeros( (length(newvert)), (length(vert)+50));
R = zeros( (length(newvert)), (length(vert)+50));
filt_trace=zeros(length(newvert),1);

disp('Start filtering ...')
for i=1:length(vert)+50;
%---------------- filter the vertical ----------------
  lag=t(i); tstar2;
  V1 = fft(newvert); 
  V1_shift=fftshift(V1); 
  V1_shift_tstr = tstr.*V1_shift; 
  V1_shift_tstr_ishift = ifftshift(V1_shift_tstr);
  filt_trace = real(ifft(V1_shift_tstr_ishift));
%-----------------------------------------------------
% now normalize the trace by its rms
  dum3=0.0;
  for ii=1:length(filt_trace)
    dum3=dum3+filt_trace(ii).^2;
  end
  rms_filt_trace = sqrt(dum3)./length(filt_trace);
  filt_trace     = filt_trace./rms_filt_trace;

% place filtered seismogram into the j-th column of the matrix
% Note that 50 samples has been hardwired in for pre-lag
  for j=i:n-50  
    toe(j,i) = filt_trace(j-i+50);
  end       % end placing filtered seismogram into column

end         % end building matrix
disp('End filtering ... set up inversion')

R=toe;

% build damping matrix
W=alpha*(eye(length(vert)+50,length(vert)+50));

% pad the radial with zeros to make it the right size
padrad=zeros(n,1);
for i=1:length(vert)
  padrad(i)=rad(i);
end

%--------------- invert ---------------
disp('multiplying')
left=(R'*R)+W; 
right=R'*padrad;
disp('inverting....')
%left_inv = inv(left);
%recf = left_inv*right;
%res = left_inv*(R'*R);
recf=bicg((left),(right),0.001,500);  

%save rf5e7 recf -ascii