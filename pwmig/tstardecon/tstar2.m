% M-file to compute t* for filtering ...

tmp_time=0.0; sum=0.0;
%tstr=zeros(n,1); 
for k=1:nlayer
  if lag <= qdz(k)
    sum=sum+(lag-tmp_time)./qs(k);
    break;
  else
    sum=sum+(qdz(k)-tmp_time)./qs(k);
    tmp_time=qdz(k);
  end
end

% Now compute the filter - call it tstr
for kk=1:n
  tstr(kk)=exp(-pi.*abs(freq(kk)).*sum);
end

