% simple m file to generate fake array of 100 stations randomly positioned
% but with a specified cutoff

cutoff=0.05;
x=zeros(100,2);
x(1,:)=rand(1,2);
nset=1;
for i=nset:100
    tryagain=1;
    while tryagain > 0
        xtest=rand(1,2);
        tryagain=0;
        for j=1:nset
            dx1=xtest(1,1)-x(j,1);
            dx2=xtest(1,2)-x(j,2);
            r=hypot(dx1,dx2);
            if r<cutoff 
                tryagain = 1;
                break;
            end
        end
    end
    nset=nset+1;
    x(nset,:)=xtest;
end
plot(x(:,1),x(:,2),'o');
scale=50.0
dx1=-0.5;
dx2=-0.5;
for i=1:100
    x(i,1)=scale*(x(i,1)+dx1);
    x(i,2)=scale*(x(i,2)+dx2);
end
figure;
plot(x(:,1),x(:,2),'o');