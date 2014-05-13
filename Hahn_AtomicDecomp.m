
function Hahn_AtomicDecomp()

clear all

load patches.mat

load dict512.mat

D=Wp';


for i = 4 : 4
    
    y=data(:,i);
    
    
    yy=reshape(y,16,16);
    
    
    figure(1)
    subplot(411)
    imagesc(yy,[0,.5])
    %   colormap(gray)
    
    a=BP(10,D,y);
    
    
    subplot(412)
    imagesc(reshape(D*a,16,16),[0,.5])
    %   colormap(gray)
    
    
    subplot(413)
    
    imagesc(yy-reshape(D*a,16,16),[0,.5])
    %   colormap(gray)
    
    
    subplot(414)
    hist(a,100)
%     pause
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    figure(2)
    colormap(jet)
    subplot(411)
    imagesc(yy,[0,.5])
    %   colormap(gray)
    
    a=MP(10,D,y);
    
    
    subplot(412)
    imagesc(reshape(D*a,16,16),[0,.5])
    %   colormap(gray)
    
    
    subplot(413)
    
    imagesc(yy-reshape(D*a,16,16),[0,.5])
    %   colormap(gray)
    
    
    subplot(414)
    hist(a,100)
%     pause



 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
    figure(3)
    colormap(jet)
    subplot(411)
    imagesc(yy,[0,.5])
    %   colormap(gray)
    
    a=LCA(y,D,0.01);
    
    
    subplot(412)
    imagesc(reshape(D*a,16,16),[0,.5])
    %   colormap(gray)
    
    
    subplot(413)
    
    imagesc(yy-reshape(D*a,16,16),[0,.5])
    %   colormap(gray)
    
    
    subplot(414)
    hist(a,100)
%     pause
    
    
    
    
end

% norm(y-reshape(D*a,16,16))

end





function [a, u] = LCA(y, D, lambda)


t=.01;
h=.0001;

d = h/t;
u = zeros(size(D,2),1);


for i=1:300
    
    
    a = ( u - sign(u).*(lambda) ) .* ( abs(u) > (lambda) );
    
    
    u =   u + d * ( D' * ( y - D*a ) - u - a  ) ;


end




end

















function [a] = MP(k,D,y)

n=size(D);
D1=zeros(n);
r=y;

for k=1:k
    
   
    [i,j] = max(abs(D'*r));
    
    D1(:,j)=D(:,j);
    
    D(:,j)=0;
    
    a = D1 \ y;
    
    r = y - D1*a;
    
    
end;




end


















function [a] = BP(k,D,y)


cvx_begin quiet;

    l=size(D);

    variable a(l(2));

    minimize( norm(D*a-y) );

    subject to;

    norm(a,1) <= k;

cvx_end;



end










