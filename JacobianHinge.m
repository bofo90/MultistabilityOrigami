function [J,t]=JacobianHinge(p0)
  J=0;
%Node coordinates
p=reshape(p0',12,1);
%Vectors
a=p(4:6)-p(1:3);    %Rotation axis
b=p(7:9)-p(1:3);    %Vector on face 1
c=p(10:12)-p(1:3);  %Vector on face 2

da=[-1 0 0
    0 -1 0
    0 0 -1
    1 0 0
    0 1 0
    0 0 1
    0 0 0
    0 0 0
    0 0 0
    0 0 0
    0 0 0
    0 0 0]';
db=[-1 0  0
     0  -1 0
     0  0 -1
     0  0  0
     0  0  0
     0  0  0
     1  0  0
     0  1  0
     0  0  1
     0  0  0
     0  0  0
     0  0  0]';
dc=[-1 0  0
     0  -1 0
     0  0 -1
     0  0  0
     0  0  0
     0  0  0
     0  0  0
     0  0  0
     0  0  0
     1  0  0
     0  1  0
     0  0  1]';
 
ka=norm(a);
na=a/ka;

nab = crossvector(a,b);
nca = crossvector(c,a);
kab = norm(nab);
kca = norm(nca);

nab = nab/kab;
nca = nca/kca;

detf=na'*(crossvector(nab,nca));
dotf=nab'*nca;
t = atan2(detf,dotf);

dna=da/ka-a*(a'*da)/ka^3;
dkab=1/(kab)*((da'*a)*(b'*b)+(a'*a)*(db'*b)-(a'*b)*((da'*b)+(db'*a)));
dkca=1/(kca)*((dc'*c)*(a'*a)+(c'*c)*(da'*a)-(c'*a)*((dc'*a)+(da'*c)));
dnab=1/kab^2*((crossvector(da,b)+crossvector(a,db))*kab-crossvector(a,b)*dkab');
dnca=1/kca^2*((crossvector(dc,a)+crossvector(c,da))*kca-crossvector(c,a)*dkca');
 
ddetf=dna'*(crossvector(nab,nca))+(na'*(crossvector(dnab,nca)+crossvector(nab,dnca)))';
ddotf=dnab'*nca+dnca'*nab;

J = -(-dotf*ddetf+detf*ddotf)/(detf^2+dotf^2);

function w=crossvector(u,v)
w=[u(2,:).*v(3,:)-u(3,:).*v(2,:);-u(1,:).*v(3,:)+u(3,:).*v(1,:);u(1,:)*v(2,:)-u(2,:)*v(1,:)];