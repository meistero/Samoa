function [A,q,rhs,sol,ev,ev_it]=  importLSE(file,nUnkn)


A=zeros(nUnkn);
q=zeros(nUnkn,1);
rhs=zeros(nUnkn,1);

fid=fopen(file)

fgetl(fid);

for i=1:nUnkn
    str=fgetl(fid);
    cellLine=textscan(str,'%f');
    y=cell2mat(cellLine);
    A(i,1:nUnkn)=y;
end

fgetl(fid);
fgetl(fid);

for i=1:nUnkn
    str=fgetl(fid);
    cell=textscan(str,'%f');
    elem=cell2mat(cell);
    q(i)=elem;
end

fgetl(fid);
fgetl(fid);

for i=1:nUnkn
    str=fgetl(fid);
    cell=textscan(str,'%f');
    elem=cell2mat(cell);
    rhs(i)=elem;
end

sol=linsolve(A,rhs);

%ev=eigs(A);
ev=zeros(nUnkn);
ev_it=zeros(nUnkn);

%ev_it=eigs(eye(nUnkn)- inv(diag(diag(A)))*A);
