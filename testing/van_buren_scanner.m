
%write a parser for van buren:

pathn='C:\Users\uqastilg\Downloads\complex_oblate_swf-master\fort.20';
isProlate=false;
% pathn='C:\Users\uqastilg\Downloads\complex_prolate_swf-master\fort.20';
% isProlate=true;

fid = fopen(pathn);
tline = fgetl(fid);
%first line different:
disp(tline)
[splits]=regexp(tline,'\s+','split');
xi=str2num(splits{2})+0;
c=str2num(splits{3})+str2num(splits{4})*1i;
m=str2num(splits{5});

tline = fgetl(fid);
n=[];
R1vb=[];
dR1vb=[];
R2vb=[];
dR2vb=[];
while ischar(tline)
    disp(tline)
    [splits]=regexp(tline,'\s+','split');
    n=[n;str2num(splits{2})];
    R1vb=[R1vb;(str2num(splits{3})+str2num(splits{4})*1i)*10^(str2num(splits{5}))];
    dR1vb=[dR1vb;(str2num(splits{6})+str2num(splits{7})*1i)*10^(str2num(splits{8}))];
    tline = fgetl(fid);

    disp(tline)
    [splits]=regexp(tline,'\s+','split');
    R2vb=[R2vb;(str2num(splits{2})+str2num(splits{3})*1i)*10^(str2num(splits{4}))];
    dR2vb=[dR2vb;(str2num(splits{5})+str2num(splits{6})*1i)*10^(str2num(splits{7}))];
    tline = fgetl(fid);
end
fclose(fid);

% now test toolbox
nmax=max(n);
[R1,dR1mod]=spheroidalR1(isProlate,nmax+mod(nmax+m,2),m,c,xi,true); %even
[R1o,dR1modo]=spheroidalR1(isProlate,nmax+mod(nmax+m,2)+1,m,c,xi,true); %odd

[R2,dR2mod,~,rnme]=spheroidalR2(isProlate,nmax+mod(nmax+m,2),m,c,xi,true);
[R2o,dR2modo,~,rnmo]=spheroidalR2(isProlate,nmax+mod(nmax+m,2)+1,m,c,xi,true); %odd

mc=min([size(rnme,1),size(rnmo,1)]);
ms=min([size(rnme,2),size(rnmo,2)]);
rnm=zeros([mc,ms]);
rnm(1:2:2*mc-1,1:ms)=rnme(1:mc,1:ms);
rnm(2:2:2*mc,1:ms)=rnmo(1:mc,1:ms);

figure(9)
imagesc(log10(abs(rnm(1:size(R2vb,1),1:min([size(rnm,2),50]))-R2vb)./abs(R2vb)))


[val,indx]=min(log10(abs(rnm(1:size(R2vb,1),:)-R2vb)./abs(R2vb)),[],2);
figure(19)
plot(indx)


R2c=zeros(mc,1);
R2c(1:2:2*mc-1)=R2(1:mc);
R2c(2:2:2*mc)=R2o(1:mc);

[R2c(1:length(R2vb)),R2vb]

