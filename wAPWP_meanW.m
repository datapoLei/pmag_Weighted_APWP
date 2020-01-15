function meanW=wAPWP_meanW(polesInputNew,polesRflag)

if nargin<2
    polesRflag=1;
end

% meanW=[lonM,latM,e95a,e95b,omega,Kx,Ky,N];
% polesInputNew=[1ageAll 2ageSpan 3lonP 4latP 5A95 6W 7ResNo 8nj 9Kj];
npts=length(polesInputNew(:,1));
meanW.ageAll=polesInputNew(:,1);
meanW.ageWin=nan(3,1);
meanW.poles=polesInputNew;
meanW.IjWjSUM=zeros(3,3);
meanW.WjNj=0;

Wj=polesInputNew(:,6);

NjTmp=polesInputNew(:,8);
polesInputNewTmp=polesInputNew(:,[3:4 9 8]);

for i=1:length(polesInputNewTmp(:,1))
    % function Ij=wAPWP_pole2I(phi, lambda, K, N)
    Ij=wAPWP_pole2I(polesInputNewTmp(i,1),polesInputNewTmp(i,2),...
        polesInputNewTmp(i,3),polesInputNewTmp(i,4));
    IjWj=Wj(i)*Ij;
    meanW.IjWjSUM=meanW.IjWjSUM+IjWj;
    
    meanW.WjNj=meanW.WjNj+Wj(i)*NjTmp(i);
end

% calculate eigenvalues and eigenvectors of vector sums
[Vtmp,Dtmp] = eig(meanW.IjWjSUM);
if isreal(Vtmp)==0
    Vtmp=real(Vtmp); Dtmp=real(Dtmp);
end
% Vtmp
% polesInputNewTmp

% sort eigenvectors by their eigenvalues
sortFlat='ascend'; sortFlatVmin=1; sortFlatVmax=3;
sortFlat='descend'; sortFlatVmin=3; sortFlatVmax=1;

[D,indEV]=sort(diag(Dtmp),sortFlat);
meanW.eigVal=diag(Dtmp(indEV,indEV));
meanW.eigVec=Vtmp(:,indEV);

[phiMIN,lambdaMIN,rMIN] = cart2sph(meanW.eigVec(1,sortFlatVmin),...
    meanW.eigVec(2,sortFlatVmin),meanW.eigVec(3,sortFlatVmin));
phiMIN=phiMIN*180/pi; lambdaMIN=lambdaMIN*180/pi;
[latA, lonA]=antipode(lambdaMIN,phiMIN); phiMINa=lonA; lambdaMINa=latA;

[phiINT,lambdaINT,rINT] = cart2sph(meanW.eigVec(1,2),...
    meanW.eigVec(2,2),meanW.eigVec(3,2));
phiINT=phiINT*180/pi; lambdaINT=lambdaINT*180/pi;
[latA, lonA]=antipode(lambdaINT,phiINT); phiINTa=lonA; lambdaINTa=latA;

[phiMAX,lambdaMAX,rMAX] = cart2sph(meanW.eigVec(1,sortFlatVmax),...
    meanW.eigVec(2,sortFlatVmax),meanW.eigVec(3,sortFlatVmax));
phiMAX=phiMAX*180/pi; lambdaMAX=lambdaMAX*180/pi;
[latA, lonA]=antipode(lambdaMAX,phiMAX); phiMAXa=lonA; lambdaMAXa=latA;

Kx=abs(sum(meanW.eigVal(1:2))/(sum(meanW.eigVal)-2*meanW.eigVal(1)));
Ky=sum(meanW.eigVal(1:2))/(sum(meanW.eigVal)-2*meanW.eigVal(2));

% meanW=[lonM,latM,e95a,e95b,omega,Kx,Ky,N];
if Kx<Ky
    e95a=140/sqrt(Kx*sum(Wj));
    e95b=140/sqrt(Ky*sum(Wj));
else
    e95a=140/sqrt(Ky*sum(Wj));
    e95b=140/sqrt(Kx*sum(Wj));
end

nptsM=npts;

% function [Dm,Im,alpha95,kappa,N]=FishM(D,I)
[DmTPM,ImTPM,alpha95TPM,kappaTPM,NTPM]=FishM(polesInputNew(:,3),polesInputNew(:,4),polesInputNew(:,5));

if distance(lambdaMIN,phiMIN,ImTPM,DmTPM)>=distance(lambdaMINa,phiMINa,ImTPM,DmTPM)
    
    omega=(phiINTa-phiMINa)/abs(phiINTa-phiMINa)*acosd(sind(lambdaINTa)/cosd(lambdaMINa));
    if isnan(omega)==1 omega=zeros(length(omega),1); end
    if isreal(omega)==0 omega=zeros(length(omega),1); end
    
    meanW.eigDir(1,:)=[phiMINa,lambdaMINa,rMIN,phiINT,lambdaINT,rINT,phiMAX,lambdaMAX,rMAX];
    meanW.eigDir(2,:)=[phiMIN,lambdaMIN,rMIN,phiINTa,lambdaINTa,rINT,phiMAXa,lambdaMAXa,rMAX];
    
    %----------------------------------------------------------------------
    
    meanW.meanW=[phiMINa,lambdaMINa,e95a,e95b,omega,Kx,Ky,nptsM];
    %----------------------------------------------------------------------
    
    if npts==1
        eigDirTmp=[phiMINa,lambdaMINa;phiINT,lambdaINT;phiMAX,lambdaMAX;...
            phiMIN,lambdaMIN;phiINTa,lambdaINTa;phiMAXa,lambdaMAXa];
        for kk=1:length(eigDirTmp(:,1))
            distTMP(kk,1)=distance(ImTPM,DmTPM,eigDirTmp(kk,2),eigDirTmp(kk,1));
        end
        indTMP=find(distTMP==min(distTMP));
        meanW.meanW(1:2)=eigDirTmp(indTMP,1:2);
    end
%     [1 distance(lambdaMIN,phiMIN,ImTPM,DmTPM) distance(lambdaMINa,phiMINa,ImTPM,DmTPM)]
else
    
    omega=(phiINT-phiMIN)/abs(phiINT-phiMIN)*acosd(sind(lambdaINT)/cosd(lambdaMIN));
    if isnan(omega)==1 omega=zeros(length(omega),1); end
    if isreal(omega)==0 omega=zeros(length(omega),1); end
    
    meanW.eigDir(1,:)=[phiMIN,lambdaMIN,rMIN,phiINT,lambdaINT,rINT,phiMAX,lambdaMAX,rMAX];
    meanW.eigDir(2,:)=[phiMINa,lambdaMINa,rMIN,phiINTa,lambdaINTa,rINT,phiMAXa,lambdaMAXa,rMAX];
    meanW.meanW=[phiMIN,lambdaMIN,e95a,e95b,omega,Kx,Ky,nptsM];
    
    if npts==1
        eigDirTmp=[phiMINa,lambdaMINa;phiINT,lambdaINT;phiMAX,lambdaMAX;...
            phiMIN,lambdaMIN;phiINTa,lambdaINTa;phiMAXa,lambdaMAXa];
        for kk=1:length(eigDirTmp(:,1))
            distTMP(kk,1)=distance(ImTPM,DmTPM,eigDirTmp(kk,2),eigDirTmp(kk,1));
        end
        indTMP=find(distTMP==min(distTMP));
        meanW.meanW(1:2)=eigDirTmp(indTMP,1:2);
    end
%     [2 distance(lambdaMIN,phiMIN,ImTPM,DmTPM) distance(lambdaMINa,phiMINa,ImTPM,DmTPM)]
end

for i=1:length(alpha95TPM)
    if isreal(alpha95TPM(i))==0 alpha95TPM(i)=0; end
end
meanW.fishM=[DmTPM,ImTPM,alpha95TPM,kappaTPM,NTPM];

% meanW=[lonM,latM,e95a,e95b,omega,Kx,Ky,N];
% ensure that when there is only one paleopole, e95a=e95b=A95
if npts==1
    meanW.meanW(3:4)=meanW.fishM(3)*ones(1,2);
    meanW.meanW(5)=0;
end

%--------------------------------------------------------------------------
% fix the err orientation issue related to the polarity of input paleopoles
if polesRflag==0  meanW.meanW(5)=-meanW.meanW(5); end
%--------------------------------------------------------------------------

% calculate fisher parameters
function [Dm,Im,alpha95,kappa,N]=FishM(D,I,a95)

D=D./(180/pi);
I=I./(180/pi);
[x,y,z]=sph2cart(D,I,1);
N=length(x);

if N>1
    R2=(sum(x))^2+(sum(y))^2+(sum(z))^2;
    R=sqrt(R2);
    
    m1=(sum(x))/R;
    m2=(sum(y))/R;
    m3=(sum(z))/R;
    
    % Fisherian Parameter: kappa(K); alpha95(A95)
    kappa=(N-1)/(N-R);
    alpha95=acosd(1-(N-R)/R*(((1/0.05)^(1/(N-1)))-1));
    
    % Convert back to (Im,Dm)
    [Dm,Im]=cart2sph(m1,m2,m3);
    Dm=Dm.*(180/pi);
    Im=Im.*(180/pi);
elseif N==1
    Dm=D.*(180/pi); Im=I.*(180/pi); alpha95=a95; kappa=0; N=1;
elseif N==0
    Dm=NaN; Im=NaN; alpha95=NaN; kappa=NaN;N=0;
end





