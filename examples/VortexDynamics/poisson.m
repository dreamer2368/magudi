function poisson(Ng,R,Ma,gamma)

%%

L = 3*R; dx = 2*L/(Ng-1);
xg = linspace(-L,L,Ng); yg=linspace(-L,L,Ng);
[Yg,Xg] = meshgrid(xg,yg);

%%

B1 = repmat([-1/60 3/20 -3/4 3/4 -3/20 1/60],Ng,1);
B2 = repmat([1/90 -3/20 3/2 -49/18 3/2 -3/20 1/90],Ng,1);

D1 = spdiags(B1,[-3 -2 -1 1 2 3],Ng,Ng);
D1(1,1:4) = [-3/2, 2, -1/2, 0];
D1(2,1:5) = [-1/2, 0, 1/2, 0, 0];
D1(3,1:6) = [1/12, -2/3, 0, 2/3, -1/12, 0];
D1(end-2,end-5:end) = [0, 1/12, -2/3, 0, 2/3, -1/12];
D1(end-1,end-4:end) = [0, 0, -1/2, 0, 1/2];
D1(end,end-3:end) = [0, 1/2, -2, 3/2];
D2 = spdiags(B2,[-3 -2 -1 0 1 2 3],Ng,Ng);
I1 = speye(Ng);

Dx = kron(I1,D1)/dx; Dy = kron(D1,I1)/dx;
Dx2 = kron(I1,D2)/dx/dx; Dy2 = kron(D2,I1)/dx/dx;
K = Dx2+Dy2;

%%

loci = importdata('loci.txt');
N = size(loci,1);
    
%%

radii = zeros(Ng,Ng);
vt = radii; vx = vt; vy = vt;
% dvtdr = vt; vtOverR = vt;
S = zeros(Ng,Ng); % S0 = S;
for k=1:N
    radii = sqrt( (Xg-loci(k,1)).^2 + (Yg-loci(k,2)).^2 );
    
    vt = 1.428./radii.*( 1-exp(-1.25*radii.^2) );
%     dvtdr(:,:,k) = 3.57*exp(-1.25*radii(:,:,k).^2)  ...
%                     - 1.428./(radii(:,:,k).^2).*( 1-exp(-1.25*radii(:,:,k).^2) );
%     vtOverR(:,:,k) = vt(:,:,k)./radii(:,:,k);
%     S0 = S0 - 2*dvtdr.*vtOverR;
    
    vx = vx - vt.*(Yg-loci(k,2))./radii;
    vy = vy + vt.*(Xg-loci(k,1))./radii;
end

S = reshape( (Dx*vx(:)).^2 + 2*(Dx*vy(:)).*(Dy*vx(:)) + (Dy*vy(:)).^2, Ng, Ng);

% figure(1)
% surf(Xg,Yg,vx,'edgecolor','none');
% colorbar;
% view([0 0 1]);
% 
% figure(2)
% surf(Xg,Yg,vy,'edgecolor','none');
% colorbar
% view([0 1 0]);
% 
% figure(3)
% surf(Xg,Yg,reshape(Dx*vx(:),Ng,Ng),'edgecolor','none');
% colorbar;
% view([0 0 1]);
% 
% figure(4)
% surf(Xg,Yg,reshape(Dy*vx(:),Ng,Ng),'edgecolor','none');
% colorbar;
% view([0 0 1]);

% figure(5)
% surf(Xg,Yg,S0,'edgecolor','none');
% colorbar;
% view([0 0 1]);

% figure(6)
% surf(Xg,Yg,log10(abs(S0-S)),'edgecolor','none');
% colorbar;
% view([0 0 1]);

%%

RHS = -(gamma-1)/gamma*Ma^2*S;
pOverRho = reshape( K\RHS(:), Ng,Ng ) + 1/gamma;

p0 = ( gamma^(1/gamma)*pOverRho ).^(gamma/(gamma-1));
rho0 = (gamma*p0).^(1/gamma);

% figure(7)
% surf(Xg,Yg,pOverRho,'edgecolor','none');
% colorbar;
% view([0 0 1]);
% 
% figure(8)
% surf(Xg,Yg,p0,'edgecolor','none');
% colorbar;
% view([0 0 1]);
% 
% figure(9)
% surf(Xg,Yg,rho0,'edgecolor','none');
% colorbar;
% view([0 0 1]);

%%

fileID = fopen('p0.bin','w');
fwrite(fileID,p0,'double');
fclose(fileID);

fileID = fopen('rho0.bin','w');
fwrite(fileID,rho0,'double');
fclose(fileID);

% fileID = fopen('rho0.bin');
% rho1 = fread(fileID,'double');
% fclose(fileID);

end