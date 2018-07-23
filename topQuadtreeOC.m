%topQuadtreeOC(512,256,0.4, 2,1,200)
function topQuadtreeOC(nelx,nely,volfrac, bSubdivision,bBlackWhite,ItMax)
% bSubdivision: 1 - subdivision
%               2 - balanced subdivision
% bBlackWhite:  1 - on
%               0 - off
clc;
close all;
MMA = 0;
 
folder = sprintf('images-%8.4f',volfrac);
mkdir(folder);
cList = zeros(ItMax,1);
vList = zeros(ItMax,1);
sList = zeros(ItMax,1);
 
penal = 3;
pNorm = 16;
 
nLevels = 4;
% Cell dimension
nc0 = 2^(nLevels+2);   % Coarset cell, e.g., nc0*nc0 = 32*32
nc = zeros(nLevels,1); 
for level = 1:nLevels 
    nc(level) = nc0/(2^(level-1));
end
 
% Elements within cell
nx0 = nelx/nc0;
ny0 = nely/nc0;
nx = zeros(nLevels,1);
ny = zeros(nLevels,1); 
for level = 1:nLevels
    nx(level) = nx0 * 2^(level-1);
    ny(level) = ny0 * 2^(level-1); 
end
 
h_min = 0;%0.00001;
 
% Design variables
h0 = ones(nx0*ny0,1);
h = {zeros(nx(1)*ny(1),1)};
for level = 2:nLevels
    h{level} = zeros(nx(level)*ny(level),1);
end
 
% A copy of design variables, for caculating the initial volume
h1 = {ones(nx(1)*ny(1),1)};
for level = 2:nLevels
    h1{level} = ones(nx(level)*ny(level),1);
end
 
% Design variables projected by the subdivision rule
hp = h; 
% Design variables converted to black-white design
hbw = h;
 
%%
num_x = nelx*nely;
[M0x] = buildSparseM0x(nx,ny,nc0,nelx,nely,num_x);
[Mx] = buildSparseMx(nLevels,nc,nx,ny,nelx,nely,num_x);  
 
[xFix] = designToX(M0x,h0,Mx,nLevels,h);
[xFix1] = designToX(M0x,h0,Mx,nLevels,h1);
initFrac = (volfrac*nelx*nely - sum(xFix(:))) / (sum(xFix1(:)) - sum(xFix(:)));
for level = 1:nLevels
    h{level} = initFrac * h1{level};
end
 
[hlinear] = linearize(h,nLevels);
[hplinear] = linearize(hp,nLevels);
 
numDesign = 0; % Number of total design variables
for level = 1:nLevels
    numDesign = numDesign + size(h{level},1);
end
 
[hde,hdeCo] = initializeDependency(nx,ny,nLevels);
 
[hdeEx,hdeExiLevels,hdeExNum,hdeExCo] = initializeExtendedDependency(nx,ny,nLevels,hde);
 
%% subdivisionProjection
if bSubdivision == 1
    [hp,hdeCo] = subdivisionProjectionPNorm(h,hp,h0,hde,hdeCo,pNorm,nLevels,nx,ny);
elseif bSubdivision == 2
    [hp,hdeExCo] = restrictedSubdivisionProjectionPNorm(h,hp,h0,hdeEx,hdeExCo,hdeExiLevels,hdeExNum,pNorm,nLevels,nx,ny);
else
    hp=h;
end 
 
beta = 1;       % beta continuation
eta = 0.5;      % projection threshold, fixed at 0.5
 
%% Adapted from 88 lines
E0 = 1;
Emin = 1e-9;
nu = 0.3;
%% PREPARE FINITE ELEMENT ANALYSIS
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
F = sparse(2*(nely+1)*(nelx+1)-nely,1,-1,2*(nely+1)*(nelx+1),1); %-nely
U = zeros(2*(nely+1)*(nelx+1),1);
fixeddofs = union([1:1:2*(nely+1)],[1]);
alldofs = [1:2*(nely+1)*(nelx+1)];
freedofs = setdiff(alldofs,fixeddofs);
 
%% INITIALIZE ITERATION
[x] = designToX(M0x,h0,Mx,nLevels,h);
x = reshape(x,nelx,nely)';
xPhys = x;
loop = 0;
loopbeta = 0;
change = 1;
hlinearold1 = hlinear;
hlinearold2 = hlinear;
low = 0;
upp = 0;
%% START ITERATION     
while change > 0.0001 && loop < ItMax
    loop = loop + 1;
    loopbeta = loopbeta + 1;
    %% FE-ANALYSIS
    sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
    K = sparse(iK,jK,sK); K = (K+K')/2;
    U(freedofs) = K(freedofs,freedofs)\F(freedofs);
    %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
    c = sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce));
    dcx = -penal*(E0-Emin)*xPhys.^(penal-1).*ce;
    dvx = ones(nely,nelx)/(nelx*nely*volfrac);
 
    cList(loop) = c;
    vList(loop) = mean(x(:));
     
    dch = {zeros(1,size(h{1},1))};
    dvh = {zeros(1,size(h{1},1))};
    for level = 2:nLevels
        dch{level} = zeros(1,size(h{level},1));
        dvh{level} = zeros(1,size(h{level},1));
    end
    dcxlinear = reshape(dcx',1,nely*nelx);
    dvxlinear = reshape(dvx',1,nely*nelx);
    for level = 1:nLevels
        dch{level} = dcxlinear*Mx{level};
        dvh{level} = dvxlinear*Mx{level};
    end
    dc = dch{1};
    dv = dvh{1};
    for level = 2:nLevels
        dc = [dc dch{level}];
        dv = [dv dvh{level}];
    end
 
    if bSubdivision > 0
        dhpdh = zeros(numDesign,numDesign);
        offset = zeros(nLevels,1);
        offset(1,1) = 0;
        for level = 2:nLevels
            offset(level,1) = offset(level-1,1) + size(h{level-1},1);
        end
         
        if bSubdivision == 1
            for level = 1:nLevels
                for i=1:nx(level)
                    for j=1:ny(level)
                        it = (j-1)*nx(level)+(i-1)+1;
                        for ii = 1:(level-1)
                            ic = hde{level}(it,ii);
                            dhpdh(offset(level,1)+it,offset(ii,1)+ic) = hdeCo{level}(it,ii);
                        end
                        dhpdh(offset(level,1)+it,offset(level,1)+it) = hdeCo{level}(it,level);
                    end
                end
            end
        elseif bSubdivision == 2        
            for level = 1:nLevels
                for i=1:nx(level)
                    for j=1:ny(level)
                        it = (j-1)*nx(level)+(i-1)+1;
                        numEx = hdeExNum{level}(it,1);
                        for ii = 1:numEx
                            ic = hdeEx{level}(it,ii);
                            iLevel = hdeExiLevels{level}(it,ii);
                            dhpdh(offset(level,1)+it,offset(iLevel,1)+ic) = hdeExCo{level}(it,ii);
                        end
                        dhpdh(offset(level,1)+it,offset(level,1)+it) = hdeExCo{level}(it,numEx+1);
                    end
                end
            end
        end
 
        dc = dc*dhpdh;
        dv = dv*dhpdh;
    end
     
    if bBlackWhite
        dx = beta * (1-tanh(beta*(hlinear-eta)).*tanh(beta*(hlinear-eta))) / (tanh(beta*eta) + tanh(beta*(1-eta)));
        dc = dc.*dx';
        dv = dv.*dx';
    end
     
    if MMA == 1
        %% MMA
        cScale = 1e-2;
        f0val = c*cScale;
        if (f0val > 100 || f0val < 1)% && maxIter > 1
            cScale = 10/c;
            f0val = c*cScale;
        end
 
        df0dx = reshape(dc,numDesign,1)*cScale;
        df0dx2 = 0*df0dx;
        dfdx = reshape(dv,1,numDesign);
        dfdx2 = 0*dfdx;
 
        iter = loopbeta;
        xval = hlinear;
        move = 0.01;
        xmin=max(h_min,xval-move);
        xmax=min(1,xval+move);
 
        fval = sum(sum(xPhys)) / (nelx*nely*volfrac) - 1;
         
        m = 1;
        n = numDesign; 
        [xmma,ymma,zmma,lam,xsi,eta_,mu,zet,s,low,upp] = ...
            mmasub(m,n,iter,xval,xmin,xmax,hlinearold1,hlinearold2,...
            f0val,df0dx,df0dx2,fval,dfdx,dfdx2,low,upp,1,0,1000,0);
        hlinearold2 = hlinearold1;
        hlinearold1 = xval;
         
        hnewlinear = xmma;
        offset = 0;
        for level = 1:nLevels
            if level == 1
                hnew = {hnewlinear((offset+1):(offset+size(h{level},1)), 1)};
            else
                hnew{level} = hnewlinear((offset+1):(offset+size(h{level},1)), 1);
            end
            offset = offset + size(h{level},1);
        end
         
        if bBlackWhite
            for level = 1:nLevels
                hbw{level} = (tanh(beta*eta) + tanh(beta*(hnew{level}-eta))) / (tanh(beta*eta) + tanh(beta*(1-eta)));
            end
        else
            hbw = hnew;
        end
 
        if bSubdivision > 0
            if bSubdivision == 1
                [hp,hdeCo] = subdivisionProjectionPNorm(hbw,hp,h0,hde,hdeCo,pNorm,nLevels,nx,ny);
            elseif bSubdivision == 2
                [hp,hdeExCo] = restrictedSubdivisionProjectionPNorm(hbw,hp,h0,hdeEx,hdeExCo,hdeExiLevels,hdeExNum,pNorm,nLevels,nx,ny);
            end
        else
            hp = hbw;
        end
        [xnew] = designToX(M0x,h0,Mx,nLevels,hp);
    else
    %% OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES
        l1 = 0; l2 = 1e9; move = 0.01;
        while (l2-l1)/(l1+l2) > 1e-3%0.5*1e-3
            lmid = 0.5*(l2+l1);
            hnewlinear = max(h_min,max(hlinear-move,min(1,min(hlinear+move,hlinear.*(-dc'./dv'/lmid).^(1/2)))));
             
            offset = 0;
            for level = 1:nLevels
                if level == 1
                    hnew = {hnewlinear((offset+1):(offset+size(h{level},1)), 1)};
                else
                    hnew{level} = hnewlinear((offset+1):(offset+size(h{level},1)), 1);
                end
                offset = offset + size(h{level},1);
            end
                 
            if bBlackWhite
                for level = 1:nLevels
                    hbw{level} = (tanh(beta*eta) + tanh(beta*(hnew{level}-eta))) / (tanh(beta*eta) + tanh(beta*(1-eta)));
                end
            else
                hbw = hnew;
            end
 
            if bSubdivision > 0                
                if bSubdivision == 1
                    [hp,hdeCo] = subdivisionProjectionPNorm(hbw,hp,h0,hde,hdeCo,pNorm,nLevels,nx,ny);
                elseif bSubdivision == 2
                    [hp,hdeExCo] = restrictedSubdivisionProjectionPNorm(hbw,hp,h0,hdeEx,hdeExCo,hdeExiLevels,hdeExNum,pNorm,nLevels,nx,ny);
                end
            else
                hp = h;
            end
            [xnew] = designToX(M0x,h0,Mx,nLevels,hp);
 
            if sum(xnew(:)) > volfrac*nelx*nely, l1 = lmid; else l2 = lmid; end
        end
    end
     
    change = max(abs(xnew(:)-x(:)));
    x = xnew;
    xPhys = reshape(x,nelx,nely)';
    h = hnew;
    [hlinear] = linearize(h,nLevels);
    [hplinear] = linearize(hp,nLevels);
 
    %% Beta-continuation for Heaviside projection
    if bBlackWhite
        if beta < 20 && (loopbeta >= 60 || change <= 0.001)
            beta = 2*beta;
            loopbeta = 0;
            change = 1;
            fprintf('Parameter beta increased to %g.\n',beta);
        end
    end
     
    tmp = (x.*(1-x)*4);
    sharp = sum(tmp(:))/(nelx*nely);
    sList(loop) = sharp;
 
    %% PRINT RESULTS
    fprintf(' It.:%5i Obj.:%11.4f Vol.:%11.4f ch.:%7.3f, Sharpness: %7.3f\n',loop,c, ...
        mean(x(:)),change,sharp);
    %% PLOT DENSITIES
    figure(1);
    colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;
     
    filename1 = sprintf('images-%8.4f\\rho-It%i.png',volfrac,loop);
    saveas(1,filename1,'png');
    if rem(loop,5) == 0
        filename1 = sprintf('images-%8.4f\\rho-It%i.fig',volfrac,loop);
        savefig(1, filename1,'compact');
    end
end
 
%% 
function [hde,hdeCo] = initializeDependency(nx,ny,nLevels)
hde = {zeros(nx(1)*ny(1),1)}; % Design variables' dependence
for level = 2:nLevels % From second level upwards, the dependence to zeroth is neglected.
    hde{level} = zeros(nx(level)*ny(level),level-1);
end
% The first level doesn't depend on the zeroth. The following is added for
% consistency
for i=1:nx(1)
    for j=1:ny(1)
        it = (j-1)*nx(1)+(i-1)+1;
        ic = it;
        hde{1}(it) = ic;
    end
end
 
hdeCo = {zeros(nx(1)*ny(1),1)}; % Design variables' dependence, coefficient of derivative
for level = 2:nLevels % From second level upwards, the dependence to zeroth is neglected.
    hdeCo{level} = zeros(nx(level)*ny(level),level);
end
 
for level = 2:nLevels
   for i=1:nx(level)
        for j=1:ny(level)
            it = (j-1)*nx(level)+(i-1)+1;
            ic = (floor((j+1)/2)-1)*nx(level-1) + (floor((i+1)/2)-1) + 1;
            hde{level}(it,level-1) = ic;
            if level > 2
                for ii = 1:(level-2)
                    hde{level}(it, ii) = hde{level-1}(ic,ii);
                end
            end
        end
    end
end
 
function [hdeEx,hdeExiLevels,hdeExNum,hdeExCo] = initializeExtendedDependency(nx,ny,nLevels,hde)
hdeEx = {zeros(nx(1)*ny(1),1)};
for level = 2:nLevels
    hdeEx{level} = zeros(nx(level)*ny(level),3^(level-1));
end
 
hdeExiLevels = {zeros(nx(1)*ny(1),1)};
for level = 2:nLevels
    hdeExiLevels{level} = zeros(nx(level)*ny(level),3^(level-1));
end
 
hdeExNum = {zeros(nx(1)*ny(1),1)};
for level = 2:nLevels
    hdeExNum{level} = zeros(nx(level)*ny(level),1);
end
for level = 2:nLevels
   for i=1:nx(level)
        for j=1:ny(level)
            it = (j-1)*nx(level)+(i-1)+1;
            ic = (floor((j+1)/2)-1)*nx(level-1) + (floor((i+1)/2)-1) + 1;
%             hde{level}(it,level-1) = ic;
%             if level > 2
%                 for ii = 1:(level-2)
%                     hde{level}(it, ii) = hde{level-1}(ic,ii);
%                 end
%             end
            if rem(i,2) == 0
                ix = i+1;
            else
                ix = i-1;
            end
            if rem(j,2) == 0
                jy = j+1;
            else
                jy = j-1;
            end
 
            icx = 0;
            icy = 0;
             
            if ix >= 1 && ix <= nx(level)
                itx = (j-1)*nx(level)+(ix-1)+1;
                icx = hde{level}(itx,level-1);
            end
            if jy >= 1 && jy <= ny(level)
                ity = (jy-1)*nx(level)+(i-1)+1;
                icy = hde{level}(ity,level-1);
            end
 
            hdeExNum{level}(it,1) = hdeExNum{level}(it,1) + 1;
            hdeEx{level}(it,hdeExNum{level}(it,1)) = ic;
            hdeExiLevels{level}(it,hdeExNum{level}(it,1)) = level-1;
            if icx ~= 0
                hdeExNum{level}(it,1) = hdeExNum{level}(it,1) + 1;
                hdeEx{level}(it,hdeExNum{level}(it,1)) = icx;
                hdeExiLevels{level}(it,hdeExNum{level}(it,1)) = level-1;
            end
            if icy ~= 0
                hdeExNum{level}(it,1) = hdeExNum{level}(it,1) + 1;
                hdeEx{level}(it,hdeExNum{level}(it,1)) = icy;
                hdeExiLevels{level}(it,hdeExNum{level}(it,1)) = level-1;
            end
             
            if level > 2
                for ii = 1:(level-2)
                    for iii = 1:hdeExNum{level-1}(ic,1)
                        bExist = 0;
                        for jj = 1:hdeExNum{level}(it,1)
                            if hdeEx{level}(it,jj) == hdeEx{level-1}(ic,iii)
                                if hdeExiLevels{level}(it,jj) == hdeExiLevels{level-1}(ic,iii)
                                    bExist = 1;
                                end
                            end
                        end
                        if bExist == 0
                            hdeExNum{level}(it,1) = hdeExNum{level}(it,1) + 1;
                            hdeEx{level}(it,hdeExNum{level}(it,1)) = hdeEx{level-1}(ic,iii);
                            hdeExiLevels{level}(it,hdeExNum{level}(it,1)) = hdeExiLevels{level-1}(ic,iii);
                        end
                    end
                     
                    if icx ~= 0
                        for iii = 1:hdeExNum{level-1}(icx,1)
                            bExist = 0;
                            for jj = 1:hdeExNum{level}(it,1)
                                if hdeEx{level}(it,jj) == hdeEx{level-1}(icx,iii)
                                    if hdeExiLevels{level}(it,jj) == hdeExiLevels{level-1}(icx,iii)
                                        bExist = 1;
                                    end
                                end
                            end
                            if bExist == 0
                                hdeExNum{level}(it,1) = hdeExNum{level}(it,1) + 1;
                                hdeEx{level}(it,hdeExNum{level}(it,1)) = hdeEx{level-1}(icx,iii);
                                hdeExiLevels{level}(it,hdeExNum{level}(it,1)) = hdeExiLevels{level-1}(icx,iii);
                            end
                        end
                    end
 
                    if icy ~= 0
                        for iii = 1:hdeExNum{level-1}(icy,1)
                            bExist = 0;
                            for jj = 1:hdeExNum{level}(it,1)
                                if hdeEx{level}(it,jj) == hdeEx{level-1}(icy,iii)
                                    if hdeExiLevels{level}(it,jj) == hdeExiLevels{level-1}(icy,iii)
                                        bExist = 1;
                                    end
                                end
                            end
                            if bExist == 0
                                hdeExNum{level}(it,1) = hdeExNum{level}(it,1) + 1;
                                hdeEx{level}(it,hdeExNum{level}(it,1)) = hdeEx{level-1}(icy,iii);
                                hdeExiLevels{level}(it,hdeExNum{level}(it,1)) = hdeExiLevels{level-1}(icy,iii);
                            end
                        end
                    end
                end
            end
        end
   end
end 
 
hdeExCo = {zeros(nx(1)*ny(1),1)}; % Design variables' dependence, coefficient of derivative
for level = 2:nLevels % From second level upwards, the dependence to zeroth is neglected.
    hdeExCo{level} = zeros(nx(level)*ny(level),3^(level-1));
end
 
function [x] = designToX(M0x,h0,Mx,nLevels,h)
x = M0x*h0;
for level = 1:nLevels
    x = x + Mx{level}*h{level};
end
 
function [hp,hdeCo] = subdivisionProjectionPNorm(h,hp,h0,hde,hdeCo,pNorm,nLevels,nx,ny)
for i=1:nx(1)
    for j=1:ny(1)
        it = (j-1)*nx(1)+(i-1)+1;
        ic = it;
        hp{1}(it) = min(h{1}(it), h0(ic));
        hdeCo{1}(it,1) = 1;
    end
end
 
for level = 2:nLevels
    de = zeros(level,1);
    for i=1:nx(level)
        for j=1:ny(level)
            it = (j-1)*nx(level)+(i-1)+1;
            de(level,1) = h{level}(it);
            for ii = 1:(level-1)
                ic = hde{level}(it, ii);
                de(ii,1) = h{ii}(ic);
            end
             
%             hp{level}(it) = 0;
%             for ii = 1:level
%                 hp{level}(it) = hp{level}(it) + (1-de(ii,1))^pNorm;
%             end
%             A = max(hp{level}(it) / level, 10^(-8));
%             hp{level}(it) = 1 - (hp{level}(it) / level)^(1/pNorm);
% 
%             for ii = 1:level
%                 hdeCo{level}(it,ii) = A^((1/pNorm)-1)*(1-de(ii,1)^(pNorm-1))/level;
%             end
            hp{level}(it) = 0;
            for ii = 1:level
                hp{level}(it) = hp{level}(it) + (de(ii,1))^(-pNorm);
            end
            A = hp{level}(it) / level;
            hp{level}(it) = A^(-1/pNorm);
 
            for ii = 1:level
                hdeCo{level}(it,ii) = A^((-1/pNorm)-1)*(de(ii,1)^(-pNorm-1))/level;
            end
        end
    end
end
 
function [hp,hdeExCo] = restrictedSubdivisionProjectionPNorm(h,hp,h0,hdeEx,hdeExCo,hdeExiLevels,hdeExNum,pNorm,nLevels,nx,ny)
for i=1:nx(1)
    for j=1:ny(1)
        it = (j-1)*nx(1)+(i-1)+1;
        ic = it;
        hp{1}(it) = min(h{1}(it), h0(ic));
        hdeExCo{1}(it,1) = 1;
    end
end
 
for level = 2:nLevels
    de = zeros(3^(level-1),1);
    for i=1:nx(level)
        for j=1:ny(level)
            it = (j-1)*nx(level)+(i-1)+1;
            numEx = hdeExNum{level}(it,1);
            de(numEx+1) = h{level}(it);
            for ii = 1:numEx
                ic = hdeEx{level}(it,ii);
                iLevel = hdeExiLevels{level}(it,ii);
                de(ii) = h{iLevel}(ic);
            end
 
%             hp{level}(it) = 0;
%             for ii = 1:(numEx+1)
%                 hp{level}(it) = hp{level}(it) + (1-de(ii))^pNorm;
%             end
%             
%             A = hp{level}(it) / (numEx+1);
%             hp{level}(it) = 1 - (hp{level}(it) / (numEx+1))^(1/pNorm);
%             
%             for ii = 1:(numEx+1)
%                 hdeExCo{level}(it,ii) = A^((1/pNorm)-1)*(1-de(ii)^(pNorm-1))/(numEx+1);
%             end
            hp{level}(it) = 0;
            for ii = 1:(numEx+1)
                hp{level}(it) = hp{level}(it) + (de(ii))^(-pNorm);
            end
             
            A = hp{level}(it) / (numEx+1);
            hp{level}(it) = (hp{level}(it) / (numEx+1))^(-1/pNorm);
             
            for ii = 1:(numEx+1)
                hdeExCo{level}(it,ii) = A^((-1/pNorm)-1)*(de(ii)^(-pNorm-1))/(numEx+1);
            end
        end
    end
end
 
function [hlinear] = linearize(h,nLevels)
hlinear = h{1};
for level = 2:nLevels
    hlinear = [hlinear; h{level}];
end
 
function [M0x] = buildSparseM0x(nx,ny,nc0,nelx,nely,num)
M0x = sparse(num,nx(1)*ny(1));
dim = nx(1)*ny(1)*(4*nc0-4);
iter = 1;
iM = zeros(dim,1);
jM = zeros(dim,1);
vM = ones(dim,1);
for i=1:nx(1)
    ii = (i-1)*nc0;
    for j=1:ny(1)
        jj = (j-1)*nc0;
        ic = (j-1)*nx(1) + (i-1)+1;   % index of the cell
         
        for k=1:nc0
            iii = ii+k;
            jjj = jj+1;
            it = (jjj-1)*nelx + (iii-1) + 1;    % index of the element
            %M0x(it,ic) = 1;
            iM(iter) = it;
            jM(iter) = ic;
            iter = iter + 1;
        end
         
        for k=1:nc0
            iii = ii+k;
            jjj = jj+nc0;
            it = (jjj-1)*nelx + (iii-1) + 1;    % index of the element
            %M0x(it,ic) = 1;
            iM(iter) = it;
            jM(iter) = ic;
            iter = iter + 1;
        end
         
        for k=2:(nc0-1)
            iii = ii+1;
            jjj = jj+k;
            it = (jjj-1)*nelx + (iii-1) + 1;    % index of the element
            %M0x(it,ic) = 1;
            iM(iter) = it;
            jM(iter) = ic;
            iter = iter + 1;
        end
         
        for k=2:(nc0-1)
            iii = ii+nc0;
            jjj = jj+k;
            it = (jjj-1)*nelx + (iii-1) + 1;    % index of the element
            %M0x(it,ic) = 1;
            iM(iter) = it;
            jM(iter) = ic;
            iter = iter + 1;
        end
    end
end
M0x = sparse(iM',jM',vM',num,nx(1)*ny(1));
 
function [Mx] = buildSparseMx(nLevels,nc,nx,ny,nelx,nely,num)
Mx = {sparse(num,nx(1)*ny(1))};
for level = 2:nLevels
    Mx{level} = sparse(num,nx(level)*ny(level));
end
for level = 1:nLevels
    nc_l = nc(level);
    dim = nx(level)*ny(level)*(4*(nc_l-2)-4);
    iter = 1;
    iM = zeros(dim,1);
    jM = zeros(dim,1);
    vM = ones(dim,1);
    for i=1:nx(level)
        ii = (i-1)*nc_l;
        for j=1:ny(level)
            jj = (j-1)*nc_l;
            ic = (j-1)*nx(level) + (i-1)+1;   % index of the corse cell
 
            for k=2:nc_l-1
                iii = ii+k;
                jjj = jj+nc_l/2;
                it = (jjj-1)*nelx + (iii-1) + 1;    % index of the element
                %Mx{level}(it,ic) = 1;
                iM(iter) = it;
                jM(iter) = ic;
                iter = iter + 1;
 
                jjj = jj+nc_l/2+1;
                it = (jjj-1)*nelx + (iii-1) + 1;    % index of the element
                %Mx{level}(it,ic) = 1;
                iM(iter) = it;
                jM(iter) = ic;
                iter = iter + 1;
            end
 
            for k=2:(nc_l/2-1)
                jjj = jj+k;
                iii = ii+nc_l/2;
                it = (jjj-1)*nelx + (iii-1) + 1;    % index of the element
                %Mx{level}(it,ic) = 1;
                iM(iter) = it;
                jM(iter) = ic;
                iter = iter + 1;
 
                iii = ii+nc_l/2+1;
                it = (jjj-1)*nelx + (iii-1) + 1;    % index of the element
                %Mx{level}(it,ic) = 1;
                iM(iter) = it;
                jM(iter) = ic;
                iter = iter + 1;
            end
             
            for k=(nc_l/2+2):(nc_l-1)
                jjj = jj+k;
                iii = ii+nc_l/2;
                it = (jjj-1)*nelx + (iii-1) + 1;    % index of the element
                %Mx{level}(it,ic) = 1;
                iM(iter) = it;
                jM(iter) = ic;
                iter = iter + 1;
 
                iii = ii+nc_l/2+1;
                it = (jjj-1)*nelx + (iii-1) + 1;    % index of the element
                %Mx{level}(it,ic) = 1;
                iM(iter) = it;
                jM(iter) = ic;
                iter = iter + 1;
            end
        end
    end
    Mx{level} = sparse(iM',jM',vM',num,nx(level)*ny(level));
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was written by Jun Wu, Department of Design Engineering,%
% Delft University of Technology, The Netherlands.                         %
% Please send your comments to: j.wu-1@@tudelft.nl                         %
%                                                                          %
% The code is intended for educational purposes and theoretical details    %
% are discussed in the paper                                               %
% "Continuous Optimization of Adaptive Quadtree Structures"                %
% Jun Wu, Computer-Aided Design, vol.102, pp.72-82, 2018                   %
%                                                                          %
% The code as well as a preprint version of the paper can be               %
% downloaded from the web-site: www.jun-wu.net                             %
%                                                                          %
% The code was developed based on the 88-line topopt code                  %
% by E. Andreassen et al (2010), Struct Multidisc Optim.                   %
%                                                                          %
% Disclaimer:                                                              %
% The authors reserves all rights but do not guaranty that the code is     %
% free from errors. Furthermore, we shall not be liable in any event       %
% caused by the use of the program.                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
