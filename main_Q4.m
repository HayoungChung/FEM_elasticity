% test script for Q4 linear FEM
clc
clear all

DEBUGFLAG__ = 0;
%% Material model (lin Elasticity)
E = 1000; v = 0.3; 
matD_voigt = E/(1-v^2)*[1 v 0;v 1 0; 0 0 (1-v)/2];

%% mesh initialization
lxy = [1 1];
exy = [10 10];
nxy = exy + 1;

totArea = prod(lxy);
nELEM = prod(exy);
nNODE = prod(nxy);

x_list = linspace(0,lxy(1),nxy(1));
y_list = linspace(0,lxy(2),nxy(2));

[x, y] = meshgrid(x_list, y_list); 
NODE = zeros(nNODE,2);
NODE(:,1) = x(:);
NODE(:,2) = y(:);

ELEM = zeros(nELEM,4);
for xx = 1:exy(1)
    for yy = 1:exy(2)
        eid = exy(2)*(xx-1)+ yy;
        ELEM(eid,1) = nxy(2)*(xx-1) + yy;
        ELEM(eid,2) = nxy(2)*(xx) + yy;
        ELEM(eid,3) = nxy(2)*(xx) + yy + 1;
        ELEM(eid,4) = nxy(2)*(xx-1) + yy +1 ;
    end
end

% plot test
if DEBUGFLAG__
    figure(1); clf; 
    hold on 
    for ee = 1:nELEM
        elem_id = ELEM(ee,:);
        plot(NODE(elem_id([1,2,3,4,1]),1),NODE(elem_id([1,2,3,4,1]),2),'r-')
    end
end

%% K matrix generation
nDOF = nNODE * 2; % dpn = 2;
dpe = 8;
K_dense = zeros(nDOF, nDOF);
K_sparse = zeros(nELEM*dpe*dpe,3);

ri = [-1 1 1 -1]/sqrt(3);
si = [-1 -1 1 1]/sqrt(3);
wi = [1 1 1 1];

for ee = 1:nELEM
    elem_id = ELEM(ee,:);
    X = NODE(elem_id,:); % [4x2]
    
    Ke = zeros(dpe,dpe);
    for gg = 1:length(wi)
        r = ri(gg); s = si(gg); w = wi(gg);
        % shape function and its derivative
        Ni   = 0.25*[(1-r)*(1-s) (1+r)*(1-s) (1+r)*(1+s) (1-r)*(1+s)];
        Ni_r = 0.25*[-1*(1-s) +1*(1-s) +1*(1+s) -1*(1+s)];
        Ni_s = 0.25*[-1*(1-r) -1*(1+r) +1*(1+r) +1*(1-r)];
        
        Ni_rs = [Ni_r;Ni_s]; % [2x4]
        matJ = Ni_rs*X;
        % matJ = [Ni_r*X(:,1) Ni_r*X(:,2);Ni_s*X(:,1) Ni_s*X(:,2)];
        detJ = det(matJ);
        
        Ni_XY = inv(matJ)*Ni_rs;
        Ni_X = Ni_XY(1,:); Ni_Y = Ni_XY(2,:);
        
        matB = zeros(3,8);
        matB(1,1:2:end) = Ni_X;
        matB(2,2:2:end) = Ni_Y;
        matB(3,1:2:end) = Ni_Y;
        matB(3,2:2:end) = Ni_X;
        
        wJ = w * detJ;
        
        Ke = Ke + matB'*matD_voigt*matB * wJ;
    end
    dof_id = [elem_id*2-1; elem_id*2];
    dof_id = dof_id(:);
    K_dense(dof_id,dof_id) = K_dense(dof_id,dof_id)+Ke;
end

%% apply BC
%clamped left
BC = zeros(nNODE,2);
xlo = abs(NODE(:,1)) <1e-3;
BC(xlo,:) = 1;
BC_dof = reshape(BC',nDOF,1);
BCid = find(BC_dof);

K_dense(BCid,:) = 0;
K_dense(:,BCid) = 0;
K_dense(BCid,BCid) = eye(length(BCid));

locF = find((abs(NODE(:,1) - lxy(1)) < 1e-3) & (abs(NODE(:,2) - lxy(2)/2) < 1e-3));
% locF = find(abs(NODE(:,1) - lxy(1)) < 1e-3);
GF = zeros(nDOF,1);
GF(locF*2-1) = 10; % y direction
% GF(locF(1)*2-1) = 5;
% GF(locF(end)*2-1) = 5;

%% solve!
u = K_dense\GF;
u2 = reshape(u,2,nNODE)';
NODE_f = NODE + u2;

figure; clf; 
hold on 
for ee = 1:nELEM
    elem_id = ELEM(ee,:);
    plot(NODE(elem_id([1,2,3,4,1]),1),NODE(elem_id([1,2,3,4,1]),2),'r--')
    plot(NODE_f(elem_id([1,2,3,4,1]),1),NODE_f(elem_id([1,2,3,4,1]),2),'k-')
end