function [K_dense,F_int] = f_Kirchhoff(u,mE,mnu)
global NODE ELEM
[mLambda, mmu] = f_lame(mE, mnu, "E_nu", "Lambda_mu");

% -----------------------------------------------
nNODE = size(NODE,1); nELEM = size(ELEM,1);
dpn = 2; dpe = 8; npe = 4;
nDOF = nNODE*dpn; 
% -----------------------------------------------
K_dense = zeros(nDOF, nDOF);
F_int   = zeros(nDOF,1);
% -----------------------------------------------
ri = [-1 1 1 -1]/sqrt(3);
si = [-1 -1 1 1]/sqrt(3);
wi = [1 1 1 1];
% -----------------------------------------------

%% K matrix generation


for ee = 1:nELEM
    elem_id = ELEM(ee,:);
    X = NODE(elem_id,:); % [4x2]
    x = [u(elem_id*2-1),u(elem_id*2)] + X;
    
    LKT  = zeros(dpe, dpe);
    LF   = zeros(dpe,1);
    
    KT_M = zeros(dpe, dpe);
    KT_G = zeros(dpe, dpe);

    for gg = 1:length(wi)
        r = ri(gg); s = si(gg); w = wi(gg);
        % shape function and its derivative ===================
        % Ni   = 0.25*[(1-r)*(1-s) (1+r)*(1-s) (1+r)*(1+s) (1-r)*(1+s)];
        Ni_r = 0.25*[-1*(1-s) +1*(1-s) +1*(1+s) -1*(1+s)];
        Ni_s = 0.25*[-1*(1-r) -1*(1+r) +1*(1+r) +1*(1-r)];
        
        Ni_rs = [Ni_r;Ni_s]; % [2x4]        
        dX_dr = Ni_rs*X;  % Jacobian [dX/dr dY/dr; dX/ds dY/ds]
        detJ = det(dX_dr);
        dx_dr = Ni_rs*x;  % Jacobian [dX/dr dY/dr; dX/ds dY/ds]
        
        % deformation gradient ================================
        % [dr/dX ds/dX; dr/dY ds/dY]*[dx/dr dy/dr; dx/ds dy/ds]
        % = [dx/dX dy/dX; dx/dY dy/dY] = F'
        matF = (dX_dr\dx_dr)'; % dxi_dXj (Holzapfel notation)
        matC = matF'*matF; 
        matE = 0.5*(matC-eye(2));
        
        % Bmatrix =============================================
        Ni_XY = dX_dr\Ni_rs;
        Ni_X  = Ni_XY(1,:);
        Ni_Y  = Ni_XY(2,:);

        matB = zeros(3,dpe);
        for kk = 1:npe
            Bk = zeros(3,2);
            Bk(1,1) = Ni_X(kk);
            Bk(2,2) = Ni_Y(kk);
            Bk(3,1) = Ni_Y(kk);
            Bk(3,2) = Ni_X(kk);
            k_col = (2*kk-1):(2*kk);
            matB(:, k_col) = Bk * matF';
        end

        % Material modelling ===================================
        % Kirchhoff material
        
        % Stress ---------------------------------------------
        PK2 = mLambda*trace(matE)*eye(2) + 2*mmu*matE; % 2nd PK
        PK1 = matF*PK2; % 1st PK
        
        % Tangent moduli -------------------------------------
        C_SE = mE/(1-mnu*mnu)*[1 mnu 0;mnu 1 0;0 0 (1-mnu)/2];
        
        % Local K matrix ========================================
        detJw = detJ*w;
        
        % geometric nonlinearity --------------------------------
        lkg = Ni_XY' * PK2 * Ni_XY * detJw;
        for ki = 1:npe
            for kj = 1:npe
                k_row = (2*ki-1):(2*ki);
                k_col = (2*kj-1):(2*kj);
                KT_G(k_row, k_col) = KT_G(k_row, k_col) + eye(2)*lkg(ki,kj);
            end
        end
        
        % material nonlinearity ---------------------------------
        lkm = matB' * C_SE * matB * detJw;
        KT_M = KT_M + lkm;
        
        LKT = KT_G + KT_M;
        
        % local internal force ===================================
        f_int = Ni_XY' * PK2 * matF' * detJw;
        for k = 1:dpn
            k_pos = k:dpn:dpe;
            LF(k_pos) = LF(k_pos) + f_int(:,k);
        end        
    end
    dof_id = [elem_id*2-1; elem_id*2];
    dof_id = dof_id(:);
    K_dense(dof_id,dof_id) = K_dense(dof_id,dof_id) + LKT;
    F_int(dof_id) = F_int(dof_id) + LF;
end
end