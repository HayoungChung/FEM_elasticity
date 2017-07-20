% test script for Q4 linear FEM
clc
clear all

DEBUGFLAG__ = 1;
ECHOFLAG__  = 1;
PLOTFLAG__  = 0;
SAVEFLAG__  = 0;

global NODE ELEM
%% Material model (HyperElasticity)
% WIP BELOW. NEED CHECKUP
% Saint-Venant Kirchhoff model 
mE = 1000; mnu = 0.3;

%% mesh initialization
lxy = [1 1];
exy = [10 10];
nxy = exy + 1;
nNODE = prod(nxy);

dpn = 2;
nDOF = nNODE * dpn; % dpn = 2;
dpe = 8;

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


%% apply BC
%clamped left
BC = zeros(nNODE,2);
xlo = abs(NODE(:,1)) <1e-3;
BC(xlo,:) = 1;
BC_dof = reshape(BC',nDOF,1);
BCid = find(BC_dof);

GF = zeros(nDOF,1);
% locF = find((abs(NODE(:,1) - lxy(1)) < 1e-3) & (abs(NODE(:,2) - lxy(2)/2) < 1e-3));
locF = find(abs(NODE(:,1) - lxy(1)) < 1e-3);
GF(locF*2-1) = 20; % y direction
GF(locF(1)*2-1) = 10;
GF(locF(end)*2-1) = 10;


%% Newton-Raphson solver initiation
MAXITER = 1000;
MAXSTEP = 1000;
XTOL = 1e-6; % 
RTOL = 1e-6; % residual
CTOL = 1e-5; % error if size of increment is smaller than CTOL

initstep = 10;
lambdaR  = 0; % increment at the present timestep
lambdaR0 = 0; % lambdaR at previous timestep
del_lambdaR = 1/initstep;

GU_u = zeros(nNODE,2); % initial displacement

for iSTEP = 1:MAXSTEP
   xm_u = reshape(GU_u',nDOF,1);
   
   % termination
   if abs(lambdaR - 1) < 1e-3
       [GKT] = f_Kirchhoff(xm_u,mE,mnu);
       break
   end
   
   % CTOL error
   if (iSTEP > 1) && (del_lambdaR < CTOL)
       error('step increment is smaller than CTOL');
   end
   
   % step increment
   if (iSTEP > 1) && (lambdaR0 + del_lambdaR >= 1) % last step
       lambdaR = 1;
       del_lambdaR = 1-lambdaR0;
   else
       lambdaR = lambdaR0 + del_lambdaR;
   end
   
   for iITER = 1:MAXITER
       F_step = GF*lambdaR;
       
       [GKT,F_int] = f_Kirchhoff(xm_u,mE,mnu);
       Res = F_step - F_int;
       Res(BCid) = 0;
       GKT(BCid,:) = 0; 
       GKT(:,BCid) = 0; 
       GKT(BCid,BCid) = eye(length(BCid)); 
       
       if iITER == 1
           R0 = Res'*Res;
       end
       
       if iITER == 2
           delx0 = norm(del_xm);
       end
       
       % iteration termination
       if (iITER > 1) && ((Res'*Res/R0 < RTOL) || (norm(del_xm)/delx0 < XTOL))
           if ECHOFLAG__
               disp('=====================================')
           end
           lambdaR0 = lambdaR;
           GU_u = reshape(xm_u,2,nNODE)'; 
           if SAVEFLAG__
               save(sprintf('./save/step_%d.mat',iSTEP),GU_u);
           end
           
           if iITER <= 4
               del_lambdaR = del_lambdaR*2;
           end
           
           if PLOTFLAG__
               NODE_f = NODE + GU_u;
               figure;
               hold on 
               for ee = 1:nELEM
                   elem_id = ELEM(ee,:);
                   plot(NODE(elem_id([1,2,3,4,1]),1),NODE(elem_id([1,2,3,4,1]),2),'r--')
                   plot(NODE_f(elem_id([1,2,3,4,1]),1),NODE_f(elem_id([1,2,3,4,1]),2),'k-')
               end
           end
           break;
       end
       
       if iITER == MAXITER
           lambdaR = lambdaR0;
           del_lambdaR = del_lambdaR/2;
           if ECHOFLAG__
               disp('-------------- MAX. ITER. NUM. --------------')
           end
           break
       end
       
       del_xm = GKT\Res;
       
       xm_u = xm_u + del_xm;
       
       if DEBUGFLAG__
          NODE_f = NODE + reshape(xm_u,2,nNODE)';
          figure(100); clf;
          hold on 
          for ee = 1:nELEM
              elem_id = ELEM(ee,:);
              plot(NODE(elem_id([1,2,3,4,1]),1),NODE(elem_id([1,2,3,4,1]),2),'r--')
              plot(NODE_f(elem_id([1,2,3,4,1]),1),NODE_f(elem_id([1,2,3,4,1]),2),'k-')
          end
          pause(0.1);
       end
       
        if ECHOFLAG__
            fprintf('%3d        %3d             %4.3e           %4.3e           %1.8f\n',...
            iSTEP,iITER,sqrt(Res'*Res),norm(del_xm),lambdaR);
        end
   end
end


%% solve!
NODE_f = NODE + GU_u;

figure; clf; 
hold on 
for ee = 1:nELEM
    elem_id = ELEM(ee,:);
    plot(NODE(elem_id([1,2,3,4,1]),1),NODE(elem_id([1,2,3,4,1]),2),'r--')
    plot(NODE_f(elem_id([1,2,3,4,1]),1),NODE_f(elem_id([1,2,3,4,1]),2),'k-')
end