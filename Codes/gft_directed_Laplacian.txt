clear all;close all;clc;
%==========================================================================
%------This Program Computes and plots GFT of a signal (defined on a directed ring Graph)
%------based on directed Laplacian----------
%==========================================================================

%-----------Creating adjacency martix for directed Ring Graph====
N = 8; %Number of nodes in the graph
W = zeros(N,N);
W1 = eye(N);
W(2:N,:) = W1(1:N-1,:);
W(1,:) = W1(N,:);

%%% OR you can load arbitrary weight matrix here====
%load W;

% G.W = sparse(double(W'));
% G = gsp_graph_default_parameters(G);
% G.plotting.vertex_size = 300;
% G.coords = G1.coords;
% G.plotting = G1.plotting
% 
% gsp_plot_graph(G);

%-----------Out and In Degree----------------------------------------------
D_out = diag(sum(W,1));
D_in = diag(sum(W,2));

%---------Out and In Laplacian---------------------------------------------
L_out = D_out - W;
L_in = D_in - W;
%==========================================================================
%---------------------------Graph Signal--------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = rand(N,1);

%==========================================================================
%----------------------Jordan Forms----------------------------------------
%==========================================================================
[V_Lin1 J_Lin1] = jordan(L_in);
[V_Lin J_Lin] = sort_jordan(V_Lin1,J_Lin1);
clear V_Lin1 J_Lin1;

for i=1:N
V_Lin(:,i) = V_Lin(:,i)/norm(V_Lin(:,i));
end


f_hat = (inv(V_Lin)) * f;

%=======================================
%plotting Spectrum
figure;
x = real(diag(J_Lin));
y = imag(diag(J_Lin));
z = abs(f_hat);
r = (x.^2 + y.^2).^(0.5);

th = (-90:1:90)*pi/180;
[TH,R] = meshgrid(th,r);
[X,Y] = pol2cart(TH,R);
zer = zeros(size(TH));

stem3(x,y,z,'filled','linewidth',2.5);
set(gca,'fontsize',25)
h = xlabel('Re$(\lambda_\ell)$','Interpreter','LaTex','FontSize',30);
set(h, 'Units', 'Normalized');
pos = get(h, 'Position'); set(h, 'Position', pos + [0.05, 0.05, 0]);
h = ylabel('Im$(\lambda_\ell)$','Interpreter','LaTex','FontSize',30);
set(h, 'Units', 'Normalized');
pos = get(h, 'Position'); set(h, 'Position', pos + [-0.03, 0.05, 0]);
h = zlabel('$| \hat{f}(\lambda_\ell)~|$','Interpreter','LaTex','FontSize',30);
set(h, 'Units', 'Normalized');
pos = get(h, 'Position'); set(h, 'Position', pos + [0, 0, -0.2]);
box off;

hold on
surf(X,Y,zer,'EdgeColor','none','LineStyle','none','FaceLighting','phong')
for i = 1: size(r)
    rho=ones(1,length(th))*r(i);
    [Xi,Yi] = pol2cart(th,rho); 
    plot(Xi,Yi,'r-');
end

Xi = linspace(0,max(r),10000);
xaxis_line = zeros(1,length(Xi));
plot(Xi,xaxis_line,'k-');

pos = get(gca, 'Position');
set(gca, 'Position', pos)
hold off;


function [U,E] = sort_jordan(Ui,Ei)
    % Sort eigenvectors and eigenvalues
    d_Ei = diag(Ei);
    [Et,inds] = sort(abs(d_Ei),'ascend');
    E = d_Ei(inds);
    Uii=Ui(:,inds);
    E = diag(E);
    U = Uii;
end
