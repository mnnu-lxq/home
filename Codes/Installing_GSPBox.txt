% Check that GSPBox is installed properly
clear;
close all;
G = gsp_logo;
% display the graph
figure;
gsp_plot_graph(G);
G1 =  gsp_minnesota;
figure;
gsp_plot_graph(G1);
title('Minnesota Road Network');
