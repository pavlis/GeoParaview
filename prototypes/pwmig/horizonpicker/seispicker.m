%
% Plots data in matrix d and enables stock matlab picker using mb1.
% User should terminate picks with return key 
% Picks are plotted on figure before exiting
function [a,b]=seispicker(d,back,xp,yp)
    figure;
    %`imagesc(back,clim);
    imagesc(back);
    colormap(hsv);
    wigb(d);
    plot(xp,yp,'-bo','LineWidth',2);
    [a,b]=ginput;
    plot(a,b,'-rs','LineWidth',3);
end
