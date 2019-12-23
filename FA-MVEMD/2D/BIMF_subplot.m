function BIMF_subplot(signal1,signal2,imf,name1,name2)

figure
subplot(2,1,1)
    %Masking wall data
    imAlpha=ones(size(signal1));
    imAlpha(isnan(signal1))=0;    
    imagesc(signal1,'AlphaData',imAlpha);
    set(gca,'color',0*[1 1 1]);
    xlabel('$x$')
    ylabel('$y$')
    axis equal;
    axis tight;
    title(sprintf('IMF %d %s ',imf,name1));
    set(gca,'TickLabelInterpreter','latex')
    colormap(parula);
    hcb = colorbar;
    colorTitleHandle = get(hcb,'Title');
    titleString = '$\frac{u}{U_{\infty}}$';
    set(colorTitleHandle ,'String',titleString,'Interpreter','latex','FontSize',14);
    set(hcb,'TickLabelInterpreter','latex');
    
subplot(2,1,2)
    %Masking wall data
    imAlpha=ones(size(signal2));
    imAlpha(isnan(signal2))=0;    
    imagesc(signal2,'AlphaData',imAlpha);
    set(gca,'color',0*[1 1 1]);
    xlabel('$x$')
    ylabel('$y$')
    axis equal;
    axis tight;
    title(sprintf('IMF %d %s ',imf,name2));
    set(gca,'TickLabelInterpreter','latex')
    colormap(parula);
    hcb = colorbar;
    colorTitleHandle = get(hcb,'Title');
    titleString = '$\frac{v}{U_{\infty}}$';
    set(colorTitleHandle ,'String',titleString,'Interpreter','latex','FontSize',14);
    set(hcb,'TickLabelInterpreter','latex');
    
end

