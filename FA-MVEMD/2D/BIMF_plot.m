function BIMF_plot(signal,Colour,imf,name1,name2)
%Masking wall data
% load('MASK_file','MASK');
% signal = MASK.*signal;

    imAlpha=ones(size(signal));
    imAlpha(isnan(signal))=0;    
    imagesc(signal,'AlphaData',imAlpha);
    set(gca,'color',0*[1 1 1]);
    xlabel('$x$')
    ylabel('$y$')
    axis equal;
    axis tight;

    switch(name1)
        case 'IMF'
            title(sprintf('%s %d %s ',name1,imf,name2));
        case 'Residue'
            title(sprintf('%s %s ',name1,name2));
        case 'Signal'
            title(sprintf('%s %s ',name1,name2));
    end
    set(gca,'TickLabelInterpreter','latex')
    colormap(Colour);
    hcb = colorbar;
    colorTitleHandle = get(hcb,'Title');
    titleString = 'u';
    set(colorTitleHandle ,'String',titleString,'Interpreter','latex','FontSize',14);
    set(hcb,'TickLabelInterpreter','latex');
end