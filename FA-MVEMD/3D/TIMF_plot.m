
function TIMF_plot(signal,Colour,nslice,imf,name1,name2)    

    [Nx, Ny, Nz] = size(signal);

    xslice = linspace(1,Nx,nslice);
    yslice = linspace(1,Ny,nslice);
    zslice = linspace(1,Nz,nslice);
    volume = slice(signal,xslice,yslice,zslice);
    axis equal;
    xlabel('x');
    ylabel('y');
    zlabel('z');
    set(gca,'TickLabelInterpreter','latex')
    switch(name1)
        case 'IMF'
            title(sprintf('%s %d %s',name1,imf,name2));
        case 'Signal'
            title(sprintf('%s %s',name1,name2));
        case 'Residue'
            title(sprintf('%s %s',name1,name2));
    end
    colorbar;
    set(volume,'EdgeColor','none',...
        'FaceColor','interp',...
        'FaceAlpha','interp')
    alpha('color')
    view(30,30);
    alphamap('rampup')
    alphamap('decrease',.1)
    colormap(Colour);
%     caxis([-3 3]); %TODO: Why commented out in all but EMD3D3V?
    hcb = colorbar;
    colorTitleHandle = get(hcb,'Title');
    titleString = '$\frac{u}{U_{\infty}}$';
    set(colorTitleHandle ,'String',titleString,'Interpreter','latex','FontSize',14);
    set(hcb,'TickLabelInterpreter','latex');
    
end