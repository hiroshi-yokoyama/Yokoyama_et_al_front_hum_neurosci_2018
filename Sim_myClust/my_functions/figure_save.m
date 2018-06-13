function figure_save(fig, fname)
    set(gca, 'FontName', 'Arial')
    print(fig, '-depsc',[fname,'.eps']);
    print(fig, '-dpng',[fname,'.png'], '-r300');
end