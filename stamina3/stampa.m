function stampa(out_dir, title_str, S_fig, S, x, y, z, step, total, dt)
    figure(S_fig);
    slice(S, x, y, z);
    shading interp
    set(gca, 'XDir', 'reverse')
    set(gca, 'ZDir', 'reverse')
    set(gcf, 'color', 'w')
    xlabel('y')
    ylabel('x')
    zlabel('z')
    axis equal
    title([title_str ' concentration profile: 0 < t = ' num2str(step*dt, '%2.2f') ' < ' num2str(total*dt, '%2.2f')])
    colorbar
    %view(90, 90)

    fname = fullfile(out_dir, sprintf("%s_%05i.png", title_str, step));
    imwrite(getframe(gcf).cdata, fname);
end