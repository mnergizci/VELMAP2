function plot_insar_fields(insar, outdir)
% Plot and save stackmap, los, and azi for each insar entry with vik colormap

    if nargin < 2
        outdir = 'insar_plots';
    end
    if ~exist(outdir, 'dir')
        mkdir(outdir);
    end

    % Load the custom colormap (vik)
    S2 = load('vik.mat');     
    vik = S2.vik;             

    for i = 1:length(insar)
        figure('Visible', 'off', 'Position', [100, 100, 1500, 500]);

        % Subplot 1: stackmap
        subplot(1, 3, 1);
        imagesc(insar(i).stackmap,  'AlphaData', ~isnan(insar(i).stackmap));
        axis image off;
        title('Stackmap');
        colormap(vik);
        c = colorbar;
        c.Label.String = '% plotting';

        % Subplot 2: los
        subplot(1, 3, 2);
        imagesc(insar(i).los, 'AlphaData', ~isnan(insar(i).los));
        axis image off;
        title('LOS');
        colormap(vik);
        c = colorbar;
        c.Label.String = '% plotting';

        % Subplot 3: azi
        subplot(1, 3, 3);
        imagesc(insar(i).azi,'AlphaData', ~isnan(insar(i).azi));
        axis image off;
        title('Azimuth');
        colormap(vik);
        c = colorbar;
        c.Label.String = '% plotting';

        % Output file naming
        if isfield(insar(i).proc, 'ID')
            filename = fullfile(outdir, [insar(i).proc.ID, '.png']);
        else
            filename = fullfile(outdir, sprintf('insar_%02d.png', i));
        end

        % Save and close
        saveas(gcf, filename);
        close(gcf);
    end

    disp(['Saved ', num2str(length(insar)), ' figures in ', outdir]);
end
