classdef visualization_tools
    %VISUALIZATION_TOOLS Small static plotting helpers for the receiver
    %   Use as:
    %       visualization_tools.plot_Y(Y)
    %       visualization_tools.plot_RVM(RVM)
    %       visualization_tools.plot_peaks(RVM, peaks)
    %
    %   This is intentionally minimal and avoids any custom style so it works
    %   in plain MATLAB installations.

    methods(Static)

        function plot_Y(Y)
            % Plot magnitude of Y (DD-domain / TF-domain depending on convention)
            figure('Name','Y magnitude','NumberTitle','off');
            imagesc(abs(Y));
            axis xy;
            colorbar;
            xlabel('M (columns)');
            ylabel('N (rows)');
            title('Magnitude of received Y');
        end

        function plot_RVM(RVM)
            figure('Name','High-res RVM','NumberTitle','off');
            imagesc(RVM);
            axis xy;
            colorbar;
            xlabel('RVM column');
            ylabel('RVM row');
            title('High-resolution Delay-Doppler (RVM)');
        end

        function plot_peaks(RVM, peaks)
            figure('Name','RVM peaks','NumberTitle','off');
            imagesc(RVM); hold on;
            % peaks is expected as L x 2 [row, col]
            if ~isempty(peaks)
                plot(peaks(:,2), peaks(:,1), 'rx', 'MarkerSize', 12, 'LineWidth', 2);
            end
            axis xy;
            colorbar;
            title('Detected Peaks on RVM');
            hold off;
        end

    end
end
