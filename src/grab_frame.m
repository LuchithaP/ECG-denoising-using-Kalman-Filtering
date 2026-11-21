function im = grab_frame(f, pngFile)
%GRAB_FRAME  Render a figure to an exact 1920x1080 RGB image (16:9).
%   Works in headless/batch MATLAB.  If pngFile is given, the image is saved.
    drawnow;
    im = print(f, '-RGBImage', '-r150');
    im = im(1:min(end,1080), 1:min(end,1920), :);
    if size(im,1) < 1080 || size(im,2) < 1920                % pad if needed
        im(end+1:1080, :, :) = 255;  im(:, end+1:1920, :) = 255;
    end
    if nargin > 1
        % Retry: on Windows a PNG that is open in a viewer can be locked.
        for attempt = 1:5
            try
                imwrite(im, pngFile); return;
            catch
                pause(1);
            end
        end
        alt = strrep(pngFile, '.png', '_new.png');
        imwrite(im, alt);
        warning('%s is locked (open in a viewer?) - saved as %s', pngFile, alt);
    end
end
