str='load psft.mat'
eval(str);

for pxsize = 5

    % irradiating backplane
    for angle = 10:5:40
        ['irradiating from ','back plane at angle = ', num2str(angle), ...
            ' and pixel size = ', num2str(pxsize) ]
        
        g4out_name = ['g4out_', num2str(angle), '.mat']
        load(g4out_name)

        tracks = cell(10000,1);
        for i = 1:10000
            try
                [H,S]=getHistory(g4out,i);
                track_temp = [H; S];
                primary_length = length(H);
                secondary_length = length(S);
                input_temp = track_temp(:,[2,3,4,6]);
                total_length = length(input_temp);


                x = track_temp(:,2);
                y = track_temp(:,3);
                z = track_temp(:,4);
                dE = track_temp(:,6);

                randx=pxsize*(rand);
                randy=pxsize*(rand);

                phi=90*rand;

                init_xy = [x, y]';
                rotation_mx = [cosd(phi), -sind(phi);...
                                sind(phi), cosd(phi)];
                rotate_xy = rotation_mx * init_xy;

                x_in = rotate_xy(1,:)+randx;
                y_in = rotate_xy(2,:)+randy;
                % convert unit to DiffuseTrack7
                x_in = x_in' ;
                y_in = y_in' ;
                z_in = (z-325);
                pos_in = [x_in, y_in, z_in];
                D = DiffuseTrack7(psft, 'manual', [pos_in, dE], 'depthcoordinates', [0, -0.65], 'pixelsize', pxsize);
                D.rotate = phi;
                D.shift = [randx, randy];
                D.primary = [pos_in(1:primary_length,:), dE(1:primary_length)];
                D.secondary = [pos_in(primary_length+1:end,:), dE(primary_length+1:end)];
                tracks{i,1} = D;

            catch exception
                continue
            end
        end
        pxsizename = num2str(floor(pxsize))
        trackName = ['ePar', num2str(angle), 'tracks_bk_', pxsizename,'.mat']
        save(trackName, 'tracks')
        clearvars g4out tracks
        g4out = [];
    end
end
