function outputImage = bwmorph_copy(varargin)
%function outputImage = bwmorph_copy(inputImage,op)
%function outputImage = bwmorph_copy(inputImage,op,numIterations)
%
% Replicating functionality of Image Processing Toolbox, for use on lawrencium.
%
% intputImage, outputImage are logical arrays or 0's and 1's (binary images)
% op may be 'dilate' or 'thin'
% numIterations for 'thin' (usually use Inf)

%% Input handling.
if nargin < 2
    error('Need at least two arguments')
%get input image from argument 1
elseif islogical(varargin{1})
    inputImage = varargin{1};
elseif isnumeric(varargin{1})
    inputImage = logical(varargin{1});
end

%get operation type from argument 2
if ischar(varargin{2})
    op = varargin{2};
else
    error('Operation type should be a string')
end

%check for argument 3
if nargin==3 && isnumeric(varargin{3})
    numIterations = varargin{3};
elseif nargin==3 && ~isnumeric(varargin{3})
    error('numIterations should be numeric')
elseif nargin > 3
    error('What do you want me to do with all these input arguments??')
else
    %default to Inf iterations for 'thin'
    numIterations = Inf;
end

%% Dilate

if strcmpi(op,'dilate')
    conn = ones(3,3);   %8-conn
    
    imgTemp = conv2(+inputImage,conn);
    imgTemp = imgTemp(2:end-1,2:end-1);
    
    outputImage = imgTemp>0;
    
%% Thin

elseif strcmpi(op,'thin')
    %start from input...
    currentImage = inputImage;
    newCurrentImage = currentImage;
    imageSize = size(inputImage);
    
    %compute checkerboard masks
%     [X,Y] = meshgrid(1:imageSize(2), 1:imageSize(1));   %switch 1 and 2 to get the right size
%     subfield{1} = (X+Y)/2 == floor((X+Y)/2);
%     subfield{2} = ~subfield{1};
    subfield{1} = true(imageSize);
    subfield{2} = true(imageSize);
    
    
    if isinf(numIterations)
        %break when image stops changing
        %   but track iterations for fun
        n=1;
        while true
            %call function
            newCurrentImage = bwmorphThinIteration(newCurrentImage, imageSize, subfield);
            
            %check for image change
            if all(newCurrentImage(:) == currentImage(:))
%                 disp(['Finished in ',num2str(n),' iterations'])
                break
            end
            
            %update image for next iteration
            %   newCurrentImage is always more updated than currentImage
            currentImage = newCurrentImage;
            
%             disp(num2str(n))
%             disp(' ')
%             newCurrentImage
            n = n+1;
        end     %iterations
    else
        %go through numIterations, or break when image stops changing
        for i=1:numIterations
            
            %call function
            newCurrentImage = bwmorphThinIteration(newCurrentImage, imageSize, subfield);
            
            %check for image change
            if all(newCurrentImage(:) == currentImage(:))
                break
            end
            
            %update image for next iteration
            %   newCurrentImage is always more updated than currentImage
            currentImage = newCurrentImage;
        end     %iterations
    end
    
    %output final image
    %   newCurrentImage is always more updated than currentImage
    outputImage = newCurrentImage;
    
%% Other (error)
    
else
    error('Unsupported operation')
end


function newCurrentImage = bwmorphThinIteration(newCurrentImage, imageSize, subfield)

        %two subiterations, which are identical except for G3 and the subfield
        for k=1:2
            %compute x_i for all pixels. x_1 is east of p; x_2 is northeast; and so on, ccw.
            %   "north" = positive x, i guess? "east" = positive y? where x,y are image subscripts.
            %   this is like how surf and pcolor plot data onto the screen.
            
            %   fixed now: "north" = negative x. so it is how matrix data looks in variable editor.
            %   that just means switching 2 and 8, 3 and 7, and 4 and 6.
            x{1} = [newCurrentImage(:,2:end), false(imageSize(1),1)];      %shift in -y direction
            
            x{8} = [newCurrentImage(2:end,2:end), false(imageSize(1)-1,1); ...
                false(1,imageSize(2))];                                 %shift in -x, -y directions
            
            x{7} = [newCurrentImage(2:end,:); ...
                false(1,imageSize(2))];                                 %shift in -x direction
            
            x{6} = [false(imageSize(1)-1,1), newCurrentImage(2:end,1:end-1); ...
                false(1,imageSize(2))];                                 %shift in -x, +y directions
            
            x{5} = [false(imageSize(1),1), newCurrentImage(:,1:end-1)];    %shift in +y direction
            
            x{4} = [false(1,imageSize(2)); ...
                false(imageSize(1)-1,1), newCurrentImage(1:end-1,1:end-1)];  %shift in +x, +y directions
            
            x{3} = [false(1,imageSize(2)); ...
                newCurrentImage(1:end-1,:)];                               %shift in +x direction
            
            x{2} = [false(1,imageSize(2)); ...
                newCurrentImage(1:end-1,2:end), false(imageSize(1)-1,1)];  %shift in +x, -y directions
            
            x{9} = x{1};    %wrap around for indexing of b_4 and n2
            
            %G1: X_H(p) = 1
            %   X_H = sum(i=1:4)(b_i), b_i = 1 iff x_(2i-1) && (x_(2i) || x_(2i+1))
            
            %compute b_i for all pixels
            b = nan([imageSize,4]);
            for j=1:4   %index of b_j
                b(:,:,j) = ~x{2*j-1} & (x{2*j} | x{2*j+1});
            end
            
            X_H = sum(b,3);
            G1 = X_H==1;
            
            %G2: 2 <= min(n1(p),n2(p)) <= 3
            %   n1(p) = sum(i=1:4)(x_(2i-1) || x_(2i))
            %   n2(p) = sum(i=1:4)(x_(2i) || x_(2i+1))
            tmp1 = nan([imageSize,4]);
            tmp2 = tmp1;
            for j=1:4
                tmp1(:,:,j) = x{2*j-1} | x{2*j};
                tmp2(:,:,j) = x{2*j} | x{2*j+1};
            end
            n1 = sum(tmp1,3);
            n2 = sum(tmp2,3);
            min_n1n2 = min(n1,n2);
            
            G2 = min_n1n2 >= 2 & min_n1n2 <= 3;
            
            if k==1
                %G3a: (x_2 | x_3 | ~x_8) & x_1 == false
                G3 = ~((x{2} | x{3} | ~x{8}) & x{1});
%                 G3 = ~x{1} | (~x{2} & ~x{3} & x{8});
                %   these two should be equivalent...
                %   the second one makes more sense as I look at the Guo paper
            elseif k==2
                %G3b: (x_6 | x_7 | ~x_4) & x_5 == false
                G3 = ~((x{6} | x{7} | ~x{4}) & x{5});
%                 G3 = ~x{5} | (~x{6} & ~x{7} & x{4});
                %   these two should be equivalent...
                %   the second one makes more sense as I look at the Guo paper
            end
            
            %debug
%             figure; 
%             subplot(2,3,1); img = [+newCurrentImage, false(imageSize(1),1); false(1,imageSize(2)+1)]; pcolor(img); axis equal; 
%             subplot(2,3,4); img = [+G1,              false(imageSize(1),1); false(1,imageSize(2)+1)]; pcolor(img); axis equal; 
%             subplot(2,3,5); img = [+G2,              false(imageSize(1),1); false(1,imageSize(2)+1)]; pcolor(img); axis equal; 
%             subplot(2,3,6); img = [+G3,              false(imageSize(1),1); false(1,imageSize(2)+1)]; pcolor(img); axis equal; 
%             subplot(2,3,2); img = [+subfield{k},     false(imageSize(1),1); false(1,imageSize(2)+1)]; pcolor(img); axis equal; 
            
            %subiteration 1
            %   set up mask of pixels to remove.
            removalMask = subfield{k} & G1 & G2 & G3;
            %   remove
            newCurrentImage(removalMask) = false;
%             subplot(2,3,3); img = [+newCurrentImage, false(imageSize(1),1); false(1,imageSize(2)+1)]; pcolor(img); axis equal; 
%             pause;
        end     %subiterations
        
        