%INSTALL_LINUX Script used to install mxTV on the Linux platform.
%
% Compiles and links the mex files for the mxTV package.
%
% Tested with the gcc compiler.
%
% See readme.txt for further instructions.

% J. Dahl^1, P.C. Hansen^2, S.H. Jensen^1 & T.L. Jensen^1
% CSI project: (1) Aalborg University, (2) Technical University of Denmark
% April 28, 2009.

any_error = false;

% Denosing.
try
    cs = sprintf('mex src/c/tv_denoise.c src/c/tools.c src/c/tv_core.c');
    eval(cs)
catch
    any_error= true;
    disp('-------------------------------------------------------------------')
    disp('You will not be able to use TVdenoise because the compilation failed.')
    disp('Follow the above instructions to locate the problem.')
end

% Inpainting.
try
    cs = sprintf('mex c/tv_inpaint.c c/tools.c c/tv_core.c');
    eval(cs)
catch
    any_error= true;
    disp('-------------------------------------------------------------------')
    disp('You will not be able to use TVinpaint because the compilation failed.')
    disp('Follow the above instructions to locate the problem.')
end

% Deblurring.
deblur_error= false;
try
    cs = sprintf('mex c/mxtrp.c c/tools.c');
    eval(cs)
    
    cs = sprintf('mex -DFFTW3 -lfftw3 c/tv_deblur.c c/tv_core.c c/tools.c');
    eval(cs);

    cs = sprintf('mex -DFFTW3 -lfftw3 c/tv_deblur_rr.c c/tv_core.c c/tools.c');
    eval(cs)
catch
    any_error= true;
    disp('-------------------------------------------------------------------')
    disp('You will not be able to use TVdeblur because the compilation failed,')
    disp('probably because the fftw3.h header file is missing. To locate the') 
    disp('problem, follow the above instructions or see the readme.txt file.')
    disp('Ignore this error if you do not need TVdeblur.')
    deblur_error = true;
end

% Check that fftw has the correct setup.
if deblur_error == false
    try	
        xtemp=TVdeblur(1,1,0.1);
        clear xtemp
    catch
        any_error = true;
        disp('----------------------------------------------------------------')
        disp('You will not be able to use TVdeblur because the setup of fftw3')
        disp('is wrong or fftw3 is missing. Follow the instructions in readme.txt')
        disp('to solve the problem if you would like to use TVdeblur.')
        disp('Ignore this error if you do not need TVdeblur.')

    end
end

if any_error == false && deblur_error == false
    disp('Install completed successfully.')
else
    disp('Installation did not complete successfully.')
end