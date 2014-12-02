function any_error = install_windows(fftwpath,libfile)
%INSTALL_WINDOWS Script used to install mxTV on the Windows platform.
%
% Compiles and links the mex files for the mxTV package.
%
% Usage: install_windows(<fftwpath>,<libfile>), where
% <fftwpath> is a string with the path to where fftw was placed, and 
% <libfile> is the path to an fftw import library file (*.lib)
% that can be used with the selected compiler.
%
% The file externlib/Lcc/libfftw3-3.lib is an import library file
% for the Lcc compiler
%
% The file externlib/msdk/libfftw3-3.lib is an import library file
% for the Microsoft Software Development Kit (SDK)
%
% Tested with the Lcc compiler bundled with older versions of Matlab.
% and Microsoft Software Development Kit (SDK)
%
% See readme.txt for further instructions.
%
% J. Dahl^1, P.C. Hansen^2, S.H. Jensen^1 & T.L. Jensen^1
% CSI project: (1) Aalborg University, (2) Technical University of Denmark
% April 28, 2009.
%
% If you want to be able to break the execution of the programs, try to
% set CTRLCBREAK = 1, which uses a non-documented MATLAB API.
% If you do not have libut.lib in your Matlab dir, try CTRLCBREAK = 2.
% Default CTRLCBREAK=0.
%

CTRLCBREAK=0;

if CTRLCBREAK==0
    sbreak = '';
elseif CTRLCBREAK==1
    sbreak = ['-DLIBUT -L' matlabroot '\extern\lib\win32\lcc -llibut'];
elseif CTRLCBREAK==2
    sbreak = ['-DLIBUT -Lexternlib -llibut'];
else
    error('Not a valid option for CTRLCBREAK')
end


ext = mexext;
if strcmp(ext(end-1:end),'64')
	arg = '-largeArrayDims';
else
	arg = '';
end

any_error = false;

% Denosing.
try
    cs = sprintf('mex %s %s c/tv_denoise.c c/tools.c c/tv_core.c',sbreak,arg);
    eval(cs)
catch
    any_error= true;
    disp('-------------------------------------------------------------------')
    disp('You will not be able to use TVdenoise because the compilation failed.')
    disp('Follow the above instructions to locate the problem.')
end

