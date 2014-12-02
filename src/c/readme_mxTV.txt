****************************************************************
*                                                              *
*                  mxTV - a Matlab package for                 *
*             Total Variation Image Reconstruction             *
*                                                              *
*             Requires Matlab version 7.5 or later             *
*                                                              *
*              This version dated August 14, 2012              *
*                                                              *
****************************************************************

This package includes Matlab and C codes for Total Variation (TV)
image reconstruction: denoising, inpainting, and deblurring.

If you use this package, please give reference to:

   J. Dahl, P.C. Hansen, S.H. Jensen and T.L. Jensen, Algorithms
   and Software for Total Variation Image Reconstruction via
   First-Order Methods, Numerical Algorithms, 53, pp. 67--92, 2010.

The code is part of the project CSI: Computational Science in
Imaging, supported by the Danish Research Council for Technology
and Production Sciences.  The work was carried out at Aalborg
University and Technical University of Denmark.

Note that the FFTW software is needed only for TV deblurring.
If the FFTW code is not installed (as described below) the codes
for denoising and in-painting will still work.


Installation guide for Windows
------------------------------

1. Go to http://www.fftw.org/install/windows.html and get the pre-
   compiled FFTW for Windows.  To determine if your computer uses
   32 or 64 bits, execute the Matlab command mexext, which returns
   a string.  If the string includes the characters '64' then your
   computer uses 64 bits, otherwise it uses 32 bits.

2. Unzip the downloaded file and put it in a new folder, for example:
      c:\Program Files\fftw-zzz-dll
   where zzz denotes the version of the files.

3. Go to Start -> Control Panel -> System, select Advanced and then
   Environmental Variables.  In the bottom half of the window, under
   System variables, select Path, and append the full path to the
   directory where you unzipped the FFTW.  You do this by adding the
   path after the last ";" and finish with a ";".

4. Go to the directory where you keep your other Matlab toolboxes,
   and unzip the mxTV.zip files to a new folder mxTV.

5. Start Matlab and go to the above mxTV folder.  If Matlab is already
   running at this point, then restart Matlab to make the change in
   the path environment variable (step 3) effective in Matlab

6. Run "mbuild -setup". For older version of Matlab,select the Lcc 
   compiler. Otherwise select, if available, the Microsoft Software
   Development Kit (SDK), or download it by following the provided
   link (then restart Matlab, run "mbuild -setup" again, and select 
   Microsoft Software Development Kit (SDK) as compiler).

7. If you use the Lcc compiler
   run install_windows(<fftwpath>,'externlib\Lcc\libfftw3-3.lib')
	 
   If you use the Microsoft Software Develoment Kit (SDK) C compiler
   run install_windows(<fftwpath>,'externlib\msdk\libfftw3-3.lib')
 	 
   <fftwpath> is a string for the path to the fftw directory in item
   2 above. Note that if the path contains spaces then use the
   format fftwpath='"C:\Program Files\fftw-zzz-dll"'.
	 
   If you use a different compiler, or if the provided import
   libraries does not work, you will see a linking error with 
   unresolved external symbols. In this case you will need to 
   generate your own import libraries, see e.g.
   http://www.fftw.org/install/windows.html. Then run install_windows
   where you provide the path to the new *.lib you generated.
	 
8. Add mxTV to Matlab's path: go to File -> Set Path -> Add Folder
   and choose the folder where mxTV is located. Then save and close.
   Alternatively, you can use the addpath command in Matlab.

9. To learn more, try the three demos:
      TVdenoise_demo, TVinpaint_demo, and TVdeblur_demo.

This installation is tested with the Lcc compiler bundled with Matlab.
More options are available in install_windows.m.


Installation guide for Linux and Unix
-------------------------------------

1. Get the fftw3 software from www.fftw.org/download.html or via your
   package manager, and install it.  Two packages are needed: libfftw3-3
   (base software) and libfftw3-dev (development software, needed for
   header files).

2. Go to the directory where you keep your other Matlab toolboxes, and
   unzip the mxTV.zip files to a new directory mxTV.

3. Start Matlab and go to the above mxTV directory.

4. Run "install_linux".  Our experience is that when installing mxTV,
   you can ignore any warnings deriving from an officially unsupported
   version of gcc.

5. Add mxTV to Matlab's path: go to File -> Set Path -> Add Folder and
   choose the folder where mxTV is located. Then save and close.
   Alternatively, you can use the addpath command in Matlab.  You may
   have to restart Matlab at this point.

6. To learn more, try the three demos:
      TVdenoise_demo, TVinpaint_demo, and TVdeblur_demo.

If you have any problems on Windows, Linus or Unix, please check the
files install_windows.m and install_linux.m.
