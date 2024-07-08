# preprocessing_pipeline_code

This code performs preprocessing of wide-field calcium imaging in mice, designed for use by the Ebner Lab at the University of Minnesota. These techniques are used in a [peer-reviewed scientific article](https://doi.org/10.1093/cercor/bhab373).  

The entire pipleine should be run from the pipeline_preprocessing.m script in the Control Scripts folder. This pipeline also requires downloading code from the code_common_to_pipelines repository.

This code calls on code created by other programmers, including [dftregistration](https://www.mathworks.com/matlabcentral/fileexchange/18401-efficient-subpixel-image-registration-by-cross-correlation), [PolyDraw](https://www.mathworks.com/matlabcentral/fileexchange/49733-subroutines-polydraw), and an altered version of [tiffread2](https://www.mathworks.com/matlabcentral/fileexchange/10298-tiffread2-m). Credit belongs to the original creators of each part of this code, and should be cited accordingly.

This respository is a work in progress. For questions, problems, or concerns, please email Sarah West at [west0883@umn.edu](west0883@umn.edu).

Please cite use of this code in any resulting projects or publications as: <br>
West, SL. (2021). MATLAB Wide-field calcium imaging preprocessing pipeline code. Accessed [date accessed]. https://github.com/west0883/preprocessing_pipeline_code
