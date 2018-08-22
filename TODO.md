## General

- Move each constant to its appropriate module.
- Standardize XYZ coordinates across submodules. (e.g., anim.nucleus.Nucleus.box_origin has ZYX where Y are the rows and X the columns, while the input of gpseq_fromfish expects X as the rows and Y as the columns).

## `gpseq_fromfish`

- Allow for FISH dot coordinates to be float (i.e., need to interpolate).
- Allow for different segmentation/analysis methods.
- Implement gaussian sum fitting for G1 selection.

## `tiff_auto3dseg`

- Clean `tools.binarization` module and allow for labeled input/output.
