# Workflow Details

This section describes the Labelmerge workflow (i.e. steps taken to produce 
intermediate and final files). Labelmerge is a 
[Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow, and thus a 
directed acyclic graph (DAG) that is automatically configured based on a set of 
rules. A simplified workflow schematic is also shown below.

## Overall workflow
[Labelmerge](https://github.com/khanlab/labelmerge) allows users to combine
segmentations from varying sources. For example, previously segmented
thalamic nuclei to be combined with an atlas of other subcortical structures. 

<img src="simple_workflow_labelmerge.png" width="800px">
_(Click on the image to enlarge)_

Snakemake workflows are organized into groups of rules, each
representing the different phases of the workflow. Each grouped set of rules 
also exist in separate rule files, which can be found under the 
[rules sub-directories](https://github.com/khanlab/labelmerge/tree/main/labelmerge/workflow/rules) 
in the workflow source. 

At a more granular level, Labelmerge reads one “base” and one “overlay” atlas image organized as BIDS 
derivatives, with the label names supplied as metadata in separate TSV files. 
An image mask is created for both the base and overlay images, mapping each unique 
label name in the parcellation scheme to their respective voxel coordinates. The two masks 
are then merged, with labels reindexed to range from one to the total number of labels across both images. 
An output image is constructed where any voxel that is labeled in both the base and overlay images takes 
the reindexed value from the overlay image. Finally, this image is written alongside a new metadata table 
where the label names are pre-appended with “base “ or “overlay “ as appropriate. 
