### HDF5 Data File Format (.h5 file)

An .h5 file is a type of the Hierarchical Data Formats (HDF) that is
used to store large, well-organized multidimensional data for quick
retrieval and analysis. `smite` extracts contents of an .h5 file into
H5Structure. The structure consists of the `Data` and `Attributes`
associated with each group. The data can be loaded and extracted. In
Figure 1 cell data is loaded.

<IMG SRC="H1.png" WIDTH=50% HEIGHT=50%><BR>Figure 1

Each level of the HDF5 file hierarchy is named in the 'Groups'
structure. For the data included here, we have two layers:
/Channel01/Zposition001/.

<IMG SRC="H2.png" WIDTH=50% HEIGHT=50%><BR>Figure 2

Zposition001/ contains various groups with their group attributes.

<IMG SRC="H3.png" WIDTH=50% HEIGHT=50%><BR>Figure 3

Zposition001/ can be further extracted.

<IMG SRC="H4.png" WIDTH=50% HEIGHT=50%><BR>Figure 4

As can be seen from the above figure, Data0001 is the 2D microscopy raw
data from the camera. A number of datasets can be written in the HDF5
file. The figure example only contain 5 datasets. The `LoadData` class in
`smite` has the functionality of loading and extracting further information
about HDF5 data.

`/Channel01/Zposition001/Data0001/Data0001` of size 5000 x 256 x 256
contains the raw data from the camera (arrays given in Analog to Digital Units
[ADU]) for 5000 frames of 256x256 (y, x) 2D images.  This dataset is the only
one absolutely needed for further analyses in ***smite***.

---

HDF5 files can be directly loaded into MATLAB.  As an example, to
access the images for Sample 1 located at

   https://datadryad.org/stash/dataset/doi:10.5061/dryad.xsj3tx9cn

run the following commands in MATLAB:

1) Find the number of sequences for a sample:  
```
   Info = h5info('LifeactApproach-HeLaCell-Sample1.h5')
```

2) Each level of the h5 file hierarchy is named in the 'Groups' structure.
For the data included here, we have two layers: /Channel01/Zposition001/.
To find the number of the sequence:
```
   Num_Dataset = numel(Info.Groups.Groups.Datasets)
```

3) To load a data set:
```
   DataSet_1 = h5read('LifeactApproach-HeLaCell-Sample1.h5', ...
                      '/Channel01/Zposition001/Data0001');
```
