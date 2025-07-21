# Relative Intensity Analysis of Axon Initial Segments

This project contains MATLAB scripts for analyzing the relative intensity of fluorescence signals along the Axon Initial Segment (AIS) from 2D microscopy images. The analysis involves identifying the AIS, measuring its properties (length, location), and quantifying fluorescence intensity profiles of different protein channels.

## Project Structure

- **`ais.m`, `ais_z3.m`**: Main scripts for single-cell analysis. These scripts guide the user through loading images, defining the axon path, and performing the analysis. `ais_z3.m` is a newer version adapted for images from a Zeiss confocal microscope.
- **`x3_channel_script.m`, `x3channel_modified_MH.m`, `mbp.m`**: Scripts for batch processing and analysis of multiple datasets. They automate the process of loading, smoothing, normalizing, and averaging data from multiple experiments.
- **`flib/`**: A library of MATLAB functions used by the main scripts for data manipulation, including loading, smoothing, normalization, interpolation, and plotting.
- **`input/`**: Directory containing the input data files. Each file typically represents a single cell's fluorescence intensity profile, with the first column being the distance along the axon and subsequent columns representing the intensity values for different channels.
- **`output/`**: Directory where the analysis results, such as averaged data and plots, are saved.

## How it Works

The analysis process can be broken down into the following steps:

1.  **Image Loading and Pre-processing**: The scripts load 2D microscopy images (e.g., TIFF files) of neurons. For multi-channel images, each channel is processed separately.
2.  **Axon Tracing**: The user manually traces the axon of interest in the image. The script records the coordinates of the traced path. These coordinates can be saved and re-used for later analysis.
3.  **Fluorescence Intensity Profiling**: The script measures the fluorescence intensity along the traced axon for each channel. It calculates both the raw intensity and a smoothed intensity profile using a 3x3 pixel averaging window and a sliding mean.
4.  **AIS Detection**: The start and end of the AIS are determined based on a fluorescence intensity threshold (by default, 33% of the maximum intensity).
5.  **Parameter Calculation**: The script calculates various parameters for each channel, including:
    *   AIS start, end, and midpoint positions
    *   AIS length
    *   Position of maximum fluorescence intensity
    *   Integrated fluorescence intensity within the AIS
6.  **Data Normalization and Averaging (for batch processing)**: For batch analysis, the scripts normalize the data (both in terms of distance along the axon and fluorescence intensity) and then average the data from multiple cells to generate a mean intensity profile with standard error.
7.  **Output Generation**: The results of the analysis, including calculated parameters and plots of the intensity profiles, are displayed to the user and can be saved to the `output/` directory.

## How to Use

1.  **Prepare your data**: Place your input data files in the `input/` directory. The data should be in a text format where the first column is the distance along the axon and subsequent columns are the fluorescence intensity values for each channel. For image analysis, you will need TIFF files.
2.  **Run the analysis**:
    *   For single-cell analysis, run `ais('your_cell_name')` or `ais_z3('your_cell_name')` in the MATLAB command window. The script will guide you through the process of tracing the axon and analyzing the data.
    *   For batch processing, configure one of the batch scripts (`x3_channel_script.m`, `x3channel_modified_MH.m`, or `mbp.m`) to point to your input directory and then run the script.
3.  **View the results**: The analysis results will be displayed in the MATLAB command window and as plots. The output data can be saved to the `output/` directory.

## Dependencies

This project requires MATLAB. No external toolboxes are required beyond the functions provided in the `flib/` directory.
