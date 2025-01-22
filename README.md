# Satellite Data Processing Challenge

## Overview
This repository addresses a challenge centered around satellite data processing. The objective is to process provided satellite data files, extract useful information, and create meaningful visual and analytical outputs.

## Challenge Tasks

### 1. File Description
- Analyze and describe the contents of the provided files:
  - `file1`
  - `file2`
  - `file3`
  - `file4`

### 2. Data Extraction
- Extract as much information as possible from the provided files. This includes:
  - Metadata (e.g., projection, resolution, timestamp, etc.)
  - Raster or vector data specifics

### 3. Visual Representation
- Generate meaningful image representations of the data, such as:
  - Heatmaps
  - Topographic maps
  - Any other visualization relevant to the data context

### 4. Reprojection
- Reproject all images to a common projection system for consistency.

### 5. Data Merging
- Merge the data from all compatible files into a single comprehensive image, ensuring alignment and consistency.

### 6. Weather Insights
- If possible, extract and provide useful weather-related information from the data.

### 7. Summary
- Summarize the work performed, including:
  - Key observations
  - Challenges encountered
  - Results achieved

## Requirements

### Programming Language
- **Python**

### Suggested Libraries
- [osgeo (gdal, osr)](https://gdal.org/): For handling geospatial data
- [numpy](https://numpy.org/): For numerical operations
- [rasterio](https://rasterio.readthedocs.io/): For raster data manipulation
- [matplotlib](https://matplotlib.org/): For creating visualizations

## Instructions
1. Clone this repository:
   ```bash
   git clone https://github.com/StormGeo/dataops-challenge-satellite
   ```
2. Push the changes to your personal GitHub account:
   - Create a new repository in your personal GitHub account.
   - Push the cloned repository to your personal repository:
     ```bash
     git remote add origin https://github.com/your-username/your-repo-name.git
     git push -u origin main
     ```
4. Add me as a participant in your repository to review the results:
   - Go to your repository on GitHub.
   - Navigate to **Settings** > **Manage Access**.
   - Add my GitHub username as a collaborator: manelceperuelo

5. You will present the results in the second-round interview. Be prepared to discuss:
   - Your approach to solving the tasks
   - Challenges encountered
   - Tools and methods used

## Notes
- The provided images may be in different projections, requiring careful handling during reprojection and merging.
- The final outputs should include both visual and textual summaries to demonstrate the results.

## Deliverable suggestions
1. Descriptions of the files
2. Image representations
3. Reprojected images
5. Merged images
6. Weather-related insights (if applicable)
7. A summary of the work

## Questions?
You are free to organize and structure the project and the code as you see fit.
If you have any questions, please contact us and we'll help you out. This is not an exam.

---

Happy coding and good luck with the challenge!
