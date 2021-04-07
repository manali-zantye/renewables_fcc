# renewables_fcc
This repository contains the optimization code for renewable-integrated flexible carbon capture problem. The optimization is demonstrated for a single power plant.

## Getting Started
The user will also have to download an excel workbook containing the complete weather data for all the stations [here](https://drive.google.com/file/d/19mJfep-_1xcI2yZVngkCkjf7cvUXjYLh/view?usp=sharing). This file should be added to [data_files](data_files).

## Overview
[main.py](main.py) implements the two-stage optimization algorithm to determine optimal system design and operation through the following 5 steps:

### Step 0: 
Obtain the time-series data for electricity price, solar and wind capacity factors for the selected power plant.
