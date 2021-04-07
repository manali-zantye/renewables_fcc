# renewables_fcc
This repository contains the optimization code for renewable-integrated flexible carbon capture problem. The optimization is demonstrated for a single power plant.

## Getting Started
The user will also have to download an excel workbook containing the complete weather data for all the stations [here](https://drive.google.com/file/d/19mJfep-_1xcI2yZVngkCkjf7cvUXjYLh/view?usp=sharing). This file should be added to [data_files](data_files).

## Overview
[main.py](main.py) implements the two-stage optimization algorithm to determine optimal system design and operation through the following 5 steps:

Step 0: Obtain the time-series data for electricity price, solar and wind capacity factors for the selected power plant.

Step 1: Calculate parameter vector and formulate original scenario tree.

Step 2: Specify the desired scenario reduction accuracy.

Step 3: Perform scenario reduction using [my_SR_code.gms](scen_red/my_SR_code.gms) and obtained reduced scenario data.

Step 4: Use reduced scenario data to obtain optimal system design using [main.gms](red_opti/main.gms).

Step 5: Fix the design decisions and solve with original dataset to obtain optimal operation using [main.gms](orig_opti/main.gms).

