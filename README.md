## Darrell Stepter  NE 204 Lab 0

## Energy Calibration with Reproducible Workflows

Data was collected with a Coaxial HPGe detector with Co-60,
Am-241, Ba-133, Cs-137 and Eu-152 sources.

A python script was written to perform an energy calibration for the
Am-241 and Cs-137 sources. Results are presented via a calibrated Ba-133 energy
spectrum and a statistical analysis of the calibration.

Instructions for using the Makefile to run the calibration script and
compiling the LaTeX lab report are posted below:

## File instructions

### Step 1: Download the data

```
make data

```
### Step 2: Ensure the data is not corrupted:

```
make validate
```

### Step 3: Test Python script is properly storing data:

```
make test
```

### Step 4: Run calibration script

```
make analysis
```

### Step 5: Compile the lab report using outputs from the calibration script

```
make
```
