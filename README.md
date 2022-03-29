# SMP-PHAT: Lightweight DoA Estimation by Merging Microphone Pairs

## Installation

Make sure FFTW is installed:

```
sudo apt-get install libfftw3-dev
```

Make sure CMake is installed:

```
sudo apt-get install cmake
```

Create a `bin` folder and select it:

```
mkdir bin
cd bin
```

Run CMake:

```
cmake ..\
```

Compile:

```
make
```

## Test speed

In the bin folder:

```
./test
```

This will take some time (few tens of seconds max), and print the performance in the terminal.

## Perform DoA estimation

Assuming the audio `in.wav` contains as many channels as the number of microphone of the arrays (choose between `respeaker_usb`, `respeaker_core`, `minidsp_uma` and `matrix_creator`), run the following command:

```
./ssl -i in.wav -a respeaker_usb -m smp -o out.csv
```

The estimated DoA will be written in CSV format to the file `out.csv`.
