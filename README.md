# ES-Project

This repository contains the code generated for the Embedded Systems Course project.
There are different mains for the different measurements we have done to characterize the behavior of the **B-L072Z-LRWAN1** Discovery kit board by STMicroelectronics.

Each code is derived from the example code available with the board on the STM32 Cube IDE.

## Boards code

### main.c

It's the original example code, available as a library for the STMCubeIDE [here](https://www.st.com/en/embedded-software/i-cube-lrwan.html).

Version 1.3.1

### main_receiver_uart.c and main_sender_uart.c

Modified code to be able to select the SF and Bandwidth on both receiver and sender boards through the *UART* interface; **main_receiver** sets the board as master by default while **main_sender** sets it to slave.
It implements the same Ping-Pong algorithm in **main.c**.

### ping.c

Code flashed only inside the sender board (battery operated). It sends only one *PING* message and then puts the LoRa radio to sleep.

## MATLAB 
This code was used to elaborate the measurements

### R_measures.m

Code used to compute the power and the energy produced by the two different buttons when connected to a resistivve load.

### SF_B_elab.m

Code used to compute and plot the energy consumed during a ping pong from the board at different values of spreading factor and bandwidth.

### Ping_elab.m

Code used to filter and elaborate the current measurements when the board send a single ping. This code also compute the energy consumption of the board and compare the single components (LoRa, MCU) consumes against the total ones.

### calcoli_capacita.m

Code used to compute the required capacity for our circuit given its energy consumption. This code also produce a rough estimate of the number of button presses required to charge the capacitor.

### Final_measures.m

Code used to elaborate the measurements done on the final prototype. Each section of this file is used to elaborate a different measurements , in some we simply filtered the signal and compared it with previous measurements and in others we conted the number of pressures required to charge the capacitor.
