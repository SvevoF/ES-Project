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

### SF_B_elab.m

### Ping_elab.m

### R_measures.m

### Final_measures.m