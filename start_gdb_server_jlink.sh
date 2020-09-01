#!/bin/sh

# PX4FLOW board:
JLinkGDBServer -device stm32f407vg --speed 1200 -if swd -port 2331
# Open MV H7 board:
#JLinkGDBServer -device stm32h743vi -speed 1200 -if swd -port 2331

