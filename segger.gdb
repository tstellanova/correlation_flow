# J-LINK GDB SERVER initialization
#
# This connects to a GDB Server listening
# for commands on localhost at tcp port 2331
target extended-remote localhost:2331

monitor speed 30
monitor reset
#
# CPU core initialization (to be done by user)
#
# Set auto JTAG speed
monitor speed auto

# Setup GDB FOR FASTER DOWNLOADS
#set remote memory-write-packet-size 1024
#set remote memory-write-packet-size fixed

# *try* to stop at the user entry point (it might be gone due to inlining)
# break main

# don't confirm when quitting debugger
define hook-quit
    set confirm off
end

#monitor semihosting enable

load

#break rust_begin_unwind
run

#stepi
