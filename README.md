# correlation_flow

Optical flow calculation using image correlation for no_std rust.

## Embedded Examples

The examples are currently designed to be used with J-Link / RTT.
In the future as tools such as probe-rs solidify, we may switch to that toolset.

- In one shell run: `./start_gdb_server_jlink.sh`
- In another shell run: `JLinkRTTClient`
- Then run your choice of examples


### PX4FLOW 

This example is intended to run on the PX4FLOW hardware.
It simply compares two image frames stored in the app binary. 

```shell script
cargo run --example px4flow --release
``` 

## Image Conversion

We used ImageMagick's `convert` command to generate raw 8-bit grayscale 
images from png files using eg: 
```shell script
convert 64sq_253_46.png -depth 8 64sq_253_46.gray
```

## Status

- [x] Detects discrete 2D image translation within a static 64x64 grid (maximum +/- 32 pixel movement)
- [x] Example that runs on embedded hardware (PX4FLOW)
- [ ] Simple test runs (currently tested under std using [glovent](https://github.com/tstellanova/glovent))