# correlation_flow

Optical flow calculation using image correlation for no_std rust.
This library will compare small grayscale image blocks (limited to 64x64 pixels currently),
and detect the phase correlation between them. This potentially enables applications such as:
- Create image panoramas by seaming together overlapping images
- Detect optical flow translation (for machine vision / robotics)
- Image registration
- Measure image similarity

## Embedded Examples
The examples are designed to be used with J-Link / RTT.
We provide a couple different ways to run these:
- With [probe-run](https://crates.io/crates/probe-run)
- With the Segger tools

#### With probe-run installed
- Simply run the example (see below) with a JLink debug probe attached to your PX4FLOW
- If you have problems, edit [config](.cargo/config) to ensure that the probe-run runner is selected

#### With segger tools installed 
- Edit [config](.cargo/config) to select the `segger.gdb` runner
- In one shell run: `./start_gdb_server_jlink.sh`
- In another shell run: `JLinkRTTClient`
- Then run your choice of examples


### PX4FLOW 
This example is intended to run on the PX4FLOW hardware.
It simply compares two image frames stored in the app binary. 

```shell script
cargo run --release --example px4flow 
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
- [ ] Support > 64x64 block sizes