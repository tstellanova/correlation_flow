/*
Copyright (c) 2020 Todd Stellanova
LICENSE: BSD3 (see LICENSE file)
*/

#![no_main]
#![no_std]

//! This example runs on the PX4FLOW machine vision board, but may also run on other
//! similar stm32h7 hardware.
//!
//! It calculates the translation between two image frames using phase correlation.
//!
use cortex_m_rt as rt;
use rt::entry;

use panic_rtt_core::{self, rprintln, rtt_init_print};
use correlation_flow::micro_rfft;

/// To keep this example simple, we pull grayscale image data into the app binary
static IMAGE0: &'static [u8] = include_bytes!("../testdata/64sq_250_30.gray");
static IMAGE1: &'static [u8] = include_bytes!("../testdata/64sq_253_46.gray");
// static IMAGE1: &'static [u8] = include_bytes!("../testdata/64sq_255_33.gray");

#[entry]
fn main() -> ! {
    rtt_init_print!(NoBlockTrim);
    rprintln!("--> MAIN --");

    // setup clock and so forth on the PX4FLOW
    let peripherals = px4flow_bsp::peripherals::setup_peripherals();

    rprintln!("start flow calcs: {}", IMAGE0.len());
    // rprintln!("buf0: {:?}", &IMAGE0[..8]);
    // rprintln!("buf1: {:?}", &IMAGE1[..8]);

    let mut correlator = micro_rfft::MicroFftContext::new();

    // let mut dwt = peripherals.2;
    let mut stopwatch_buf = [0u32; 2];
    let mut stopwatch = peripherals.2.stopwatch(&mut stopwatch_buf);
    const FRAME_COUNT: u32 = 10;
    let mut last_flow = (0i16, 0i16);

    stopwatch.reset();
    for _ in 0..FRAME_COUNT {
        last_flow = correlator.measure_translation(&IMAGE1, &IMAGE0);
        //rprintln!("flow: {:?}", last_flow);
    }
    stopwatch.lap();
    if let Some(lap_time) = stopwatch.lap_time(1) {
        let frac_secs = lap_time.as_secs_f32();
        rprintln!("{} frames in {:.2} seconds: {:.2} fps", FRAME_COUNT, frac_secs, (FRAME_COUNT as f32)/frac_secs);
    }
    rprintln!("forward flow: {:?}", last_flow);
    rprintln!("ground truth: (3, 16)");

    for _ in 0..FRAME_COUNT {
        last_flow = correlator.measure_translation(&IMAGE0, &IMAGE1);
        // rprintln!("flow: {:?}", last_flow);
    }
    rprintln!("reverse flow: {:?}", last_flow);

    rprintln!("<--- DONE --");

    loop {}
}
