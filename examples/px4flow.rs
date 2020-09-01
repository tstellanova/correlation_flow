/*
Copyright (c) 2020 Todd Stellanova
LICENSE: BSD3 (see LICENSE file)
*/

#![no_main]
#![no_std]


use cortex_m_rt as rt;
use rt::entry;


use panic_rtt_core::{self, rprintln, rtt_init_print};

use correlation_flow::fwht;

static IMAGE0: &'static [u8] = include_bytes!("../testdata/64sq_253_46.gray");
static IMAGE1: &'static [u8] = include_bytes!("../testdata/64sq_250_30.gray");


#[entry]
fn main() -> ! {
    rtt_init_print!(NoBlockTrim);
    rprintln!("--> MAIN --");

    // setup clock and so forth on the PX4FLOW
    let _peripherals = px4flow_bsp::peripherals::setup_peripherals();

    rprintln!("start flow calcs: {}", IMAGE0.len());
    rprintln!("buf0: {:?}", &IMAGE0[0..8]);
    rprintln!("buf1: {:?}", &IMAGE1[0..8]);

    const COLS: usize = 64;
    const ROWS: usize = 64;

    let mut correlator = fwht::HadamardCorrelator::new(COLS, ROWS);

    const FRAME_COUNT: u32 = 10;
    let mut last_flow = (0i16, 0i16);
    for _ in 0..FRAME_COUNT {
        last_flow = correlator.measure_translation(&IMAGE1, &IMAGE0);
        rprintln!("flow: {:?}", last_flow);
    }
    rprintln!("estimated flow: {:?}", last_flow);
    //TODO use eg DWT::get_cycle_count() to measure elapsed time

    rprintln!("<--- DONE --");

    loop {

    }

}
