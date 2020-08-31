/*
Copyright (c) 2020 Todd Stellanova
LICENSE: BSD3 (see LICENSE file)
*/

#![no_main]
#![no_std]

use stm32f4xx_hal as p_hal;
use p_hal::stm32 as pac;

use cortex_m_rt as rt;
use rt::entry;

use pac::interrupt;

use panic_rtt_core::{self, rprint, rprintln, rtt_init_print};

use embedded_hal::digital::v2::OutputPin;
use embedded_hal::digital::v2::ToggleableOutputPin;

// use correlation_flow::micro_rfft;
use correlation_flow::fwht;

const GYRO_REPORTING_RATE_HZ: u16 = 95;
const GYRO_REPORTING_INTERVAL_MS: u16 = 1000 / GYRO_REPORTING_RATE_HZ;


use core::sync::atomic::{AtomicPtr, Ordering};
use px4flow_bsp::board::Board;
use px4flow_bsp::dcmi::{
    ImageFrameBuf,  FRAME_BUF_LEN,
};

static mut BOARD_PTR: AtomicPtr<Board> = AtomicPtr::new(core::ptr::null_mut());
/// should be called whenever DMA2 completes a transfer
#[interrupt]
fn DMA2_STREAM1() {
    // forward to DCMI's interrupt handler
    unsafe {
        (*BOARD_PTR.load(Ordering::SeqCst)).handle_dma2_stream1_interrupt();
    }
}

/// should be called whenever DCMI completes a frame
#[interrupt]
fn DCMI() {
    // forward to DCMI's interrupt handler
    unsafe {
        (*BOARD_PTR.load(Ordering::SeqCst)).handle_dcmi_interrupt();
    }
}

/// Setup core-coupled RAM buffers for faster image manipulation
#[link_section = ".ccmram.IMG_BUFS"]
static mut FAST_IMG0: ImageFrameBuf = [0u8; FRAME_BUF_LEN];

#[link_section = ".ccmram.IMG_BUFS"]
static mut FAST_IMG1: ImageFrameBuf = [0u8; FRAME_BUF_LEN];

// #[link_section = ".ccmram.IMG_BUFS"]
// static mut  FAST_IMG_BUFS: [ImageFrameBuf; 2] = [[0u8; FRAME_BUF_LEN]; 2];


static IMAGE0: &'static [u8] = include_bytes!("../testdata/64sq_253_46.gray");
static IMAGE1: &'static [u8] = include_bytes!("../testdata/64sq_250_30.gray");


#[entry]
fn main() -> ! {
    rtt_init_print!(NoBlockTrim);
    rprintln!("--> MAIN --");

    // just setup clock and so forth
    let _periphs = px4flow_bsp::peripherals::setup_peripherals();

    rprintln!("start flow calcs: {}", IMAGE0.len());
    rprintln!("buf0: {:?}", &IMAGE0[0..8]);
    rprintln!("buf1: {:?}", &IMAGE1[0..8]);

    const COLS: usize = 64;
    const ROWS: usize = 64;

    let mut correlator = fwht::HadamardCorrelator::new(COLS, ROWS);

    const FRAME_COUNT: u32 = 10;
    let mut last_flow = (0i16, 0i16);
    for _ in 0..FRAME_COUNT {
        last_flow = correlator.calculate_flow(&IMAGE1, &IMAGE0);
        rprintln!("flow: {:?}", last_flow);
    }
    rprintln!("estimated flow: {:?}", last_flow);
    //TODO use eg DWT::get_cycle_count() to measure elapsed time

    rprintln!("<--- DONE --");

    loop {

    }

}

// #[entry]
// fn main() -> ! {
//     rtt_init_print!(BlockIfFull);
//     rprintln!("-- > MAIN --");
//
//     let mut board = Board::default();
//     // this provides the interrupt handler access to the shared Board struct
//     unsafe {
//         BOARD_PTR.store(&mut board, Ordering::SeqCst);
//     }
//
//     let loop_interval = GYRO_REPORTING_INTERVAL_MS as u8;
//     rprintln!("loop_interval: {}", loop_interval);
//
//     let _ = board.activity_led.set_high();
//     let _ = board.comms_led.set_high();
//     let _ = board.error_led.set_high();
//
//     // let mut fast_img_bufs: [_; 2] = unsafe { FAST_IMG_BUFS };
//     let fast_img_bufs: [_; 2] = unsafe { [&mut FAST_IMG0, &mut FAST_IMG1] };
//     // This is how we can enable a grayscale test pattern on the MT9V034
//     // let _ = board.camera_config.as_mut().unwrap().
//     //     enable_pixel_test_pattern(true, PixelTestPattern::DiagonalShade);
//
//     // let mut _mr_fft = micro_rfft::MicroFftContext::new();
//
//     if let Some(dcmi_wrap) = board.dcmi_wrap.as_mut() {
//         dcmi_wrap.enable_capture(&board.dma2);
//     }
//     let mut img_count: u32 = 0;
//     let mut flow_img_idx = 0;
//     loop {
//         for _ in 0..10 {
//             // for _ in 0..10 {
//             //     // read the 6dof frequently
//             //     if let Some(six_dof) = board.gyro.as_mut() {
//             //         if let Ok(_sample) = six_dof.gyro() {
//             //             //rprintln!("gyro {}, {}, {}", _sample.x, _sample.y, _sample.z );
//             //         }
//             //     }
//             // }
//             if let Some(dcmi_wrap) = board.dcmi_wrap.as_mut() {
//                 let dst = fast_img_bufs[flow_img_idx].as_mut();
//                 if let Ok(read_len) = dcmi_wrap.read_available(dst) {
//                     if read_len > 0 {
//                         flow_img_idx = (flow_img_idx + 1) % 2;
//                         // in this example we calculate flow from two most recent frames
//                         // let flow = mr_fft.calculate_flow_fft(
//                         //     &new_frame,
//                         //     &prev_frame,
//                         //     64usize,
//                         //     64usize,
//                         // );
//                         // rprintln!("{} flow: {}, {}", img_count, flow.0, flow.1);
//
//                         let _ = board.activity_led.toggle();
//                         img_count += 1;
//                     }
//                 }
//             }
//             rprintln!("{} old {:?} new {:?}",img_count, &fast_img_bufs[0][0..32], &fast_img_bufs[1][0..32]);
//
//         }
//         let _ = board.comms_led.toggle();
//     }
// }

// /// output image data as 8-bit raw pixels in base64 encoded format, to RTT
// fn dump_pixels(image_count: u32, buf: &[u8]) {
//     rprintln!("\n--- {}", image_count);
//
//     //process input chunks that are multiples of 12 bytes (for base64 continuity)
//     const CHUNK_SIZE: usize = 24;
//     let total_len = buf.len();
//     let mut read_idx = 0;
//     while read_idx < total_len {
//         let max_idx = total_len.min(read_idx + CHUNK_SIZE);
//         let wrapper = Base64Display::with_config(
//             &buf[read_idx..max_idx],
//             base64::STANDARD,
//         );
//         rprint!("{}", wrapper);
//         read_idx += CHUNK_SIZE;
//     }
// }
