/*
Copyright (c) 2020 Todd Stellanova
LICENSE: BSD3 (see LICENSE file)
*/
#![no_std]

use microfft::Complex32;
// use num_complex::{Complex32};

#[cfg(feature = "complex_fft")]
pub mod micro_cfft;
#[cfg(feature = "real_fft")]
pub mod micro_rfft;
#[cfg(feature = "fwht")]
pub mod fwht;

/// Convert 8-bit sample data to normalized f32
fn fill_u8_samples_to_f32(input: &[u8], output: &mut [f32]) {
    let n = input.len();
    for i in 0..n {
        output[i] = (input[i] as f32) / 255.0;
    }
}

pub fn fill_u8_samples_to_i16(input: &[u8], output: &mut [i16]) {
    let n = input.len();
    for i in 0..n {
        output[i] = input[i] as i16;
    }
}

/// Find the peak real value in an array of Complex32
/// - returns the (x,y) position of the peak
fn find_peak_real(
    input: &[Complex32],
    cols: usize,
) -> (usize, usize, f32)
{
    let mut peak_val = 0f32;
    let mut peak_idx: usize = 0;
    for i in 0..input.len() {
        if input[i].re > peak_val {
            peak_idx = i;
            peak_val = input[i].re;
        }
    }

    let x = peak_idx % cols;
    let y = peak_idx / cols;
    (x, y, peak_val)
}