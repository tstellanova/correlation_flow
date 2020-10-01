/*
Copyright (c) 2020 Todd Stellanova
LICENSE: BSD3 (see LICENSE file)
*/
#![no_std]

//! Calculates 2D image translation using image correlation
//!
//! Could be used for:
//! - stitching image panoramas together (by calculating the point where two images overlap)
//! - image registration
//! - measuring optical flow (by measuring the 2D translation between image frames)
//!
//! This library to operate in very low memory complexity on no_std rust
//! for embedded systems.
//!

use num_complex::Complex32;
pub mod micro_rfft;


/// Convert 8-bit mono image data to normalized f32
pub(crate) fn fill_u8_samples_to_f32(input: &[u8], output: &mut [f32]) {
    let n = input.len();
    for i in 0..n {
        output[i] = (input[i] as f32) / 255.0;
    }
}

pub(crate) fn find_peak_real(
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
