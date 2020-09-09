/*
Copyright (c) 2020 Todd Stellanova
LICENSE: BSD3 (see LICENSE file)
*/

//! Implements image phase correlation using the Fast Walsh-Hadamard Transform
//! with sign-only correlation. Could be used to:
//! - stitch images together (by calculating the point where images overlap)
//! - measure optical flow (by measuring the 2D translation between images)
//!
//! References:
//! - Izumi ITO and Hitoshi KIYA, "DCT SIGN-ONLY CORRELATION WITH APPLICATION TO IMAGE MATCHING
//! AND THE RELATIONSHIP WITH PHASE-ONLY CORRELATION" [DOI 10.1109/ICASSP.2007.366138](https://doi.org/10.1109/ICASSP.2007.366138)
//! Authors: Izumi ITO and Hitoshi KIYA (2007)
//! - Hitoshi KIYA et al, "DCT Sign Only Correlation and Its Application to Image Registration"
//! [DOI 10.1109/APCCAS.2006.342490](https://doi.org/10.1109/APCCAS.2006.342490)
//!
//!

/// Currently the number of allowed samples is statically fixed at NSAMPLES
const FRAME_DIM: usize = 64;
pub const NSAMPLES: usize = FRAME_DIM * FRAME_DIM;

/// Length of an edge of the SAD blocks
const SAD_BLOCK_DIM: usize = 2;
const SAD_BLOCK_LEN: usize = SAD_BLOCK_DIM * SAD_BLOCK_DIM;


// use generic_array::typenum::U5;
// use generic_array::{ArrayLength, GenericArray};

/// Buffer size for processing a full frame
type FwhtFrameBufI16 = [i16; NSAMPLES];

/// Fast Walsh Hadamard transform correlator:
/// Uses FWHT sign-only transformation to measure displacement vector between two images.
/// Can be used to measure optical flow by comparing a series of images.
///
/// Minimum static memory used is (3 * 2) = 6 * NSAMPLES bytes
pub struct HadamardCorrelator {
    /// The number of columns in the image.
    /// This and NSAMPLES defines the "shape" of the image matrix.
    cols: usize,
    rows: usize,
    /// the center x (column) of a frame
    frame_center_x: usize,
    /// the center y (row) of a frame
    frame_center_y: usize,

    /// Scratch buffers used for processing
    i16_scratch0: FwhtFrameBufI16,
    i16_scratch1: FwhtFrameBufI16,
    i16_scratch2: FwhtFrameBufI16,
}

impl HadamardCorrelator {
    pub fn new(columns: usize, rows: usize) -> Self {
        let nsamples = columns * rows;
        assert_eq!(nsamples, NSAMPLES, "size restricted to NSAMPLES");
        let frame_center_x = (columns / 2) + (columns % 2);
        let frame_center_y = (rows / 2) + (rows % 2);

        Self {
            cols: columns,
            rows,
            frame_center_x,
            frame_center_y,
            i16_scratch0: [0i16; NSAMPLES],
            i16_scratch1: [0i16; NSAMPLES],
            i16_scratch2: [0i16; NSAMPLES],
        }
    }

    /// Fill working buffer from input samples,
    /// then transform, then sign-reduce
    fn fill_transform_reduce(input: &[u8], output: &mut [i16]) {
        Self::fill_u8_samples_to_i16(&input, output);
        // Self::fwht_transform_1d_i16(output);
        Self::fwht_transform_2d_i16(output, FRAME_DIM);
        Self::sign_reduce(output);
    }

    /// Copy u8 sample bytes to i16 working buffer
    fn fill_u8_samples_to_i16(input: &[u8], output: &mut [i16]) {
        let n = input.len().min(output.len());
        for i in 0..n {
            output[i] = input[i] as i16;
        }
    }

    /// Fast filter to eliminate low-change frames
    fn center_block_change_check(
        &self,
        new_frame: &[u8],
        old_frame: &[u8],
        sad_threshold: u16,
    ) -> bool {
        let mut subframe1 = [0u8; 16];
        let mut subframe0 = [0u8; 16];
        Self::fill_block_from_frame(
            &old_frame,
            &mut subframe0,
            self.frame_center_x - 2,
            self.frame_center_y - 2,
            self.cols,
            4,
        );
        Self::fill_block_from_frame(
            &new_frame,
            &mut subframe1,
            self.frame_center_x - 2,
            self.frame_center_y - 2,
            self.cols,
            4,
        );
        let quick_sad = Self::sum_abs_diffs(&subframe0, &subframe1);
        if quick_sad < (sad_threshold * (subframe0.len() as u16)) {
            return false;
        }
        true
    }

    /// Measure 2D translation, the movement of image or camera
    /// between two image frames.  This is the basis of optical flow measurement.
    /// - Returns (dx, dy)
    /// - `old_frame` is an image sample frame that was obtained earlier in time
    /// - `new_frame` is an image sample frame that was obtained later in time
    pub fn measure_translation(
        &mut self,
        new_frame: &[u8],
        old_frame: &[u8],
    ) -> (i16, i16) {
        let max_col_idx = self.frame_center_x;
        let max_row_idx = self.frame_center_y;
        if !self.center_block_change_check(new_frame, old_frame, 2) {
            return (0, 0);
        }

        Self::fill_transform_reduce(&old_frame, &mut self.i16_scratch0);
        Self::fill_transform_reduce(&new_frame, &mut self.i16_scratch1);
        // scratch0 and scratch1 now contain sets of -1/+1 which are the sign-reduced
        // form of their FWHT-transformed values

        Self::hadamard_product_i16(
            &self.i16_scratch0,
            &self.i16_scratch1,
            &mut self.i16_scratch2,
        );
        // scratch2 is now { -1.0, 0.0, 1.0 } for cross-spectral power

        // inverse of FWHT is itself
        Self::fwht_transform_1d_i16(&mut self.i16_scratch2);

        // find the magnitude of dx, dy
        let (max_x, max_y) =
            Self::find_averaged_peak(&self.i16_scratch2, self.cols);
        if max_x > max_col_idx || max_y > max_row_idx {
            // bogus motion
            return (0, 0);
        }
        // we now know the magnitude of the movement vector but not the direction
        let (sdx, sdy) = Self::sad_block_search(
            &old_frame, &new_frame, self.cols, max_x, max_y,
        );
        let dx = sdx * (max_x as i16);
        let dy = sdy * (max_y as i16);
        (dx, dy)
    }

    /// Apply Fast Walsh Hadamard Transform (FWHT) in-place to buffer
    fn fwht_transform_1d_i16(buf: &mut [i16]) {
        let mut h = 1;
        while h < buf.len() {
            let h_leap = h << 1;
            for i in (0..buf.len()).step_by(h_leap) {
                for j in i..i + h {
                    let j_leap = j + h;
                    let x = buf[j];
                    let y = buf[j_leap];
                    //println!("j x,y: {} {},{}", j, x, y);
                    // we use saturating arithmetic here because, for our application,
                    // this is sufficient to search for peaks
                    buf[j] = x.saturating_add(y);
                    buf[j_leap] = x.saturating_sub(y);
                }
            }
            h = h_leap;
        }
    }

    fn fwht_transform_2d_i16(buf: &mut [i16], ncols: usize) {
        let nrows = buf.len() / ncols;

        //transform rows first
        for row in 0..nrows {
            let idx = row*ncols;
            Self::fwht_transform_1d_i16(&mut buf[idx..idx+ncols]);
        }

        let mut col_buf = [0i16; FRAME_DIM]; //nrows
        //  transform columns
        // TODO use iterator instead?
        for col in 0..ncols {
            // fill column buffer
            for row in 0..nrows {
                col_buf[row] = buf[row*ncols + col];
            }
            Self::fwht_transform_1d_i16(&mut col_buf);
            for row in 0..nrows {
                //copy back
                buf[row*ncols + col] = col_buf[row];
            }
        }

    }


        /// Reduce input values to just their sign: +1 or -1
    /// For our application we disallow zero
    fn sign_reduce(buf: &mut [i16]) {
        for i in 0..buf.len() {
            buf[i] = match buf[i] {
                n if n > 0 => 1,
                _ => -1,
            };
        }
    }

    /// Find the SAD between the search block and a block from the old_frame
    fn calc_sad_at_offset(
        old_frame: &[u8],
        search_block: &[u8],
        col: usize,
        row: usize,
        frame_cols: usize,
    ) -> u16 {
        let mut test_block = [0u8; SAD_BLOCK_LEN];
        Self::fill_block_from_frame(
            &old_frame,
            &mut test_block,
            col,
            row,
            frame_cols,
            SAD_BLOCK_DIM,
        );
        Self::sum_abs_diffs(&search_block, &test_block)
    }

    /// The sequence of +/- dx, +/- dy used for direction determination
    const QUADRANT_SIGN_SEQUENCE: [(i16, i16); 4] =
        [(-1, -1), (1, -1), (-1, 1), (1, 1)];

    /// - Compare a sample block centered in the new frame to four blocks at:
    /// (-dx, -dy), (+dx, -dy), (-dx, +dy), (+dx, +dy) -- see `QUADRANT_SIGN_SEQUENCE`
    /// - Determine which comparison has the lowest Sum of Absolute Differences (SAD),
    /// - Returns the corresponding sign of (dx, dy)
    /// - `old_frame` and `new_frame` are a sequence of equal-sized image frames
    /// - `frame_cols` is the number of columns per row in an image frame
    /// - (`dx`, `dy`) comprise the previously measured magnitude of the translation between frames
    fn sad_block_search(
        old_frame: &[u8],
        new_frame: &[u8],
        frame_cols: usize,
        dx: usize,
        dy: usize,
    ) -> (i16, i16) {
        const SAD_BLOCK_HALF_DIM: usize = SAD_BLOCK_DIM / 2;
        assert_eq!(
            old_frame.len(),
            new_frame.len(),
            "old and new frames must be same size"
        );
        let frame_rows = new_frame.len() / frame_cols;
        let max_col = frame_cols - SAD_BLOCK_HALF_DIM;
        let max_row = frame_rows - SAD_BLOCK_HALF_DIM;
        let ctr_y = (frame_rows / 2) + (frame_rows % 2);
        let ctr_x = (frame_cols / 2) + (frame_cols % 2);

        // ctr_x_lo , ctr_y_lo are the upper-left corner of the new frame sample block
        let ctr_x_lo = ctr_x - SAD_BLOCK_HALF_DIM;
        let ctr_y_lo = ctr_y - SAD_BLOCK_HALF_DIM;
        // next we calculate the edges of the blocks to compare with (from the old frame)
        let x_lo = ctr_x_lo.wrapping_sub(dx); // we use wrap to detect out-of-bounds
        let x_hi = ctr_x_lo + dx; // we use > max_col to detect out-of-bounds
        let y_lo = ctr_y_lo.wrapping_sub(dy); // we use wrap to detect out-of-bounds
        let y_hi = ctr_y_lo + dy; // we use  > max_row to detect out-of-bounds

        let mut search_block = [0u8; SAD_BLOCK_LEN];
        Self::fill_block_from_frame(
            &new_frame,
            &mut search_block,
            ctr_x_lo,
            ctr_y_lo,
            frame_cols,
            SAD_BLOCK_DIM,
        );

        // find the lowest Sum of Absolute Differences between search block and surrounding blocks
        let mut min_sad = u16::MAX;
        let mut min_idx: usize = 0;
        // NOTE ths order needs to match QUADRANT_SIGN_SEQUENCE:
        if x_lo < ctr_x_lo && y_lo < ctr_y_lo {
            let cur_sad = Self::calc_sad_at_offset(
                &old_frame,
                &search_block,
                x_lo,
                y_lo,
                frame_cols,
            );
            if cur_sad < min_sad {
                min_sad = cur_sad;
                min_idx = 0;
                // break early if SAD is zero (maximum exact match)
                if cur_sad == 0 {
                    return Self::QUADRANT_SIGN_SEQUENCE[min_idx];
                }
            }
        }

        if x_hi < max_col && y_lo < ctr_y_lo {
            let cur_sad = Self::calc_sad_at_offset(
                &old_frame,
                &search_block,
                x_hi,
                y_lo,
                frame_cols,
            );
            if cur_sad < min_sad {
                min_sad = cur_sad;
                min_idx = 1;
                // break early if SAD is zero (maximum exact match)
                if cur_sad == 0 {
                    return Self::QUADRANT_SIGN_SEQUENCE[min_idx];
                }
            }
        }

        if x_lo < ctr_x_lo && y_hi < max_row {
            let cur_sad = Self::calc_sad_at_offset(
                &old_frame,
                &search_block,
                x_lo,
                y_hi,
                frame_cols,
            );
            if cur_sad < min_sad {
                min_sad = cur_sad;
                min_idx = 2;
                // break early if SAD is zero (maximum exact match)
                if cur_sad == 0 {
                    return Self::QUADRANT_SIGN_SEQUENCE[min_idx];
                }
            }
        }

        if x_hi < max_col && y_hi < max_row {
            let cur_sad = Self::calc_sad_at_offset(
                &old_frame,
                &search_block,
                x_hi,
                y_hi,
                frame_cols,
            );
            if cur_sad < min_sad {
                // min_sad = cur_sad;
                min_idx = 3;
                // break early if SAD is zero (maximum exact match)
                if cur_sad == 0 {
                    return Self::QUADRANT_SIGN_SEQUENCE[min_idx];
                }
            }
        }

        Self::QUADRANT_SIGN_SEQUENCE[min_idx]
    }

    /// Take a block-size bite out of a larger image sample frame.
    /// Currently guaranteed to explode on out-of-bounds:
    /// relies on frames being much bigger than blocks.
    pub fn fill_block_from_frame(
        frame: &[u8],
        block: &mut [u8],
        frame_start_x: usize,
        frame_start_y: usize,
        frame_cols: usize,
        block_dim: usize,
    ) {
        let frame_len = frame.len();
        //memcpy one row at a time
        for block_row in 0..block_dim {
            let block_idx = block_row * block_dim;
            let frame_y = frame_start_y + block_row;
            let frame_idx = frame_y * frame_cols + frame_start_x;
            if (frame_idx + block_dim) < frame_len {
                block[block_idx..block_idx + block_dim]
                    .copy_from_slice(&frame[frame_idx..frame_idx + block_dim]);
            }
        }
    }

    /// Sum of Absolute Differences between two sample blocks
    fn sum_abs_diffs(block0: &[u8], block1: &[u8]) -> u16 {
        let mut sum_diff: u16 = 0;
        for i in 0..block0.len() {
            let val_b = block1[i];
            let val_a = block0[i];
            //TODO someday use eg core::arch::arm::vsubq_u32()
            sum_diff += (if val_a > val_b {
                val_a - val_b
            } else {
                val_b - val_a
            }) as u16;
        }
        sum_diff
    }

    /// Multiply two arrays element-wise ("Hadamard product")
    fn hadamard_product_i16(f0: &[i16], f1: &[i16], result: &mut [i16]) {
        for i in 0..f0.len() {
            result[i] = f0[i] * f1[i]
        }
    }

    /// Find the average position of top two values in the given slice.
    /// - `columns` is the number of columns per row in the input data
    /// - Returns (dx, dy) of the peak position in the image
    fn find_averaged_peak(input: &[i16], columns: usize) -> (usize, usize) {
        let peaks = Self::find_top_two_peaks(&input, columns);
        //println!("peaks: {:?}", peaks);

        let mut avg_x = 0;
        let mut avg_y = 0;
        for i in 0..2 {
            let peak = peaks[i];
            avg_x += peak.0;
            avg_y += peak.1;
        }
        avg_x = avg_x / 2;
        avg_y = avg_y / 2;
        (avg_x, avg_y)
    }

    /// Find the top two peaks in an cols x rows matrix
    fn find_top_two_peaks(
        input: &[i16],
        cols: usize,
    ) -> [(usize, usize, i16); 2] {
        let mut working = [(0, 0i16); 2];
        let sample_count = input.len();
        for i in 1..sample_count {
            let val = input[i];
            if val > working[0].1 {
                if val > working[1].1 {
                    // swap max into min peak
                    working[0] = working[1];
                    working[1] = (i, val);
                } else {
                    //replace min peak
                    working[0] = (i, val);
                }
            }
        }

        let min_peak = (working[0].0 % cols, working[0].0 / cols);
        let max_peak = (working[1].0 % cols, working[1].0 / cols);
        [
            (min_peak.0, min_peak.1, working[0].1),
            (max_peak.0, max_peak.1, working[1].1),
        ]
    }
}

// Note these tests only run on non-embedded by commenting out dev-dependencies in Cargo.toml
#[cfg(test)]
mod tests {
    use super::*;

    const FRAME_25_DIM: usize = 5;
    #[rustfmt::skip]
    const FRAME_25: [u8; FRAME_25_DIM * FRAME_25_DIM] = [
        10, 20, 30, 40, 50,
        11, 21, 31, 41, 51,
        12, 22, 32, 42, 52,
        13, 23, 33, 43, 53,
        14, 24, 34, 44, 54,
    ];

    const FRAME_64_DIM: usize = 8;
    #[rustfmt::skip]
    const FRAME_64: [u8; FRAME_64_DIM * FRAME_64_DIM] = [
        10, 20, 30, 40, 50, 60, 70, 80,
        11, 21, 31, 41, 51, 61, 71, 81,
        12, 22, 32, 42, 52, 62, 72, 82,
        13, 23, 33, 43, 53, 63, 73, 83,
        14, 24, 34, 44, 54, 64, 74, 84,
        15, 25, 35, 45, 55, 65, 75, 85,
        16, 26, 36, 46, 56, 66, 76, 86,
        17, 27, 37, 47, 57, 67, 77, 87,
    ];
    #[rustfmt::skip]
    const FRAME_64_SHIFT_P1P1: [u8; FRAME_64_DIM * FRAME_64_DIM] = [
        21, 31, 41, 51, 61, 71, 81, 91,
        22, 32, 42, 52, 62, 72, 82, 92,
        23, 33, 43, 53, 63, 73, 83, 93,
        24, 34, 44, 54, 64, 74, 84, 94,
        25, 35, 45, 55, 65, 75, 85, 95,
        26, 36, 46, 56, 66, 76, 86, 96,
        27, 37, 47, 57, 67, 77, 87, 97,
        28, 38, 48, 58, 68, 78, 88, 98,
    ];
    #[rustfmt::skip]
    const FRAME_64_SHIFT_N1_P1: [u8; FRAME_64_DIM * FRAME_64_DIM] = [
        1, 11, 21, 31, 41, 51, 61, 71,
        2, 12, 22, 32, 42, 52, 62, 72,
        3, 13, 23, 33, 43, 53, 63, 73,
        4, 14, 24, 34, 44, 54, 64, 74,
        5, 15, 25, 35, 45, 55, 65, 75,
        6, 16, 26, 36, 46, 56, 66, 76,
        7, 17, 27, 37, 47, 57, 67, 77,
        8, 18, 28, 38, 48, 58, 68, 78,
    ];

    #[test]
    fn block_loading() {
        let mut block = [0u8; SAD_BLOCK_LEN];
        let frame = FRAME_25;
        const FRAME_COLS: usize = FRAME_25_DIM;
        let start_x: usize = 2;
        let start_y: usize = 1;
        HadamardCorrelator::fill_block_from_frame(
            &frame,
            &mut block,
            start_x,
            start_y,
            FRAME_COLS,
            SAD_BLOCK_DIM,
        );

        assert_eq!(block, [31, 41, 32, 42]);
    }

    #[test]
    fn sad_search() {
        let (sdx, sdy) = HadamardCorrelator::sad_block_search(
            &FRAME_64,
            &FRAME_64_SHIFT_P1P1,
            FRAME_64_DIM,
            1,
            1,
        );

        assert_eq!(sdx, 1i16, "sign dx should be P1");
        assert_eq!(sdy, 1i16, "sign dy should be P1");

        let (sdx, sdy) = HadamardCorrelator::sad_block_search(
            &FRAME_64,
            &FRAME_64_SHIFT_N1_P1,
            FRAME_64_DIM,
            1,
            1,
        );
        assert_eq!(sdx, -1i16, "sign dx should be N1");
        assert_eq!(sdy, 1i16, "sign dy should be P1");
    }
}
