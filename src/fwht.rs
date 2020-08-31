/*
Copyright (c) 2020 Todd Stellanova
LICENSE: BSD3 (see LICENSE file)
*/


//! Implement image phase correlation using
//! the Fast Walsh-Hadamard Transform
//! with sign-only correlation
//!
//! References:
//! - "DCT SIGN-ONLY CORRELATION WITH APPLICATION TO IMAGE MATCHING
//! AND THE RELATIONSHIP WITH PHASE-ONLY CORRELATION":
//! Authors: Izumi ITO and Hitoshi KIYA (2007)
//! - "DCT Sign Only Correlation and Its Application to Image Registration":
//! Authors: H. Kiya et al (2006)
//!
use num_integer::Integer;
const NSAMPLES: usize = 64*64;


const BLOCK_DIM: usize = 2;
const BLOCK_LEN: usize = BLOCK_DIM*BLOCK_DIM;
type SadBlock = [u8; BLOCK_LEN];

/// Fast Walsh Hadamard transform correlator:
/// Use FWHT sign-only transformation to measure optical flow
/// via image displacement between two images.
/// Static memory used is (3 * 2) = 6 * NSAMPLES bytes
pub struct HadamardCorrelator {
    cols: usize,
    i16_scratch0: [i16; NSAMPLES],
    i16_scratch1: [i16; NSAMPLES],
    i16_scratch2: [i16; NSAMPLES],
}

impl HadamardCorrelator {
    pub fn new(columns: usize, rows: usize) -> Self {
        let nsamples = columns * rows;
        assert_eq!(nsamples, NSAMPLES, "nsamples restricted to NSAMPLES");

        Self {
            cols: columns,
            i16_scratch0: [0i16; NSAMPLES],
            i16_scratch1: [0i16; NSAMPLES],
            i16_scratch2: [0i16; NSAMPLES],
        }
    }

    /// Fill working buffer from input samples and transform them
    fn fill_and_transform(input: &[u8], output: &mut [i16] )
    {
        super::fill_u8_samples_to_i16(&input, output);
        Self::fwht_transform_i16(output);
        Self::sign_reduce(output);
    }

    /// Estimate translational optical flow using DCT sign-only correlation
    pub fn calculate_flow(
        &mut self,
        new_frame: &[u8],
        old_frame: &[u8],
    ) -> (i16, i16) {
        Self::fill_and_transform(&old_frame, &mut self.i16_scratch0);
        Self::fill_and_transform(&new_frame, &mut self.i16_scratch1);

        Self::hadamard_product_i16(
            &self.i16_scratch0,
            &self.i16_scratch1,
            &mut self.i16_scratch2,
        );
        // scratch2 is now { -1.0, 0.0, 1.0 } for cross-spectral power

        // inverse of FWHT is itself
        Self::fwht_transform_i16(&mut self.i16_scratch2);

        // find the magnitude of dx, dy
        let (max_x, max_y) = Self::find_averaged_peak(&self.i16_scratch2, self.cols);
        // we now know the magnitude of the movement vector but not the direction
        let (sdx, sdy) = Self::sad_block_search(&old_frame, &new_frame,self.cols, max_x, max_y);
        let dx = sdx * (max_x as i16);
        let dy = sdy * (max_y as i16);
        (dx, dy)
    }

    /// Apply Fast Walsh Hadamard Transform in-place to buffer
    fn fwht_transform_i16(buf: &mut [i16]) {
        let mut h = 1;
        while h < buf.len() {
            let h_leap = h << 1;
            for i in (0..buf.len()).step_by(h_leap) {
                for j in i..i+h {
                    let j_leap = j+h;
                    let x = buf[j];
                    let y = buf[j_leap];
                    //println!("j x,y: {} {},{}", j, x, y);
                    buf[j] = x.saturating_add(y);
                    buf[j_leap] = x.saturating_sub(y);
                }
            }
            h = h_leap;
        }
    }

    /// Reduce input values to just their sign
    fn sign_reduce(buf: &mut [i16]) {
        for i in 0..buf.len() {
            buf[i] = match buf[i] {
                n if n > 0 => 1i16,
                _ => -1i16
            };
        }
    }

    /// The sequence of +/- dx, +/- dy used for direction determination
    const SIGN_SEQUENCE: [(i16,i16); 4] = [ (-1, -1), (1, -1), (-1, 1), (1, 1)];

    /// compare a block centered in the new frame to four blocks at:
    /// (-dx, -dy), (-dx, +dy), (+dx, -dy), (+dx, +dy)
    /// determine which has the lowest Sum of Absolute Differences (SAD),
    /// return the corresponding sign of (dx, dy) to make that true
    fn sad_block_search(old_frame: &[u8], new_frame: &[u8], frame_cols: usize, dx: usize, dy: usize)
    -> (i16, i16)
    {
        const BLOCK_HALF_DIM: usize = BLOCK_DIM /2;
        let frame_rows = new_frame.len() / frame_cols;
        let ctr_y = (frame_rows / 2) + (frame_rows % 2) ;
        let ctr_x = (frame_cols / 2) + (frame_cols % 2) ;
        // println!("ctr_x {} ctr_y {} BLOCK_HALF_DIM {}", ctr_x, ctr_y, BLOCK_HALF_DIM);

        assert_eq!(old_frame.len(), new_frame.len(),"old and new frames must be same size");
        let ctr_x_lo = ctr_x - BLOCK_HALF_DIM;
        let ctr_y_lo = ctr_y - BLOCK_HALF_DIM;
        let x_lo = ctr_x_lo - dx;
        let x_hi = ctr_x_lo + dx;
        let y_lo = ctr_y_lo - dy;
        let y_hi = ctr_y_lo + dy;

        // println!("lo {},{} hi {},{}", x_lo, y_lo, x_hi, y_hi);

        let mut search_block = [0u8; BLOCK_LEN];
        Self::fill_block_from_frame(&new_frame, &mut search_block, ctr_x_lo, ctr_y_lo, frame_cols, BLOCK_DIM);
        let mut oldies: [SadBlock; 4] = [[0u8; BLOCK_LEN]; 4];
        // NOTE ths order needs to match SIGN_SEQUENCE
        Self::fill_block_from_frame(&old_frame, &mut oldies[0], x_lo, y_lo, frame_cols, BLOCK_DIM);
        Self::fill_block_from_frame(&old_frame, &mut oldies[1], x_hi, y_lo, frame_cols, BLOCK_DIM);
        Self::fill_block_from_frame(&old_frame, &mut oldies[2], x_lo, y_hi, frame_cols, BLOCK_DIM);
        Self::fill_block_from_frame(&old_frame, &mut oldies[3], x_hi, y_hi, frame_cols, BLOCK_DIM);

        // find the lowest SAD
        let mut min_sad = u16::MAX;
        let mut min_idx: usize = 0;
        for i in 0..Self::SIGN_SEQUENCE.len() {
            let cur_sad = Self::sum_abs_diffs(&search_block, &oldies[i]);
            if cur_sad < min_sad {
                min_sad = cur_sad;
                min_idx = i;
                if cur_sad == 0 { break; }
            }
        }

        Self::SIGN_SEQUENCE[min_idx]
    }

    /// Take a block-size bite out of a larger frame
    /// Currently guaranteed to explode on out-of-bounds
    fn fill_block_from_frame(frame: &[u8],
                             block: &mut [u8],
                             frame_start_x: usize,
                             frame_start_y: usize,
                             frame_cols: usize,
                             block_dim: usize) {
        // let start_x = start_idx % frame_cols;
        // let start_y = start_idx / frame_cols;
        for i in 0..block.len() {
            let block_y = i / block_dim;
            let block_x = i % block_dim;
            let frame_x = frame_start_x + block_x;
            let frame_y = frame_start_y + block_y;
            let frame_idx = frame_y * frame_cols + frame_x;
            block[i] = frame[frame_idx];
        }
    }

    fn sum_abs_diffs(block0: &[u8], block1: &[u8]) -> u16 {
        let mut sum_diff: u16 = 0;
        for i in 0..block0.len() {
            sum_diff += block0[i].wrapping_sub(block1[i]) as u16;
        }
        //println!("diff {:?} {:?} -> {}", block0, block1, sum_diff);
        sum_diff
    }

    /// Multiply two arrays element-wise
    fn hadamard_product_i16(f0: &[i16], f1: &[i16], result: &mut [i16]) {
        //TODO could use SIMD here ?
        for i in 0..f0.len() {
            // NOTE number of zeros seems vanishingly small
            result[i] = f0[i] * f1[i]
        }
    }


    /// Find the average peak of top two
    fn find_averaged_peak(
        result: &[i16],
        columns: usize,
    ) -> (usize, usize) {
        let peaks = Self::find_top_two_peaks(&result, columns);
        //println!("peaks: {:?}", peaks);

        let mut avg_x = 0;
        let mut avg_y = 0;
        for i in 0..2 {
            let peak = peaks[i];
            avg_x += peak.0;
            avg_y += peak.1;
        }
        avg_x = avg_x.div_ceil(&2);
        avg_y = avg_y.div_ceil(&2);
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



#[cfg(test)]
mod tests {
    use super::*;

    const FRAME_25_DIM: usize = 5;
    const FRAME_25: [u8; FRAME_25_DIM * FRAME_25_DIM] = [
        10, 20, 30, 40, 50,
        11, 21, 31, 41, 51,
        12, 22, 32, 42, 52,
        13, 23, 33, 43, 53,
        14, 24, 34, 44, 54 ];

    const FRAME_64_DIM: usize = 8;
    const FRAME_64: [u8; FRAME_64_DIM * FRAME_64_DIM] = [
        10, 20, 30, 40, 50, 60, 70, 80,
        11, 21, 31, 41, 51, 61, 71, 81,
        12, 22, 32, 42, 52, 62, 72, 82,
        13, 23, 33, 43, 53, 63, 73, 83,
        14, 24, 34, 44, 54, 64, 74, 84,
        15, 25, 35, 45, 55, 65, 75, 85,
        16, 26, 36, 46, 56, 66, 76, 86,
        17, 27, 37, 47, 57, 67, 77, 87 ];

    const FRAME_64_SHIFT_P1P1: [u8; FRAME_64_DIM * FRAME_64_DIM] = [
        21, 31, 41, 51, 61, 71, 81, 91,
        22, 32, 42, 52, 62, 72, 82, 92,
        23, 33, 43, 53, 63, 73, 83, 93,
        24, 34, 44, 54, 64, 74, 84, 94,
        25, 35, 45, 55, 65, 75, 85, 95,
        26, 36, 46, 56, 66, 76, 86, 96,
        27, 37, 47, 57, 67, 77, 87, 97,
        28, 38, 48, 58, 68, 78, 88, 98 ];


    const FRAME_64_SHIFT_N1_P1: [u8; FRAME_64_DIM * FRAME_64_DIM] = [
        1, 11, 21, 31, 41, 51, 61, 71,
        2, 12, 22, 32, 42, 52, 62, 72,
        3, 13, 23, 33, 43, 53, 63, 73,
        4, 14, 24, 34, 44, 54, 64, 74,
        5, 15, 25, 35, 45, 55, 65, 75,
        6, 16, 26, 36, 46, 56, 66, 76,
        7, 17, 27, 37, 47, 57, 67, 77,
        8, 18, 28, 38, 48, 58, 68, 78 ];


    #[test]
    fn block_loading() {
        let mut block = [0u8; BLOCK_LEN];
        let frame = FRAME_25;
        const FRAME_COLS: usize = FRAME_25_DIM;
        let start_x: usize = 2;
        let start_y: usize = 1;
        HadamardCorrelator::fill_block_from_frame(
            &frame,
            &mut block,
            start_x, start_y,
            FRAME_COLS, BLOCK_DIM);

        assert_eq!(block, [31,41,32,42]);
    }


    #[test]
    fn sad_search() {
        let (sdx, sdy) =
            HadamardCorrelator::sad_block_search(
            &FRAME_64, &FRAME_64_SHIFT_P1P1, FRAME_64_DIM,1, 1);

        assert_eq!(sdx, 1i16, "sign dx should be P1");
        assert_eq!(sdy, 1i16, "sign dy should be P1");


        let (sdx, sdy) =
            HadamardCorrelator::sad_block_search(
                &FRAME_64, &FRAME_64_SHIFT_N1_P1, FRAME_64_DIM,1, 1);
        assert_eq!(sdx, -1i16, "sign dx should be N1");
        assert_eq!(sdy, 1i16, "sign dy should be P1");
    }

}