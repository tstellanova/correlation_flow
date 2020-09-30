/*
Copyright (c) 2020 Todd Stellanova
LICENSE: BSD3 (see LICENSE file)
*/

//! Image phase correlation using the embedded-friendly microfft crate.

use microfft::{complex::cfft_4096, real::rfft_4096};
use num_complex::Complex32;
use num_traits::Zero;
use num_traits::real::Real;

/// The number of columns supported
pub const COL_DIM: usize = 64;
/// The number of rows supported
pub const ROW_DIM: usize = 64;

pub const BLOCK_LEN: usize = COL_DIM * ROW_DIM;
type BlockBuffer = [f32; BLOCK_LEN];

/// Performs flow correlation using (2*4 + 8) = 16 * BLOCK_LEN bytes
pub struct MicroFftContext {
    rscratch0: BlockBuffer,
    rscratch1: BlockBuffer,
    cscratch0: [Complex32; BLOCK_LEN],
}

impl MicroFftContext {
    pub fn new() -> Self {
        Self {
            rscratch0: [0f32; BLOCK_LEN],
            rscratch1: [0f32; BLOCK_LEN],
            cscratch0: [Complex32::zero(); BLOCK_LEN],
        }
    }

    /// Calculate the translation (flow) between two 8-bit mono image frames
    /// - old_frame and new_frame image frames with pixels in row-major-order.
    /// - Returns (x, y) of translation
    pub fn measure_translation(
        &mut self,
        new_frame: &[u8],
        old_frame: &[u8],
    ) -> (i16, i16) {
        self.calculate_flow_fft(new_frame, old_frame, COL_DIM, ROW_DIM)
    }

    /// Output in self.cscratch0
    fn cross_power_spectrum_norm(&mut self, buf0: &[u8], buf1: &[u8]) {
        super::fill_u8_samples_to_f32(&buf0, &mut self.rscratch0);
        super::fill_u8_samples_to_f32(&buf1, &mut self.rscratch1);

        //rfft_4096() doesn't allocate Complex32 buffers-- it reuses the provided scratch buffers
        let fft0 = rfft_4096(&mut self.rscratch0);
        let mut fft1 = rfft_4096(&mut self.rscratch1);

        // f0 * f1' / | f0 * f1'|
        Self::complex_conjugate(&mut fft1);
        Self::clear_complex(&mut self.cscratch0);
        Self::complex_hadamard_product_norm(&fft0, &fft1, &mut self.cscratch0);
    }

    ///
    /// Calculate the translation (flow) between two 8-bit mono image frames
    /// - old_frame and new_frame image frames with pixels in row-major-order.
    /// - Returns (x, y) of translation
    ///
    pub fn calculate_flow_fft(
        &mut self,
        new_frame: &[u8],
        old_frame: &[u8],
        columns: usize,
        rows: usize,
    ) -> (i16, i16) {
        assert_eq!(columns, ROW_DIM, "Only {} rows supported", COL_DIM);
        assert_eq!(rows, ROW_DIM, "Only {} rows supported", ROW_DIM);
        let nsamples = columns * rows;
        let x_limit = columns / 2;
        let y_limit = rows / 2;
        assert_eq!(nsamples, BLOCK_LEN, "nsamples restricted to {}", BLOCK_LEN);

        self.cross_power_spectrum_norm(&old_frame, &new_frame);
        Self::invert_fft(&mut self.cscratch0);

        let (max_x, max_y, _max_val) =
            super::find_peak_real(&self.cscratch0, columns);

        // calculate wrapping as negative displacement:
        let dx: i16 = if max_x > x_limit {
            max_x as i16 - columns as i16
        } else {
            max_x as i16
        };
        let dy: i16 = if max_y > y_limit {
            max_y as i16 - rows as i16
        } else {
            max_y as i16
        };
        (dx, dy)
    }

    /// zero out the complex slice
    fn clear_complex(data: &mut [Complex32]) {
        for i in 0..data.len() {
            data[i].set_zero();
        }
    }

    /// Inverse FFT
    fn invert_fft(data: &mut [Complex32]) {
        Self::neg_im_components(data);

        //now forward-transform
        cfft_4096(data);

        Self::neg_im_components(data);
    }

    /// negate the imaginary component of each element of the slice
    fn neg_im_components(data: &mut [Complex32]) {
        for i in 0..data.len() {
            data[i].im = -data[i].im;
        }
    }

    /// Perform complex conjugate on the entire slice
    fn complex_conjugate(data: &mut [Complex32]) {
        for i in 0..data.len() {
            data[i] = data[i].conj();
        }
    }

    /// Normalized complex Hadamard product
    fn complex_hadamard_product_norm(
        f0: &[Complex32],
        f1: &[Complex32],
        result: &mut [Complex32],
    ) {
        let n = f0.len();
        for i in 0..n {
            let prod = f0[i] * f1[i];
            //let prod_norm = (prod.re*prod*re + prod.im*prod.im).sqrt();
            //let prod_norm = prod.re.hypot(prod.im);
            let prod_norm = prod.norm_sqr().sqrt();
            result[i] = prod / prod_norm;
        }
    }
}
