/*
Copyright (c) 2020 Todd Stellanova
LICENSE: BSD3 (see LICENSE file)
*/

use num_traits::Zero;
use microfft::Complex32;
// use num_complex::Complex32;
use micromath::F32Ext;

const NSAMPLES: usize = 64*64;
type RealSampleBuf = [f32; NSAMPLES];

/// Performs flow correlation using (2*4 + 8) = 16 * NSAMPLES bytes
pub struct MicroFftContext {
    rscratch0: RealSampleBuf,
    rscratch1: RealSampleBuf,
    cscratch0: [Complex32; NSAMPLES],
}

impl MicroFftContext {

    pub fn new() -> Self {
        Self {
            rscratch0: [0f32; NSAMPLES],
            rscratch1: [0f32; NSAMPLES],
            cscratch0: [Complex32::zero(); NSAMPLES],
        }
    }

    /// Calculate the cross-power spectrum of two sets of samples
    /// Normalize the result.
    /// Output in self.cscratch0
    fn cross_power_spectrum_norm(
        &mut self,
        buf0: &[u8],
        buf1: &[u8],
    )
    {
        super::fill_u8_samples_to_f32(&buf0, &mut self.rscratch0);
        super::fill_u8_samples_to_f32(&buf1, &mut self.rscratch1);

        //rfft_4096() doesn't allocate Complex32 buffers-- it reuses the provided scratch buffers
        let fft0 = microfft::real::rfft_4096(&mut self.rscratch0);
        let mut fft1 = microfft::real::rfft_4096(&mut self.rscratch1);

        // f0 * f1' / | f0 * f1'|
        Self::complex_conjugate(&mut fft1);
        Self::clear_complex(&mut self.cscratch0);
        Self::complex_hadamard_product_norm(&fft0, &fft1, &mut self.cscratch0);
    }

    /// Measure optical flow by calculating the displacement between
    /// and old image frame and a new image frame.
    /// - Input is a pair of 8-bit grayscale images
    /// - Output is a (dx, dy) estimate of the pixel displacement between images
    ///
    pub fn calculate_flow_fft(
        &mut self,
        new_frame: &[u8],
        old_frame: &[u8],
        columns: usize,
        rows: usize,
    ) -> (i16, i16) {
        let nsamples = columns * rows;
        let x_limit = columns / 2;
        let y_limit = rows / 2;
        assert_eq!(nsamples, NSAMPLES, "nsamples restricted to 4096");
        self.cross_power_spectrum_norm(&old_frame, &new_frame);
        Self::invert_fft(&mut self.cscratch0);

        let (max_x, max_y, _max_val) =
            super::find_peak_real(&self.cscratch0, columns);
        //println!("{} {} {}", max_x, max_y, _max_val);

        // calculate wrapping as negative displacement:
        let dx: i16 = if max_x > x_limit { max_x as i16 - columns as i16 } else { max_x as i16 };
        let dy: i16 = if max_y > y_limit { max_y as i16 - rows as i16 } else { max_y as i16 };
        (dx, dy)
    }

    /// Zero out all complex values
    fn clear_complex(data: &mut [Complex32]) {
        for i in 0..data.len() {
            data[i].set_zero();
        }
    }

    /// Perform an inverse FFT in-place on the data provided
    fn invert_fft(data: &mut [Complex32]) {
        //swap real and imaginary components
        Self::swap_complex_components(data);

        //now forward-transform
        microfft::complex::cfft_4096(data);

        //reswap real and imaginary components
        Self::swap_complex_components(data);
    }

    /// Swap real and imaginary components
    fn swap_complex_components(data: &mut [Complex32]) {
        for i in 0..data.len() {
            let tmp = data[i].re;
            data[i].re = data[i].im;
            data[i].im = tmp;
        }
    }

    /// Conjugate all complex values in the sequence
    fn complex_conjugate(data: &mut [Complex32]) {
        for i in 0..data.len() {
            data[i] = data[i].conj();
        }
    }


    /// Normalized complex Hadamard product of two arrays
    fn complex_hadamard_product_norm(
        f0: &[Complex32],
        f1: &[Complex32],
        result: &mut [Complex32]
    )
    {
        let n = f0.len();
        for i in 0..n {
            let prod = f0[i] * f1[i];
            let mag = prod.re.hypot(prod.im);//TODO use prod.norm when available?
            //let mag = prod.norm_sqr().sqrt();
            result[i] = prod / mag;
        }

    }

}
