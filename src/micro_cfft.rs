/*
Copyright (c) 2020 Todd Stellanova
LICENSE: BSD3 (see LICENSE file)
*/

use num_traits::Zero;
use microfft::Complex32;
use micromath::F32Ext;

const NSAMPLES: usize = 64*64;
type SampleBuf = [Complex32; NSAMPLES];

/// Performs flow correlation using (3*8) = 24 * NSAMPLES bytes
pub struct MicroFftContext {
    scratch0: SampleBuf,
    scratch1: SampleBuf,
    scratch2: SampleBuf,
}

impl MicroFftContext {

    pub fn new() -> Self {
        Self {
            scratch0: [Complex32::zero(); NSAMPLES],
            scratch1: [Complex32::zero(); NSAMPLES],
            scratch2: [Complex32::zero(); NSAMPLES],
        }
    }

    pub fn cross_power_spectrum_norm(
        &mut self,
        buf0: &[u8],
        buf1: &[u8],
    )
    {
        Self::fill_u8_samples_to_complex(&buf0, &mut self.scratch0);
        Self::fill_u8_samples_to_complex(&buf1, &mut self.scratch1);

        microfft::complex::cfft_4096(&mut self.scratch0);
        microfft::complex::cfft_4096(&mut self.scratch1);

        // f0 * f1' / | f0 * f1'|
        Self::complex_conjugate(&mut self.scratch1);
        Self::complex_hadamard_product_norm(&self.scratch0, &self.scratch1, &mut self.scratch2);
    }

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
        Self::invert_fft(&mut self.scratch2);
        let (max_x, max_y, _max_val) =
            super::find_peak_real(&self.scratch2, columns);
        //println!("{} {} {}", max_x, max_y, _max_val);

        // calculate wrapping as negative displacement:
        let dx: i16 = if max_x > x_limit { max_x as i16 - columns as i16 } else { max_x as i16 };
        let dy: i16 = if max_y > y_limit { max_y as i16 - rows as i16 } else { max_y as i16 };
        (dx, dy)
    }

    fn invert_fft(data: &mut SampleBuf) {
        //swap real and imaginary components
        Self::swap_complex_components(data);

        //now forward-transform
        microfft::complex::cfft_4096(data);

        //reswap real and imaginary components
        Self::swap_complex_components(data);
    }

    fn swap_complex_components(data: &mut [Complex32]) {
        //swap real and imaginary components
        for i in 0..data.len() {
            let tmp = data[i].re;
            data[i].re = data[i].im;
            data[i].im = tmp;
        }
    }

    fn complex_conjugate(data: &mut SampleBuf) {
        for i in 0..data.len() {
            data[i] = data[i].conj();
        }
    }


    fn fill_u8_samples_to_complex(frame: &[u8], output: &mut SampleBuf) {
        let n = frame.len();
        // normalize the 8-bit samples to 0..1 range
        for i in 0..n {
            output[i].re = (frame[i] as f32) / 255.0;
            output[i].im = 0f32;
        }
    }

    /// Normalized complex Hadamard product
    fn complex_hadamard_product_norm(
        f0: &SampleBuf,
        f1: &SampleBuf,
        result: &mut SampleBuf
    )
    {
        let n = f0.len();
        for i in 0..n {
            let prod = f0[i] * f1[i];
            let mag = prod.norm_sqr().sqrt(); //TODO use prod.norm when available?
            result[i] = prod / mag;
            // println!("{:?} * {:?} = {:?}", f0[i], f1[i], prod);
        }

    }

}
