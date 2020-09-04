/*
Copyright (c) 2020 Todd Stellanova
LICENSE: BSD3 (see LICENSE file)
*/
#![cfg_attr(not(test), no_std)]

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

pub mod fwht;
