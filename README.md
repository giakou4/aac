<p align="center">
  <img src="https://www.provideocoalition.com/wp-content/uploads/AAC-recortado.jpg">
</p>

[![made-with-matlab](https://img.shields.io/badge/Made%20with-MATLAB-cb6015)](https://www.mathworks.com/products/matlab.html)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/giakou4/game_theory_spatial_PD/LICENSE)

# MULTIMEDIA SYSTEMS

## 1. Advanced Audio Encoder 

The work aims to implement an Advanced Audio Coding (AAC) encoder/decoder. Variations of the AAC are used by many international standards such as MPEG-2, MPEG-4, H.264 etc. The version presented in here is more like the 3GPP TS 26.403 specification where some processing steps are missing. An exception is the psychoacoustic model, which is a slightly simplified version of MPEG AAC. Despite the simplifications, this version leads to quite good results. AAC encoding and decoding belongs to the waveform compression category and attempts to represent the original signal in such a way that its decoded version sounds as similar as possible to the original signal. The psychoacoustic model that allows the introduction of distortions of the signal (noise due to quantization) which is below the audibility threshold is used as a fidelity criterion. For this reason the Psychoacoustic Model mechanism that guides the Quantizer mechanism plays a dominant role. To reduce the excess information, AAC basically uses the transformation coding approach implemented using the so-called Modified Discrete Cosine Transform (MDCT) at the Filterbank stage, while for entropy coding it uses the Huffman coding implemented at the homonymous stage. More specifically, during coding the original audio signal (for us stereo with sampling 48000 samples/sec) is divided into 50% overlapping sections (frames) of 2048 samples. Then each frame is encoded autonomously and therefore the finally encoded bitstream consists of the sequence of bits that correspond to the successive frames. To validate the functionality, AAC is applied to the song Licor De Calandraca.

## 2. Contents   

### Level 1
* AACoder1.m
* iAACoder1.m
* demoAAC1.m
* SSC.m
* filterbank.m
* iFilterbank.m

### Level 2
* AACoder2.m
* iAACoder2.m
* demoAAC2.m
* TNS.m
* iTNS.m
* TableB219.mat

### Level 3
* AACoder3.m
* iAACoder3.m
* demoAAC3.m
* psycho.m
* AACquantizer.m
* iAACquantizer.m
* TableB219.mat
* loadLUT.m
* encodeHuff.m
* decodeHuff.m
* huffCodebooks.mat
* huffCodebookSF.mat

### Report
* report.pdf

### Notes
* Note 1: For Level 2 and 3, the folder are included in code with with command: addpath.
* Note 2: Scripts work only in MATLAB 2020b

## Support

Reach out to me:
- [giakou4's email](mailto:giakonick98@gmail.com "giakonick98@gmail.com")

## License
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/giakou4/multimedia/LICENSE)
