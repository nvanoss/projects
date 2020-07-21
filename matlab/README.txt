Updated 07/21/2020                                                                                 |

Hello and welcome to my MATLAB directory!

SAW_Streaming_Simulation: I simulated Surface Acoustic Wave streaming including Leaky SAW
	propagation and 3D attenuation.

VNA_S11_data_processing: I wrote this script to extract scattering parameter data at key
	frequencies, then plot and compare values across data from multiple files to distinguish
	more optimal Interdigital Transducer designs from less promising devices.

nvofitting_centroid: this code calculates polynomial fittings of spectra up to the fifth order.
	It then uses these fittings to calculate the centroid of spectra, comparing and tracking
	changes across time, sample concentration, and between control groups. This code was used to
	quantify the results presented in an HIV paper to be published soon in a reputable journal.
	This file was purposefully left separate from the following similar file to limit the
	lengthy runtime	of the script dealing with large dimensions in a five-dimension array.

nvofitting_specint: this code calculates polynomial fittings of spectra up to the fifth order.
	It then uses these fittings to calculate the integral of spectra, comparing and tracking
	changes	across time, sample concentration, and between control groups. This code was used to
	quantify the results presented in an HIV paper to be published soon in a reputable journal.

nvo_image_intensity_extraction: this simple script converts the intensity of fluorescence images
	into numerical values and compares them, tracking the shift over time. It additionally
	verifies that the image is entirely composed of green light, detecting any leakage in the
	optical setup.

nvo_image_intensity_extraction_exposure_time_test: I altered the above script to help PhD students
	test whether an increase in the exposure time of the camera corresponded linearly to the
	increase in intensity. Note that many odd jobs such as this one were requested of me which
	resulted in additional MATLAB files that are not contained in this repository for brevity
	but can be supplied to an employer by request.

