Updated 07/20/2020                                                                                 |

Hello and welcome to my LabVIEW directory!

Continuous Voltage Acquisition and Logging: This code automates an electro-polling experiment
	designed to simulate the electrical signals emitted by neurons. The code commands the
	potentiostat to stimulate the sample in a pattern based on user input, recording the
	resulting behavior by automating measurements from a DAQ and photodiode.

Spectral Integration Threshold Ranging: This code automates spectrometer measurements, displaying
	the spectrum in real time and calculating the integral. The program was made robust to
	accomodate varying initial spectra with a startup state that lets the user check the
	spectra and adjust parameters before starting the experiment.

Spectral Integration Fast Version: This code was adapted from a previous version of the above code,
	removing all unessential components to test the mimimum spectrometer integration time
	supported by the LabVIEW interface. The result was that the photodiode is significantly
	faster and a better candidate for some experiments despite its inability to measure
	wavelength.