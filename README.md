# Lorenz Attractor ANN

### Description

Create a ANN model to predict the Lorenz system. It is nonlinear, non-periodic, three-dimensional and deterministic system given by:

dX/dt &= \sigma(Y-X)
	
dY/dt &= -XZ+\rho X-Y
	
dZ/dt &= XY-\beta Z


A perceptron network and a RNN network was created. The RNN suprosingly performed worse and took alot more epochs for training than the RNN did for a similar prediction accuracy.

### Next steps

There are problems and next steps this project could be taken.
* Exploring the RNN for faster training times and implimenting a scikit learn hyperparameter choosing workflow could be implimented. 
* Implimenting a multiprocessing framework for training the RNN.
* Evaluating RNN on different initial states and physical parameters just the the perceptron model.
* **Generating a bunch of different attractors and giving the model a much bigger training set with more variation.**
* General Code cleanup for more intuitive flow and efficiency

### References

[1] M. Madondo, T. Gibbons, Learning and Modeling Chaos Using LSTM Recurrent Neural Networks, MICS 2018 Proc. (2018) 1–14. http://micsymposium.org/mics2018/proceedings/MICS_2018_paper_26.pdf.

[2] C. Hernandez, A. Martinez, L.F. Mingo, J. Castellanos, Controlling Lorenz Chaos with Neural Networks, (n.d.) 6–9.

[3] L. Zhang, Artificial neural networks model design of Lorenz chaotic system for EEG pattern recognition and prediction, 2017 IEEE Life Sci. Conf. LSC 2017. 2018-January (2018) 39–42. https://doi.org/10.1109/LSC.2017.8268138.

[4] N.Q. Ann, D. Pebrianti, M.F. Abas, L. Bayuaji, M. Syafrullah, Parameter prediction for lorenz attractor by using deep neural network, Indones. J. Electr. Eng. Informatics. 8 (2020) 532–540. https://doi.org/10.11591/ijeei.v8i3.1272.

[5] R.J. Frank, N. Davey, S.P. Hunt, Time series prediction and neural networks, J. Intell. Robot. Syst. Theory Appl. 31 (2001) 91–103. https://doi.org/10.1023/A:1012074215150.