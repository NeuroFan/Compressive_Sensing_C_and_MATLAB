Three well-known Sparse Recovery Algorithms implemented in C and MATLAB:


      I. Orthogonal Matching Persuits (OMP),
      II. Iterative Hard Thresholding (IHT) and 
      III. Approximage Message Passing (AMP)  

If you used the code please cite our paper [1].

Important excerpts from the paper:


      For fair comparison of the algorithms, one must take into
      consideration the window length of the signal, the sparsity
      degree of the test signal, the hyperparameters of the algorithms,
      such as termination criteria and desired reconstruction quality.
      From Fig. 4 [following figure], we observed that the OMP algorithm is fastest in
      reconstruction, whereas the AMP and IHT algorithms that are
      known to be computationally cheaper, appear to be slower. This
      is due to the low sparsity degree and short signal length. The
      OMP algorithm gives better performance for less sparse signals
      [27], and here the experiments were done with signals with less
      than 10% occupancy. The IHT algorithm’s performance is
      relatively independent from the sparsity degree, and the
      performance of AMP is less sensitive to sparsity degree than the
      OMP [27], [28]. These issues are rather strong practical
      arguments for flexible designs. 

  ![alt text]( https://github.com/NeuroFan/Compressive_Sensing/blob/master/performance_comparison.png)



*Reference* 

[1] M. Safarpour, I. Hautala and O. Silvén, "An Embedded Programmable Processor for Compressive Sensing Applications," 2018 IEEE Nordic Circuits and Systems Conference (NORCAS): NORCHIP and International Symposium of System-on-Chip (SoC), Tallinn, Estonia, 2018, pp. 1-5.

doi: 10.1109/NORCHIP.2018.8573494
