Three well-known Sparse Recovery Algorithms implemented in C/C++ and MATLAB:


      I. Orthogonal Matching Pursuit ( OMP ),
      II. Iterative Hard Thresholding ( IHT ),  
      III. Approximage Message Passing ( AMP )  
# Introduction 

      OMP is a greedy algorithm introduced as an extension to the
      well-established Matching Pursuit algorithm [17]. The OMP
      algorithm iteratively finds the best matrix columns that
      correspond to the non-zero coefficients of the sparse signal, and
      then performs a least squares (LS) optimization in the subspace
      formed from current and previously selected columns.
      AMP and the IHT algorithms do not require LS in each
      iteration, and instead perform simple vector truncation, which
      results in an iterative completion of the sparse signal. The
      parameters, such as step size and threshold, are critical in the
      performance of AMP and IHT algorithms. The optimum
      parameters of AMP are chosen based on the experimental results
      of [13] . In the case of the IHT algorithm, a different flavor of
      algorithm called Normalized-IHT was implemented, where the
      step size is automatically determined in each iteration.
      Important excerpts from the paper:

# Results

      For fair comparison of the algorithms, one must take into
      consideration the window length of the signal, the sparsity
      degree of the test signal, the hyperparameters of the algorithms,
      such as termination criteria and desired reconstruction quality.
      From Fig. 4 [following figure], we observed that the OMP algorithm is fastest in
      reconstruction, whereas the AMP and IHT algorithms that are
      known to be computationally cheaper, appear to be slower. 
      
 ![alt text]( https://github.com/NeuroFan/Compressive_Sensing/blob/master/performance_comparison.png)

      
      This is due to the low sparsity degree and short signal length. The
      OMP algorithm gives better performance for less sparse signals
      [27], and here the experiments were done with signals with less
      than 10% occupancy. The IHT algorithm’s performance is
      relatively independent from the sparsity degree, and the
      performance of AMP is less sensitive to sparsity degree than the
      OMP [27], [28]. These issues are rather strong practical
      arguments for flexible designs. 


#Citation 

If you used the code please cite our paper [1].

You can find the paper in this repository or through following link:
https://github.com/NeuroFan/Compressive_Sensing/blob/master/An%20Embedded%20Programmable%20Processor%20for%20Compressive.pdf

*Reference* 

[1] M. Safarpour, I. Hautala and O. Silvén, "An Embedded Programmable Processor for Compressive Sensing Applications," 2018 IEEE Nordic Circuits and Systems Conference (NORCAS): NORCHIP and International Symposium of System-on-Chip (SoC), Tallinn, Estonia, 2018, pp. 1-5.

doi: 10.1109/NORCHIP.2018.8573494
