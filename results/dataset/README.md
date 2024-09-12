## About
Dataset accompanying the paper "_Adaptive connectivity control in networked multi-agent systems: A distributed approach_" by M. Krizmancic and S. Bogdan submitted to PLOS ONE journal on April 30, 2024.

## Content
- Vector images of the figures presented in the paper.
- Data files containing the values used to build the figures.

### List of figures and associated files
| **Figure num.** | **Vector image**                                             | **Data file**                                  |
|-----------------|--------------------------------------------------------------|------------------------------------------------|
| Fig 3.          | 6line_adapt_sigma_Lambda.svg<br>6line_fixed_sigma_Lambda.svg | 6line_adapt_sigma.mat<br>6line_fixed_sigma.mat |
| Fig 4.          | new_initial_topology.svg                                     | new_*.mat                                      |
| Fig 5.          | new_estimate_real_Lambda.svg                                 | new_estimate_real.mat                          |
| Fig 6.          | new_follow_Lambda.svg                                        | new_follow.mat                                 |
| Fig 7.          | new_follow_Snapshots.svg                                     | new_follow.mat                                 |
| Fig 8.          | new_fail_Lambda.svg                                          | new_fail.mat                                   |
| Fig 9.          | new_energy_Lambda.svg                                        | new_energy.mat                                 |
| Fig 10.         | 15_initial_topology.svg                                      | 15_fail.mat                                    |
| Fig 11.         | 15_fail_Lambda.svg                                           | 15_fail.mat                                    |

### List of variables in the data files
| **Variable**      | **Type**          | **Description**                           |
|-------------------|-------------------|-------------------------------------------|
| A                 | 4-D double        | History of A matrices for each agent      |
| A0                | n-by-n double     | A matrix of the initial topology          |
| K_l2              | n-by-steps double | History of K_l2 limits for each agent     |
| estimation_period | double            | First steps without control               |
| lambda            | n-by-steps double | History of lambda_2 values for each agent |
| l2_refs           | 2-D double        | Lambda_2 reference signal                 |
| params            | struct            | Simulation parameters                     |
| rand_stream       | <variable>        | Type of random stream used in simulation  |
| signal_refs       | cell              | Signals for simulating signal strength    |
| steps             | double            | Total number of steps in the simulation   |

## Usage
1. Unzip the folder to your desired location.
2. Open the data files in MATLAB and explore the data.
    ```matlab
    load('<data_file_name>.mat');
    disp(lambda);
    ```
    -- OR --   
    Open the data files in Python and explore the data.
    ```python
    import scipy.io as sio
    mat_contents = sio.loadmat('<data_file_name>.mat')
    print(mat_contents['lambda'])
    ```
3. For more information, refer to the paper and the documentation of your preferred tool ([Matlab](https://www.mathworks.com/help/matlab/ref/load.html), [Python](https://docs.scipy.org/doc/scipy/reference/generated/scipy.io.loadmat.html))