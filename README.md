# LQR-integral-Q-learning

This repository provides the open source code for reproducing the simulation results (Fig 1(a),(d)) presented in the publication:

"Lee, J.Y., Park, J.B., and Choi, Y.H., Integral Q-learning and explorized policy iteration for adaptive optimal control of continuous-time linear systems, Automatica 11(48), pp. 2850~2859, 2012."

To reproduce the results in the paper, please run the code as follows (tested in MATLAB R2012a (32bit) Edition).

1. Set the MATLAB working directory to the cloned local repository path in your machine;
    
2. Clear the environment using the following commands:
``` octave-workspace 	
	close all
	clear all
	clc
```

3. Run 
``` octave-workspace 	
	main.m
```

4. If necessary, change the hyper-parameters in `main.m`.

For a bug report or any issue, please send an e-mail to [Dr. Jaeyoung Lee](mailto:jyounglee@yonsei.ac.kr?subject=[GitHub:%20LQR-integral-Q-learning]%20Bug%20Report%20or%20Any%20Issues).
